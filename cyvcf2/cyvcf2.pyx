#cython: profile=False
from __future__ import print_function
import os
import os.path as op
import sys
from collections import defaultdict
import atexit
import tempfile
import numpy as np
from array import array
import math

from libc cimport stdlib
cimport numpy as np
np.seterr(invalid='ignore')
np.import_array()
import locale

ENC = locale.getpreferredencoding()

from cython cimport view

from cpython.version cimport PY_MAJOR_VERSION

# overcome lack of __file__ in cython
import inspect
if not hasattr(sys.modules[__name__], '__file__'):
    __file__ = inspect.getfile(inspect.currentframe())



def par_relatedness(vcf_path, samples, ncpus, sites, min_depth=5, each=1):
    from multiprocessing import Pool
    p = Pool(ncpus)
    cdef np.ndarray aibs, an, ahets

    for i, fname in enumerate(p.imap_unordered(_par_relatedness, [
        (vcf_path, samples, min_depth, i, ncpus, each, sites) for i in range(ncpus)])):

        arrays = np.load(fname)
        os.unlink(fname)
        if i == 0:
            aibs, an, ahets = arrays['ibs'], arrays['n'], arrays['hets']
        else:
            aibs += arrays['ibs']
            an += arrays['n']
            ahets += arrays['hets']

    return VCF(vcf_path, samples=samples)._relatedness_finish(aibs, an, ahets)


def par_het(vcf_path, samples, ncpus, qsites, min_depth=8, percentiles=(10, 90),
                  int each=1, int offset=0):
    from multiprocessing import Pool
    p = Pool(ncpus)

    any_counts, sum_counts, het_counts = 0, 0, 0
    all_gt_types, mean_depths, sites = [], [], []
    maf_lists = defaultdict(list)
    for ret in p.imap_unordered(_par_het, [(vcf_path, samples, qsites, min_depth, i, ncpus, each) for i in range(ncpus)]):
        (mean_depths_, maf_lists_, het_counts_, sum_counts_,
                all_gt_types_, sites_, any_counts_) = ret
        mean_depths.extend(mean_depths_)
        any_counts += any_counts_
        sites.extend(sites_)
        all_gt_types.extend(all_gt_types_)
        sum_counts += sum_counts_ # an array
        het_counts += het_counts_ # an array
        for k in maf_lists_:
            li = maf_lists_[k]
            maf_lists[k].extend(li)
    mean_depths = np.array(mean_depths, dtype=np.int32).T
    vcf = VCF(vcf_path, samples=samples, gts012=True)
    return vcf._finish_het(mean_depths, maf_lists,
                           percentiles,
                           het_counts, sum_counts, all_gt_types,
                           sites, any_counts)


def _par_het(args):
    vcf_path, samples, sites, min_depth, offset, ncpus, each = args
    each *= ncpus
    vcf = VCF(vcf_path, samples=samples, gts012=True)
    return vcf.het_check(min_depth=min_depth, sites=sites, each=each, offset=offset, _finish=False)


def _par_relatedness(args):
    vcf_path, samples, min_depth, offset, ncpus, each, sites = args
    vcf = VCF(vcf_path, samples=samples, gts012=True)
    each = each * ncpus
    vibs, vn, vhet = vcf._site_relatedness(min_depth=min_depth, offset=offset, each=each, sites=sites)
    # to get around limits of multiprocessing size of transmitted data, we save
    # the arrays to disk and return the file
    fname = tempfile.mktemp(suffix=".npz")
    atexit.register(os.unlink, fname)
    np.savez_compressed(fname, ibs=np.asarray(vibs), hets=np.asarray(vhet), n=np.asarray(vn))
    return fname

cdef unicode xstr(s):
    if type(s) is unicode:
        # fast path for most common case(s)
        return <unicode>s
    elif PY_MAJOR_VERSION < 3 and isinstance(s, bytes):
        # only accept byte strings in Python 2.x, not in Py3
        return (<bytes>s).decode('ascii')
    elif isinstance(s, unicode):
        # an evil cast to <unicode> might work here in some(!) cases,
        # depending on what the further processing does.  to be safe,
        # we can always create a copy instead
        return unicode(s)
    else:
        raise TypeError(...)

def r_(int[::view.contiguous] a_gts, int[::view.contiguous] b_gts, float f, int32_t n_samples):
    return r_unphased(&a_gts[0], &b_gts[0], f, n_samples)

cdef set_constants(VCF v):
    v.HOM_REF = 0
    v.HET = 1
    if v.gts012:
        v.HOM_ALT = 2
        v.UNKNOWN = 3
    else:
        v.UNKNOWN = 2
        v.HOM_ALT = 3

cdef class VCF:

    cdef htsFile *hts
    cdef const bcf_hdr_t *hdr
    cdef tbx_t *idx
    cdef hts_idx_t *hidx
    cdef int n_samples
    cdef int PASS
    cdef bytes fname
    cdef bint gts012
    cdef bint lazy
    cdef list _seqnames
    # holds a lookup of format field -> type.
    cdef dict format_types

    cdef readonly int HOM_REF
    cdef readonly int HET
    cdef readonly int HOM_ALT
    cdef readonly int UNKNOWN

    def __init__(self, fname, mode="r", gts012=False, lazy=False, samples=None):
        if fname == b"-" or fname == "-":
            fname = b"/dev/stdin"
        if not op.exists(fname):
            raise Exception("bad path: %s" % fname)
        fname, mode = to_bytes(fname), to_bytes(mode)
        self.hts = hts_open(fname, mode)
        if self.hts == NULL:
            raise IOError("Error opening %s" % fname)
        if self.hts.format.format != vcf and self.hts.format.format != bcf:
            raise IOError("%s if not valid bcf or vcf" % fname)

        cdef bcf_hdr_t *hdr
        hdr = self.hdr = bcf_hdr_read(self.hts)
        if samples is not None:
            self.set_samples(samples)
        self.n_samples = bcf_hdr_nsamples(self.hdr)
        self.PASS = -1
        self.fname = to_bytes(fname)
        self.gts012 = gts012
        self.lazy = lazy
        self._seqnames = []
        set_constants(self)
        self.format_types = {}

    cdef get_type(self, fmt):
        fmt = from_bytes(fmt)
        if not fmt in self.format_types:
            s = self[fmt]
            self.format_types[fmt] = s["Type"]

        return from_bytes(self.format_types[fmt])

    def add_to_header(self, line):
        ret = bcf_hdr_append(self.hdr, to_bytes(line))
        if ret != 0:
            raise Exception("couldn't add '%s' to header")
        ret = bcf_hdr_sync(self.hdr)
        if ret != 0:
            raise Exception("couldn't add '%s' to header")
        return ret

    def add_info_to_header(self, adict):
        return self.add_to_header("##INFO=<ID={ID},Number={Number},Type={Type},Description=\"{Description}\">".format(**adict))

    def add_format_to_header(self, adict):
        return self.add_to_header("##FORMAT=<ID={ID},Number={Number},Type={Type},Description=\"{Description}\">".format(**adict))

    def add_filter_to_header(self, adict):
        return self.add_to_header("##FILTER=<ID={ID},Description=\"{Description}\">".format(**adict))

    def set_samples(self, samples):
        if samples is None:
            samples = "-".encode()
        if isinstance(samples, list):
            samples = to_bytes(",".join(samples))
        else:
            samples = to_bytes(samples)

        ret = bcf_hdr_set_samples(self.hdr, <const char *>samples, 0)
        assert ret >= 0, ("error setting samples", ret)
        if ret != 0 and samples != "-":
            s = samples.split(",")
            if ret < len(s):
                sys.stderr.write("warning: not all samples in PED found in VCF\n")

    def update(self, id, type, number, description):
        ret = bcf_hdr_append(self.hdr, "##INFO=<ID={id},Number={number},Type={type},Description=\"{description}\">".format(id=id, type=type, number=number, description=description))
        if ret != 0:
            raise Exception("unable to update to header: %d", ret)
        ret = bcf_hdr_sync(self.hdr)
        if ret != 0:
            raise Exception("unable to update to header")

    def _bcf_region(VCF self, region):
        if self.hidx == NULL:
            self.hidx = bcf_index_load(self.fname)
        assert self.hidx != NULL, ("error loading .csi index for %s" % self.fname)
        cdef bcf1_t *b
        cdef int ret
        cdef hts_itr_t *itr

        itr = bcf_itr_querys(self.hidx, self.hdr, to_bytes(region))
        if itr == NULL:
            sys.stderr.write("no intervals found for %s at %s\n" % (self.fname, region))
            raise StopIteration
        try:
            while True:
                b = bcf_init()
                ret = bcf_itr_next(self.hts, itr, b)
                if ret < 0:
                    bcf_destroy(b)
                    break
                yield newVariant(b, self)
        finally:
            if itr != NULL:
                hts_itr_destroy(itr)


    def __call__(VCF self, region=None):
        if not region:
            yield from self
            raise StopIteration

        if self.fname.decode(ENC).endswith('.bcf'):
            yield from self._bcf_region(region)
            raise StopIteration

        if self.idx == NULL:
            if not (op.exists(from_bytes(self.fname)+ ".tbi") or
                    op.exists(from_bytes(self.fname) + ".csi")):
                raise Exception("can't extract region without tabix or csi index for %s" % self.fname)


            self.idx = tbx_index_load(to_bytes(self.fname))
            assert self.idx != NULL, "error loading tabix index for %s" % self.fname

        cdef hts_itr_t *itr
        cdef kstring_t s
        cdef bcf1_t *b
        cdef int slen, ret

        itr = tbx_itr_querys(self.idx, to_bytes(region))

        if itr == NULL:
            sys.stderr.write("no intervals found for %s at %s\n" % (self.fname, region))
            raise StopIteration

        try:
            slen = tbx_itr_next(self.hts, self.idx, itr, &s)
            while slen > 0:
                b = bcf_init()
                ret = vcf_parse(&s, self.hdr, b)
                if ret > 0:
                    bcf_destroy(b)
                    raise Exception("error parsing")
                yield newVariant(b, self)
                slen = tbx_itr_next(self.hts, self.idx, itr, &s)
        finally:
            stdlib.free(s.s)
            hts_itr_destroy(itr)

    def header_iter(self):
        cdef int i
        for i in range(self.hdr.nhrec):
            yield newHREC(self.hdr.hrec[i], self.hdr)

    def ibd(self, int nmax=-1):
        assert self.gts012
        import itertools

        cdef int i, rl, n_bins = 16

        samples = self.samples
        sample_to_idx = {s: samples.index(s) for s in samples}
        sample_pairs = list(itertools.combinations(samples, 2))
        # values of bins, run_length

        cdef int n = 0
        cdef float pi
        cdef int[:] b
        cdef int[:] gts
        cdef int idx0, idx1
        bins = np.zeros((len(sample_pairs), n_bins), dtype=np.int32)
        rls = np.zeros(len(sample_pairs), dtype=np.int32)

        for v in self:
            if n == nmax: break
            n += 1
            gts = v.gt_types
            pi = v.aaf
            for i, (s0, s1) in enumerate(sample_pairs):
                b = bins[i, :]
                idx0, idx1 = sample_to_idx[s0], sample_to_idx[s1]
                rls[i] = ibd(gts[idx0], gts[idx1], rls[i], pi, &b[0], n_bins)

        return {sample_pairs[i]: bins[i, :] for i in range(len(sample_pairs))}

    # pull something out of the HEADER, e.g. CSQ
    def __getitem__(self, key):
        key = to_bytes(key)
        cdef bcf_hrec_t *b = bcf_hdr_get_hrec(self.hdr, BCF_HL_INFO, b"ID", key, NULL);
        cdef int i
        if b == NULL:
            b = bcf_hdr_get_hrec(self.hdr, BCF_HL_FMT, b"ID", key, NULL);
        if b == NULL:
            b = bcf_hdr_get_hrec(self.hdr, BCF_HL_GEN, key, NULL, NULL);
            if b == NULL:
                raise KeyError(key)
            d = {from_bytes(b.key): from_bytes(b.value)}
        else:
            d =  {from_bytes(b.keys[i]): from_bytes(b.vals[i]) for i in range(b.nkeys)}
        #bcf_hrec_destroy(b)
        return d

    def __contains__(self, key):
        try:
            self[key]
            return True
        except KeyError:
            return False

    contains = __contains__


    def __dealloc__(self):
        if self.hdr != NULL:
            bcf_hdr_destroy(self.hdr)
            self.hdr = NULL
        if self.hts != NULL:
            hts_close(self.hts)
            self.hts = NULL
        if self.idx != NULL:
            tbx_destroy(self.idx)
        if self.hidx != NULL:
            hts_idx_destroy(self.hidx)

    def __iter__(self):
        return self

    def __next__(self):

        cdef bcf1_t *b = bcf_init()
        cdef int ret
        with nogil:
            ret = bcf_read(self.hts, self.hdr, b)
        if ret >= 0:
            return newVariant(b, self)
        else:
            bcf_destroy(b)
        raise StopIteration

    property samples:
        def __get__(self):
            cdef int i
            return [str(self.hdr.samples[i].decode('utf-8')) for i in range(self.n_samples)]

    property raw_header:
        def __get__(self):
            cdef int hlen
            s = bcf_hdr_fmt_text(self.hdr, 0, &hlen)
            return from_bytes(s)

    property seqnames:
        def __get__(self):
            if len(self._seqnames) > 0: return self._seqnames
            cdef char **cnames
            cdef int i, n = 0
            if self.fname.decode(ENC).endswith('.bcf'):
                if self.hidx == NULL:
                    self.hidx = bcf_index_load(self.fname)
                if self.hidx != NULL:
                    cnames = bcf_index_seqnames(self.hidx, self.hdr, &n)
            else:
                if self.idx == NULL and (op.exists(from_bytes(self.fname)+ ".tbi") or
                        op.exists(from_bytes(self.fname) + ".csi")):
                    self.idx = tbx_index_load(to_bytes(self.fname))
                if self.idx !=NULL:
                    cnames = tbx_seqnames(self.idx, &n)
            if n == 0:
                cnames = bcf_hdr_seqnames(self.hdr, &n)

            self._seqnames = [cnames[i].decode() for i in range(n)]
            stdlib.free(cnames)
            return self._seqnames

    def plot_relatedness(self, riter):
        import pandas as pd
        from matplotlib import pyplot as plt
        from matplotlib import gridspec
        import seaborn as sns
        sns.set_style("ticks")

        df = []
        for row in riter:
          row['jtags'] = '|'.join(row['tags'])
          df.append(row)


        df = pd.DataFrame(df)
        fig = plt.figure(figsize=(9, 9))

        gs = gridspec.GridSpec(2, 1, height_ratios=[3.5, 1])

        ax0, ax1 = plt.subplot(gs[0]), plt.subplot(gs[1])

        if "error" in df.columns:
            # plot all gray except points that don't match our expectation.
            import matplotlib
            matplotlib.rcParams['pdf.fonttype'] = 42
            import matplotlib.colors as mc
            colors = [mc.hex2color(h) for h in ('#b6b6b6', '#ff3333')]
            for i, err in enumerate(("ok", "error")):
                subset = df[df.error == err]
                subset.plot(kind='scatter', x='rel', y='ibs0', c=colors[i],
                          edgecolor=colors[0],
                          label=err, ax=ax0, s=17 if i == 0 else 35)
            sub = df[df.error == "error"]
            for i, row in sub.iterrows():
                ax0.annotate(row['sample_a'] + "\n" + row['sample_b'],
                        (row['rel'], row['ibs0']), fontsize=8)
        else:
            # color by the relation derived from the genotypes.
            colors = sns.color_palette("Set1", len(set(df.jtags)))
            for i, tag in enumerate(set(df.jtags)):
                subset = df[df.jtags == tag]
                subset.plot(kind='scatter', x='rel', y='ibs0', c=colors[i],
                          label=tag, ax=ax0)

            ax0.legend()

        ax0.set_ylim(ymin=0)
        ax0.set_xlim(xmin=df.rel.min())

        ax1.set_xlim(*ax0.get_xlim())
        ax1.hist(df.rel, 40)
        ax1.set_yscale('log', nonposy='clip')
        return fig

    def gen_variants(self, sites,
                    offset=0, each=1, call_rate=0.8):

        seqnames = set(self.seqnames)

        if sites is not None:
            if isinstance(sites, basestring):
                isites = []
                for i in (x.strip().split(":") for x in open(sites)):
                    # handle 'chr' prefix on incoming.
                    if not i[0] in seqnames and 'chr' + i[0] in seqnames:
                        i[0] = 'chr' + i[0]
                    i[1] = int(i[1])
                    isites.append(i)
            else:
                isites = sites

        cdef Variant v
        cdef int k, last_pos
        if sites:
            isites = isites[offset::each]
            def gen():
                ref, alt = None, None
                j = 0
                for i, osite in enumerate(isites):
                    if len(osite) >= 4:
                        chrom, pos, ref, alt = osite[:4]
                    else:
                        chrom, pos = osite[:2]
                    for v in self("%s:%s-%s" % (chrom, pos, pos)):
                        if len(v.ALT) != 1: continue
                        if ref is not None:
                            if v.REF != ref: continue
                            if alt is not None:
                                if v.ALT[0] != alt: continue
                        if v.call_rate < call_rate: continue
                        yield i, v
                        j += 1
                        break
        else:
            def gen():
                last_pos, k = -10000, 0
                for v in self:
                    if abs(v.POS - last_pos) < 5000: continue
                    if len(v.REF) != 1: continue
                    if len(v.ALT) != 1: continue
                    if v.call_rate < 0.5: continue
                    if not 0.03 < v.aaf < 0.6: continue
                    if np.mean(v.gt_depths > 7) < 0.5: continue
                    last_pos = v.POS
                    if k >= offset and k % each == 0:
                        yield k, v
                    k += 1
                    if k > 20000: break
        return gen

    def het_check(self, min_depth=8, percentiles=(10, 90), _finish=True,
                  int each=1, int offset=0,
                  sites=None):

        cdef int i, k, n_samples = len(self.samples)
        cdef Variant v
        cdef np.ndarray het_counts = np.zeros((n_samples,), dtype=np.int32)

        cdef np.ndarray sum_depths = np.zeros((n_samples,), dtype=np.int32)
        cdef np.ndarray sum_counts = np.zeros((n_samples,), dtype=np.int32)
        cdef int any_counts = 0

        # keep the sites and gts that we used for PCA
        used_sites, all_gt_types = [], []

        mean_depths = []

        gen = self.gen_variants(sites, each=each, offset=offset)
        maf_lists = defaultdict(list)
        idxs = np.arange(n_samples)
        for i, v in gen():
            if v.CHROM in ('X', 'chrX'): break
            if v.aaf < 0.01: continue
            if v.call_rate < 0.5: continue
            used_sites.append("%s:%d:%s:%s" % (v.CHROM, v.start + 1, v.REF, v.ALT[0]))
            alts = v.gt_alt_depths
            assert len(alts) == n_samples
            depths = (alts + v.gt_ref_depths).astype(np.int32)
            sum_depths += depths
            sum_counts += (depths > min_depth)
            any_counts += 1
            mean_depths.append(depths)

            mafs = alts / depths.astype(float)
            gt_types = v.gt_types
            hets = gt_types == 1
            het_counts[hets] += 1
            for k in idxs[hets]:
                if depths[k] <= min_depth: continue
                maf_lists[k].append(mafs[k])
            all_gt_types.append(np.array(gt_types, dtype=np.uint8))

        if _finish:
            mean_depths = np.array(mean_depths, dtype=np.int32).T
            return self._finish_het(mean_depths, maf_lists,
                                    percentiles,
                                    het_counts, sum_counts, all_gt_types,
                                    used_sites, any_counts)
        return (mean_depths, maf_lists, het_counts, sum_counts,
                all_gt_types, used_sites, any_counts)


    def _finish_het(self, mean_depths, maf_lists, percentiles, het_counts,
            sum_counts, all_gt_types, sites, any_counts):

            sample_ranges = {}
            for i, sample in enumerate(self.samples):
                qs = np.asarray(np.percentile(maf_lists[i] or [0], percentiles))
                sample_ranges[sample] = dict(zip(['p' + str(p) for p in percentiles], qs))
                sample_ranges[sample]['range'] = qs.max() - qs.min()
                sample_ranges[sample]['het_ratio'] = het_counts[i] / float(any_counts)
                sample_ranges[sample]['het_count'] = het_counts[i]
                sample_ranges[sample]['sampled_sites'] = sum_counts[i]
                sample_ranges[sample]['mean_depth'] = np.mean(mean_depths[i])
                sample_ranges[sample]['median_depth'] = np.median(mean_depths[i])
            # used for the peddy paper.
            if os.environ.get('PEDDY_MAF_DUMP'):
                path = os.environ['PEDDY_MAF_DUMP']
                import cPickle
                cPickle.dump({s: maf_lists[i] for i, s in enumerate(self.samples)}, open(path, 'wb', -1))


            return sample_ranges, sites, np.transpose(all_gt_types)


    def site_relatedness(self, sites=None,
                         min_depth=5, each=1):

        vibs, vn, vhet = self._site_relatedness(sites=sites, min_depth=min_depth, each=each)
        return self._relatedness_finish(vibs, vn, vhet)


    cdef _site_relatedness(self, sites=None,
            int min_depth=5, int each=1, int offset=0):
        """
        sites must be an file of format: chrom:pos1:ref:alt where
        we match on all parts.
        it must have a matching file with a suffix of .bin.gz that is the binary
        genotype data. with 0 == hom_ref, 1 == het, 2 == hom_alt, 3 == unknown.
        min_depth applies per-sample
        """
        cdef int n_samples = len(self.samples)
        cdef int k, i
        assert each >= 0

        gen = self.gen_variants(sites, offset=offset, each=each)

        cdef int32_t[:, ::view.contiguous] ibs = np.zeros((n_samples, n_samples), np.int32)
        cdef int32_t[:, ::view.contiguous] n = np.zeros((n_samples, n_samples), np.int32)
        cdef int32_t[:] hets = np.zeros((n_samples, ), np.int32)
        cdef int32_t[:] gt_types = np.zeros((n_samples, ), np.int32)
        cdef int32_t[:] depths = np.zeros((n_samples, ), np.int32)

        cdef Variant v

        for j, (i, v) in enumerate(gen()):
            gt_types = v.gt_types
            krelated(&gt_types[0], &ibs[0, 0], &n[0, 0], &hets[0], n_samples)

        return ibs, n, hets

    def relatedness(self, int n_variants=35000, int gap=30000, float min_af=0.04,
                    float max_af=0.8, float linkage_max=0.2, min_depth=8):

        cdef Variant v

        cdef int last = -gap, nv = 0, nvt=0
        cdef int *last_gts = NULL
        samples = self.samples
        cdef int n_samples = len(samples)
        cdef float aaf = 0.0
        cdef int n_unlinked = 0

        cdef int32_t[:, ::view.contiguous] ibs = np.zeros((n_samples, n_samples), np.int32)
        cdef int32_t[:, ::view.contiguous] n = np.zeros((n_samples, n_samples), np.int32)
        cdef int32_t[:] hets = np.zeros((n_samples, ), np.int32)
        cdef int32_t[:] gt_types = np.zeros((n_samples, ), np.int32)

        for v in self:
            nvt += 1
            if last_gts == NULL:
                if v._gt_types == NULL:
                    v.gt_types
                last_gts = v._gt_types
            if v.POS - last < gap and v.POS > last:
                continue
            if v.call_rate < 0.5: continue
            # require half of the samples to meet the min depth
            if np.mean(v.gt_depths > min_depth) < 0.5: continue
            aaf = v.aaf
            if aaf < min_af: continue
            if aaf > max_af: continue
            if linkage_max < 1 and v.POS - last < 40000:
                if v._gt_types == NULL:
                    v.gt_types
                # require 5 unlinked variants
                if r_unphased(last_gts, v._gt_types, 1e-5, n_samples) > linkage_max:
                    continue
                n_unlinked += 1
                if n_unlinked < 5:
                    continue

            n_unlinked = 0

            if v._gt_types == NULL:
                gt_types = v.gt_types
            last, last_gts = v.POS, v._gt_types

            v.relatedness(ibs, n, hets)
            nv += 1
            if nv == n_variants:
                break
        sys.stderr.write("tested: %d variants out of %d\n" % (nv, nvt))
        return self._relatedness_finish(ibs, n, hets)

    cdef dict _relatedness_finish(self,
                                  int32_t[:, ::view.contiguous] _ibs,
                                  int32_t[:, ::view.contiguous] _n,
                                  int32_t[:] _hets):
        samples = self.samples
        cdef int sj, sk, ns = len(samples)

        res = {'sample_a': [], 'sample_b': [],
                'rel': array('f'),
                'hets_a': array('I'),
                'hets_b': array('I'),
               'shared_hets': array('I'),
               'ibs0': array('I'),
               'ibs2': array('I'),
               'n': array('I')}

        cdef float bot

        for sj in range(ns):
            sample_j = samples[sj]
            for sk in range(sj, ns):
                if sj == sk: continue
                sample_k = samples[sk]

                # calculate relatedness. we use the geometric mean.
                bot = math.exp(0.5 * (math.log(1 + _hets[sk]) + math.log(1 + _hets[sj])))
                #bot = (_hets[sk] + _hets[sj])
                phi = (_ibs[sk, sj] - 2.0 * _ibs[sj, sk]) / (bot)

                res['sample_a'].append(sample_j)
                res['sample_b'].append(sample_k)
                res['hets_a'].append(_hets[sj])
                res['hets_b'].append(_hets[sk])
                res['rel'].append(phi) # rel is 2*kinship
                res['ibs0'].append(_ibs[sj, sk])
                res['shared_hets'].append(_ibs[sk, sj])
                res['ibs2'].append(_n[sk, sj])
                res['n'].append(_n[sj, sk])
        return res

cdef class Variant(object):
    cdef bcf1_t *b
    cdef VCF vcf
    cdef int *_gt_types
    cdef int *_gt_ref_depths
    cdef int *_gt_alt_depths
    cdef void *fmt_buffer
    cdef int *_gt_phased
    cdef float *_gt_quals
    cdef int *_int_gt_quals
    cdef int *_gt_idxs
    cdef int _gt_nper
    cdef int *_gt_pls
    cdef float *_gt_gls
    cdef readonly INFO INFO
    cdef int _ploidy

    cdef readonly int POS

    def __cinit__(self):
        self.b = NULL
        self._gt_types = NULL
        self._gt_phased = NULL
        self._gt_pls = NULL
        self._ploidy = -1

    def __repr__(self):
        return "Variant(%s:%d %s/%s)" % (self.CHROM, self.POS, self.REF, ",".join(self.ALT))

    def __str__(self):
        cdef kstring_t s
        s.s, s.l, s.m = NULL, 0, 0
        vcf_format(self.vcf.hdr, self.b, &s)
        try:
            return s.s[:s.l].decode()
        finally:
            stdlib.free(s.s)

    def __dealloc__(self):
        if self.b is not NULL:
            bcf_destroy(self.b)
            self.b = NULL
        if self._gt_types != NULL:
            stdlib.free(self._gt_types)
        if self._gt_ref_depths != NULL:
            stdlib.free(self._gt_ref_depths)
        if self._gt_alt_depths != NULL:
            stdlib.free(self._gt_alt_depths)
        if self._gt_phased != NULL:
            stdlib.free(self._gt_phased)
        if self._gt_quals != NULL:
            stdlib.free(self._gt_quals)
        if self._int_gt_quals != NULL:
            stdlib.free(self._int_gt_quals)
        if self._gt_idxs != NULL:
            stdlib.free(self._gt_idxs)
        if self._gt_pls != NULL:
            stdlib.free(self._gt_pls)
        if self._gt_gls != NULL:
            stdlib.free(self._gt_gls)

    property gt_bases:
        def __get__(self):
            cdef np.ndarray gt_types = self.gt_types
            cdef int i, n = self.ploidy, j=0, k
            cdef char **alleles = self.b.d.allele
            #cdef dict d = {i:alleles[i] for i in range(self.b.n_allele)}
            cdef list d = [from_bytes(alleles[i]) for i in range(self.b.n_allele)]
            d.append(".") # -1 gives .
            cdef list a = []
            cdef list phased = list(self.gt_phases)
            cdef list lookup = ["/", "|"]
            cdef int unknown = 3 if self.vcf.gts012 else 2
            for i in range(0, n * self.vcf.n_samples, n):
                if n == 2:
                    if gt_types[j] == unknown:
                        a.append("./.")
                    else:
                        try:
                            d[self._gt_idxs[i+1]]
                            a.append(d[self._gt_idxs[i]] + lookup[phased[j]] + d[self._gt_idxs[i+1]])
                        except IndexError:
                            a.append(d[self._gt_idxs[i]])
                elif n == 1:
                    a.append(d[self._gt_idxs[i]])
                else:
                    raise Exception("gt_bases not implemented for ploidy > 2")

                j += 1
            return np.array(a, np.str)



    def relatedness(self,
                    int32_t[:, ::view.contiguous] ibs,
                    int32_t[:, ::view.contiguous] n,
                    int32_t[:] hets):
        if not self.vcf.gts012:
            raise Exception("must call relatedness with gts012")
        if self._gt_types == NULL:
            self.gt_types
        return krelated(<int32_t *>self._gt_types, &ibs[0, 0], &n[0, 0], &hets[0], self.vcf.n_samples)

    property num_called:
        def __get__(self):
            if self._gt_types == NULL:
                self.gt_types
            cdef int n = 0, i = 0
            if self.vcf.gts012:
                for i in range(self.vcf.n_samples):
                    if self._gt_types[i] != 3:
                        n+=1
            else:
                for i in range(self.vcf.n_samples):
                    if self._gt_types[i] != 2:
                        n+=1
            return n

    property call_rate:
        def __get__(self):
            if self.vcf.n_samples > 0:
                return float(self.num_called) / self.vcf.n_samples

    property aaf:
        def __get__(self):
            num_chroms = 2.0 * self.num_called
            if num_chroms == 0.0:
                return 0.0
            return float(self.num_het + 2 * self.num_hom_alt) / num_chroms

    property nucl_diversity:
        def __get__(self):
            num_chroms = 2.0 * self.num_called
            p = self.aaf
            return (num_chroms / (num_chroms - 1.0)) * 2 * p * (1 - p)

    property num_hom_ref:
        def __get__(self):
            if self._gt_types == NULL:
                self.gt_types
            cdef int n = 0, i = 0
            for i in range(self.vcf.n_samples):
                if self._gt_types[i] == 0:
                    n+=1
            return n

    property num_het:
        def __get__(self):
            if self._gt_types == NULL:
                self.gt_types
            cdef int n = 0, i = 0
            for i in range(self.vcf.n_samples):
                if self._gt_types[i] == 1:
                    n+=1
            return n

    property num_hom_alt:
        def __get__(self):
            if self._gt_types == NULL:
                self.gt_types
            cdef int n = 0, i = 0
            if self.vcf.gts012:
                for i in range(self.vcf.n_samples):
                    if self._gt_types[i] == 2:
                        n+=1
            else:
                for i in range(self.vcf.n_samples):
                    if self._gt_types[i] == 3:
                        n+=1
            return n

    property num_unknown:
        def __get__(self):
            if self._gt_types == NULL:
                self.gt_types
            cdef int n = 0, i = 0
            for i in range(self.vcf.n_samples):
                if self._gt_types[i] == 2:
                    n+=1
            return n

    def format(self, itag, vtype=None):
        """
        type is one of [int, float, str]
        returns None if the key isn't found.
        """
        if vtype is None:
            vtype = self.vcf.get_type(itag)

        cdef bytes tag = to_bytes(itag)
        cdef bcf_fmt_t *fmt = bcf_get_fmt(self.vcf.hdr, self.b, tag)
        cdef int n = 0, nret
        cdef void *buf = NULL;
        cdef int typenum = 0
        if vtype == "Integer" or vtype == int:
            nret = bcf_get_format_int32(self.vcf.hdr, self.b, tag, <int **>&buf, &n)
            typenum = np.NPY_INT32
        elif vtype == "Float" or vtype == float:
            nret = bcf_get_format_float(self.vcf.hdr, self.b, tag, <float **>&buf, &n)
            typenum = np.NPY_FLOAT32
        elif vtype == "String" or vtype == str or vtype == "Character":
            vtype = str
            nret = bcf_get_format_string(self.vcf.hdr, self.b, tag, <char ***>&buf, &n)
            typenum = np.NPY_STRING
        else:
            raise Exception("type %s not supported to format()" % vtype)
        if nret < 0:
            return None

        cdef char **dst
        cdef int i
        cdef np.npy_intp shape[2]
        shape[0] = <np.npy_intp> self.vcf.n_samples
        shape[1] = fmt.n # values per sample

        if vtype == str:
            dst = <char **>buf
            v = [dst[i] for i in range(self.vcf.n_samples)]
            xret = np.array(v, dtype=str)
            stdlib.free(dst[0])
            stdlib.free(dst)
            return xret

        iv = np.PyArray_SimpleNewFromData(2, shape, typenum, buf)
        iret = np.array(iv)
        stdlib.free(buf)
        return iret

    property gt_types:
        def __get__(self):
            cdef int ndst, ngts, n, i, nper, j = 0, k = 0
            cdef int a
            if self.vcf.n_samples == 0:
                return []
            if self._gt_types == NULL:
                self._gt_phased = <int *>stdlib.malloc(sizeof(int) * self.vcf.n_samples)
                ndst = 0
                ngts = bcf_get_genotypes(self.vcf.hdr, self.b, &self._gt_types, &ndst)
                nper = ndst / self.vcf.n_samples
                self._ploidy = nper
                self._gt_idxs = <int *>stdlib.malloc(sizeof(int) * self.vcf.n_samples * nper)
                if ndst == 0 or nper == 0:
                    return []
                for i in range(0, ndst, nper):
                    for k in range(i, i + nper):
                        a = self._gt_types[k]
                        if a >= 0:
                            self._gt_idxs[k] = bcf_gt_allele(a)
                        else:
                            self._gt_idxs[k] = a

                    self._gt_phased[j] = self._gt_types[i] > 0 and <int>bcf_gt_is_phased(self._gt_types[i+1])
                    j += 1

                if self.vcf.gts012:
                    n = as_gts012(self._gt_types, self.vcf.n_samples, nper)
                else:
                    n = as_gts(self._gt_types, self.vcf.n_samples, nper)
            cdef np.npy_intp shape[1]
            shape[0] = <np.npy_intp> self.vcf.n_samples
            return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT32, self._gt_types)

    property ploidy:
        def __get__(self):
            if self._ploidy == -1:
                self.gt_types
            return self._ploidy

    property gt_phred_ll_homref:
        def __get__(self):
            if self.vcf.n_samples == 0:
                return []
            cdef int ndst = 0, nret=0, n, i, j, nper

            cdef int imax = np.iinfo(np.int32(0)).max

            if self._gt_pls == NULL and self._gt_gls == NULL:
                nret = bcf_get_format_int32(self.vcf.hdr, self.b, "PL", &self._gt_pls, &ndst)
                if nret < 0:
                    nret = bcf_get_format_float(self.vcf.hdr, self.b, "GL", &self._gt_gls, &ndst)
                    if nret < 0:
                        return []
                    else:
                        for i in range(nret):
                            if self._gt_gls[i] <= -2147483646:
                                # this gets translated on conversion to PL
                                self._gt_gls[i] = imax / -10.0
                else:
                    for i in range(nret):
                        if self._gt_pls[i] < 0:
                            self._gt_pls[i] = imax

                self._gt_nper = nret / self.vcf.n_samples
            cdef np.npy_intp shape[1]
            shape[0] = <np.npy_intp> self._gt_nper * self.vcf.n_samples
            if self._gt_pls != NULL:
                pls = np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT32,
                        self._gt_pls)[::self._gt_nper]
                return pls
            else:
                gls = np.PyArray_SimpleNewFromData(1, shape, np.NPY_FLOAT32,
                        self._gt_gls)[::self._gt_nper]
                gls = (-10 * gls).round().astype(np.int32)
                return gls

    property gt_phred_ll_het:
        def __get__(self):
            if self.vcf.n_samples == 0:
                return []
            if self._gt_pls == NULL and self._gt_gls == NULL:
                # NOTE: the missing values for all homref, het, homalt are set
                # by this call.
                self.gt_phred_ll_homref
            cdef np.npy_intp shape[1]
            shape[0] = <np.npy_intp> self._gt_nper * self.vcf.n_samples
            if self._gt_pls != NULL:
                if self._gt_nper > 1:
                    ret = np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT32, self._gt_pls)[1::self._gt_nper]
                    return ret

                return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT32, self._gt_pls)
            else:
                if self._gt_nper > 1:
                    gls = np.PyArray_SimpleNewFromData(1, shape, np.NPY_FLOAT32,
                            self._gt_gls)[1::self._gt_nper]
                else:
                    gls = np.PyArray_SimpleNewFromData(1, shape, np.NPY_FLOAT32, self._gt_gls)
                gls = (-10 * gls).round().astype(np.int32)
                return gls

    property gt_phred_ll_homalt:
        def __get__(self):
            if self.vcf.n_samples == 0:
                return []
            if self._gt_pls == NULL and self._gt_gls == NULL:
                self.gt_phred_ll_homref
            cdef np.npy_intp shape[1]
            shape[0] = <np.npy_intp> self._gt_nper * self.vcf.n_samples
            if self._gt_pls != NULL:
                if self._gt_nper > 1:
                    return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT32,
                            self._gt_pls)[2::self._gt_nper]
                return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT32, self._gt_pls)
            else:
                if self._gt_nper > 1:
                    gls = np.PyArray_SimpleNewFromData(1, shape, np.NPY_FLOAT32,
                            self._gt_gls)[2::self._gt_nper]
                else:
                    gls = np.PyArray_SimpleNewFromData(1, shape, np.NPY_FLOAT32,
                            self._gt_gls)
                gls = (-10 * gls).round().astype(np.int32)
                return gls

    property gt_ref_depths:
        def __get__(self):
            cdef int ndst, nret = 0, n, i, j = 0, nper = 0
            if self.vcf.n_samples == 0:
                return []
            if self._gt_ref_depths == NULL:
                ndst = 0
                # GATK
                nret = bcf_get_format_int32(self.vcf.hdr, self.b, "AD", &self._gt_ref_depths, &ndst)
                if nret > 0:
                    nper = nret / self.vcf.n_samples
                    if nper == 1:
                        stdlib.free(self._gt_ref_depths); self._gt_ref_depths = NULL
                        return -1 + np.zeros(self.vcf.n_samples, np.int32)

                    for i in range(0, nret, nper):
                        self._gt_ref_depths[j] = self._gt_ref_depths[i]
                        j += 1
                elif nret == -1:
                    # Freebayes
                    # RO has to be 1:1
                    nret = bcf_get_format_int32(self.vcf.hdr, self.b, "RO", &self._gt_ref_depths, &ndst)
                    if nret < 0:
                        stdlib.free(self._gt_ref_depths); self._gt_ref_depths = NULL
                        return -1 + np.zeros(self.vcf.n_samples, np.int32)
                # TODO: add new vcf standard.
                else:
                    stdlib.free(self._gt_ref_depths); self._gt_ref_depths = NULL
                    return -1 + np.zeros(self.vcf.n_samples, np.int32)

                for i in range(self.vcf.n_samples):
                    if self._gt_ref_depths[i] < 0:
                        self._gt_ref_depths[i] = -1
            else:
                pass

            cdef np.npy_intp shape[1]
            shape[0] = <np.npy_intp> self.vcf.n_samples
            return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT32, self._gt_ref_depths)


    property gt_alt_depths:
        def __get__(self):
            cdef int ndst, nret = 0, n, i, j = 0, k = 0, nper = 0
            if self.vcf.n_samples == 0:
                return []
            if self._gt_alt_depths == NULL:
                ndst = 0
                # GATK
                nret = bcf_get_format_int32(self.vcf.hdr, self.b, "AD", &self._gt_alt_depths, &ndst)
                if nret > 0:
                    nper = nret / self.vcf.n_samples
                    if nper == 1:
                        stdlib.free(self._gt_alt_depths); self._gt_alt_depths = NULL
                        return (-1 + np.zeros(self.vcf.n_samples, np.int32))

                    for i in range(0, nret, nper):
                        self._gt_alt_depths[j] = self._gt_alt_depths[i+1]
                        # add up all the alt alleles
                        for k in range(2, nper):
                            self._gt_alt_depths[j] += self._gt_alt_depths[i+k]
                        j += 1

                elif nret == -1:
                    # Freebayes
                    nret = bcf_get_format_int32(self.vcf.hdr, self.b, "AO", &self._gt_alt_depths, &ndst)
                    nper = nret / self.vcf.n_samples
                    if nret < 0:
                        stdlib.free(self._gt_alt_depths); self._gt_alt_depths = NULL
                        return -1 + np.zeros(self.vcf.n_samples, np.int32)
                    for i in range(0, nret, nper):
                        self._gt_alt_depths[j] = self._gt_alt_depths[i]
                        for k in range(1, nper):
                            self._gt_alt_depths[j] += self._gt_alt_depths[i+k]
                        j += 1
                else:
                    stdlib.free(self._gt_alt_depths); self._gt_alt_depths = NULL
                    return -1 + np.zeros(self.vcf.n_samples, np.int32)

                # TODO: add new vcf standard.
            for i in range(self.vcf.n_samples):
                if self._gt_alt_depths[i] < 0:
                    self._gt_alt_depths[i] = -1

            cdef np.npy_intp shape[1]
            shape[0] = <np.npy_intp> self.vcf.n_samples
            return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT32, self._gt_alt_depths)

    property gt_quals:
        def __get__(self):
            if self.vcf.n_samples == 0:
                return []
            cdef int ndst = 0, nret, n, i
            cdef int *gq
            cdef np.ndarray[np.float32_t, ndim=1] a
            if self._gt_quals == NULL and self._int_gt_quals == NULL:
                nret = bcf_get_format_int32(self.vcf.hdr, self.b, "GQ", &self._int_gt_quals, &ndst)
                if nret == -2: # defined as int
                    ndst = 0
                    nret = bcf_get_format_float(self.vcf.hdr, self.b, "GQ", &self._gt_quals, &ndst)
                if nret < 0 and nret != -2:
                    return -1.0 + np.zeros(self.vcf.n_samples, np.float32)
            cdef np.npy_intp shape[1]
            shape[0] = <np.npy_intp> self.vcf.n_samples
            if self._int_gt_quals != NULL:
                a = np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT32, self._int_gt_quals).astype(np.float32)
                a[a < 0] = -1
            else:
                a = np.PyArray_SimpleNewFromData(1, shape, np.NPY_FLOAT32, self._gt_quals)
                # this take up 10% of the total vcf parsing time. fix!!
                a[np.isnan(a)] = -1
            return a

    property gt_depths:
        def __get__(self):
            if self.vcf.n_samples == 0:
                return []
            # unfortunately need to create a new array here since we're modifying.
            r = np.array(self.gt_ref_depths, np.int32)
            a = np.array(self.gt_alt_depths, np.int32)
            # keep the -1 for empty.
            rl0 = r < 0
            al0 = a < 0
            r[rl0] = 0
            a[al0] = 0
            depth = r + a
            depth[rl0 & al0] = -1
            return depth

    property gt_phases:
        def __get__(self):
            # run for side-effect
            if self._gt_phased == NULL:
                self.gt_types
            cdef np.npy_intp shape[1]
            shape[0] = <np.npy_intp> self.vcf.n_samples

            return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT32, self._gt_phased).astype(bool)


    property REF:
        def __get__(self):
            return self.b.d.allele[0].decode()

    property ALT:
        def __get__(self):
            cdef int i
            return [self.b.d.allele[i].decode() for i in range(1, self.b.n_allele)]

    property is_snp:
        def __get__(self):
            cdef int i
            if len(self.b.d.allele[0]) > 1: return False
            for i in range(1, self.b.n_allele):
                if not self.b.d.allele[i] in (b"A", b"C", b"G", b"T"):
                    return False
            return True

    property is_indel:
        def __get__(self):
            cdef int i
            is_sv = self.is_sv
            if len(self.b.d.allele[0]) > 1 and not is_sv: return True

            if len(self.REF) > 1 and not is_sv: return True

            for i in range(1, self.b.n_allele):
                alt = self.b.d.allele[i]
                if alt == b".":
                    return True
                if len(alt) != len(self.REF):
                    if not is_sv:
                        return True
            return False

    property is_transition:
        def __get__(self):
            if len(self.ALT) > 1: return False

            if not self.is_snp: return False
            ref = self.REF
                # just one alt allele
            alt_allele = self.ALT[0]
            if ((ref == b'A' and alt_allele == b'G') or
                (ref == b'G' and alt_allele == b'A') or
                (ref == b'C' and alt_allele == b'T') or
                (ref == b'T' and alt_allele == b'C')):
                    return True
            return False

    property is_deletion:
        def __get__(self):
            if len(self.ALT) > 1: return False

            if not self.is_indel: return False
            alt = self.ALT[0]
            if alt is None or alt == ".":
                return True

            if len(self.REF) > len(alt):
                return True
            return False

    property is_sv:
        def __get__(self):
            return self.INFO.get(b'SVTYPE') is not None

    property CHROM:
        def __get__(self):
            return bcf_hdr_id2name(self.vcf.hdr, self.b.rid).decode()

    property var_type:
        def __get__(self):
           if self.is_snp:
               return "snp"
           elif self.is_indel:
               return "indel"
           elif self.is_sv:
               return "sv"
           else:
               return "unknown"

    property var_subtype:
        def __get__(self):
            if self.is_snp:
                if self.is_transition:
                    return "ts"
                if len(self.ALT) == 1:
                    return "tv"
                return "unknown"

            elif self.is_indel:
                if self.is_deletion:
                    return "del"
                if len(self.ALT) == 1:
                    return "ins"
                else:
                    return "unknown"

            svt = self.INFO.get("SVTYPE")
            if svt is None:
                return "unknown"
            if svt == "BND":
                return "complex"
            if self.INFO.get('IMPRECISE') is None:
                return svt
            return self.ALT[0].strip('<>')

    property start:
        def __get__(self):
            return self.b.pos

    property end:
        def __get__(self):
            return self.b.pos + self.b.rlen

    property ID:
        def __get__(self):
            cdef char *id = self.b.d.id
            if id == b".": return None
            return id

    property FILTER:
        def __get__(self):
            cdef int i
            cdef int n = self.b.d.n_flt
            if n == 1:
                if self.vcf.PASS != -1:
                    if self.b.d.flt[0] == self.vcf.PASS:
                        return None
                else:
                    v = bcf_hdr_int2id(self.vcf.hdr, BCF_DT_ID, self.b.d.flt[0])
                    if v == b"PASS":
                        self.vcf.PASS = self.b.d.flt[0]
                        return None
                    return v
            if n == 0:
                return None
            return b';'.join(bcf_hdr_int2id(self.vcf.hdr, BCF_DT_ID, self.b.d.flt[i]) for i in range(n))

        def __set__(self, filters):
            if isinstance(filters, basestring):
                filters = filters.split(";")
            cdef bcf_hdr_t *h = self.vcf.hdr
            cdef int *flt_ids = <int *>stdlib.malloc(sizeof(int) * len(filters))
            for i, fname in enumerate(filters):
                flt_ids[i] = bcf_hdr_id2int(h, BCF_DT_ID, to_bytes(fname))
            ret = bcf_update_filter(h, self.b, flt_ids, len(filters))
            stdlib.free(flt_ids)
            if ret != 0:
                raise Exception("not able to set filter: %s", filters)

    property QUAL:
        def __get__(self):
            cdef float q = self.b.qual
            if bcf_float_is_missing(q):
                return None
            return q

cdef inline HREC newHREC(bcf_hrec_t *hrec, bcf_hdr_t *hdr):
    cdef HREC h = HREC.__new__(HREC)
    h.hdr = hdr
    h.hrec = hrec
    return h

cdef class HREC(object):
    cdef bcf_hdr_t *hdr
    cdef bcf_hrec_t *hrec

    def __cinit__(HREC self):
        pass

    def __dealloc__(self):
        #bcf_hrec_destroy(self.hrec)
        self.hrec = NULL
        self.hdr = NULL

    @property
    def type(self):
        return ["FILTER", "INFO", "FORMAT", "CONTIG", "STR", "GENERIC"][self.hrec.type]

    def __getitem__(self, key):
        for i in range(self.hrec.nkeys):
            if self.hrec.keys[i] == key:
                return self.hrec.vals[i]
        raise KeyError

    def info(self, extra=False):
        """
        return a dict with commonly used stuffs
        """
        d = {}
        for k in ('Type', 'Number', 'ID', 'Description'):
            try:
                d[k] = self[k]
            except KeyError:
                continue
        d['HeaderType'] = self.type
        if extra:
            for i in range(self.hrec.nkeys):
                k = self.hrec.keys[i]
                if k in d: continue
                d[k] = self.hrec.vals[i]
        return d

    def __repr__(self):
        return str(self.info())

cdef class INFO(object):
    cdef bcf_hdr_t *hdr
    cdef bcf1_t *b
    cdef int _i

    def __cinit__(INFO self):
        self._i = 0

    def __setitem__(self, key, value):
        # only support strings for now.
        if value is True or value is False:

            ret = bcf_update_info_flag(self.hdr, self.b, to_bytes(key), b"", int(value))
            if ret != 0:
                raise Exception("not able to set flag", key, value, ret)
            return

        ret = bcf_update_info_string(self.hdr, self.b, to_bytes(key), to_bytes(value))
        if ret != 0:
            raise Exception("not able to set: %s -> %s (%d)", key, value, ret)

    cdef _getval(INFO self, bcf_info_t * info, char *key):

        if info.len == 1:
            if info.type == BCF_BT_INT8:
                if info.v1.i == INT8_MIN:
                    return None
                return <int>(info.v1.i)

            if info.type == BCF_BT_INT16:
                if info.v1.i == INT16_MIN:
                    return None
                return <int>(info.v1.i)

            if info.type == BCF_BT_INT32:
                if info.v1.i == INT32_MIN:
                    return None
                return <int>(info.v1.i)

            if info.type == BCF_BT_FLOAT:
                if bcf_float_is_missing(info.v1.f):
                    return None
                return info.v1.f

        if info.type == BCF_BT_CHAR:
            v = info.vptr[:info.vptr_len]
            if len(v) > 0 and v[0] == 0x7:
                return None
            return from_bytes(v)

        return bcf_array_to_object(info.vptr, info.type, info.len)

    def __getitem__(self, okey):
        okey = to_bytes(okey)
        cdef char *key = okey
        cdef bcf_info_t *info = bcf_get_info(self.hdr, self.b, key)
        if info == NULL:
            raise KeyError(key)
        return self._getval(info, key)

    def get(self, key, default=None):
        try:
            return self.__getitem__(key)
        except KeyError:
            return default

    def __iter__(self):
        self._i = 0
        return self

    def __next__(self):
        cdef bcf_info_t *info = NULL
        cdef char *name
        while info == NULL:
            if self._i >= self.b.n_info:
                raise StopIteration
            info = &(self.b.d.info[self._i])
            self._i += 1
        name = bcf_hdr_int2id(self.hdr, BCF_DT_ID, info.key)
        return name.decode(), self._getval(info, name)


# this function is copied verbatim from pysam/cbcf.pyx
cdef bcf_array_to_object(void *data, int type, int n, int scalar=0):
    cdef char    *datac
    cdef int8_t  *data8
    cdef int16_t *data16
    cdef int32_t *data32
    cdef float   *dataf
    cdef int      i

    if not data or n <= 0:
        return None

    if type == BCF_BT_CHAR:
        datac = <char *>data
        value = datac[:n].decode() if datac[0] != bcf_str_missing else None
    else:
        value = []
        if type == BCF_BT_INT8:
            data8 = <int8_t *>data
            for i in range(n):
                if data8[i] == bcf_int8_vector_end:
                    break
                value.append(data8[i] if data8[i] != bcf_int8_missing else None)
        elif type == BCF_BT_INT16:
            data16 = <int16_t *>data
            for i in range(n):
                if data16[i] == bcf_int16_vector_end:
                    break
                value.append(data16[i] if data16[i] != bcf_int16_missing else None)
        elif type == BCF_BT_INT32:
            data32 = <int32_t *>data
            for i in range(n):
                if data32[i] == bcf_int32_vector_end:
                    break
                value.append(data32[i] if data32[i] != bcf_int32_missing else None)
        elif type == BCF_BT_FLOAT:
            dataf = <float *>data
            for i in range(n):
                if bcf_float_is_vector_end(dataf[i]):
                    break
                value.append(dataf[i] if not bcf_float_is_missing(dataf[i]) else None)
        else:
            raise TypeError('unsupported info type code')

        if not value:
            value = None
        elif scalar and len(value) == 1:
            value = value[0]
        else:
            value = tuple(value)

    return value

cdef inline Variant newVariant(bcf1_t *b, VCF vcf):
    cdef Variant v = Variant.__new__(Variant)
    v.b = b
    if not vcf.lazy:
        with nogil:
            bcf_unpack(v.b, 15)
    else:
        with nogil:
            bcf_unpack(v.b, 1|2|4)

    v.vcf = vcf
    v.POS = v.b.pos + 1
    cdef INFO i = INFO.__new__(INFO)
    i.b, i.hdr = b, vcf.hdr
    v.INFO = i
    return v

cdef to_bytes(s, enc=ENC):
    if not isinstance(s, bytes):
        return s.encode(enc)
    return s

cdef from_bytes(s):
    if isinstance(s, bytes):
        return s.decode(ENC)
    return s


cdef class Writer(object):
    cdef htsFile *hts
    cdef bcf_hdr_t *hdr
    cdef public bytes name
    cdef bint header_written

    def __init__(self, fname, VCF tmpl):
        self.name = to_bytes(fname)
        self.hts = hts_open(self.name, "w")

        cdef bcf_hdr_t *h = tmpl.hdr
        cdef bcf_hdr_t *hdup = bcf_hdr_dup(h)
        self.hdr = hdup
        self.header_written = False
        samples = to_bytes(",".join(tmpl.samples))
        bcf_hdr_set_samples(self.hdr, samples, 0)

    def write_record(self, Variant var):
        if not self.header_written:
            bcf_hdr_write(self.hts, self.hdr)
            self.header_written = True
        return bcf_write(self.hts, self.hdr, var.b)

    def close(self):
        if self.hts != NULL:
            hts_close(self.hts)
            self.hts = NULL

    def __dealloc__(self):
        bcf_hdr_destroy(self.hdr)
        self.hdr = NULL
        self.close()

