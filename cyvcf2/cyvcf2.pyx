#cython: profile=False
#cython: embedsignature=True
from __future__ import print_function
import os
import sys
from collections import defaultdict
import atexit
import tempfile
import numpy as np
from array import array
import math
import ctypes


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


def r_(int32_t[::view.contiguous] a_gts, int32_t[::view.contiguous] b_gts, float f, int32_t n_samples):
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
    """
    VCF class holds methods to iterate over and query a VCF.

    Parameters
    ----------
    fname: str
        path to file
    gts012: bool
        if True, then gt_types will be 0=HOM_REF, 1=HET, 2=HOM_ALT, 3=UNKNOWN. If False, 3, 2 are flipped.
    lazy: bool
        if True, then don't unpack (parse) the underlying record until needed.
    strict_gt: bool
        if True, then any '.' present in a genotype will classify the corresponding element in the gt_types array as UNKNOWN.
    samples: list
        list of samples to extract from full set in file.


    Returns
    -------
    VCF object for iterating and querying.
    """

    cdef htsFile *hts
    cdef const bcf_hdr_t *hdr
    cdef tbx_t *idx
    cdef hts_idx_t *hidx
    cdef int n_samples
    cdef int PASS
    cdef bytes fname
    cdef bint gts012
    cdef bint lazy
    cdef bint strict_gt
    cdef list _seqnames
    cdef list _seqlens
    # holds a lookup of format field -> type.
    cdef dict format_types

    #: The constant used to indicate the the genotype is HOM_REF.
    cdef readonly int HOM_REF
    #: The constant used to indicate the the genotype is HET.
    cdef readonly int HET
    #: The constant used to indicate the the genotype is HOM_ALT.
    cdef readonly int HOM_ALT
    #: The constant used to indicate the the genotype is UNKNOWN.
    cdef readonly int UNKNOWN

    def __init__(self, fname, mode="r", gts012=False, lazy=False, strict_gt=False, samples=None, threads=None):
        cdef hFILE *hf

        if isinstance(fname, basestring):
            if fname == b"-" or fname == "-":
                fname = b"/dev/stdin"
            fname, mode = to_bytes(fname), to_bytes(mode)
            self.hts = hts_open(fname, mode)
            self.fname = fname
        else:
            mode = to_bytes(mode)
            hf = hdopen(int(fname), mode)
            self.hts = hts_hopen(hf, "<file>", mode)

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
        self.gts012 = gts012
        self.lazy = lazy
        self.strict_gt = strict_gt
        self._seqnames = []
        self._seqlens = []
        set_constants(self)
        self.format_types = {}
        if threads is not None:
            self.set_threads(threads)

    def set_threads(self, int n):
        v = hts_set_threads(self.hts, n)
        if v < 0:
            raise Exception("error setting number of threads: %d" % v)

    cdef get_type(self, fmt):
        fmt = from_bytes(fmt)
        if not fmt in self.format_types:
            s = self.get_header_type(fmt, order=[BCF_HL_FMT])
            self.format_types[fmt] = s["Type"]
        return from_bytes(self.format_types[fmt])

    def add_to_header(self, line):
        """Add a new line to the VCF header.

        Parameters
        ----------
        line: str
            full vcf header line.
        """
        ret = bcf_hdr_append(self.hdr, to_bytes(line))
        if ret != 0:
            raise Exception("couldn't add '%s' to header")
        ret = bcf_hdr_sync(self.hdr)
        if ret != 0:
            raise Exception("couldn't add '%s' to header")
        return ret

    def add_info_to_header(self, adict):
        """Add a INFO line to the VCF header.

        Parameters
        ----------
        adict: dict
            dict containing keys for ID, Number, Type, Description.
        """
        return self.add_to_header("##INFO=<ID={ID},Number={Number},Type={Type},Description=\"{Description}\">".format(**adict))

    def add_format_to_header(self, adict):
        """Add a FORMAT line to the VCF header.

        Parameters
        ----------
        adict: dict
            dict containing keys for ID, Number, Type, Description.
        """
        return self.add_to_header("##FORMAT=<ID={ID},Number={Number},Type={Type},Description=\"{Description}\">".format(**adict))

    def add_filter_to_header(self, adict):
        """Add a FILTER line to the VCF header.

        Parameters
        ----------
        adict: dict
            dict containing keys for ID, Description.
        """
        return self.add_to_header("##FILTER=<ID={ID},Description=\"{Description}\">".format(**adict))

    def set_samples(self, samples):
        """Set the samples to be pulled from the VCF; this must be called before any iteration.

        Parameters
        ----------
        samples: list
            list of samples to extract.
        """
        if samples is None:
            samples = "-".encode()
        if isinstance(samples, list):
            samples = to_bytes(",".join(samples))
        else:
            samples = to_bytes(samples)

        ret = bcf_hdr_set_samples(self.hdr, <const char *>samples, 0)
        assert ret >= 0, ("error setting samples", ret)
        if ret != 0 and samples != "-":
            s = from_bytes(samples).split(",")
            if ret < len(s):
                sys.stderr.write("warning: not all requested samples found in VCF\n")

        self.n_samples = bcf_hdr_nsamples(self.hdr)

    def update(self, id, type, number, description):
        """Update the header with an INFO field of the given parameters.

        Parameters
        ----------
        id: str
            ID
        type: str
            valid VCF type
        number: str
             valid VCF number
        description: str
             description of added line.
        """
        ret = bcf_hdr_append(self.hdr, "##INFO=<ID={id},Number={number},Type={type},Description=\"{description}\">".format(id=id, type=type, number=number, description=description))
        if ret != 0:
            raise Exception("unable to update to header: %d", ret)
        ret = bcf_hdr_sync(self.hdr)
        if ret != 0:
            raise Exception("unable to update to header")

    def set_index(self, index_path=""):
        self.hidx = hts_idx_load2(to_bytes(self.fname), to_bytes(index_path))
        if self.hidx == NULL:
            self.idx = tbx_index_load2(to_bytes(self.fname), to_bytes(index_path))
        if self.hidx == NULL and self.idx == NULL:
          raise OSError("unable to open index:'%s' for '%s'" % (index_path, self.fname))

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
                if bcf_subset_format(self.hdr, b) != 0:
                    sys.stderr.write("could not subset variant")
                    bcf_destroy(b)
                    break
                yield newVariant(b, self)
        finally:
            if itr != NULL:
                hts_itr_destroy(itr)

    def __call__(VCF self, region=None):
        """
        Extract the region from the VCF.

        Parameters
        ----------
        region: str
           region string like chr1:1234-34566 or 'chr7

        Returns
        -------
        An Iterator over the requested region.
        """
        if not region:
            yield from self
            raise StopIteration

        if self.fname.decode(ENC).endswith('.bcf'):
            yield from self._bcf_region(region)
            raise StopIteration

        if self.idx == NULL:
            self.idx = tbx_index_load(to_bytes(self.fname))
            assert self.idx != NULL, "error loading tabix index for %s" % self.fname

        cdef hts_itr_t *itr
        cdef kstring_t s
        cdef bcf1_t *b
        cdef int slen = 1, ret = 0
        cdef bytes bregion = to_bytes(region)
        cdef char *cregion = bregion

        with nogil:
            itr = tbx_itr_querys(self.idx, cregion)

        if itr == NULL:
            sys.stderr.write("no intervals found for %s at %s\n" % (self.fname, region))
            raise StopIteration

        try:
            while 1:
                with nogil:
                    slen = tbx_itr_next(self.hts, self.idx, itr, &s)
                    if slen > 0:
                            b = bcf_init()
                            ret = vcf_parse(&s, self.hdr, b)
                if slen <= 0:
                    break
                if ret > 0:
                    bcf_destroy(b)
                    stdlib.free(s.s)
                    hts_itr_destroy(itr)
                    raise Exception("error parsing")
                yield newVariant(b, self)
        finally:
            stdlib.free(s.s)
            hts_itr_destroy(itr)

    def variant_from_string(self, variant_string):
        cdef bcf1_t *b = bcf_init()
        cdef kstring_t s
        tmp = to_bytes(variant_string)
        s.s = tmp
        s.l = len(variant_string)
        s.m = len(variant_string)
        ret = vcf_parse(&s, self.hdr, b)
        if ret > 0:
            bcf_destroy(b)
            raise Exception("error parsing:" + variant_string + " return value:" + ret)

        v = newVariant(b, self)
        return v

    def header_iter(self):
        """
        Iterate over fields in the HEADER
        """
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
    cpdef get_header_type(self, key, order=[BCF_HL_INFO, BCF_HL_FMT]):
        """Extract a field from the VCF header by id.

        Parameters
        ----------
        key: str
           ID to pull from the header.
        Returns
        -------
        rec: dict
            dictionary containing header information.
        """
        key = to_bytes(key)
        cdef bcf_hrec_t *b
        cdef int i
        for typ in order:
            b = bcf_hdr_get_hrec(self.hdr, typ, b"ID", key, NULL);
            if b != NULL:
                break
        if b == NULL:
            b = bcf_hdr_get_hrec(self.hdr, BCF_HL_GEN, key, NULL, NULL);
            if b == NULL:
                raise KeyError(key)
            d = {from_bytes(b.key): from_bytes(b.value)}
        else:
            d =  {from_bytes(b.keys[i]): from_bytes(b.vals[i]) for i in range(b.nkeys)}
        #bcf_hrec_destroy(b)
        return d

    def __getitem__(self, key):
        return self.get_header_type(key)

    def __contains__(self, key):
        """Check if the given ID is in the header."""
        try:
            self[key]
            return True
        except KeyError:
            return False

    contains = __contains__

    def close(self):
        if self.hts != NULL:
            if self.fname != "-":
                # TODO flush
                hts_close(self.hts)
            self.hts = NULL

    def __dealloc__(self):
        if self.hts != NULL and self.hdr != NULL:
            bcf_hdr_destroy(self.hdr)
            self.hdr = NULL
        self.close()
        if self.idx != NULL:
            tbx_destroy(self.idx)
        if self.hidx != NULL:
            hts_idx_destroy(self.hidx)

    def __iter__(self):
        return self

    def __next__(self):

        cdef bcf1_t *b
        cdef int ret
        if self.hts == NULL:
            raise Exception("attempt to iterate over closed/invalid VCF")
        with nogil:
            b = bcf_init()
            ret = bcf_read(self.hts, self.hdr, b)
        if ret >= 0:
            return newVariant(b, self)
        else:
            bcf_destroy(b)
        raise StopIteration

    property samples:
        "list of samples pulled from the VCF."
        def __get__(self):
            cdef int i
            return [str(self.hdr.samples[i].decode('utf-8')) for i in range(self.n_samples)]

    property raw_header:
        "string of the raw header from the VCF"
        def __get__(self):
            cdef int hlen
            s = bcf_hdr_fmt_text(self.hdr, 0, &hlen)
            return from_bytes(s)

    property seqlens:
        def __get__(self):
            if len(self._seqlens) > 0: return self._seqlens
            cdef int32_t nseq;
            cdef int32_t* sls = bcf_hdr_seqlen(self.hdr, &nseq)
            self._seqlens = [sls[i] for i in range(nseq)]
            stdlib.free(sls)
            return self._seqlens

    property seqnames:
        "list of chromosomes in the VCF"
        def __get__(self):
            if len(self._seqnames) > 0: return self._seqnames
            cdef char **cnames
            cdef int i, n = 0
            cnames = bcf_hdr_seqnames(self.hdr, &n)
            if n == 0 and self.fname.decode(ENC).endswith('.bcf'):
                if self.hidx == NULL:
                    self.hidx = bcf_index_load(self.fname)
                if self.hidx != NULL:
                    cnames = bcf_index_seqnames(self.hidx, self.hdr, &n)
            elif n == 0:
                if self.idx == NULL:
                    self.idx = tbx_index_load(to_bytes(self.fname))
                if self.idx !=NULL:
                    cnames = tbx_seqnames(self.idx, &n)

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
                    if not i[0] in seqnames and 'chr' + i[0] in seqnames:
                        i[0] = 'chr' + i[0]
                    i[1] = int(i[1])
                    isites.append(i)
            else:
                isites = sites
                for i in isites:
                    if not i[0] in seqnames and 'chr' + i[0] in seqnames:
                        i[0] = 'chr' + i[0]

        cdef Variant v
        cdef int k, last_pos
        if sites:
            isites = isites[offset::each]
            ref, alt = None, None
            j = 0
            for i, osite in enumerate(isites):
                if len(osite) >= 4:
                    chrom, pos, ref, alt = osite[:4]
                else:
                    chrom, pos = osite[:2]
                for v in self("%s:%s-%s" % (chrom, pos, pos)):
                    if len(v.ALT) != 1: continue
                    if ref is not None and v.REF != ref: continue
                    if alt is not None and v.ALT[0] != alt: continue
                    if v.call_rate < call_rate: continue
                    yield i, v
                    j += 1
                    break
        else:
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

        maf_lists = defaultdict(list)
        idxs = np.arange(n_samples)
        for i, v in self.gen_variants(sites, each=each, offset=offset):
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


        cdef int32_t[:, ::view.contiguous] ibs = np.zeros((n_samples, n_samples), np.int32)
        cdef int32_t[:, ::view.contiguous] n = np.zeros((n_samples, n_samples), np.int32)
        cdef int32_t[:] hets = np.zeros((n_samples, ), np.int32)
        cdef int32_t[:] gt_types = np.zeros((n_samples, ), np.int32)
        cdef int32_t[:] depths = np.zeros((n_samples, ), np.int32)
        cdef double[:] alt_freqs = np.zeros((n_samples,), np.double)

        cdef Variant v

        for j, (i, v) in enumerate(self.gen_variants(sites, offset=offset, each=each)):
            gt_types = v.gt_types
            alt_freqs = v.gt_alt_freqs
            krelated(&gt_types[0], &ibs[0, 0], &n[0, 0], &hets[0], n_samples,
                    &alt_freqs[0])

        return ibs, n, hets

    def relatedness(self, int n_variants=35000, int gap=30000, float min_af=0.04,
                    float max_af=0.8, float linkage_max=0.2, min_depth=8):
        cdef Variant v

        cdef int last = -gap, nv = 0, nvt=0
        cdef int32_t *last_gts = NULL
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
            if _hets[sj] == 0:
                print("peddy: no hets found for sample %s\n" % sample_j, file=sys.stderr)
            for sk in range(sj, ns):
                if sj == sk: continue
                sample_k = samples[sk]

                # calculate relatedness. we use the geometric mean.
                #bot = math.exp(0.5 * (math.log(1 + _hets[sk]) + math.log(1 + _hets[sj])))
                #bot = (_hets[sk] + _hets[sj])/2.0
                bot = min(_hets[sk], _hets[sj])
                if bot == 0:
                    bot = max(_hets[sk], _hets[sj])
                    if bot == 0:
                        # set to negative value if we are unable to calculate it.
                        bot = -1

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

cdef class Allele(object):
    cdef int32_t *_raw
    cdef int i

    cdef int _value(self):
        if self._raw[self.i] < 0: return self._raw[self.i]
        return (self._raw[self.i] >> 1) - 1

    @property
    def phased(self):
        return self._raw[self.i] & 1 == 1

    @phased.setter
    def phased(self, bint ph):
        if ph:
            self._raw[self.i] = (self._value() + 1)<<1|1
        else:
            self._raw[self.i] = (self._value() + 1)<<1

    @property
    def value(self):
        if self._raw[self.i] < 0: return self._raw[self.i]
        return (self._raw[self.i] >> 1) - 1

    @value.setter
    def value(self, int value):
        if value < 0:
            self._raw[self.i] = value
            return
        if self.phased:
            self._raw[self.i] = (value + 1)<<1|1
        else:
            self._raw[self.i] = (value + 1)<<1

    def __repr__(self):
        if self.value < 0: return "."
        return str(self.value) + ("|" if self.phased else "/")

cdef inline Allele newAllele(int32_t *raw, int i):
    cdef Allele a = Allele.__new__(Allele)
    a._raw = raw
    a.i = i
    return a

cdef class Genotypes(object):
    cdef int32_t *_raw
    cdef readonly int n_samples
    cdef readonly int ploidy
    def __cinit__(self):
        self.ploidy = 0
        self.n_samples = 0
        self._raw = NULL
    def __dealloc__(self):
        if self._raw != NULL:
            stdlib.free(self._raw)

    def phased(self, int i):
        """
        a boolean indicating that the ith sample is phased.
        """
        return (self._raw[i * self.ploidy + 1] & 1) == 1

    def alleles(self, int i):
        cdef list result = []
        cdef int32_t v
        for j in range(self.ploidy):
            v = self._raw[i * self.ploidy + j]
            result.append((v >> 1) - 1)
        return result

    def array(Genotypes self):
        """
        array returns an int16 numpy array  of shape n_samples, (ploidy + 1).
        The last column indicates phased (1 is phased, 0 is unphased).
        The other columns indicate the alleles, e.g. [0, 1, 1] is 0|1.
        """
        cdef np.int16_t* to_return = <np.int16_t *>stdlib.malloc(sizeof(np.int16_t)
                                                   * self.n_samples
                                                   * (self.ploidy+1))

        cdef int ind
        cdef int allele
        cdef int p = self.ploidy + 1

        for ind in range(self.n_samples):
            for allele in range(self.ploidy):
                to_return[ind * p + allele] = (self._raw[ind * self.ploidy + allele] >> 1) - 1
            to_return[ind * p + self.ploidy] = (self._raw[ind * self.ploidy + 1] & 1) == 1

        cdef np.npy_intp shape[2]
        shape[0] = self.n_samples
        shape[1] = self.ploidy + 1
        return np.PyArray_SimpleNewFromData(
            2,
            shape,
            np.NPY_INT16,
            to_return
        )

    def __getitem__(self, int i):
        ## return the Allele objects for the i'th sample.
        cdef int k
        return [newAllele(self._raw, k) for k in range(i*self.ploidy,(i+1)*self.ploidy)]

cdef inline Genotypes newGenotypes(int32_t *raw, int ploidy, int n_samples):
    cdef Genotypes gs = Genotypes.__new__(Genotypes)
    gs._raw = raw
    gs.ploidy = ploidy
    gs.n_samples = n_samples
    return gs

cdef class Variant(object):
    """
    Variant represents a single VCF Record.

    It is created internally by iterating over a VCF.

    Attributes
    ----------

    INFO: `INFO`
       a dictionary-like field that provides access to the VCF INFO field.

    """
    cdef bcf1_t *b
    cdef VCF vcf
    cdef int32_t *_gt_types
    cdef int *_gt_ref_depths
    cdef int *_gt_alt_depths
    cdef float *_gt_alt_freqs
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
    cdef list _genotypes

    cdef readonly int POS

    def __init__(self, *args, **kwargs):
        raise TypeError("Variant object cannot be instantiated directly.")

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
        if self._gt_alt_freqs != NULL:
            stdlib.free(self._gt_alt_freqs)
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
        "numpy array indicating the alleles in each sample."
        def __get__(self):
            cdef np.ndarray gt_types = self.gt_types
            cdef int i, n = self.ploidy, j=-1, a, b
            cdef char **alleles = self.b.d.allele
            #cdef dict d = {i:alleles[i] for i in range(self.b.n_allele)}
            cdef list d = [from_bytes(alleles[i]) for i in range(self.b.n_allele)]
            d.append(".") # -1 gives .
            cdef list bases = ["./." for _ in range(self.vcf.n_samples)]
            cdef np.ndarray phased = self.gt_phases
            cdef list lookup = ["/", "|"]
            cdef int unknown = 3 if self.vcf.gts012 else 2
            for i in range(0, n * self.vcf.n_samples, n):
                j += 1
                if n == 2:
                    if (gt_types[j] == unknown) and (not self.vcf.strict_gt):
                        continue
                    else:
                        a = self._gt_idxs[i]
                        b = self._gt_idxs[i + 1]
                        if a >= -1 and b >= -1:
                          bases[j] = d[a] + lookup[phased[j]] + d[b]
                        else:
                          bases[j] = d[a]
                elif n == 1:
                    bases[j] = d[self._gt_idxs[i]]
                else:
                    raise Exception("gt_bases not implemented for ploidy > 2")

            return np.array(bases, np.str)

    def relatedness(self,
                    int32_t[:, ::view.contiguous] ibs,
                    int32_t[:, ::view.contiguous] n,
                    int32_t[:] hets):
        if not self.vcf.gts012:
            raise Exception("must call relatedness with gts012")
        if self._gt_types == NULL:
            self.gt_types
        cdef double[:] alt_freqs = self.gt_alt_freqs
        return krelated(<int32_t *>self._gt_types, &ibs[0, 0], &n[0, 0],
                &hets[0], self.vcf.n_samples, &alt_freqs[0])

    property num_called:
        "number of samples that were not UKNOWN."
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
        "proprtion of samples that were not UKNOWN."
        def __get__(self):
            if self.vcf.n_samples > 0:
                return float(self.num_called) / self.vcf.n_samples

    property aaf:
        "alternate allele frequency across samples in this VCF."
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
        "number homozygous reference samples at this variant."
        def __get__(self):
            if self._gt_types == NULL:
                self.gt_types
            cdef int n = 0, i = 0
            for i in range(self.vcf.n_samples):
                if self._gt_types[i] == 0:
                    n+=1
            return n

    property num_het:
        "number heterozygous samples at this variant."
        def __get__(self):
            if self._gt_types == NULL:
                self.gt_types
            cdef int n = 0, i = 0
            for i in range(self.vcf.n_samples):
                if self._gt_types[i] == 1:
                    n+=1
            return n

    property num_hom_alt:
        "number homozygous alternate samples at this variant."
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
        "number unknown samples at this variant."
        def __get__(self):
            if self._gt_types == NULL:
                self.gt_types
            cdef int n = 0, i = 0
            for i in range(self.vcf.n_samples):
                if self._gt_types[i] == 2:
                    n+=1
            return n

    property FORMAT:
        "VCF FORMAT field for this variant."
        def __get__(self):
            cdef int i
            cdef bcf_fmt_t fmt
            cdef char *key
            keys = []
            for i in range(self.b.n_fmt):
                fmt = self.b.d.fmt[i];
                key = bcf_hdr_int2id(self.vcf.hdr, BCF_DT_ID, fmt.id)
                keys.append(from_bytes(key))
            return keys

    def format(self, field, vtype=None):
        """format returns a numpy array for the requested field.

        The numpy array shape will match the requested field. E.g. if the fields
        has number=3, then the shape will be (n_samples, 3).

        Parameters
        ----------
        field: str
            FORMAT field to get the values.

        Returns
        -------
        numpy array.
        """
        if vtype is None:
            vtype = self.vcf.get_type(field)

        cdef bytes tag = to_bytes(field)
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

    @property
    def genotype(self):
        if self.vcf.n_samples == 0: return None
        cdef int32_t *gts = NULL
        cdef int ndst = 0
        if bcf_get_genotypes(self.vcf.hdr, self.b, &gts, &ndst) <= 0:
            raise Exception("couldn't get genotypes for variant")
        return newGenotypes(gts, ndst/self.vcf.n_samples, self.vcf.n_samples)

    @genotype.setter
    def genotype(self, Genotypes g):
        cdef int ret = bcf_update_genotypes(self.vcf.hdr, self.b, g._raw, self.vcf.n_samples * g.ploidy)
        if ret < 0:
            raise Exception("error setting genotypes with: %s" % g)

    property genotypes:
        """genotypes returns a list for each sample Indicating the allele and phasing.

        e.g. [0, 1, True] corresponds to 0|1
        while [1, 2, False] corresponds to 1/2
        """
        def __get__(self):
            if self.vcf.n_samples == 0: return None
            if self._genotypes is not None:
              return self._genotypes
            cdef int32_t *gts = NULL
            cdef int i, j, nret, off = 0, ndst = 0, k = 0
            cdef int n_samples = self.vcf.n_samples
            #self._genotypes = []
            nret = bcf_get_genotypes(self.vcf.hdr, self.b, &gts, &ndst)
            if nret < 0:
                raise Exception("error parsing genotypes")
            nret /= n_samples
            self._genotypes = [[] for _ in range(n_samples)]

            for i in range(n_samples):
              k = i * nret
              for j in range(nret):
                  #assert k + j < ndst
                  if bcf_gt_is_missing(gts[k + j]):
                      self._genotypes[i].append(-1)
                      continue
                  if gts[k + j] == bcf_int32_vector_end:
                      break
                  self._genotypes[i].append(bcf_gt_allele(gts[k + j]))
              self._genotypes[i].append(
                    bool(bcf_gt_is_phased(gts[k+1 if k+1 < ndst else k])))

            stdlib.free(gts)
            return self._genotypes

        def __set__(self, gts):
            cdef int n_samples = self.vcf.n_samples
            if len(gts) != n_samples:
                raise Exception("genotypes: must set with a number of gts equal the number of samples in the vcf")
            elif len(gts) == 0:
                nret = 0
            else:
                nret = max(len(gt)-1 for gt in gts)
            cdef int * cgts = <int *>stdlib.malloc(sizeof(int) * nret * n_samples)
            cdef int i, j, k
            self._genotypes = None

            for i in range(n_samples):
                k = i * nret
                for j in range(nret):
                    if j == len(gts[i]) - 1:
                        cgts[k + j] = bcf_int32_vector_end #bcf_gt_phased(-1)
                        break
                    else:
                        cgts[k + j] = bcf_gt_phased(gts[i][j]) if gts[i][-1] else bcf_gt_unphased(gts[i][j])
            ret = bcf_update_genotypes(self.vcf.hdr, self.b, cgts, n_samples * nret)
            if ret < 0:
                raise Exception("error setting genotypes with: %s" % gts)
            stdlib.free(cgts)

    def set_pos(self, int pos0):
        """
        set the POS to the given 0-based position
        """
        self.b.pos = pos0
        self.POS = self.b.pos + 1

    def set_format(self, name, np.ndarray data not None):
        """
        set the format field given by name..
        data must be a numpy array of type float or int
        """
        cdef int n_samples = self.vcf.n_samples
        if len(data) % n_samples != 0:
            raise Exception("format: len(data) must be a multiple of number of samples in vcf.")

        cdef np.ndarray[np.float32_t, mode="c"] afloat
        cdef np.ndarray[np.int32_t, mode="c"] aint

        cdef int size = data.shape[0]
        if len((<object>data).shape) > 1:
            size *= data.shape[1]

        cdef int ret
        if np.issubdtype(data.dtype, np.int):
            aint = data.astype(np.int32).reshape((size,))
            ret = bcf_update_format_int32(self.vcf.hdr, self.b, to_bytes(name), &aint[0], size)
        elif np.issubdtype(data.dtype, np.float):
            afloat = data.astype(np.float32).reshape((size,))
            ret = bcf_update_format_float(self.vcf.hdr, self.b, to_bytes(name), &afloat[0], size)
        else:
            raise Exception("format: currently only float and int numpy arrays are supported. got %s", data.dtype)
        if ret < 0:
            raise Exception("error (%d) setting format with: %s" % (ret, data[:100]))

    property gt_types:
        """gt_types returns a numpy array indicating the type of each sample.

        HOM_REF=0, HET=1. For `gts012=True` HOM_ALT=2, UKNOWN=3
        """
        def __get__(self):
            cdef int ndst = 0, ngts, n, i, nper, j = 0, k = 0
            cdef int a
            if self.vcf.n_samples == 0:
                return np.array([])
            if self._gt_types == NULL:
                self._gt_phased = <int *>stdlib.malloc(sizeof(int) * self.vcf.n_samples)
                ngts = bcf_get_genotypes(self.vcf.hdr, self.b, &self._gt_types, &ndst)
                nper = ndst / self.vcf.n_samples
                self._ploidy = nper
                self._gt_idxs = <int *>stdlib.malloc(sizeof(int) * self.vcf.n_samples * nper)
                if ndst == 0 or nper == 0:
                    return np.array([])
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
                    n = as_gts012(self._gt_types, self.vcf.n_samples, nper, self.vcf.strict_gt)
                else:
                    n = as_gts(self._gt_types, self.vcf.n_samples, nper, self.vcf.strict_gt)
            cdef np.npy_intp shape[1]
            shape[0] = <np.npy_intp> self.vcf.n_samples
            return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT32, self._gt_types)

    property ploidy:
        """get the ploidy of each sample for the given record."""
        def __get__(self):
            if self._ploidy == -1:
                self.gt_types
            return self._ploidy

    property gt_phred_ll_homref:
        """get the PL of Hom ref for each sample as a numpy array."""
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
        """get the PL of het for each sample as a numpy array."""
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
        """get the PL of hom_alt for each sample as a numpy array."""
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
        """get the count of reference reads as a numpy array."""
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
        """get the count of alternate reads as a numpy array."""
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

    property gt_alt_freqs:
        """get the freq of alternate reads as a numpy array."""
        def __get__(self):
            if self.vcf.n_samples == 0:
                return []
            t = np.array(self.gt_depths, np.float)
            a = np.array(self.gt_alt_depths, np.float)

            # for which samples are the alt or total depths unknown?
            tU = t < 0
            aU = a < 0
            # for which samples is the total depth 0?
            t0 = t == 0

            ## initialize
            alt_freq = t.astype(float)

            # default
            alt_freq[t0] = 0
            alt_freq[aU] = 0
            alt_freq[tU] = -1

            # compute the alt_freq when not unknown and no div0 error 
            clean = ~tU & ~aU & ~t0
            alt_freq[clean] = (a[clean] / t[clean])

            return alt_freq

    property gt_quals:
        """get the GQ for each sample as a numpy array."""
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
        """get the read-depth for each sample as a numpy array."""
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
        """get a boolean indicating wether each sample is phased as a numpy array."""
        def __get__(self):
            # run for side-effect
            if self._gt_phased == NULL:
                self.gt_types
            cdef np.npy_intp shape[1]
            shape[0] = <np.npy_intp> self.vcf.n_samples

            return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT32, self._gt_phased).astype(bool)

    property REF:
        "the reference allele."
        def __get__(self):
            return self.b.d.allele[0].decode()

        def __set__(self, ref):
          # int bcf_update_alleles_str(const bcf_hdr_t *hdr, bcf1_t *line, const char *alleles_string);
          alleles = (ref + "," + ",".join(self.ALT)).encode()
          if bcf_update_alleles_str(self.vcf.hdr, self.b, alleles) != 0:
            raise ValueError("couldn't set reference to:" + str(ref))

    property ALT:
        "the list of alternate alleles."
        def __get__(self):
            cdef int i
            return [self.b.d.allele[i].decode() for i in range(1, self.b.n_allele)]

        def __set__(self, alts):
          # int bcf_update_alleles_str(const bcf_hdr_t *hdr, bcf1_t *line, const char *alleles_string);
          if not isinstance(alts, list):
            alts = [alts]
          alleles = (self.REF + "," + ",".join(alts)).encode()
          if bcf_update_alleles_str(self.vcf.hdr, self.b, alleles) != 0:
            raise ValueError("couldn't set alternates to:" + str(alts))

    property is_snp:
        "boolean indicating if the variant is a SNP."
        def __get__(self):
            cdef int i
            if len(self.b.d.allele[0]) > 1: return False
            for i in range(1, self.b.n_allele):
                if not self.b.d.allele[i] in (b"A", b"C", b"G", b"T"):
                    return False
            return self.b.n_allele > 1

    property is_indel:
        "boolean indicating if the variant is an indel."
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
        "boolean indicating if the variant is a transition."
        def __get__(self):
            if len(self.ALT) != 1: return False

            if not self.is_snp: return False
            ref = self.REF
                # just one alt allele
            alt_allele = self.ALT[0]
            if ((ref == 'A' and alt_allele == 'G') or
                (ref == 'G' and alt_allele == 'A') or
                (ref == 'C' and alt_allele == 'T') or
                (ref == 'T' and alt_allele == 'C')):
                    return True
            return False

    property is_deletion:
        "boolean indicating if the variant is a deletion."
        def __get__(self):
            if len(self.ALT) > 1: return False

            if not self.is_indel: return False
            if len(self.ALT) == 0:
                return True
            alt = self.ALT[0]
            if alt is None or alt == ".":
                return True

            if len(self.REF) > len(alt):
                return True
            return False

    property is_sv:
        "boolean indicating if the variant is an SV."
        def __get__(self):
            return self.INFO.get(b'SVTYPE') is not None

    property CHROM:
        "chromosome of the variant."
        def __get__(self):
            return bcf_hdr_id2name(self.vcf.hdr, self.b.rid).decode()

    property var_type:
        "type of variant (snp/indel/sv)"
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
        "0-based start of the variant."
        def __get__(self):
            return self.b.pos

    property end:
        "end of the variant. the INFO field is parsed for SVs."
        def __get__(self):
            return self.b.pos + self.b.rlen

    property ID:
        "the value of ID from the VCF field."
        def __get__(self):
            cdef char *id = self.b.d.id
            if id == b".": return None
            return from_bytes(id)

        def __set__(self, value):
            sanitized = str(value) if value is not None else '.'
            ret = bcf_update_id(self.vcf.hdr, self.b, to_bytes(sanitized))
            if ret != 0:
                raise Exception("not able to set ID: %s", value)

    property FILTER:
        """the value of FILTER from the VCF field.

        a value of PASS in the VCF will give None for this function
        """
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
                    return from_bytes(v)
            if n == 0:
                return None
            return from_bytes(b';'.join(bcf_hdr_int2id(self.vcf.hdr, BCF_DT_ID, self.b.d.flt[i]) for i in range(n)))

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
        "the float value of QUAL from the VCF field."
        def __get__(self):
            cdef float q = self.b.qual
            if bcf_float_is_missing(q):
                return None
            return q

        def __set__(self, value):
            if value is None:
                bcf_float_set(&self.b.qual, bcf_float_missing)
            else:
                self.b.qual = value


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
        key = from_bytes(key)
        if key == "HeaderType":
            return self.type
        cdef int i
        for i in range(self.hrec.nkeys):
            if from_bytes(self.hrec.keys[i]) == key:
                return from_bytes(self.hrec.vals[i])
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
    """
    INFO is create internally by accessing `Variant.INFO`

    is acts like a dictionary where keys are expected to be in the INFO field of the Variant
    and values are typed according to what is specified in the VCF header

    Items can be deleted with del v.INFO[key] and accessed with v.INFO[key] or v.INFO.get(key)
    """
    cdef bcf_hdr_t *hdr
    cdef bcf1_t *b
    cdef int _i

    def __cinit__(INFO self):
        self._i = 0

    def __delitem__(self, okey):
        okey = to_bytes(okey)
        cdef char *key = okey
        cdef bcf_info_t *info = bcf_get_info(self.hdr, self.b, key)
        if info == NULL:
            raise KeyError(key)
        cdef int htype = bcf_hdr_id2type(self.hdr, BCF_HL_INFO, info.key)
        cdef int ret = bcf_update_info(self.hdr, self.b, key,NULL,0,htype)
        if ret != 0:
            raise Exception("error deleting %s" % key)

    def __setitem__(self, key, value):
        # only support strings for now.
        if value is True or value is False:

            ret = bcf_update_info_flag(self.hdr, self.b, to_bytes(key), b"", int(value))
            if ret != 0:
                raise Exception("not able to set flag", key, value, ret)
            return
        cdef int32_t iint
        cdef float ifloat
        if isinstance(value, int):
            iint = value
            ret = bcf_update_info_int32(self.hdr, self.b, to_bytes(key), &iint, 1)
        elif isinstance(value, float):
            ifloat = value
            ret = bcf_update_info_float(self.hdr, self.b, to_bytes(key), &ifloat, 1)
        else:
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

        # FLAG.
        if info.len == 0:
            return bcf_hdr_id2type(self.hdr, BCF_HL_INFO, info.key) == BCF_HT_FLAG

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
        try:
            return s.decode(ENC)
        except UnicodeDecodeError:
            return s.decode('utf8', errors='replace')
    return s


cdef class Writer(VCF):
    """
    Writer class makes a VCF Writer.

    Parameters
    ----------
    fname: str
        path to file
    tmpl: VCF
        a template to use to create the output header.

    Returns
    -------
    VCF object for iterating and querying.
    """
    #cdef htsFile *hts
    #cdef bcf_hdr_t *hdr
    cdef public bytes name
    cdef bint header_written
    cdef const bcf_hdr_t *ohdr

    def __init__(Writer self, fname, VCF tmpl, mode="w"):
        self.name = to_bytes(fname)
        if fname.endswith(".gz") and mode == "w":
            mode = "wz"
        if fname.endswith(".bcf") and mode == "w":
            mode = "wb"
        self.hts = hts_open(self.name, to_bytes(mode))
        if self.hts == NULL:
            raise Exception("error opening file: %s" % self.name)

        bcf_hdr_sync(tmpl.hdr)
        self.ohdr = tmpl.hdr
        self.hdr = bcf_hdr_dup(tmpl.hdr)
        bcf_hdr_sync(self.hdr)
        self.header_written = False

    @classmethod
    def from_string(Writer cls, fname, header_string, mode="w"):
        cdef Writer self = Writer.__new__(Writer)

        self.name = to_bytes(fname)
        if fname.endswith(".gz") and mode == "w":
            mode = "wz"
        if fname.endswith(".bcf") and mode == "w":
            mode = "wb"
        self.hts = hts_open(self.name, to_bytes(mode))
        cdef char *hmode = "w"
        self.hdr = bcf_hdr_init(hmode)
        if bcf_hdr_parse(self.hdr, to_bytes(header_string)) != 0:
            raise Exception("error parsing header:" + header_string)
        if bcf_hdr_sync(self.hdr) != 0:
            raise Exception("error syncing header:" + header_string)
        self.header_written = False
        self.n_samples = bcf_hdr_nsamples(self.hdr)
        return self


    def write_header(Writer self):
        bcf_hdr_write(self.hts, self.hdr)
        self.header_written = True

    def write_record(Writer self, Variant var):
        "Write the variant to the writer."
        cdef bcf_hrec_t *h
        if not self.header_written:
            self.write_header()
        if var.b.errcode == BCF_ERR_CTG_UNDEF:
            h = bcf_hdr_id2hrec(self.ohdr, BCF_DT_CTG, 0, var.b.rid)
            if h == NULL:
                raise Exception("contig %d unknown and not found in header" % var.b.rid)
            if bcf_hdr_add_hrec(self.hdr, h) < 0:
                raise Exception("error adding contig %d to header" % var.b.rid)
            bcf_hdr_sync(self.hdr)
        elif var.b.errcode != 0:
            raise Exception("variant to be written has errorcode: %d" % var.b.errcode)
        return bcf_write(self.hts, self.hdr, var.b)

    def close(Writer self):
        if not self.header_written:
            self.write_header()
        super().close()
