#cython: profile=False
import os.path
import sys
from libc cimport stdlib
import numpy as np
cimport numpy as np
np.import_array()
np.seterr(invalid='ignore')

from cython cimport view

def r_(int[::view.contiguous] a_gts, int[::view.contiguous] b_gts, float f, int32_t n_samples):
    return r_unphased(&a_gts[0], &b_gts[0], f, n_samples)

cdef class VCF(object):

    cdef htsFile *hts
    cdef const bcf_hdr_t *hdr
    cdef tbx_t *idx
    cdef int n_samples
    cdef int PASS
    cdef char *fname
    cdef bint gts012
    cdef bint lazy


    def __init__(self, fname, mode="r", gts012=False, lazy=False):
        if fname == "-":
            fname = "/dev/stdin"
        if not os.path.exists(fname):
            raise Exception("bad path: %s" % fname)
        self.hts = hts_open(fname, mode)
        cdef bcf_hdr_t *hdr
        hdr = self.hdr = bcf_hdr_read(self.hts)
        assert bcf_hdr_set_samples(self.hdr, "-", 0) == 0, ("error setting samples")
        self.n_samples = bcf_hdr_nsamples(self.hdr)
        self.PASS = -1
        self.fname = fname
        self.gts012 = gts012
        self.lazy = lazy

    def __call__(VCF self, char *region):
        if self.idx == NULL:
            # we load the index on first use if possible and re-use
            if not os.path.exists(self.fname + ".tbi"):
                raise Exception("cant extraction region without tabix index")
            self.idx = tbx_index_load(self.fname + ".tbi")
            assert self.idx != NULL, "error loading index for %s" % self.fname

        cdef hts_itr_t *itr = tbx_itr_querys(self.idx, region)
        assert itr != NULL, "error starting query for %s at %s" % (self.name, region)
        cdef kstring_t s
        cdef bcf1_t *b
        cdef int slen, ret

        try:
            slen = tbx_itr_next(self.hts, self.idx, itr, &s)
            while slen > 0:
                b = bcf_init()
                ret = vcf_parse(&s, self.hdr, b)
                if ret > 0:
                    raise Exception("error parsing")
                yield newVariant(b, self)
                slen = tbx_itr_next(self.hts, self.idx, itr, &s)
        finally:
            stdlib.free(s.s)
            hts_itr_destroy(itr)

    # pull something out of the HEADER, e.g. CSQ
    def __getitem__(self, key):
        cdef bcf_hrec_t *b = bcf_hdr_get_hrec(self.hdr, BCF_HL_INFO, "ID", key, NULL);
        cdef int i
        if b == NULL:
            b = bcf_hdr_get_hrec(self.hdr, BCF_HL_GEN, key, NULL, NULL);
            if b == NULL:
                raise KeyError
            d = {b.key: b.value}
        else:
            d =  {b.keys[i]: b.vals[i] for i in range(b.nkeys)}
        #bcf_hrec_destroy(b)
        return d

    def __dealloc__(self):
        if self.hdr != NULL:
            bcf_hdr_destroy(self.hdr)
            self.hdr = NULL
        if self.hts != NULL:
            hts_close(self.hts)
            self.hts = NULL
        if self.idx != NULL:
            tbx_destroy(self.idx)

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
            return [self.hdr.samples[i] for i in range(self.n_samples)]

    property raw_header:
        def __get__(self):
            cdef int hlen
            s = bcf_hdr_fmt_text(self.hdr, 0, &hlen)
            return s

    def plot_relatedness(self, riter):
        import pandas as pd
        from matplotlib import pyplot as plt
        from matplotlib import gridspec
        import seaborn as sns

        df = []
        for row in riter:
          row['jtags'] = '|'.join(row['tags'])
          df.append(row)

        df = pd.DataFrame(df)
        fig = plt.figure(figsize=(9, 9))

        gs = gridspec.GridSpec(2, 1, height_ratios=[3.5, 1])
        colors = sns.color_palette("Set1", len(set(df.jtags)))

        ax0, ax1 = plt.subplot(gs[0]), plt.subplot(gs[1])

        for i, tag in enumerate(set(df.jtags)):
            subset = df[df.jtags == tag]
            subset.plot(kind='scatter', x='rel', y='ibs', c=colors[i],
                      label=tag, ax=ax0)

        ax0.legend()
        ax0.set_ylim(ymin=0)
        ax0.set_xlim(xmin=-0.2)

        ax1.set_xlim(*ax0.get_xlim())
        ax1.hist(df.rel, 60)
        ax1.set_yscale('log', nonposy='clip')
        return fig

    def relatedness(self, int n_variants=3000, int gap=30000, float min_af=0.02,
                    float max_af=0.6, float linkage_max=0.05):

        cdef Variant v

        cdef int last = -gap, nv = 0, nvt=0
        cdef int *last_gts
        samples = self.samples
        cdef int n_samples = len(samples)
        cdef float aaf
        cdef int n_unlinked = 0

        a = np.zeros((n_samples, n_samples), np.float64)
        n = np.zeros((n_samples, n_samples), np.int32)
        ibs0 = np.zeros((n_samples, n_samples), np.int32)

        for v in self:
            nvt += 1
            if last_gts == NULL:
                if v._gt_types == NULL:
                    v.gt_types
                last_gts = v._gt_types
            if v.POS - last < gap and v.POS > last:
                continue
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
                v.gt_types
            last, last_gts = v.POS, v._gt_types

            v.relatedness(a, n, ibs0)
            nv += 1
            if nv == n_variants:
                break
        sys.stderr.write("tested: %d variants out of %d\n" % (nv, nvt))

        a /= n
        ibs0 = ibs0.astype(float) / n.astype(float)
        for sj, sample_j in enumerate(samples):
            for sk, sample_k in enumerate(samples[sj:], start=sj):
                if sj == sk: continue

                rel, ibs = a[sj, sk], ibs0[sj, sk]
                pair = sample_j, sample_k

                d = {'pair': pair, 'rel': rel, 'ibs': ibs}

                # self or twin
                if rel > 0.8:
                    d['tags'] = ['identical twins', 'self']
                elif rel > 0.7:
                    opts = ['identical twins', 'self']
                    if ibs < 0.012:
                        opts += ['parent-child']
                    else:
                        opts += ['full-siblings']
                    d['tags'] = opts

                # sib or parent-child
                elif 0.3 < rel < 0.7:
                    if ibs > 0.018:

                        d['tags'] = ['full siblings']
                    elif ibs < 0.012:
                        d['tags'] = ['parent-child']
                    else:
                        d['tags'] = ['parent-child', 'full siblings']
                elif 0.15 < rel < 0.3:
                    d['tags'] = ['related level 2']
                elif rel < 0.04:
                    d['tags'] = ['unrelated']
                elif rel < 0.15:
                    d['tags'] = ['distant relations', 'unrelated']
                else:
                    raise Exception('impossible')
                yield d

cdef class Variant(object):
    cdef bcf1_t *b
    cdef VCF vcf
    cdef int *_gt_types
    cdef int *_gt_ref_depths
    cdef int *_gt_alt_depths
    cdef int *_gt_phased
    cdef float *_gt_quals
    cdef int *_int_gt_quals
    cdef int *_gt_idxs
    cdef int _gt_nper
    cdef int *_gt_pls
    cdef float *_gt_gls
    cdef readonly INFO INFO

    cdef readonly int POS

    def __cinit__(self):
        self.b = NULL
        self._gt_types = NULL
        self._gt_phased = NULL
        self._gt_pls

    def __repr__(self):
        return "Variant(%s:%d %s/%s)" % (self.CHROM, self.POS, self.REF,
                ",".join(self.ALT))

    def __str__(self):
        cdef kstring_t s
        s.s, s.l, s.m = NULL, 0, 0
        vcf_format(self.vcf.hdr, self.b, &s)
        return ks_release(&s)

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
            if self._gt_idxs == NULL:
                self.gt_types
            cdef int i, n = self.b.n_allele, j=0
            cdef char **alleles = self.b.d.allele
            cdef dict d = {i:alleles[i] for i in range(n)}
            cdef list a = []
            cdef np.ndarray phased = self.gt_phases
            cdef char **lookup = ["/", "|"]
            for i in range(0, n * self.vcf.n_samples, n):
                a.append(d.get(self._gt_idxs[i], ".")
                        + lookup[phased[j]] +
                        d.get(self._gt_idxs[i+1], "."))
                j += 1
            return np.array(a, np.str)

    cpdef relatedness(self, double[:, ::view.contiguous] asum,
                          int32_t[:, ::view.contiguous] n,
                          int32_t[:, ::view.contiguous] ibs0):
        if not self.vcf.gts012:
            raise Exception("must call relatedness with gts012")
        if self._gt_types == NULL:
            self.gt_types
        cdef int n_samples = self.vcf.n_samples
        return related(self._gt_types, &asum[0, 0], &n[0, 0], &ibs0[0, 0], n_samples)

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

    property gt_types:
        def __get__(self):
            cdef int ndst, ngts, n, i, nper, j = 0, k = 0
            if self.vcf.n_samples == 0:
                return []
            if self._gt_types == NULL:
                self._gt_phased = <int *>stdlib.malloc(sizeof(int) * self.vcf.n_samples)
                ndst = 0
                ngts = bcf_get_genotypes(self.vcf.hdr, self.b, &self._gt_types, &ndst)
                nper = ndst / self.vcf.n_samples
                self._gt_idxs = <int *>stdlib.malloc(sizeof(int) * self.vcf.n_samples * nper)
                for i in range(0, ndst, nper):
                    self._gt_idxs[i] = bcf_gt_allele(self._gt_types[i])
                    for k in range(i + 1, i + nper):
                        self._gt_idxs[k] = bcf_gt_allele(self._gt_types[k])
                    self._gt_phased[j] = <int>bcf_gt_is_phased(self._gt_types[i+1])
                    j += 1

                if self.vcf.gts012:
                    n = as_gts012(self._gt_types, self.vcf.n_samples)
                else:
                    n = as_gts(self._gt_types, self.vcf.n_samples)
            #print self.vcf.fname, self.POS, [self._gt_phased[i] for i in range(self.vcf.n_samples)]
            cdef np.npy_intp shape[1]
            shape[0] = <np.npy_intp> self.vcf.n_samples
            return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT32, self._gt_types)

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
            return self.b.d.allele[0]

    property ALT:
        def __get__(self):
            cdef int i
            return [self.b.d.allele[i] for i in range(1, self.b.n_allele)]

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
            return self.INFO.get('SVTYPE') is not None

    property CHROM:
        def __get__(self):
            return bcf_hdr_id2name(self.vcf.hdr, self.b.rid)

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
            cdef bcf_hdr_t *h = self.vcf.hdr
            cdef int n = self.b.d.n_flt
            if n == 0:
                return None
            if n == 1:
                if self.vcf.PASS != -1:
                    if self.b.d.flt[0] == self.vcf.PASS:
                        return None
                else:
                    v = bcf_hdr_int2id(h, BCF_DT_ID, self.b.d.flt[0])
                    if v == b"PASS":
                        self.vcf.PASS = self.b.d.flt[0]
                        return None
                    return v
            return ';'.join(bcf_hdr_int2id(h, BCF_DT_ID, self.b.d.flt[i]) for i in range(n))

    property QUAL:
        def __get__(self):
            cdef float q = self.b.qual
            if bcf_float_is_missing(q):
                return None
            return q

cdef class INFO(object):
    cdef bcf_hdr_t *hdr
    cdef bcf1_t *b
    cdef int _i

    def __cinit__(INFO self):
        self._i = 0

    cdef _getval(INFO self, bcf_info_t * info):

        if info.len == 1:
            if info.type == BCF_BT_INT8:
                if info.v1.i == INT8_MIN:
                    raise KeyError
                return <int>(info.v1.i)

            if info.type == BCF_BT_INT16:
                if info.v1.i == INT16_MIN:
                    raise KeyError
                return <int>(info.v1.i)

            if info.type == BCF_BT_INT32:
                if info.v1.i == INT32_MIN:
                    raise KeyError
                return <int>(info.v1.i)

            if info.type == BCF_BT_FLOAT:
                if bcf_float_is_missing(info.v1.f):
                    raise KeyError
                return info.v1.f

        if info.type == BCF_BT_CHAR:
            v = info.vptr[:info.vptr_len]
            if v[0] == 0x7:
                raise KeyError
            return v

        return bcf_array_to_object(info.vptr, info.type, info.len)

    def __getitem__(self, okey):
        okey = str(okey)
        cdef char *key = okey
        cdef bcf_info_t *info = bcf_get_info(self.hdr, self.b, key)
        if info == NULL:
            raise KeyError
        return self._getval(info)

    def get(self, char *key, default=None):
        try:
            return self.__getitem__(key)
        except KeyError:
            return None

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
        return name, self._getval(info)

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
        value = datac[:n] if datac[0] != bcf_str_missing else None
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

cdef class Writer(object):
    cdef htsFile *hts
    cdef bcf_hdr_t *hdr
    cdef public str name
    cdef bint header_written

    # TODO: see cgotabix AddHeaderInfo for how to add stuff to header.

    def __init__(self, fname, VCF tmpl):
        self.name = fname
        self.hts = hts_open(fname, "w")
        cdef bcf_hdr_t *h = tmpl.hdr
        cdef bcf_hdr_t *hdup = bcf_hdr_dup(h)
        self.hdr = hdup
        self.header_written = False

    def write_record(self, Variant var):
        if not self.header_written:
            bcf_hdr_write(self.hts, self.hdr)
            self.header_written = True
        return bcf_write(self.hts, self.hdr, var.b)

    def close(self):
        if self.hts != NULL:
            hts_close(self.hts)

    def __dealloc__(self):
        bcf_hdr_destroy(self.hdr)
        self.hdr = NULL
        self.close()
        self.hts = NULL
