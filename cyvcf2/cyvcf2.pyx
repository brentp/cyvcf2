#cython: profile=False
import os.path
from libc cimport stdlib
import numpy as np
cimport numpy as np
np.import_array()
np.seterr(invalid='ignore')

import atexit

cdef class VCF(object):

    cdef htsFile *hts
    cdef const bcf_hdr_t *hdr
    cdef int n_samples
    cdef int PASS
    cdef char *fname

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

    def __init__(self, fname, mode="r"):
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

    def __dealloc__(self):
        if self.hdr != NULL:
            bcf_hdr_destroy(self.hdr)
            self.hdr = NULL
        if self.hts != NULL:
            hts_close(self.hts)
            self.hts = NULL

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

            #cdef np.npy_intp shape[1]
            #shape[0] = <np.npy_intp> self.vcf.n_samples * 2
            #return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT32, self._gt_idxs)
    def relatedness(self, asum, n, vmax=3):
        self._relatedness(asum, n, vmax)

    cdef inline void _relatedness(self, double[:, ::1] asum, int[:, ::1] N, double vmax):
        """
        N track the number of snps that were used for each sample. The original
        used a fixed N for all samples due to GWAS, but for sequencing, we must
        allow for missing N.
        """
        cdef int j, k
        cdef float denom, numer, val
        cdef np.ndarray[dtype=np.float_t, ndim=1] x = np.array(self.gt_types, dtype=np.float64)
        # flip the alleles as they did in the paper. so homalt is 0
        # also flip the aaf or it doesn't work out.
        cdef float pi = self.aaf
        if pi == 0 or pi == 1: return
        x[x == 2] = -1
        x[x == 3] = 2

        denom = 2.0 * pi * (1.0 - pi)
        for j in range(self.vcf.n_samples):
            if x[j] < 0:
                continue
            for k in range(j, self.vcf.n_samples):
                if x[k] < 0:
                    continue
                N[j, k] += 1
                if j != k:
                    numer = (x[j] - 2.0 * pi) * (x[k] - 2.0 * pi)
                    val = numer / denom
                else:
                    numer = (x[j]**2.0 - (1.0 + 2.0 * pi) * x[j] + 2.0 * pi**2.0)
                    val = numer / denom

                if val < vmax:
                    asum[j, k] += val
                else:
                    asum[j, k] += vmax

    property num_called:
        def __get__(self):
            if self._gt_types == NULL:
                self.gt_types
            cdef int n = 0, i = 0
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
                import sys
                for i in range(0, ndst, nper):
                    self._gt_idxs[i] = bcf_gt_allele(self._gt_types[i])
                    for k in range(i + 1, i + nper):
                        self._gt_idxs[k] = bcf_gt_allele(self._gt_types[k])
                    self._gt_phased[j] = <int>bcf_gt_is_phased(self._gt_types[i+1])
                    j += 1

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
            if self._gt_pls == NULL and self._gt_gls == NULL:
                nret = bcf_get_format_int32(self.vcf.hdr, self.b, "PL", &self._gt_pls, &ndst)
                if nret < 0:
                    nret = bcf_get_format_float(self.vcf.hdr, self.b, "GL", &self._gt_gls, &ndst)
                    if nret < 0:
                        return []
                    else:
                        for i in range(nret):
                            if self._gt_gls[i] < -2147483646:
                                # this gets translated to -1 when converted to
                                # pls
                                self._gt_gls[i] = 0.1
                else:
                    for i in range(nret):
                        if self._gt_pls[i] < 0:
                            self._gt_pls[i] = -1

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
                if not self.b.d.allele[i] in "ACGT":
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
        import sys
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
    with nogil:
        bcf_unpack(v.b, 15)
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

    def __init__(self, fname, VCF tmpl):
        self.name = fname
        self.hts = hts_open(fname, "w")
        cdef bcf_hdr_t *h = tmpl.hdr
        self.hdr = h
        bcf_hdr_write(self.hts, h)

    def write_record(self, Variant var):
        return bcf_write(self.hts, self.hdr, var.b)

    def close(self):
        if self.hts != NULL:
            hts_close(self.hts)
