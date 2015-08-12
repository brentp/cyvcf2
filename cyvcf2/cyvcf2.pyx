#cython: profile=True
import os.path
from libc cimport stdlib
import numpy as np
cimport numpy as np
np.import_array()


cdef class VCF(object):

    cdef htsFile *hts
    cdef const bcf_hdr_t *hdr
    cdef int n_samples
    cdef int PASS

    def __init__(self, fname, mode="r"):
        if not os.path.exists(fname):
            raise Exception("bad path: %s" % fname)
        self.hts = hts_open(fname, mode)
        cdef bcf_hdr_t *hdr
        hdr = self.hdr = bcf_hdr_read(self.hts)
        assert bcf_hdr_set_samples(self.hdr, "-", 0) == 0, ("error setting samples")
        self.n_samples = bcf_hdr_nsamples(self.hdr)
        self.PASS = -1

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
            cdef int t
            return [self.hdr.samples[i] for i in range(self.n_samples)]

cdef class Variant(object):
    cdef bcf1_t *b
    cdef VCF vcf
    cdef int *_gt_types
    cdef int *_gt_ref_depths
    cdef int *_gt_alt_depths
    cdef bint *_gt_phased
    cdef float *_gt_quals
    cdef int *_int_gt_quals
    cdef int *_gt_idxs
    cdef int _gt_nper
    cdef int *_gt_pls
    cdef float *_gt_gls

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
            d[-1] = "."
            cdef list a = []
            cdef np.ndarray phased = self.gt_phases
            cdef char **lookup = ["/", "|"]
            for i in range(0, n * self.vcf.n_samples, n):
                a.append(d[self._gt_idxs[i]] + lookup[phased[j]] +
                        d[self._gt_idxs[i+1]])
            return np.array(a, np.str)

            #cdef np.npy_intp shape[1]
            #shape[0] = <np.npy_intp> self.vcf.n_samples * 2
            #return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT32, self._gt_idxs)

    property gt_types:
        def __get__(self):
            cdef int ndst, ngts, n, i, nper, j = 0
            if self._gt_types == NULL:
                self._gt_phased = <bint *>stdlib.malloc(sizeof(bint) * self.vcf.n_samples)
                ndst = 0
                ngts = bcf_get_genotypes(self.vcf.hdr, self.b, &self._gt_types, &ndst)
                nper = ndst / self.vcf.n_samples
                self._gt_idxs = <int *>stdlib.malloc(sizeof(int) * self.vcf.n_samples * nper)
                for i in range(0, ndst, nper):
                    self._gt_phased[j] = bcf_gt_is_phased(self._gt_types[i])
                    self._gt_idxs[i] = bcf_gt_allele(self._gt_types[i])
                    for k in range(i + 1, i + nper):
                        self._gt_idxs[k] = bcf_gt_allele(self._gt_types[k])
                    j += 1

                n = as_gts(self._gt_types, self.vcf.n_samples)
            cdef np.npy_intp shape[1]
            shape[0] = <np.npy_intp> self.vcf.n_samples
            return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT32, self._gt_types)

    property gt_phred_ll_homref:
        def __get__(self):
            cdef int ndst = 0, nret=0, n, i, j, nper
            if self._gt_pls == NULL and self._gt_gls == NULL:
                nret = bcf_get_format_int32(self.vcf.hdr, self.b, "PL", &self._gt_pls, &ndst)
                if nret < 0:
                    nret = bcf_get_format_float(self.vcf.hdr, self.b, "GL", &self._gt_gls, &ndst)
                    if nret < 0:
                        return []
                self._gt_nper = nret / self.vcf.n_samples
            cdef np.npy_intp shape[1]
            shape[0] = <np.npy_intp> self.vcf.n_samples
            if self._gt_pls != NULL:
                return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT32, self._gt_pls)[::3]
            else:
                gls = np.PyArray_SimpleNewFromData(1, shape, np.NPY_FLOAT32, self._gt_gls)[::3]
                gls = (-10 * gls).round().astype(np.int32)
                return gls

    property gt_phred_ll_het:
        def __get__(self):
            if self._gt_pls == NULL and self._gt_gls == NULL:
                self.gt_phred_ll_homref
            cdef np.npy_intp shape[1]
            shape[0] = <np.npy_intp> self.vcf.n_samples
            if self._gt_pls != NULL:
                return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT32, self._gt_pls)[1::3]
            else:
                gls = np.PyArray_SimpleNewFromData(1, shape, np.NPY_FLOAT32, self._gt_gls)[1::3]
                gls = (-10 * gls).round().astype(np.int32)
                return gls

    property gt_phred_ll_homalt:
        def __get__(self):
            if self._gt_pls == NULL and self._gt_gls == NULL:
                self.gt_phred_ll_homref
            cdef np.npy_intp shape[1]
            shape[0] = <np.npy_intp> self.vcf.n_samples
            if self._gt_pls != NULL:
                return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT32, self._gt_pls)[2::3]
            else:
                gls = np.PyArray_SimpleNewFromData(1, shape, np.NPY_FLOAT32, self._gt_gls)[2::3]
                gls = (-10 * gls).round().astype(np.int32)
                return gls

    property gt_ref_depths:    
        def __get__(self):
            cdef int ndst, nret = 0, n, i, j = 0
            if self._gt_ref_depths == NULL:
                ndst = 0
                # GATK
                nret = bcf_get_format_int32(self.vcf.hdr, self.b, "AD", &self._gt_ref_depths, &ndst)
                if nret > 0:
                    nper = nret / self.vcf.n_samples

                    for i in range(0, nret, nper):
                        self._gt_ref_depths[j] = self._gt_ref_depths[i]
                        j += 1
                else:
                    # Freebayes
                    # RO has to be 1:1
                    nret = bcf_get_format_int32(self.vcf.hdr, self.b, "RO", &self._gt_ref_depths, &ndst)
                    assert nret > 0
                # TODO: add new vcf standard.

            cdef np.npy_intp shape[1]
            shape[0] = <np.npy_intp> self.vcf.n_samples
            return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT32, self._gt_ref_depths)

    property gt_alt_depths:    
        def __get__(self):
            cdef int ndst, nret = 0, n, i, j = 0, k = 0
            if self._gt_alt_depths == NULL:
                ndst = 0
                # GATK
                nret = bcf_get_format_int32(self.vcf.hdr, self.b, "AD", &self._gt_alt_depths, &ndst)
                if nret > 0:
                    nper = nret / self.vcf.n_samples

                    for i in range(0, nret, nper):
                        self._gt_alt_depths[j] = self._gt_alt_depths[i+1]
                        # add up all the alt alleles
                        for k in range(2, nper):
                            self._gt_alt_depths[j] += self._gt_alt_depths[i+k]
                        j += 1
                else:
                    # Freebayes
                    nret = bcf_get_format_int32(self.vcf.hdr, self.b, "AO", &self._gt_alt_depths, &ndst)
                    assert nret > 0
                    nper = nret / self.vcf.n_samples
                    for i in range(0, nret, nper):
                        self._gt_alt_depths[j] = 0
                        for k in range(nper):
                            self._gt_alt_depths[j] = self._gt_alt_depths[i+k]
                        j += 1
                # TODO: add new vcf standard.

            cdef np.npy_intp shape[1]
            shape[0] = <np.npy_intp> self.vcf.n_samples
            return np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT32, self._gt_alt_depths)

    property gt_quals:
        def __get__(self):
            cdef int ndst = 0, nret, n, i
            cdef int *gq
            cdef np.ndarray[np.float32_t, ndim=1] a
            if self._gt_quals == NULL and self._int_gt_quals == NULL:
                nret = bcf_get_format_int32(self.vcf.hdr, self.b, "GQ", &self._int_gt_quals, &ndst)
                if nret == -2: # defined as int
                    ndst = 0
                    nret = bcf_get_format_float(self.vcf.hdr, self.b, "GQ", &self._gt_quals, &ndst)
                #    print("OK")
                if nret < 0 and nret != -2:
                    # TODO: how to fill quals? nan?
                    return np.zeros(self.vcf.n_samples, np.float32)
            cdef np.npy_intp shape[1]
            shape[0] = <np.npy_intp> self.vcf.n_samples
            if self._int_gt_quals != NULL:
                a = np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT32, self._int_gt_quals).astype(np.float32)
            else:
                a = np.PyArray_SimpleNewFromData(1, shape, np.NPY_FLOAT32, self._gt_quals)
            # this take up 10% of the total vcf parsing time. fix!!
            #a[a == -2147483648] = np.nan
            return a

    property gt_depths:
        def __get__(self):
            return self.gt_ref_depths + self.gt_alt_depths

    property gt_phases:
        def __get__(self):
            # run for side-effect 
            if self._gt_phased == NULL:
                self.gt_types

            cdef np.npy_intp shape[1]
            shape[0] = <np.npy_intp> self.vcf.n_samples
            return np.PyArray_SimpleNewFromData(1, shape, np.NPY_BOOL, self._gt_phased)


    property REF:
        def __get__(self):
            return self.b.d.allele[0]

    property ALT:
        def __get__(self):
            cdef int i
            return [self.b.d.allele[i] for i in range(1, self.b.n_allele)]

    property CHROM:
        def __get__(self):
            return bcf_hdr_id2name(self.vcf.hdr, self.b.rid)

    property var_type:
        def __get__(self):
            return "snp"

    property var_subtype:
        def __get__(self):
            return "snp"


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
            # TODO: convert "PASS" to None?
            return ';'.join(bcf_hdr_int2id(h, BCF_DT_ID, self.b.d.flt[i]) for i in range(n))

    property QUAL:
        def __get__(self):
            cdef float q = self.b.qual
            if bcf_float_is_missing(q):
                return None
            return q

    def INFO_get(self, char *key, default=None):
        cdef bcf_info_t *info = bcf_get_info(self.vcf.hdr, self.b, key)
        if info == NULL:
            return default
            raise KeyError

        if info.len == 1:
            if info.type == BCF_BT_INT8:
                if info.v1.i == INT8_MIN:
                    return default
                return <int>(info.v1.i)

            if info.type == BCF_BT_INT16:
                if info.v1.i == INT16_MIN:
                    return default
                return <int>(info.v1.i)

            if info.type == BCF_BT_INT32:
                if info.v1.i == INT32_MIN:
                    return default
                return <int>(info.v1.i)

            if info.type == BCF_BT_FLOAT:
                if bcf_float_is_missing(info.v1.f):
                    return default
                return info.v1.f

        if info.type == BCF_BT_CHAR:
            v = info.vptr[:info.vptr_len]
            if v[0] == 0x7:
                return default
            return v

        return bcf_array_to_object(info.vptr, info.type, info.len)


#    property INFO:
#        def __get__(self):
#            if self._INFO is NULL:
#                self._INFO = {}
#            if

    # var_type and var_sub_type


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


cdef Variant newVariant(bcf1_t *b, VCF vcf):
    cdef Variant v = Variant.__new__(Variant)
    v.b = b
    with nogil:
        bcf_unpack(v.b, 15)
    v.vcf = vcf
    v.POS = v.b.pos + 1
    return v
