from libc.stdint cimport int32_t, uint32_t, int8_t, int16_t, uint8_t
import numpy as np
cimport numpy as np
np.import_array()

cdef extern from "relatedness.h":
    int related(int *gt_types, double *asum, int32_t *N, int32_t *ibs0,
                int32_t *ibs2, int32_t n_samples)
    float r_unphased(int *a_gts, int *b_gts, float f, int32_t n_samples)
    int ibd(int agt, int bgt, int run_length, float pi, int *bins, int32_t n_bins)

    int krelated(int32_t *gt_types, int32_t *ibs, int32_t *n, int32_t *hets, int32_t n_samples)

cdef extern from "helpers.h":
    int as_gts(int *gts, int num_samples, int ploidy);
    int as_gts012(int *gts, int num_samples, int ploidy);

cdef extern from "htslib/kstring.h":

    ctypedef struct kstring_t:
        size_t l, m;
        char *s;

    inline char *ks_release(kstring_t *s)

cdef extern from "htslib/hts.h":
    ctypedef struct hFILE:
        pass

    cdef union ufp:
        hFILE *hfile;

    cdef enum htsExactFormat:
        unknown_format,
        binary_format, text_format,
        sam, bam, bai, cram, crai, vcf, bcf, csi, gzi, tbi, bed

    ctypedef struct htsFormat:
        htsExactFormat format
    
    ctypedef struct htsFile:
        ufp fp
        htsFormat format





    int hts_detect_format(hFILE *fp, htsFormat *fmt);


    htsFile *hts_open(char *fn, char *mode);

    cdef int hts_verbose = 1

    ctypedef struct hts_itr_t:
        pass

    ctypedef struct hts_idx_t:
        pass

    hts_idx_t *bcf_index_load(char *fn)

    #int hts_itr_next(BGZF *fp, hts_itr_t *iter, void *r, void *data);
    void hts_itr_destroy(hts_itr_t *iter);
    void hts_idx_destroy(hts_idx_t *idx);

cdef extern from "htslib/tbx.h":

    ctypedef struct tbx_t:
        pass

    tbx_t *tbx_index_load(const char *fn);
    hts_itr_t *tbx_itr_queryi(tbx_t *tbx, int tid, int beg, int end)
    hts_itr_t *tbx_itr_querys(tbx_t *tbx, char *reg)
    int tbx_itr_next(htsFile *fp, tbx_t *tbx, hts_itr_t *iter, void *data);
    void tbx_destroy(tbx_t *tbx);

cdef extern from "htslib/vcf.h":

    ctypedef struct hts_idx_t:
        pass

    int bcf_itr_next(htsFile *, hts_itr_t* iter, bcf1_t*)
    hts_itr_t *bcf_itr_querys(hts_idx_t *, void *, char *);


    const int BCF_DT_ID = 0;
    const int BCF_DT_SAMPLE = 2;

    uint32_t bcf_float_missing = 0x7F800001;


    const int BCF_BT_NULL   = 0
    const int BCF_BT_INT8   = 1
    const int BCF_BT_INT16  = 2
    const int BCF_BT_INT32  = 3
    const int BCF_BT_FLOAT  = 5
    const int BCF_BT_CHAR   = 7

    const int bcf_str_missing = 0x07
    const int bcf_str_vector_end = 0

    const int INT8_MIN = -128
    const int INT16_MIN = -32768
    const int INT32_MIN = -2147483648

    const int bcf_int8_vector_end  = -127
    const int bcf_int16_vector_end  = -32767
    const int bcf_int32_vector_end  = -2147483647

    const int bcf_int8_missing  = -127
    const int bcf_int16_missing  = -32767
    const int bcf_int32_missing  = -2147483647


    ctypedef union uv1:
        int32_t i; # integer value
        float f;   # float value

    ctypedef struct variant_t:
        pass

    ctypedef struct bcf_fmt_t:
        int n; # n

    ctypedef struct bcf_info_t:
        int key;        # key: numeric tag id, the corresponding string is bcf_hdr_t::id[BCF_DT_ID][$key].key
        int type, len;  # type: one of BCF_BT_* types; len: vector length, 1 for scalars
        #} v1; # only set if $len==1; for easier access
        uv1 v1
        uint8_t *vptr;          # pointer to data array in bcf1_t->shared.s, excluding the size+type and tag id bytes
        uint32_t vptr_len;      # length of the vptr block or, when set, of the vptr_mod block, excluding offset
        uint32_t vptr_off;
        uint32_t vptr_free;   # vptr offset, i.e., the size of the INFO key plus size+type bytes
               # indicates that vptr-vptr_off must be freed; set only when modified and the new


    ctypedef struct bcf_dec_t:
        int m_fmt, m_info, m_id, m_als, m_allele, m_flt; # allocated size (high-water mark); do not change
        int n_flt;  # Number of FILTER fields
        int *flt;   # FILTER keys in the dictionary
        char *id;      # ID block (\0-seperated)
        char *als;     # REF+ALT block (\0-seperated)
        char **allele;      # allele[0] is the REF (allele[] pointers to the als block); all null terminated
        bcf_info_t *info;   # INFO
        bcf_fmt_t *fmt;     # FORMAT and individual sample
        variant_t *var;     # $var and $var_type set only when set_variant_types called
        int n_var, var_type;
        int shared_dirty;   # if set, shared.s must be recreated on BCF output
        int indiv_dirty;    # if set, indiv.s must be recreated on BCF output

    ctypedef struct bcf1_t:
        int32_t rid;  #// CHROM
        int32_t pos;  #// POS
        int32_t rlen; #// length of REF
        float qual;   #// QUAL
        uint32_t n_info, n_allele;
        #uint32_t n_fmt:8, n_sample:24;
        #kstring_t shared, indiv;
        bcf_dec_t d; #// lazy evaluation: $d is not generated by bcf_read(), but by explicitly calling bcf_unpack()
        int max_unpack;        # // Set to BCF_UN_STR, BCF_UN_FLT, or BCF_UN_INFO to boost performance of vcf_parse when some of the fields wont be needed
        int unpacked;          # // remember what has been unpacked to allow calling bcf_unpack() repeatedly without redoing the work
        int unpack_size[3];    # // the original block size of ID, REF+ALT and FILTER
        int errcode;   # // one of BCF_ERR_* codes

    ctypedef struct bcf_idpair_t:
        pass

    const int BCF_HL_FLT  = 0 # header line
    const int BCF_HL_INFO = 1
    const int BCF_HL_FMT  = 2
    const int BCF_HL_CTG  = 3
    const int BCF_HL_STR  = 4 # structured header line TAG=<A=..,B=..>
    const int BCF_HL_GEN  = 5 # generic header line

    ctypedef struct bcf_hrec_t:
        int type;       # One of the BCF_HL_* type
        char *key;      # The part before '=', i.e. FILTER/INFO/FORMAT/contig/fileformat etc.
        char *value;    # Set only for generic lines, NULL for FILTER/INFO, etc.
        int nkeys;              # Number of structured fields
        char **keys;    # The key=value pairs
        char **vals;    # The key=value pairs

    ctypedef struct kstring_t:
        pass

    ctypedef struct bcf_hdr_t:
        int32_t n[3];
        bcf_idpair_t *id[3];
        void *dict[3];         # ID dictionary, contig dict and sample dict
        char **samples;
        bcf_hrec_t **hrec;
        int nhrec, dirty;
        int ntransl;    # for bcf_translate()
        int *transl[2]; # for bcf_translate()
        int nsamples_ori;        # for bcf_hdr_set_samples()
        uint8_t *keep_samples;
        kstring_t mem;


    bint bcf_float_is_missing(float f)
    bint bcf_float_is_vector_end(float f)

    void bcf_destroy(bcf1_t *v);
    bcf1_t * bcf_init();
    int vcf_parse(kstring_t *s, const bcf_hdr_t *h, bcf1_t *v);

    bcf_hdr_t *bcf_hdr_read(htsFile *fp);

    int bcf_hdr_set_samples(bcf_hdr_t *hdr, const char *samples, int is_file);
    int bcf_hdr_nsamples(const bcf_hdr_t *hdr);
    void bcf_hdr_destroy(const bcf_hdr_t *hdr)
    char *bcf_hdr_fmt_text(const bcf_hdr_t *hdr, int is_bcf, int *len);

    int bcf_write(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v);
    int bcf_hdr_write(htsFile *fp, bcf_hdr_t *h);
    int vcf_format(const bcf_hdr_t *h, const bcf1_t *v, kstring_t *s);

    bcf_hrec_t *bcf_hdr_get_hrec(const bcf_hdr_t *hdr, int type, const char *key, const char *value, const char *str_class);
    void bcf_hrec_destroy(bcf_hrec_t *)

    int hts_close(htsFile *fp);

    int bcf_read(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v) nogil;

    const char *bcf_hdr_id2name(const bcf_hdr_t *hdr, int rid);
    const char *bcf_hdr_int2id(const bcf_hdr_t *hdr, int type, int int_id)
    int bcf_hdr_id2int(const bcf_hdr_t *hdr, int type, const char *id);

    int bcf_unpack(bcf1_t *b, int which) nogil;


    bcf_fmt_t *bcf_get_fmt(const bcf_hdr_t *hdr, bcf1_t *line, const char *key);

    int bcf_get_genotypes(const bcf_hdr_t *hdr, bcf1_t *line, int **dst, int *ndst);
    int bcf_get_format_int32(const bcf_hdr_t *hdr, bcf1_t *line, char * tag, int **dst, int *ndst);
    int bcf_get_format_float(const bcf_hdr_t *hdr, bcf1_t *line, char * tag, float **dst, int *ndst)
    int bcf_get_format_string(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, char ***dst, int *ndst);

    int bcf_get_format_values(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, void **dst, int *ndst, int type);
    int bcf_gt_is_phased(int);
    int bcf_gt_is_missing(int);
    int bcf_gt_allele(int);
    bint bcf_float_is_missing(float);
    bcf_info_t *bcf_get_info(const bcf_hdr_t *hdr, bcf1_t *line, const char *key);

    int bcf_update_info(const bcf_hdr_t *hdr, bcf1_t *line, const char *key, const void *values, int n, int type);
    int bcf_update_filter(const bcf_hdr_t *hdr, bcf1_t *line, int *flt_ids, int n);


    int bcf_add_id(const bcf_hdr_t *hdr, bcf1_t *line, const char *id);
    int bcf_update_id(const bcf_hdr_t *hdr, bcf1_t *line, const char *id);
    int bcf_hdr_append(bcf_hdr_t * hdr, char *);
    int bcf_hdr_sync(bcf_hdr_t *h);
    bcf_hdr_t *bcf_hdr_dup(bcf_hdr_t *h);

    int bcf_update_info_int32(const bcf_hdr_t *hdr, bcf1_t * line, const char *key, const int32_t *values, int n)
    int bcf_update_info_float(const bcf_hdr_t *hdr, bcf1_t * line, const char *key, const float *values, int n)
    int bcf_update_info_flag(const bcf_hdr_t *hdr, bcf1_t * line, const char
            *key, const char *value, int n)
    int bcf_update_info_string(const bcf_hdr_t *hdr, bcf1_t * line, const char *key, const char *values)
    #define bcf_update_info_flag(hdr,line,key,string,n)    bcf_update_info((hdr),(line),(key),(string),(n),BCF_HT_FLAG)
    #define bcf_update_info_float(hdr,line,key,values,n)   bcf_update_info((hdr),(line),(key),(values),(n),BCF_HT_REAL)
    #define bcf_update_info_flag(hdr,line,key,string,n)    bcf_update_info((hdr),(line),(key),(string),(n),BCF_HT_FLAG)
    #define bcf_update_info_string(hdr,line,key,string)    bcf_update_info((hdr),(line),(key),(string),1,BCF_HT_STR)

    # free the array, not the values.
    char **bcf_index_seqnames(hts_idx_t *idx, bcf_hdr_t *hdr, int *n);
    char **tbx_seqnames(tbx_t *tbx, int *n)
    char **bcf_hdr_seqnames(bcf_hdr_t *hdr, int *n);



