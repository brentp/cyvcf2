# 1 "htslib/sam.h"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "htslib/sam.h"
# 30 "htslib/sam.h"
# 1 "/usr/lib/gcc/x86_64-linux-gnu/5/include/stdint.h" 1 3 4
# 9 "/usr/lib/gcc/x86_64-linux-gnu/5/include/stdint.h" 3 4
# 1 "/usr/include/stdint.h" 1 3 4
# 25 "/usr/include/stdint.h" 3 4
# 1 "/usr/include/features.h" 1 3 4
# 367 "/usr/include/features.h" 3 4
# 1 "/usr/include/x86_64-linux-gnu/sys/cdefs.h" 1 3 4
# 410 "/usr/include/x86_64-linux-gnu/sys/cdefs.h" 3 4
# 1 "/usr/include/x86_64-linux-gnu/bits/wordsize.h" 1 3 4
# 411 "/usr/include/x86_64-linux-gnu/sys/cdefs.h" 2 3 4
# 368 "/usr/include/features.h" 2 3 4
# 391 "/usr/include/features.h" 3 4
# 1 "/usr/include/x86_64-linux-gnu/gnu/stubs.h" 1 3 4
# 10 "/usr/include/x86_64-linux-gnu/gnu/stubs.h" 3 4
# 1 "/usr/include/x86_64-linux-gnu/gnu/stubs-64.h" 1 3 4
# 11 "/usr/include/x86_64-linux-gnu/gnu/stubs.h" 2 3 4
# 392 "/usr/include/features.h" 2 3 4
# 26 "/usr/include/stdint.h" 2 3 4
# 1 "/usr/include/x86_64-linux-gnu/bits/wchar.h" 1 3 4
# 27 "/usr/include/stdint.h" 2 3 4
# 1 "/usr/include/x86_64-linux-gnu/bits/wordsize.h" 1 3 4
# 28 "/usr/include/stdint.h" 2 3 4
# 36 "/usr/include/stdint.h" 3 4

# 36 "/usr/include/stdint.h" 3 4
typedef signed char int8_t;
typedef short int int16_t;
typedef int int32_t;

typedef long int int64_t;







typedef unsigned char uint8_t;
typedef unsigned short int uint16_t;

typedef unsigned int uint32_t;



typedef unsigned long int uint64_t;
# 65 "/usr/include/stdint.h" 3 4
typedef signed char int_least8_t;
typedef short int int_least16_t;
typedef int int_least32_t;

typedef long int int_least64_t;






typedef unsigned char uint_least8_t;
typedef unsigned short int uint_least16_t;
typedef unsigned int uint_least32_t;

typedef unsigned long int uint_least64_t;
# 90 "/usr/include/stdint.h" 3 4
typedef signed char int_fast8_t;

typedef long int int_fast16_t;
typedef long int int_fast32_t;
typedef long int int_fast64_t;
# 103 "/usr/include/stdint.h" 3 4
typedef unsigned char uint_fast8_t;

typedef unsigned long int uint_fast16_t;
typedef unsigned long int uint_fast32_t;
typedef unsigned long int uint_fast64_t;
# 119 "/usr/include/stdint.h" 3 4
typedef long int intptr_t;


typedef unsigned long int uintptr_t;
# 134 "/usr/include/stdint.h" 3 4
typedef long int intmax_t;
typedef unsigned long int uintmax_t;
# 10 "/usr/lib/gcc/x86_64-linux-gnu/5/include/stdint.h" 2 3 4
# 31 "htslib/sam.h" 2
# 1 "htslib/hts.h" 1
# 31 "htslib/hts.h"
# 1 "/usr/lib/gcc/x86_64-linux-gnu/5/include/stddef.h" 1 3 4
# 149 "/usr/lib/gcc/x86_64-linux-gnu/5/include/stddef.h" 3 4
typedef long int ptrdiff_t;
# 216 "/usr/lib/gcc/x86_64-linux-gnu/5/include/stddef.h" 3 4
#typedef long unsigned int size_t;
# 328 "/usr/lib/gcc/x86_64-linux-gnu/5/include/stddef.h" 3 4
typedef int wchar_t;
# 426 "/usr/lib/gcc/x86_64-linux-gnu/5/include/stddef.h" 3 4
typedef struct {
  long long __max_align_ll __attribute__((__aligned__(__alignof__(long long))));
  long double __max_align_ld __attribute__((__aligned__(__alignof__(long double))));
} max_align_t;
# 32 "htslib/hts.h" 2


# 1 "htslib/hts_defs.h" 1
# 35 "htslib/hts.h" 2
# 1 "htslib/hts_log.h" 1
# 37 "htslib/hts_log.h"

# 37 "htslib/hts_log.h"
enum htsLogLevel {
    HTS_LOG_OFF,
    HTS_LOG_ERROR,
    HTS_LOG_WARNING = 3,
    HTS_LOG_INFO,
    HTS_LOG_DEBUG,
    HTS_LOG_TRACE
};


void hts_set_log_level(enum htsLogLevel level);


enum htsLogLevel hts_get_log_level();






extern int hts_verbose;
# 69 "htslib/hts_log.h"
void hts_log(enum htsLogLevel severity, const char *context, const char *format, ...)
__attribute__((__format__ (printf, 3, 4)));
# 36 "htslib/hts.h" 2






typedef struct BGZF BGZF;


struct cram_fd;
struct hFILE;
struct hts_tpool;



typedef struct __kstring_t {
    size_t l, m;
    char *s;
} kstring_t;
# 128 "htslib/hts.h"
enum htsFormatCategory {
    unknown_category,
    sequence_data,
    variant_data,
    index_file,
    region_list,
    category_maximum = 32767
};

enum htsExactFormat {
    unknown_format,
    binary_format, text_format,
    sam, bam, bai, cram, crai, vcf, bcf, csi, gzi, tbi, bed,
    json,
    format_maximum = 32767
};

enum htsCompression {
    no_compression, gzip, bgzf, custom,
    compression_maximum = 32767
};

typedef struct htsFormat {
    enum htsFormatCategory category;
    enum htsExactFormat format;
    struct { short major, minor; } version;
    enum htsCompression compression;
    short compression_level;
    void *specific;
} htsFormat;
# 166 "htslib/hts.h"
typedef struct {
    uint32_t is_bin:1, is_write:1, is_be:1, is_cram:1, is_bgzf:1, dummy:27;
    int64_t lineno;
    kstring_t line;
    char *fn, *fn_aux;
    union {
        BGZF *bgzf;
        struct cram_fd *cram;
        struct hFILE *hfile;
    } fp;
    htsFormat format;
} htsFile;
# 186 "htslib/hts.h"
typedef struct {
    struct hts_tpool *pool;
    int qsize;
} htsThreadPool;


enum sam_fields {
    SAM_QNAME = 0x00000001,
    SAM_FLAG = 0x00000002,
    SAM_RNAME = 0x00000004,
    SAM_POS = 0x00000008,
    SAM_MAPQ = 0x00000010,
    SAM_CIGAR = 0x00000020,
    SAM_RNEXT = 0x00000040,
    SAM_PNEXT = 0x00000080,
    SAM_TLEN = 0x00000100,
    SAM_SEQ = 0x00000200,
    SAM_QUAL = 0x00000400,
    SAM_AUX = 0x00000800,
    SAM_RGAUX = 0x00001000,
};


enum hts_fmt_option {

    CRAM_OPT_DECODE_MD,
    CRAM_OPT_PREFIX,
    CRAM_OPT_VERBOSITY,
    CRAM_OPT_SEQS_PER_SLICE,
    CRAM_OPT_SLICES_PER_CONTAINER,
    CRAM_OPT_RANGE,
    CRAM_OPT_VERSION,
    CRAM_OPT_EMBED_REF,
    CRAM_OPT_IGNORE_MD5,
    CRAM_OPT_REFERENCE,
    CRAM_OPT_MULTI_SEQ_PER_SLICE,
    CRAM_OPT_NO_REF,
    CRAM_OPT_USE_BZIP2,
    CRAM_OPT_SHARED_REF,
    CRAM_OPT_NTHREADS,
    CRAM_OPT_THREAD_POOL,
    CRAM_OPT_USE_LZMA,
    CRAM_OPT_USE_RANS,
    CRAM_OPT_REQUIRED_FIELDS,
    CRAM_OPT_LOSSY_NAMES,
    CRAM_OPT_BASES_PER_SLICE,


    HTS_OPT_COMPRESSION_LEVEL = 100,
    HTS_OPT_NTHREADS,
    HTS_OPT_THREAD_POOL,
    HTS_OPT_CACHE_SIZE,
};




typedef struct hts_opt {
    char *arg;
    enum hts_fmt_option opt;
    union {
        int i;
        char *s;
    } val;
    struct hts_opt *next;
} hts_opt;
# 265 "htslib/hts.h"
int hts_opt_add(hts_opt **opts, const char *c_arg);







int hts_opt_apply(htsFile *fp, hts_opt *opts);




void hts_opt_free(hts_opt *opts);
# 288 "htslib/hts.h"
int hts_parse_format(htsFormat *opt, const char *str);
# 301 "htslib/hts.h"
int hts_parse_opt_list(htsFormat *opt, const char *str);






extern const unsigned char seq_nt16_table[256];




extern const char seq_nt16_str[];




extern const int seq_nt16_int[];






const char *hts_version(void);







int hts_detect_format(struct hFILE *fp, htsFormat *fmt);






char *hts_format_description(const htsFormat *format);
# 369 "htslib/hts.h"
htsFile *hts_open(const char *fn, const char *mode);
# 385 "htslib/hts.h"
htsFile *hts_open_format(const char *fn, const char *mode, const htsFormat *fmt);






htsFile *hts_hopen(struct hFILE *fp, const char *fn, const char *mode);






int hts_close(htsFile *fp);






const htsFormat *hts_get_format(htsFile *fp);






const char *hts_format_file_extension(const htsFormat *format);
# 422 "htslib/hts.h"
int hts_set_opt(htsFile *fp, enum hts_fmt_option opt, ...);

int hts_getline(htsFile *fp, int delimiter, kstring_t *str);
char **hts_readlines(const char *fn, int *_n);
# 434 "htslib/hts.h"
char **hts_readlist(const char *fn, int is_file, int *_n);
# 444 "htslib/hts.h"
int hts_set_threads(htsFile *fp, int n);







int hts_set_thread_pool(htsFile *fp, htsThreadPool *p);







void hts_set_cache_size(htsFile *fp, int n);
# 469 "htslib/hts.h"
int hts_set_fai_filename(htsFile *fp, const char *fn_aux);
# 482 "htslib/hts.h"
int hts_check_EOF(htsFile *fp);
# 508 "htslib/hts.h"
struct __hts_idx_t;
typedef struct __hts_idx_t hts_idx_t;

typedef struct {
    uint64_t u, v;
} hts_pair64_t;

typedef int hts_readrec_func(BGZF *fp, void *data, void *r, int *tid, int *beg, int *end);

typedef struct {
    uint32_t read_rest:1, finished:1, is_cram:1, dummy:29;
    int tid, beg, end, n_off, i;
    int curr_tid, curr_beg, curr_end;
    uint64_t curr_off;
    hts_pair64_t *off;
    hts_readrec_func *readrec;
    struct {
        int n, m;
        int *a;
    } bins;
} hts_itr_t;




    hts_idx_t *hts_idx_init(int n, int fmt, uint64_t offset0, int min_shift, int n_lvls);
    void hts_idx_destroy(hts_idx_t *idx);
    int hts_idx_push(hts_idx_t *idx, int tid, int beg, int end, uint64_t offset, int is_mapped);
    void hts_idx_finish(hts_idx_t *idx, uint64_t final_offset);







int hts_idx_save(const hts_idx_t *idx, const char *fn, int fmt) __attribute__ ((__warn_unused_result__));
# 553 "htslib/hts.h"
int hts_idx_save_as(const hts_idx_t *idx, const char *fn, const char *fnidx, int fmt) __attribute__ ((__warn_unused_result__));







hts_idx_t *hts_idx_load(const char *fn, int fmt);






hts_idx_t *hts_idx_load2(const char *fn, const char *fnidx);
# 582 "htslib/hts.h"
uint8_t *hts_idx_get_meta(hts_idx_t *idx, uint32_t *l_meta);
# 596 "htslib/hts.h"
int hts_idx_set_meta(hts_idx_t *idx, uint32_t l_meta, uint8_t *meta, int is_copy);

    int hts_idx_get_stat(const hts_idx_t* idx, int tid, uint64_t* mapped, uint64_t* unmapped);
    uint64_t hts_idx_get_n_no_coor(const hts_idx_t* idx);
# 616 "htslib/hts.h"
long long hts_parse_decimal(const char *str, char **strend, int flags);
# 625 "htslib/hts.h"
const char *hts_parse_reg(const char *str, int *beg, int *end);

    hts_itr_t *hts_itr_query(const hts_idx_t *idx, int tid, int beg, int end, hts_readrec_func *readrec);
    void hts_itr_destroy(hts_itr_t *iter);

    typedef int (*hts_name2id_f)(void*, const char*);
    typedef const char *(*hts_id2name_f)(void*, int);
    typedef hts_itr_t *hts_itr_query_func(const hts_idx_t *idx, int tid, int beg, int end, hts_readrec_func *readrec);

    hts_itr_t *hts_itr_querys(const hts_idx_t *idx, const char *reg, hts_name2id_f getid, void *hdr, hts_itr_query_func *itr_query, hts_readrec_func *readrec);
    int hts_itr_next(BGZF *fp, hts_itr_t *iter, void *r, void *data) __attribute__ ((__warn_unused_result__));
    const char **hts_idx_seqnames(const hts_idx_t *idx, int *n, hts_id2name_f getid, void *hdr);
# 650 "htslib/hts.h"
    int hts_file_type(const char *fname);






struct errmod_t;
typedef struct errmod_t errmod_t;

errmod_t *errmod_init(double depcorr);
void errmod_destroy(errmod_t *em);







int errmod_cal(const errmod_t *em, int n, int m, uint16_t *bases, float *q);






typedef struct probaln_par_t {
    float d, e;
    int bw;
} probaln_par_t;

int probaln_glocal(const uint8_t *ref, int l_ref, const uint8_t *query, int l_query, const uint8_t *iqual, const probaln_par_t *c, int *state, uint8_t *q);






    struct hts_md5_context;
    typedef struct hts_md5_context hts_md5_context;
# 705 "htslib/hts.h"
    hts_md5_context *hts_md5_init(void);


    void hts_md5_update(hts_md5_context *ctx, const void *data, unsigned long size);


    void hts_md5_final(unsigned char *digest, hts_md5_context *ctx);




    void hts_md5_reset(hts_md5_context *ctx);




    void hts_md5_hex(char *hex, const unsigned char *digest);


    void hts_md5_destroy(hts_md5_context *ctx);


static inline int hts_reg2bin(int64_t beg, int64_t end, int min_shift, int n_lvls)
{
    int l, s = min_shift, t = ((1<<((n_lvls<<1) + n_lvls)) - 1) / 7;
    for (--end, l = n_lvls; l > 0; --l, s += 3, t -= 1<<((l<<1)+l))
        if (beg>>s == end>>s) return t + (beg>>s);
    return 0;
}

static inline int hts_bin_bot(int bin, int n_lvls)
{
    int l, b;
    for (l = 0, b = bin; b; ++l, b = (((b) - 1) >> 3));
    return (bin - (((1<<(((l)<<1) + (l))) - 1) / 7)) << (n_lvls - l) * 3;
}





static inline int ed_is_big(void)
{
    long one= 1;
    return !(*((char *)(&one)));
}
static inline uint16_t ed_swap_2(uint16_t v)
{
    return (uint16_t)(((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8));
}
static inline void *ed_swap_2p(void *x)
{
    *(uint16_t*)x = ed_swap_2(*(uint16_t*)x);
    return x;
}
static inline uint32_t ed_swap_4(uint32_t v)
{
    v = ((v & 0x0000FFFFU) << 16) | (v >> 16);
    return ((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8);
}
static inline void *ed_swap_4p(void *x)
{
    *(uint32_t*)x = ed_swap_4(*(uint32_t*)x);
    return x;
}
static inline uint64_t ed_swap_8(uint64_t v)
{
    v = ((v & 0x00000000FFFFFFFFLLU) << 32) | (v >> 32);
    v = ((v & 0x0000FFFF0000FFFFLLU) << 16) | ((v & 0xFFFF0000FFFF0000LLU) >> 16);
    return ((v & 0x00FF00FF00FF00FFLLU) << 8) | ((v & 0xFF00FF00FF00FF00LLU) >> 8);
}
static inline void *ed_swap_8p(void *x)
{
    *(uint64_t*)x = ed_swap_8(*(uint64_t*)x);
    return x;
}
# 32 "htslib/sam.h" 2
# 51 "htslib/sam.h"
typedef struct {
    int32_t n_targets, ignore_sam_err;
    uint32_t l_text;
    uint32_t *target_len;
    int8_t *cigar_tab;
    char **target_name;
    char *text;
    void *sdict;
} bam_hdr_t;
# 154 "htslib/sam.h"
typedef struct {
    int32_t tid;
    int32_t pos;
    uint16_t bin;
    uint8_t qual;
    uint8_t l_qname;
    uint16_t flag;
    uint8_t unused1;
    uint8_t l_extranul;
    uint32_t n_cigar;
    int32_t l_qseq;
    int32_t mtid;
    int32_t mpos;
    int32_t isize;
} bam1_core_t;
# 187 "htslib/sam.h"
typedef struct {
    bam1_core_t core;
    int l_data;
    uint32_t m_data;
    uint8_t *data;

    uint64_t id;

} bam1_t;
# 270 "htslib/sam.h"
    bam_hdr_t *bam_hdr_init(void);
    bam_hdr_t *bam_hdr_read(BGZF *fp);
    int bam_hdr_write(BGZF *fp, const bam_hdr_t *h) __attribute__ ((__warn_unused_result__));
    void bam_hdr_destroy(bam_hdr_t *h);
    int bam_name2id(bam_hdr_t *h, const char *ref);
    bam_hdr_t* bam_hdr_dup(const bam_hdr_t *h0);

    bam1_t *bam_init1(void);
    void bam_destroy1(bam1_t *b);
    int bam_read1(BGZF *fp, bam1_t *b) __attribute__ ((__warn_unused_result__));
    int bam_write1(BGZF *fp, const bam1_t *b) __attribute__ ((__warn_unused_result__));
    bam1_t *bam_copy1(bam1_t *bdst, const bam1_t *bsrc);
    bam1_t *bam_dup1(const bam1_t *bsrc);

    int bam_cigar2qlen(int n_cigar, const uint32_t *cigar);
    int bam_cigar2rlen(int n_cigar, const uint32_t *cigar);
# 298 "htslib/sam.h"
    int32_t bam_endpos(const bam1_t *b);

    int bam_str2flag(const char *str);
    char *bam_flag2str(int flag);
# 324 "htslib/sam.h"
hts_idx_t *sam_index_load(htsFile *fp, const char *fn);







hts_idx_t *sam_index_load2(htsFile *fp, const char *fn, const char *fnidx);
# 341 "htslib/sam.h"
int sam_index_build(const char *fn, int min_shift) __attribute__ ((__warn_unused_result__));
# 350 "htslib/sam.h"
int sam_index_build2(const char *fn, const char *fnidx, int min_shift) __attribute__ ((__warn_unused_result__));
int sam_index_build3(const char *fn, const char *fnidx, int min_shift, int nthreads) __attribute__ ((__warn_unused_result__));


    hts_itr_t *sam_itr_queryi(const hts_idx_t *idx, int tid, int beg, int end);
    hts_itr_t *sam_itr_querys(const hts_idx_t *idx, bam_hdr_t *hdr, const char *region);
# 366 "htslib/sam.h"
    int sam_open_mode(char *mode, const char *fn, const char *format);




    char *sam_open_mode_opts(const char *fn,
                             const char *mode,
                             const char *format);

    typedef htsFile samFile;
    bam_hdr_t *sam_hdr_parse(int l_text, const char *text);
    bam_hdr_t *sam_hdr_read(samFile *fp);
    int sam_hdr_write(samFile *fp, const bam_hdr_t *h) __attribute__ ((__warn_unused_result__));

    int sam_parse1(kstring_t *s, bam_hdr_t *h, bam1_t *b) __attribute__ ((__warn_unused_result__));
    int sam_format1(const bam_hdr_t *h, const bam1_t *b, kstring_t *str) __attribute__ ((__warn_unused_result__));
    int sam_read1(samFile *fp, bam_hdr_t *h, bam1_t *b) __attribute__ ((__warn_unused_result__));
    int sam_write1(samFile *fp, const bam_hdr_t *h, const bam1_t *b) __attribute__ ((__warn_unused_result__));
# 398 "htslib/sam.h"
uint8_t *bam_aux_get(const bam1_t *b, const char tag[2]);







int64_t bam_aux2i(const uint8_t *s);







double bam_aux2f(const uint8_t *s);






char bam_aux2A(const uint8_t *s);






char *bam_aux2Z(const uint8_t *s);






uint32_t bam_auxB_len(const uint8_t *s);
# 445 "htslib/sam.h"
int64_t bam_auxB2i(const uint8_t *s, uint32_t idx);
# 456 "htslib/sam.h"
double bam_auxB2f(const uint8_t *s, uint32_t idx);
# 469 "htslib/sam.h"
int bam_aux_append(bam1_t *b, const char tag[2], char type, int len, const uint8_t *data);
# 478 "htslib/sam.h"
int bam_aux_del(bam1_t *b, uint8_t *s);
# 487 "htslib/sam.h"
int bam_aux_update_str(bam1_t *b, const char tag[2], int len, const char *data);
# 504 "htslib/sam.h"
typedef union {
    void *p;
    int64_t i;
    double f;
} bam_pileup_cd;
# 528 "htslib/sam.h"
typedef struct {
    bam1_t *b;
    int32_t qpos;
    int indel, level;
    uint32_t is_del:1, is_head:1, is_tail:1, is_refskip:1, aux:28;
    bam_pileup_cd cd;
} bam_pileup1_t;

typedef int (*bam_plp_auto_f)(void *data, bam1_t *b);

struct __bam_plp_t;
typedef struct __bam_plp_t *bam_plp_t;

struct __bam_mplp_t;
typedef struct __bam_mplp_t *bam_mplp_t;







    bam_plp_t bam_plp_init(bam_plp_auto_f func, void *data);
    void bam_plp_destroy(bam_plp_t iter);
    int bam_plp_push(bam_plp_t iter, const bam1_t *b);
    const bam_pileup1_t *bam_plp_next(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp);
    const bam_pileup1_t *bam_plp_auto(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp);
    void bam_plp_set_maxcnt(bam_plp_t iter, int maxcnt);
    void bam_plp_reset(bam_plp_t iter);
# 566 "htslib/sam.h"
    void bam_plp_constructor(bam_plp_t plp,
                             int (*func)(void *data, const bam1_t *b, bam_pileup_cd *cd));
    void bam_plp_destructor(bam_plp_t plp,
                            int (*func)(void *data, const bam1_t *b, bam_pileup_cd *cd));

    bam_mplp_t bam_mplp_init(int n, bam_plp_auto_f func, void **data);
# 580 "htslib/sam.h"
    void bam_mplp_init_overlaps(bam_mplp_t iter);
    void bam_mplp_destroy(bam_mplp_t iter);
    void bam_mplp_set_maxcnt(bam_mplp_t iter, int maxcnt);
    int bam_mplp_auto(bam_mplp_t iter, int *_tid, int *_pos, int *n_plp, const bam_pileup1_t **plp);
    void bam_mplp_reset(bam_mplp_t iter);
    void bam_mplp_constructor(bam_mplp_t iter,
                              int (*func)(void *data, const bam1_t *b, bam_pileup_cd *cd));
    void bam_mplp_destructor(bam_mplp_t iter,
                             int (*func)(void *data, const bam1_t *b, bam_pileup_cd *cd));
# 597 "htslib/sam.h"
int sam_cap_mapq(bam1_t *b, const char *ref, int ref_len, int thres);
int sam_prob_realn(bam1_t *b, const char *ref, int ref_len, int flag);
