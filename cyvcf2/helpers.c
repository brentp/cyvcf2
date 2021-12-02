#include <helpers.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))

int as_gts(int32_t *gts, int num_samples, int ploidy, int strict_gt, int HOM_ALT, int UNKNOWN) {
    int j = 0, i, k;
    int missing= 0, found=0;
    for (i = 0; i < ploidy * num_samples; i += ploidy){
        missing = 0;
    found = 0;
        for (k = 0; k < ploidy; k++) {
            if bcf_gt_is_missing(gts[i+k])  {
                missing += 1;
            }
        }
        if (missing == ploidy) {
            gts[j++] = UNKNOWN; // unknown
            continue;
        } else if ( (missing != 0) && (strict_gt == 1) ) {
            gts[j++] = UNKNOWN; // unknown
            continue;
        }

        if(ploidy == 1 || gts[i+1] == bcf_int32_vector_end) {
            int a = bcf_gt_allele(gts[i]);
            if (a == 0) {
                   gts[j++] = 0;
            } else if (a == 1) {
                gts[j++] = HOM_ALT;
            } else {
                gts[j++] = UNKNOWN;
            }
            continue;
        }

        int a = bcf_gt_allele(gts[i]);
        int b = bcf_gt_allele(gts[i+1]);

        if((a == 0) && (b == 0)) {
            gts[j++] = 0; //  HOM_REF
            continue;
        }
        //fprintf(stderr, "i: %d\tmissing:%d\ta:%d\tb:%d\n", i/ploidy, missing, a, b);
        if ((missing > 0) && ((a == 0) || (b == 0))) {
            // if a single allele is missing e.g 0/. it's still encoded as hom ref because it has no alts
            gts[j++] = 0; // HOM_REF
            continue;
        }
        else if((a == 1) && (b == 1)) {
            gts[j] = HOM_ALT; //  HOM_ALT
        }
        else if((a != b)) {
            gts[j] = 1; //  HET
        }
        else if((a == b)) {
            gts[j] = HOM_ALT; //  HOM_ALT
        } else {
            gts[j] = UNKNOWN; // unknown
        }
        j++;
    }
    return j;
}

KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)
typedef khash_t(vdict) vdict_t;

// this is taken directly from atks/vt
int32_t* bcf_hdr_seqlen(const bcf_hdr_t *hdr, int32_t *nseq)
{
    vdict_t *d = (vdict_t*)hdr->dict[BCF_DT_CTG];
    int tid, m = kh_size(d);
    int32_t *lens = (int32_t*) malloc(m*sizeof(int32_t));
    khint_t k;
    int found = 0;

    for (k=kh_begin(d); k<kh_end(d); k++)
    {
        if ( !kh_exist(d,k) ) continue;
        tid = kh_val(d,k).id;
        lens[tid] = bcf_hrec_find_key(kh_val(d, k).hrec[0],"length");
        int j;
        if (lens[tid] > 0 && sscanf(kh_val(d, k).hrec[0]->vals[lens[tid]],"%d",&j) )
            lens[tid] = j;
    if(lens[tid] > 0){
      found++;
    }
    }
    *nseq = m;
    // found is used to check that we actually got the lengths.
    if(found == 0){
      *nseq = -1;
    }
    return lens;
}


