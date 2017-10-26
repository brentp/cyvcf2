#include <helpers.h>

int as_gts(int32_t *gts, int num_samples, int ploidy, int strict_gt) {
    int j = 0, i, k;
	int missing= 0;
    for (i = 0; i < ploidy * num_samples; i += ploidy){
		missing = 0;
		for (k = 0; k < ploidy; k++) {
			if (gts[i+k] <= 0) {
				missing += 1;
			}
			//fprintf(stderr, "%d\n", gts[i + k]);
        }
		if (missing == ploidy) {
			gts[j++] = 2; // unknown
			continue;
		} else if ( (missing != 0) && (strict_gt == 1) ) {
			gts[j++] = 2; // unknown
			continue;
		}

		if(ploidy == 1) {
			int a = bcf_gt_allele(gts[i]);
			if (a == 0) {
			   	gts[j++] = 0;
			} else if (a == 1) {
				gts[j++] = 3;
			} else {
				gts[j++] = 2;
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
            gts[j] = 3; //  HOM_ALT
        }
        else if((a  != b)) {
            gts[j] = 1; //  HET
        }
        else if((a == b)) {
            gts[j] = 3; //  HOM_ALT
        } else {
            gts[j] = 2; // unknown
        }
        j++;
    }
    return j;
}

int as_gts012(int32_t *gts, int num_samples, int ploidy, int strict_gt) {
    int j = 0, i, k;
	int missing;
    for (i = 0; i < ploidy * num_samples; i += ploidy){
		missing = 0;
		for (k = 0; k < ploidy; k++) {
			if (gts[i+k] <= 0) {
				missing += 1;
			}
			//fprintf(stderr, "%d\n", gts[i + k]);
        }
		if (missing == ploidy) {
			gts[j++] = 3; // unknown
			continue;
		} else if ( (missing != 0) && (strict_gt == 1) ) {
			gts[j++] = 3; // unknown
			continue;
		}

		if(ploidy == 1) {
			int a = bcf_gt_allele(gts[i]);
			if (a == 0) {
			   	gts[j++] = 0;
			} else if (a == 1) {
				gts[j++] = 2;
			} else {
				gts[j++] = 3;
			}
			continue;
		}

        int a = bcf_gt_allele(gts[i]);
        int b = bcf_gt_allele(gts[i+1]);

        if((a == 0) && (b == 0)) {
            gts[j++] = 0; //  HOM_REF
            continue;
        }
        if ((missing > 0) && ((a == 0)  || b == 0)) {
            gts[j] = 0; // HOM_REF
        }
        else if((a == 1) && (b == 1)) {
            gts[j] = 2; //  HOM_ALT
        }
        else if((a  != b)) {
            gts[j] = 1; //  HET
        }
        else if ((a == b)) {
            gts[j] = 2; //  HOM_ALT
        } else {
            gts[j] = 3; // unknown
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

    for (k=kh_begin(d); k<kh_end(d); k++)
    {
        if ( !kh_exist(d,k) ) continue;
        tid = kh_val(d,k).id;
        lens[tid] = bcf_hrec_find_key(kh_val(d, k).hrec[0],"length");
        int j;
        if ( sscanf(kh_val(d, k).hrec[0]->vals[lens[tid]],"%d",&j) )
            lens[tid] = j;
    }
	*nseq = m;
    return lens;
}


