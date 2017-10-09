#include <htslib/vcf.h>

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
