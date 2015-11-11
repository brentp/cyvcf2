#include <htslib/vcf.h>

int as_gts(int *gts, int num_samples, int ploidy) {
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
			gts[j++] = 2; // unknown
			continue;
		} else if (missing != 0) {
			gts[j++] = 1; // HET
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
            gts[j] = 0; //  HOM_REF
        }
        else if((a == 1) && (b == 1)) {
            gts[j] = 3; //  HOM_ALT
        }
        else if((a  != b)) {
            gts[j] = 1; //  HET
        } else {
            gts[j] = 2; // unknown
        }
        j++;
    }
    // free the latter part of the array that we don't need.
    //free((void *)(&(gts[j])));
    return j;
}

int as_gts012(int *gts, int num_samples, int ploidy) {
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
		} else if (missing != 0) {
			gts[j++] = 1; // HET
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
            gts[j] = 0; //  HOM_REF
        }
        else if((a == 1) && (b == 1)) {
            gts[j] = 2; //  HOM_ALT
        }
        else if((a  != b)) {
            gts[j] = 1; //  HET
        } else {
            gts[j] = 3; // unknown
        }
        j++;
    }
    // free the latter part of the array that we don't need.
    //free((void *)(&(gts[j])));
    return j;
}
