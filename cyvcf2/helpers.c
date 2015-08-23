#include <htslib/vcf.h>

int as_gts(int *gts, int num_samples) {
    int j = 0, i;
    for (i = 0; i < 2 * num_samples; i += 2){
        if (bcf_gt_is_missing(gts[i]) && bcf_gt_is_missing(gts[i+1])){
            gts[j++] = 2; // unknown
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

int as_gts012(int *gts, int num_samples) {
    int j = 0, i;
    for (i = 0; i < 2 * num_samples; i += 2){
        if (bcf_gt_is_missing(gts[i]) && bcf_gt_is_missing(gts[i+1])){
            gts[j++] = 3; // unknown
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
