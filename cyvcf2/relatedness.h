#include <stdint.h>

// internal calculate of alternate allele frequency.
inline float aaf(int *gt_types, int32_t n_samples){
	float af = 0;
	int i, n_called = 0;

	for (i = 0; i < n_samples; i++){
		if(gt_types[i] == 3){
			continue;
		}
		af += gt_types[i];
		n_called += 1;
	}
	return af / (float)(2 * n_called);
}

// related takes an array of genotypes (0=HOM_REF, 1=HET, 2=HOMALT, 3=UNKNOWN) and updates asum and N
// which are used to calculate relatedness between samples j, k as asum[j, k] / N[j, k].
// This should be called on few hundred to a few thousand variants that are
// not in linkage and have an aaf > 1 / n_samples (or so).
// asum and N are of length n_samples * n_samples and assumed to be in C order.
int related(int *gt_types, double *asum, int32_t *N, int32_t n_samples) {

	int n_used = 0;
	int32_t j, k;
	float pi = aaf(gt_types, n_samples);
	float numer, val;
	float valj, valk;
	float denom = 2.0 * pi * (1.0 - pi);

	for(j=0; j <n_samples; j++){
		// skip unknown
		if(gt_types[j] == 3){
			continue;
		}
		valj = (float)gt_types[j];
		n_used++;
		for(k=j; k<n_samples; k++){
			if(gt_types[k] == 3){
				continue;
			}
			valk = (float)gt_types[k];
			if(j != k){
				// multiply by 2 here to get the correct scale. differs from
				// original paper.
				numer = 2.0 * (valj - 2.0 * pi) * (valk - 2.0 * pi);
			}else {
				numer = (valj * valj) - (1.0 + 2.0 * pi) * valj + 2.0 * pi * pi;
				asum[j * n_samples + k]+= 1;
			}
			val = numer / denom;
			// heuristic to avoid too-large values
			if((j != k) && (val > 2.9)) {
				val = 2.9;
			} else if (val < -2.5){
				continue;
			}
			asum[j * n_samples + k] += val;
			N[j * n_samples + k]+= 1;
		}
	}
	return n_used;
}
