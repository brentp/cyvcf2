#include <stdint.h>

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

int related(int *gt_types, double *asum, int32_t *N, int32_t n_samples) {

	int n_used = 0, nn;
	int32_t j, k;
	float pi = aaf(gt_types, n_samples);
	float numer, val;
	float valj, valk;
	float denom = 2.0 * pi * (1.0 - pi);

	for(j=0; j <n_samples; j++){
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
				numer = 2.0 * (valj - 2.0 * pi) * (valk - 2.0 * pi);
			}else {
				numer = (valj * valj) - (1.0 + 2.0 * pi) * valj + 2.0 * pi * pi;
				asum[j * n_samples + k]+= 1;
			}
			val = numer / denom;
			if((j != k) && (val > 2.5)) {
				val = 2.5;
			} else if (val < -2.5){
				continue;
			}
			asum[j * n_samples + k] += val;
			N[j * n_samples + k]+= 1;
			//#N[j * n_samples + k] = nn + 1;
			//fprintf(stderr, "after:%d\n", N[j * n_samples + k]);
			//fprintf(stderr, "asum:%.1f\n", asum[j * n_samples + k]);
		}
	}
	return n_used;
}
