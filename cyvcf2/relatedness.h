#include <stdint.h>
#include <math.h>

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
// The result value should be ~1 for self and idential twins, ~0.5 for sibs and parent off-spring
// though that usually seems to be ~0.4 in practice.
// This should be called on few hundred to a few thousand variants that are
// not in linkage and have an aaf > 1 / n_samples (or so).
// asum and N are of length n_samples * n_samples and assumed to be in C order.
int related(int *gt_types, double *asum, int32_t *N, int32_t *ibs0, int32_t n_samples) {

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
		valj = gt_types[j];
		n_used++;
		for(k=j; k<n_samples; k++){
			if(gt_types[k] == 3){
				continue;
			}
			valk = gt_types[k];
			if(j != k){
				// multiply by 2 here to get the correct scale. differs from
				// original paper.
				numer = 2.0 * (valj - 2.0 * pi) * (valk - 2.0 * pi);
				ibs0[j * n_samples + k] += (valj != 1 && valk != 1 && valj != valk);
			} else {
				numer = (valj * valj) - (1.0 + 2.0 * pi) * valj + 2.0 * pi * pi;
				asum[j * n_samples + k]+= 1;
			}
			val = numer / denom;
			// heuristic to avoid too-large values
			if((j != k) && ((val > 4.5))) {
				val = 4.5;
			} else if (val < -4.5){
				continue;
			}

			asum[j * n_samples + k] += val;
			N[j * n_samples + k]+= 1;
		}
	}
	return n_used;
}

float r_unphased(int *a_gts, int *b_gts, float f, int32_t n_samples) {
	// http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2710162/pdf/GEN1823839.pdf
	// https://github.com/alanrogers/covld/blob/master/estimate_ld.c
    float vA, vB, vAB, cov, nsqr;
    int suma = 0, sumb = 0, sumaa = 0, sumbb = 0, sumab = 0;
    int i, n=0, a, b;

    for(i=0; i<n_samples; i++) {
        a = a_gts[i];
		if (a == 3) continue;
        b = b_gts[i];
		if (b == 3) continue;

        n += 1;
        suma += a;
        sumb += b;
        sumaa += a*a;
        sumbb += b*b;
        sumab += a*b;
    }

    nsqr = (double) n*(n-1);
    cov = (n*sumab - suma*sumb)/nsqr;
    vA = (n*sumaa - suma*suma)/nsqr;
    vB = (n*sumbb - sumb*sumb)/nsqr;

    vAB = vA*vB;
    if(vAB > 0) {
        return cov/sqrt(vAB);
    }
	return 0.0;
}
