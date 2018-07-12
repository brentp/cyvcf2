#include "htslib/hts.h"
#include "htslib/vcf.h"

int main() {
	htsFile *hts = hts_open("seg.vcf.gz", "r");
	bcf_hdr_t *hdr = bcf_hdr_read(hts);

    htsFile *ohts = hts_open("-", "w");
	bcf_hdr_t *ohdr = bcf_hdr_dup(hdr);

	bcf_hdr_write(ohts, ohdr);

	bcf1_t *b = bcf_init();
    if(bcf_read(hts, hdr, b) < 0) {
		exit(1);
	}

	bcf_write(ohts, ohdr, b);
}
