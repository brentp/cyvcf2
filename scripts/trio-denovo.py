import numpy as np
import sys
from scipy.stats import binom_test
path = sys.argv[1]

import cyvcf2
vcf = cyvcf2.VCF(path)
PRO, MOM, DAD = range(3)

for v in vcf:
    if v.QUAL < 10: continue
    if np.any(v.gt_depths < 10): continue
    refs, alts = v.gt_ref_depths, v.gt_alt_depths
    if not all(v.gt_types == [vcf.HET, vcf.HOM_REF, vcf.HOM_REF]): continue
    if alts[MOM] > 1 or alts[DAD] > 1: continue
    print("%s\t%d\t%s\t%s\t%s" % (v.CHROM, v.start, v.REF,
        ",".join(v.ALT), ",".join(v.gt_bases)))
