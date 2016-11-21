import sys
from pysam import VariantFile

vcf = VariantFile(sys.argv[1])
n = 0

for v in vcf:
    if len(v.alts) > 1: continue
    if v.qual < 20: continue
    gts = [s['GT'] for s in v.samples.values()]
    an = sum(len(gt) for gt in gts)
    ac = sum(sum(gt) for gt in gts)
    aaf = (float(ac) / float(an))
    if aaf > 0.05: continue
    n += 1
print(n)
