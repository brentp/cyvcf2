import sys
import cyvcf2
n = 0
for v in cyvcf2.VCF(sys.argv[1]):
    if len(v.ALT) > 1: continue
    if v.QUAL < 20: continue
    if v.aaf > 0.05: continue
    n += 1
print(n)
