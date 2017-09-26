import sys
from vcf import Reader
import gzip

vcf = Reader(open(sys.argv[1], 'rb'))
n = 0

for v in vcf:
    if len(v.ALT) > 1: continue
    if v.QUAL < 20: continue
    if v.aaf[0] > 0.05: continue
    n += 1
print(n)
