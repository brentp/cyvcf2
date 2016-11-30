from tabulate import tabulate
from collections import OrderedDict
vals = """
#bcftools
986798

real    3m49.024s
user    3m48.620s
sys 0m7.848s
#cyvcf2
984214

real    4m5.190s
user    4m5.048s
sys 0m0.112s
#pysam
984214

real    33m12.417s
user    33m11.920s
sys 0m0.244s
""".split("#")
vals = [x.strip() for x in vals if x.strip()]

key = 'time (seconds)'

def parse_group(g):
    lines = g.split("\n")
    d = OrderedDict([
            ('name', lines[0].strip()),
            ('time', next(x for x in lines if x.startswith('user')).split()[1].rstrip('s'))
            ])
    minutes = float(d['time'].split('m')[0])
    seconds = float(d['time'].split('m')[1])
    d[key] = 60 * minutes + seconds
    del d['time']
    return d



tbl = [parse_group(g) for g in vals]
base = next(d for d in tbl if d['name'] == 'cyvcf2')[key]
for d in tbl:
    d['ratio'] = "%.2f" % (d[key] / base)
    d[key] = "%.1f" % d[key]

print(r"\begin{table}[h]")
print(r"\caption{Timing VCF filtering}")
print(tabulate(tbl, headers="keys", tablefmt="latex"))
print(r"\end{table}")
