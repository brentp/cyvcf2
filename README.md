cyvcf2
======

Fast python **(2 and 3)** parsing of VCF and BCF including region-queries.

[![Build Status](https://travis-ci.org/brentp/cyvcf2.svg?branch=master)](https://travis-ci.org/brentp/cyvcf2)

cyvcf2 is a cython wrapper around [htslib](https://github.com/samtools/htslib) built for fast parsing of [Variant Call Format](https://en.m.wikipedia.org/wiki/Variant_Call_Format) (VCF) files.

Attributes like `variant.gt_ref_depths` return a numpy array directly so they are immediately ready for downstream use.
**note** that the array is backed by the underlying C data, so, once `variant` goes out of scope. The array will contain nonsense.
To persist a copy, use: `cpy = np.array(variant.gt_ref_depths)` instead of just `arr = variant.gt_ref_depths`.

Example
=======

The example below shows much of the use of cyvcf2.

```Python
from cyvcf2 import VCF

for variant in VCF('some.vcf.gz'): # or VCF('some.bcf')
	variant.REF, variant.ALT # e.g. REF='A', ALT=['C', 'T']


	variant.CHROM, variant.start, variant.end, variant.ID, \
				variant.FILTER, variant.QUAL

	# numpy arrays of specific things we pull from the sample fields.
	# gt_types is array of 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
	variant.gt_types, variant.gt_ref_depths, variant.gt_alt_depths # numpy arrays
	variant.gt_phases, variant.gt_quals, variant.gt_bases # numpy array


	## INFO Field.
	## extract from the info field by it's name:
	variant.INFO.get('DP') # int
	variant.INFO.get('FS') # float
	variant.INFO.get('AC') # float

	# convert back to a string.
	str(variant)


	## sample info...

	# Get a numpy array of the depth per sample:
    dp = variant.format('DP')
    # or of any other format field:
    sb = variant.format('SB')
    assert sb.shape == (n_samples, 4) # 4-values per

# to do a region-query:

vcf = VCF('some.vcf.gz')
for v in vcf('11:435345-556565'):
    if v.INFO["AF"] > 0.1: continue
    print(str(v))

```

Installation
============

```
pip install cyvcf2
```

Testing
=======

Tests can be run with:

```
python setup.py test
```

See Also
========

Pysam also [has a cython wrapper to htslib](https://github.com/pysam-developers/pysam/blob/master/pysam/cbcf.pyx) and one block of code here is taken directly from that library. But, the optimizations that we want for gemini are very specific so we have chosen to create a separate project.
