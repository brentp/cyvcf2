cyvcf2
======

Fast python **(2 and 3)** parsing of VCF and BCF including region-queries.

[![Build Status](https://travis-ci.org/brentp/cyvcf2.svg?branch=master)](https://travis-ci.org/brentp/cyvcf2)

cyvcf2 is a cython wrapper around [htslib](https://github.com/samtools/htslib) built for fast parsing of [Variant Call Format](https://en.m.wikipedia.org/wiki/Variant_Call_Format) (VCF) files.
It is targetted toward our use-case in [gemini](http://gemini.rtfd.org) but should also be of general utility.

On a file with 189 samples that takes [cyvcf](https://github.com/arq5x/cyvcf) **21 seconds** to parse and extract all sample information, it takes `cyvcf2` **1.4 seconds**.

Attributes like `variant.gt_ref_depths` return a numpy array directly so they are immediately ready for downstream use.
**note** that the array is backed by the underlying C data, so, once `variant` goes out of scope. The array will contain nonsense.
To persist a copy, use: `cpy = np.array(variant.gt_ref_depths)` instead of just `arr = variant.gt_ref_depths`.

Example
=======

The example below shows much of the use of cyvcf2.

```Python
from cyvcf2 import VCF

for variant in VCF('some.vcf.gz'): # or VCF('some.bcf')

	variant.gt_types # numpy array
	variant.gt_ref_depths, variant.gt_alt_depths # numpy arrays
	variant.gt_phases, variant.gt_quals # numpy arrays
	variant.gt_bases # numpy array
	variant.CHROM, variant.start, variant.end, variant.ID, \
				variant.REF, variant.ALT, variant.FILTER, variant.QUAL
	variant.INFO.get('DP') # int
	variant.INFO.get('FS') # float
	variant.INFO.get('AC') # float
    a = variant.gt_phred_ll_homref # numpy array
    b = variant.gt_phred_ll_het # numpy array
    c = variant.gt_phred_ll_homalt # numpy array

	str(variant)
	# Get a numpy array of the depth per sample:
    dp = variant.format('DP', int)
    # or of any other format field:
    sb = variant.format('SB', float)
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
