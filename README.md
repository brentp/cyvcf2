cyvcf2
======

<!-- ghp-import -p docs/build/html/ -->
[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](http://brentp.github.io/cyvcf2/)

If you use cyvcf2, please cite the [paper](https://academic.oup.com/bioinformatics/article/2971439/cyvcf2)


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

## pip
```
pip install cyvcf2
```

## github (building htslib and cyvcf2 from source)

```
git clone --recursive https://github.com/brentp/cyvcf2
cd cyvcf2/htslib
autoheader
autoconf
./configure --enable-libcurl
make

cd ..
pip install -e .
```

On **OSX**, using brew, you may have to set the following as indicated by the brew install:

```
For compilers to find openssl you may need to set:
  export LDFLAGS="-L/usr/local/opt/openssl/lib"
  export CPPFLAGS="-I/usr/local/opt/openssl/include"

For pkg-config to find openssl you may need to set:
  export PKG_CONFIG_PATH="/usr/local/opt/openssl/lib/pkgconfig"
```

Testing
=======

Tests can be run with:

```
python setup.py test
```

CLI
=======
Run with `cyvcf2 path_to_vcf`

```
$ cyvcf2 --help
Usage: cyvcf2 [OPTIONS] <vcf_file> or -

  fast vcf parsing with cython + htslib

Options:
  -c, --chrom TEXT                Specify what chromosome to include.
  -s, --start INTEGER             Specify the start of region.
  -e, --end INTEGER               Specify the end of the region.
  --include TEXT                  Specify what info field to include.
  --exclude TEXT                  Specify what info field to exclude.
  --loglevel [DEBUG|INFO|WARNING|ERROR|CRITICAL]
                                  Set the level of log output.  [default:
                                  INFO]
  --silent                        Skip printing of vcf.
  --help                          Show this message and exit.
```


See Also
========

Pysam also [has a cython wrapper to htslib](https://github.com/pysam-developers/pysam/blob/master/pysam/cbcf.pyx) and one block of code here is taken directly from that library. But, the optimizations that we want for gemini are very specific so we have chosen to create a separate project.

Performance
===========

For the performance comparison in the paper, we used [thousand genomes chromosome 22](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz)
With the full comparison runner [here](https://github.com/brentp/cyvcf2/blob/master/scripts/compare.sh).
