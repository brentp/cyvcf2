cyvcf2
======

Note: cyvcf2 versions < 0.20.0 require htslib < 1.10. cyvcf2 versions >= 0.20.0 require htslib >= 1.10

<!-- ghp-import -p docs/build/html/ -->
The latest documentation for cyvcf2 can be found here:

[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](http://brentp.github.io/cyvcf2/)

If you use cyvcf2, please cite the [paper](https://academic.oup.com/bioinformatics/article/2971439/cyvcf2)


Fast Python **(3.8+)** parsing of VCF and BCF including region-queries.


[![Build](https://github.com/brentp/cyvcf2/actions/workflows/build.yml/badge.svg)](https://github.com/brentp/cyvcf2/actions/workflows/build.yml)

cyvcf2 is a cython wrapper around [htslib](https://github.com/samtools/htslib) built for fast parsing of [Variant Call Format](https://en.m.wikipedia.org/wiki/Variant_Call_Format) (VCF) files.

Attributes like `variant.gt_ref_depths` work for diploid samples and return a numpy array directly so they are immediately ready for downstream use.
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

# to query "all records" via __call__:
# this uses the index (HTS_IDX_START), so an index is required.
# if no index is available, this yields zero records.
all_vars = list(vcf())
```

Installation
============

## pip or uv with bundled htslib

If a binary wheel is available for your platform, pip and uv will install it
without building from source.

```
pip install cyvcf2
uv pip install cyvcf2
```

## pip or uv with system htslib

Assuming you have already built and installed htslib version 1.12 or higher.
```
CYVCF2_HTSLIB_MODE=EXTERNAL pip install --no-binary cyvcf2 cyvcf2
CYVCF2_HTSLIB_MODE=EXTERNAL uv pip install --no-binary cyvcf2 cyvcf2
```

## windows (experimental, only tested on MSYS2)

Assuming you have already built and installed htslib.
```
CYVCF2_HTSLIB_MODE=EXTERNAL pip install cyvcf2
```

## github (building htslib and cyvcf2 from source)

```
git clone --recursive https://github.com/brentp/cyvcf2
cd cyvcf2
CYVCF2_HTSLIB_MODE=BUILTIN python -m pip install .
CYVCF2_HTSLIB_MODE=BUILTIN uv pip install .
# or to use a system htslib.so
CYVCF2_HTSLIB_MODE=EXTERNAL python -m pip install .
CYVCF2_HTSLIB_MODE=EXTERNAL uv pip install .
```

Source builds use scikit-build-core and CMake. Python build dependencies such
as Cython, NumPy, CMake, and Ninja are installed by PEP 517 build isolation, but
you still need a C compiler plus the htslib native dependencies for your
platform.

Most users can leave `CYVCF2_CYTHONIZE` unset. It is mainly for developers
building from a Git checkout who need to control when CMake regenerates
`cyvcf2.c` from `cyvcf2/cyvcf2.pyx`:

- `AUTO` (the default) uses an existing `cyvcf2/cyvcf2.c` when present and runs
  Cython only when the file is missing.
- `ON` regenerates the C file after edits to `cyvcf2.pyx` or `cyvcf2.pxd`, or
  when testing Cython-generated output.
- `OFF` requires an existing C file.

For example:

```
CYVCF2_CYTHONIZE=ON python -m pip install .
python -m pip install . --config-settings=cmake.define.CYVCF2_CYTHONIZE=ON
```

The legacy `CYTHONIZE=1` environment variable is still accepted.

To build source distributions and wheels locally:

```
python -m build
uv build
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

Install `pytest`, then tests can be run with:

```
pytest
uv run --with pytest pytest
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

Pysam also [has a cython wrapper to htslib](https://github.com/pysam-developers/pysam/blob/master/pysam/libcbcf.pyx) and one block of code here is taken directly from that library. But, the optimizations that we want for gemini are very specific so we have chosen to create a separate project.

Performance
===========

For the performance comparison in the paper, we used [thousand genomes chromosome 22](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz)
With the full comparison runner [here](https://github.com/brentp/cyvcf2/blob/main/scripts/compare.sh).
