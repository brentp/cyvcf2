.. cyvcf2 documentation master file, created by
   sphinx-quickstart on Mon Nov 14 08:40:54 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

cyvcf2
======

cyvcf2 is a fast VCF parser for python (2 and 3) released under the MIT license.
It is a `cython` wrapper for `htslib <https://github.com/samtools/htslib/>`_ developed in 
the `Quinlan lab <http://quinlanlab.org/>`_.

See the :ref:`api` for detailed documentation, but the most common usage is summarized in the snippet below:

.. code-block:: python

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


        ## per-sample info...
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


Installation
============

.. code-block:: bash

    pip install cyvcf2

or via bioconda.

Testing
=======

Tests can be run with:

.. code-block:: bash

    python setup.py test

See Also
========

Pysam also [has a cython wrapper to htslib](https://github.com/pysam-developers/pysam/blob/master/pysam/cbcf.pyx) and one block of code here is taken directly from that library. But, the optimizations that we want for gemini are very specific so we have chosen to create a separate project.

.. toctree::
   :maxdepth: 2

   docstrings

