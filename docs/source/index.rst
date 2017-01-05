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

        # single sample of 0|1 in vcf becomes [[0, 1, True]]
        # 2 samples of 0/0 and 1|1 would be [[0, 0, False], [1, 1, True]]
        print v.genotypes 

Modifying Existing Records
==========================

`cyvcf2` is optimized for fast reading and extraction from existing files.
However, it also offers some means of modifying existing VCFs. Here, wrapper
show an example of how to annotate variants with the genes that they overlap.


.. code-block:: python

    from cyvcf2 import VCF, Writer
    vcf = VCF(VCF_PATH)
    # adjust the header to contain the new field
    # the keys 'ID', 'Description', 'Type', and 'Number' are required.
    vcf.add_info_to_header({'ID': 'gene', 'Description': 'overlapping gene',
        'Type':'Character', 'Number': '1'})

    # create a new vcf Writer using the input vcf as a template.
    w = Writer(f, vcf)

    for v in vcf:
        # The get_gene_intersections function is not shown.
        # This could be any operation to find intersections
        # or any manipulation required by the user.
        genes = get_gene_intersections(v)
        if genes is not None:
            v.INFO["gene"] = ",".join(genes)
        w.write_record(v)

    w.close(); vcf.close()


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

