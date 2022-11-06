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

.. this example is also in the writing.rst page

`cyvcf2` is optimized for fast reading and extraction from existing files.
However, it also offers some means of modifying existing VCFs. Here, we
show an example of how to modify the INFO field at a locus to annotate
variants with the genes that they overlap.

.. code-block:: python

    from cyvcf2 import VCF, Writer
    vcf = VCF(VCF_PATH)
    # adjust the header to contain the new field
    # the keys 'ID', 'Description', 'Type', and 'Number' are required.
    vcf.add_info_to_header({'ID': 'gene', 'Description': 'overlapping gene',
        'Type':'Character', 'Number': '1'})

    # create a new vcf Writer using the input vcf as a template.
    fname = "out.vcf"
    w = Writer(fname, vcf)

    for v in vcf:
        # The get_gene_intersections function is not shown.
        # This could be any operation to find intersections
        # or any manipulation required by the user.
        genes = get_gene_intersections(v)
        if genes is not None:
            v.INFO["gene"] = ",".join(genes)
        w.write_record(v)

    w.close(); vcf.close()

More info on writing vcfs can be found :doc:`here <writing>`

Setting Genotyping Strictness
=============================

By default, genotypes containing a single missing allele (e.g. genotypes such
as ``0/.``, ``./0``, ``1/.``, or ``./1``) are classified as heterozygous ("HET"
or id: 1) instead of "UNKNOWN" (or id: 2).  A case can be made for either
classification on these partially missing genotypes depending on the situation.

A better way to explicitly control the genotype classification for these cases,
is to use the "strict genotype", or ``strict_gt``, feature.  Here's an example
usage:

.. code-block:: python

    from cyvcf2 import VCF
    vcf = VCF("/path/to/vcf/file", strict_gt=True)
    for variant in vcf:
        # do something

When the ``strict_gt`` flag is enabled, `cyvcf2` will treat any genotype
containing a missing allele (containing a '.') as an UNKNOWN genotype;
otherwise, genotypes like ``0/.``, ``./0``, ``1/.``, or ``./1`` will be
classified as heterozygous ("HET").

The "strict genotype" feature is not enabled by default (i.e. ``strict_gt`` is
set to `False` by default) to preserve backwards compatibility.

Installation
============

.. code-block:: bash

    pip install cyvcf2

or via bioconda.

Testing
=======

Install `pytest`, then tests can be run with:

.. code-block:: bash

    pytest

Known Limitations
=================
* `cyvcf2` currently does not support reading VCFs encoded with UTF-8 with non-ASCII characters in the contents of string-typed FORMAT fields.

For limitations on writing VCFs, see :ref:`here <Limitations with writing>`

See Also
========

Pysam also `has a cython wrapper to htslib <https://github.com/pysam-developers/pysam/blob/master/pysam/cbcf.pyx>`_ and one block of code here is taken directly from that library. But, the optimizations that we want for gemini are very specific so we have chosen to create a separate project.

.. toctree::
   :maxdepth: 1

   docstrings
   writing

