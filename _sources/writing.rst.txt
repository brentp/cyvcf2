Modifying Existing VCFs
=======================

`cyvcf2` is optimized for fast reading and extraction from existing files.
However, it also offers some means of modifying existing VCFs.

Modifying the INFO field
------------------------

.. this is the same example as on the index.rst page

Here, we show an example of how to modify the INFO field at a locus to annotate
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

Modifying genotypes and the FORMAT field
----------------------------------------

Here is an example where we filter some calls from records in a VCF. This demonstrates
how to add to the FORMAT field and how to modify genotypes at a locus.

.. code-block:: python

    from cyvcf2 import VCF, Writer
    vcf = VCF(VCF_PATH)
    # adjust the header to contain the new field
    # the keys 'ID', 'Description', 'Type', and 'Number' are required.
    vcf.add_format_to_header({
        'ID': 'FILTER_CODE',
        'Description': 'Numeric code for filtering reason',
        'Type': 'Integer',
        'Number': '1'
    })

    # create a new vcf Writer using the input vcf as a template.
    fname = "out.vcf"
    w = Writer(fname, vcf)

    for v in vcf:
        # The filter_samples function is not shown.
        # This could be any manipulation required by the user.
        # Since we specified the FILTER_CODE format field is an Integer
        # with Number=1, reasons must be an n x 1 numpy array of integers
        # where n is the number of samples.
        indicies, reasons = filter_samples(v)
        if indicies:
            # add the reasons array to the format dictionary at this locus
            v.set_format('FILTER_CODE', reasons)
            for index in indicies:
                # overwrite the genotypes of each filtered locus to be nocalls
                # Note: until the reassignment to v.genotypes below, this
                # leaves the v Variant object in an inconsistent state
                v.genotypes[index] = [-1]*v.ploidy + [False]
            # it is necessary to reassign the genotypes field
            # so that the v Variant object reprocess it and future calls
            # to functions like v.genotype.array() properly reflect
            # the changed genotypes
            v.genotypes = v.genotypes

        w.write_record(v)

    w.close(); vcf.close()

.. _Limitations with writing:

Known Limitations with Writing VCFs
-----------------------------------
* `cyvcf2` currently does not support writing VCFs encoded with UTF-8 with non-ASCII characters in the contents of string-typed FORMAT fields.
* `cyvcf2` currently does not support writing string type format fields with number>1.

