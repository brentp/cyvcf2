.. _api:

cyvcf2 API
==========

.. _genotype-phasing:

Genotype phasing
----------------

``Variant.gt_phases``, the final item in each ``Variant.genotypes`` sample,
``Variant.genotype.phased(i)``, and the final column of
``Variant.genotype.array()`` summarize phasing with one flag per sample.
For samples with more than one actual allele, this flag uses the phase bit
of the second allele. For diploid calls this corresponds to the
single separator: ``0/1`` is unphased and ``0|1`` is phased.

This convention is retained for backwards compatibility, but does not
describe all phasing in polyploid calls. For example, ``0|1/2`` reports
``True`` while ``0/1|2`` reports ``False``; neither is fully phased.
Do not interpret the flag as meaning that every allele is phased.
Here, polyploid means more than two allele copies in a sample's genotype,
not more than two possible alleles at a variant site.

For haploid samples, the flag uses the only allele's phase bit, including
when the sample is padded to match other samples in a mixed-ploidy record.
Padding slots are not additional alleles and their bits are not phasing
information. htslib 1.22 and later treat plain called haploid genotypes such
as ``1`` as implicitly phased; an unphased missing haploid call ``.`` is
not made phased by padding.

To inspect individual allele phase bits, use
``variant.genotype[sample_index][allele_index].phased`` for actual alleles
(not padding), or inspect the original GT separators in ``str(variant)``.

.. autoclass:: cyvcf2.cyvcf2.VCF
    :members:
    :no-undoc-members:
    :exclude-members: het_check,ibd,next,plot_relatedness,set_threads,site_relatedness,update
    
.. autoclass:: cyvcf2.cyvcf2.Variant
    :members:
    :exclude-members: relatedness

.. autoclass:: cyvcf2.cyvcf2.Writer
    :members:
    :no-undoc-members:

.. autoclass:: cyvcf2.cyvcf2.INFO
    :no-undoc-members:
