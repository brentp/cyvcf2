import pytest

from ..cyvcf2 import VCF


@pytest.mark.parametrize("gts012, HOM_REF, HOM_ALT, UNKNOWN", [
    (False, 0, 3, 2),
    (True, 0, 2, 3),
])
def test_haploid_second_alt_is_hom_alt(tmp_path, gts012, HOM_REF, HOM_ALT, UNKNOWN):
    # A haploid call of the second ALT allele ("2") is a called genotype, as
    # the diploid "2/2" is, so it should not be counted as UNKNOWN (#204).
    path = tmp_path / "haploid-multiallelic.vcf"
    path.write_text(
        '##fileformat=VCFv4.2\n'
        '##contig=<ID=chrX>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\n'
        'chrX\t1\t.\tA\tC,G\t.\tPASS\t.\tGT\t0\t1\t2\t.\t2/2\n'
    )
    with VCF(path, gts012=gts012) as vcf:
        v = next(vcf)
        assert v.gt_types.tolist() == [HOM_REF, HOM_ALT, HOM_ALT, UNKNOWN, HOM_ALT]
        assert v.num_called == 4
        assert v.num_unknown == 1
        assert v.num_hom_alt == 3
