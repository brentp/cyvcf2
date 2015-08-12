from cyvcf2 import VCF, Variant
import os.path
from nose.tools import assert_raises

HERE = os.path.dirname(__file__)
VCF_PATH = os.path.join(HERE, "test.vcf.gz")

def test_init():
    v = VCF(VCF_PATH)
    assert v

def test_bad_init():
    assert_raises(Exception, VCF, "XXXXX")

def test_samples():
    v = VCF(VCF_PATH)
    assert len(v.samples) == 189

def test_next():
    v = VCF(VCF_PATH)
    variant = next(v)
    assert isinstance(variant, Variant)

def test_attrs():
    # 1 10172   .       CCCTAA  C       92.0    PASS
    v = VCF(VCF_PATH)
    variant = next(v)
    assert variant.POS == 10172
    assert variant.CHROM == '1'
    assert variant.ID is None, variant.ID
    assert variant.start == 10171
    assert variant.end == 10177, variant.end
    assert variant.FILTER is None
    assert variant.QUAL == 92.0

    assert variant.REF == "CCCTAA"
    assert variant.ALT == ["C"]

def test_vartype():
    v = VCF(VCF_PATH)
    variant = next(v)
    assert variant.var_type == "indel", variant.var_type

def test_iterate():

    for i, v in enumerate(VCF(VCF_PATH), start=1):
        pass
    assert i == 115, i

