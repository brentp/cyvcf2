from cyvcf2 import VCF
import os.path
from nose.tools import assert_raises

HERE = os.path.dirname(__file__)
VCF_PATH = os.path.join(HERE, "test.vcf.gz")

def test_init():
    v = VCF(VCF_PATH)
    assert v

def test_bad_init():

    assert_raises(Exception, VCF, "XXXXX")
