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


def test_var_type():
    v = VCF(VCF_PATH)
    variant = next(v)
    assert variant.var_type == "indel", variant.var_type
    # 1 10172   .       CCCTAA  C       92.0    PASS
    for variant in v:
        if variant.POS == 10478: break
    else:
        raise Exception
    assert variant.var_type == "snp", variant.var_type

def _get_samples(v):
    import gzip
    import numpy as np
    import sys

    def _get_gt(s):
        if not ":" in s:
            return 2
        s = s.split(":", 1)[0]
        if s in ("0/0", "0|0"):
            return 0
        if s in ("0/1", "0|1", "0/.", "1/.", "./1", "1/."):
            return 1
        if s in ("1/1", "1|1"):
            return 3
        return 2

    for i, line in enumerate(gzip.open(VCF_PATH), start=1):
        if line[0] == "#": continue
        toks = line.strip().split("\t")
        if not (toks[0] == v.CHROM and int(toks[1]) == v.POS): continue
        if toks[3] != v.REF: continue
        if toks[4] not in v.ALT: continue

        samples = toks[9:]
        return  np.array([_get_gt(s) for s in samples], np.int)
    else:
        raise Exception("not found")


def test_gt_types():
    v = VCF(VCF_PATH)
    for variant in v:

        gt_types = variant.gt_types
        o = _get_samples(variant)
        assert (gt_types == o).all(), (variant, variant.CHROM, variant.POS, zip(gt_types, o))

def test_iterate():

    for i, v in enumerate(VCF(VCF_PATH), start=1):
        pass
    assert i == 115, i

