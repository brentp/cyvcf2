from cyvcf2 import VCF, Variant
import os.path
from nose.tools import assert_raises

HERE = os.path.dirname(__file__)
VCF_PATH = os.path.join(HERE, "test.vcf.gz")
VCF_PATH2 = os.path.join(HERE, "test.snpeff.vcf")
VCF_PHASE_PATH = os.path.join(HERE, "test.comp_het.3.vcf")

def test_init():
    v = VCF(VCF_PATH)
    assert v

def test_type():
    vcf = VCF(VCF_PATH)
    for v in vcf:
        if len(v.REF) == 1 and len(v.ALT[0]) == 1:
            assert v.var_type == 'snp'
        elif v.ALT[0][0] != "<":
            assert v.var_type == 'indel'
        else:
            print v.var_type, v.REF, v.ALT

def test_pls():
    vcf = VCF(VCF_PATH)
    v = next(vcf)

    assert v.gt_phred_ll_homref[0] == 0, v.gt_phred_ll_homref[0]
    assert v.gt_phred_ll_het[0] == 7, v.gt_phred_ll_het[0]
    assert v.gt_phred_ll_homalt[0] == 922, v.gt_phred_ll_homalt[0]

    import numpy as np
    imax = np.iinfo(np.int32(0)).max
    # missing
    assert v.gt_phred_ll_homref[1] == imax, v.gt_phred_ll_homref[1]
    assert v.gt_phred_ll_het[1] == imax, v.gt_phred_ll_het[1]
    assert v.gt_phred_ll_homalt[1] == imax, v.gt_phred_ll_homalt[1]

def test_str():
    vcf = VCF(VCF_PATH, lazy=True)
    v = next(vcf)
    s = str(v)

    assert "10172\t.\tCCCTAA\t" in s
    assert "fitcons_float=0.1266" in s


def test_region():
    vcf = VCF(VCF_PATH)

    start = 12783
    end = 13783
    k = 0
    reg = '1:%d-%d' % (start, end)
    for var in vcf(reg):
        k += 1
        assert var.start <= end, var
        assert var.end >= start, var
        assert isinstance(var.REF, basestring)
        assert isinstance(var.ALT, list)
    assert k == 28, k

def test_phases():
    vcf = VCF(VCF_PHASE_PATH)

    v = next(vcf)
    assert all(v.gt_phases), v.gt_phases

    v = next(vcf)
    assert all(v.gt_phases)

    v = next(vcf)
    assert all(v.gt_phases[1::2])
    assert not any(v.gt_phases[0::2])

    v = next(vcf)
    assert not any(v.gt_phases[:-1])
    assert v.gt_phases[-1]

    v = next(vcf)
    assert not any(v.gt_phases)

def test_bad_init():
    assert_raises(Exception, VCF, "XXXXX")

def test_samples():
    v = VCF(VCF_PATH)
    assert len(v.samples) == 189

def test_next():
    v = VCF(VCF_PATH)
    variant = next(v)
    assert isinstance(variant, Variant)

def test_info_dict():
    v = VCF(VCF_PATH)
    variant = next(v)
    d = dict(variant.INFO)
    assert d != {}, d
    toks = _get_line_for(variant)

    info = toks[7].split(";")
    keys = [x.split('=')[0] for x in info]
    for k in keys:
        assert k in d, (k, info)


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

def _get_line_for(v):
    import gzip

    for i, line in enumerate(gzip.open(VCF_PATH), start=1):
        if line[0] == "#": continue
        toks = line.strip().split("\t")
        if not (toks[0] == v.CHROM and int(toks[1]) == v.POS): continue
        if toks[3] != v.REF: continue
        if toks[4] not in v.ALT: continue
        return toks
    else:
        raise Exception("not found")


def _get_samples(v):
    import numpy as np

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
    toks = _get_line_for(v)
    samples = toks[9:]
    return np.array([_get_gt(s) for s in samples], np.int)

def test_header_info():
    v = VCF(VCF_PATH)
    csq = v['CSQ']
    assert csq['ID'] == "CSQ"
    assert "Description" in csq


    assert_raises(KeyError, v.__getitem__, 'XXXXX')

def test_snpeff_header():
    v = VCF(VCF_PATH2)

    f = v['SnpEffVersion']
    assert f != {}, f
    assert 'SnpEffVersion' in f

"""
def test_info_update():
    vcf = VCF(VCF_PATH)
    v = next(vcf)
    ret = v.INFO.update({'k': 22})
    #print ret
    #assert v.INFO['k'] == 22
    #assert ret == 0, ret
"""

def test_gt_types():
    v = VCF(VCF_PATH)
    for variant in v:

        gt_types = variant.gt_types
        o = _get_samples(variant)
        assert (gt_types == o).all(), (variant, variant.CHROM, variant.POS, zip(gt_types, o))

def test_raw_header():
    v = VCF(VCF_PATH)
    h = v.raw_header.strip().split("\n")
    s = h[0]
    assert s == "##fileformat=VCFv4.1", s
    assert len(h) == 185, len(h)



def test_iterate():

    for i, v in enumerate(VCF(VCF_PATH), start=1):
        pass
    assert i == 115, i

