from cyvcf2 import VCF, Variant, Writer
import os.path
from nose.tools import assert_raises
import tempfile
import os
import atexit

HERE = os.path.dirname(__file__)
VCF_PATH = os.path.join(HERE, "test.vcf.gz")
VCF_PATH2 = os.path.join(HERE, "test.snpeff.vcf")
VCF_PHASE_PATH = os.path.join(HERE, "test.comp_het.3.vcf")

try:
    basestring
except NameError:
    basestring = (str, bytes)

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
            print(v.var_type, v.REF, v.ALT)

def test_ibd():
    samples = ['101976-101976', '100920-100920', '100231-100231']
    vcf = VCF(VCF_PATH, gts012=True, samples=samples)
    res = vcf.ibd()
    assert len(res) == 3, (len(res))
    arr = res[(b'101976-101976', b'100920-100920')]
    assert len(arr) > 0

def test_relatedness():
    vcf = VCF(VCF_PATH, gts012=True)
    viter = iter(vcf.relatedness(gap=0, linkage_max=2))
    res = next(viter)
    assert "ibs0" in res
    assert "ibs2" in res
    assert "ibs2*" in res
    #vcf = VCF(VCF_PATH, gts012=True)
    #for r in viter:
    #    print r['pair'], r['ibs0'], r['ibs2'], r['ibs2*']



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
    for var in vcf(reg.encode()):
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

    info = toks[7].split(b";")
    keys = [x.split(b'=')[0] for x in info]
    for k in keys:
        assert k in d, (k, info)


def test_attrs():
    # 1 10172   .       CCCTAA  C       92.0    PASS
    v = VCF(VCF_PATH)
    variant = next(v)
    assert variant.POS == 10172
    assert variant.CHROM == b"1"
    assert variant.ID is None, variant.ID
    assert variant.start == 10171
    assert variant.end == 10177, variant.end
    assert variant.FILTER is None
    assert variant.QUAL == 92.0

    assert variant.REF == b"CCCTAA"
    assert variant.ALT == [b"C"]

def test_writer():

    v = VCF(VCF_PATH)
    f = tempfile.mktemp(suffix=".vcf")
    atexit.register(os.unlink, f)
    o = Writer(f, v)
    rec = v.next()
    rec.INFO["AC"] = "3"
    rec.FILTER = ["LowQual"]
    o.write_record(rec)

    rec.FILTER = ["LowQual", "VQSRTrancheSNP99.90to100.00"]
    o.write_record(rec)


    rec.FILTER = "PASS"
    o.write_record(rec)

    o.close()

    expected = ["LowQual", "LowQual;VQSRTrancheSNP99.90to100.00", None]

    for i, variant in enumerate(VCF(f)):
        assert variant.FILTER == expected[i], (variant.FILTER, expected[i])




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
        if line[0] == b"#": continue
        toks = line.strip().split(b"\t")
        if not (toks[0] == v.CHROM and int(toks[1]) == v.POS): continue
        if toks[3] != v.REF: continue
        if toks[4] not in v.ALT: continue
        return toks
    else:
        raise Exception("not found")


def _get_samples(v):
    import numpy as np

    def _get_gt(s):
        if not b":" in s:
            return 2
        s = s.split(b":", 1)[0]
        if s in (b"0/0", b"0|0"):
            return 0
        if s in (b"0/1", b"0|1", b"0/.", b"1/.", b"./1", b"1/."):
            return 1
        if s in (b"1/1", b"1|1"):
            return 3
        return 2
    toks = _get_line_for(v)
    samples = toks[9:]
    return np.array([_get_gt(s) for s in samples], np.int)

def test_header_info():
    v = VCF(VCF_PATH)
    csq = v[b'CSQ']
    assert csq[b'ID'] == b"CSQ"
    assert b"Description" in csq


    assert_raises(KeyError, v.__getitem__, b'XXXXX')

def test_snpeff_header():
    v = VCF(VCF_PATH2)

    f = v[b'SnpEffVersion']
    assert f != {}, f
    assert b'SnpEffVersion' in f

#def test_info_update():
#    vcf = VCF(VCF_PATH)
#    v = next(vcf)
#    ret = v.INFO.update({'k': 22})
    #print ret
    #assert v.INFO['k'] == 22
    #assert ret == 0, ret

def test_gt_types():
    v = VCF(VCF_PATH)
    for variant in v:

        gt_types = variant.gt_types
        o = _get_samples(variant)
        assert (gt_types == o).all(), (variant, variant.CHROM, variant.POS, zip(gt_types, o))

def test_raw_header():
    v = VCF(VCF_PATH)
    h = v.raw_header.strip().split(b"\n")
    s = h[0]
    assert s == b"##fileformat=VCFv4.1", s
    assert len(h) == 185, len(h)



def test_iterate():

    for i, v in enumerate(VCF(VCF_PATH), start=1):
        pass
    assert i == 115, i


def test_haploid():

    for (gts012, (HOM_REF, HOM_ALT, UNKNOWN)) in ((False, [0, 3, 2]), (True, [0, 2, 3])):
        vcf = VCF("%s/test-haploidX.vcf" % HERE, gts012=gts012)
        for i, v in enumerate(vcf):
            assert not any("/" in b for b in v.gt_bases), (v.start + 1, v.gt_bases)
            if i == 0:
                assert (v.gt_types == [HOM_ALT, HOM_ALT, HOM_ALT]).all(), v.gt_types
            elif v.start == 2800676:
                assert (v.gt_types == [UNKNOWN, HOM_ALT, UNKNOWN]).all(), v.gt_types
            elif v.start == 2832771:
                assert (v.gt_types == [HOM_REF, HOM_ALT, HOM_ALT]).all(), v.gt_types
                break

            if v.start == 2700156:
                assert (v.gt_bases == ['A', 'A', 'A']).all(), v.gt_bases
        break


def test_diploid():
    vcf = VCF("%s/test.comp_het.3.vcf" % HERE)
    for v in vcf:
        if v.start == 17362:
            assert (v.gt_bases == ['TTCT|TTCT', 'TTCT|TTCT', 'TTCT|TTCT',
                                   'TTCT|TTCT', 'TTCT|TTCT', 'TTCT|TTCT', 'TTCT|T', 'TTCT|T',
                                   'TTCT|TTCT', 'TTCT|T', 'TTCT|T',
                                   'TTCT|TTCT']).all(), v.gt_bases

