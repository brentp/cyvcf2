from __future__ import print_function
from cyvcf2 import VCF, Variant, Writer
import numpy as np
import os.path
from nose.tools import assert_raises
import tempfile
import sys
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

def test_format_str():
    vcf = VCF(os.path.join(HERE, "test-format-string.vcf"))

    f = next(vcf).format("RULE")
    assert list(f) == ['F', 'G']
    f = next(vcf).format("RULE")
    assert list(f) == ['F2,F3,F4', 'G2,G3,G4']


def test_ibd():
    samples = ['101976-101976', '100920-100920', '100231-100231']
    vcf = VCF(VCF_PATH, gts012=True, samples=samples)
    res = vcf.ibd()
    assert len(res) == 3, (len(res))
    arr = res[('101976-101976', '100920-100920')]
    assert len(arr) > 0

def test_relatedness():
    vcf = VCF(VCF_PATH, gts012=True)
    df = vcf.relatedness(gap=0, linkage_max=2)
    assert "ibs0" in df, df
    assert "rel" in df
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
    for var in vcf(reg):
        k += 1
        assert var.start <= end, var
        assert var.end >= start, var
        assert isinstance(var.REF, basestring)
        assert isinstance(var.ALT, list)
    assert k == 28, k

def test_empty_info():
    for v in VCF(VCF_PHASE_PATH):
        dict(v.INFO)

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
    assert variant.CHROM == "1"
    assert variant.ID is None, variant.ID
    assert variant.start == 10171
    assert variant.end == 10177, variant.end
    assert variant.FILTER is None
    assert variant.QUAL == 92.0

    assert variant.REF == "CCCTAA"
    assert variant.ALT == ["C"]

def test_empty():
    p = os.path.join(HERE, "empty.vcf")
    assert os.path.exists(p)
    assert_raises(IOError, VCF, p)

def test_writer():

    v = VCF(VCF_PATH)
    f = tempfile.mktemp(suffix=".vcf")
    atexit.register(os.unlink, f)

    o = Writer(f, v)
    rec = next(v)
    rec.INFO["AC"] = "3"
    rec.FILTER = ["LowQual"]
    o.write_record(rec)

    rec.FILTER = ["LowQual", "VQSRTrancheSNP99.90to100.00"]
    o.write_record(rec)


    rec.FILTER = "PASS"
    o.write_record(rec)

    o.close()

    expected = ["LowQual".encode(), "LowQual;VQSRTrancheSNP99.90to100.00".encode(), None]

    for i, variant in enumerate(VCF(f)):
        assert variant.FILTER == expected[i], (variant.FILTER, expected[i])

def test_add_info_to_header():
    v = VCF(VCF_PATH)
    v.add_info_to_header({'ID': 'abcdefg', 'Description': 'abcdefg',
        'Type':'Character', 'Number': '1'})
    # NOTE that we have to add the info to the header of the reader,
    # not the writer because the record will be associated with the reader
    f = tempfile.mktemp(suffix=".vcf")
    atexit.register(os.unlink, f)
    w = Writer(f, v)
    import sys
    rec = next(v)

    rec.INFO["abcdefg"] = "XXX"
    w.write_record(rec)
    w.close()

    v = next(VCF(f))
    ret = v.INFO["abcdefg"]
    if isinstance(ret, bytes):
        ret = ret.decode()
    assert ret == "XXX", (dict(v.INFO), v.INFO["abcdefg"])

def test_add_flag():
    vcf = VCF(VCF_PATH)
    vcf.add_info_to_header({'ID': 'myflag', 'Description': 'myflag',
        'Type':'Flag', 'Number': '0'})
    # NOTE that we have to add the info to the header of the reader,
    # not the writer because the record will be associated with the reader
    f = tempfile.mktemp(suffix=".vcf")
    atexit.register(os.unlink, f)
    w = Writer(f, vcf)
    rec = next(vcf)

    rec.INFO["myflag"] = True
    w.write_record(rec)
    w.close()

    v = next(VCF(f))
    assert v.INFO["myflag"] is None, dict(v.INFO)

    f = tempfile.mktemp(suffix=".vcf")
    atexit.register(os.unlink, f)
    w = Writer(f, vcf)
    rec.INFO["myflag"] = False
    w.write_record(rec)
    v = next(VCF(f))
    assert_raises(KeyError, v.INFO.__getitem__, "myflag")


def test_constants():
    v = VCF(VCF_PATH)
    assert v.HOM_REF == 0
    assert v.HET == 1
    assert v.UNKNOWN == 2
    assert v.HOM_ALT == 3

    v = VCF(VCF_PATH, gts012=True)
    assert v.HOM_REF == 0
    assert v.HET == 1
    assert v.HOM_ALT == 2
    assert v.UNKNOWN == 3


def test_add_filter_to_header():
    v = VCF(VCF_PATH)
    # NOTE that we have to add the filter to the header of the reader,
    # not the writer because the record will be associated with the reader
    v.add_filter_to_header({'ID': 'abcdefg', 'Description': 'abcdefg'})

    f = tempfile.mktemp(suffix=".vcf")
    atexit.register(os.unlink, f)
    w = Writer(f, v)
    rec = next(v)

    rec.FILTER = ["abcdefg"]
    w.write_record(rec)
    w.close()

    v = next(VCF(f))
    ret = v.FILTER
    if isinstance(ret, bytes):
        ret = ret.decode()

    assert ret == "abcdefg", v.FILTER

def test_seqnames():
    v = VCF(VCF_PATH)
    assert v.seqnames == ['1'], v.seqnames

    b = VCF('{}/test.snpeff.bcf'.format(HERE))
    assert b.seqnames[0] == 'chr1', b.seqnames
    assert b.seqnames[-1] == 'chrY', b.seqnames


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
        line = line.decode()
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


    assert_raises(KeyError, v.__getitem__, b'XXXXX')

def test_snpeff_header():
    v = VCF(VCF_PATH2)

    f = v['SnpEffVersion']
    assert f != {}, f
    assert 'SnpEffVersion' in f

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
    h = v.raw_header.strip().split("\n")
    s = h[0]
    assert s == "##fileformat=VCFv4.1", s
    assert len(h) == 185, len(h)



def test_iterate():

    for i, v in enumerate(VCF(VCF_PATH), start=1):
        pass
    assert i == 115, i


def test_empty_call():

    for i, v in enumerate(VCF(VCF_PATH)(), start=1):
        pass
    assert i == 115, i
    del i

    for i, v in enumerate(VCF(VCF_PATH)(''), start=1):
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


def test_format():

    vcf = VCF('{}/test.vcf.gz'.format(HERE))
    for v in vcf:
        a = v.format('PL', int)
        assert a.shape == (189, 3)

        a = v.format('SB', int)
        assert a is None

        a = v.format('DP', int)
        assert a.shape == (189, 1)

def test_header_stuff():
    vcf = VCF('{}/test.vcf.gz'.format(HERE))
    import sys
    seen_formats, seen_infos = 0, 0
    for h in vcf.header_iter():
        i = h.info(extra=True)
        assert isinstance(i, dict)
        seen_formats += i['HeaderType'] == 'FORMAT'
        seen_infos += i['HeaderType'] == 'INFO'
    assert seen_formats == 9, seen_formats
    assert seen_infos == 73, seen_infos


def test_bcf():
    vcf = VCF('{}/test.snpeff.bcf'.format(HERE))
    l = sum(1 for _ in vcf)
    assert l == 10, l

    # NOTE: this is 0 becuase we don't SEEK.
    l = sum(1 for _ in vcf())
    assert l == 0, l


    viter = vcf("1:69260-69438")
    sys.stderr.write("\nOK\n")
    sys.stderr.flush()
    l = list(viter)
    assert len(l) == 0, len(l)

    iter = vcf("chr1:69260-69438")
    l = list(iter)
    assert len(l) == 2, len(l)


def test_issue12():
    fields = "ADP_ALL ADPD ADPO ADP_PASS ADPR AFR AMBIG BMF_PASS BMF_QUANT AF_FAILED FA_FAILED FM_FAILED FP_FAILED FR_FAILED MD_FAILED IMPROPER MQ_FAILED OVERLAP PV_FAILED QSS".split()

    vcf = VCF('{}/bug.vcf.gz'.format(HERE))
    for v in vcf:
        for f in fields:
            vals = v.format(f)
            if vals is not None:
                assert vals.dtype in (np.int32, np.float32), (f, vals.dtype)

        vals = v.format("RVF")
        assert vals.dtype in (np.float32, np.float64)

    assert_raises(KeyError, v.format, "RULE")

def test_gt_bases_nondiploid():
    """Ensure gt_bases works with more complex base representations.
    """
    vcf = VCF('{}/test_gt_bases.vcf.gz'.format(HERE))
    expected = {0: ['C/C', 'C/C'], 1: ['C/T', 'C'], 2: ['C', 'C/T'], 3: ['C', 'C']}
    for i, v in enumerate(vcf):
        assert v.gt_bases.tolist() == expected[i], (v.gt_bases.tolist(), expected[i])
