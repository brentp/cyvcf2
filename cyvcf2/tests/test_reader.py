from __future__ import print_function
import os.path
import tempfile
import sys
import os
import atexit
try:
  from pathlib import Path
except ImportError:
  from pathlib2 import Path  # python 2 backport
import warnings

import numpy as np
import pytest

from ..cyvcf2 import VCF, Variant, Writer


HERE = os.path.dirname(__file__)
VCF_PATH = os.path.join(HERE, "test.vcf.gz")
VCF_PATH2 = os.path.join(HERE, "test.snpeff.vcf")
VCF_PATH3 = os.path.join(HERE, "test.mnp.vcf")
VCF_PHASE_PATH = os.path.join(HERE, "test.comp_het.3.vcf")
VCF_ALTFREQ_PATH = os.path.join(HERE, "test_gt_alt_freqs.vcf")

try:
    basestring
except NameError:
    basestring = (str, bytes)

def test_init():
    # string
    v = VCF(VCF_PATH)
    assert v
    expected_count = sum(1 for _ in v)
    v.close()

    # Path
    v = VCF(Path(VCF_PATH))
    value = sum(1 for _ in v)
    assert value == expected_count

    # file descriptor
    with open(VCF_PATH) as fp:
        fd = fp.fileno()
        v = VCF(fd)
        assert sum(1 for _ in v) == expected_count
        v.close()  # this should not close the file descriptor originally opened

    # file-like object
    with open(VCF_PATH) as fp:
        v = VCF(fp)
        assert sum(1 for _ in v) == expected_count
        v.close()  # this should not close the file descriptor originally opened

def test_type():
    vcf = VCF(VCF_PATH)
    for v in vcf:
        if len(v.REF) == 1 and len(v.ALT[0]) == 1:
            assert v.var_type == 'snp'
        elif v.ALT[0][0] != "<":
            assert v.var_type == 'indel'
        else:
            print(v.var_type, v.REF, v.ALT)

def test_type_mnp():
    vcf = VCF(VCF_PATH3)
    for v in vcf:
        if len(v.ALT) == 1:
            if (v.REF, v.ALT[0]) in [("CGT","CGG"), ("AGG","CGA")]:
                assert v.var_type == 'mnp'
            if (v.REF, v.ALT[0]) in [("GCA","GA")]:
                assert v.var_type == 'indel'
        if len(v.ALT) == 2:
            if (v.REF, v.ALT[0], v.ALT[1]) in [("TCGGT","GCGGG","GCGGT")]:
                assert v.var_type == 'mnp'
            if (v.REF, v.ALT[0], v.ALT[1]) in [("ATTAC","ATAC","AATAC")]:
                assert v.var_type == 'indel'
        if len(v.ALT) == 3:
            if (v.REF, v.ALT[0], v.ALT[1], v.ALT[2]) in [("GCC","TCC","GCA","GCG")]:
                assert v.var_type == 'mnp'

def test_format_str():
    vcf = VCF(os.path.join(HERE, "test-format-string.vcf"))

    f = next(vcf).format("RULE")
    assert list(f) == ['F', 'G']
    f = next(vcf).format("RULE")
    assert list(f) == ['F2,F3,F4', 'G2,G3,G4']

def test_missing_samples():
    samples = ['101976-101976', 'sample_not_in_vcf']
    with warnings.catch_warnings(record=True) as w:
        vcf = VCF(VCF_PATH, gts012=True, samples=samples)
        assert len(w) == 1
        assert "not all requested samples found" in str(w[-1].message)
    assert len(vcf.samples) == 1
    vcf.close()
    samples = '101976-101976,sample_not_in_vcf'
    with warnings.catch_warnings(record=True) as w:
        vcf = VCF(VCF_PATH, gts012=True, samples=samples)
        assert len(w) == 1
        assert "not all requested samples found" in str(w[-1].message)
    assert len(vcf.samples) == 1

def test_ibd():
    samples = ['101976-101976', '100920-100920', '100231-100231']
    vcf = VCF(VCF_PATH, gts012=True, samples=samples)
    res = vcf.ibd()
    assert len(res) == 3, (len(res))
    arr = res[('101976-101976', '100920-100920')]
    assert len(arr) > 0

def test_relatedness():
    vcf = VCF(VCF_PATH, gts012=True)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("ignore", category=DeprecationWarning)
        df = vcf.relatedness(gap=0, linkage_max=2)
        assert len(w) == 1
        assert "tested:" in str(w[-1].message)
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

def test_gt_alt_freqs():
    vcf = VCF(VCF_ALTFREQ_PATH)

    v = next(vcf)
    assert v.gt_alt_freqs[0] == 0.2
    assert v.gt_alt_freqs[1] == 1.0

    v = next(vcf)
    assert v.gt_alt_freqs[0] == 0.5
    assert v.gt_alt_freqs[1] == 0.9

    v = next(vcf)
    assert v.gt_alt_freqs[0] == 0.0
    assert v.gt_alt_freqs[1] == 0.0

    v = next(vcf)
    assert v.gt_alt_freqs[0] == 0.0
    assert v.gt_alt_freqs[1] == 1.0

    v = next(vcf)
    assert v.gt_alt_freqs[0] == 0.0
    assert v.gt_alt_freqs[1] == 0.0

    v = next(vcf)
    assert v.gt_alt_freqs[0] == -1
    assert v.gt_alt_freqs[1] == -1

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
    with pytest.raises(Exception):
        VCF("XXXXX")

def test_samples():
    v = VCF(VCF_PATH)
    assert len(v.samples) == 189

def test_next():
    v = VCF(VCF_PATH)
    variant = next(v)
    assert isinstance(variant, Variant)

def test_variant():
    with pytest.raises(TypeError):
        Variant()

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
    with pytest.raises(IOError):
        VCF(p)

def test_format_field():
    vcf = VCF(VCF_PATH)
    for v in vcf:
        assert isinstance(v.FORMAT, list)

def test_writer_from_string():

    header = """##fileformat=VCFv4.1
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr2,length=249250621,assembly=hg19>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samplea
"""

    w = Writer.from_string("out.vcf", header)
    w.write_header()
    v = w.variant_from_string("chr1\t234\t.\tA\tC\t40\tPASS\t.\tGT\t0/0")
    w.write_record(v)
    w.close()

def test_isa():
    vcf = VCF(os.path.join(HERE, "test.isa.vcf"))
    for i, v in enumerate(vcf):
        if i in {0, 1, 2, 3}:
            assert v.is_indel
            assert not v.is_snp
        if i in {4}:
            assert v.is_snp




def run_writer(writer, filename, rec):
    rec.INFO["AC"] = "3"
    rec.FILTER = ["LowQual"]
    writer.write_record(rec)

    rec.FILTER = ["LowQual", "VQSRTrancheSNP99.90to100.00"]
    writer.write_record(rec)

    rec.FILTER = "PASS"
    writer.write_record(rec)

    writer.close()

    expected = ["LowQual", "LowQual;VQSRTrancheSNP99.90to100.00", None]

    for i, variant in enumerate(VCF(filename)):
        assert variant.FILTER == expected[i], (variant.FILTER, expected[i])

def test_writer():
    v = VCF(VCF_PATH)
    f = tempfile.mktemp(suffix=".vcf")
    atexit.register(os.unlink, f)
    rec = next(v)

    # string
    run_writer(Writer(f, v), f, rec)

    # Path
    path = Path(f)
    run_writer(Writer(path, v), f, rec)

    # file descriptor
    with open(VCF_PATH) as fp:
        fd = fp.fileno()
        run_writer(Writer(fd, v), f, rec)

    # file-like object
    with open(VCF_PATH) as fp:
        run_writer(Writer(fp, v), f, rec)

def test_filters():
    v = VCF(VCF_PATH)
    f = tempfile.mktemp(suffix=".vcf")
    atexit.register(os.unlink, f)
    rec = next(v)

    writer = Writer(f, v)

    rec.INFO["AC"] = "3"
    rec.FILTER = ["LowQual"]
    writer.write_record(rec)

    rec.FILTER = ["LowQual", "VQSRTrancheSNP99.90to100.00"]
    writer.write_record(rec)

    rec.FILTER = "PASS"
    writer.write_record(rec)

    rec.FILTER = []
    writer.write_record(rec)

    writer.close()

    expected_filter = ["LowQual", "LowQual;VQSRTrancheSNP99.90to100.00", None, None]
    expected_filters = [["LowQual"], ["LowQual", "VQSRTrancheSNP99.90to100.00"], ["PASS"], []]

    for i, variant in enumerate(VCF(f)):
        assert variant.FILTER == expected_filter[i], (variant.FILTER, expected_filter[i])
        assert variant.FILTERS == expected_filters[i], (variant.FILTERS, expected_filters[i])

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

def test_read_flag():
    vcf = VCF(VCF_PATH)
    for v in vcf:
        assert ("in_exac_flag" in str(v)) == v.INFO.get('in_exac_flag', False)

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

    fh = VCF(f)
    v = next(fh)
    fh.close()
    assert v.INFO["myflag"] is True, dict(v.INFO)

    f = tempfile.mktemp(suffix=".vcf")
    atexit.register(os.unlink, f)
    w = Writer(f, vcf)
    rec.INFO["myflag"] = False
    w.write_record(rec)
    w.close()
    fh = VCF(f)
    v = next(fh)
    fh.close()
    with pytest.raises(KeyError):
        v.INFO["myflag"]

def test_issue198():
    vcf = VCF(os.path.join(HERE, "issue_198.vcf"), strict_gt=True)
    for v in vcf:
        assert all(v.gt_types == [2]), v.gt_types

    vcf = VCF(os.path.join(HERE, "issue_198.vcf"), strict_gt=False)
    for v in vcf:
        assert all(v.gt_types == [0]), v.gt_types

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
    assert v.seqnames ==  [u'1', u'2', u'3', u'4', u'5', u'6', u'7', u'8', u'9', u'10', u'11', u'12', u'13', u'14', u'15', u'16', u'17', u'18', u'19', u'20', u'21', u'22', u'X', u'Y', u'MT', u'GL000207.1', u'GL000226.1', u'GL000229.1', u'GL000231.1', u'GL000210.1', u'GL000239.1', u'GL000235.1', u'GL000201.1', u'GL000247.1', u'GL000245.1', u'GL000197.1', u'GL000203.1', u'GL000246.1', u'GL000249.1', u'GL000196.1', u'GL000248.1', u'GL000244.1', u'GL000238.1', u'GL000202.1', u'GL000234.1', u'GL000232.1', u'GL000206.1', u'GL000240.1', u'GL000236.1', u'GL000241.1', u'GL000243.1', u'GL000242.1', u'GL000230.1', u'GL000237.1', u'GL000233.1', u'GL000204.1', u'GL000198.1', u'GL000208.1', u'GL000191.1', u'GL000227.1', u'GL000228.1', u'GL000214.1', u'GL000221.1', u'GL000209.1', u'GL000218.1', u'GL000220.1', u'GL000213.1', u'GL000211.1', u'GL000199.1', u'GL000217.1', u'GL000216.1', u'GL000215.1', u'GL000205.1', u'GL000219.1', u'GL000224.1', u'GL000223.1', u'GL000195.1', u'GL000212.1', u'GL000222.1', u'GL000200.1', u'GL000193.1', u'GL000194.1', u'GL000225.1', u'GL000192.1', u'NC_007605', u'hs37d5', u'phix'], v.seqnames

    b = VCF('{}/test.snpeff.bcf'.format(HERE), threads=3)
    assert b.seqnames[0] == 'chr1', b.seqnames
    assert b.seqnames[-1] == 'chrY', b.seqnames

def test_different_index():
    b = VCF('{}/test.snpeff.bcf'.format(HERE), threads=3)
    b.set_index("{}/test-diff.csi".format(HERE))
    s = 0
    for r in b("chr1:69427-69429"):
        s += 1
    assert s == 1

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

    with gzip.open(VCF_PATH) as f:
        for i, line in enumerate(f, start=1):
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
        if s in ("0/0", "0|0", "0/.", "./0", "0|."):
            return 0
        if s in ("0/1", "0|1", "1/.", "./1", "1/."):
            return 1
        if s in ("1/1", "1|1"):
            return 3
        return 2
    toks = _get_line_for(v)
    samples = toks[9:]
    return np.array([_get_gt(s) for s in samples], np.int32)

def test_header_info():
    v = VCF(VCF_PATH)
    csq = v['CSQ']
    assert csq['ID'] == "CSQ"
    assert "Description" in csq


    with pytest.raises(KeyError):
        v[b'XXXXX']

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
    with warnings.catch_warnings(record=True) as w:
        l = list(viter)
        assert len(w) == 1
        assert "no intervals found" in str(w[-1].message)
    assert len(l) == 0, len(l)

    iter = vcf("chr1:69260-69438")
    l = list(iter)
    assert len(l) == 2, len(l)

def test_vcf_no_intervals():
    vcf = VCF('{}/test.vcf.gz'.format(HERE))
    viter = vcf("not_a_chrom")
    with warnings.catch_warnings(record=True) as w:
        l = list(viter)
        assert len(w) == 1
        assert "no intervals found" in str(w[-1].message)
    assert len(l) == 0, len(l)

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

    with pytest.raises(KeyError):
        v.format("RULE")

def test_gt_bases_nondiploid():
    """Ensure gt_bases works with more complex base representations.
    """
    vcf = VCF('{}/test_gt_bases.vcf.gz'.format(HERE))
    expected = {0: ['C/C', 'C/C'], 1: ['C/T', 'C'], 2: ['C', 'C/T'], 3: ['C', 'C']}
    for i, v in enumerate(vcf):
        assert v.gt_bases.tolist() == expected[i], (v.gt_bases.tolist(), expected[i])

def fmap(fn, a):
    return list(map(fn, a))

def allclose(a, b):
    return np.allclose(np.array(a, dtype=float), np.array(b, dtype=float))

def test_set_format_float():
    vcf = VCF('{}/test-format-string.vcf'.format(HERE))
    assert vcf.add_format_to_header(dict(ID="PS", Number=1, Type="Float", Description="PS example")) == 0
    v = next(vcf)
    v.set_format("PS", np.array([0.555, 1.111], dtype=np.float32))
    assert allclose(fmap(float, get_gt_str(v, "PS")), np.array([0.555, 1.111]))

    v.set_format("PS", np.array([8.555, 11.111], dtype=np.float64))
    assert allclose(fmap(float, get_gt_str(v, "PS")), [8.555, 11.111])

    v.set_format("PS", np.array([9998.555, 99911.111], dtype=np.float32))
    obs = fmap(float, get_gt_str(v, "PS"))
    assert allclose(obs, [9998.555, 99911.111]), obs

def test_set_format_int_a():
    vcf = VCF('{}/test-format-string.vcf'.format(HERE))
    assert vcf.add_format_to_header(dict(ID="PI", Number=1, Type="Integer", Description="Int example")) == 0
    v = next(vcf)
    v.set_format("PI", np.array([5, 1], dtype=int))
    assert allclose(fmap(float, get_gt_str(v, "PI")), [5, 1])

def test_set_format_int_b():
    vcf = VCF('{}/test-format-string.vcf'.format(HERE))
    assert vcf.add_format_to_header(dict(ID="PI", Number=1, Type="Integer", Description="Int example")) == 0
    v = next(vcf)

    v.set_format("PI", np.array([855, 11], dtype=np.int64))
    assert allclose(fmap(float, get_gt_str(v, "PI")), [855, 11])

def test_set_format_int_c():
    vcf = VCF('{}/test-format-string.vcf'.format(HERE))
    assert vcf.add_format_to_header(dict(ID="PI", Number=1, Type="Integer", Description="Int example")) == 0
    v = next(vcf)

    v.set_format("PI", np.array([9998, 99911], dtype=np.int32))
    obs = fmap(float, get_gt_str(v, "PI"))
    assert allclose(obs, [9998, 99911]), obs

def test_set_format_int3():
    "test that we can handle multiple (in this case 3) values per sample"
    vcf = VCF('{}/test-format-string.vcf'.format(HERE))
    assert vcf.add_format_to_header(dict(ID="P3", Number=3, Type="Integer", Description="Int example")) == 0
    v = next(vcf)
    exp = np.array([[1, 11, 111], [2, 22, 222]], dtype=int)
    v.set_format("P3", exp)
    res = get_gt_str(v, "P3")
    assert res == ["1,11,111", "2,22,222"], (res, str(v))

    assert np.allclose(v.format("P3"), exp)

def test_set_format_str_bytes_second_longer():
    vcf = VCF('{}/test-format-string.vcf'.format(HERE))
    assert vcf.add_format_to_header(dict(ID="STR", Number=1, Type="String", Description="String example")) == 0
    v = next(vcf)

    v.set_format("STR", np.array([b'foo', b'barbaz']))
    assert np.all(v.format('STR') == np.array(['foo', 'barbaz']))

def test_set_format_str_bytes_first_longer():
    vcf = VCF('{}/test-format-string.vcf'.format(HERE))
    assert vcf.add_format_to_header(dict(ID="STR", Number=1, Type="String", Description="String example")) == 0
    v = next(vcf)

    v.set_format("STR", np.array([b'foobar', b'baz']))
    assert np.all(v.format('STR') == np.array(['foobar', 'baz']))

def test_set_format_str_bytes_number3():
    # Confirm currently not supported
    vcf = VCF('{}/test-format-string.vcf'.format(HERE))
    assert vcf.add_format_to_header(dict(ID="STR", Number=3, Type="String", Description="String example")) == 0
    v = next(vcf)

    contents = np.array([[b'foo', b'barbaz', b'biz'], [b'blub', b'bloop', b'blop']])
    with pytest.raises(Exception):
        v.set_format("STR", contents)

def test_set_gts():
    vcf = VCF('{}/test-format-string.vcf'.format(HERE))
    v = next(vcf)

    v.genotypes = [[1, 1, True], [0, 0, False]]
    assert get_gt_str(v) == ["1|1", "0/0"]

    v.genotypes = [[-1, 1, False], [-1, 0, False]]
    assert get_gt_str(v) == ["./1", "./0"]

    v.genotypes = [[-1, -1, False], [0, 0, True]]
    assert get_gt_str(v) == ["./.", "0|0"]

    v.genotypes = [[2, 2, True], [0, 2, True]]
    assert get_gt_str(v) == ["2|2", "0|2"]

    v.genotypes = [[0, 1, 2, False], [1, 2, True]]
    s = get_gt_str(v)
    assert s == ["0/1/2", "1|2"]

    v.genotypes = [[1, 2, False], [0, 1, 2, True]]
    assert get_gt_str(v) == ["1/2", "0|1|2"]

    v.genotypes = [[0, 1, 2, False], [0, 1, 2, True]]
    assert get_gt_str(v) == ["0/1/2", "0|1|2"]

def test_info_del():
    vcf = VCF(os.path.join(HERE, "test-hemi.vcf"))
    v = next(vcf)
    d = str(v)
    assert ";DP=" in d
    del v.INFO["DP"]
    d = str(v)
    assert not ";DP=" in d

def test_filter_id():
    vcf = VCF(os.path.join(HERE, "test-hemi.vcf"))
    v = next(vcf)
    assert v.ID == "ID1"
    assert v.FILTER == "FAIL"


def get_gt_str(variant, key="GT"):
    idx = variant.FORMAT.index(key)
    return [x.split(":")[idx] for x in str(variant).strip().split("\t")[9:]]

def test_access_gts():
    vcf = VCF('{}/test-format-string.vcf'.format(HERE))
    """
7	55086956	.	C	G	0	.	.	GT:ADP_ALL:RULE	0/0:6728,1:F	1|1:22,1:G
7	55086957	.	T	A,C,G	0	.	.	GT:ADP_ALL:RULE	1/2:6768,2,2,1:F2,F3,F4	2|3:1,2,3,4:G2,G3,G4
7	55086958	.	T	G	0	.	.	GT:ADP_ALL:RULE	0/1/.:6768,2,2,1:F2,F3,F4	0:1,2,3,4:G2,G3,G4
7	55086959	.	T	G,T	0	.	.	GT:ADP_ALL:RULE	.	0|2:1,2,3,4:G2,G3,G4
    """

    v = next(vcf)
    gts = v.genotypes
    assert gts == [[0, 0, False], [1, 1, True]], gts

    v = next(vcf)
    assert v.genotypes == [[1, 2, False], [2, 3, True]], v.genotypes

    v = next(vcf)
    assert v.genotypes == [[0, 1, -1, False], [0, True]], v.genotypes

    v = next(vcf)
    assert v.genotypes == [[-1, True], [0, 2, True]], v.genotypes

def test_access_genotype():
    vcf = VCF('{}/test-format-string.vcf'.format(HERE))
    v = next(vcf)
    gts = v.genotype
    #print(str(v), file=sys.stderr)

    # indexing directly gives a list of Allele objects (for diploid, the list
    # has length 2)
    alleles = gts[0]
    assert alleles[0].value == 0
    assert alleles[1].value == 0
    assert alleles[0].phased == False

    alleles = gts[1]
    assert alleles[0].value == 1
    assert alleles[1].value == 1
    assert alleles[1].phased == True

    assert np.all(gts.array()[:, :2] == np.array([[0, 0], [1, 1]]))
    assert np.all(gts.array()[:, 2] == np.array([False, True]))

    assert alleles[1].phased == True
    alleles[1].phased = False
    assert alleles[1].phased == False
    alleles = gts[1]
    assert alleles[1].phased == False

    assert np.all(gts.array()[:, 2] == np.array([False, False]))

    # can also just get the phased stats of the nth sample:
    assert gts.phased(0) == False
    # note this got updated above
    assert gts.phased(1) == False

    # and the alleles of the nth sample.
    assert gts.alleles(0) == [0, 0]
    #print(gts.alleles(1), file=sys.stderr)
    assert gts.alleles(1) == [1, 1]


    alleles = gts[0]
    assert alleles[0].value == 0
    alleles[0].value = 1
    assert alleles[0].value == 1
    alleles = gts[0]
    assert alleles[0].value == 1
    assert alleles[0].phased == False

    assert np.all(gts.array()[:, :2] == np.array([[1, 0], [1, 1]]))

    gts[1][0].value = 0

    assert np.all(gts.array()[:, :2] == np.array([[1, 0], [0, 1]]))

    # update the varint
    v.genotype = gts
    assert "1/0:6728,1:F	0/1:22,1:G" in str(v)

def test_access_genotype_array():
    vcf = VCF('{}/test-format-string.vcf'.format(HERE))
    """
7	55086956	.	C	G	0	.	.	GT:ADP_ALL:RULE	0/0:6728,1:F	1|1:22,1:G
7	55086957	.	T	A,C,G	0	.	.	GT:ADP_ALL:RULE	1/2:6768,2,2,1:F2,F3,F4	2|3:1,2,3,4:G2,G3,G4
7	55086958	.	T	G	0	.	.	GT:ADP_ALL:RULE	0/1/.:6768,2,2,1:F2,F3,F4	0:1,2,3,4:G2,G3,G4
7	55086959	.	T	G,T	0	.	.	GT:ADP_ALL:RULE	.	0|2:1,2,3,4:G2,G3,G4
    """

    v = next(vcf)
    np.testing.assert_array_equal(
        v.genotype.array(),
        np.array([[0, 0, 0], [1, 1, 1]], dtype=np.int16)
    )

    v = next(vcf)
    np.testing.assert_array_equal(
        v.genotype.array(),
        np.array([[1, 2, 0], [2, 3, 1]], dtype=np.int16)
    )

    v = next(vcf)
    np.testing.assert_array_equal(
        v.genotype.array(),
        np.array([[0, 1, -1, 0], [0, -2, -2, 1]], dtype=np.int16)
    )

    v = next(vcf)
    np.testing.assert_array_equal(
        v.genotype.array(),
        np.array([[-1, -2, 1], [0, 2, 1]], dtype=np.int16)
    )

    # test fill value
    np.testing.assert_array_equal(
        v.genotype.array(fill=-9),
        np.array([[-1, -9, 1], [0, 2, 1]], dtype=np.int16)
    )

def test_alt_homozygous_gt():
    vcf = VCF(os.path.join(HERE, "test-multiallelic-homozygous-alt.vcf.gz"))
    assert vcf is not None
    v = next(vcf)
    assert v
    assert v.gt_bases[0] == '<*:DEL>/<*:DEL>'

    vcf = VCF(os.path.join(HERE, "test-multiallelic-homozygous-alt.vcf.gz"), gts012=True)
    assert vcf is not None
    v = next(vcf)
    assert v
    assert v.gt_bases[0] == '<*:DEL>/<*:DEL>'

def test_write_missing_contig():
    input_vcf = VCF('{}/seg.vcf.gz'.format(HERE))
    output_vcf = Writer('/dev/null', input_vcf)
    for v in input_vcf:
        v.genotypes = [[1,1,False]]
        output_vcf.write_record(v)
    output_vcf.close()

def test_set_samples():
    vcf = VCF(VCF_PATH)
    assert len(vcf.samples) == 189, len(vcf.samples)
    vcf.set_samples([vcf.samples[2]])
    assert len(vcf.samples) == 1
    v = next(vcf)
    assert len(v.gt_types) == 1

def test_hrec():

    vcf = VCF(VCF_PATH)
    for item in vcf.header_iter():
        info = item.info()
        if info['HeaderType'] != 'GENERIC':
            assert 'ID' in info

def test_issue44():
    vcf = VCF('{}/issue_44.vcf'.format(HERE))
    w = Writer('__o.vcf', vcf)
    for v in vcf:
        tmp = v.genotypes
        #print(tmp, file=sys.stderr)
        v.genotypes = tmp
        w.write_record(v)
    w.close()
    #           "./."            "."          ".|."           "0|0"
    expected = [[-1, -1, False], [-1, False], [-1, -1, True], [0, 0, True]]
    #print("", file=sys.stderr)
    for i, v in enumerate(VCF('__o.vcf')):
        #print(v.genotypes, file=sys.stderr)
        assert v.genotypes == [expected[i]], (i, v.genotypes, expected[i])
    os.unlink("__o.vcf")

def test_id_field_updates():
    # 1 10172   .       CCCTAA  C       92.0    PASS
    v = VCF(VCF_PATH)
    variant = next(v)
    assert variant.ID is None, variant.ID

    variant.ID = 'foo'
    assert variant.ID == 'foo', variant.ID

    variant.ID = 100
    assert variant.ID == '100', variant.ID

    variant.ID = 100.1
    assert variant.ID == '100.1', variant.ID

    variant.ID = '.'
    assert variant.ID is None, variant.ID

    variant.ID = None
    assert variant.ID is None, variant.ID

def test_set_pos():
    test_vcf = '{}/test-strict-gt-option-flag.vcf.gz'.format(HERE)
    vcf = VCF(test_vcf, gts012=False)
    v = next(vcf)

    orig_pos, orig_start = v.POS, v.start
    v.set_pos(22)
    assert v.start == 22
    assert v.POS == 23

def test_set_chrom_when_contig_not_in_header():
    test_vcf = '{}/test-strict-gt-option-flag.vcf.gz'.format(HERE)
    new_chrom = "NEW"
    vcf = VCF(test_vcf, gts012=False)
    original_seqnames = vcf.seqnames
    assert new_chrom not in original_seqnames
    v = next(vcf)

    with warnings.catch_warnings(record=True) as w:
        v.CHROM = new_chrom
        assert len(w) == 1
        assert "added new contig" in str(w[-1].message)
    assert v.CHROM == new_chrom
    expected_seqnames = sorted(original_seqnames + [new_chrom])
    assert vcf.seqnames == expected_seqnames

def test_set_chrom_after_contig_is_added_to_header():
    test_vcf = '{}/test-strict-gt-option-flag.vcf.gz'.format(HERE)
    new_chrom = "NEW"
    vcf = VCF(test_vcf, gts012=False)
    original_seqnames = vcf.seqnames
    vcf.add_to_header("##contig=<ID={},length=15>".format(new_chrom))
    expected_seqnames = sorted(original_seqnames + [new_chrom])
    assert vcf.seqnames == expected_seqnames
    v = next(vcf)

    v.CHROM = new_chrom
    assert v.CHROM == new_chrom

def test_set_qual():
    v = VCF(VCF_PATH)
    variant = next(v)
    assert variant.QUAL == 92.0

    variant.QUAL = 30.0
    assert variant.QUAL == 30.0

    with pytest.raises(TypeError):
        variant.QUAL = "30.0"

    variant.QUAL = None
    assert variant.QUAL is None, 'variant.QUAL is {}'.format(variant.QUAL)

def test_strict_gt_option_flag():
    test_vcf = '{}/test-strict-gt-option-flag.vcf.gz'.format(HERE)

     #T, C
     #0/0     1/1     0/1     1/0     1/.     ./1     0/.     ./0     ./.

    truth_gt_bases = ('T/T', 'C/C', 'T/C', 'C/T', 'C/.', './C', 'T/.', './T', './.')
    truth_genotypes = (
        [0, 0, False],
        [1, 1, False],
        [0, 1, False],
        [1, 0, False],
        [1, -1, False],
        [-1, 1, False],
        [0, -1, False],
        [-1, 0, False],
        [-1, -1, False],
    )

    vcf = VCF(test_vcf, gts012=False)
    variant = next(vcf)

    msg = "VCF(gts012=False, strict_gt=False) not working"
    truth_gt_types = (0, 3, 1, 1, 1, 1, 0, 0, 2)
    print(variant.gt_bases.tolist(), file=sys.stderr)
    print(variant.gt_types.tolist(), file=sys.stderr)
    assert bool(tuple(variant.gt_bases.tolist()) == truth_gt_bases), '{} [gt_bases]'.format(msg)

    """
    assert bool(tuple(variant.gt_types.tolist()) == truth_gt_types), '{} [gt_types]'.format(msg)
    assert bool(tuple(variant.genotypes) == truth_genotypes), '{} (genotypes)'.format(msg)

    vcf = VCF(test_vcf, gts012=False, strict_gt=True)
    variant = next(vcf)

    msg = "VCF(gts012=False, strict_gt=True) not working"
    truth_gt_types = (0, 3, 1, 1, 2, 2, 2, 2, 2)
    assert bool(tuple(variant.gt_bases.tolist()) == truth_gt_bases), '{} [gt_bases]'.format(msg)
    assert bool(tuple(variant.gt_types.tolist()) == truth_gt_types), '{} [gt_types]'.format(msg)
    assert bool(tuple(variant.genotypes) == truth_genotypes), '{} (genotypes)'.format(msg)


    vcf = VCF(test_vcf, gts012=True)
    variant = next(vcf)

    msg = "VCF(gts012=True, strict_gt=False) not working"
    truth_gt_types = (0, 2, 1, 1, 1, 1, 0, 0, 3)
    assert tuple(variant.gt_bases.tolist()) == truth_gt_bases, '{} [gt_bases]'.format(msg)
    #sys.stderr.write("\nobs:%s\n" % variant.gt_types.tolist())
    #sys.stderr.write("exp:%s\n" % list(truth_gt_types))
    assert tuple(variant.gt_types.tolist()) == truth_gt_types, '{} [gt_types]'.format(msg)
    assert tuple(variant.genotypes) == truth_genotypes, '{} (genotypes)'.format(msg)

    vcf = VCF(test_vcf, gts012=True, strict_gt=True)
    variant = next(vcf)

    msg = "VCF(gts012=True, strict_gt=True) not working"
    truth_gt_types = (0, 2, 1, 1, 3, 3, 3, 3, 3)
    assert tuple(variant.gt_bases.tolist()) == truth_gt_bases, '{} [gt_bases]'.format(msg)
    assert tuple(variant.gt_types.tolist()) == truth_gt_types, '{} [gt_types]'.format(msg)
    assert tuple(variant.genotypes) == truth_genotypes, '{} (genotypes)'.format(msg)
    """

def test_alt_repr():
    v = os.path.join(HERE, "test-alt-repr.vcf")
    vcf = VCF(v, gts012=True, strict_gt=False)
    v = next(vcf)
    assert np.all(v.gt_types == np.array([0, 1, 2, 0, 1, 3]))

    v = os.path.join(HERE, "test-alt-repr.vcf")
    vcf = VCF(v, gts012=False, strict_gt=False)
    v = next(vcf)
    assert np.all(v.gt_types == np.array([0, 1, 3, 0, 1, 2]))


"""
def test_seqlens():
    v = VCF(VCF_PATH)
    assert v.seqlens == [249250621, 243199373, 198022430, 191154276,
        180915260, 171115067, 159138663, 146364022, 141213431, 135534747,
        135006516, 133851895, 115169878, 107349540, 102531392, 90354753,
        81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560,
        59373566, 16569, 4262, 15008, 19913, 27386, 27682, 33824, 34474, 36148,
        36422, 36651, 37175, 37498, 38154, 38502, 38914, 39786, 39929, 39939,
        40103, 40531, 40652, 41001, 41933, 41934, 42152, 43341, 43523, 43691,
        45867, 45941, 81310, 90085, 92689, 106433, 128374, 129120, 137718,
        155397, 159169, 161147, 161802, 164239, 166566, 169874, 172149, 172294,
        172545, 174588, 179198, 179693, 180455, 182896, 186858, 186861, 187035,
        189789, 191469, 211173, 547496, 171823, 35477943, 5386], v.seqlens
        """

def test_closed_iter():
    path = os.path.join(HERE, "test-alt-repr.vcf")
    vcf = VCF(path, gts012=True, strict_gt=False)
    vcf.close()

    with pytest.raises(Exception):
        next(vcf)

def test_issue72():
    path = os.path.join(HERE, "test-alt-repr.vcf")
    vcf = VCF(path, gts012=True, strict_gt=False)

    v = next(vcf)
    assert v.INFO['DQ'] == 1

    assert v.format('DQ') is not None

def test_is_transition():
    vcf = VCF(VCF_ALTFREQ_PATH)

    for r in vcf:
        assert r.is_transition

def test_decomposed():
    vcf = VCF(os.path.join(HERE, "decomposed.vcf"))
    v = next(vcf)
    #0/.	./0	1/.	./1	./.
    assert np.all(v.gt_types == np.array([vcf.HOM_REF, vcf.HOM_REF, vcf.HET, vcf.HET, vcf.UNKNOWN]))


def test_fd():

    fh = open(os.path.join(HERE, "decomposed.vcf"))
    fn = fh.fileno()

    vcf = VCF(fn)
    v = next(vcf)
    assert np.all(v.gt_types == np.array([vcf.HOM_REF, vcf.HOM_REF, vcf.HET, vcf.HET, vcf.UNKNOWN]))
    fh.close()
    vcf.close()


def test_set_reference():
    fh = open(os.path.join(HERE, "decomposed.vcf"))
    fn = fh.fileno()

    vcf = VCF(fn)
    for v in vcf:
      v.REF = "CCCCA"
      assert  "CCCCA" in str(v)

def test_set_alternates():
    fh = open(os.path.join(HERE, "decomposed.vcf"))
    fn = fh.fileno()

    vcf = VCF(fn)
    for v in vcf:
      v.ALT = "TTT,GGG"
      assert  "TTT,GGG" in str(v)

      v.ALT = ["AAAC", "CCCA"]
      assert "AAAC,CCCA" in str(v)

def test_no_seqlen():

    vcf_path = os.path.join(HERE, "no-seq-len.vcf")
    vcf = VCF(vcf_path)
    assert vcf.seqnames == ["3"]
    with pytest.raises(AttributeError):
        vcf.seqlens

def test_set_unknown_format():
    vcf = VCF(VCF_PATH)
    vcf.add_format_to_header({'ID':'NEW', 'Type':'Float', 'Number':1, 'Description':'...'})

    v = next(vcf)
    arr = np.array([[1.1] for s in vcf.samples])
    arr[-1][0] = np.nan
    v.set_format('NEW', arr)
    record = str(v)
    parts = record.split()
    assert parts[-1][-1] == '.'
    assert parts[-2][-3:] == '1.1'

def test_invalid_header():
    # htslib produces the error "Empty sample name: trailing spaces/tabs in the header line?"
    p = os.path.join(HERE, "test-invalid-header.vcf")
    assert os.path.exists(p)
    with pytest.raises(Exception):
        VCF(p)


def test_genotypes():
    """
. 1
./. ./1
./. 1
0/1 0
0 0/1
     """
    exp_array = [[-1,  0],
       [-1, -1,  0],
       [-1, -1,  0],
       [0, 1,  0],
       [0, -2, 1],
       ]

    non_strict_exp_num = [
    [0, 0, 1, 1],
    [0, 1, 1, 0],
    [0, 0, 1, 1],
    [1, 1, 0, 0],
    [1, 1, 0, 0],
    ] 

    strict_exp_num = [x[:] for x in non_strict_exp_num]
    strict_exp_num[1] = [0, 0, 2, 0] # both unknown

    for strict_gt in (False, True):
        vcf = VCF(os.path.join(HERE, "test-genotypes.vcf"), strict_gt=strict_gt)

        exp_num = strict_exp_num if strict_gt else non_strict_exp_num


        for i, v in enumerate(vcf):
            #if i != 3: continue
            obs = [v.num_hom_ref, v.num_het, v.num_unknown, v.num_hom_alt]
            assert obs == exp_num[i], ("error with num_*")

            a = v.genotype.array()[0] # only 0'th item
            print("i:", i, " a:", v.genotype.array()[0], " exp:", exp_array[i])
            assert (a == exp_array[i]).all(), " error with array"

def test_issue236():


    for strict_gt in (False, True):
        res = [[], []]

        for gts in (False, True):
            vcf = VCF(os.path.join(HERE, "test-genotypes.vcf"), gts012=gts, strict_gt=strict_gt)
            for v in vcf:
                res[int(gts)].append(v.num_unknown)

        assert res[0] == res[1]
        assert len(res[0]) > 0


def test_issue17_no_gt():
    vcf = VCF(os.path.join(HERE, "test-no-genotypes.vcf"))
    with pytest.raises(Exception):
        for v in vcf:
            v.num_called  # Used to give segmentation fault
