import numpy as np
from cyvcf2 import VCF, Variant, Writer
import os.path
from nose.tools import assert_raises

HERE = os.path.dirname(__file__)
HEM_PATH = os.path.join(HERE, "test-hemi.vcf")
VCF_PATH = os.path.join(HERE, "test.vcf.gz")



def check_var(v):
    s = [x.split(":")[0] for x in str(v).split("\t")[9:]]
    lookup = {'0/0': 0, '0/1': 1, './1': 1, '1/.': 1, '0/.': 1, './0': 1, '1/1': 3, '.': 2, './.': 2}
    expected = np.array([lookup[ss] for ss in s])
    obs = v.gt_types
    assert np.all(expected == obs), zip(expected, obs)


def test_hemi():
    """
    make sure that we are getting the correct gt_types
    for hemizygous variants
    """

    for p in (HEM_PATH, VCF_PATH):
        vcf = VCF(p)

        for v in vcf:
            check_var(v)
