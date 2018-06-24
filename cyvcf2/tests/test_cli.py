from cyvcf2 import VCF
import os.path

from click.testing import CliRunner

from cyvcf2.cli import cyvcf2

HERE = os.path.dirname(__file__)
VCF_PATH = os.path.join(HERE, "test.vcf.gz")


#Sanity check for cli
def test_version():
    runner = CliRunner()
    result = runner.invoke(cyvcf2, [VCF_PATH])
    assert result.exit_code == 0
