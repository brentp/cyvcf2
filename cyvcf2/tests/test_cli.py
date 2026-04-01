from cyvcf2 import VCF
import os.path

from click.testing import CliRunner

from cyvcf2.cli import cyvcf2
from cyvcf2.cyvcf2 import set_htslib_log_level

HERE = os.path.dirname(__file__)
VCF_PATH = os.path.join(HERE, "test.vcf.gz")
VCF_WITH_WARNINGS_PATH = os.path.join(HERE, "test.warnings.vcf.gz")


#Sanity check for cli
def test_version():
    runner = CliRunner()
    result = runner.invoke(cyvcf2, [VCF_PATH])
    assert result.exit_code == 0


def test_htslib_log_level(capfd):
    runner = CliRunner()

    for log_level in [4, 5]:
        set_htslib_log_level(log_level)
        runner = CliRunner()
        _ = runner.invoke(cyvcf2, [VCF_WITH_WARNINGS_PATH])
        captured = capfd.readouterr()
        assert "[W::vcf_parse_format" in captured.err
        assert "[I::hts_idx_check_local" in captured.err

    set_htslib_log_level(3)
    runner = CliRunner()
    _ = runner.invoke(cyvcf2, [VCF_WITH_WARNINGS_PATH])
    captured = capfd.readouterr()
    assert "[W::vcf_parse_format" in captured.err
    assert "[I::hts_idx_check_local" not in captured.err

    for log_level in [0, 1, 2]:
        set_htslib_log_level(log_level)
        runner = CliRunner()
        _ = runner.invoke(cyvcf2, [VCF_WITH_WARNINGS_PATH])
        captured = capfd.readouterr()
        assert "[W::vcf_parse_format" not in captured.err
        assert "[I::hts_idx_check_local" not in captured.err
