import os
from cyvcf2 import VCF
from click.testing import CliRunner
from cyvcf2.cli import cyvcf2
from cyvcf2.cyvcf2 import set_htslib_log_level, _init_logging_from_env

HERE = os.path.dirname(__file__)
VCF_PATH = os.path.join(HERE, "test.vcf.gz")
VCF_WITH_WARNINGS_PATH = os.path.join(HERE, "test.warnings.vcf.gz")


# Sanity check for cli
def test_version():
    runner = CliRunner()
    result = runner.invoke(cyvcf2, [VCF_PATH])
    assert result.exit_code == 0


def test_htslib_log_level(capfd):
    runner = CliRunner()

    for log_level in [4, 5]:
        set_htslib_log_level(log_level)
        _ = runner.invoke(cyvcf2, [VCF_WITH_WARNINGS_PATH])
        captured = capfd.readouterr()
        assert "[W::vcf_parse_format" in captured.err
        assert "[I::hts_idx_check_local" in captured.err

    set_htslib_log_level(3)
    _ = runner.invoke(cyvcf2, [VCF_WITH_WARNINGS_PATH])
    captured = capfd.readouterr()
    assert "[W::vcf_parse_format" in captured.err
    assert "[I::hts_idx_check_local" not in captured.err

    for log_level in [0, 1, 2]:
        set_htslib_log_level(log_level)
        _ = runner.invoke(cyvcf2, [VCF_WITH_WARNINGS_PATH])
        captured = capfd.readouterr()
        assert "[W::vcf_parse_format" not in captured.err
        assert "[I::hts_idx_check_local" not in captured.err


def test_htslib_log_level_env(monkeypatch, capfd):
    runner = CliRunner()

    monkeypatch.setenv("CYVCF2_HTSLIB_LOG_LEVEL", "5")
    _init_logging_from_env()
    _ = runner.invoke(cyvcf2, [VCF_WITH_WARNINGS_PATH])
    captured = capfd.readouterr()
    assert "[W::vcf_parse_format" in captured.err
    assert "[I::hts_idx_check_local" in captured.err

    monkeypatch.setenv("CYVCF2_HTSLIB_LOG_LEVEL", "3")
    _init_logging_from_env()
    _ = runner.invoke(cyvcf2, [VCF_WITH_WARNINGS_PATH])
    captured = capfd.readouterr()
    assert "[W::vcf_parse_format" in captured.err
    assert "[I::hts_idx_check_local" not in captured.err

    monkeypatch.setenv("CYVCF2_HTSLIB_LOG_LEVEL", "0")
    _init_logging_from_env()
    _ = runner.invoke(cyvcf2, [VCF_WITH_WARNINGS_PATH])
    captured = capfd.readouterr()
    assert "[W::vcf_parse_format" not in captured.err
    assert "[I::hts_idx_check_local" not in captured.err

    monkeypatch.delenv("CYVCF2_HTSLIB_LOG_LEVEL", raising=False)
    _init_logging_from_env()
    _ = runner.invoke(cyvcf2, [VCF_WITH_WARNINGS_PATH])
    captured = capfd.readouterr()
    assert "[W::vcf_parse_format" in captured.err
    assert "[I::hts_idx_check_local" not in captured.err
