from io import StringIO

from ..cyvcf2 import Writer

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path  # python 2 backport


class TestFileModeInference:
    def test_defaultModeWithVcfFname_returnsUncompressedVcf(self):
        fname = "test.vcf"
        mode = None

        actual = Writer._infer_file_mode(fname, mode)
        expected = "w"

        assert actual == expected

    def test_defaultModeWithBcfFname_returnsCompressedBcf(self):
        fname = "test.bcf"
        mode = None

        actual = Writer._infer_file_mode(fname, mode)
        expected = "wb"

        assert actual == expected

    def test_defaultModeWithBcfFnameAndUncompressedMode_returnsUncompressedBcf(self):
        fname = "test.bcf"
        mode = "wbu"

        actual = Writer._infer_file_mode(fname, mode)
        expected = mode

        assert actual == expected

    def test_defaultModeWithCompressedBcfFname_returnsCompressedBcf(self):
        fname = "test.bcf.gz"
        mode = None

        actual = Writer._infer_file_mode(fname, mode)
        expected = "wb"

        assert actual == expected

    def test_defaultModeWithCompressedVcfFname_returnsCompressedVcf(self):
        fname = "test.vcf.gz"
        mode = None

        actual = Writer._infer_file_mode(fname, mode)
        expected = "wz"

        assert actual == expected

    def test_defaultModeWithStdOut_returnsUncompressedVcf(self):
        fname = "-"
        mode = None

        actual = Writer._infer_file_mode(fname, mode)
        expected = "w"

        assert actual == expected

    def test_defaultModeWithNonVcfName_returnsUncompressedVcf(self):
        fname = "foo"
        mode = None

        actual = Writer._infer_file_mode(fname, mode)
        expected = "w"

        assert actual == expected

    def test_defaultModeWithIntFname_returnsUncompressedVcf(self):
        fname = 1
        mode = None

        actual = Writer._infer_file_mode(fname, mode)
        expected = "w"

        assert actual == expected

    def test_defaultModeWithHandle_returnsUncompressedVcf(self):
        fname = StringIO()
        mode = None

        actual = Writer._infer_file_mode(fname, mode)
        expected = "w"

        assert actual == expected

    def test_defaultModeWithPosixPath_returnsUncompressedVcf(self):
        fname = Path("test.vcf")
        mode = None

        actual = Writer._infer_file_mode(fname, mode)
        expected = "w"

        assert actual == expected

    def test_explicitMode_doesNotInferFromFname(self):
        fname = "test.vcf.gz"
        mode = "wb4"

        actual = Writer._infer_file_mode(fname, mode)
        expected = "wb4"

        assert actual == expected
