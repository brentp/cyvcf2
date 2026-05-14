from ._version import __version__
from .cyvcf2 import (VCF, Variant, Writer, r_ as r_unphased, par_relatedness,
                     par_het)
Reader = VCFReader = VCF
