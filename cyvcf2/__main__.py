import sys

from .cli import cyvcf2 as cli

"""
cyvcf2.__main__
~~~~~~~~~~~~~~~~~~~~~
The main entry point for the command line interface.
Invoke as ``cyvcf2`` (if installed)
or ``python -m cyvcf2`` (no install required).
"""

if __name__ == "__main__":
    # exit using whatever exit code the CLI returned
    sys.exit(cli())