import os
import glob
import sys
import subprocess
import platform

import pkg_resources
from setuptools import setup, Extension, dist

if sys.version_info.major == 2 and sys.version_info.minor != 7:
    sys.stderr.write("ERROR: cyvcf2 is only for python 2.7 or greater you are running %d.%d\n", (sys.version_info.major, sys.version_info.minor))
    sys.exit(1)

# Install numpy right now
dist.Distribution().fetch_build_eggs(['numpy'])
import numpy as np


def get_version():
    """Get the version info from the mpld3 package without importing it"""
    import ast

    with open(os.path.join("cyvcf2", "__init__.py"), "r") as init_file:
      module = ast.parse(init_file.read())

    version = (ast.literal_eval(node.value) for node in ast.walk(module)
         if isinstance(node, ast.Assign)
         and node.targets[0].id == "__version__")
    try:
      return next(version)
    except StopIteration:
          raise ValueError("version could not be located")


def no_cythonize(extensions, **_ignore):
    for extension in extensions:
        sources = []
        for sfile in extension.sources:
            path, ext = os.path.splitext(sfile)
            if ext in (".pyx", ".py"):
                sfile = path + ".c"
            sources.append(sfile)
        extension.sources[:] = sources
    return extensions


# Build the Cython extension by statically linking to the bundled htslib
sources = [
    x for x in glob.glob('htslib/*.c') 
    if not any(e in x for e in ['irods', 'plugin'])
]
sources += glob.glob('htslib/cram/*.c')
# Exclude the htslib sources containing main()'s
sources = [x for x in sources if not x.endswith(('htsfile.c', 'tabix.c', 'bgzip.c'))]
sources.append('cyvcf2/helpers.c')

extra_libs = []
if platform.system() != 'Darwin':
    extra_libs.append('crypt')
if bool(int(os.getenv("LIBDEFLATE", 0))):
    extra_libs.append('deflate')

extensions = [Extension("cyvcf2.cyvcf2",
                        ["cyvcf2/cyvcf2.pyx"] + sources,
                        libraries=['z', 'bz2', 'lzma', 'curl', 'ssl'] + extra_libs,
                        extra_compile_args=["-Wno-sign-compare", "-Wno-unused-function",
                            "-Wno-strict-prototypes",
                            "-Wno-unused-result", "-Wno-discarded-qualifiers"],
                        include_dirs=['htslib', 'cyvcf2', np.get_include()])]


CYTHONIZE = bool(int(os.getenv("CYTHONIZE", 0)))
if CYTHONIZE:
    try:
        from Cython.Build import cythonize
    except ImportError:
        sys.stderr.write(
            "Cannot find Cython. Have you installed all the requirements?\n"
            "Try pip install -r requirements.txt\n"
        )
        sys.exit(1)
    compiler_directives = {"language_level": 2, "embedsignature": True}
    extensions = cythonize(extensions, compiler_directives=compiler_directives)
else:
    extensions = no_cythonize(extensions)


setup(
    name="cyvcf2",
    description="fast vcf parsing with cython + htslib",
    long_description_content_type="text/markdown",
    url="https://github.com/brentp/cyvcf2/",
    long_description=open("README.md").read(),
    license="MIT",
    author="Brent Pedersen",
    author_email="bpederse@gmail.com",
    version=get_version(),
    ext_modules=extensions,
    packages=['cyvcf2', 'cyvcf2.tests'],
    entry_points=dict(
        console_scripts=[
            'cyvcf2 = cyvcf2.__main__:cli',
        ],
    ),
    install_requires=['numpy', 'coloredlogs', 'click'],
    include_package_data=True,
    zip_safe=False,
)
