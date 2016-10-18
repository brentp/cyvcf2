from setuptools import setup
from setuptools import Extension
import os
import glob

import sys

if sys.version_info.major == 2 and sys.version_info.minor != 7:
    sys.stderr.write("ERROR: cyvcf2 is only for python 2.7 or greater you are running %d.%d\n", (sys.version_info.major, sys.version_info.minor))
    sys.exit(1)

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

excludes = ['irods', 'plugin']

sources = [x for x in glob.glob('htslib/*.c') if not any(e in x for e in excludes)] + glob.glob('htslib/cram/*.c')
# these have main()'s
sources = [x for x in sources if not x.endswith(('htsfile.c', 'tabix.c', 'bgzip.c'))]
sources.append('cyvcf2/helpers.c')

install_requires = []
try:
    import numpy as np
except ImportError:
    install_requires.append("numpy")

try:
    from Cython.Distutils import build_ext
except ImportError:
    install_requires.append("cython>=0.22.1")

if not install_requires:
    cmdclass = {'build_ext': build_ext}
    extension = [Extension("cyvcf2.cyvcf2",
                           ["cyvcf2/cyvcf2.pyx"] + sources,
                           libraries=['z'],
                           include_dirs=['htslib', 'cyvcf2', np.get_include()])]
else:
    cmdclass = {}
    extension = []

setup(
    name="cyvcf2",
    description="fast vcf parsing with cython + htslib",
    url="https://github.com/brentp/cyvcf2/",
    long_description=open("README.md").read(),
    license="MIT",
    author="Brent Pedersen",
    author_email="bpederse@gmail.com",
    version=get_version(),
    cmdclass=cmdclass,
    ext_modules=extension,
    packages=['cyvcf2', 'cyvcf2.tests'],
    test_suite='nose.collector',
    tests_require='nose',
    install_requires=install_requires,
    include_package_data=True,
    zip_safe=False,
)
