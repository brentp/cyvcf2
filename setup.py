from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np
import os
import glob

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


setup(
    name="cyvcf2",
    description="fast vcf parsing with cython + htslib",
    url="https://github.com/brentp/cyvcf2/",
    long_description=open("README.md").read(),
    license="MIT",
    author="Brent Pedersen",
    author_email="bpederse@gmail.com",
    version=get_version(),
    ext_modules=cythonize([
        Extension("cyvcf2.cyvcf2", ["cyvcf2/cyvcf2.pyx"],
                  libraries=["hts"],
                  library_dirs= os.environ.get("LD_LIBRARY_PATH", "").split(":"),
                  include_dirs=os.environ.get('C_INCLUDE_PATH', '').split(":") + ["cyvcf2", np.get_include()] ,
                  extra_objects=["cyvcf2/helpers.c"])
        ], include_path=["cyvcf2"]),
    packages=['cyvcf2', 'cyvcf2.tests'],
    test_suite='nose.collector',
    tests_require='nose',
    install_requires=['cython>=0.22.1'],
    include_package_data=True,
)
