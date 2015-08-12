from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import glob

setup(
    name="cyvcf2",
    ext_modules=cythonize([
        Extension("cyvcf2.cyvcf2", ["cyvcf2/cyvcf2.pyx"], libraries=["hts"],
            include_dirs=["cyvcf2"], extra_objects=["cyvcf2/helpers.c",
                "cyvcf2/helpers.h"] + glob.glob("cyvcf2/blosc*") +
            glob.glob("cyvcf2/shuff*"))
     ], include_path=["cyvcf2"]),
    packages=['cyvcf2', 'cyvcf2.tests'],
    #test_suite='nose.collector',
    #tests_require='nose',
    include_package_data=True,
)
