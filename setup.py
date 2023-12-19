import ctypes
import glob
import os
import platform
import shutil
import sys
import subprocess

from setuptools import setup, Extension, Command
from setuptools.command.build_ext import build_ext
from setuptools.command.sdist import sdist

if sys.version_info.major == 2 and sys.version_info.minor != 7:
    sys.stderr.write(
        "ERROR: cyvcf2 is only for python 2.7 or greater you are running %d.%d\n",
        (sys.version_info.major, sys.version_info.minor),
    )
    sys.exit(1)

import numpy as np


def get_version():
    """Get the version info from the mpld3 package without importing it"""
    import ast

    with open(os.path.join("cyvcf2", "__init__.py"), "r") as init_file:
        module = ast.parse(init_file.read())

    version = (
        ast.literal_eval(node.value)
        for node in ast.walk(module)
        if isinstance(node, ast.Assign) and node.targets[0].id == "__version__"
    )
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


def check_libhts():
    os_type = platform.system()
    if os_type == "Linux":
        lib_name = "libhts.so"
    elif os_type == "Darwin":  # macOS
        lib_name = "libhts.dylib"
    elif os_type == "Windows":
        lib_name = "hts-3.dll"
    else:
        return False  # Unsupported OS

    try:
        ctypes.CDLL(lib_name)
        return True
    except Exception:
        return False


def build_htslib(htslib_configure_options, static_mode):
    current_directory = os.getcwd()
    os.chdir(os.path.join(current_directory, "htslib"))

    if os.path.exists("config.status"):
        print("# cyvcf2: config.status exists, skip configure htslib")
    else:
        subprocess.run(["autoreconf", "-i"], check=True)

        configure_args = ["./configure"]
        if static_mode:
            configure_args.append("CFLAGS=-fPIC")
        if htslib_configure_options:
            configure_args.extend(htslib_configure_options.split())

        subprocess.run(configure_args, check=True)
    subprocess.run(["make"], check=True)

    os.chdir(current_directory)


def pre_sdist():
    current_directory = os.getcwd()
    os.chdir(os.path.join(current_directory, "htslib"))

    # generate version.h
    subprocess.run(["make", "htscodecs/htscodecs/version.h"], check=True)

    # remove redundant file
    redudant_files = ["htscodecs.mk"]
    for redudant_file in redudant_files:
        if os.path.exists(redudant_file):
            os.remove(redudant_file)

    os.chdir(current_directory)


class cyvcf2_build_ext(build_ext):
    def run(self):
        print("# cyvcf2: htslib mode is {}".format(CYVCF2_HTSLIB_MODE))
        if CYVCF2_HTSLIB_MODE == "BUILTIN" or not check_libhts():
            if platform.system() == "Windows":
                # Windows htslib can't be built internall
                raise RuntimeError(
                    "Required library 'hts-3.dll' not found. Please install htslib first."
                )

            print(
                "# cyvcf2: htslib configure options is {}".format(
                    CYVCF2_HTSLIB_CONFIGURE_OPTIONS
                )
            )
            build_htslib(
                CYVCF2_HTSLIB_CONFIGURE_OPTIONS, CYVCF2_HTSLIB_MODE == "BUILTIN"
            )

        if CYVCF2_HTSLIB_MODE == "BUILTIN":
            # add htslib linked libraries
            # libcrypto is bundled with libssl, capabale to replace libssl
            all_dynamic_libs = {
                "z",
                "bz2",
                "lzma",
                "curl",
                "deflate",
                "crypto",
            }

            extra_libs = []
            # read the htslib config.status file to get the linked libraries
            with open(os.path.join("htslib", "config.status"), "r") as f:
                for line in f:
                    if 'S["static_LIBS"]' in line:
                        linked_libs_str = line.split("=")[1].strip()[1:-1]
                        print("# cyvcf2: htslib librarys is {}".format(linked_libs_str))
                        linked_libs = linked_libs_str.split()
                        for lib in linked_libs:
                            if lib[2:] in all_dynamic_libs:  # remove -l prefix
                                extra_libs.append(lib[2:])

            self.extensions[0].libraries = self.extensions[0].libraries + extra_libs

        super().run()


class cyvcf2_sdist(sdist):
    def run(self):
        if platform.system() == "Windows":
            raise RuntimeError("cyvcf2 Can't build sdist on Windows")

        pre_sdist()

        super().run()


class clean_ext(Command):
    description = "clean up Cython temporary files and htslib build files"
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        print("cleaning build files")
        if os.path.exists("build"):
            shutil.rmtree("build")

        if os.path.exists("cyvcf2.egg-info"):
            shutil.rmtree("cyvcf2.egg-info")

        print("cleaning Cython temporary files")
        cyvcf2_c_path = os.path.join("cyvcf2", "cyvcf2.c")
        if os.path.exists(cyvcf2_c_path):
            os.remove(cyvcf2_c_path)

        lib_files = glob.glob("cyvcf2/cyvcf2.cpython-*")
        for file in lib_files:
            os.remove(file)

        if platform.system() != "Windows":
            print("cleaning htslib build files")
            current_directory = os.getcwd()
            os.chdir(os.path.join(current_directory, "htslib"))
            subprocess.run(["make", "distclean"], check=True)
            os.chdir(current_directory)


# How to link against HTSLIB
# BUILTIN:  build and static link against htslib from
#           builtin htslib code (default)
# EXTERNAL: use shared libhts.so compiled outside of
#           cyvcf2
if platform.system() == "Windows":
    # can't static link to htslib on Windows
    htslib_mode_default = "EXTERNAL"
else:
    htslib_mode_default = "BUILTIN"

CYVCF2_HTSLIB_MODE = os.environ.get("CYVCF2_HTSLIB_MODE", htslib_mode_default)

if platform.system() == "Windows" and CYVCF2_HTSLIB_MODE == "BUILTIN":
    print(
        "# cyvcf2 WARNING: The use of cyvcf2 on Windows is experimental. It will not work when statically linked to htslib. Fallback to htslib EXTERNAL mode"
    )
    CYVCF2_HTSLIB_MODE = "EXTERNAL"

CYVCF2_HTSLIB_CONFIGURE_OPTIONS = os.environ.get(
    "CYVCF2_HTSLIB_CONFIGURE_OPTIONS", None
)

htslib_objects = []
htslib_librarys = []
htslib_library_dirs = []

if CYVCF2_HTSLIB_MODE == "BUILTIN":
    htslib_objects = ["htslib/libhts.a"]
else:
    htslib_librarys = ["hts"]
    if not check_libhts():
        htslib_library_dirs = ["htslib"]

htslib_include_dirs = ["htslib", "htslib/htslib"]

# Build the Cython extension by statically linking to the bundled htslib
sources = ["cyvcf2/cyvcf2.pyx", "cyvcf2/helpers.c"]

extensions = [
    Extension(
        "cyvcf2.cyvcf2",
        sources,
        extra_objects=htslib_objects,
        libraries=htslib_librarys,
        extra_compile_args=[
            "-Wno-sign-compare",
            "-Wno-unused-function",
            "-Wno-strict-prototypes",
            "-Wno-unused-result",
            "-Wno-discarded-qualifiers",
        ],
        include_dirs=["cyvcf2", np.get_include()] + htslib_include_dirs,
        library_dirs=htslib_library_dirs,
    )
]


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
    if os.path.exists('htslib/config.status') and CYVCF2_HTSLIB_MODE == "BUILTIN":
        os.unlink('htslib/config.status')
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
    packages=["cyvcf2"],
    entry_points=dict(
        console_scripts=[
            "cyvcf2 = cyvcf2.__main__:cli",
        ],
    ),
    python_requires=">=3.7",
    install_requires=["numpy", "coloredlogs", "click"],
    include_package_data=True,
    zip_safe=False,
    cmdclass={
        "clean_ext": clean_ext,
        "build_ext": cyvcf2_build_ext,
        "sdist": cyvcf2_sdist,
    },
)
