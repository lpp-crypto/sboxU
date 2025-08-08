import os
import sys


from setuptools import find_packages, setup
from setuptools.extension import Extension
from Cython.Build import cythonize

if sys.platform == 'darwin':	#macOs
    os.environ["CC"] = "clang"
    os.environ["CXX"] = "clang"
else:
    os.environ["CC"] = "g++"
    os.environ["CXX"] = "g++"
extra_compile_args = ["-O3", "-march=native", "-std=c++17", "-pthread", "-Wno-narrowing"]	#narrowing warnings in fp_lat when calling shape_t{p}
extra_link_args=[]

if sys.platform == 'darwin':
    extra_compile_args += ['-lomp', '-I/usr/local/opt/libomp/include']
    extra_link_args += ['-lomp', '-L/usr/local/opt/libomp/include']
else:
    extra_compile_args += ['-fopenmp']
    extra_link_args += ['-fopenmp']



# !SECTION! Building modules

# core functions over F_2, used throughout the library
f2functions_module = Extension(
        "sboxUv2.f2functions.cython_functions",
        sources = ["sboxUv2/f2functions/cython_functions.pyx" ], #  /!\ name of pyx file must match extension name
        language = "c++",
        extra_link_args = extra_link_args,
        extra_compile_args = extra_compile_args,
)

# core module (S_box and friends)
sbox_module = Extension(
        "sboxUv2.sbox.cython_functions",
        sources = ["sboxUv2/sbox/cython_functions.pyx"], #  /!\ name of pyx file must match extension name
        language = "c++",
        extra_link_args = extra_link_args,
        extra_compile_args = extra_compile_args,
)


# dealing with statistical properties
statistics_module = Extension(
        "sboxUv2.statistics.cython_functions",
        sources = ["sboxUv2/statistics/cython_functions.pyx"], #  /!\ name of pyx file must match extension name
        language = "c++",
        extra_link_args = extra_link_args,
        extra_compile_args = extra_compile_args,
)


# APN functions an friends
apn_module = Extension(
        "sboxUv2.apn.cython_functions",
        sources = ["sboxUv2/apn/cython_functions.pyx"], #  /!\ name of pyx file must match extension name
        language = "c++",
        extra_link_args = extra_link_args,
        extra_compile_args = extra_compile_args,
)

    
# !SECTION! Final setup 
    
setup( # names and others are specified in the pyproject.toml file
    packages = find_packages(),
    ext_modules=cythonize(
        [
                f2functions_module,
                sbox_module,
                statistics_module,
                apn_module,
        ],
        language_level = "3",
    ),
)


