import os
import sys
import re


from setuptools import find_packages, setup
from setuptools.extension import Extension
from Cython.Build import cythonize


# !SECTION! Setting up the compilation of the C++ part

if sys.platform == 'darwin':	#macOs
    os.environ["CC"] = "clang"
    os.environ["CXX"] = "clang"
else:
    os.environ["CC"] = "g++"
    os.environ["CXX"] = "g++"
extra_compile_args = ["-O3", "-march=native", "-std=c++20", "-pthread", "-Wno-narrowing"]	#narrowing warnings in fp_lat when calling shape_t{p}
extra_link_args=[]

if sys.platform == 'darwin':
    extra_compile_args += ['-lomp', '-I/opt/homerew/opt/libomp/include']
    extra_link_args += ['-lomp', '-L/opt/homebrew/opt/libomp/lib']
else:
    extra_compile_args += ['-fopenmp']
    extra_link_args += ['-fopenmp']



# !SECTION! Declaring cython extensions


def declare_cython(full_module_name):
    src = re.sub(r"\.", r"/", full_module_name) + ".pyx"
    return Extension(
        full_module_name,
        sources = [ src ],
        language = "c++",
        extra_link_args = extra_link_args,
        extra_compile_args = extra_compile_args,
    )


all_cython_extensions = [ declare_cython(name) for name in [
    
    "sboxUv2.core.f2functions.cython_functions",
    "sboxUv2.core.sbox.cython_functions",
    "sboxUv2.algorithms.cython_functions",
    "sboxUv2.statistics.cython_functions",
    "sboxUv2.ccz.cython_functions",
    "sboxUv2.ccz.affine_equivalence.cython_functions",
    "sboxUv2.core.spectrum.cython_functions",
    "sboxUv2.core.anf.cython_functions"
]]

    
# !SECTION! Final setup 
    
setup( # names and others are specified in the pyproject.toml file
    packages = find_packages(),
    ext_modules=cythonize(
        all_cython_extensions,
        language_level = "3",
    ),
)


