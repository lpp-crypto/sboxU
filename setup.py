from setuptools import setup
from distutils.core import Extension
from Cython.Build import cythonize
import os
from sys import platform

if platform == 'darwin':	#macOs
	os.environ["CC"] = "clang"
	os.environ["CXX"] = "clang"
else:
	os.environ["CC"] = "g++"
	os.environ["CXX"] = "g++"
extra_compile_args = ["-O3", "-march=native", "-std=c++17", "-pthread", "-Wno-narrowing"]	#narrowing warnings in fp_lat when calling shape_t{p}
extra_link_args=[]

if platform == 'darwin':
	extra_compile_args += ['-lomp', '-I/usr/local/opt/libomp/include']
	extra_link_args += ['-lomp', '-L/usr/local/opt/libomp/include']
else:
	extra_compile_args += ['-fopenmp']
	extra_link_args += ['-fopenmp']

setup(
    name = "sboxU",
    version = "1.3.0",
    description = "SAGE/Python functions useful for studying S-boxes and Boolean functions such as computing the DDT, computing the Walsh spectrum, affine equivalence testing...",
    packages = ['sboxU', 'sboxU.known_functions', 'sboxU.sboxU_cython'],
    ext_modules=cythonize(
        [
            Extension(
                f"sboxU.sboxU_cython.{name}",
                sources = [ f"sboxU/sboxU_cython/{name}.pyx" ],
                include_dirs = [ "sboxU/sboxU_cython/" ],
                language = "c++",
                extra_link_args = extra_link_args,
                extra_compile_args = extra_compile_args,
            )
                for name in [ "cpp_diff_lin", "cpp_utils", "cpp_equiv", "cpp_equiv_approx", "cpp_ccz" ]
        ],
        language_level = "3",
    ),
)
