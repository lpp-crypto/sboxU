from setuptools import setup
from distutils.core import Extension
from Cython.Build import cythonize
import os

os.environ["CC"] = "g++"
os.environ["CXX"] = "g++"

module_utils = Extension("cpp_utils",
		sources=["cpp_utils.pyx"], 
		include_dirs=['.'], 
		language='c++')
module_diff_lin = Extension("cpp_diff_lin",
		sources=["cpp_diff_lin.pyx"], 
		include_dirs=['.'], 
		language='c++')
module_equiv = Extension("cpp_equiv",
		sources=["cpp_equiv.pyx"], 
		include_dirs=['.'], 
		language='c++')
module_ccz = Extension("cpp_ccz",
		sources=["cpp_ccz.pyx"], 
		include_dirs=['.'], 
		language='c++')
module_fp = Extension("cpp_fp",
		sources=["cpp_fp.pyx"], 
		include_dirs=['.'], 
		language='c++')

setup(name='cpp_utils',    ext_modules=cythonize([module_utils],    language_level = "3"))
setup(name='cpp_diff_lin', ext_modules=cythonize([module_diff_lin], language_level = "3"))
setup(name='cpp_equiv',    ext_modules=cythonize([module_equiv],    language_level = "3"))
setup(name='cpp_ccz',      ext_modules=cythonize([module_ccz],      language_level = "3"))
setup(name='cpp_fp',       ext_modules=cythonize([module_fp],       language_level = "3"))
