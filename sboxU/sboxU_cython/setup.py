from setuptools import setup
from distutils.core import Extension
from Cython.Build import cythonize
import os

os.environ["CC"] = "g++"
os.environ["CXX"] = "g++"

module_utils = Extension("utils",
		sources=["utils.pyx"], 
		include_dirs=['.'], 
		language='c++')
module_diff_lin = Extension("diff_lin",
		sources=["diff_lin.pyx"], 
		include_dirs=['.'], 
		language='c++')
module_equiv = Extension("equiv",
		sources=["equiv.pyx"], 
		include_dirs=['.'], 
		language='c++')
module_ccz = Extension("ccz",
		sources=["ccz.pyx"], 
		include_dirs=['.'], 
		language='c++')
module_fp = Extension("fp",
		sources=["fp.pyx"], 
		include_dirs=['.'], 
		language='c++')

setup(name='utils', ext_modules=cythonize([module_utils], language_level = "3"))
setup(name='diff_lin', ext_modules=cythonize([module_diff_lin], language_level = "3"))
setup(name='equiv', ext_modules=cythonize([module_equiv], language_level = "3"))
setup(name='ccz', ext_modules=cythonize([module_ccz], language_level = "3"))
setup(name='fp', ext_modules=cythonize([module_fp], language_level = "3"))
