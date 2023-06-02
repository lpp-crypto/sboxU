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
extra_compile_args = ["-O3", "-march=native", "-std=c++17"]
extra_link_args=[]

HOME = os.path.expanduser('~')
if os.path.exists(HOME + '/.fftw'):
	extra_compile_args += ["-I" + os.path.expanduser("~") + "/.fftw/include", "-L" + os.path.expanduser("~") + "/.fftw/lib", "-lfftw3_omp", "-lfftw3", "-lm"]
	extra_link_args += ["-L" + os.path.expanduser("~") + "/.fftw/lib", "-lfftw3_omp", "-lfftw3", "-lm"]
	if platform == 'darwin':
		extra_compile_args += ['-lomp', '-I/usr/local/opt/libomp/include']
		extra_link_args += ['-lomp', '-L/usr/local/opt/libomp/include']
	else:
		extra_compile_args += ['-fopenmp']
		extra_link_args += ['-fopenmp']



if os.path.exists(HOME + '/.fftw'):
	module_diff_lin = Extension("cpp_diff_lin",
		                    sources=["cpp_diff_lin.pyx"],
                                    libraries=['m', 'fftw3'],
		                    include_dirs=['.'], 
		                    language='c++',
 			   	    extra_link_args=extra_link_args,
				    extra_compile_args=extra_compile_args)
else:
	module_diff_lin = Extension("cpp_diff_lin_no_fp_lat",
		                    sources=["cpp_diff_lin_no_fp_lat.pyx"],
		                    include_dirs=['.'], 
		                    language='c++',
 			   	    extra_link_args=extra_link_args,
				    extra_compile_args=extra_compile_args)						

module_utils = Extension("cpp_utils",
		         sources=["cpp_utils.pyx"], 
		         include_dirs=['.'],
		         language='c++',
			 extra_link_args=extra_link_args,
                         extra_compile_args=extra_compile_args)

module_equiv = Extension("cpp_equiv",
		         sources=["cpp_equiv.pyx"], 
		         include_dirs=['.'], 
		         language='c++',
			 extra_link_args=extra_link_args,
                         extra_compile_args=extra_compile_args)

module_equiv_approx = Extension("cpp_equiv_approx",
		         sources=["cpp_equiv_approx.pyx"], 
		         include_dirs=['.'], 
		         language='c++',
			 extra_link_args=extra_link_args,
                         extra_compile_args=extra_compile_args)

module_ccz = Extension("cpp_ccz",
		       sources=["cpp_ccz.pyx"], 
		       include_dirs=['.'], 
		       language='c++',
 		       extra_link_args=extra_link_args,
                       extra_compile_args=extra_compile_args)

setup(name='cpp_utils',    ext_modules=cythonize([module_utils],    language_level = "3"))
setup(name='cpp_diff_lin', ext_modules=cythonize([module_diff_lin], language_level = "3"))
setup(name='cpp_equiv',    ext_modules=cythonize([module_equiv],    language_level = "3"))
setup(name='cpp_equiv_approx',    ext_modules=cythonize([module_equiv_approx],    language_level = "3"))
setup(name='cpp_ccz',      ext_modules=cythonize([module_ccz],      language_level = "3"))
