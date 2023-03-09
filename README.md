# SboxU v1.1

## Description

SAGE/Python functions useful for studying S-boxes and Boolean
functions such as computing the DDT, computing the Walsh spectrum,
affine equivalence testing... Some of them are implemented in C++
and/or in parallel to improve their performances.


## Dependencies

Most functions in `sboxU` only depend on a recent version of SAGE (it was tested with
version 9.0).

**/!\ Functions related to linear approximations of sboxes defined in fields of base prime p != 2 require the installation of the `fftw` library found at http://www.fftw.org/ and `openmp`. They must be downloaded and installed prior to running `setup.py` using the instructions below, in particular with the same directory prefix for fftw.** LAT-related functions now accept an additional argument p (assumed to be equal to 2 by default).  `fftw3` must be configured with the flags `--enable-shared --enable-openmp --prefix $HOME/.fftw` (as seen below).


### **Optional**: Installing openmp and fftw on Ubuntu

Installing openmp can be done with:

    sudo apt-get install libomp-dev

To download fftw, first go into your working directory, then input the following commands:

    wget http://www.fftw.org/fftw-3.3.10.tar.gz
    tar -xf fftw-3.3.10.tar.gz
    cd fftw-3.3.10/

The commands to install fftw with the right options are, once in the `fftw-*/` directory:

    ./configure --enable-shared --enable-openmp --prefix $HOME/.fftw
    make
    make install

### **Optional**: Installing openmp and fftw on macOS

Installing openmp can be done with:

    brew install libomp

To download fftw, first go into your working directory, then input the following commands:

    curl -O http://www.fftw.org/fftw-3.3.10.tar.gz
    tar -xf fftw-3.3.10.tar.gz
    cd fftw-3.3.10/

Installing fftw requires `gfortran` (`brew install gfortran`), as `clang` can have trouble linking with openmp:

    ./configure CC="gfortran" --enable-shared --enable-openmp --prefix $HOME/.fftw
    make
    make install


## Usage

To retrieve this module, use the following command:

    git clone https://github.com/lpp-crypto/sboxU/

Then, you need to compile the content of the `sboxU/sboxU_cython` folder.
This process only relies on SAGE To compile it:

    cd sboxU/sboxU_cython
    sage setup.py build_ext --inplace

This compiles the C++ part of `sboxU` and allows it to be called from
a SAGE script. To use it in your project, simply move the `sboxU`
folder to your project's directory. You can then import `sboxU` like
any other python module.  As an example of the functions provided by
`sboxU`, the SAGE script `example.py` stored alongside the folder
`sboxU` generates random permutations and tests their affine
equivalence.
    
Then, in order to use the functions provided by sboxU, simply
paste/make a link to the folder sboxU from the directory in which you
are working. For example, if you have a script called `test.py` in
which you want to call the functions in sboxU then your directory
should look like this:

    $ ls
    test.py
    README.md
    sboxu/
    
and the file `test.py` should contain `from sboxU import *`. You can
also import sboxU in a sage notebook session provided that you started
sage from the folder containing it.


## Troubleshooting

- `undefined symbol: _Py_ZeroStruct` If you get such an error, there
  are probably mismatches between the python versions used. This
  typically happens if you have 2 versions of SAGE installed, used the
  most recent one to compile `sboxU`, but then use an older one in
  your scripts.

## Authors

[Mathias Joly](https://github.com/MathiasJoly)

[Léo Perrin](https://who.paris.inria.fr/Leo.Perrin/)

[Aurélien Boeuf](https://who.paris.inria.fr/Aurelien.Boeuf/)
