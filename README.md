# SboxU v1.2

## Description

SAGE/Python functions useful for studying S-boxes and Boolean
functions such as computing the DDT, computing the Walsh spectrum,
affine equivalence testing... Some of them are implemented in C++
and/or in parallel to improve their performances.


## Dependencies

Most functions in `sboxU` only depend on a recent version of SAGE (it was tested with
version 9.0). Some use `openmp` for multithreading, and you may need to install it in order to successfully compile. 

**UPDATE (v1.2)**: `fftw` is no longer required. It was replaced with `pocketfft`, which is contained entirely within `pocketfft_hdronly.h`.

Installing openmp can be done with:

    sudo apt-get install libomp-dev

or, on macOS:

    brew install libomp


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

## Acknowledgements

When operating over fields of odd characteristic, we use the [pocketfft] (https://gitlab.mpcdf.mpg.de/mtr/pocketfft) library (`pocketfft_hdronly.h` was taken from the `cpp` branch) and we would like to thank their authors.
