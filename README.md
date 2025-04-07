# SboxU v1.3

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


## Install

### Direct Installation via SAGE

To use this module, use the following command:

    sage --pip install git+https://github.com/lpp-crypto/sboxU

This compiles the C++ part of `sboxU` and installs the full `sboxU` module
in your SAGE's python environment. You can then import `sboxU` like
any other python module.

### Installation via `nix`

To create a development shell that contains SAGE with `sboxU` installed, use the following command:

    nix develop github:lpp-crypto/sboxU

The nix package and development shell are currently only implemented for `x86_64-linux` systems.

## Example Usage

As an example of the functions provided by `sboxU`, you can download and run the SAGE script [`example.py`](https://raw.githubusercontent.com/lpp-crypto/sboxU/refs/heads/master/example.py).
There `sboxU` generates random permutations and tests their affine
equivalence.

On a unix based system (e.g. Linux or macOS), use the following commands:

    wget https://raw.githubusercontent.com/lpp-crypto/sboxU/refs/heads/master/example.py
    sage example.py


To write your own script using `sboxU`, this file should contain `from sboxU import *`.
You can also import sboxU in a sage notebook session.


## Troubleshooting

- `undefined symbol: _Py_ZeroStruct` If you get such an error, there
  are probably mismatches between the python versions used. This
  typically happens if you have 2 versions of SAGE installed, used the
  most recent one to compile `sboxU`, but then use an older one in
  your scripts.

## Authors (in alphabetical order)

[Jens Alich](https://informatik.rub.de/ac-personen/alich/)

[Jules Baudrin](https://who.paris.inria.fr/Jules.Baudrin/)

[Aurélien Boeuf](https://who.paris.inria.fr/Aurelien.Boeuf/)

[Xavier Bonnetain](https://bonneta.in/)

[Alain Couvreur](http://www.lix.polytechnique.fr/Labo/Alain.Couvreur/)

[Mathias Joly](https://github.com/MathiasJoly)

[Léo Perrin](https://who.paris.inria.fr/Leo.Perrin/)

## Acknowledgements

When operating over fields of odd characteristic, we use the [pocketfft] (https://gitlab.mpcdf.mpg.de/mtr/pocketfft) library (`pocketfft_hdronly.h` was taken from the `cpp` branch) and we would like to thank their authors.
