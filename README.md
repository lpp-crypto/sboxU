# SboxU v0.1


## Disclaimer

`sboxU` is at a very early stage. Its API might change over time!

## Description

SAGE/Python functions useful for studying S-boxes and Boolean
functions such as computing the DDT, computing the Walsh spectrum,
affine equivalence testing... Some of them are implemented in C++
and/or in parallel to improve their performances.


## Dependencies

The `SboxU` was only tested on Linux (Ubuntu 16.04). To install it,
you need the following packages.

```
libboost-python-dev
libpython-dev
sage
cmake
```

## Usage

To retrieve this module, use the following command:

git clone https://github.com/lpp-crypto/sboxU/

Then, move to the `sbox-utils/sboxU` directory and run:

```
cmake .
make
```

This compiles the C++ part of `sboxU` and allows it to be called from
a SAGE script. To use it in your project, simply move the `sboxU`
folder to your project's directory. You can then import `sboxU` like
any other python module.  As an example of the functions provided by
`sboxU`, the SAGE script `example.py` stored alongside the folder
`sboxU` generates random permutations and tests their affine
equivalence.


# For MAC Users

`sboxU` should work out of the box on a MAC. If Sage complains about
`ValueError: unknown locale: UTF-8 in Python`, try adding the
following lines to your `~/.bash_profile`:

```
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
```

To install the boost library needed, you can use the following
command:

```
brew install boost-python
```

