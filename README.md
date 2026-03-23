# sboxU

![The logo of sboxU, showing a box being bult or disassembled.](./docs/logo-v2-5.png)

## Description


`sboxU` is a SAGE/Python library that is intended to systematize knowledge about the black-box analysis of vectorial Boolean functions, p-ary fuctions, and in particular about the algorithms relevant to their study. To this end, it provides a wide variety of fuctions performing for instance the following tasks:
- generating S-boxes from univariate polynomials,
- multi-threaded computation of the Walsh spectrum,
- display of the "Pollock representation" on a DDT,
- computation of the automorphisms of the graph of an APN funotion,
- identify all the bijections in the CCZ-equivalence class of a function,
- compuptation iof the non-linear invariants of a function,
- ... and much, much more!

At its core, `sboxU` is a C++ library providing convenient abstractions for S-boxe, affine maps, etc.; as well as the algorithms operating on them. Then, a cython layer exposes these functions to SAGE.

If you use `sboxU` in a published paper, please cite it using the following bibtex entry:

```
    @misc{sboxU,
    authors={Léo Perrin,
    Jules Baudrin,
    Xavier Bonnetain,
    Alain Couvreur,
    Merlin Fruchon,
    Mathias Joly,
    Pierre Galissant,
    Lukas Stennes
    },
    year=2026,
    title={sbox{U}: black-box analysis of discrete functions},
    howpublished={Available online at \url{https://github.com/lpp-crypto/sboxU/}}
    }
```



### Documentation 

Some tests/examples are provided in the `tests` folder. You must compile and install `sboxU` as explained below in order for them to work. 

- The SAGE API is documented [here](https://who.paris.inria.fr/Leo.Perrin/code/sboxU/sage/) (note in particular the search box!).
- If you need direct access to its C++ internals, a (still quite incomplete) API documentation is available [here](https://who.paris.inria.fr/Leo.Perrin/code/sboxU/cpp/).



### sboxU_CPP

The `C++` component of `sboxU` can be used on its own, see [the relevant folder](./sboxUv2/cpp/README.md).

### Executable Scripts

During its installation, `sboxU` creates some executables that are installed in your path. To see their interfacce, simply run them with the `-h` argument.

- **sboxU_apn_db_generation** uses known lists of APN fuctions to generate a local `tinySQL` database containing exactly one representstive per EA-class of functions 


## Installing SboxU


### Dependencies

Most functions in `sboxU` only depend on a recent version of SAGE (it was tested for example with version 10.5). Some use `openmp` for multithreading, and you may need to install it in order to successfully compile. 

Installing openmp can be done with:

    sudo apt-get install libomp-dev

or, on macOS:

    brew install libomp


### Downloading and Compiling


#### Straight from Github 

To be able to use `sboxU` in your scripts, simply run the following command. It will download and compile `sboxU` (and also create some executable scripts, see above).

```
    sage --pip install git+https://github.com/lpp-crypto/sboxU
```

#### From a Local Copy of the Repository

Alternatively, especially if you want to work on `sboxU`, you can first clone this repository and then install `sboxU` from its content. To do this, simply `cd` to the directory containing this makefile and run

```
sage -pip install -e .
```


## Contributing

### How to

TODO

### Contributors

- [Jens Alich](https://informatik.rub.de/ac-personen/alich/)
- [Jules Baudrin](https://who.paris.inria.fr/Jules.Baudrin/)
- [Aurélien Boeuf](https://who.paris.inria.fr/Aurelien.Boeuf/)
- [Xavier Bonnetain](https://bonneta.in/)
- [Alain Couvreur](http://www.lix.polytechnique.fr/Labo/Alain.Couvreur/)
- [Mathias Joly](https://github.com/MathiasJoly)
- [Merlin Fruchon](https://who.paris.inria.fr/Merlin.Fruchon/)
- Pierre Galissant
- [Léo Perrin](https://who.paris.inria.fr/Leo.Perrin/)
- [Lukas Stennes](https://informatik.rub.de/symcrypt/personen/stennes/)





