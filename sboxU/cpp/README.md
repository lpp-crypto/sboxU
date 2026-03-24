# sboxU/C++


## What is this?

`sboxU` is ostensibly a Python/SAGE library, however, its speed comes thanks to the content of this directory: the `sboxU/C++` module.

Its content is compiled in a way that enables its call from Python when you run `sage -pip install .` in the root directory of this project, but it cannot be called from pure C++ programs. In order to do this, you must follow the steps below.

## Compiling

Create a directory called `build` in this folder (if it doesn't exist yet), and, from there, run

```sh
cmake ..
make
```

This will in particular create some executables corresponding to tests in the folder `build`. These have names of the form `test_<something>`. If you want to add tests (which is always a very good idea!), then add them in the `tests` folder, and update the end of the `CMakeLists.txt` to add a new executable target corresponding to your program.
