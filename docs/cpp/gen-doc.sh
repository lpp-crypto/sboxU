#!/usr/bin/env sh


# In order to generate the documentation of the C++ API, simply run
# this script in the folder in which it lives. It requires `doxygen`,
# which is for instance available in the repositories of ubuntu.

doxygen ./doxygen.conf
open ./html/index.html
