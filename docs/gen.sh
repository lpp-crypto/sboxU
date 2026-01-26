#!/usr/bin/env sh


# SAGE code documentation
# =======================

# -- cleaning up
rm source/sboxUv2*.rst
rm -rf ./build/*


PYTHONPATH="$PYTHONPATH:.."
# -- parsing the docstrings
sphinx-apidoc ../sboxUv2 -o source -M

# -- generating the bibliography
#sage ./biblio.py "gen"

# -- generating the HTML
sphinx-build -M html source build

# -- opening it
xdg-open ./build/html/index.html 


# C++ code documentation
# ======================

# -- parsing source code
doxygen ./cpp/doxygen.conf

# -- opening it
xdg-open ./html/index.html
