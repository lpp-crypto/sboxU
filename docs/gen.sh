#!/usr/bin/env sh

# cleaning up
rm source/sboxUv2*.rst
rm -rf ./build/*


# parsing the docstrings
sphinx-apidoc ../sboxUv2 -o source -M

# generating the bibliography
sage ./biblio.py "gen"

# generating the HTML
sphinx-build -M html source build

# opening it
open ./build/html/index.html 
