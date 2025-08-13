#!/usr/bin/env sh

# cleaning up
rm source/sboxUv2.*.rst
rm -rf ./build/*


# parsing the docstrings
sphinx-apidoc ../sboxUv2 -o source -M

# generating the HTML
make html
