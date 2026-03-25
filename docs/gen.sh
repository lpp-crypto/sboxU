#!/usr/bin/env sh


# SAGE code documentation
# =======================

# -- cleaning up
rm -f source/sboxU*.rst
rm -rf ./build/*

echo "building additional doc"
# grabbing tutorials/tests
pandoc ../tests/statistics/anomalies.md -o source/anomalies.rst


# building index.rst
pandoc ../README.md -o source/index.rst
sed "s/docs\/source\/logo-v2-5.png/logo-v2-5.png/g" -i source/index.rst
echo "
.. toctree::
   :maxdepth: 5
   :caption: Contents:

   ./modules.rst
   ./apn-study.rst
   ./bibliography.rst
   ./anomalies.rst
" >> source/index.rst 



PYTHONPATH="$PYTHONPATH:.."
# -- parsing the docstrings
sphinx-apidoc -M -o source ../sboxU "../sboxU/scripts/apnDB/*QAM.py" "../sboxU/scripts/apnDB/BeierleLeander.py" "../sboxU/scripts/apnDB/apn_8bit_BLP22.py" 

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
