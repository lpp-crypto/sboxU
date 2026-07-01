#!/usr/bin/env sh
# Invoke inside $sboxU_root/docs

# Dependencies : doxygen, graphviz, pandoc

# SAGE code documentation
# =======================

# Environnment

platform="$(uname)"

if [ "$platform" = "Darwin" ]; then
    viewer="open"
else
    viewer="xdg-open"
fi


# -- cleaning up
rm -f source/sboxU*.rst
rm -rf ./build/*

echo "building additional doc"
# grabbing tutorials/tests
find ../tests -name "*.md" | while read -r md; do
    base="$(basename "$md" .md)"
    pandoc "$md" -o "source/${base}.rst"
done


# building index.rst
pandoc ../README.md -o source/index.rst
sed "s/docs\/source\/logo-v2-5.png/logo-v2-5.png/g" -i source/index.rst
printf "\n.. toctree::\n   :maxdepth: 1\n   :caption: Tutorials\n\n" >> source/index.rst
find ../tests -name "*.md" | sort | while read -r md; do
    base="$(basename "$md" .md)"
    printf "   ./%s.rst\n" "$base" >> source/index.rst
done
echo "
.. toctree::
   :maxdepth: 2
   :caption: API documentation

   ./sboxU.rst
" >> source/index.rst



PYTHONPATH="$PYTHONPATH:.."
# -- parsing the docstrings - we choose to use sage python since the sphinx theme is installed via sage (cf pyproject.toml)
sage -python -m sphinx.cmd.apidoc -M -o source ../sboxU "../sboxU/scripts/apnDB/*QAM.py" "../sboxU/scripts/apnDB/BeierleLeander.py" "../sboxU/scripts/apnDB/apn_8bit_BLP22.py" 

# -- generating the bibliography
#sage ./biblio.py "gen"

# -- generating the HTML
sage -python -m sphinx.cmd.build -M html source build

# -- opening it
"$viewer" ./build/html/index.html


# C++ code documentation
# ======================

# -- parsing source code
doxygen ./cpp/doxygen.conf

# -- opening it
"$viewer" ./html/index.html
