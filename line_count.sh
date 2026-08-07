#!/usr/bin/env sh


for extension in cpp hpp pyx pxd py; do
    echo $extension
    echo "=========\n"
    wc -l $(find . -name "*.$extension" | grep -v build | grep -v apnDB |grep -v "cython_functions.cpp")
    echo "\n\n"
done
    
