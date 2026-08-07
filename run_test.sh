#!/usr/bin/env sh

tangle_and_run() {
    if [ $1.md -nt $1.py ]; then # we tangle only when the md file is
                                 # more recent than the script
        echo "tangling "$1.md
        sboxU_test_tangle $1.md
    fi
    echo "running "$1.py
    sage $1.py
}

if sage setup.py build_ext --inplace -j 8 > /dev/null; then
    echo "sboxU compilation succeeded"
    if [ $# -eq 0 ]; then # case without any argument
        files=$(find tests/** -name "*.md" | sed "s/\.md$//")
    else # case where some tests are specified
        files="$@"
    fi
    for f in $files; do
        tangle_and_run $f
        if [ $? -ne 0 ]; then
            echo "Test $f failed, interrupting"
            return 1
        fi    
    done
    return 0
else
    echo "sboxU compilation failed"
    return 1
fi

