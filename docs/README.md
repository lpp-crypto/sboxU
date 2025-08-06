Documentation
==============


In order to generate the documentation, you need `sphinx`, and the theme "readthedocs", which you can install with `pip`[^1]

Then, **assuming you have already installed `sboxUv2`**, run the following commands from this folder:

```
sphinx-apidoc -f -o source ../sboxUv2 -M
make html
```

The documentation can then be browsed starting from `./build/html/index.html`.


The long term idea is to use the github integration with `readthedocs` to automatically generate the documentation and store it there. Stay tuned.

[^1]: Careful if you have installed SAGE with conda: expect to fight with virtual environments. 
