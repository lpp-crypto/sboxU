Documentation
==============


In order to generate the documentation, you need `sphinx`, and the theme "readthedocs", which you can install with `pip`[^1]

Then, run the following commands from this folder:

```
sphinx-apidoc -f -o source ../sboxUv2 -M
make html
```

The documentation can then be browed starting from `./build/html/index.html`.

[^1]: Careful if you have installed SAGE with conda: expect to fight with virtual environments. 
