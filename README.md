# sboxuV2

This branch is used for the development of `sboxUv2`, the 2.0 version of `sboxU` that is essentially a full re-write intended to:
1. Allow easier maintenance. This is a big one as more and more people work not just *with* but *on* `sboxU`. The current code-base does not really allow that. In particular, tests will be provided, and the structure of the code itself will be clearer (more modules, each with submodules if necessary).l
2. Support S-boxes with different input and output sizes.
3. Make sure the `C++` code can talk to itself, instead of consisting of a myriad of independent functions that are just called from `python`. This will allow a cleanup of said `C++` code, and enable moving larger and larger parts of the logic to pure `C++`.


To install it, you need to have the code of this branch somewhere in your file system. Then, simply run:

```
sage -pip install .
```

Then, you will be able to `import sboxUv2` from other sage scripts!
