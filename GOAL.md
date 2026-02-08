# Project Goal
Create optimized working versions of one- and two-electron integral code written in Rust.

## What exists already
The [pyquante](https://github.com/rpmuller/pyquante2) project has integrals written in pure python (https://github.com/rpmuller/pyquante2/tree/master/src/pyquante2/ints) and rewritten in C for speed (https://github.com/rpmuller/pyquante2/tree/master/src/pyquante2/cints)

There are typically different versions of the code.
- One electron integrals. These cover overlap, kinetic energy, and potential energy integrals. They are in the files https://github.com/rpmuller/pyquante2/blob/master/src/pyquante2/ints/one.py and https://github.com/rpmuller/pyquante2/blob/master/src/pyquante2/cints/cints.c.
- Electron repulsion integrals have three different versions, typically:
    - A "slow" version in https://github.com/rpmuller/pyquante2/blob/master/src/pyquante2/ints/two.py and https://github.com/rpmuller/pyquante2/blob/master/src/pyquante2/ints/one.py and https://github.com/rpmuller/pyquante2/blob/master/src/pyquante2/cints/cints.c
    - A version using Rys quadrature, in https://github.com/rpmuller/pyquante2/blob/master/src/pyquante2/ints/rys.py and https://github.com/rpmuller/pyquante2/blob/master/src/pyquante2/cints/crys.c
    - A version using the method of Head-Gordon and Pople, in https://github.com/rpmuller/pyquante2/blob/master/src/pyquante2/ints/hgp.py and https://github.com/rpmuller/pyquante2/blob/master/src/pyquante2/cints/chgp.c

## Testing
The Python code typically has doctests that show how the code is called and determine ground truth for the values.

## Boys integrals
- There is a working Rust implementation of the Boys integral, which is used in the Head-Gordon and Pople implementations, in the https://github.com/rpmuller/boys repository.
