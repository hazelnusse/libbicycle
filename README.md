## Bicycle Model
This project is a work in progress on bicycle modelling.  The primary goal is
to derive symbolic equations of motion for a bicycle model and output these
equations to C/C++ code that is can be used in numerical studies.  A secondary
goal is to experiment with code output techniques that may eventually be
included as part of sympy.physics.mechanics.

## Requirements

### Python code
* sympy
* numpy

### C++ code
* cmake 2.8 or later
* Eigen 3.1 or later

## Instructions

### Running the Python code

    $ cd derivation
    $ python derivation.py

### Building the C++ code

    $ mkdir build
    $ cd build
    $ cmake -D <PATH TO EIGEN3 FOLDER> ..
    $ make
