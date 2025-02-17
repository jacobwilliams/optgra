OPTGRA is an optimization algorithm developed and implemented by Johannes Schoenmaekers, it is specifically designed for near-linear constrained problems, which commonly occur in trajectory optimization.

[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)
[![CI Status](https://github.com/jacobwilliams/optgra/actions/workflows/CI.yml/badge.svg)](https://github.com/jacobwilliams/optgra/actions)
[![last-commit](https://img.shields.io/github/last-commit/jacobwilliams/optgra)](https://github.com/jacobwilliams/optgra/commits/master)
<!-- [![GitHub release](https://img.shields.io/github/release/jacobwilliams/optgra.svg)](https://github.com/jacobwilliams/optgra/releases/latest) -->
<!-- [![codecov](https://codecov.io/gh/jacobwilliams/optgra/branch/master/graph/badge.svg)](https://codecov.io/gh/jacobwilliams/optgra) -->



This is a Modern Fortran refactoring. It is a work in progress.

### Notes

It is specifically designed for near-linear optimization problems with many constraints. When optimizing a problem, Optgra will first move towards satisfying the constraints, then move along the feasible region boundary to optimize the merit function, fixing constraint violations as they occur.

For this, constraints and the merit function are linearized. Optgra will perform less well on very non-linear merit functions or constraints.

### Documentation

The latest API documentation for the `master` branch can be found [here](https://jacobwilliams.github.io/optgra/). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford).


### See also

  * [pyoptgra](https://github.com/esa/pyoptgra) -- Python interface to OPTGRA.
