OPTGRA is an optimization algorithm developed and implemented by Johannes Schoenmaekers, it is specifically designed for near-linear constrained problems, which commonly occur in trajectory optimization.

This is a Modern Fortran refactoring. It is a work in progress.

### Notes

It is specifically designed for near-linear optimization problems with many constraints. When optimizing a problem, Optgra will first move towards satisfying the constraints, then move along the feasible region boundary to optimize the merit function, fixing constraint violations as they occur.

For this, constraints and the merit function are linearized. Optgra will perform less well on very non-linear merit functions or constraints.

### See also

  * [pyoptgra](https://github.com/esa/pyoptgra) -- Python interface to OPTGRA.
