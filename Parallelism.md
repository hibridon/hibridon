The major computational effort in calculations with the Hibridon code is associated with a series of matrix-matrix (BLAS3) operations (matrix inversions, diagonizations, and multiplications), which are carried out by calls to the [LAPACK](/www.netlib.org/lapack) routines `DGETRF`, `DGETRI`, `DSYEVR`, and `DGEMM`.

The latest version of Intel's ([MKL](https://www.intel.com/content/www/us/en/develop/documentation/get-started-with-mkl-for-dpcpp/top.html) math library as well as of Apple's [vecLib](https://developer.apple.com/documentation/accelerate/veclib) include, automatically, parallel versions of these LAPACK routines, so that considerable performance enhancement can be achieved when running on multi-cpu, multi-core workstations.


At the end of every job, the [output](output) file shows the total CPU and wall clock times for the calculation.

The wall time is the actual elapsed time for the job, while the CPU time is the time the job would take on a single cpu.

The ratio of *t*<sub>wall</sub> to *t*<sub>CPU</sub> is a measure of the use of parallelism.