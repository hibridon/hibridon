### Matrix-Matrix opérations

The major computational effort in calculations with the Hibridon code is associated with a series of matrix-matrix (BLAS3) operations (matrix inversions, diagonizations, and multiplications), which are carried out by calls to the [`LAPACK`](http://www.netlib.org/lapack) routines `DGETRF`, `DGETRI`, `DSYEVR`, and `DGEMM`.

The latest version of Intel's ([`MKL`](https://www.intel.com/content/www/us/en/develop/documentation/get-started-with-mkl-for-dpcpp/top.html) math library as well as of Apple's [vecLib](https://developer.apple.com/documentation/accelerate/veclib) include, automatically, parallel versions of these LAPACK routines, so that considerable performance enhancement can be achieved when running on multi-cpu, multi-core workstations.


At the end of every job, the [output](./com/out.md) file shows the total CPU and wall clock times for the calculation.

The wall time is the actual elapsed time for the job, while the CPU time is the time the job would take on a single cpu.

The ratio of *t*<sub>wall</sub> to *t*<sub>CPU</sub> is a measure of the use of parallelism.


### Hyperfine cross sections
Some of the loops in the [`HYPXSC`](./com/hypxsc.md) module of Hibridon's code are also explicitly parallel using OpenMP directives.

### Setting the maximum number of threads
The number of threads used for the computation can be tuned using the environment variable `OMP_NUM_THREADS`:
```bash
# Set maximum number of threads to 8:
export OMP_NUM_THREADS=8
```