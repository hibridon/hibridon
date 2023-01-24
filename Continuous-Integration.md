Each time Hibridon source code is modified merging a Pull Request or directly pushing to the master branch, the new source code is automatically built and tested. This continuous integration ensures that the new modifications added to the source code does not unexpectedly affect the results of Hibridon calculations.

On the main page of the [Hibridon GitHub repository](https://github.com/hibridon/hibridon/), the two labels at the top of the build instructions allow you to check if the last bunch of tests were successful (passing):

<img width="400px" src="https://user-images.githubusercontent.com/20931653/214405225-1f61f2e1-af8a-4210-a2bb-37fb8e1c55a5.jpg">

 

Currently, these tests are done on two different platforms (macOS and Linux), using different build modes (Release and Debug), different compilers, and different BLAS/LAPACK libraries.
The full list of the different configurations used for the continuous integration is the following:

| Platform  | OS version  | Build mode  | Compiler  | BLAS/LAPACK Library |
|---|---|---|---|---|
| macOS  | 11.2  | Debug  | gfortran 9  |  Intel MKL |
| macOS  | 11.2  | Debug  | gfortran 10 |  Intel MKL |
| macOS  | 11.2  | Debug  | gfortran 11 |  Intel MKL |
| macOS  | 11.2  | Debug  | Intel ifort 2022  |  Intel MKL |
| macOS  | 11.2  | Release  | gfortran 9 |  Intel MKL |
| macOS  | 11.2  | Release  | gfortran 10 |  Intel MKL |
| macOS  | 11.2  | Release  | gfortran 11 |  Intel MKL |
| macOS  | 11.2  | Release  | Intel ifort 2022 |  Intel MKL |
|---|---|---|---|---|
| Linux  | Debian 9  | Debug | gfortran 6 |  Atlas |
| Linux  | Debian 9  | Debug | Intel ifort 2021 |  Atlas |
| Linux  | Debian 9  | Release | gfortran 6 |  Atlas |
| Linux  | Debian 9  | Release | Intel ifort 2021 |  Atlas |

