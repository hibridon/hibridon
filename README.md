# Hibridon v4.4 beta

[![Build Status](https://jenkins.ipr.univ-rennes1.fr/buildStatus/icon?job=hibridon-build)](https://jenkins.ipr.univ-rennes1.fr/job/hibridon-build/)

Computer Software for
Molecular Inelastic Scattering and Photodissociation

Documentation is available on
- [http://www2.chem.umd.edu/groups/alexander/hibridon]
- [<hibridon_root_path>/doc](doc)

Changes:
- [ReleaseNotes](ReleaseNotes)
- [README_4.2.1](README_4.2.1) : changes from `hibridon` 4.1.5 to 4.2.1
- [README_4.3](README_4.3) : changes from `hibridon` 4.2.1 to 4.3 
- [README_4.4](README_4.4) : changes since `hibridon` 4.3

## Build instructions

Required tools:
* CMake >= 2.6
* Git
* Fortran compiler with support for fortran 90 and fpp preprocessing (eg: gfortran or ifort).

Required libraries:

* Blas compatible library (eg Intel's Math Kernel Library)
* Lapack compatible library (eg Intel's Math Kernel Library)

### 1. Create a directory to store hibridon source

```bash
mkdir -p /tmp/hib_src
```
### 2. Get hibridon source code

```bash
cd /tmp/hib_src
git clone https://github.com/hibridon/hibridon.git
```
This will create a directory /tmp/hib_src/hibridon, which is a clone of https://github.com/hibridon/hibridon.git 
### 3. Create a directory to store hibridon's build

```bash
mkdir /tmp/hib_build
```
### 4. Configure hibridon's build

```bash
cd /tmp/hib_build
cmake /tmp/hib_src
```
This will automatically find the required libraries and create a Makefile to build hibridon. 

#### to use mkl

make sure your environment variable `MKL_ROOT` is set to use the mkl lbrary you want. Depending on the system you are using, This can be achieved with one of the following methods:
- `source
- `module load mkl/latest`

Once `MKL_ROOT` is set properly, you just have to tell cmake that you want to use mkl (see [https://cmake.org/cmake/help/v3.14/module/FindBLAS.html]):

```bash
cd /tmp/hib_build
cmake -DBLA_VENDOR=Intel10_64lp /tmp/hib_src
```


### 5. Build hibridon

```bash
make
```
### 6. Test hibridon

The following command will run all hibridon's tests:

```bash
make test
```

You can also run a test suite (a group of tests). For example, the following command will run the test suite `short`:

```bash
make testsuite_short
```

You can also run a single test. For example, the following command will run the test `arn2`:

```bash
ctest -L '^arn2$'
```

### 7. one liner example

This one line command configures, builds and tests hibridon from a directory containing the source code.

```
graffy@graffy-ws2:/tmp/hibridon.build$ rm -R ./* ; script -q /dev/null --command  "cmake -DCMAKE_BUILD_TYPE=Debug /home/graffy/work/hibridon" && script -q /dev/null --command "make 2>&1" | tee "/home/graffy/work/hibridon/refactor_notes/make_$(date).stdout" && ctest
```
