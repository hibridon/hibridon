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

1. create a directory to store hibridon source
```bash
mkdir -p /tmp/hib_src
```
2. get hibridon source code
```bash
cd /tmp/hib_src
git clone https://github.com/hibridon/hibridon.git
```
This will create a directory /tmp/hib_src/hibridon, which is a clone of https://github.com/hibridon/hibridon.git 
3. create a directory to store hibridon's build
```bash
mkdir /tmp/hib_build
```

4. configure hibridon's build

```bash
cd /tmp/hib_build
cmake /tmp/hib_src
```
This will automatically find the required libraries and create a Makefile to build hibridon. 

5. build hibridon

```bash
make
```
6. test hibridon

The following command will run hibridon's unit tests.

```bash
make test
```
