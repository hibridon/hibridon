# Hibridon v5.0 alpha

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
* CMake >= 3.3
* Git (optional)
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

### 4. Create a configuration file for your potential energy surface

Create a new file with the `.user.cmake` extension at the root of the build directory (/tmp/hib_build) (e.g. `nh3h2.user.cmake`).

Paste the following content and edit to suit your needs:
```
# NH3-H2 sample cmake script

set(EXE_NAME "nh3h2.exe")
set(POT_SRC_FILE "/home/NH3H2/pot_nh3h2_2009.F90")
set(p_T_MATRIX_SIZE "500")
```
where 
* `EXE_NAME` will be the name of the hibridon executable
* `POT_SRC_FILE` is the full path to your potential fortran source code*
* `p_T_MATRIX_SIZE` is the T_MATRIX_SIZE parameter

\*Simply the filename (e.g. `pot_nh3h2_2009.F90` or `./pot_nh3h2_2009.F90`) if your potential is located in the build directory.


### 5. Configure hibridon's build

```bash
cd /tmp/hib_build
cmake /tmp/hib_src/hibridon
```
This will automatically find the required libraries and create a Makefile to build hibridon. 

#### to use mkl

make sure your environment variable `MKL_ROOT` is set to use the mkl lbrary you want. Depending on the system you are using, This can be achieved with one of the following methods:
- `source <intel_mkl_root_dir>/bin/mklvars.sh`
- `module load mkl/latest`

Once `MKL_ROOT` is set properly, you just have to tell cmake that you want to use mkl (see [https://cmake.org/cmake/help/v3.14/module/FindBLAS.html]):

```bash
cd /tmp/hib_build
cmake -DBLA_VENDOR=Intel10_64lp /tmp/hib_src
```

#### to configure the build with intel fortran...

make sure you have set the environment variables such that the `ifort` command works. Depending on the system you are using, This can be achieved with one of the following methods:
- `source <intel_mkl_root_dir>/bin/compilervars.sh`
- `module load mkl/latest`

Once `ifort` is in your path, you just have to tell cmake that you want to use ifort as the fortran compiler using the option `-DCMAKE_Fortran_COMPILER=ifort`, eg:

```bash
cd /tmp/hib_build
cmake -DCMAKE_Fortran_COMPILER=ifort /tmp/hib_src
```

### 6. Test hibridon (OPTIONAL)

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

### 7. Build hibridon

```bash
make <EXEC_NAME>
```
where `EXEC_NAME` is the executable name defined in your `.user.cmake` file.

### 8. one liner example

This one line command configures, builds and tests hibridon from a directory containing the source code.

```
graffy@graffy-ws2:/tmp/hibridon.build$ rm -R ./* ; script -q /dev/null --command  "cmake -DCMAKE_BUILD_TYPE=Debug /home/graffy/work/hibridon" && script -q /dev/null --command "make 2>&1" | tee "/home/graffy/work/hibridon/refactor_notes/make_$(date).stdout" && ctest
```

## For code contributors

### Code coverage

Code coverage option `ENABLE_CODE_COVERAGE` allows the delvelopers to identify the portions of hibridon source code that are not yet covered by the tests.

To activate code coverage, add `-DENABLE_CODE_COVERAGE=ON` to the cmake command. This option will generate code coverage info files when running tests.

Then, `make html_coverages`, will convert these coverage files into html reports:
- `<hibridon_build_dir>/coverage/<test_id>/index.html`: a report that shows the code covered by the test `<test_id>`
- `<hibridon_build_dir>/coverage/total/index.html`: a report that shows the code covered by all tests

### Performance profiling

To activate profiling, add `-DENABLE_PROFILING=ON` to the cmake command. This will build and run hibridon with profiling option. When run, each test will additionnaly create a `call_graph.pdf` file which shows where time was spent during the test.
