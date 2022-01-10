# Hibridon v5.0 alpha

[![CI_macpro6](https://github.com/hibridon/hibridon/actions/workflows/CI_macpro6.yml/badge.svg?branch=master)](https://github.com/hibridon/hibridon/actions/workflows/CI_macpro6.yml)

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

### 4. Create a configuration file for your potential energy surfaces

In this section we assume that you only have one hibridon potential energy surface code, which is located in `~/hibridon/my_pots/pot_nh3h2.F90`. Please note that a ready-to-use sample is available in https://github.com/hibridon/hibridon.git/doc/samples/arno.

Create a `CMake` project file `~/hibridon/my_pots/CMakeLists.txt` that will describe the way to build `nh3h2.exe` from `~/hibridon/my_pots/pot_nh3h2.F90` and hibridon source code. Here's an example of `~/hibridon/my_pots/CMakeLists.txt`'s contents:

```cmake
# this is a minimal cmake example for creating a hibridon executable with a user-defined potential
cmake_minimum_required (VERSION 3.3)
project (my_pots)
enable_language (Fortran)

# add hibridon library
add_subdirectory("/tmp/hib_src" hibridon)

# declare a new executable using hibridon's add_hibexe cmake function, where:
# - the 1st argument (here nh3h2.exe) is the name of the resulting executable
# - the 2nd argument (here "pot_nh3h2.F90") is the file path of the user provided potential file
# - the 3rd argument (here "kmax") is the size of the t matrix
#    - kmax : for normal cases
#    - kbig : for special cases (only arn2_big test uses it)
add_hibexe(nh3h2.exe "pot_nh3h2.F90" "kmax")
```

Please note that this project `my_pots` only describes how to build one hibridon executable (`nh3h2.exe`). However, more potential executables can be added in the same project: to do that, just add more `add_hibexe` statements.

Also please note that this project includes `hibridon` through the use of `add_subdirectory` statement. This means that building your project will build `hibridon`.

### 5. Configure your project's build

```bash
mkdir -p /tmp/my_pots_build
cd /tmp/my_pots_build
cmake ~/hibridon/my_pots
```
This will automatically find the required libraries and create a Makefile to build hibridon. 

#### to use mkl

make sure your environment variable `MKL_ROOT` is set to use the mkl lbrary you want. Depending on the system you are using, This can be achieved with one of the following methods:
- `source <intel_mkl_root_dir>/bin/mklvars.sh`
- `module load mkl/latest`

Once `MKL_ROOT` is set properly, you just have to tell cmake that you want to use mkl (see [https://cmake.org/cmake/help/v3.14/module/FindBLAS.html]):

```bash
cd /tmp/my_pots_build
cmake -DBLA_VENDOR=Intel10_64lp ~/hibridon/my_pots
```

#### to configure the build with intel fortran...

make sure you have set the environment variables such that the `ifort` command works. Depending on the system you are using, This can be achieved with one of the following methods:
- `source <intel_mkl_root_dir>/bin/compilervars.sh`
- `module load mkl/latest`

Once `ifort` is in your path, you just have to tell cmake that you want to use ifort as the fortran compiler using the option `-DCMAKE_Fortran_COMPILER=ifort`, eg:

```bash
cd /tmp/my_pots_build
cmake -DCMAKE_Fortran_COMPILER=ifort ~/hibridon/my_pots
```

### 6. Test hibridon (OPTIONAL)

The following command will run all hibridon's tests:

```bash
cd /tmp/my_pots_build
make test
```

You can also run a test suite (a group of tests). For example, the following command will run the test suite `short`:

```bash
cd /tmp/my_pots_build
make testsuite_short
```

You can also run a single test. For example, the following command will run the test `arn2`:

```bash
cd /tmp/my_pots_build
ctest -L '^arn2$'
```

### 7. Build your potential executables

```bash
cd /tmp/my_pots_build
make nh3h2.exe
```
This should create the file `/tmp/my_pots_build/nh3h2.exe`

### 8. one liner example

This one line command configures, builds and tests hibridon from a directory `/home/graffy/work/hibridon` containing the source code (please note that this only builds and tests hibridon library; it doesn't build any user-provided potential energy surface).

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
