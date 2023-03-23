
[![Full CI on macOS-11.2](https://github.com/hibridon/hibridon/actions/workflows/full_macOS-11.2.yml/badge.svg)](https://github.com/hibridon/hibridon/actions/workflows/full_macOS-11.2.yml)
[![Long CI on macOS-11.2](https://github.com/hibridon/hibridon/actions/workflows/long_macOS-11.2.yml/badge.svg)](https://github.com/hibridon/hibridon/actions/workflows/long_macOS-11.2.yml)

[![Full CI on Debian 9](https://github.com/hibridon/hibridon/actions/workflows/full_Debian-9.yml/badge.svg)](https://github.com/hibridon/hibridon/actions/workflows/full_Debian-9.yml)
[![Long CI on Debian 9](https://github.com/hibridon/hibridon/actions/workflows/long_Debian-9.yml/badge.svg)](https://github.com/hibridon/hibridon/actions/workflows/long_Debian-9.yml)

[What does that mean ?](https://github.com/hibridon/hibridon/wiki/Continuous-Integration)

---
# Hibridon v5.0

Hibridon is a program package to solve the close-coupled equations which occur in the quantum treatment of inelastic atomic and molecular collisions. Gas-phase scattering, photodissociation, collisions of atoms and/or molecules with flat surfaces, and bound states of weakly-bound complexes can be treated

The full documentation is available on [https://github.com/hibridon/hibridon/wiki](https://github.com/hibridon/hibridon/wiki)

Changes:
- [CHANGELOG](CHANGELOG.md)


#  Prerequisites

## Required tools:
* [CMake](https://cmake.org/install/) >= 3.3
* Fortran compiler with support for fortran 2008 and fpp preprocessing (e.g. [GNU Fortran](https://fortran-lang.org/learn/os_setup/install_gfortran)).
* [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) (optional)


## Required libraries:

* BLAS compatible library (e.g. [Netlib BLAS](http://www.netlib.org/blas/))
* LAPACK compatible library (e.g. [Netlib LAPACK](http://www.netlib.org/lapack/))
  

# Build instructions

## 1. Get Hibridon source code
### Using GIT

```bash
git clone https://github.com/hibridon/hibridon.git ~/hib_src
```
This will create a directory `~/hib_src/` which is a clone of https://github.com/hibridon/hibridon.git 

### Without GIT
- Download Hibridon source code as a zip archive: [hibridon-master.zip](https://github.com/hibridon/hibridon/archive/refs/heads/master.zip)

- Extract the content of the archive into the directory of your choice, e.g.: `~/hib_src/` 

## 2. Create a configuration file for your potential energy surfaces (PESs)

Create the directories that will contain Hibridon build files and store your project configuration file:
```bash
mkdir -p ~/hib_build/project/
```

Create a `CMake` project configuration file:
```bash
touch ~/hib_build/project/CMakeLists.txt
```

Copy and past the following example into your `CMake` project configuration file:

```cmake
# this is a minimal cmake example for creating a hibridon executable with an user-supplied PES
cmake_minimum_required (VERSION 3.3)
project (my_pots)
enable_language (Fortran)

# add hibridon library by specifying the directory in which the source files are
add_subdirectory("~/hib_src/" hibridon)

# declare new executables using hibridon's add_hibexe cmake function, where:
#  - the 1st argument is the name of the resulting executable
#  - the 2nd argument is the file path of the user provided potential file
#  - the 3rd argument is the size of the t matrix:
#    - "kmax": for normal cases
#    - "kbig": for special cases (only arn2_big test uses it)

add_hibexe(NH3-H2.exe "~/my_pots/pot_nh3h2.F90" "kmax") # NH3-H2

# You can add as many executables as you want by using the add_hibexe function:
#add_hibexe(OH-H2.exe "~/my_pots/pot_ohh2.F90" "kmax") # OH-H2
#add_hibexe(He-CO.exe "~/my_pots/pot_heco.F90" "kmax") # He-CO
```

Adapt this example to suit your needs.



## 3. Configure your project's build

```bash
cd ~/hib_build
cmake ./project/
```
This will automatically find the required libraries and compiler and create a Makefile to build Hibridon. 


### Useful CMake options and special cases:



- **Use a specific compiler**
- 
    ```bash
    cd ~/hib_build
    cmake ./project/ -DCMAKE_Fortran_COMPILER=<compiler>
    ```
    Where `<compiler>` is your compiler executable e.g.:
    - `gfortran`
    - `gfortran-8`
    - `ifort` 

- **Use a specific BLAS/LAPACK library**

    ```bash
    cd ~/hib_build
    cmake ./project/ -DBLA_VENDOR=<BLAS_LIB> 
    ```
    Where `<BLAS_LIB>` is your BLAS library e.g.:
    - `OpenBLAS`
    - `Intel10_64lp`
    - `Apple`
  
    see [CMake BLAS/LAPACK VENDORS](https://cmake.org/cmake/help/latest/module/FindBLAS.html#blas-lapack-vendors) for a full list of supported libraries


- **Enable Hibridon testing**
 
    ```bash
    cd ~/hib_build
    cmake ./project/ -DBUILD_TESTING=ON
    ```
- **Build in debug mode**

  By default, the source files are compiled in "Release" mode.
  You can build Hibridon in a "Debug" mode that will disable all compiler optimizations and produce debugging symbols table and traceback informations:
    ```bash
    cd ~/hib_build
    cmake ./project/ -DCMAKE_BUILD_TYPE=Debug
    ```

- **Change the default install prefix**
  By default, the executables will be installed in `/usr/local/bin/`.
  You can change this by setting the `CMAKE_INSTALL_PREFIX` variable:
    ```bash
    cd ~/hib_build
    cmake ./project/ -DCMAKE_INSTALL_PREFIX=<path>
    ```
    Where `<path>` is the path to the directory where you want to install the executables.

## 4. Make your project
- **Make all executables**
    ```bash
    make
    ```
  The executable files will be put in the current directory (`~/hib_build`).
- **Make a specific executable** 
  ```bash
    make <executable>
  ```
  Where `<executable>` is one of the executable you defined in the CMakeLists.txt configuration file, e.g.: `NH3-H2.exe`.

- **Install executables in /usr/local/bin/ (Optional)** 

  You can install all the generatted executables in the `<prefix>/bin/` directory by running:
    ```bash
    make install
    ```
  By default, `<prefix>` is set to `/usr/local/`, which is the standard location for user-installed software on macOS and most Linux distributions.

  Note that this will overwrite any existing executable with the same name in the destination folder.

  The default <prefix> directory (`/usr/local/`) can be changed by setting the `CMAKE_INSTALL_PREFIX` (see previous section).

## 5. Test Hibridon (Optional)
Hibridon testing must be enabled (see Section 3: Configure your project's build).
The following commands need to be executed within the Hibridon build directory (e.g. ~/hib_build/project/hibridon/).

- **Run all tests**
    ```bash
    ctest
    ```
- **Run a test suite (group of tests)** 
    ```bash
    ctest -L <testsuite>
    ```
    Where testsuite is the name of the test suite among:
    - `coverage` (covers most of the source code)
    - `quick` (runs only quick tests)
    - `benchmark` (runs long tests for the purpose of measuring `hibridon`'s performance)
    
- **Run tests for only one PES** 
    ```bash
    ctest -L <label>
    ```

    Where `<label>` is the name of a group of tests associated to one PES.
    
    The full list of avaible labels can be printed using: 
    ```bash
    ctest --print-labels
    ```

### 6. Run your project
Once your executables are created (step 4), you can interactively run hibridon using:
```bash
./<executable> -k <kmax>
# or 
./<executable> --kmax <kmax>
```
Where `<kmax>` is the maximum number of channels (see [Memory Requirements](https://github.com/hibridon/hibridon/wiki/Memory-requirements)).


You can also run your Hibridon executable using a command file:
```bash
./<executable> -k <kmax> < <command file>
# or
./<executable> -k <kmax> -c <command file>
# or
./<executable> -k <kmax> --com <command file>
```
Where `<commands file>` is a file containing the input commands you want to execute (see [Batch-mode-and-background-execution](https://github.com/hibridon/hibridon/wiki/Batch-mode-and-background-execution))

### 7. One liner example

This one line command configures, builds and tests Hibridon from a directory `~/hib_build` containing the source code.

```bash
 git clone https://github.com/hibridon/hibridon.git ~/hib_src && mkdir -p ~/hib_build && cmake -DCMAKE_Fortran_COMPILER=gfortran -DBUILD_TESTING=ON -S ~/hib_src/ -B ~/hib_build && cd ~/hib_build/ && make && ctest -L coverage
```

Please note that this only builds and tests Hibridon library; it doesn't build any user-provided PES.

## For code contributors

<!---
### Code coverage

Code coverage option `ENABLE_CODE_COVERAGE` allows the delvelopers to identify the portions of hibridon source code that are not yet covered by the tests.

To activate code coverage, add `-DENABLE_CODE_COVERAGE=ON` to the cmake command. This option will generate code coverage info files when running tests.

Then, `make html_coverages`, will convert these coverage files into html reports:
- `<hibridon_build_dir>/coverage/<test_id>/index.html`: a report that shows the code covered by the test `<test_id>`
- `<hibridon_build_dir>/coverage/total/index.html`: a report that shows the code covered by all tests
-->

### Performance profiling

To activate profiling, add `-DENABLE_PROFILING=ON` to the cmake command. This will build and run hibridon with profiling option. When run, each test will additionnaly create a `<test_id>_call_graph.pdf` file which shows where time was spent during the test.
