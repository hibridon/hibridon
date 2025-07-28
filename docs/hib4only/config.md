# CONFIG file

To install, compile and link the Hibridon code it is essential to have a `CONFIG` file which corresponds properly to your machine architecture, and to your choice
of compiler and math library.

The script `../bin/makefirst` attempts to create a `CONFIG` file appropriate to your machines architecture. As a reference, we have supplied with version 4.3.4 several `CONFIG` files, ennumerated below

We have confirmed that Hibridon v.4.3.4 installs, compiles, and runs correctly under the architectures listed below. If you have a slightly different architecture, or compiler/math library suite, then use the following `CONFIG` files as a guide.

## CONFIG files

| Hardware                     | OS       | compiler                                                  | math library                                                | CONFIG file                              |
| ---------------------------- | -------- | --------------------------------------------------------- | ----------------------------------------------------------- | ---------------------------------------- |
| Apple G5                     | OSX 10.3 | xlf 8.1                                                   | /Library/Frameworks/vecLib.framework/vecLib                 | CONFIG_APPLE_G4_G5                       |
| Apple Intel Core2duo or Xeon | OSX 10.5 | ifort 10.1.xxx 32 bit                                     | /Library/Frameworks/Intel_MKL.framework/Versions/10.0.2.018 | CONFIG_APPLE_INTEL_10.1.012              |
| Apple Intel Core2duo or Xeon | OSX 10.5 | ifort 10.1.xxx 64 bit                                     | /Library/Frameworks/Intel_MKL.framework/Versions/10.0.2.018 | CONFIG_APPLE_INTEL_10.1.012_64           |
| Apple Intel Core2duo or Xeon | OSX 10.5 | ifort 11.0.xxx 64 bit                                     | .../intel/Compiler/11.0/xxx/Frameworks/mkl                  | CONFIG_APPLE_INTEL_11.0.064              |
| Apple Intel Core2duo or Xeon | OSX 10.6 | ifort 11.1.xxx 64 bit                                     | .../intel/Compiler/11.1/xxx/Frameworks/mkl                  | CONFIG_APPLE_INTEL_11.1.076              |
| Intel X86 Xeon               | Linux    | /opt/intel/fce/10.1.xxx/bin/ifort[^b]                     | .../intel/mkl/10.0.xxx                                      | CONFIG_X86_INTEL_10.1.013_MKL_10.0.011   |
| Intel X86 Xeon               | Linux    | /opt/intel/3.2.2/Compiler/11.1/xxx/bin/intel64/ifortt[^b] | .../intel/3.2.2/mkl/10.2.2.025                              | CONFIG_X86_INTEL_11.1.059_MKL_10.2.2.025 |
| Intel X86 Xeon               | Linux    | .../pgi/linux86-64/7.1-3/7.1/bin/pgf95[^c]                | .../pgi/linux86-64/7.1-3/lib[^c]                            | CONFIG_X86_PGF95_ACML                    |

- [^a]: On Intel X86 platforms, to run jobs with total number of channels (`kmax` greater than 3941, see <a href=linking.html>linking</a>), you must change the `FI64` line in the `CONFIG` file to read `FI64="-mcmodel=medium -shared-intel"`
- [^b]: On your installation, the Intel compiler and mkl may be located elsewhere than /opt.  Check with your system administrator.
- [^c]: Check with your system administrator for the exact location of the pgi directory tree.


