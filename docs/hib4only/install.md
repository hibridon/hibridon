# Downloading and Installing the Hibridon Package (Hibridon < 5.0 only)

Version 4.3 of the Hibridon ([copyright](copyright.md)) package
is now available for the following machines:
- Apple G5 workstations running under OSX 10.3 (Panther) or OSX 10.4 (Tiger)
  - Fortran Compiler: [xlf 8.1](http://www.absoft.com/xlf/index.htm)
- Apple Intel workstations 
  - OSX 10.5 (Leopard); Fortran Compiler: [v. 10.1, and v. 11.0](http://www.intel.com/cd/software/products/asmo-na/eng/compilers/267426.htm)
  - OSX 10.6 (Snow Leopard); Fortran Compiler: [v. 10.1, v. 11.1 and v. 12](http://www.intel.com/cd/software/products/asmo-na/eng/compilers/267426.htm) 

- Xeon x86 workstations running under Linux (Redhat 4 has been tested)
  - Fortran Compilers:  [Intel ifort v. 10.0](http://www.intel.com/support/performancetools/fortran/linux/)  and v 11.1, as well as [Portland Group pgf95 v.7.1](http://www.pgroup.com/products/workpghpf.htm)

- SGI Altix (Itanium ia64) running under Linux (Suse has been tested)
  - Fortran Compilers: [Portland Group pgf95 v.7.1](http://www.pgroup.com/products/workpghpf.htm)
  
<!--
- IBM RS/6000 workstations running under AIX 4.3
  - Fortran Compiler:  xlf 7.1
- HP 9000 series workstations running under HPUX 11.0
-->

Only the Apple (10.5 and 10.6) and Linux x86 versions will be actively maintained.

You can obtain a copy of Hibridon 4.3 as a `g'zipped tar` [archive](https://www2.chem.umd.edu/groups/alexander/hibridon/hib43/hib43.tar.gz) in binary format (~1 MB)

## Prerequisites

If you are installing on a platform that uses the Intel fortan compiler, make sure that the compiler environment
is properly set up. 

- For ifort v. 10, this is explained in the directory `...intel/fce/{version}/doc/INSTALL.htm#settingup`.
Specifically, the scripts 
  - `...intel/fce/{version}/bin/ifortvars.sh(.csh)`
  - `...intel/fce/{version}/bin/idbvars.sh(.csh)`

should be executed. The "fce" directories contain the 64 bit versions of the compiler. As can be seen in [timing](./timing.md), the 64 bit code executed ~ 1.5 times faster than the 32 bit code (invoked by sourceing the "fc" rather than "fce" files).

- For ifort v. 11 (testing has been done with v. 11.0.064, 11.1.058,  11.1.076, 11.1.089) this is explained in the file  `/opt/intel/Compiler/11.1/xxx/Documentation/en_US/documentation_f.htm`.<br>
Specifically, the scripts
  - ` source /opt/intel/Compiler/11.1/xxx/bin/ifortvars.csh intel64`
  - `source /opt/intel/Compiler/11.1/xxx/bin/ifortvars.sh intel64`

should be executed (by source'ing these commands) in your `.profile (.login)` file.

- For ifort v. 12 (testing has been done with v. 12.0.0) this is explained in the file  `/opt/intel/composerxe-2011.0.085/Documentation/en_US/documentation_f.htm`. Specifically, the scripts
  - ` source /opt/intel/bin/compilervars.csh intel64`
  - ` source /opt/intel/bin/ifortvars.csh intel64`
  - ` source /opt/intel/bin/compilervars.sh intel64`
  - ` source /opt/intel/bin/ifortvars.sh intel64`
should be executed (by source'ing these commands) in your `.profile (.login)` file.

- ⚠️ For Apple OSX 10.5 (Leopard), Hibridon works with version 10.1 of the Intel Compiler under both 32 bit (` .../intel/fc`) and 64 bit (` .../intel/fce`) modes, as well as with version 11.0 under 64 bit mode.

- ⚠️ For Apple OSX 10.6 (Snow Leopard), Hibridon works with version 10.1 of the Intel Compiler under  32 bit (` .../intel/fc`) mode as well as with version 11.1 and 12 under 64 bit mode.

- You will need to have a local BLAS and LAPACK library
  - If you are planning to use Intel's MKL library, make sure that the required environment variables are properly set, 
as explained in the file `userguide.pdf`.  
    - On Apple machines and ifort version 10 or 11 this is located in `/Library/Frameworks/Intel_MKL.framework/Documentation`.  If you are using the 64 bit compiler, you should use the 64 bit (`em64t`) MKL libraries, rather than the 32 bit libraries.
    - On Apple machines and ifort version 12 this is located in `/opt/intel/composerxe-2011.0.085/Documentation/en_US/documentation_f.htm`

    - On Intel Linux platforms, you may have to ask your system administrator where this file is located.
  - On Apple Intel platforms running ifort 10.x or 11.x, the command to set the required MKL environment variables is
    - `source /opt/intel/Compiler/1X.y/zzz/Frameworks/mkl/tools/environment/mklvarsem64t.csh`
    - `source /opt/intel/Compiler/11.0/064/Frameworks/mkl/tools/environment/mklvars32.csh` (for 32 bit)
  - On Apple Intel platforms running ifort 12.x, the command to set the required MKL environment variables is
```
source /opt/intel/mkl/bin/intel64/mklvars_intel64.csh
```

    - If you are using the Portland Group compiler and planning to use the ACML Library, you may have to ask your system administrator where the library is located.

## End of prerequisites

After you have obtained the g'zipped archive, follow these steps:

1. Create a main hibridon directory. We shall designate this here as `hib43x`, although you can give it any name you want. 
2. Move the g'zipped tar archive into this directory
3. Uncompress and unpack the tar archive
```
gunzip hib4.3.tar.gz |tar xvf -
```
4. The Hibridon hypertext help facility (a copy of the help facility you are reading here) is contained  in the directory `hib43x/doc/hib_html`. The initial target file for your web browser should be `hib43x/doc/hib_html/index.html`
5. Extend your PATH environmental variable to include the `bin` subdirectory of your main hibridon directory, that is `hib43x/bin`  
6. ⚠️ Read the [license](./license.md) agreement.  If you **do not accept** all provisions of this agreement, **stop here** and destroy the software. If you **accept** the license agreement, please read [the instructions for acknowledgement](../Home.md) of the code in any future publications from your group that use the Hibridon package.
7. `cd` to the subdirectory `hib43x/bin`.
8. Execute the command
```
makefirst
```
This creates a number of local configuration files.  

  - You will be asked if you accept the license. Respond "yes" if you accept, otherwise destroy the software. 
  - You will also be asked to supply
    - The name of your FORTRAN compiler (the install script will suggest a default name, which you can change)
    - The name of your C compiler (the install script will suggest a default name, which you can change)
    - The path of your BLAS and LAPACK libraries (the install script will suggest a default name, which you can change)

The script `makefirst` will create the file `CONFIG` in your top-level Hibridon directory. You should
inspect this file to make sure you agree with the choice of compilers, compiler options, and BLAS/LAPACK libraries. Several default `CONFIG` files are supplied for reference, or in case the `makefirst` script does not execute properly at your installation.

The `makefirst` script will create a `CONFIG` file suitable for either a 64-bit or 32-bit.  This choice will depend on what your default compiler is (which you can verify by the command `which ifort`).

9. Execute the command

```
hib_makeobj
```
This compiles all the hibridon source code and creates the on-line executable `hib43x/bin/hib_help`.
The file `hib43x/obj.log` contains a log of the compilation, in case of failure. Check this log to search for any compiler errors.

To obtain help at any time you can either
  - Execute the command `hib_help` and follow instructions.
  - Use a Web Browser to read the Hibridon HTML [help library](http://www.chem.umd.edu/groups/alexander/hibridon/hibhelp.html).

10. Test your code by executing the command
```
hibtest
```

which runs a series of [test](../tests.md) calculations. After the tests are run, the results are saved
in the directory `hib43x/testnew`. All files in this directory should be compared with the canonical results, which are stored in the directory `hib43x/tests`.  
  - On Apple platforms this comparison can be done graphically, using Apple's `FileMerge` application, by executing the command
  
```
compare_OSX
```
  - On non-Apple platforms, this comparison can be done by excuting the command 
```
compare.com
```

The comparison between the files in `../tests` and those in `../testnew` should agree to within at least 6 significant figures.

11. You are now ready to [link](./linking.md) and [run](../com/run.md) your own code for your particular scattering problem.  Use one of the [potential](../potlist.md) subroutines in the directory `hib43x/src/pot` as a guide for writing your own. You will certainly find it helpful to consult the [Basis](../basis.md) help files.

12. You can run several timing benchmarks by the command `hib_timer` and compare with the results obtained on other platforms (as discussed in more detail on the [timing](./timing.md) page).

13. *Good luck!*

---

[📧](mailto:mha@umd.edu) Please report all bugs that you may find. They will (hopefully) be eliminated in a timely manner.
