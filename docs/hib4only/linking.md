# Linking an Executable (the `makehib` command) (for Hibridon < 5.0)

To link an executable module you must execute the command

```sh
makehib [-b -d -c] {potname} [kmax]
```

where
- `{potname}` is the name of the [potential](../potlist.md) subroutine `pot_{potname}.f`, which
describes the system under investigation and which must be located in the Hibridon subdirectory `src/pot`
- `kmax` is the numerical value of the maximum number of channels. If this field is omitted, a default value of `kmax = 151` is used.
  - ⚠️  On Intel X86 platforms, to link and run jobs with kmax greater than 3941, you must change the `FI64` line in the [`CONFIG`](./config.md) file to read:
    > `FI64="-mcmodel=medium -shared-intel"`
- `-b` (optional) forces creation of a code requires fewer matrices and can thus handle larger number of channels (for more details see [Memory Requirements](../memory.md))
- `-d` (optional) forces creation of a code in which the [potential](../potlist.md) subroutine `pot_{potname}.f` is compiled with the debug (f77 -g) option.
- `-c` (optional) forces creation of a code in which direct access I/O is handled by C rather than FORTRAN routines.

the `makehib` command creates an executable module `bin/progs/hib_{potname}_kmax`.

If the `-b` or `-d` `-f` options (or both) are present, the 
executable module is named, respectively,
- `hibb_{potname}_kmax`
- `hibd_{potname}_kmax`
- `hibc_{potname}_kmax`
- `hibbd_{potname}_kmax`
- and so forth

The  file `h_{potname}.log` contains a log of this link, which may be examined in case of failure.

For example,

```sh
makehib arno 251
```

will create an executable module `bin/progs/hib_arno_251`

Here is the terminal output from this command:

```
potname is arn2
current directory is /Users/mha/hib43/bin/progs
compiling pot subroutine for hib_arn2_151 ...
potname is pot_arn2.f

making fortran converter ...

running fortran converter ...

compiling pot_arn2.f ...
pot_arn2 successfully compiled ...

hib_arn2_151 being created at date 12/26/07 07:30:01 EST with version 4.3

building  Hibridon v.4.3 with Fortran direct i/o routines ...

system parameters are:
System is OS X 10.4.11
Machine type is unix-darwin unix-ifort
C compiler is /usr/bin/cc -DNAME_LU
Fortran compiler is /opt/intel/fc/10.1.007/bin/ifort -O3 -save
Mathematical Libraries are:
-L/Library/Frameworks/Intel_MKL.framework/Versions/10.0.1.014/Libraries/32 -lmkl_intel -lmkl_intel_thread -lmkl_core -lmkl_intel_thread -lmkl_core -lmkl_intel_thread -lmkl_core -lguide -lpthread
compiling hiversion ...

converting himain with kmax = 151

linking hib_arn2_151
  c-based direct i/o routines used

removing hiiolib*.o and (if present) hiunix.o from directory src ...

successful link; moving hib_arn2_151 to ../bin/progs ...

examine /Users/mha/hib43/h_arn2.log for diagnostic and/or error messages
```

⚠️ Your link will be unsuccessful unless all the requisite `.o` and `common` files have been previously created.  See the [installation](./install.md)
instructions.
