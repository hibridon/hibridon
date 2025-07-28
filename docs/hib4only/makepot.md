
# Testing your potential subroutine (the `makepot` command) (only for Hibridon < 5.0)

To test your [potential subroutine](../potlist.md) you must execute the command

```
makepot {potname}
```

where
- `{potname}` is the name of the [potential](../potlist.md) subroutine `pot_{potname}.f`, which
describes the system under investigation and which must be located in the Hibridon subdirectory `src/pot`

the `makepot` command creates an executable module `bin/progs/runpot_{potname}_kmax`.

The  file `p_{potname}.log` contains a log of this link, which may be examined in case of failure.

---

For example,

```
makepot arno_ccsdt
```

will create an executable module `bin/progs/runpot_arno_ccsdt`

Here is the terminal output from this command

```
/Users/mha/hib43/bin/progs> makepot arno_ccsdt

making fortran converter ...

running fortran converter ...

compiling runpot.f ...

linking runpot_arno_ccsdt

successful link; moving runpot_arno_ccsdt to ../bin/progs ...

examine /Users/mha/hib43/p_arno_ccsdt.log for diagnostic and/or error messages
```

Running the resulting executable code `runpot_{potname}` will allow you to print out the expansion coefficients
of the potential at various values of $R$, as controlled in the subroutine "driver" in the file `../src/pot/{potname}.f`.

For example, for the `arno_ccsdt` potential subroutine created as described above

```
/Users/mha/hib43/bin/progs> runpot_arno_ccsdt
  r (bohr)
6
 vsum
  2.03110015E-03 -1.58292849E-03  4.61031455E-03 -1.00506297E-03  9.77742918E-04 -1.80847723E-04  1.05529707E-04 -2.34828441E-05  9.56131165E-06
  vdif
 -5.94863282E-04  6.42939848E-04 -1.14560854E-03  2.24197698E-04 -1.92397409E-04  2.63754376E-05 -1.38655839E-05
  r (bohr)
```

