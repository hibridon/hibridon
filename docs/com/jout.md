# JOUT

`JOUT` is an array containing the values of the rotational quantum number for which
- S-matrix elements will be printed out if the flag [`WRSMAT`]../param/wrsmat) is set .TRUE.
- Integral cross sections will be printed out in response to the command [`INTCRS`]../com/intcrs.md)
- Partial cross sections will be printed out in response to the command [`PARTC`]../com/partc.md)

The length of the array `JOUT` is given by the absolute value of the parameter [`NNOUT`]../param/nnout.md).
If [`WRSMAT`]../param/wrsmat.md) is.TRUE., then
  - If [`NNOUT`]../param/nnout.md) is POSITIVE, then the only S-matrix elements which will be saved are those for which both the initial and final rotational quantum numbers correspond to one of the values in the array JOUT.
  - If [`NNOUT`]../param/nnout.md) is NEGATIVE, any column of the S-matrix for which the initial rotational quantum number is included in the array JOUT will be printed.

For the printing of integral or partial cross sections, which is controlled by the commands [`INTCRS`]../com/intcrs.md) and [`PARTC`]../com/partc.md), the sign of [`NNOUT`]../param/nnout.md) is of no importance

Changing or entering the JOUT array is done by the command
```
JOUT,NNOUT,jout(1),jout(2).....jout(NNOUT)
```
as, for example,
```
JOUT,3,jout(1),jout(2),jout(3)
```
If the number of values of jout entered is less than [`NNOUT`]../param/nnout.md), then [`NNOUT`]../param/nnout.md) is reset automatically. The string of input values should be terminated with a backslash (\) if other commands follow on the same input line.