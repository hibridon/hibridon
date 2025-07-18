The parameter NIOUT sets the length of the array [INDOUT](INDOUT) of quantum number(s) (in addition to the array [JOUT](JOUT)) for which
  - S-matrix elements will be printed out if the flag [WRSMAT](PRSMAT-and-WRSMAT) is set .TRUE.
  - Integral cross sections will be printed out in response to the command [INTCRS](INTCRS)
  - Partial cross sections will be printed out in response to the command [PARTC](PARTC)

The length of the array INDOUT is given by the absolute value of NIOUT. If [WRSMAT](PRSMAT-and-WRSMAT) is.TRUE., then

  - If NIOUT is POSITIVE, then the only S-matrix elements which will be saved are those for which both the initial and final values of the additional quantum numbers correspond to one of the values in the array INDOUT.
  - If NIOUT is NEGATIVE, any column of the S-matrix for which the initial value of the additional quantum number is included in the array INDOUT will be printed.

Changing or setting the value of NIOUT can be done only by setting the [INDOUT](INDOUT) parameter as follows:
```
INDOUT,NIOUT,indout(1),indout(2).....indout(NIOUT)
```
as, for example,
```
INDOUT,3,indout(1),indout(2),indout(3)
```
If the number of values of indout entered is less than NIOUT, then NIOUT is reset automatically. The string of input values should be terminated with a backslash (\\) if other commands follow on the same input line.

Click [here](Parameters) for more information on reading in, changing, or saving the values of any parameter.