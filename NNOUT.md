The parameter NNOUT sets the length of the array [JOUT](JOUT) which contains the values of the rotational quantum number for which
  - S-matrix elements will be printed out if the flag [WRSMAT](PRSMAT-and-WRSMAT) is set .TRUE.
  - Integral cross sections will be printed out in response to the command [INTCRS](INTCRS)
  - Partial cross sections will be printed out in response to the command [PARTC](PARTC)

The length of the [JOUT](JOUT) array is given by the absolute value of the parameter NNOUT.

If NNOUT is POSITIVE, then the only S-matrix elements which will be saved are those for which both the initial and final rotational quantum numbers correspond to one of the values in the array JOUT.

If NNOUT is NEGATIVE, any column of the S-matrix for which the initial rotational quantum number is included in the array JOUT will be printed.
For the printing of integral or partial cross sections, which is controlled by the commands [INTCRS](INTCRS) and [PARTC](PARTC), the sign of NNOUT is of no importance

Changing or setting the value of NNOUT can be done only by setting the [JOUT](JOUT) parameter as follows:
```
JOUT,NNOUT,jout(1),jout(2).....jout(NNOUT)
```
as, for example,
```
JOUT,3,jout(1),jout(2),jout(3)
```
If the number of values of jout entered is less than NNOUT, then NNOUT is reset automatically. The string of input values should be terminated with a backslash (\\) if other commands follow on the same input line.

Click [here](Parameters) for more information on reading in, changing, or saving the values of any parameter.