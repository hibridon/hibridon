# DIFFER

The command `DIFFER` requests a comparison of two S matrices, or of the moduli of two S matrices. The command line syntax is
```
differ,{jobnam1},{jobnam2},iprint,ienerg,thresh
```
where
- {jobnam1},{jobnam2}:  Names of the two S-matrix files to be compared, the first in the file {jobnam1}n.smt and the second in file {jobnam2}n.smt. Here n denotes the value of the parameter ienerg (see below). In order to create these S-matrix files, you must have carried out the prior calculation with the flag [`WRSMAT`](../param/wrsmat.md) = .TRUE.

- iprint:  Controls the amount of output, as follows:
  - iprint = 0:  Only the values of the average relative and largest relative differences are printed
  - iprint = 1:  The values of the average relative, average absolute, largest relative, and largest absolute differences are printed
  - iprint = 2:  All the S-matrix elements are printed followed by the values of the average relative, average absolute, largest relative, and largest absolute differences

- ienerg:  The cardinal value of the [energy](./energ.md) to which the S matrices refer. If ienerg = 2, say, then the second energy S matrices [{jobnam1}2.smt and {jobnam2}2.smt] are compared

- thresh:  All S-matrix elements smaller than thresh are ignored in the comparison. If thresh is a positive number (or zero), the comparison and the output involves the **moduli** of the S-matrix elements. If thresh is negative, the comparison and the output involve both the real and the imaginary parts of the S matrix.

The default values of these variables are:
- iprint = 1
- ienerg = 1
- thresh = 1.E-5
- ⚠️Both {jobnam1} and {jobnam2} must be specified

As an example, the command
```
differ,oldjob,newjob,0,2
```
would request a comparison between the moduli of the S matrices in files oldjob2.smt and newjob2.smt