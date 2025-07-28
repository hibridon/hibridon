# S Matrix parameters

## PRSMAT

If PRSMAT = .TRUE., then the real and the imaginary S-matrix elements will be printed out at the end of the regular output file. To print the S matrix to a separate file, set the flag WRSMAT = .TRUE.

⚠️ S matrices can be quite large, so be wary of setting PRSMAT = .TRUE.

## WRSMAT

If WRSMAT = .TRUE., then selected (or all) S matrix elements will be printed out to a separate file {jobname}n.smt, where n designates the value of the energy, e.g. at the first energy the S-matrix elements will be in the file {jobname}n.smt, etc. The subset of S-matrix matrix elements which will be printed are defined by the set of rotational quantum numbers contained in the array [`JOUT`](../com/jout.md) as follows:
- If NNOUT is positive, then the only S-matrix elements which will be saved are those for which BOTH the initial and final rotational quantum numbers correspond to one of the value in the array [`JOUT`](../com/jout.md)
- If NNOUT is negative, any column of the S matrix for which the initial rotational quantum number is included in the array [`JOUT`](../com/jout.md) will be printed

# notes

The [command](../Commands.md) [`PRINTS`](../com/prints.md) allows the printing of an S matrix which has been previously written to the file {jobname}n.smt.