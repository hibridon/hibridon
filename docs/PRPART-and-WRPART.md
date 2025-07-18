If PRPART = .TRUE., then partial cross sections will be printed out at the end of the regular [output](Output) file. In the case of [CS](Coupled-States) calculations the contribution of each value of nu is also printed. The columns of the printed matrix of partial cross sections correspond to the initial state and the rows, to the final state.

⚠️ Partial cross sections will be written only for transitions for which both the intial and final J are specified in the array [JOUT](JOUT) and for which both the initial and final index are specified in the array [JOUT](JOUT)

The printing of partial cross sections to a separate file is controlled by the flag WRPART.
If WRPART = .TRUE., then integral cross sections will be printed out to a separate [file](Files) {jobname}n.pcs where n designates the value of the energy, e.g. for the first energy the integral cross sections will be in the file {jobname}1.pcs, etc. Partial cross sections are computed exactly as described earlier on this page.

The contents of the file {jobname}n.pcs can be examined or printed using the command [PARTC](PARTC).