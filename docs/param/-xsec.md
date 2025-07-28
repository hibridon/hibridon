# integral cross section parameters

## PRXSEC

If `PRXSEC` = .TRUE., then integral cross sections will be printed out at the end of the regular output file. Integral cross sections are obtained by summing over all partial cross sections with

[`JTOT1`](../param/jtot1.md)  &leq; Jtot  &leq; [`JTOT2`](../param/jtot2.md)

If [`JTOTD`](../param/jtotd.md) > 1, a trapezoidal interpolation scheme is used to approximate this sum. The columns of the printed matrix of integral cross sections correspond to the initial state and the rows, to the final state.


In the case of [`CS`](../Coupled-states.md) calculations the integral cross sections represent a sum over all the values of the [`CS`](../Coupled-states.md) projection index nu lying between [`NUMIN`](../param/numin.md) and [`NUMAX`](../param/numin.md), inclusive. Only transitions into and/or out of rotational levels with J .le. [`NUMAX`](../param/numin.md) will be converged. In general, cross sections for only these transitions are printed. However, setting `IPRINT` = 2 forces printing of all the cross sections, even those which are not converged with respect to the sum over the projection index nu.

In the case of [`CC`](../Close-coupled-equations.md) calculations, to obtain integral cross sections summed over both values of the parity index [`JLPAR`](../param/jlpar.md), set `JLPAR` = 0.

## WRXSEC

The printing of integral cross sections to a separate file is controlled by the flag `WRXSEC`, as follows:
If `WRXSEC` = .TRUE., then integral cross sections will be printed out to a separate file *{jobname}n.ics* where *n* designates the cardinal value of the energy, e.g. for the first energy the integral cross sections will be in the file *{jobname}1.ics*, etc. Integral cross sections are computed exactly as described earlier on this page.


# notes

The contents of the file *{jobname}n.ics* can be examined or printed using the command [`PRINTC`](../com/printc.md).