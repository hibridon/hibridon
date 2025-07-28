1. [`JTOT1`](#jtot1)
2. [`JTOT2`](#jtot2)
3. [`JTOTD`](#jtotd)


## JTOT1
`JTOT1` is the initial value of the total angular momentum in a [`CC`](../Close-coupled-equations.md) calculation. For collisions involving a molecule with half- integer spin (flag [`FLAGHF`](./flaghf.md)) then the true total angular momentum is equal to Jtot+1/2.

For collisions of a molecule with an uncorrugated surface, you should set `JTOT1` = 0.

## JTOT2

`JTOT2` is the final value of the total angular momentum in a [`CC`](../Close-coupled-equations.md) calculation or the the orbital angular momentum (**L**) in a [`CS`](Coupled states) calculation. In the calculation of integral and/or differential cross sections, you must ensure that the parameter `JTOT2` has been set large enough so that increasing its value will not significantly change the cross section(s) of interest.

For collisions of a molecule with an uncorrugated surface, you should set `JTOT2`=0.

## JTOTD

`JTOTD` is the step size for the total angular momentum in a [`CC`](../Close-coupled-equations.md) calculation or the the orbital angular momentum (**L**) in a [`CS`](../Coupled-states.md) calculation. Scattering calculations are performed at every `JTOTD` value of JTOT (or **L**) lying between (and including) `JTOT1` and `JTOT2`. For collisions involving a molecule with half-integer spin (flag [`FLAGHF`](./flaghf.md)) then the total angular momentum is half integer (although  remains integer). In this case the true values of Jtot range from `JTOT1`+1/2 to `JTOT2`+1/2 in steps of `JTOTD`.

If the flags [`WRXSEC`](./wrxsec.md) and/or [`PRXSEC`](./prxsec.md) are .true., the program accumulates the contribution to the integral cross sections due to all partial waves lying between `JTOT1` and `JTOT2`. If `JTOTD` = 1, this is done by adding the degeneracy weighted sums of the squares of the T matrix at each partial wave. If, however, `JTOTD` > 1, then a Simpson's rule scheme is used to interpolate the missing T-matrix elements in this sum.

For collisions of a molecule with an uncorrugated surface, you should set `JTOTD` = 0.