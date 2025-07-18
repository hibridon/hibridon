In the [CC](Close-coupled-equations) formulation, there is no coupling between channels whose internal wavefunction have opposite total parity. JLPAR is the parameter which selects the parity of the channels which are included in a CC calculation.

If JLPAR=+1 or -1, then for molecules in electronic states of odd multiplicity (integer spin) only those channels for which
```
e(-1) J - S - s + L - Jtot = JLPAR
```
are included in the channel basis. For molecules in electronic states of even multiplicity (half-integer spin) only those channels for which
```
e(-1)J - S - s + L - Jtot - 1/2 = JLPAR
```
are included in the channel basis. Here
  - e = symmetry index of molecular state in Hund's case (a) (e = +1 or -1)
  - S = total spin
  - s = 1 for sigma-minus states, 0 otherwise
  - L = orbital angular momentum

❗ In general, in the [CC](Close-coupled-equations) formulation, calculations for both values of JLPAR must be performed to determine integral and/or differential cross sections. Setting JLPAR = 0 in a CC calculation will ensure that the calculation is done first for JLPAR = 1, and subsequently for JLPAR = -1.

❗ The parameter JLPAR is ignored in [coupled states](Coupled-States) calculations.