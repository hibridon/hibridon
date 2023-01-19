The values of the total energy at which the scattering calculation will be run are specified by the input line
```
ENERG = {energy1},{energy2}...{energyN}
```

If other input data occurs on the same line, then the string of energies must be terminated by a backslash.

Note that if the length of the array of energies is greater than the current value of [NERG](NERG), then [NERG](NERG) is reset accordingly. The maximum number of energies cannot exceed 25.

If the first energy is negative, then the calculations are to be preformed on a regular grid of energies, specified by
```
ENERG = {negative number},{energy1},{energ2},{nerg}
```
where
- energy1: The initial energy
- energy2: The final energy
- nerg: The number of energies (less than or equal to 25)

The spacing for the grid of energies is determined as

  >delta_e=(energ2-energ1)/nerg

This spacing is then rounded to the nearest 0.01 cm-1, so that the final energy becomes

  >e_final=e_initial + (nerg-1)*delta_e

which may not be exactly equal to the final energy which you initially specified.


In all cases, the energy values are reordered in decreasing order, starting with the highest value.
