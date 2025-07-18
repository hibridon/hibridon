NERG is the number of total energies at which the calculation is to be performed. Both the log-derivative and AIRY propagation are significantly faster for the second and subsequent energies. Thus if cross sections are desired over a range of energies, it is advantageous to set NERG > 1.

The maximum value of NERG is controlled by a PARAMETER statement in the subroutine LOGAIR (himain.f) and can not exceed 25 (because of restrictions imposed by the length of the names of the logical units reserved for storing the partial and integral cross sections and selected S-matrix elements). If NERG>1, then three scratch files (units 10-12) are opened to store the required energy independent matrices as well as the channel quantum numbers.

The value of NERG can be changed with the command
```
NERG={newvalue}
```
If {newvalue} is larger than the value of NERG currently stored, then you will be warned to increase accordingly the array of energies at which the calculations are to be performed. This can be done with the command
```
ENERG={energy1},{energy2}...{energyn}\
```
where n is the new value of NERG. The string of input values must be terminated with a backslash if other input occurs on the same line.

Alternatively, you can specify a grid of collision energies b the command
```
ENERG = {negative number},{energy1},{energ2},{nerg}
```
where
  - energy1: The initial energy
  - energy2: The final energy
  - nerg: The number of energies (less than or equal to 25)

The spacing for the grid of energies is determined as
```
delta_e=(energ2-energ1)/nerg
```
This spacing is then rounded to the nearest 0.01 $\text{cm}^{-1}$, so that the final energy becomes
```
e_final=e_initial + (nerg-1)*delta_e
```
which may not be exactly equal to the final energy which you initially specified.

In all cases, the energy values are reordered in decreasing order, starting with the highest value.