For each of the twelve collision [systems](systems) that can currently be treated by the Hibridon package, the channel list, channel energies, and coupling matrix elements of the potential are determined by an appropriate BASIS subroutine. This is done before the propagation is initiated. A specific description of these particular BASIS subroutines can be obtained by clicking the button for the following cases:

The appropriate BASIS (and [SYSDAT](sysdat) subroutines are selected by the input parameter `BASISTYP`, whose value is the index of the following ordered list.

To print out the list of system dependent parameters, enter:
```
BASISTY=i
```

where  `i`  is an integer ranging from 0 to 15, and 99. Then, enter:
```
SHOW
```

The subsequent command:
```
save=test.inp
```
will print out a sample input file called  `Test.inp` , which you can then examine.


## Basis Routines

### [1. A<sup>1</sup>&Sigma; molecule and a spherical atom - ba1sg](Basis-subroutines:-ba1sg)

### [2. A<sup>2</sup>&Sigma; molecule and a spherical atom - ba2sg](Basis-subroutines:-ba2sg)

### [3. A<sup>2</sup>&Pi; molecule and a spherical atom - ba2pi](Basis-subroutines:-ba2pi)

### [4. Collisional Mixing of <sup>2</sup>&Pi; and <sup>2</sup>&Sigma; states of a diatomic molecule induced by a spherical atom](Basis-subroutines:-basgpi)

### [5. Collisions of a molecule in a  <sup>1</sup>&Pi;, <sup>2</sup>&Pi;, or <sup>3</sup>&Pi; state with a spherical atom - bapi](Basis-subroutines:-bapi)

### [6. A singlet symmetric top molecule and a spherical atom - bastp](Basis-subroutines:-bastp)

### [7. An atom in a  <sup>1</sup>P and/or  <sup>3</sup>P state and a spherical atom - bas13p](Basis-subroutines:-ba13p)

### [8. Two  <sup>1</sup>&Sigma; diatomics - ba2mol](Basis-subroutines:-ba2mol)

### [9. A symmetric top and a linear molecule - bastpln](Basis-subroutines:-bastpln)

### [10. Collisions of an atom in a  <sup>2</sup>P state and an atom in a <sup>2</sup>S electronic state - ba22p](Basis-subroutines:-ba22p)

### [11. Collisions between a  <sup>1</sup>&Delta; molecule and a spherical atom](Basis-subroutines:-ba1del)

### [12. An atom in a  <sup>2</sup>P electronic state and a homonuclear diatomic - bah2p](Basis-subroutines:-bah2p)

### [13. An atom in a <sup>3</sup>P electronic state and a homonuclear diatomic - bah3p](Basis-subroutines:-bah3p)

### [14. Collisions between a <sup>2</sup>&Delta; molecule and a spherical atom](Basis-subroutines:-ba2del)

### [15. Collisions between an atom in a <sup>2</sup>P electronic state and a heteronuclear diatomic](Basis-subroutines:-badiat2p)

### [16. Collisions between of an atom with a symmetric top](Basis-subroutines:-bastp)

### [17. Collision of an Atom with a CH<sub>2</sub>(X <sup>3</sup>B<sub>1</sub>) (0,v<sub>2</sub>,0) bender vibrational level](Basis-subroutines:-bach2x)

### [18. Collisions between of an atom with a symmetric top, with no inversion doubling](Basis-subroutines:-bastp1)

### [19. Collisional mixing of <sup>2</sup>&Pi; and <sup>2</sup>&Sigma; states of a diatomic molecule induced by a spherical atom (no spectroscopic perturbations](Basis-subroutines:-basgpi1)

### [20. Collisions between a <sup>2</sup>&Pi; molecule and a <sup>1</sup>&Sigma; molecule](Basis-subroutines:-ba2pi1sg)

### [21. Collision of a symmetric top with a linear molecule in a <sup>1</sup>&Sigma; state, using coupled representation of the PES](Basis-subroutines:-bastp1sg)

### [22. Collisions involving the <sup>1</sup>D and <sup>3</sup>P states of an atom with electron configuration p<sup>2</sup> or p<sup>4</sup> with a structureless atom](Basis-subroutines:-ba1d3p)

### [23. Collisions of a <sup>3</sup>P atom  with a <sup>4</sup>p electron configuration in a <sup>2</sup>S state](Basis-subroutines:-ba3p2s)

### [24. Collision of a spherical top molecule with an atom ](Basis-subroutines:-basphtp)

### [25. Collisions of two unlike a closed-shell (<sup>1</sup>&Sigma;) diatomic molecule](Basis-subroutines:-ba1sg1sg)

### [26. Collisions of a <sup>2</sup>&Sigma; diatomic molecule with a closed-shell (<sup>1</sup>&Sigma;) diatomic molecule](Basis-subroutines:-ba2sg1sg)

### [27. Collisions of an atom with an asymmetric top](Basis-subroutines:-baastp1)

### [28. Collision of a <sup>3</sup>&Sigma; diatomic molecule with a closed-shell (<sup>1</sup>&Sigma;) diatomic molecule](Basis-subroutines:-ba3sg1sg)

### [29. Collision of an atom with an asymmetric top](Basis-subroutines:-baastp2)

### [30. Collisions of an asymmetric top with a closed-shell linear molecule](Basis-subroutines:-baastp3)


### [99. User-defined basis routine - bausr](Basis-subroutines:-bausr)

Setting the parameter `BASISTYPE = 99` allows the user to define his own basis routine, in addition to the fifteen types enumerated above.  For an example see the file `tests/ch3i/pot_ch3i.F90`