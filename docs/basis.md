# Basis subroutines
In addition to the general parameters, input parameters must be supplied to define the particular collision system to be treated and to control the choice of channels and internal energies.

Currently the Hibridon package can accomodate collision systems listed at the bottom of this page. Consult the page dedicated to each basis for a further description of its system-dependent variables.

note:

> For systems (1) - (7), (9), and (11) below, by setting the [flag](./Flags.md) [`FLAGSU`](./param/flagsu.md) = `.TRUE.`, it is possible to treat the collision of the
partner (atom or molecule) with internal structure with a flat, rigid surface.


For each of the collision systems that can currently be treated by the Hibridon package, the channel list, channel energies, and coupling matrix elements of the potential are determined by an appropriate BASIS subroutine. This is done before the propagation is initiated. A specific description of these particular BASIS subroutines can be obtained by clicking the button for the following cases:

The appropriate BASIS (and [`SYSDAT`](sysdat)) subroutines are selected by the input parameter `BASISTYP`, whose value is the index of the following ordered list.

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
will print out a sample input file called `Test.inp` , which you can then examine.


## Basis Routines

### [1. A<sup>1</sup>&Sigma; molecule and a spherical atom - ba1sg](./bases/ba1sg.md)

### [2. A<sup>2</sup>&Sigma; molecule and a spherical atom - ba2sg](./bases/ba2sg.md)

### [3. A<sup>2</sup>&Pi; molecule and a spherical atom - ba2pi](./bases/ba2pi.md)

### [4. Collisional Mixing of <sup>2</sup>&Pi; and <sup>2</sup>&Sigma; states of a diatomic molecule induced by a spherical atom](./bases/basgpi.md)

### [5. Collisions of a molecule in a  <sup>1</sup>&Pi;, <sup>2</sup>&Pi;, or <sup>3</sup>&Pi; state with a spherical atom - bapi](./bases/bapi.md)

### [6. A singlet symmetric top molecule and a spherical atom - bastp](./bases/bastp.md)

### [7. An atom in a  <sup>1</sup>P and/or  <sup>3</sup>P state and a spherical atom - ba13p](./bases/ba13p.md)

### [8. Two  <sup>1</sup>&Sigma; diatomics - ba2mol](./bases/ba2mol.md)

### [9. A symmetric top and a linear molecule - bastpln](./bases/bastpln.md)

### [10. Collisions of an atom in a  <sup>2</sup>P state and an atom in a <sup>2</sup>S electronic state - ba22p](./bases/ba22p.md)

### [11. Collisions between a  <sup>1</sup>&Delta; molecule and a spherical atom](./bases/ba1del.md)

### [12. An atom in a  <sup>2</sup>P electronic state and a homonuclear diatomic - bah2p](./bases/bah2p.md)

### [13. An atom in a <sup>3</sup>P electronic state and a homonuclear diatomic - bah3p](./bases/bah3p.md)

### [14. Collisions between a <sup>2</sup>&Delta; molecule and a spherical atom](./bases/ba2del.md)

### [15. Collisions between an atom in a <sup>2</sup>P electronic state and a heteronuclear diatomic](./bases/badiat2p.md)

### [16. Collisions between of an atom with a symmetric top](./bases/baastp.md)

### [17. Collision of an Atom with a CH<sub>2</sub>(X <sup>3</sup>B<sub>1</sub>) (0,v<sub>2</sub>,0) bender vibrational level](./bases/bach2x.md)

### [18. Collisions between of an atom with a symmetric top, with no inversion doubling](./bases/bastp1.md)

### [19. Collisional mixing of <sup>2</sup>&Pi; and <sup>2</sup>&Sigma; states of a diatomic molecule induced by a spherical atom (no spectroscopic perturbations](./bases/basgpi1.md)

### [20. Collisions between a <sup>2</sup>&Pi; molecule and a <sup>1</sup>&Sigma; molecule](./bases/ba2pi1sg.md)

### [21. Collision of a symmetric top with a linear molecule in a <sup>1</sup>&Sigma; state, using coupled representation of the PES](./bases/bastp1sg.md)

### [22. Collisions involving the <sup>1</sup>D and <sup>3</sup>P states of an atom with electron configuration p<sup>2</sup> or p<sup>4</sup> with a structureless atom](./bases/ba1d3p.md)

### [23. Collisions of a <sup>3</sup>P atom  with a <sup>4</sup>p electron configuration in a <sup>2</sup>S state](./bases/ba3p2s.md)

### [24. Collision of a spherical top molecule with an atom ](./bases/basphtp.md)

### [25. Collisions of two unlike a closed-shell (<sup>1</sup>&Sigma;) diatomic molecule](./bases/ba1sg1sg.md)

### [26. Collisions of a <sup>2</sup>&Sigma; diatomic molecule with a closed-shell (<sup>1</sup>&Sigma;) diatomic molecule](./bases/ba2sg1sg.md)

### [27. Collisions of an atom with an asymmetric top](./bases/baastp1.md)

### [28. Collision of a <sup>3</sup>&Sigma; diatomic molecule with a closed-shell (<sup>1</sup>&Sigma;) diatomic molecule](./bases/ba3sg1sg.md)

### [29. Collision of an atom with an asymmetric top](./bases/baastp2.md)

### [30. Collisions of an asymmetric top with a closed-shell linear molecule](./bases/baastp3.md)


### 99. User-defined basis routine - bausr

Setting the parameter `BASISTYPE = 99` allows the user to define his own basis routine, in addition to the types enumerated above. For an example see the file `tests/ch3i/pot_ch3i.F90`