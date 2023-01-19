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

### [1. A<sup>1</sup>&Sigma; molecule and a spherical atom - ba1sg](Basis-subroutines---ba1sg)
Ref:  W. A. Lester, Jr., Meth. Comput. Phys.  10 , 211 (1971).


### [2. A<sup>2</sup>&Sigma; molecule and a spherical Atom - ba2sg](Basis-subroutines---ba2sg)
Ref:  M. H. Alexander, J. Chem. Phys.  76 , 3637 (1982).


### [3. A<sup>2</sup>&Pi; molecule and a spherical atom - ba2pi](Basis-subroutines---ba2pi)
Ref:  M. H. Alexander, J. Chem. Phys.  76 , 5974 (1982);\
M. H. Alexander, Chem. Phys.  92 , 337 (1985);\
G. C. Corey and M. H. Alexander, J. Chem. Phys.  85 , 5652 (1986).


### [4. A<sup>1</sup>&Delta; molecule and a spherical atom - ba1del](Basis-subroutines---ba1del)
Ref: D. G. Sauder, D. Patel-Misra, and P. J. Dagdigian, J. Chem. Phys.  91 , 5316 (1989).

[](ba2del.html) [](rightarrowsmall.gif)

### [5. A<sup>2</sup>&Delta; molecule and a spherical atom - ba2del](Basis-subroutines---ba2del)
Ref:  B. Nizamov, P. J. Dagdigian, Y.-R. Tzeng, and M. H. Alexander, J. Chem. Phys.  115 , 800 (2001).


### [6. Collisional mixing of  <sup>2</sup>&Pi;  and  <sup>2</sup>&Sigma;  states of a diatomic molecule induced by a spherical atom - basgpi](Basis-subroutines---basgpi)
Ref:  M. H. Alexander and G. C. Corey, J. Chem. Phys.  84 , 100 (1986);\
H.-J. Werner, B. Follmeg, M. H. Alexander, and D. Lemoine,  ibid.   91 , 5425 (1989).


### [7. Collisions of a molecule in a  <sup>1</sup>&Pi;, <sup>2</sup>&Pi;, or <sup>3</sup>&Pi; state with a spherical atom - bapi](Basis-subroutines---bapi)
Ref:  M. H. Alexander, J. Chem. Phys.  76 , 5974 (1982);\
M. H. Alexander, Chem. Phys.  92 , 337 (1985);\
G. C. Corey and M. H. Alexander, J. Chem. Phys.  85 , 5652 (1986);\
D. Lemoine, G. C. Corey, M. H. Alexander, and J. Derouard, Chem. Phys.  118 , 357 (1987).


### [8. A singlet symmetric top molecule and a spherical atom - bastp](Basis-subroutines---bastp)
Ref:  S. Green, J. Chem. Phys.  64 , 3463 (1976);  67 , 816 (1979).


### [9. A singlet asymmetric top and a spherical atom - baastp](Basis-subroutines---baastp)
Ref:  B. J. Garrison  et al. , J. Chem. Phys.  65 , 2193 (1976);\
S. Green, J. Chem. Phys.  64 , 3463 (1976);  67 , 816 (1979).


### [10. An atom in a  <sup>1</sup>P and/or  <sup>3</sup>P state and a spherical atom - bas13p](Basis-subroutines---ba13p)
Ref:  B. Pouilly, T. Orlikowski, and M. H. Alexander, J. Phys. B  18 , 1953 (1985);\
B. Pouilly and M. H. Alexander, J. Chem. Phys.  86 , 4790 (1987);\
B. Pouilly, J.-M. Robbe, and M. H. Alexander,  ibid.   91 , 1658 (1989).


[](ba2mol.html) [](rightarrowsmall.gif)

### [11. Two  <sup>1</sup>&Sigma; diatomics - ba2mol](Basis-subroutines---ba2mol)
Ref:  S. Green, J. Chem. Phys.  62 , 2271 (1975);\
A. E. DePristo and M. H. Alexander,   ibid.   66 , 1334 (1977);\
M. H. Alexander,   ibid.   73 , 5135 (1980).


### [12. A symmetric top and a linear molecule - bastpln](Basis-subroutines---bastpln)
Ref:  C. Rist, M. H. Alexander, and P. Valiron, J. Chem. Phys.  98 , 4662 (1993).

[](ba22p.html) [](rightarrowsmall.gif)

### [13. Collisions of an atom in a  <sup>2</sup>P state and an atom in a <sup>2</sup>S electronic state - ba22p](Basis-subroutines---ba22p)
Ref:  M. H. Alexander, B. Pouilly, and T. Duhoo, J. Chem. Phys.  99 , 1752 (1993).

### [14. An atom in a  <sup>2</sup>P electronic state and a homonuclear diatomic - bah2p](Basis-subroutines---bah2p)
Ref:  M.-L. Dubernet and J. M. Hutson, J. Chem. Phys.  101, 1939 (1994).


### [16. An atom in a <sup>3</sup>P electronic state and a homonuclear diatomic - bah3p](Basis-subroutines---bah3p)
Ref:  M.-L. Dubernet and J. M. Hutson, J. Chem. Phys.  101, 1939 (1994).

### [17. Collisions between a  <sup>2</sup>&Delta; Molecule and a Spherical Atom](Basis-subroutines---ba2del)

### [18. Collisions between of an atom with a symmetric top, with no inversion doubling](Basis-subroutines---bastp1)

### [19. Collisional mixing of <sup>2</sup>&Pi; and <sup>2</sup>&Sigma; states of a diatomic molecule induced by a spherical atom (no spectroscopic perturbations](Basis-subroutines---basgpi1)

### [20. Collisions between a <sup>2</sup>&Pi; molecule and a <sup>1</sup>&Sigma; molecule](Basis-subroutines---ba2pi1sg)

### [21. Collision of a symmetric top with a linear molecule in a <sup>1</sup>&Sigma; state, using coupled representation of the PES](Basis-subroutines---bastp1sg)

### [22. Collisions involving the <sup>1</sup>D and <sup>3</sup>P states of an atom with electron configuration p<sup>2</sup> or p<sup>4</sup> with a structureless atom](Basis-subroutines---ba1d3p)

### [23. Collisions of a <sup>3</sup>P atom  with a <sup>4</sup>p electron configuration in a <sup>2</sup>S state](Basis-subroutines---ba3p2s)

### [24. Collision of a spherical top molecule with an atom ](Basis-subroutines---basphtp)

### [25. Collisions of two unlike a closed-shell (<sup>1</sup>&Sigma;) diatomic molecule](Basis-subroutines---ba1sg1sg)

### [26. Collisions of a <sup>2</sup>&Sigma; diatomic molecule with a closed-shell (<sup>1</sup>&Sigma;) diatomic molecule](Basis-subroutines---ba2sg1sg)

### [27. Collisions of an atom with an asymmetric top](Basis-subroutines---baastp1)

### [28. Collision of a <sup>3</sup>&Sigma; diatomic molecule with a closed-shell (<sup>1</sup>&Sigma;) diatomic molecule](Basis-subroutines---ba3sg1sg)

### [29. Collision of an atom with an asymmetric top](Basis-subroutines---baastp2)

### [30. Collisions of an asymmetric top with a closed-shell linear molecule](Basis-subroutines---baastp3)


### [99. User-defined basis routine - bausr](Basis-subroutines---bausr)

Setting the parameter `BASISTYPE = 99` allows the user to define his own basis routine, in addition to the fifteen types enumerated above.  For an example see the file `tests/ch3i/pot_ch3i.F90`