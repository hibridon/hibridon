The Hibridon code can treat photodissociation (in a time-independent manner) as well as scattering. The former involves only a slight change in the asymptotic [boundary conditions](BOUNDARY). In addition a separate subroutine GROUND must be included which determines the ground state (bound) wavefunction mutiplied by the transition dipole matrix projected onto the same channel basis used to expand the photofragment (scattering) wavefunction.

For a good reference see

D. E. Manolopoulos and M. H. Alexander, J. Chem. Phys. **97**, 2527 (1992); M. H. Alexander, C. Rist, and D. E. Manolopoulos, ibid **97**, 4836 (1992).

For more information examine the potential subroutine for the photodissociation of the CH<sub>3</sub>I molecule used in automated testing of Hibridon:
```
tests/ch3i/pot_ch3i.F90
```

If PHOTOF = .FALSE., then inelastic scattering [boundary conditions](BOUNDARY) are assumed.
If PHOTOF = .TRUE., then photodissociation [boundary conditions](BOUNDARY) are assumed.