### Collisions involving an atom in a <sup>1,3</sup>P electronic state


> BASISTYPE = 7
>
> System subroutine:  sy13p
>
> Basis subroutine:  hiba07_13p.F90
>
> Ref. B. Pouilly, T. Orlikowski, and M. H. Alexander, J. Phys. B **18**, 1953 (1985); B. Pouilly and M. H. Alexander, J. Chem. Phys. **86**, 4790 (1987); B. Pouilly, J.-M. Robbe, and M. H. Alexander, ibid. **91**, 1658 (1989).


#### The definition of the integer system dependent parameters is as follows:

- NTERM (the number of potential surfaces involved: this should be 1 in this case). This parameter can not be changed
- NSTATE: number of electronic states included
  - 0: just singlet state
  - 1: just triplet state
  - 2: both singlet and triplet state
- IPOL: selects whether $^1P_{1e}$ states are combined to form proper $^1\Pi_e$ and $^1\Sigma_e$ combinations [see B. Pouilly and M. H. Alexander, Chem. Phys. **145**, 191 (1990)]
  - 0: normal case, no linear combinations made
  - 1: linear combinations made (just for [JLPAR](JLPAR) = -1)

#### The definition of the real system dependent parameters is as follows:

- E-3P0, E-3P1, E-3P2: energies of $^3P$ spin-orbit levels
- E-1P1: energy of $^1P$ level relative to lower $^3P$ level than follow the parameters for the morse-spline-van der waals description of potential curves [see B. Pouilly and M. H. Alexander, J. Chem. Phys. **86**, 4790 (1987)]
- DE(3PI), DE(3SG), DE(1PI), DE(1SG): Morse De
- RE(3PI), RE(3SG), RE(1PI), RE(1SG): Morse Re
- BE(3PI), BE(3SG), BE(1PI), BE(1SG): Morse ß
- RL(3PI), RL(3SG), RL(1PI), RL(1SG): Onset of Long-range behavior
- C(3PI), C(3SG), C(1PI), C(1SG): Coefficent of R-6
- CMIX: extent of mixing between $^1P_1$ and $^3P_1$ levels; the mixing angle in a 2x2 rotation is theta=acos(CMIX)

#### The system dependent input paramaters are as follows:

- line 1: NSTATE, IPOL, NPOT
where NPOT is an integer which designates a particular choice of potential within the POT subroutine, and the other parameters are defined above
- line 2: E-3P0, E-3P1, E-3P2
- line 3: E-1P1, CMIX
- line 4: DE, RE, BE, RL, CL for 3 state
- line 5: DE, RE, BE, RL, CL for 3 state
- line 6: DE, RE, BE, RL, CL for 1 state
- line 7: DE, RE, BE, RL, CL for 1 state
- line 8: DEMOR, REMOR, BEMOR, DISSMOR (additional parameters, see the src/pot/pot_cahe.f subroutine)
- line 9: RGAUS, AGAUS, ALPHG (additional parameters, see the src/pot/pot_cahe.f subroutine)