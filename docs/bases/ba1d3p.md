### Collisions involving the <sup>1</sup>D and <sup>3</sup>P states of an atom with electron configuration p<sup>2</sup> or p<sup>4</sup> with a structureless atom


> BASISTYPE = 22
>
> System subroutine:  sy1d3p
>
> Basis subroutine:  hiba22_1d3p.F90
>
> Ref. P. J. Dagdigian, M. H. Alexander, and J. Klos J. Chem. Phys. **143**, 054306 (2015).

#### The definition of the integer system dependent parameters is as follows:

- NTERM: the number of potential surfaces involved. This parameter can not be changed

- NSTATE: number of electronic states included:
  - set NSTATE=0 for just the <sup>1</sup>D state
  - set NSTATE=1 for just the <sup>3</sup>P state
  - set NSTATE=2 for both <sup>1</sup>D and <sup>3</sup>P states

#### The definition of the real system dependent parameters is as follows:

- EN1D: energy of the <sup>1</sup>D state relative to that of ths <sup>3</sup>P state.