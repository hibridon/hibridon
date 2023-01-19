### Collisional Mixing of <sup>2</sup>&Pi; and <sup>2</sup>&Sigma; states of a diatomic molecule induced by a spherical atom


> BASISTYPE = 4
>
> System subroutine:  sysgpi
>
> Basis subroutine:  hiba04_sysgpi.F90
>
> Ref. M. H. Alexander and G. C. Corey, J. Chem. Phys. **84**, 100 (1986); H.-J. Werner, B. Follmeg, M. H. Alexander, and D. Lemoine, ibid. **91**, 5425 (1989).

#### The definition of the integer system dependent parameters is as follows:

- `NTERM`: (the number of potential surfaces involved: this should be 1 in this case). This parameter can not be changed
- `NVMINS`: lowest vibrational level for sigma state
- `NVMAXS`: higest vibrational level for sigma state
- `NVMINP`: lowest vibrational level for pi state
- `NVMAXP`: higest vibrational level for pi state
- `JMIN`: the minimum rotational angular momenta for each <sup>2</sup>&Pi; channel
- `JMAX`: the maximum rotational angular momenta for each <sup>2</sup>&Pi; channel in each spin-orbit manifold with convention omega .le. j .le. jmax+0.5
- `NMIN`: the minimum Hund's case (b) rotational angular momenta for the <sup>2</sup>&Sigma; state
- `NMAX`: the maximum Hund's case (b) rotational angular momenta for the <sup>2</sup>&Sigma; state

  ⚠️ `JMIN`, `JMAX`, `NMIN`, `NMAX` are defined separately for each vibrational level

- `IPERT`: Each vibrational Pi level ivp may be perturbed by one sigma vibrational level IVS. This level is given by ivs=ipert(ivp) (see hisysgpi)
- `IGUPI`: permutation inversion symmetry of <sup>2</sup>&Pi; electronic state
  - `IGU` = 1 for gerade states
  - `IGU` = -1 for ungerade states
- `IGUSG`: permutation inversion symmetry of <sup>2</sup>&Sigma; electronic state
  - `IGU` = 1 for gerade states
  - `IGU` = -1 for ungerade states

  ⚠️ for heteronuclear molecules both `IGUPI` and `IGUSG` should be +1

- `NPARPI`: number of <sup>2</sup>&Pi; symmetry doublets included
  - `NPAR` = 2 will ensure both lambda doublets
  - `NPARPI` = 1 just &epsilon; = 1 doublets
  - `NPARPI` = -1, just &epsilon; = -1 doublets
  - `NPARSG`: number of spin doublets in <sup>2</sup>&Sigma; state included (NPAR = 2 will ensure both spin doublets)
- `ISYMSG`: if ISYM =+1, then the electronic symmetry of the <sup>2</sup>&Sigma; state is sigma-plus if ISYM = -1, then the electronic symmetry is sigma-minus ISA: s/a symmetry index, if the molecule is homonuclear (ihomo=t) then, if isa=+1 then only the s-levels (both for sigma and pi) are included in the basis, if isa=-1, then only the a-levels are included
- `ISG`: if isg=1 and ipi=0 then <sup>2</sup>&Sigma; + atom scattering

- `IPI`: if ipi=1 and isg=0 then <sup>2</sup>&Pi; + atom scattering. If isg=1 and ipi=1 then <sup>2</sup>&Pi;-<sup>2</sup>&Sigma; + atom scattering