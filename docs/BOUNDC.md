The Hibridon code can determine energies and eigenvectors of weakly-bound states as well as scattering and photodissociation cross sections. This involves only a slight change in the [methodology](Bound-States). To determine bound-state energies, set:

`BOUNDC` = .TRUE.

This will invoke a variational calculation of all bound-state levels using the distributed Gaussian method of I. P. Hamilton and J. C. Light [J. Chem. Phys. **84**, 306 (1986)]. The R-dependence of the wavefunction will be expanded in a set of Gaussian functions:

$\chi_m(R) = e^{-\alpha(R-R_m)^2}$

This set of functions is defined by several of the Hibridon input [parameters](Parameters), as follows:
- `R1`: smallest value of *R<sub>m</sub>*
- `R2`: largest value of *R<sub>m</sub>*
- `SPAC`: spacing between successive values of *R<sub>m</sub>*
- `C`: parameter which determines the exponential scale factor of the distributed Gaussian functions, with *&alpha;*=(`C`/`SPAC`)<sup>2</sup>.
- `EIGMIN`: lower limit on the minimum allowed eigenvalue of the overlap matrix. If the minimum eigenvalue is less than this value, the parameter `C` should be increased.
- `DELR`, `HSIMP`: The matrix elements of W(R) are evaluated by a Simpson's rule integration extending from `R1`–`DELR` to `R2`+`DELR` in steps of `HSIMP`.

Reasonable starting values for these parameters are (see the file [bound_parameters](https://www2.chem.umd.edu/groups/alexander/hibridon/hib43/bound_parameters.html) for a detailed discussion of setting these parameters):
  - `SPAC`: 0.4
  - `C`: 0.5
  - `EIGMIN`: 1.e-6
  - `DELR`: 1
  - `DR`: 0.1

If the flag `CSFLAG` is .TRUE., then a bound-state calculation will be carried out within the centrifugal decoupling approximation including only channels for which the projection of the total angular momentum along `R` is equal to
- [`NUMIN`](NUMIN), (if [`FLAGHF`](FLAGHF) = .FALSE.)
- [`NUMIN`](NUMIN)+0.5, (if [`FLAGHF`](FLAGHF) = .TRUE.)

If `BOUNDC` = .FALSE. (the default), then a scattering or photodissociation calculation is performed.