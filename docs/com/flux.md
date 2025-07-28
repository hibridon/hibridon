# FLUX

The command FLUX requests the determination of the collision or photodissociation flux as a function of R using using data stored during a previous calculation. For more details on the definition and determination of the flux see:

M. H. Alexander, J. Chem. Phys. **95**, 8931 (1991); **96**, 6672 (1992). D. E. Manolopoulos and M. H. Alexander, J. Chem. Phys. **97**, 2527 (1992); M. H. Alexander, C. Rist, and D. E. Manolopoulos, J. Chem. Phys. **97**, 4836 (1992).

The command line syntax is
```
FLUX,{jobnam},IFLUX,{additional variables}
```
where
  - {jobname}:  the jobname under which the wavefunction information from the previous calculation is stored, as the file {jobname}.wfu

  ⚠️ In order to create the {jobname}.wfu file, the previous calculation must have been carried out with the flags [`WAVEFL`]../param/wavefl.md) and [`AIRYFL`]AIRYFL-and-PRAIRY) both set .TRUE.

  - `IFLUX` =  
    - 1:  for a determination of the fluxes in the diabatic (asymptotic) basis
    - −1:  for a determination of the fluxes in the locally adiabatic basis
    - 2:  for a determination of the adiabatic energies (see also the [`EADIAB`]./eadiab.md) command)
    - −2:  for fluxes summed over the sign of the additional channel index
    - 3:  for fluxes in coordinate space
    - 4:  for determination of the transformation matrix C(R), defined in Eq. (6) of the help file describing the [close coupling method](../Close-coupled-equations.md)


If `IFLUX = −1,1, or −2` then the {additional variables} are `THRESH`, `IPRINT`, `DAMP`, `INITIAL-J`, `INITIAL-L`, and `INITIAL-IND`. These are defined as follows:
  - `THRESH`:  If the local wavevector in any channel is less than THRESH, then the wavefunction component associated with this channel is killed off. The default value of THRESH is 0 for scattering calculations, and -1.e+9 for photodissociation calculations

  - `IPRINT`:  Normally, fluxes are not printed inside of any energetically closed region. To print these out set IPRINT not equal to zero. The default value of IPRINT is 0 for scattering and 1 for photodissociation

  - `DAMP`:  the damping factor for closed channel components (the default value is one)

  - `INITIAL-J`:  the rotational angular momentum of the initial channel

  - `INITIAL-L`:  the orbital angular momentum of the initial channel

  - `INITIAL-IND`:  the value of the additional index of the initial channel

  - In addition, channel fluxes are printed out ONLY for those channels which are specified in the arrays [`JOUT`]../param/jtot.md) and [`INDOUT`]../param/indout.md).

❗For a photodissociation calculation it is not nessary to specify `THRESH`, `IPRINT`, `DAMP`, `INITIAL-J`, `INITIAL-L`, or `INITIAL-IND`. In this case the command line should be
```
FLUX,JOB,IFLUX
```

If `IFLUX = 3` then the {additional variables} are `NR`, `RMIN`, `DR`, `THRESH`, `IPRINT`, and `DAMP`, where

  - `NR`:  | NR | is the number of coordinate space points at which the flux is to be determined
    - if NR > 0, then the calculated flux is a sum over the flux associated with all channels with positive values of the additional index
    - if NR < 0, then the calculated flux is a sum over the flux associated with all channels with negative values of the additional index

  - `RMIN`:  Minimum value of the internal coordinate r

  - `DR`:  Step size in the internal coordinate r

  - `THRESH`:  as defined above

  - `IPRINT`:  as defined above

  - `DAMP`:  as defined above

❗Coordinate state fluxes will be determined for all values of the internal coordinate r ranging from RMIN to RMIN+(NR -1) DR


If `IFLUX = 4` then {additional variables} designates

  - `nstate`:  The number of states for which the asymptotic → locally adiabatic transformation matrix is to be printed out
  - `Rmin`:  The minimum value of R at which the asymptotic → locally adiabatic transformation matrix is to be printed

  - `DR`:  The step size

  - `Rmax`:  The maximum value of R at which the asymptotic --> locally adiabatic transformation matrix is to be printed

If the variable DR is ≤ 0, then the transformation matrix is printed out just for R=Rmin. Otherwise, the transformation matrix is printed out for R = Rmin : DR : Rmax.

❗nstate cannot exceed 12, it will be truncated if the input value exceeds 12.

⚠️ Fluxes are determined at all values of R lying between [`RENDLD`]RSTART,-RENDLD,-and-RENDAI) and [`RENDAI`]RSTART,-RENDLD,-and-RENDAI), with the grid size and spacing controlled by the same parameters which govern the AIRY integration, namely [`SPAC`]../param/spac.md), [`FSTFAC`]../param/fstfac.md), [`TOLAI`]TOLAI-and-RINCR), and [`RINCR`]TOLAI-and-RINCR).

⚠️ The determination of collision and photodissociation fluxes should be attempted only if you have a sophisticated understanding of the time-independent quantum description of molecular collisions.