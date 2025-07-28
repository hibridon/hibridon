# DIFCRS

The command DIFCRS allows you to compute from previously determined S-matrix elements the
- differential cross sections
- diagonal quadrupole and octupole alignment moments $A_0^{(2)}$ and $A_0^{(4)}$
- off-diagonal alignment moment $A_2^{(2)}$
- m-dependent cross sections
- steric (oriented) cross sections.

The units of the cross setions are square Angstroms. The moments are dimensionless. The definition of the moments $A_0^{(2)}$ and $A_0^{(4)}$ is given in C. H. Greene and R. N. Zare, J. Chem. Phys. **78**, 6741 (1983). The off-diagonal alignment moment $A_2^{(2)}$ is defined [here](https://www2.chem.umd.edu/groups/alexander/hibridon/hib43/A22plus_short.pdf).

```
difcrs,{jobname},j1,in1,j2,in2,ang1,ang2,dang,ienerg,jtotend,ipr,mflag,stflag,$\alpha$,$\beta$
```
where
- {jobnam} the jobname under which the S matrices have been stored {jobname}n.smt. Here n denotes the value of the parameter ienerg (see below). These S matrices must have been previously generated using  
  - [`JLPAR`](../param/jlpar.md) = 0 (this generates S matrices for both parities)
  - [`JTOT1`](../param/jtot1.md) = 1 (ensuring the determination of S-matrices at every partial wave)
  - [`WRSMAT`](../param/wrsmat.md) = .true. (ensuring that the S-matrices are written to file {jobname}.smt).
The default value of {jobname} is the value you have set with the command [`JOB`]Files), or, if no value has been set, {jobname}=JOB

- j1,in1  rotational quantum number and additional index for the initial state

- j2,in2   rotational quantum number and additional index for the final state

- ang1,ang2  initial and final angle (in degrees)

- dang   step size (in degrees) for scan through angles.

- ienerg  the cardinal value of the [energy](./energ.md) for which the differential cross section is computed. i.e. if ienerg = 2, then the second energy S matrices [{jobname}2.smt] are used

- jtotend   the maximum value of $\mathbf{J}_{\rm tot}$ included in determining the scattering amplitude. The value of jtotend can not be larger than the variable [`JTOT2`]../param/jtot.md) used in the initial calculation

- ipr: 
  - If ipr = 0 (the default) the degeneracy-averaged differential cross sections and product rotational alignment and hexapole moments are NOT printed to the normal output file but only to the file {jobnam}n.dcs
  - if ipr .ne. 0, the differential cross sections and the alignment (A0(2)) and hexapole moments (A0(4)) of the products are also printed to the normal output file and to stdout.

- mflag:
  - If mflag = 0 (the default) only the degeneracy-averaged differential cross section and the alignment (A0(2)) and hexapole moments (A0(4)) of the products are calculated.
  - If mflag .ne. 0, then
    - All m→m' differential cross sections are calculated and printed. Quantization is in the collision frame, where the initial relative velocity vector defines the z axis.
    - The integral m→ m' cross sections are determined, by integration from ang1 to ang2 in steps of dang. 
    - The diagonal (ρm,m) and the real part of the 2nd supra-diagonal (ρm,m+2) elements of the rotational density matrix of the scattered products are output into file `<jobnam>n.rho`. These quantities are defined here in terms of the scattering amplitudes.

- stflag:
  - stflag = 0 is the default. If stflag .ne. 0, then "heads" and "tails" steric cross sections are calculated (see M. H. Alexander, Faraday Discuss. **113**, 437 (1999). This is only allowed if [`FLAGHF`](../param/flaghf.md) = .true. and [`BASISTY`](../basis.md) = 3 (doublet pi). If stflag .ne. 0, then the following two parameters must be defined:
    - $\alpha, \beta$   parameters which define the mixture of e and f lamda-doublet states in the "heads" or "tails" orientation.

NB: the "heads" state for m-initial negative is defined as:

| heads > = α | jme > - β | jmf >

and for m-initial positive as
| heads > = α | jme > + β | jmf >