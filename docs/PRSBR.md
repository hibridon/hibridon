Ref. P. J. Dagdigian, M. H. Alexander, and J. Klos J. Chem. Phys. **143**, 054306 (2015)

The command `PRSBR` allows you to calculate pressure-broadening cross sections.

For details, see R. Shafer and R. G. Gordon, J. Chem. Phys. **58**, 5422 (1973); S. Green, J. Boissoles, and C. Boulet, J. Quant. Spectrosc. Radiat. Transfer **39**, 33 (1988).

With this module, the pressure broadening cross section for an isolated spectral line, or collisionally coupled lines, can be computed. The kinetic energies for the initial and final levels should be the same in the two scattering calculations. In addition, [`JTOT1`](JTOT.md) and [`JTOT2`](JTOT.md) should be the same for both calculations.

The command line is
```
PRSBR,<job1>,<ienerg1>,<job2>,<ienerg2>,<k>,<j1>,<in1>,<j2>,<in2>,<diag>,<j1p>,<in1p>,<j2p>,<in2p>
```
where the input parameters are:

- `<job1>`,`<ienerg1>`: These parameters locate the S-matrix for a rotational level in the initial vibrational or electronic state. The program searches for an S-matrix in the file `<job1><ienerg1}.smt`, where `<ienerg1>` is an integer (`<ienerg1>` = 1, 2, ... etc.)

- `<job2>`,`<ienerg2>`: These parameters locate the S-matrix for a rotational level in the final vibrational or electronic state. The program searches for an S-matrix in the file `<job2><ienerg2>.smt`, where `<ienerg2>` is an integer (`<ienerg2>` = 1, 2, ... etc.)

- `<k>`: This is the tensor order of the cross section:
  - `<k>` = 0 for isotropic Raman scattering,
  - `<k>` = 1 for dipole-allowed transitions,
  - `<k>` = 2 for anisotropic Raman scattering.

- `<j1>`: The rotational angular momentum of the rotational level in the initial state.

- `<in1>`: The state index for the rotational level in the initial state.

- `<j2>`: The rotational angular momentum of the rotational level in the final state.

- `<in2>`: The state index for the rotational level in the final state.

- `<diag>`: This parameter equals one if a cross section for an isolated spectral line is to be computed. `<diag>` should be set to zero to consider the collisional coupling of two spectral lines.

The parameters below are required only if `<diag>` = 0.

- `<j1p>`: The rotational angular momentum of the rotational level in the  initial state of the coupled spectral line.

- `<in1p>`: The state index for the rotational level in the initial state of the coupled spectral line.

- `<j2p>`: The rotational angular momentum of the rotational level in the final state of the coupled spectral line.

- `<in2p>`: The state index for the rotational level in the final state of the coupled spectral line.

The real and imaginary parts of the pressure broadening cross section are printed in the terminal output, as well as in an auxiliary file, `<job1><ienerg1>.ppb`. Also outputted to this file are the partial cross sections vs. `JTOT`.

⚠️ For accurate calculation of the pressure broadening cross section, especially the imaginary part, Logd integration should be used over the full range of $R$.
