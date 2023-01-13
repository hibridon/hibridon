# Hibridon ChangeLog

## [v5.0.0] - 2023-01-31

### Changed

- Modified licence for Hibridon 5.0: [`GPLv3`](https://www.gnu.org/licenses/gpl-3.0.en.html) licence
- Added version control using [`GIT`](https://git-scm.com), source code hosted on [GitHub](www.github.com/hibridon/hibridon)
- Replaced custom build system with `CMake` ([#8](https://github.com/hibridon/hibridon/issues/8))
- Changed `fortran` source code from fixed source form to free source form
- Code is now using `fortran` features up to the 2008 standard
- Removed previous custom preprocessing (`ftconv`) system and replaced with `fpp` (Fortran PreProcessor)
- Added interactive command `testpot` that replaces the `testpot` program
- Made the `v2` array (angular coupling matrices) growable to replace its arbitrary hardcoded size - This reduces Hibridon memory footprint ([#51](https://github.com/hibridon/hibridon/issues/51))
- Replaced hardcoded `kmax` value with a runtime argument (`-k <KMAX>`) - This allows to control Hibridon memory footprint without recompiling ([#44](https://github.com/hibridon/hibridon/issues/44))
- Updated most common blocks into modules ([#7](https://github.com/hibridon/hibridon/issues/7))
- Changed informations displayed in headers of the log files
- Made quantum state variables local instead of global ([#201](https://github.com/hibridon/hibridon/issues/201))

### Added

- Added automatic testing of Hibridon ([#160](https://github.com/hibridon/hibridon/issues/160))
- Added continuous integration - This validates the code for MacOS and Linux using either gfortran or ifort ([`c6d06c0`](https://github.com/hibridon/hibridon/commit/c6d06c0), [`b56c086`](https://github.com/hibridon/hibridon/commit/b56c086))
- Added the `-c <COMFILE>` runtime argument to specify a command file ([`32b6ada`](https://github.com/hibridon/hibridon/commit/32b6ada))

### Fixed

- Fixed bug to allow proper handling of the `WFU` files bigger than 2GB (!!! WFU file format has changed !!!) ([#110](https://github.com/hibridon/hibridon/issues/110))

## [v4.4.0] - 2021-03-19

- **2013-11-01** : Added routines ofr collision of `1D/3P` atom with atom
- **2013-10-01** : Symmetric top basis routines (`bastp`, `bastp1`) now have capability to molecules with equivalent spin 1 nuclei, e.g. CD3 and ND3
- **2013-06-01** : Added built-in 2pi--1sigma molecule-molecule basis, `ibasty=20` (`hiba1sg`, `hiba2pi1sg`, `hisystem` also modified)
- **2013-06-01** : Two functions, `is_twomol()` and `is_j12()` added in `hiba1sg`, to check if a basis type (`ibasty`) is a molecule-molecule basis, and if it requires a `j12` array to complete the description of the channel basis.  Statements checking these in various files modified accordingly. (So that when adding a new basis various `I/O` and cross sections routines do not need to be changed.)
- **2013-06-01** : `Hibound` routine rewritten, now uses dynamic allocation, and  provides a simple wave vector analysis when `wavefl=true`.  It calculates only (-`tolai`) bound states if `tolai<0`.
- **2013-06-01** : `Twomol` is now passed to the basis routine.
- **2013-06-01** : Routines `tf3jm0`, `tf3j`, `tj6j` and `tf9j` routines added to `hituil`. These functions take (integer) `2j`, `2m` as parameters.  The `3j` and `6j` functions are faster than the `xf` versions (because conversion between real  number and integers is not fast and not so reliable).
- **2013-06-01** : Eliminate superfluous `nlami` array in basis routines `baastp`, `bach2x`, `basyp`, `bastp1`, `bastpln`
- **2013-06-01** : Dynamically allocate memory for arrays etc. in `hitrnprt` module
- **2013-06-01** : Dynamically allocate memory for arrays etc. in hiprsbr module and added sinqr routine to `hiiolib_f` (to determine max `jtot` and max channel length)
- **2013-04-01** : Added `prsbr` module to compute pressure broadening cross sections
- **2013-03-01** : Revised `bah3p` basis routine for the calculation of cross sections. Shifted zero of energy so that lowest level has `eint = 0`.
- **2013-03-01** : Revised `bah2p` basis routine for the calculation of cross sections. Shifted zero of energy so that lowest level has `eint = 0`.
- **2012-11-01** : Added common block `conlamp` to avoid duplicate use of `lamnum`
- **2012-11-01** : Increased size of `jmx` to 90m and sizes of `lbufs` and `lbuflb`, in partens
- **2012-11-01** : Rewrote `hitrnprt.f` so that all transport cross sections computed with a single command
- **2012-09-01** : Increased size of `jout` array to 50
- **2012-09-01** : Increased maximum size of asymmetric top basis expansion
- **2012-09-01** : Took care of problem with `ibasty=4` and 19 basis types in  `sigk` subroutine in `hitensor`
- **2012-09-01** : Increased lengths of buffers in `partens` common block
- **2012-04-01** : Renomed duplicate variable definitions in `hisystem`
- **2012-04-01** : `pot_nh3h2_2009` no longer writes to units that are not opened.
- **2012-04-01** : Corrected the subroutine in `hisystem.f` that writes input files for  `bastpln`
- **2012-04-01** : New interface to `eadiab`.  The original interface can be accessed  from the flux command.
- **2012-04-01** : Allow writing only adiabatic energies to `wfu` files when `WRSMAT` set `false`. Backward compatibility is lost since `WRSMAT` has to be set `true` for flux/psi calculations. A couple of input files in the test have been modified.
- **2012-04-01** : Simplify the `del` function in `hiutil.f` (used by `f3j0` and `f6j`).   This significantly improves the performance of the bastpln  basis subroutine.
- **2012-04-01** : Removed `cf9j` in `hibastpln.f` which calculates 9-j symbols for  integers but take double precision numbers as input.   The function is never called and can be substituted with  existing f9j and xf9j functions.
- **2012-04-01** : Moved `f9j` routine from `hibastpln.f` to `hiutil.f`
- **2012-04-01** : Fixed bug in `wrhead` subroutine (corrected for alignment)
- **2012-04-01** : Revised for stream I/O of wfu files
- **2012-02-01** : Implemented calculation of differential cross sections for molecule-molecule collisions (`ibasty=9`)
- **2012-01-01** : Moved `dlmo` and `pm0` routines from `hibrid3.f` to `hiutil.f`
- **2012-01-01** : Fixed problem of calculation of elastic cross sections for molecule-molecule collisions with `INTCRS`
- **2012-01-01** : `J12` now printed out for molecule-molecule collisions when `PRINTS` called
- **2012-01-01** : Replaced buffering s-matrix i/o routines with streaming i/o routines (hibridon now requires ifort compiler)
- **2012-01-01** : Time given to nearest second, but still works if date changes
- **2012-01-01** : Last line "quit/exit" now optional in a .com file
- **2012-01-01** : Added Q. Ma as a contributor
- **2011-12-01** : Promoted P. Dagdigian from contributor to author
- **2011-12-01** : Increased `nchtop` to 1000 in `wavewr`
- **2011-12-01** : Increase `kmxbas` to 19 to accommodate more basis routines
- **2011-12-01** : Added `basgpi1` basis routine - to determine angular coupling potential for collisional transfer between a 2sigma state  and one or more vibrational levels in a 2pi state (no isolated-molecule mixing of the sigma and pi states)
- **2011-11-01** : Half-integral j's now printed out with `intcrs` command
- **2011-11-01** : Added `bastp1` basis routine - to determine angular coupling potential for collision of symmetric top molecule with no inversion splitting (like CH3)
- **2011-10-01** : Added calculation of tensor differential cross sections in the geometric apse frame
- **2011-08-01** : Added module to compute transport cross sections
- **2011-08-01** : Included inversion splitting is symmetric top energies and changed mothod of selecting rotational levels (`jmax`, `emax`) in `bastpln` basis routine
- **2011-08-01** : Included inversion splitting is symmetric top energies and changed mothod of selecting rotational levels (`jmax`, `emax`) in `bastp` basis routine
- **2011-05-01** : Added module to compute cross sections for mixed singlet-triple states, e.g. CH3 X3B1 - a1A1 collisions
- **2011-03-01** : Added calculation of tensor differential cross sections in the helicity frame
- **2011-01-01** : Modified method of calculation of hyperfine-resolvec cross sections
- **2010-12-01** : Commented out diagnostic print statement in `diffcrs`
- **2010-12-01** : Added `bach2x` basis routine - to determine angular coupling potential for collision of CH2(X3B1) in a (0,v2,0) bender level (non-rigid molecule)
- **2010-07-01** : Increased size of `jout` array in common block `cosout`
- **2009-09-01** : Added asymmetric top basis


## [v4.3.7] - 2009-08-01

- **2009-09-01** : Added asymmetric top basis
- **2009-08-01** : Added capability of calculating hyperfine-resolved integral cross sections
- **2009-08-01** : Added F. Lique to contributor list
- **2009-01-01** : Updated and corrected `hitensor.f`, `himain.t` to allow correct determination of tensor cross sections for 2Pi molecules
- **2007-12-17** : allow one more significant figures in printed integral cross sections in `hibrid2.f`
- **2007-12-16** : 
  - initialize `izero` in `wavewr` and `waverd`
  - replace calls to `mxma` in `hibrid3` with `dgemm` under `unix-darwin`
- **2007-12-15** : 
  - let `abstol=1.d-16` instead of `0d0` in call to `dsyevr` in subroutine `potent` 
  - change to `character*12` time strings in `propag`
  - initialize `izero` in `smatrx`
  - replace call to `mxma` with call to `dgemm` in subroutine `steppr`
  - set `abstol = 1.d-16` instead of `0d0` in to call `dsyever` in subroutine `wavevc`
- **2007-12-08** : change format in `hibrid2` to use floating rather than fixed format printing `ics` files
- **2007-12-04** : in `hibrid4` truncate `oldlab` and `oldpot` to 40 characters in `waverd` in `hiiolib_c.f` and `hiiolib_f.f` dimension `iword(32)` in restrt `himain.t` ifort `mtime` changed
- **2007-12-03** :
  - change `hibrid2` to append integral cross section file in `prsg`, `prsgpi` for ifort 
  - change the name of runtest to hibtest removed inclusion of `hiunix.c` in `makepot` scripts
  - change `hibrid5` to allow accurate determination of sum of partial cross sections after restart.
  - add `compare_OSX.com` in tests to allow use of `FileMerge.app` under OSX 
- **2007-12-02** : change `hinput` to truncate `jobname` to 8 characters change `hiiolib` to append sequential, formatted files in `openf` for ifort
- **2007-12-01** : change `hiiolib.f` (both c and fortran) to allow longer output for `JOUT`
- **2007-11-30** : use `hiunix.c` from v 2006.0 from molpro
- **2007-11-30** : use `machines.h` (v. 2006.3 Patch(2006.3)) from molpro
- **2007-11-29** : in `hiutil.f`, change function second, dater for fortran, rather than c, calls
- **2007-11-26** : add `cc_hib` procedure in `~/bin`
- **2007-11-26** : add dummy fehler subroutine to `hiutil.f` for compatibility with molpro2006.3 `c_routines`
- **2007-11-26** :   
  - remove nooptimization flags from `hiversion.t`
  - change `maketar` to remove `/src/himain.f`
  - remove nooptimization flags from `himain.t`.
  - `himain.f` is no longer explicitly compiled in `makeobj`
  - explicitly compile `hiinput.f` and `hiversion.f` with -n in `makeobj`
  - eliminate compilation of `himain.f` in `makeobj`
  - replace `hiunix.c` with combination of `rdabsf.c` and `timing_molpro.c` from molpro2006.3
  - replace `machines.h` with molpro2006.3 version

- **2007-11-25** : eliminate makemachine, `makeconfig`, replace `makefirst`

## [v4.2.1] - 2006-09-01

- **2006-09-01** : added nelib `spline.f` and `seval.f` to `hipotutil`
- **2006-07-25** : `hiba2pi.f` changed so that `ntop=n` always in bound state calculations
- **2006-07-21** : `hiba2pi.f` comment at the beginning of `vlm2pi` changed to indicate reference to correction given by G.C. Corey and M. H. Alexander, JCP 85, 5652 (1986)
- **2006-06-10** :
  - changed to correct cs centrifugal barrier for bound state 1sig (`hiba1sg.f`) and 2sig (`hiba2sg.f`)
  - changed `hibound.f` to write out size of total basis changed `hiiolib` to reset `numax=numin` in bound state problems if CS `flag=true` and if `numin>numax` in input
- **2006-06-09** : added `-lSystemStubs` to `makeconfig` for Darwin OS X 10.4
- **2006-04-11** :
  - changed `makeconfig` to allow for `aix5/xlf8.1` on powerpc_4
  - changed `makeobj` to use cc for `hiunix` rather than xlc
- **2005-10-25** :
  - changed call in `hibound` from `dsygv` to `dsygvd` - much faster.
  - change `kaux3` in `himain.t`
- **2005-10-24** : changed `himain.t`, `hibound`:  error in eigenvalues for bound state problem corrected
- **2005-10-20** : changed `hinput`, `hisystem`, `hiba1sg` to allow a 15th basis type:  2P atom with a possibly heteronuclear diatomic
- **2005-03-15** : wavevectors in angstrom in `hibrid1.f` - conversion factor corrected
- **2004-04-11** : `mnscript` & `vscript` changed to allow four-element machine type from `machine.exe`
- **2004-04-11** : `makeconfig` changed for PA8000/HPUX11.0
- **2004-04-07** : `hibrid3.f` changed so that rles works on aix machines
- **2004-04-06** : small changes to `makehib` and `makepot`
- **2004-04-06** : small format change in `prsgpi` (`hibrid2.f`)
- **2004-04-04** : change `makepot` to eliminate adding `hiiolib.o` in link
- **2004-02-29** : `maketar` changed 
- **2004-02-29** : copyright updated
- **2004-02-24** : replace "call rs" with `dsyev` in `hibound` for `unix-darwin`
- **2004-02-24** : replace remaining `call abort` with `stop`
- **2004-02-24** : change `ftconv.exe` to `ftconv_hib.exe` to avoid conflict with molpro
- **2004-02-24** : all `cstart mac`'s eliminated
- **2004-02-24** : rles changed to call lapack routines `dgetrf` and `dgetrs`
- **2004-02-24** : call to `tqlrat` in `wavevc` replaced by call to `dsyevr`
- **2004-02-23** : eigenvalue/eigenvector call in potent replaced by call to `dsyevr`
- **2004-02-23** : `maketar`, `makeconfig` updated
- **2004-02-23** : `hibound` corrected so that call to lapack `dsygv` is made
- **2004-02-23** : minor typo corrected in `/bin/runtest`
- **2004-02-22** : calls to `smxinv` replaced with call to `syminv` for matrix inversion, this calls lapack `dsytrf` and `dsytri`
- **2003-12-28** : `himain` and himatrix changed to allow call to lapack `dspev` for matrix diagonalization in MacOSX version
- **2003-12-28** : `hibrid1.f` changed to allow call to `dgemm` in dtrans instead of `mxma` in MacOSX version
- **2003-12-28** : `hiiolib.f` changed to ensure that sratch files for cs calculations are closed before they are opened
- **2003-12-27** : `makeconfig` and `makemachine` changed for MacOSX 10.3.2 with xlf 8.1beta


## [v4.1.5] - 2003-12-26

- **2003-04-07** :  `hibrid1` changed to allow calculation of hexapole moment in difcrs
- **2003-04-07** : formatting slightly changed in basis routines
- **2003-04-07** : some html help files updated
- **2003-04-05** : added aroh tests, revised `compare.com` and runtest
- **2003-04-04** : formatting change in hiversion
- **2003-04-03** : makemachine and makeconfig updated for sun-sparc
- **2003-04-03** : indexing error corrected in `hiba2pi.f`, which caused unpredictable errors on sunsparc machines
- **2003-04-03** : several small formatting errors corrected in subroutine intcrs in `hibrid5`
- **2003-02-24** : `hiunix.c` updated to include molpro2002.3 `unix.c` and machines.h

## [v4.1.4] - 2003-02-21

- **2003-02-21** : `hibrid5.f` modified to allow calculation of integral cross sections for `jtot_max` less than `jtot2` used in calculating smt file
- **2002-06-03** : `hibah2p.f` modified to correct errors in CC calculations
- **2001-10-10** : `hibrid2` changed to allow appending for `.xsc` files on hp and sgi routines
- **2001-10-09** : `hibrid2` and `hinput` changed to allow minimum threshold  for printing of integral cross sections
- **2001-10-08** : makeconfig changed to include opt1 for all machines
- **2001-10-08** : maketar brought up-to-date with latest help file locations
- **2001-10-08** : `hibhelp.html`, `copyright.html`, `install.htmli`, `acknow.html` modified; `bah3p.html` and `ba2de.html` added
- **2001-10-07** : `makemachine` updated
- **2001-10-07** : `hkey` and `hkey.aix` changed to reflect mha's current email address
- **2001-10-07** : `makeobj` changed so that `hiiolib` is compiled only at optimization level `FOPT1` and `hinput`, `himain,` and `hiversion` are  compiled with optimization level `NOOPT`
- **2001-10-07** : `maketar` changed to include option for inclusion of html files and to keep `include.list` in `/src/pot`
- **2001-10-07** : `makeconfig` changed to include `-liblapack` on hpmachines
- **2001-10-03** : changes in `hiiolib` and `hinput` to restrict `jobnam` to 6 characters
- **2001-10-01** : `hivector` changes so that calls to lapack routines not made if blas3
- **2001-10-01** : `himatrix` changed so that calls to all blas3 and blas2 subroutines are eliminated if unix-blas3 is set in CONFIG
- **2001-10-01** : `himatrix.f` changed so that in mxma and mxmb cstart unix-blas3 changed to cstart unix-blas3 .and.not.unix-hp800
- **2001-10-01** : prsgpi in `hibrid2.f` changed to re-initialize insize at every passage through subroutine
- **2001-10-01** : `himain.t` changed to force skip to next partial wave when `numax` lowered in cs calculations
- **2001-10-01** : change inline optimization directives for unix-hp in `hiversion.t`, `himain.t` and `hinput.f` to be compatible with f90
- **2001-09-19** : formatting changes in `hibrid5` again to allow output of larger channel indices
- **2001-03-25** : formatting changed in `hibrid5` to allow output of channel indexes >100 in size
- **2000-06-03** : potscript changed to incorporate consistently changes made on ibm (2000-05-18) and hp (1999-02-15)
- **2000-05-31** : `makeconfig` changed to integrate latest changes on ibms and hp
- **2000-05-18** : `ftn_hib` changed so the "compilation done" is echoed at the end
- **2000-05-18** : `runhscript`, `potscript`, `ftscript` changed so that alias rm is correct and other small errors corrected
- **2000-05-18** : tests changed to include differential and steric differential calculations for Ar-NO
- **2000-05-18** : `maketar_full` changed so that repeated common directories are not copied
- **2000-05-17** : `makeconfig` changed so that correct essl library is used for POWER, POWER2, and POWER3 machines correct sgi CONFIG is created for machinetype "unix unix-iris" and "unix unix-i4 unix-iris"
- **2000-05-17** : `hiunix.c` changed to newer molpro version

## [v4.1.3] - 2000-04-13

- **2000-04-13** : `makeconfig` added so that aix CONFIG is created for machinetype "unix unix-ibm" and "unix unix-i4 unix-ibm"
- **1999-10-13** : `hisystem.f`, `hinput.f`, `himain.f`, `hiba1sg.f` changed and `hiba2del.f` added to allow additional system: doublet-delta molecule plus atom
- **1999-09-08** : `src/common/parhlp.t` changed to allow longer names for $basedir
- **1999-09-08** : `makemachine` changed so that -qextname added on aix machines
- **1999-09-07** : `makehib` changed so that `-d` option adds `/opt/langtools/lib/end.o` on HP-UX
- **1999-04-01** : `boundc` added to `hiba2pi` to allow correct cent decoupling diagonal potential for bound state calculation
- **1999-03-31** : flush6 added to hp-unix in `himain.t`
- **1999-02-17** : `hibrid1` changed to complex*16 in ampli and difcrs
- **1999-02-15** : `hscript`, `ftscript`, `vscript`, `potscript`, `mnscript` revised to allow 4 components in machine variable
- **1999-02-15** : `ftconv` updated following molpro98.1, compiles to create `ftconv.exe` and `ftconv_nobatch.exe`
- **1999-02-05** : unix-blas3 added to unix-ibm essl configuration

## [v4.1.2] - 1999-02-04

- **1999-02-04** : -qextname compiler option added to `makeconfig` molpro 98.1 elements added to `makeconfig`
- **1999-02-04** : `arno` steric effect added to `hibrid4`, tests upgraded dxsc test for arn2
- **1997-09-01** : -12/30/98 various undocumented changes

## [v4.1.1] - 1997-08-18

- **1997-08-18** : `tests/bench.com` corrected (`inp=Arn2_dxsec.inp`)
- **1997-08-13** : `hiunix.c` replaced by MOLPRO 96.4 unix.c also `macdecls.h` and `machines.h` added from MOLPRO 96.4
- **1997-08-13** : aix ftn -qextname compiler option added to makeconfig


## [v4.1.0] - 1997-05-18

- **1997-05-18** : `hibrid5` and `hiiolib` changed to allow `nucros`to be transferred in `intchk`, `intpol`, `intpl2`, and `intpl3`
- **1997-05-16** : makehib changed to allow -b and -d options
- **1997-05-16** : tests expanded to include test of restart option and -b
- **1997-05-15** : `himain.t`, `hibrid5.f`, `hiiolib.f` changed to enable -b option in  compiling `tq1`, `tq2`, `tq3` matrices not used unless `wavefl`, `photof = .t.`
- **1997-05-13** : all `hiba*.f` codes changed so that rcut is not used if `boundc=.t.`
- **1997-05-08** : `hibah2p.f` corrected in case 1C calculation of energies
- **1997-05-06** : `himain.t` and `hibrid3.f` changed so that nstep is set for logd integration at first partial wave and then not changed
- **1997-05-06** : `no_opt` made compatible for all platforms in makeconfig
- **1997-05-04** : added `compare.com` to tests, changed runtest
- **1997-05-04** : converted all single-precision constants to dp in hiutil, hitensor, `hibrid1`, `hibrid2`, `hibrid3`, `hibrid4`, `pot_vfit.f`, hiamp.f
- **1997-05-04** : corrected `hivector.f` slightly
- **1907-04-30** : `hivector.f` changed to include idamin and idmin
- **1997-04-30** : release 4.1 created
- **1997-04-29** : `hinput.f` changed to include correct pcod's for bound state
- **1997-04-29** : `hiiolib.f` changed to include gendat adapted for bound state
- **1997-04-29** : runtest changed to include new bound state test
- **1997-04-28** : runtest changed to include correct vfit,n2phetest
- **1997-04-28** : `hicommon.all` changed to include new partens
- **1997-04-24** : `hibrid3` changed so that no loop if nlam=0
- **1997-04-23** : `hibrid3` and flow changed  so that return from bound is  correct
- **1997-04-23** : himain changed so that nchtop = nmax in bound calculation
- **1997-04-23** : `hibound.f` changed completely
- **1997-04-23** : `himatrix.f` expanded to incorporate subroutines `rsg`, `rebak`, and `reduce` from `EISPACK`
- **1997-04-18** : runtest changed to copy `FOLLMEG.BIN` from `tests` to `testnew`
- **1997-04-18** : hinput updated to allow 8 input variables for `tenxsc`
- **1997-04-18** : `maketar` updated to remove `src/*.hold` files
- **1997-04-18** : `hibound.f` changed to give dummy return if not unix-aix machines
- **1997-04-17** : `hitensor` updated
- **1997-04-16** : `pot_vfit` changed to give correct comparison with follmeg
- **1997-04-16** : `hibrid.hlp` changed to reflect following change in flux command
- **1997-04-16** : `hibrid4.f` changed (transmt changed) so that if rout < 0 transformation matrices will be printed out at all R
- **1997-04-08** : `pot_vfit.f` changed (sdot -> ddot, and `vv0`)
- **1997-04-08** : `hibrid3.f` changed to correct `testptn` and `testpt`
- **1997-04-08** : `hisystem.f` changed to eliminate variables n0max and interp from sy2sg
- **1997-04-08** : `hiba2sg.f` changed to accept `lammin=0,` and to eliminate variables `n0max` and `interp`
- **1997-03-17** : `makecommon` updated to automatically create soft link in subdirectory `src/pot`
- **1997-03-17** : updated `himain.t` to include `boundc` and `boundf` as logical parameters
- **1997-03-17** : updated `common/parcode` to include boundc as logical parameter
- **1997-03-17** : updated `common/logdef` and  `common/lpar` to include boundc as logical  parameter
- **1997-03-17** : updated `hiiolib.f` to include `boundc` as logical parameter
- **1997-03-17** : updated `hinput.f` to include `boundc` as logical parameter
- **1997-03-17** : updated `hibrid3.f` to include moonbong's changes (additional print to `resenergy.dat` in `smatop`, now disabled, for use in cnne calculations, susan gregurick's bound state programs (`bound`, `boundwavfn,i` `h_basis`)
- **1997-03-17** : updated `hibrid1.f` to include moonbong's changes (increased dimensions, also change to difcrs)
- **1997-03-14** : updated `hiba1sg.f` to expand format for > 999 channels
- **1997-03-14** : updated `hisystem.f` to include check of system parameters  after they have been changed
- **1997-03-09** :  remove bug in `dcopy` in subroutine potent (this does not affect any computed results)
- **1997-03-05** :  change format of output in 22P calculations in `smatop`

## [v4.0.2] - 1996-11-19

- **1996-11-19** : release 4.0.2 created
- **1996-11-19** : create `hkey.irix64` to accomodate `resolv.conf` moving to `/usr/etc`
- **1996-11-19** : update `/bin/makefirst` to 96.3 to include IRIX64 machine type
- **1996-05-27** : release 4.0.1 created
- **1996-05-27** : update `/bin/makeconfig` to include `-fast unix-noblas` for SunOS
- **1996-05-27** : update `doc/install.html` to include Sun workstations
- **1996-05-27** : update `doc/hibhelp.html` to include Sun workstations
- **1996-05-27** : update `src/hiblas.f` src/`himain.f` src/`hiiolib.f` src/`hivector.f`  `src/himatrix.f` for SunOS
- **1996-05-24** : update `bin/makemachine` to 96.1 for sites where license already accepted
- **1996-05-24** : update `bin/makefirst` to 96.2 to prompt for licence only on first installation
- **1996-05-24** : patch1 tar archive created
- **1996-05-23** : all tests rerun with code updated to include compiler option -qdpc=e
- **1996-05-23** : updated `bin/makeconfig` to 96.3 to include compiler option -qdpc=e for AIX
- **1996-05-22** : `tests/Ch3itest.flx` and `tests/Ch3itest.psi` updated with corrected pot_ch3i.f
- **1996-05-22** : `src/pot/pot_ch3i.f` corrected
- **1996-05-21** : `bin/potcopy` and `bin/maketar` updated to 96.2 to use `/scratch/mha/htar`
- **1996-05-20** : `bin/makepot` updated to 96.2 to include `hiunix.o` (for sgi)
- **1996-05-20** : updated `pot_ch3i.f` to include generic driver
- **1996-05-16** : `bin/makeconfig` upgraded to 96.2 to restrict `-O2`





[v5.0.0]: https://github.com/hibridon/hibridon/releases/tag/v5.0.0
[v4.4.0]: https://github.com/hibridon/hibridon/releases/tag/pdagdigian-20210319
[v4.3.7]: https://github.com/hibridon/hibridon/releases/tag/v4.3.7
[v4.2.1]: https://github.com/hibridon/hibridon/releases/tag/v4.2
[v4.1.5]: https://github.com/hibridon/hibridon/releases/tag/v4.1.5
[v4.1.4]: https://github.com/hibridon/hibridon/releases/tag/missing
[v4.1.3]: https://github.com/hibridon/hibridon/releases/tag/missing
[v4.1.2]: https://github.com/hibridon/hibridon/releases/tag/missing
[v4.1.1]: https://github.com/hibridon/hibridon/releases/tag/missing
[v4.1.0]: https://github.com/hibridon/hibridon/releases/tag/missing
[v4.0.2]: https://github.com/hibridon/hibridon/releases/tag/missing

