# Hibridon ChangeLog

- 31/03/2022 Fixed a bug that caused the `isa` parameter written in the `ics` file to be wrong in most cases

## The following changes have been made in bringing hibridon 4.3 to release level 4.4

- 11/13     Added routines ofr collision of 1D/3P atom with atom

- 10/13     Symmetric top basis routines (bastp, bastp1) now have capability
          to molecules with equivalent spin 1 nuclei, e.g. CD3 and ND3

- 6/13      Added built-in 2pi--1sigma molecule-molecule basis, ibasty=20
          (hiba1sg, hiba2pi1sg, hisystem also modified)

- 6/13      Two functions, is_twomol() and is_j12() added in hiba1sg, to check 
          if a basis type (ibasty) is a moleule-molecule basis, and if it 
          requires a j12 array to complete the description of the channel
          basis.  Statements checking these in various files modified accordingly.  
          (So that when adding a new basis various I/O and cross sections 
          routines do not need to be changed.)

- 6/13      Hibound routine rewritten, now uses dynamic allocation, and 
          provides a simple wave vector analysis when wavefl=true.  It
          calculates only (-tolai) bound states if tolai<0.

- 6/13      Twomol is now passed to the basis routine.

- 6/13      Routines tf3jm0, tf3j, tj6j and tf9j routines added to hituil. These 
          functions take (integer) 2j, 2m as parameters.  The 3j and 6j functions 
          are faster than the xf versions (because conversion between real 
          number and integers is not fast and not so reliable).

- 6/13      Eliminate superfluous nlami array in basis routines baastp, bach2x,
          basyp, bastp1, bastpln

- 6/13      Dynamically allocate memory for arrays etc. in hitrnprt module

- 6/13      Dynamically allocate memory for arrays etc. in hiprsbr module
          and added sinqr routine to hiiolib_f (to determine max jtot and
          max channel length)

- 4/13      Added prsbr module to compute pressure broadening cross sections

- 3/13      Revised bah3p basis routine for the calculation of cross sections.
          Shifted zero of energy so that lowest level has eint = 0.

- 3/13      Revised bah2p basis routine for the calculation of cross sections.
          Shifted zero of energy so that lowest level has eint = 0.

- 11/12     Added common block conlamp to avoid duplicate use of lamnum

- 11/12     Increased size of jmx to 90m and sizes of lbufs and lbuflb, in partens

- 11/12     Rewrote hitrnprt.f so that all transport cross sections computed
          with a single command

- 9/12      Increased size of jout array to 50

- 9/12      Increased maximum size of asymmetric top basis expansion

- 9/12      Took care of problem with ibasty=4 and 19 basis types in 
          sigk subroutine in hitensor

- 9/12      Increased lengths of buffers in partens common block

- 4/12      Renomed duplicate variable definitions in hisystem

- 4/12      pot_nh3h2_2009 no longer writes to units that are not opened.

- 4/12      Corrected the subroutine in hisystem.f that writes input files for 
          bastpln

- 4/12      New interface to eadiab.  The original interface can be accessed 
          from the flux command.

- 4/12      Allow writing only adiabatic energies to wfu files when WRSMAT set 
          false.  Backward compatibility is lost since WRSMAT has to be set 
          true for flux/psi calculations.  A couple of input files in the 
          test have been modified.

- 4/12      Simplify the del function in hiutil.f (used by f3j0 and f6j).  
          This significantly improves the performance of the bastpln 
          basis subroutine.

- 4/12      Removed cf9j in hibastpln.f which calculates 9-j symbols for 
          integers but take double precision numbers as input.  
          The function is never called and can be substituted with 
          existing f9j and xf9j functions.

- 4/12      Moved f9j routine from hibastpln.f to hiutil.f

- 4/12      Fixed bug in wrhead subroutine (corrected for alignment)

- 4/12      Revised for stream I/O of wfu files

- 2/12      Implemented calculation of differential cross sections
          for molecule-molecule collisions (ibasty=9)

- 1/12      Moved dlmo and pm0 routines from hibrid3.f to hiutil.f

- 1/12      Fixed problem of calculation of elastic cross sections
          for molecule-molecule collisions with INTCRS

- 1/12      J12 now printed out for molecule-molecule collisions
          when PRINTS called

- 1/12      Replaced buffering s-matrix i/o routines with
          streaming i/o routines (hibridon now requires ifort compiler)

- 1/12      Time given to nearest second, but still works if date changes

- 1/12      Last line "quit/exit" now optional in a .com file

- 1/12      Added Q. Ma aa a contributor

- 12/11     Promoted P. Dagdigian from contributor to author

- 12/11     Increased nchtop to 1000 in wavewr

- 12/11     Increase kmxbas to 19 to accommodate more basis routines

- 12/11     Added basgpi1 basis routine - to determine angular coupling
          potential for collisional transfer between a 2sigma state 
          and one or more vibrational levels in a 2pi state (no
          isolated-molecule mixing of the sigma and pi states)

- 11/11     Half-integral j's now printed out with intcrs command

- 11/11     Added bastp1 basis routine - to determine angular coupling
          potential for collision of symmetric top molecule with no
          inversion splitting (like CH3)

- 10/11     Added calculation of tensor differential cross sections
          in the geometric apse frame

- 08/11     Added module to compute transport cross sections

- 08/11     Included inversion splitting is symmetric top energies and
          changed mothod of selecting rotational levels (jmax, emax)
          in bastpln basis routine

- 08/11     Included inversion splitting is symmetric top energies and
          changed mothod of selecting rotational levels (jmax, emax)
          in bastp basis routine

- 05/11     Added module to compute cross sections for mixed singlet-triple
          states, e.g. CH3 X3B1 - a1A1 collisions

- 03/11     Added calculation of tensor differential cross sections
          in the helicity frame

- 01/11     Modified method of calculation of hyperfine-resolvec cross sections

- 12/10     Commented out diagnostic print statement in diffcrs

- 12/10     Added bach2x basis routine - to determine angular coupling
          potential for collision of CH2(X3B1) in a (0,v2,0) bender
          level (non-rigid molecule)

- 07/10     Increased size of jout array in common block cosout

- 09/09     Added asymmetric top basis




## to migrate to release 4.4, you should download the upgraded tar file for the complete package from the website

[http://www2.chem.umd.edu/groups/alexander/hibridon/hib44](http://www2.chem.umd.edu/groups/alexander/hibridon/hib44)

and follow the installation instructions contained on the webpage
     www2.chem.umd.edu/groups/alexander/hibridon/hib44/install.html

to retain a copy of any previous versions, you may wish to
carry out the intall in a new main hibridon directory.

the current release level has been checked for the following
operating system and compiler environments:

- MacOSX v. 10.7.x
  - ifort compiler  composer_xe_2011
  - MKL:  at least 10.2


## The following changes have been made in bringing hibridon 4.2.1 to release level 4.3

- 09/09     Added asymmetric top basis

- 08/09     Added capability of calculating hyperfine-resolved integral cross sections
          Added F. Lique to contributor list

- 12/08-1/09 Updated and corrected hitensor.f, himain.t to allow correct determination of tensor
           cross sections for 2Pi molecules

- 12/17/07  allow one more significant figures in printed integral cross sections in hibrid2.f

- 12/16/07  initialize izero in wavewr and waverd
          replace calls to mxma in hibrid3 with dgemm under unix-darwin

- 12/15/07  let abstol=1.d-16 instead of 0d0 in call to dsyevr in subroutine potent
          change to character*12 time strings in propag
          initialize izero in smatrx
          replace call to mxma with call to dgemm in subroutine steppr
          set abstol = 1.d-16 instead of 0d0 in to call dsyever in subroutine wavevc

- 12/8/07   change format in hibrid2 to use floating rather than fixed format printing ics files

- 12/4/07   in hibrid4 truncate oldlab and oldpot to 40 characters in waverd
          in hiiolib_c.f and hiiolib_f.f dimension iword(32) in restrt
          himain.t ifort mtime changed

- 12/3/07   change hibrid2 to append integral cross section file in prsg, prsgpi for ifort
          change the name of runtest to hibtest
          removed inclusion of hiunix.c in makepot scripts
          change hibrid5 to allow accurate determination of sum of partial cross sections after restart
          add compare_OSX.com in tests to allow use of FileMerge.app under OSX 

- 12/2/07   change hinput to truncate jobname to 8 characters
          change hiiolib to append sequential, formatted files in openf for ifort

- 12/1/07   change hiiolib.f (both c and fortran) to allow longer output for JOUT

- 11/30/07  use hiunix.c from v 2006.0 from molpro

- 11/30/07  use machines.h (v. 2006.3 Patch(2006.3)) from molpro

- 11/29/07  in hiutil.f, change function second, dater for fortran, rather than c, calls

- 11/26/07  add cc_hib procedure in ~/bin

- 11/26/07  add dummy fehler subroutine to hiutil.f for compatibility with molpro2006.3 c_routines

- 11/26/07:  
    - remove nooptimization flags from hiversion.t
    - change maketar to remove /src/himain.f
    - remove nooptimization flags from himain.t.
    - himain.f is no longer explicitly compiled in makeobj
    - explicitly compile hiinput.f and hiversion.f with -n in makeobj
    - eliminate compilation of himain.f in makeobj
    - replace hiunix.c with combination of rdabsf.c and timing_molpro.c from molpro2006.3
    - replace machines.h with molpro2006.3 version

- 11/25/07  eliminate makemachine, makeconfig, replace makefirst

## to migrate to release 4.3, you should download the upgraded tar file for the omplete package from the website

[www.chem.umd.edu/groups/alexander/hibridon/hib43](www.chem.umd.edu/groups/alexander/hibridon/hib43)

and follow the installation instructions contained on the webpage
     www.chem.umd.edu/groups/alexander/hibridon/hib43/install.html

to retain a copy of any previous versions, you may wish to
carry out the intall in a new main hibridon directory.

the current release level has been checked for the following operating system and compiler environments:
- MacOSX v. 10.5.x
  - ifort compiler 10.1.012, 11.0.064, 11.1.058
  - MKL:  10.0.2.018, 11.0.064, 11.1.058
- MacOSX v. 10.4.11
  - ifort compiler 9.1.041
  - xlf compiler 8.1

not yet fully operational on the following machines
- AIX v. 4.2, 4.3
  - xlf compiler 3.2.3, 3.2.4, 7.1
- HP-UX v. 9.01, 9.02, 9.03, 10.01, 10.20, 11.0
- IRIX 6.5 f77v7.1
- OSF v3.2
- SunOS v2 rel. 4.1.4, 5.5.1, 5.7 (f77 v6.2)

## The following changes have been made in bringing hibridon 4.1.5 to relase level 4.2.1:

- 9/1/06    added nelib spline.f and seval.f to hipotutil

- 7/25/06   hiba2pi.f changed so that ntop=n always in bound state calculations

- 7/21/06   hiba2pi.f comment at the beginning of vlm2pi changed to indicate
          reference to correction given by G.C. Corey and M. H. Alexander, JCP 85, 5652 (1986)

- 6/10/06   changed to correct cs centrifugal barrier for bound state
          1sig (hiba1sg.f) and 2sig (hiba2sg.f)
          changed hibound.f to write out size of total basis
          changed hiiolib to reset numax=numin in bound state problems
             if CS flag=true and if numin>numax in input

- 6/9/06    added -lSystemStubs to makeconfig for Darwin OS X 10.4

- 4/11/06   changed makeconfig to allow for aix5/xlf8.1 on powerpc_4
          changed makeobj to use cc for hiunix rather than xlc

- 10/25/05  changed call in hibound from dsygv to dsygvd - much faster.
          change kaux3 in himain.t

- 10/24/05  changed himain.t, hibound:  error in eigenvalues for bound state problem
          corrected

- 10/20/05  changed hinput, hisystem, hiba1sg to allow a 15th basis type:  
          2P atom with a possibly heteronuclear diatomic

- 3/15/05   wavevectors in angstrom in hibrid1.f - conversion factor corrected

- 4/11/04   mnscript & vscript changed to allow four-element machine type from 
          machine.exe

- 4/11/04   makeconfig changed for PA8000/HPUX11.0

- 4/7/04    hibrid3.f changed so that rles works on aix machines

- 4/6/04    small changes to makehib and makepot

- 4/6/04    small format change in prsgpi (hibrid2.f)

- 4/4/04    change makepot to eliminate adding hiiolib.o in link

- 2/29/04   maketar changed 

- 2/29/04   copyright updated

- 2/24/04   replace "call rs" with dsyev in hibound for unix-darwin

- 2/24/04   replace remaining "call abort" with "stop"

- 2/24/04   change ftconv.exe to ftconv_hib.exe to avoid conflict with
          molpro

- 2/24/04   all "cstart mac"'s eliminated

- 2/24/04   rles changed to call lapack routines dgetrf and dgetrs

- 2/24/04   call to tqlrat in wavevc replaced by call to dsyevr

- 2/23/04   eigenvalue/eigenvector call in potent replaced by call to dsyevr

- 2/23/04   maketar, makeconfig updated

- 2/23/04   hibound corrected so that call to lapack dsygv is made

- 2/23/04   minor typo corrected in /bin/runtest

- 2/22/04   calls to smxinv replaced with call to syminv for matrix
          inversion, this calls lapack dsytrf and dsytri

- 12/28/03  himain and himatrix changed to allow call to lapack dspev for
          matrix diagonalization in MacOSX version

- 12/28/03  hibrid1.f changed to allow call to dgemm in dtrans instead
          of mxma in MacOSX version

- 12/28/03  hiiolib.f changed to ensure that sratch files for cs
          calculations are closed before they are opened

- 12/27/03  makeconfig and makemachine changed for MacOSX 10.3.2 with xlf 8.1beta

## to migrate to release 4.2, you should download the upgraded tar file for the complete package from the website
[www.chem.umd.edu/physical/alexander/hibridon](www.chem.umd.edu/physical/alexander/hibridon)    

and follow the installation instructions contained on the webpage [www.chem.umd.edu/physical/alexander/hibridon/install.html](www.chem.umd.edu/physical/alexander/hibridon/install.html)

to retain a copy of any previous versions, you may wish to
carry out the intall in a new main hibridon directory.

the current release level has been checked for the following
operating system and compiler environments:

- MacOSX v. 10.3.2
  - xlf compiler 8.1

not yet fully operational on the following machines
- AIX v. 4.2, 4.3
  - xlf compiler 3.2.3, 3.2.4, 7.1
- HP-UX v. 9.01, 9.02, 9.03, 10.01, 10.20, 11.0
- IRIX 6.5 f77v7.1
- OSF v3.2
- SunOS v2 rel. 4.1.4, 5.5.1, 5.7 (f77 v6.2)