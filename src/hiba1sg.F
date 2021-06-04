*************************************************************************
*                                                                       *
*                            basis  library                             *
*                                                                       *
*************************************************************************
*                          routines included:                           *
*                                                                       *
*  1. basis      dispatcher routine to select basis                     *
*                see in header for basis subroutine for list of basis   *
*                routines currently available                           *
*  2. ba1sg      basis for singlet sigma scattering                     *
*  3. vlm1sg     calculates v-lamda matrices for above                  *
*  4. is_twomol check if a basis is for mol-mol collision (j=10j1+j2)   *
*  5. is_j12    check if j12 is used in a basis                         *
*                                                                       *
*     current revision:  24-jul-2019 (p.dagdigian)                       *
*                                                                       *
*************************************************************************
* --------------------------------------------------------------------
      subroutine basis (j, l, is, jhold, ehold, ishold, nlevel, nlevop,
     :                  sc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu,
     :                  csflag, clist, bastst, ihomo, nu, numin, jlpar,
     :                  twomol, n, nmax, ntop)
* --------------------------------------------------------------------
*
*  dispatcher to select the correct basis routine for given problem
*  the variable in common block /coselb/ ibasty determines the
*  basis used in the current calculation
*
*  basis routines currently available are:
*  ibasty         basis routine          kind of problem
*    1              ba1sg             singlet sigma scattering
*    2              ba2sg             doublet sigma scattering
*    3              ba2pi             doublet pi scattering
*    4              basgpi            sigma/pi scattering
*    5              bapi              general pi scattering
*    6              bastp             symmetric top scattering
*    7              ba13p             1/3 P atom scattering
*    8              ba2mol            1sigma + 1sigma
*    9              bastpln           symmetric top + linear molecule
*    10             ba22p             2/2 P atom scattering
*    11             ba1del            singlet delta scattering
*    12             bah2p             homonuclear + 2P atom
*    13             bah3p             homonuclear + 3P atom
*    14             ba2del            doublet delta scattering
*    15             badiat2p          heteronuclear + 2P atom  **to check
*    16             baastp            asymmetric top scattering
*    17             bach2x            CH2(X B1) (0,v2,0) bender levels
*    18             bastp1            symmetric top - no inversion doubling
*    19             basgpi1           2sigma | 2pi + atom (no pertubations)
*    20             ba2pi1sg          2pi molecule + 1sigma molecules
*    21             bastp1sg          symmetric top + 1sigma molecules
*    22             ba1d3p            1D/3P atom + closed-shell atom
*    23             ba3p2s            3P atom + 2S atom
*    24             basphtp           spherical top + atom scattering
*    25             ba1s1sg           1sigma + 1sigma (different molecules)
*    26             ba2sg1sg          2sigma + 1sigma molecules
*    27             baastp1           C2v asymmetric top scattering, w body-frame
*                                     quant axis along C2 axis (compatible w MOLSCAT)
*    28             ba3sg1sg          3sigma + 1sigma molecules
*    29             baastp2           chiral asymmetric top + atom scattering
*    30             baastp3           C2v asymmetric top + linear molecule
*    99 or higher   bausr             user defined basis
*  author: b. follmeg
*  current revision of list:  20-jun-2019 (p.dagdigian)
* --------------------------------------------------------------------
*  variables in call list:
*    j:        on return contains rotational quantum numbers for each
*              channel
*    l:        on return contains orbital angular momentum for each
*              channel
*    is:       on return contains symmetry index for each channel
*    jhold:    on return contains rotational quantum numbers for each
*              rotational level
*    ehold:    on return contains energy in hartrees of each rotational
*              level
*    ishold:   on return contains symmetry index of each rotational level
*    nlevel:   on return contains number of energetically distinct
*              rotational levels used in channel basis
*    nlevop:   on return contains number of energetically distinct
*              rotational levels used in channel basis which are open
*              asymptotically
*    sc1,sc2:  scratch vectors of length at least nmax
*    sc3,sc4:  scratch vectors of length at least nmax
*              these scratch vectors are not used here
*    rcut:     cut-off point for keeping higher energy channels
*              if any open channel is still closed at r=rcut, then
*              all closed channels as well any open channels which are
*              still closed at r=rcut are dropped from basis
*              note that this is ignored in molecule-surface collisions!!!
*    jtot:     total angular momentum
*    flaghf:   if .true., then system with half-integer spin
*              if .false., then system with integer spin (this is the case
*              here)
*    flagsu:   if .true., then molecule-surface collisons
*    csflag:   if .true., then coupled-states calculation
*              if .false., the close-coupled calculation
*    clist:    if .true., then quantum numbers and energies listed for
*              each channel
*    bastst:   if .true., then execution terminates after the first call
*              to basis
*    ihomo:    if .true. , then homonuclear molecule
*              if .false., then heteronuclear molecule
*              if the molecule is homonuclear (ihomo = .true.), the
*              rotational levels included go from jmin to jmax in steps
*              of 2 and only even lambda terms in the anisotropy are
*              included
*    nu:       coupled-states projection index
*    numin:    minimum coupled states projection index
*              for cc calculations nu and numin are both set = 0 by calling
*              program
*    jlpar:    total parity of included channels in cc calculation
*              parity=(-1)**jlpar (by definition parity=(-1)**(j+l+jtot) )
*    n:        on return equals number of channels
*    nmax:     maximum dimension of arrays
*    ntop:     maximum row dimension of all matrices passed to subroutines
*              propag and soutpt.  ntop is set in basis only if nu = numin
*              otherwise it is unchanged from the value supplied by the
*              calling program
*    note!!!   if flaghf = .true., then the true values of the rotational
*    quantum numbers, the total angular momentum, and the coupled-states
*    projection index are equal to the values stored in j, jtot, and nu
*    plus 1/2
      implicit double precision (a-h,o-z)
      integer j, l, is, jhold, ishold, nlevel, nlevop, jtot, nu,
     :        jlpar, n, nmax
*      real ehold, sc1, sc2, sc3, sc4, rcut
      logical flaghf, flagsu, csflag, clist, bastst, ihomo, twomol
      integer ibasty
      common /coselb/ ibasty
      dimension j(1), l(1), is(1), jhold(1), ehold(1), ishold(1),
     :         sc1(1), sc2(1), sc3(1), sc4(1)
*
*  select basis routine according to value of ibasty
      if (ibasty .ge. 99) then
*  user supplied routine
        call bausr(j, l, is, jhold, ehold, ishold, nlevel, nlevop,
     :                  sc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu,
     :                  csflag, clist, bastst, ihomo, nu, numin, jlpar,
     :                  n, nmax, ntop)
        return
      endif
      goto (100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,
     :      1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,
     :      2500,2600,2700,2800,2900,3000)
     :      ibasty
*  singlet sigma basis
100   call ba1sg (j, l, is, jhold, ehold, ishold, nlevel, nlevop,
     :                  sc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu,
     :                  csflag, clist, bastst, ihomo, nu, numin, jlpar,
     :                  n, nmax, ntop)
      return
*  doublet sigma basis
200   call ba2sg (j, l, is, jhold, ehold, ishold, nlevel,
     :                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :                  flaghf, flagsu, csflag, clist, bastst,
     :                  ihomo, nu, numin, jlpar, n, nmax, ntop)
      return
*  doublet pi basis
300   call ba2pi (j, l, is, jhold, ehold, ishold, nlevel,
     :                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :                  flaghf, flagsu, csflag, clist, bastst,
     :                  ihomo, nu, numin, jlpar, n, nmax, ntop)
      return
*  sigma/pi basis
400   call basgpi (j, l, is, jhold, ehold, ishold, nlevel,
     :                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :                  flaghf, flagsu, csflag, clist, bastst,
     :                  ihomo, nu, numin, jlpar, n, nmax, ntop)
      return
*   general pi basis
500    call bapi (j, l, is, jhold, ehold, ishold, nlevel,
     :                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :                  flaghf, flagsu, csflag, clist, bastst,
     :                  ihomo, nu, numin, jlpar, n, nmax, ntop)
      return
*  symmetric top basis, with inversion doubling
600   call bastp (j, l, is, jhold, ehold, ishold, nlevel,
     :                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :                  flaghf, flagsu, csflag, clist, bastst,
     :                  ihomo, nu, numin, jlpar, n, nmax, ntop)
      return
*  1/3 P atom basis
700   call ba13p (j, l, is, jhold, ehold, ishold, nlevel,
     :                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :                  flaghf, flagsu, csflag, clist, bastst,
     :                  ihomo, nu, numin, jlpar, n, nmax, ntop)
      return
* 1sigma + 1sigma basis
800   call ba2mol (j, l, is, jhold, ehold, ishold, nlevel,
     :                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :                  flaghf, flagsu, csflag, clist, bastst,
     :                  ihomo, nu, numin, jlpar, n, nmax, ntop)
      return
* symmetric top + 1 sigma basis
900   call bastpln (j, l, is, jhold, ehold, ishold, nlevel,
     :                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :                  flaghf, flagsu, csflag, clist, bastst,
     :                  ihomo, nu, numin, jlpar, n, nmax, ntop)
      return
* 2/2 P atom basis
1000   call ba22p (j, l, is, jhold, ehold, ishold, nlevel,
     :                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :                  flaghf, flagsu, csflag, clist, bastst,
     :                  ihomo, nu, numin, jlpar, n, nmax, ntop)
      return
*  singlet delta basis
1100  call ba1del (j, l, is, jhold, ehold, ishold, nlevel,
     :                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :                  flaghf, flagsu, csflag, clist, bastst,
     :                  ihomo, nu, numin, jlpar, n, nmax, ntop)
      return
*  homonuclear + 2P atom basis
1200  call bah2p (j, l, is, jhold, ehold, ishold, nlevel,
     :                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :                  flaghf, flagsu, csflag, clist, bastst,
     :                  ihomo, nu, numin, jlpar, n, nmax, ntop)
      return
*  homonuclear + 3P atom basis
1300  call bah3p (j, l, is, jhold, ehold, ishold, nlevel,
     :                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :                  flaghf, flagsu, csflag, clist, bastst,
     :                  ihomo, nu, numin, jlpar, n, nmax, ntop)
      return
*  doublet delta basis
1400  call ba2del (j, l, is, jhold, ehold, ishold, nlevel,
     :                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :                  flaghf, flagsu, csflag, clist, bastst,
     :                  ihomo, nu, numin, jlpar, n, nmax, ntop)
      return
*  heteronuclear + 2P atom basis
c1500  call bah2p (j, l, is, jhold, ehold, ishold, nlevel,
c    :                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
c    :                  flaghf, flagsu, csflag, clist, bastst,
c    :                  ihomo, nu, numin, jlpar, n, nmax, ntop)
1500  call badiat2p (j, l, is, jhold, ehold, ishold, nlevel,
     :                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :                  flaghf, flagsu, csflag, clist, bastst,
     :                  ihomo, nu, numin, jlpar, n, nmax, ntop)
*  asymmetric top basis
1600   call baastp (j, l, is, jhold, ehold, ishold, nlevel,
     :                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :                  flaghf, flagsu, csflag, clist, bastst,
     :                  ihomo, nu, numin, jlpar, n, nmax, ntop)
      return
*  CH2(X 3B1) (0,v2,0) bender levels
1700   call bach2x (j, l, is, jhold, ehold, ishold, nlevel,
     :                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :                  flaghf, flagsu, csflag, clist, bastst,
     :                  ihomo, nu, numin, jlpar, n, nmax, ntop)
      return
*  symmetric top basis, with no inversion doubling
1800   call bastp1 (j, l, is, jhold, ehold, ishold, nlevel,
     :                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :                  flaghf, flagsu, csflag, clist, bastst,
     :                  ihomo, nu, numin, jlpar, n, nmax, ntop)
      return
*  2sig-2pi + atom scattering (one 2sigma state and one or more 2pi
*   vibrational levels, no sigma-pi spectroscopic perturbations)
1900   call basgpi1 (j, l, is, jhold, ehold, ishold, nlevel,
     :                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :                  flaghf, flagsu, csflag, clist, bastst,
     :                  ihomo, nu, numin, jlpar, n, nmax, ntop)
      return
*  2pi + 1sigma molecules
 2000 call ba2pi1sg(j, l, is, jhold, ehold, ishold, nlevel,
     $     nlevop, rcut, jtot, flaghf, flagsu, csflag, clist,
     $     bastst, ihomo, nu, numin, jlpar, twomol, n, nmax, ntop)
      return
*  symmetric top + 1sigma molecules
 2100 call bastp1sg(j, l, is, jhold, ehold, ishold, nlevel,
     $     nlevop, rcut, jtot, flaghf, flagsu, csflag, clist,
     $     bastst, ihomo, nu, numin, jlpar, twomol, n, nmax, ntop)
      return
*  1D/3P atom + closed-shell atom
 2200 call ba1d3p(j, l, is, jhold, ehold, ishold, nlevel,
     $     nlevop, rcut, jtot, flaghf, flagsu, csflag, clist,
     $     bastst, ihomo, nu, numin, jlpar, n, nmax, ntop)
      return
*   3P atom + 2S atom
 2300 call ba3p2s (j, l, is, jhold, ehold, ishold, nlevel,
     :     nlevop, rcut, jtot, flaghf, flagsu, csflag, clist,
     :     bastst, ihomo, nu, numin, jlpar, n, nmax, ntop)
      return
*   spherical top + atom
 2400 call basphtp (j, l, is, jhold, ehold, ishold, nlevel,
     :     nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :     flaghf, flagsu, csflag, clist, bastst, ihomo,
     :     nu, numin, jlpar, n, nmax, ntop)
      return
*   1sigma + 1sigma basis (different molecules)
 2500 call ba1s1sg (j, l, is, jhold, ehold, ishold, nlevel,
     :     nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :     flaghf, flagsu, csflag, clist, bastst, ihomo,
     :     nu, numin, jlpar, n, nmax, ntop)
      return
*   2sigma + 1sigma basis
 2600 call ba2s1sg (j, l, is, jhold, ehold, ishold, nlevel,
     :     nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :     flaghf, flagsu, csflag, clist, bastst, ihomo,
     :     nu, numin, jlpar, n, nmax, ntop)
      return
*   C2v asymmetric top + atom scattering
 2700 call baastp1 (j, l, is, jhold, ehold, ishold, nlevel,
     :     nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :     flaghf, flagsu, csflag, clist, bastst, ihomo,
     :     nu, numin, jlpar, n, nmax, ntop)
      return
*   3sigma + 1sigma basis
 2800 call ba3s1sg (j, l, is, jhold, ehold, ishold, nlevel,
     :     nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :     flaghf, flagsu, csflag, clist, bastst, ihomo,
     :     nu, numin, jlpar, n, nmax, ntop)
      return
*   chiral asymmetric top + atom scattering
 2900 call baastp2 (j, l, is, jhold, ehold, ishold, nlevel,
     :     nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :     flaghf, flagsu, csflag, clist, bastst, ihomo,
     :     nu, numin, jlpar, n, nmax, ntop)
      return
*   C2v asymmetric top + linear molecule scattering
 3000 call baastp3 (j, l, is, jhold, ehold, ishold, nlevel,
     :     nlevop, sc1, sc2, sc3, sc4, rcut, jtot,
     :     flaghf, flagsu, csflag, clist, bastst, ihomo,
     :     nu, numin, jlpar, n, nmax, ntop)
      return
      end
*--------------------------------------------------------------------
      subroutine ba1sg (j, l, is, jhold, ehold, ishold, nlevel, nlevop,
     :                  sc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu,
     :                  csflag, clist, bastst, ihomo, nu, numin, jlpar,
     :                  n, nmax, ntop)
* --------------------------------------------------------------------
*  subroutine to determine angular coupling potential
*  for collision of singlet-sigma molecule with a structureless atom
*  authors:  millard alexander and hans-joachim werner
*  current revision date:  18-jun-2006 by mha
* --------------------------------------------------------------------
*  variables in call list:
*    j:        on return contains rotational quantum numbers for each
*              channel
*    l:        on return contains orbital angular momentum for each
*              channel
*    is:       on return contains symmetry index for each channel
*    jhold:    on return contains rotational quantum numbers for each
*              rotational level
*    ehold:    on return contains energy in hartrees of each rotational
*              level
*    ishold:   on return contains symmetry index of each rotational level
*    nlevel:   on return contains number of energetically distinct
*              rotational levels used in channel basis
*    nlevop:   on return contains number of energetically distinct
*              rotational levels used in channel basis which are open
*              asymptotically
*    sc1,sc2:  scratch vectors of length at least nmax
*    sc3,sc4:  scratch vectors of length at least nmax
*              these scratch vectors are not used here
*    rcut:     cut-off point for keeping higher energy channels
*              if any open channel is still closed at r=rcut, then
*              all closed channels as well any open channels which are
*              still closed at r=rcut are dropped from basis
*              note that this is ignored in molecule-surface collisions!!!
*    jtot:     total angular momentum
*    flaghf:   if .true., then system with half-integer spin
*              if .false., then system with integer spin (this is the case
*              here)
*    flagsu:   if .true., then molecule-surface collisons
*    csflag:   if .true., then coupled-states calculation
*              if .false., the close-coupled calculation
*    clist:    if .true., then quantum numbers and energies listed for
*              each channel
*    bastst:   if .true., then execution terminates after the first call
*              to basis
*    ihomo:    if .true. , then homonuclear molecule
*              if .false., then heteronuclear molecule
*              if the molecule is homonuclear (ihomo = .true.), the
*              rotational levels included go from jmin to jmax in steps
*              of 2 and only even lambda terms in the anisotropy are
*              included
*    nu:       coupled-states projection index
*    numin:    minimum coupled states projection index
*              for cc calculations nu and numin are both set = 0 by calling
*              program
*    jlpar:    total parity of included channels in cc calculation
*              only those channels are included for which
*                  (-1)**(j+l-jtot)=jlpar
*    n:        on return equals number of channels
*    nmax:     maximum dimension of arrays
*    ntop:     maximum row dimension of all matrices passed to subroutines
*              propag and soutpt.  ntop is set in basis only if nu = numin
*              otherwise it is unchanged from the value supplied by the
*              calling program
*    note!!!   if flaghf = .true., then the true values of the rotational
*    quantum numbers, the total angular momentum, and the coupled-states
*    projection index are equal to the values stored in j, jtot, and nu
*    plus 1/2
*  variables in common block /cosysr/
*    isrcod:   number of real system dependent variables
*    brot:     rotational constant in cm-1
*  variables in common block /cosysi/
*    nscode:   total number of system dependent variables
*    isicod:   numbr of integer system dependent variables
*    nterm:    number of different associated legendre terms in
*              expansion of potential
*    jmin:     the minimum rotational angular momenta
*    jmax:     the maximum rotational angular momenta
*  variable in common block /cocent/
*    cent:      array containing centrifugal barrier of each channel
*  variable in common block /coeint/
*    eint:      array containing channel energies (in hartree)
*  variables in common block /coered/
*    ered:      collision energy in atomic units (hartrees)
*    rmu:       collision reduced mass in atomic units
*               (mass of electron = 1)
*  variable in common block /conlam/
*    nlam:      the number of angular coupling terms actually used
*    nlammx:    the maximum number of angular coupling terms
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx
*  variable in common block /cov2/
*    nv2max:    maximum core memory allocated for the v2 matrix
*    v2:        lower triangle of nonzero angular coupling matrix elements
*               only nonzero elements are stored
*  variable in common block /coiv2/
*   iv2:        matrix address of v2 matrix for each non-zero element
*  variable in common block /coskip/
*   nskip  for a homonuclear molecule lamda is running in steps of nskip=2
*          for a heteronuclear molecule nskip=1
*
*   iskip   same as nskip, used for consistency check
*
*  subroutines called:
*   vlm1sg:    returns singlet-sigma angular coupling coefficient for
*              particular choice of initial and final channel quantum numbers
*
* --------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical ihomo, flaghf, csflag, clist, flagsu, bastst
      include "common/parbas"
      include "common/parbasl"
      common /covib/ nvib,ivib(maxvib)
      common /coipar/ iiipar(9), iprint
      common /cosysi/ nscode, isicod, nterm, nvmin, nvmax,
     :                iscod(2,maxvib)
      common /cosysr/ isrcod, junkr, rpar(4,maxvib)
      common /cov2/ nv2max, junkv, v2(1)
      common /coiv2/ iv2(1)
      common /conlam/ nlam, nlammx, lamnum(1)
      common /cocent/ cent(1)
      common /coeint/ eint(1)
      common /coered/ ered, rmu
      common /coskip/ nskip, iskip
      common /coconv/ econv, xmconv
      dimension j(5), l(5), jhold(5), ehold(5), sc1(5), sc2(5), sc3(5),
     :          sc4(5), ishold(5), is(5)
*   econv is conversion factor from cm-1 to hartrees
*   xmconv is converson factor from amu to atomic units
      zero = 0.d0
      two = 2.d0
*  check for consistency in the values of flaghf and csflag
      if (flaghf) then
        write (6, 1)
        write (9, 1)
1       format (' *** FLAGHF = .TRUE. FOR SINGLET SYSTEM; ABORT ***' )
        if (bastst) then
          return
        else
          call exit
        end if
      end if
      if (flagsu .and. .not. csflag) then
        write (6, 8)
        write (9, 8)
8      format
     :   ('  *** CSFLAG = .FALSE. FOR SURFACE CALCULATION; ABORT ***')
        if (bastst) then
          return
        else
          call exit
        end if
      end if
* check that jmin .ge. nu for bound state calculations (this shouldn't be a bu
* but is?)
      nsum = 0
        iva=0
        ive=0
*  check that requested vib levels have been defined
        do 10 i=1,nvib
        if(ivib(i).ne.ivib(1)+i-1) then
          write(6,5) (ivib(k),k=1,nvib)
5         format(/' INPUT ERROR: NON-SEQUENTIAL VIBRATIONAL',
     :            ' STATES: ',10i5)
          call exit
        end if
        if(ivib(i).eq.nvmin) iva=i
        if(ivib(i).eq.nvmax) ive=i
10      continue
        if(iva.eq.0.or.ive.eq.0) then
          write(6,15) nvmin,nvmax,(ivib(i),i=1,nvib)
15        format(/' PARAMETERS UNDEFINED FOR VIBRATIONAL STATE'/
     :    1x,' REQUESTED STATES:',i3,'-',i2,' DEFINED STATES:',10i5)
          call exit
        end if

      if (clist) then
        if (flagsu) then
          write (6,16) rmu * xmconv, ered * econv, jtot, jlpar
          write (9,16) rmu * xmconv, ered * econv, jtot, jlpar
16        format(/,' **  1-SIGMA ON UNCORRUGATED SURFACE ** RMU=', f9.4,
     :      '             E=', f7.2, /, ' JTOT=', i5, 2x,' JLPAR=',i2)
        else
          if (.not. csflag) then
            write (6,20) rmu * xmconv, ered * econv, jtot, jlpar
            write (9,20) rmu * xmconv, ered * econv, jtot, jlpar
20          format(/,' **  CC SINGLET SIGMA DIATOMIC ** RMU=', f9.4,
     :           '       E=',f7.2,'   JTOT=', i5, 2x,' JLPAR=',i2)
          else
            write (6,25) rmu * xmconv, ered * econv, jtot, nu
            write (9,25) rmu * xmconv, ered * econv, jtot, nu
25          format(/,' **  CS SINGLET SIGMA DIATOMIC ** RMU=', f9.4,
     :           '       E=',f7.2,'   JTOT=', i5, 2x,' NU=',i2)
          end if
        end if
        if (.not. flagsu) write (9,30) rcut
30      format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=',
     :          f6.2)
*  assign quantum numbers and energies for rotational levels
        write(6,31) ' State    B(v)',' D(v)','H(v)','E(v)'
        write(9,31) ' State    B(v)',' D(v)','H(v)','E(v)'
31      format(/2(a,8x),a,10x,a)
        do 40 iv=iva,ive
        write(6,35) ivib(iv),(rpar(jj,iv),jj=1,4)
40      write(9,35) ivib(iv),(rpar(jj,iv),jj=1,4)
35      format(1x,i3,2x,3g12.5,f15.8)
      end if
      n=0
      nskip = 1
      if (ihomo) nskip = 2
      do 120 iv=iva,ive
      jmin=iscod(1,iv)
      if (boundc.and.csflag) then
          if (jmin.lt.nu) then
             write(6, 7) jmin, nu
             write(9, 7) jmin, nu
7            format
     +    (/' JMIN = ',i2,', .LT. NU = ',i2,' IN BASIS; JMIN RESET')
             jmin=nu
          endif
      endif
      jmax=iscod(2,iv)
      brot=rpar(1,iv)/econv
      drot=rpar(2,iv)/econv
      hrot=rpar(3,iv)/econv
      evib=rpar(4,iv)/econv
      do 120  ji = jmin, jmax, nskip
        jj1=ji*(ji+1)
        ee=brot*jj1 - drot*jj1**2 + hrot*jj1**3 + evib
        if (.not. csflag) then
*  here for cc calculations
        lmax = jtot + ji
        lmin = iabs (jtot - ji)
        do 110  li = lmin, lmax
          ix = (-1) ** (ji + li - jtot)
          if (ix .eq. jlpar) then
*  here for correct combination of j and l
            n = n + 1
            if (n .gt. nmax) go to 130
            l(n) = li
            cent(n) = li * (li + 1.)
            is(n) = ivib(iv)
            j(n) = ji
            eint(n) = ee
          end if
110     continue
        else
*  here for cs calculations
          if (ji .lt. nu) go to 120
          n = n + 1
          if (n .gt. nmax) go to 130
          l(n) = jtot
          if (.not.boundc) then
            cent(n) = jtot * (jtot + 1)
          else
            xjtot=jtot
            xj=j(n)
            xnu=nu
            cent(n)=xjtot*(xjtot+1)+xj*(xj+1)-2*xnu*xnu
          endif
          is(n) = ivib(iv)
          j(n) = ji
          eint(n) = ee
        end if
120   continue
130   if (n .gt. nmax) then
        write (9, 140) n, nmax
        write (6, 140) n, nmax
140     format(/' *** NCHANNELS=', i4,' .GT. MAX DIMENSION OF',
     :         i4,'; ABORT')
        if (bastst) then
          return
        else
          call exit
        end if
      end if
*  now check to see if any of the open channels are closed at r=rcut
*  this is not done for molecule-surface collisions or for rcut < 0
*  and for bound state calculations
      if (.not.flagsu .and. rcut .gt. 0.d0 .and. .not.boundc) then
        emin = 1.e+7
        do 145  i = 1, n
          if (eint(i) .le. ered) then
*  here if channel is
            if ( jtot * (jtot + 1) / (2. * rmu * rcut * rcut)
     :          .gt. (ered - eint(i)) ) then
*  here if channel is open asymptotically but closed at r = rcut
              if (eint(i) .lt. emin) emin = eint(i)
*  emin now contains the lowest channel energy for which this
*  condition is met
            end if
          end if
 145    continue
*  now eliminate all channels with eint .ge. emin if any of the channels
*  are open asymptotically but closed at r = rcut
        if (emin.lt.ered) then
          nn=n
          n = 0
          do 150 i = 1, nn
            if (eint(i) .lt. emin) then
*  here if this channel is to be included
              n = n + 1
              eint(n) = eint(i)
              is(n) = is(i)
              j(n) = j(i)
              cent(n) = cent(i)
              l(n) = l(i)
            end if
150       continue
*  reset number of channels
        end if
      end if
*  return if no channels
      if (n .eq. 0) return
      if (nu .eq. numin) then
        ntop = max(n, nlevop)
*  ntop is the maximum row dimension of all matrices passed in the
*  call list of subroutines propag and soutpt.
*  for fps make sure this is an odd number, for faster bank access.
*  this has no effect on vax or cray
        if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax)
     :       ntop = ntop + 1
      else
        if (n.gt.ntop) then
          write (6, 160) nu, n, ntop
          write (9, 160) nu, n, ntop
160       format (' *** NCH = ',i3, ' AT NU = ',i2, ' .GT. NTOP = ',i3,
     :            '; ABORT **',/,
     :    '     CHECK RCUT')
          call exit
        endif
      end if
      nlevel = 0
      nlevop = 0
*  form list of all energetically open rotational levels included in the
*  calculations and their energies
*  if homonuclear diatomic, skip space is two
      do 200 iv=iva,ive
      jmin=iscod(1,iv)
      jmax=iscod(2,iv)
      brot=rpar(1,iv)/econv
      drot=rpar(2,iv)/econv
      hrot=rpar(3,iv)/econv
      evib=rpar(4,iv)/econv
      do 200  ji = jmin, jmax, nskip
        jj1=ji*(ji+1)
        ee=brot*jj1 - drot*jj1**2 + hrot*jj1**3 + evib
        nlevel = nlevel + 1
        ehold(nlevel) = ee
        jhold(nlevel) = ji
        ishold(nlevel) = ivib(iv)
200   continue
*  now sort this list to put closed levels at end
*  also determine number of levels which are open
      nlevop = 0
      if (nlevel .gt. 1) then
        do 80  i = 1, nlevel - 1
          if (ehold(i) .le. ered) then
             nlevop = nlevop + 1
          else
            do 75 ii = i + 1, nlevel
              if (ehold(ii) .le. ered) then
                nlevop = nlevop + 1
                ikeep = jhold(i)
                jhold(i) = jhold(ii)
                jhold(ii) = ikeep
                ikeep = ishold(i)
                ishold(i) = ishold(ii)
                ishold(ii) = ikeep
                ekeep = ehold(i)
                ehold(i) = ehold(ii)
                ehold(ii) = ekeep
                go to 80
              end if
75          continue
          end if
80      continue
      else
* here for only one level
        if (ehold(1) .le. ered) then
          nlevop=1
        else
          nlevop=0
        endif
      endif
      if (nlevop .le. 0) then
        write (9, 85)
        write (6, 85)
85      format('*** NO OPEN LEVELS IN BA1SG; ABORT')
        if (bastst) return
        call exit
      endif
      if (ehold(nlevel) .le. ered) nlevop = nlevop + 1
      if (nu .eq. numin) then
        ntop = max(n, nlevop)
*  ntop is the maximum row dimension of all matrices passed in the
*  call list of subroutines propag and soutpt.
*  for fps make sure this is an odd number, for faster bank access.
*  this has no effect on vax or cray
        if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax)
     :       ntop = ntop + 1
      end if
*  now list channels if requested
      if (clist) then
        write (6, 255)
        write (9, 255)
255     format(/'   N   V   J    L      EINT(CM-1)',/)
        do 265  i = 1, n
          write (6, 260) i, is(i), j(i), l(i), eint(i) * econv
          write (9, 260) i, is(i), j(i), l(i), eint(i) * econv
260       format (3i4, i5, f13.3)
265     continue
        write (6, 256)
256     format(/' OPEN LEVELS:'//
     1          '   N   V   J      EINT(CM-1)',/)
        do 266  i = 1, nlevop
          write (6, 261) i, ishold(i), jhold(i),  ehold(i) * econv
261       format (3i4, f13.3)
266     continue
      end if
*  now calculate coupling matrix elements
      if (bastst .and. iprint.gt.1) then
        write (6, 280)
        write (9, 280)
280     format (/'  VR VC  LAMBDA ILAM    I   ICOL  IROW',
     :           '  IV2       VEE')
      end if
*
*  the following only for safety. Checks that pot routine
*  contains all required vib. levels in correct order
*
      do 295 irow=1,n
      ivr=is(irow)
      do 290 icol=1,irow
      ivc=is(icol)
      do 290 ivb=1,ntv(1)
      if(ivrow(ivb,1).eq.ivr.and.
     :   ivcol(ivb,1).eq.ivc) goto 295
290   continue
      write(6,291) ivr,ivc
      write(9,291) ivr,ivc
291   format(/' VIBRATIONAL STATE NOT DEFINED IN POT',2i5)
      call exit
295   continue
* i counts v2 elements
* inum counts v2 elements for given lambda
* ilam counts numver of v2 matrices
* ij is address of given v2 element in present v2 matrix
      i = 0
      ilam=0
      do 320 iv=1,ntv(1)
          ivr=ivrow(iv,1)
          ivc=ivcol(iv,1)
      do 320 il = lammin(1), lammax(1), nskip
        lb=il
        ilam=ilam+1
        inum = 0
        ij=0
        do 310  icol= 1, n
          do 300  irow = icol, n
            if(is(irow).ne.ivr.or.is(icol).ne.ivc) goto 300
            ij = ntop * (icol - 1) +irow
*  here for coupling between molecular rotational levels
            call vlm1sg (j(irow), l(irow), j(icol), l(icol), jtot,
     :                   nu, lb, vee, csflag)
            if (vee .eq. 0) goto 300
              i = i + 1
              inum = inum + 1
               if (bastst .and. iprint.gt.1) then
                 write (6, 340) ivr,ivc,il,ilam,i,irow,icol,
     :                          ij,vee
                 write (9, 340) ivr,ivc,il,ilam,i,irow,icol,
     :                          ij,vee
340              format (1x,2i3,6i6, g17.8)
               end if
              if (i .le. nv2max) then
                v2(i) = vee
                iv2(i) = ij
              end if
300       continue
310     continue
        if(ilam.gt.nlammx) then
          write(6,311) ilam, nlammx
311       format(/' ILAM =',i3,' .GT. NLAMMX =',i3,' IN BA1SG; ABORT')
          call exit
        end if
        lamnum(ilam) = inum
        if (bastst .and. iprint.ge.1) then
          write (6, 420) ivr,ivc,il,inum
          write (9, 420) ivr,ivc,il,inum
420       format(' IVR=',i2,'  IVC=',i2,'  LAMBDA=',i2,
     :           ' NUMBER OF NONZERO MATRIX ELEMENTS',i5)
        end if
        if(inum.ne.0) nlam=ilam
320   continue
      if ( i.gt. nv2max) then
        write (6, 350) i, nv2max
        write (9, 350) i, nv2max
350     format (' *** NUMBER OF NONZERO V2 ELEMENTS = ',i6,
     :           ' .GT. NV2MAX=',i6,'; ABORT ***')
        if (bastst) then
          return
        else
          call exit
        end if
      end if
      if (clist) then
        write (6, 360) i
        write (9, 360) i
360     format (' ** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS',
     :           i6)
      end if
      return
      end
* --------------------------------------------------------------------
      subroutine vlm1sg (j1, l1, j2, l2, jtot, nu, lb, v, csflag)
* --------------------------------------------------------------------
*  subroutine to evaluate the angular coupling matrix element for rotationally
*  inelastic collisions of a 1sigma molecule with a structureless target
*  current revision date: 5-sep-88
*  variables in call list:
*  j1,l1:    initial rotational and orbital angular momenta
*  j2,l2:    final rotational and orbital angular momenta
*  jtot:     total angular momentum
*  nu:       coupled states projection index
*  lb:       value of legendre expansion index
*  v:        on return:  contains desired coupling matrix element
*  csflag:   if .true., then coupled states calculation
*  subroutines called:
*  xf3j:     evaluates 3j symbol
*  xf6j:     evaluates 6j symbol
* --------------------------------------------------------------------
      implicit double precision (a-h,o-z)
*      real v, xj1, xj2, xjtot, xl1, xl2, xlb, xnu, zero, xf3j, xf6j
      integer ij, il, j1, j2, jtot, l1, l2, lb, nu
      logical csflag
      zero = 0.
      v = 0.
      xj1 = j1
      xl1 = l1
      xj2 = j2
      xl2 = l2
      xjtot = jtot
      xlb = lb
      xnu = nu
      if (.not. csflag) then
*  here for close-coupling, for the matrix element see w.lester, meth. comp.
*  phys. 10, 211 (1971)
        ij = (-1) ** (j1 + j2 + lb)
        il = (-1) ** (l1 + l2 + lb)
        if (ij .ne. -1  .and. il .ne. -1) then
          v = xf3j (xj1, xlb, xj2, zero, zero, zero) *
     :        xf3j (xl1, xlb, xl2, zero, zero, zero) *
     :        xf6j (xj1, xj2, xlb, xl2, xl1, xjtot) *
     :        sqrt ( (2 * j1 + 1.) * (2 * j2 + 1.) * (2 * l1 + 1.)
     :             * (2 * l2 + 1.) ) * (-1) ** ( j1 + j2 - jtot )
        end if
      else
*  here for coupled-states, see p. mcguire and d.j. kouri, j. chem. phys.
*  60, 2488 (1974)
        v = xf3j (xj1, xlb, xj2, zero, zero, zero) *
     :      xf3j (xj1, xlb, xj2, -xnu, zero, xnu) *
     :      sqrt ( (2 * j1 + 1.) * (2 * j2 + 1.) ) * (-1) ** nu
      end if
      return
      end
c     ------------------------------------------------------------------
      logical function is_twomol(ibasty)
c     ------------------------------------------------------------------
c
c     checks if a basis is for molecule-molecule collision (j=10j1+j2)
c
c     written by q. ma
c     current revision:  24-jul-2019 (p.dagdigian)
c     ------------------------------------------------------------------
      implicit none
      integer, intent(in) :: ibasty
      if ((ibasty .eq. 9) .or. (ibasty .eq. 20) .or. (ibasty .eq. 21)
     $      .or. (ibasty .eq. 25) .or. (ibasty .eq. 26)
     $      .or. (ibasty .eq. 28) .or. (ibasty .eq. 30) 
     $      .or. (ibasty .eq. 100))
     $     then
         is_twomol = .true.
      else
         is_twomol = .false.
      end if
      return
      end function is_twomol
c     ------------------------------------------------------------------
      logical function is_j12(ibasty)
c     ------------------------------------------------------------------
c
c     checks if j12 is used in a basis
c
c     written by q. ma
c     current revision:  17-oct-2018 (p.dagdigian)
c     ------------------------------------------------------------------
      implicit none
      integer, intent(in) :: ibasty
      logical :: is_twomol
      if (is_twomol(ibasty) .or. (ibasty .eq. 12)
     $     .or. (ibasty .eq. 13) .or. (ibasty .eq. 15)
     $     .or. (ibasty .eq. 23)) then
         is_j12 = .true.
      else
         is_j12 = .false.
      end if
      return
      end function is_j12
