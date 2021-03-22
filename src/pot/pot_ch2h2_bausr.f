c     pot_ch2h2.f
c
c     based on pot_ohh2_bausr.f, written by Q. Ma 
c     author: Paul Dagdigian
c     current revision:  18-jan-2021
c
c     Pot routine for a collision between CH2 and H2.
c
c     PES to be read from data file specified in the input file.
c
c     This basis routine should have ibasty = 100, so hibridon will
c     treat it as a mol-mol problem with j12.
c
c     DUMMY SUBROUTINES FOR ALL USER-DEFINED BASIS/POT.
      include "common/ground"
c
c     ------------------------------------------------------------------
c     This module contains (explictly) the number of terms and their
c     indices in the expansion of the PES.  Its contents should be
c     filled in the pot routine.
c
c     This module replaces lammin, lammax, mproj in hibridon.
c
      module mod_lams
      implicit none
c
      type lm_type
      integer :: l1, m1, l2, ltot
      end type lm_type
c
      type(lm_type), dimension(:), allocatable :: lms
      end module mod_lams
c     ----------------------------------------------------------------- 
C     Module containing the shared arrays for this pot routine
c
      module mod_ch2h2
      implicit none
      integer :: nr1
      real(8), dimension(:), allocatable :: rr1
      real(8), dimension(:, :), allocatable :: coef, spl_b, 
     $  spl_c, spl_d
      real(8) econv, xmconv, sq4pi
      parameter (econv=219474.6315343234d0,
     $  xmconv=0.0005485799094979479d0, sq4pi=3.544907701811032d0)
      end module mod_ch2h2
*  -----------------------------------------------------------------------
c     THE FOLLOWING SOUBROUTINE WILL BE THE MAIN FUNCTION FOR MAKEPOT.
c     NOTE THAT VVL IS IN HARTREES.
      subroutine driver
      use mod_lams
      use mod_ch2h2
      implicit none
c
*      include 'pot_ch2h2_bausr_common.f'
      common /conlam/ nlam
      common /covvl/ vvl
      character*40 filenm
      real(8), dimension(1) :: vvl
      double precision r, vv0
      integer i, nlam
c
      filenm = 'pot_ch2h2.dat'
*      print *, 'Input the name of potdata file:'
*      read (5, *, end=99) filenm
      call loapot(10, filenm)
 10   print *, 'R (bohr), Ctrl+D to exit:'
      read (5, *, end=99) r
      call pot(vv0, r)
      write (6, *) 'angular coefficients'
      write (6, 20) (lms(i)%l1, lms(i)%m1, lms(i)%l2, lms(i)%ltot,
     $     vvl(i) * econv, i=1, nlam)
 20   format (3(4(i2, 1x), 1x, 1pe16.8, 2x))
      goto 10
 99   end
c     ------------------------------------------------------------------
c     DATA FILE, IF REQUIRED, CAN BE LOADED WITH THIS SOUBROUTINE.
      subroutine loapot(iunit, filnam)
      use mod_lams
      use mod_ch2h2
      implicit none
c
*      include 'pot_ch2h2_bausr_common.f'
      include 'common/parpot'
c
      character*(*) filnam
      integer :: iunit, ir, iv, nv, ifrsts, junks, nlam, nlammx, lamnum
      integer MAX_NR, MAX_NVB
      real(8), dimension(1) :: vvl
      common /covvl/ vvl
      common /conlam/ nlam
      character*255 datfl
      parameter (MAX_NR=400, MAX_NVB=800)
C     
c     WHEN HIBRIDON LOADS, A STRING CONTAINING ONLY ONE SPACE WILL BE
c     PASSED TO THIS SUBROUTINE.  MAKE SURE NOTHING WILL BE DONE.
      if (filnam .eq. ' ') return
C
      potnam = 'Dagdigian AVTZ-F12a CH2-H2 PES'
C
c     LOAD DATA TO A COMMON BLOCK ONLY USED IN THIS FILE, AND DO PRE-
c     PROCESSING IF NECESSARY.
C     
C     Read  angular coefficients
      call datfln(trim(filnam), datfl)
      open (unit=iunit, file=datfl, status='old')
      read (iunit, *) nr1
      if (nr1 .gt. MAX_NR) call raise('too many R in data file.')
      if (allocated(rr1)) deallocate(rr1)
      allocate(rr1(nr1))
      read (iunit, *) rr1
      read (iunit, *) nv
      nlam = nv
      if (nv .gt. MAX_NVB)
     $     call raise('too many vlm terms in data file.')
      nlam = nv
      if (allocated(lms)) deallocate(lms)
      allocate(lms(nv))
      if (allocated(spl_b)) deallocate(spl_b)
      allocate(spl_b(nr1, nv))
      if (allocated(spl_c)) deallocate(spl_c)
      allocate(spl_c(nr1, nv))
      if (allocated(spl_d)) deallocate(spl_d)
      allocate(spl_d(nr1, nv))
      if (allocated(coef)) deallocate(coef)
      allocate(coef(nr1, nv))
C     
      do iv = 1, nv
         read (iunit, *) lms(iv)%l1, lms(iv)%m1, lms(iv)%l2, 
     :     lms(iv)%ltot
         read (iunit, *) (coef(ir, iv), ir = 1, nr1)
      end do
c  convert to a.u.
      coef = coef / econv
C     
      close(unit=iunit)
C     Spline parameters prepared here
      do iv = 1, nv
         call spline(nr1, rr1, coef(1, iv), spl_b(1, iv), 
     $        spl_c(1, iv), spl_d(1, iv))
      end do
      return
      end subroutine loapot
c     ------------------------------------------------------------------
c     A POT ROUTINE RETURNS VVL ARRAY (COVVL BLOCK) FOR A GIVEN R.
c     ALWAYS SET VV0 = 0 TO AVOID CONFUSION.  THE ISOTROPIC TERM CAN BE
c     EASILY TREAT AS ONE TERM IN THE POTENTIAL EXPANSION.  NOTE THAT
c     VVL SHOULD BE IN HARTREES.
      subroutine pot(vv0, r_raw)
      use mod_ch2h2
      implicit none
c
*      include 'pot_ch2h2_bausr_common.f'
      real(8), intent(in) :: r_raw
      real(8), intent(out) :: vv0
      real(8) seval, r, vvl
      integer :: iv, nv, nlam
      common /covvl/ vvl(1)
      common /conlam/ nlam
c
      if (r_raw .lt. 3.5d0) then
         r = 3.5d0
      else
         r = r_raw
      end if
      vv0 = 0d0
      do iv = 1, nlam
        vvl(iv) = seval(nr1, r, rr1, coef(1, iv),
     $    spl_b(1, iv), spl_c(1, iv), spl_d(1, iv))
      end do
      return
      end
c     ------------------------------------------------------------------
c     THIS SUBROUTINE GOVERNS THE INPUT/OUTPUT OF THE BASIS ROUTINE.
c     ONLY IREAD IS USED: RETURN DIRECTLY IF ZERO.
      subroutine syusr(irpot, readpt, iread)
*  subroutine to read in system dependent parameters for CH2(X 3B1) (0,v2,0)
*      bender + atom scattering
*  if iread = 1 read data from input file
*  if iread = 0 return after defining variable names
*  NOTE THAT THIS VERSION DOES NOT USE DATA FROM VFIT (FOLLMEG FORM)
*  modified from syatp subroutine
*  
*  based on sysch2x, writted by p.dagdigian
*  current revision date: 16-jan-2021 by p.dagdigian
*  -----------------------------------------------------------------------
*  variables in common bloc /cosysr/
*    isrcod:   total number of real system dependent variables
*    emax:     the maximum rotational energy (in cm-1) for a channel to be
*              included in the basis
*  variables in common block /cosysi/
*    nscode:   total number of variable names which are passed to HINPUT
*              nscod must equal isrcod + isicod + 3
*    isicod:   total number of integer system dependent variables
*    nterm:    number of different associated legendre terms in
*              expansion of potential
*    numpot:   the number of the potential used, this variable is passed
*              to the pot subroutine
*    ipotsy:   cylindrical symmetry of potential.  Use this parm to
*              distinguish between AB2 and ABC type molecules.  for the
*              former, only even lambda allowed (set delta_lambda = ipotsy)
*    iop:      ortho/para label for molecular states. If ihomo=.true. then onl
*              para states will be included if iop=1 and only ortho states if
*              iop=-1
*    ivbend:   bend vibrational quantum number (can equal 0 to 3)
*    jmax:     maximum rotational quantum number for the molecule
*    j2min:    minimum rotational quantum number for H2 collision partner
*    j2max:    maximum rotational quantum number for H2 collision partner
*    ipotsy2:  symmetry of the H2 collision partner (set to 2)
*  variable in common /cosys/
*    scod:    character*8 array of dimension nscode, which contains names
*             of all system dependent parameters.  Note that the ordering
*             of the variable names in scod must correspond to the ordering
*             of the variable names in cosysi followed by the ordering of
*             variable names in cosysr followed by LAMMIN, LAMMAX, and MPROJ
*  -----------------------------------------------------------------------
      logical readpt, existf
      double precision emax, b2rot
      integer numpot
      integer icod, ircod, lencod
      integer i, iop, iread, irpot, isicod, isrcod, ipotsy, jmax,
     :        nscode, nterm, ivbend
      character*8 scod
      character*1 dot
      character*(*) fname
      character*60 line, filnam, potfil, filnm1
      parameter (icod=9, ircod=2)
      parameter (lencod = icod + ircod + 3)
      include "common/parbas"
      common /coiout/ niout, indout(20)
      common /cosys/ scod(lencod)
      common /cosysi/ nscode, isicod, nterm, numpot, ipotsy, iop,
     :  ivbend, jmax, j2min, j2max, ipotsy2
      common /cosysr/ isrcod, junkr, emax, b2rot
      common /conlam/ nlam
      save potfil
*  number and names of system dependent parameters
*  first all the system dependent integer variables
*  in the same order as in the common block /cosysi/
*  variable names must be in uppercase of maximum length 6 characters
*  NTERM must be the first variable
*  followed by the system dependent real variables
*  in the same order as in the common block /cosysr/
*  then the three variable names LAMMIN, LAMMAX, MPROJ, in that order
      include "common/comdot"
      scod(1)='NTERM'
      scod(2)='NUMPOT'
      scod(3)='IPOTSY'
      scod(4)='IOP'
      scod(5)='IVBEND'
      scod(6)='JMAX'
      scod(7)='J2MIN'
      scod(8)='J2MAX'
      scod(9)='IPOTSY2'
      scod(10)='EMAX'
      scod(11)='B2ROT'
      scod(12)='LAMMIN'
      scod(13)='LAMMAX'
      scod(14)='MPROJ'
      nscode = lencod
      isicod = icod
      isrcod = ircod
      irpot = 1
*  set default values for asymmetric top scattering
*  (symmetric CH2 molecule)
       nterm = 7
      if (iread .eq. 0) then
        mproj(1) = 0
        mproj(2) = 1
        mproj(3) = 2
        mproj(4) = 3
        mproj(5) = 4
        mproj(6) = 5
        mproj(7) = 6

        lammin(1) = 2
        lammin(2) = 1
        lammin(3) = 2
        lammin(4) = 3
        lammin(5) = 4
        lammin(6) = 5
        lammin(7) = 6
        lammax(1) = 6
        lammax(2) = 5
        lammax(3) = 6
        lammax(4) = 7
        lammax(5) = 8
        lammax(6) = 7
        lammax(7) = 8


        ipotsy = 2
        iop = 1
        jmax = 3
        niout=5
*
*  INDOUT VALUES TO BE SET ***
        niout=11
        indout(1)=0
        indout(2)=1
        indout(3)=-1
        indout(4)=2
        indout(5)=-2
        indout(6)=3
        indout(7)=-3
        indout(8)=4
        indout(9)=-4
        indout(10)=5
        indout(11)=-5
      endif
      potfil=' '
      if (iread .eq. 0) return
*  line 18
      read (8, *, err=80) ipotsy, iop, ivbend
*  line 19
      read (8, *, err=80) jmax, emax
*  line 20
      read (8, *, err=80) j2min, j2max, ipotsy2, b2rot
      read (8, *, err=80) potfil
      call loapot(10, potfil)
      close(8)
      return
* here if read error occurs
80    write(6,90)
90    format(/'   *** ERROR DURING READ FROM INPUT FILE ***')
      return
*  -----------------------------------------------------------------------
      entry ptrusr (fname, readpt)
      line = fname
      readpt = .true.
100   if (readpt) then
        l=1
        call parse(line,l,filnam,lc)
        if(lc.eq.0) then
          write(6,102)
102       format(' FILENAME MISSING FOR POTENTIAL INPUT')
        end if
        j=index(filnam(1:lc),dot)
        if(j.eq.0) then
           call gennam(potfil,filnam,0,'BIN',lc)
           filnam = potfil
        end if
        potfil=filnam
        filnm1 = 'potdata/'//filnam
        inquire(file=filnm1,exist=existf)
        if(.not.existf) then
          write(6,105) filnam(1:lc)
105      format(' FILE ',(a),' NOT FOUND')
         return
        end if
* now call loapot(iunit,filnam) routine to read potential parameters
        call loapot(1,filnam)
      end if
      close (8)
      return
c     ------------------------------------------------------------------
      entry savusr(readpt)
c     WRITE THE LAST FEW LINES OF THE INPUT FILE.
*  save input parameters for CH2(X 3B1) (0,v2,0) bender + H2 scattering
*  the order of the write statements should be identical to the read statement
*  above. for consistency with the data file written by gendat, format
*  statements should reserve the first 30 spaces for data, spaces 31-33 should
*  be left blank, and the names of the variables should be printed in spaces
*  34-80
*  line 18:
      write (8, 220) ipotsy, iop, ivbend
220   format (3i4, 18x,'   ipotsy, iop, ivbend')
*  line 20
      write (8, 230) jmax, emax
230   format (i4, g12.5,14x, '   jmax, emax')
      write (8, 240) j2min, j2max, ipotsy2, b2rot
240   format(3i4,14x, '   j2min, j2max, ipotsy2, b2rot')
      write (8, *) potfil
      return
      end
* --------------------------------------------------------------------
c     THE BASIS ROUTINE.  JHOLD, EHOLD, ISHOLD, NLEVEL, NLEVOP ARE TO
c     BE RETURNED FOR LEVEL INFORMATION; J, L, IS, EINT, CENT, N,
c     NTOP, AND POSSIBLY J12 ARE TO BE RETURNED FOR CHANNEL INFORMATION;
c     ALSO CALCULATE NON-ZERO COUPLING MATRIX ELEMENTS AND RETURN IN V2,
c     IV2, LAMNUM.
c
* ----------------------------------------------------------------------
      subroutine bausr (j, l, is, jhold, ehold, ishold, nlevel, nlevop,
     :                  etemp, fjtemp, fktemp, fistmp,
     :                  rcut, jtot, flaghf, flagsu,
     :                  csflag, clist, bastst, ihomo, nu, numin, jlpar,
     :                  n, nmax, ntop)
* --------------------------------------------------------------------
*  subroutine to determine angular coupling potential for collision
*  of CH2(X 3B1) in a (0,v2,0) bender vibrational level with a structureless
*  atom or with an uncorrugated surface
*  this electronic state can be described as a non-rigid bender, and
*  its rotational energies are poorly described by the standard expressions.
*  the rotational energies of this electronic state are read from a table of
*  energies from a MORBID calculation by P. Jensen [P. R. Bunker and P. Jensen,
*  JCP 85, 3724 (1986) and P. Jensen, private communication]
*  because the effective A rotational constant is very large, especially for
*  the (020) and (030) vibrational levels, the rotational wave functions
*  are well approximated by symmetric top wave functions
*  the scattering calculations are carried out with a spin-free basis
*  since CH2(X 3B1) is well described by Hund's case (b) coupling
*
*  the MORBID energies go up to j = 14.  the energies for higher j (<= 25),
*  for IS = 0,-1,1,...,5, were estimated by extrapolation, fitting the MORBID
*  energies to the rigid rotor expression x0 + x1*j*(j+1) for each IS stack.
*  the extrapolated energies are appended to the MORBID energy list below.
*
*  prmatp, which is in hibaastp.f, is called
*
*  author:  paul dagdigian
*  based on hibach2x.f
*  current revision: 18-jan-2021 (p.dagdigian)
* --------------------------------------------------------------------
*  variables in call list:
*    j:        on return contains rotational quantum number for each
*              channel
*    l:        on return contains orbital angular momentum for each
*              channel
*    is:       on return contains symmetry index (ieps * kp) for each
*              channel
*  Note that we have adopted the following convention for the symmetry
*  index "is" so that on return is = 0/1, with the symmetric top basis
*  functions written as [|j,kp> + (-1)^is*|j,-kp)\, respectively.
*    jhold:    on return contains rotational quantum numbers for each
*              rotational level
*    ehold:    on return contains energy in hartrees of each rotational
*              level
*    ishold:   on return contains symmetry index of each energetically
*              distinct level
*    nlevel:   on return contains number of energetically distinct
*              rotational levels used in channel basis
*    nlevop:   on return contains number of energetically distinct
*              rotational levels used in channel basis which are open
*              asymptotically
*    etemp:    scratch array used to create channel list
*    fjtemp:   scratch array used to create channel list
*    fktemp:   scratch array used to create channel list
*    fistmp:   scratch array used to create channel list
*    rcut:     cut-off point for keeping higher energy channels
*              if any open channel is still closed at r=rcut, then
*              all closed channels as well any open channels which are
*              still closed at r=rcut are dropped from basis
*              note that this is ignored in molecule-surface collisions!!!
*    jtot:     total angular momentum
*              in cc calculation jtot is the total angular momentum
*              in cs calculation jtot is the l-bar quantum number
*    flaghf:   if .true., then system has half-integer spin
*              if .false., then system has integer spin
*    flagsu:   if .true., then molecule-surface collisons
*    csflag:   if .true., then coupled-states calculation
*              if .false., the close-coupled calculation
*              CSFLAG IS IGNORED IN THIS BASIS ROUTINE.  ONLY CLOSE
*              COUPLED CALCULATIONS CARRIED OUT
*    clist:    if .true., then quantum numbers and energies listed for
*              each channel
*    bastst:   if .true., then execution terminates after the first call
*              to basis
*    ihomo:    if .true., then the molecule posseses interchange symmetry
*              (e.g. H2O, CH2, H2CO), so that only the ortho or para levels will be
*              included depending on the value of the parameter iop in common
*              /cosysi/ (see below)
*              subroutine currently to allow only BA2 or A2BC type molecules,
*              where the a inertial axis is an axis of C2 symmetry.
*    nu:       coupled-states projection index
*    numin:    minimum coupled states projection index
*              for cc calculations nu and numin are both set = 0 by calling
*              program
*    jlpar:    total parity of included channels in cc calculation
*              only those channels are included for which
*                  -is*(-1)**(j+kp+l-jtot)=jlpar
*              where parity designates the parity of the molecular state
*              in cs calculation jlpar is set equal to 1 in calling program.
*              note that the initial minus sign comes from the electronic
*              symmetry of the CH2(X 3B1) state.
*    n:        on return equals number of channels
*    nmax:     maximum dimension of arrays
*    ntop:     maximum row dimension of all matrices passed to subroutines
*              propag and soutpt.  ntop is set in basis only if nu = numin
*              otherwise it is unchanged from the value supplied by the
*              calling program
*  variables in common bloc /cosysr/
*    isrcod:   total number of real system dependent variables
*    emax:     the maximum rotational energy (in cm^-1) for a channel to be
*              included in the basis
*    b2rot:    rotational constant of collision partner
*  variables in common block /cosysi/
*    nscode:   total number of variable names which are passed to HINPUT
*              nscod must equal isrcod + isicod + 3
*    isicod:   total number of integer system dependent variables
*    nterm:    number of different associated legendre terms in
*              expansion of potential
*    numpot:   the number of the potential used, this variable is passed
*              to the pot subroutine
*    ipotsy:   cylindrical symmetry of potential.  Use this parm to
*              distinguish between BA2 and ABC type triatomics
*              set ipotsy = 2 and 1, respectively, for these types of
*              molecules
*    iop:      ortho/para label for molecular states. If ihomo=.true. then only
*              para states will be included if iop=1 and only ortho states if
*              iop=-1
*    ivbend:   bend vibrational quantum number (can equal 0 to 3)
*    jmax:     the maximum rotational angular momentum for the asymmetric top
*    j2min:    minimum rotational quantum number for H2 collision partner
*    j2max:    maximum rotational quantum number for H2 collision partner
*    ipotsy2:  symmetry of the H2 collision partner (set to 2)
*  variable in common block /cocent/
*    cent:      array containing centrifugal barrier of each channel
*  variable in common block /coeint/
*    eint:      array containing channel energies (in hartree)
*               the zero of energy is assumed to be the 0(0,0) level
*  variable in common /coj12/ 
*    j12:       array containing the j12 quantum number of each channel  
*  variables in common block /coered/
*    ered:      collision energy in atomic units (hartrees)
*    rmu:       collision reduced mass in atomic units
*               (mass of electron = 1)
*  variables in common block /conlam/
*    nlam:      the number of angular coupling terms actually used
*    nlammx:    the maximum number of angular coupling terms
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx
*  variables in common block /cov2/
*    nv2max:    maximum core memory allocated for the v2 matrix
*    v2:        lower triangle of nonzero angular coupling matrix elements
*               stored in packed column form that is (1,1), (2,1), (3,1) ...
*               (n,1), (2,2), (3,2) ... (n,2), etc.
*  variable in common block /coiv2/
*   iv2:        row+column index of v2 matrix for each non-zero element
*  variables in common block /coconv/
*   econv:      conversion factor from cm-1 to hartrees
*   xmconv:     converson factor from amu to atomic units
*  subroutines called:
* --------------------------------------------------------------------
      use mod_lams
      implicit double precision (a-h,o-z)
      logical flaghf, csflag, clist, flagsu, ihomo, bastst
*      include 'pot_ch2h2_bausr_common.f'
      character*1 slab
      include "common/parbas"
      include "common/parbasl"
      common /cosysi/ nscode, isicod, nterm, numpot, ipotsy, iop,
     :  ivbend, jmax, j2min, j2max, ipotsy2
      common /coipar/ iiipar(9), iprint
      common /cosysr/ isrcod, junkr, emax, b2rot
      common /cov2/ nv2max, junkv, v2(1)
      common /coiv2/ iv2(1)
      common /conlam/ nlam, nlammx, lamnum(1)
      common /cocent/ cent(1)
      common /coeint/ eint(1)
      common /coj12/  j12(1)
      common /coered/ ered, rmu
      common /coconv/ econv, xmconv
      dimension j(1), l(1), is(1), jhold(1), ehold(1),
     :          ishold(1), etemp(1), fjtemp(1), fktemp(1),
     :          fistmp(1)
      dimension j2rot(100), e2rot(100) 
*  table of rotational energies from Jensen's MORBID calculation
      dimension ch2x_e(6,282)
      data num_x /282/
      data ch2x_e /
*         N    IS     (000)     (010)     (020)     (030)
     :   0, 0, 0.d0, 0.d0, 0.d0, 0.d0,
*
     :   1, 0, 15.63553d0, 15.45891d0, 15.43323d0, 15.64797d0,
     :   1,-1, 78.32436d0, 169.09658d0, 366.82807d0, 530.70204d0,
     :   1, 1, 79.51666d0, 170.36361d0, 368.29006d0, 523.48901d0,
*
     :   2, 0, 46.88059d0, 46.36126d0, 46.29544d0, 46.95978d0,
     :   2,-1, 108.47063d0, 199.02816d0, 396.62786d0, 560.42062d0,
     :   2, 1, 112.04307d0, 202.82491d0, 401.00989d0, 565.77708d0,
     :   2,-2, 276.2217d0, 467.70371d0, 794.44085d0, 1068.0774d0,
     :   2, 2, 276.23868d0, 467.70641d0, 794.41897d0, 1068.01296d0,
*
     :   3, 0, 93.68375d0, 92.67639d0, 92.57833d0, 93.96709d0,
     :   3,-1, 153.66321d0, 243.90508d0, 441.31633d0, 605.00397d0,
     :   3, 1, 160.79454d0, 251.48552d0, 450.06853d0, 615.70287d0,
     :   3,-2, 323.40787d0, 514.84382d0, 841.61174d0, 1115.61821d0,
     :   3, 2, 323.49228d0, 514.85703d0, 841.61174d0, 1115.29725d0,
     :   3,-3, 566.95251d0, 848.80601d0, 1272.07249d0, 1628.01346d0,
     :   3, 3, 566.95251d0, 848.8062d0, 1272.07338d0, 1628.01346d0,
*
     :   4, 0, 155.96925d0, 154.35917d0, 154.2698d0, 156.71712d0,
     :   4,-1, 213.87056d0, 303.70319d0, 500.88041d0, 664.45975d0,
     :   4, 1, 225.72555d0, 316.30813d0, 515.4407d0, 682.25776d0,
     :   4,-2, 386.2685d0, 577.651d0, 904.72857d0, 1178.99664d0,
     :   4, 2, 386.51967d0, 577.6894d0, 904.40064d0, 1178.04029d0,
     :   4,-3, 630.16553d0, 912.02902d0, 1335.4209d0, 1691.45712d0,
     :   4, 3, 630.16695d0, 912.03037d0, 1335.42717d0, 1691.47595d0,
     :   4,-4, 933.10362d0, 1293.31254d0, 1793.22753d0, 2208.02691d0,
     :   4, 4, 933.10362d0, 1293.31253d0, 1793.2275d0, 2208.02638d0,
*
     :   5, 0, 233.63883d0, 231.3509d0, 231.35444d0, 235.27185d0,
     :   5,-1, 289.05169d0, 378.39092d0, 575.30329d0, 738.79987d0,
     :   5, 1, 306.77579d0, 397.24348d0, 597.09287d0, 765.4286d0,
     :   5,-2, 464.7578d0, 656.08615d0, 983.43531d0, 1258.20461d0,
     :   5, 2, 465.33765d0, 656.17216d0, 982.67145d0, 1255.99832d0,
     :   5,-3, 709.12703d0, 990.9801d0, 1414.48228d0, 1770.57379d0,
     :   5, 3, 709.13263d0, 990.98545d0, 1414.50789d0, 1770.64856d0,
     :   5,-4, 1012.43898d0, 1372.72078d0, 1872.74664d0, 2284.37526d0,
     :   5, 4, 1012.43899d0, 1372.72076d0, 1872.74638d0, 2284.3724d0,
     :   5,-5, 1361.9649d0, 1789.89377d0, 2352.74712d0, 2807.83902d0,
     :   5, 5, 1361.96492d0, 1789.89377d0, 2352.74712d0, 2807.83903d0,
*
     :   6, 0, 326.5744d0, 323.58061d0, 323.81411d0, 329.70628d0,
     :   6,-1, 379.15694d0, 467.92977d0, 664.5649d0, 828.04206d0,
     :   6, 1, 403.87055d0, 494.22983d0, 694.98328d0, 865.19618d0,
     :   6,-2, 558.81944d0, 750.10111d0, 1077.807d0, 1353.23188d0,
     :   6, 2, 559.96366d0, 750.26462d0, 1076.28417d0, 1348.89565d0,
     :   6,-3, 803.799d0, 1085.60662d0, 1509.16645d0, 1865.22073d0,
     :   6, 3, 803.81555d0, 1085.62257d0, 1509.24594d0, 1865.44285d0,
     :   6,-4, 1107.54901d0, 1467.89166d0, 1968.00782d0, 2578.48378d0,
     :   6, 4, 1107.54909d0, 1467.89152d0, 1968.00654d0, 2578.4808d0,
     :   6,-5, 1457.52339d0, 1885.56692d0, 2449.33205d0, 2903.35367d0,
     :   6, 5, 1457.52437d0, 1885.56693d0, 2449.33204d0, 2903.35373d0,
     :   6,-6, 1844.25956d0, 2330.97403d0, 2930.25208d0, 3513.85456d0,
     :   6, 6, 1844.25956d0, 2330.97403d0, 2930.25208d0, 3513.85456d0,
*
     :   7, 0, 434.64152d0, 430.96652d0, 431.62851d0, 440.10517d0,
     :   7,-1, 484.12896d0, 572.27492d0, 768.64259d0, 932.21172d0,
     :   7, 1, 516.92109d0, 607.19927d0, 809.06212d0, 981.53327d0,
     :   7,-2, 668.38708d0, 859.63909d0, 1187.80302d0, 1464.06581d0,
     :   7, 2, 670.41296d0, 859.91582d0, 1185.07641d0, 1456.45324d0,
     :   7,-3, 914.13513d0, 1195.84461d0, 1619.35019d0, 1975.21798d0,
     :   7, 3, 914.17574d0, 1195.88419d0, 1619.56233d0, 1975.76609d0,
     :   7,-4, 1218.38366d0, 1578.76105d0, 2078.92492d0, 2509.27024d0,
     :   7, 4, 1218.38394d0, 1578.76056d0, 2078.92026d0, 2509.24075d0,
     :   7,-5, 1568.87085d0, 1997.01464d0, 2561.46578d0, 3014.50462d0,
     :   7, 5, 1568.86981d0, 1997.01465d0, 2561.4658d0, 3014.50499d0,
     :   7,-6, 1956.12001d0, 2442.97072d0, 3041.52118d0, 3626.07942d0,
     :   7, 6, 1956.12001d0, 2442.97072d0, 3041.52118d0, 3626.07941d0,
     :   7,-7, 2372.97755d0, 2911.02087d0, 3588.65795d0, 4171.26331d0,
     :   7, 7, 2372.97755d0, 2911.02087d0, 3588.65795d0, 4171.26331d0,
*
     :   8, 0, 557.69379d0, 553.41775d0, 554.77588d0, 566.55792d0,
     :   8,-1, 603.90369d0, 691.37582d0, 887.51184d0, 1051.34409d0,
     :   8, 1, 645.82526d0, 736.06494d0, 939.27162d0, 1114.40176d0,
     :   8,-2, 793.38508d0, 984.6352d0, 1313.37676d0, 1590.69094d0,
     :   8, 2, 796.69505d0, 985.06352d0, 1308.8688d0, 1578.4224d0,
     :   8,-3, 1040.08018d0, 1321.61833d0, 1744.82517d0, 2100.34092d0,
     :   8, 3, 1040.16769d0, 1321.70472d0, 1745.37068d0, 2101.52639d0,
     :   8,-4, 1344.88476d0, 1705.25488d0, 2205.39935d0, 2636.86682d0,
     :   8, 4, 1344.88559d0, 1705.2534d0, 2205.3855d0, 2636.79319d0,
     :   8,-5, 1695.94482d0, 2124.15956d0, 2689.13456d0, 3141.17035d0,
     :   8, 5, 1695.9433d0, 2124.15961d0, 2689.13475d0, 3141.17188d0,
     :   8,-6, 2083.77011d0, 2570.73389d0, 3168.37941d0, 3753.99427d0,
     :   8, 6, 2083.77011d0, 2570.73389d0, 3168.3794d0, 3753.99419d0,
     :   8,-7, 2501.19786d0, 3039.38026d0, 3717.14641d0, 4300.028d0,
     :   8, 7, 2501.19786d0, 3039.38026d0, 3717.14641d0, 4300.02803d0,
     :   8,-8, 2942.65216d0, 3525.76312d0, 4240.43496d0, 4859.15267d0,
     :   8, 8, 2942.65216d0, 3525.76312d0, 4240.43496d0, 4859.15267d0,
*
     :   9, 0, 695.5779d0, 690.83615d0, 693.23362d0, 709.15131d0,
     :   9,-1, 738.41142d0, 825.17693d0, 1021.14713d0, 1185.48662d0,
     :   9, 1, 790.46812d0, 880.73686d0, 1085.5462d0, 1263.74715d0,
     :   9,-2, 933.72919d0, 1125.01698d0, 1454.47595d0, 1733.08841d0,
     :   9, 2, 938.8095d0, 1125.63318d0, 1447.47245d0, 1714.61434d0,
     :   9,-3, 1181.56969d0, 1462.83985d0, 1884.60871d0, 2240.31775d0,
     :   9, 3, 1181.74061d0, 1463.01113d0, 1886.5794d0, 2242.63688d0,
     :   9,-4, 1486.9862d0, 1847.28937d0, 2347.32119d0, 2779.77075d0,
     :   9, 4, 1486.98828d0, 1847.28553d0, 2347.28555d0, 2779.60233d0,
     :   9,-5, 1838.6746d0, 2266.91431d0, 2832.26772d0, 3283.21735d0,
     :   9, 5, 1838.67198d0, 2266.91451d0, 2832.26859d0, 3283.22247d0,
     :   9,-6, 2227.13546d0, 2714.173d0, 3310.7158d0, 3897.46203d0,
     :   9, 6, 2227.13546d0, 2714.17301d0, 3310.71578d0, 3897.46161d0,
     :   9,-7, 2645.19368d0, 3183.4785d0, 3861.33062d0, 4444.35389d0,
     :   9, 7, 2645.19368d0, 3183.4785d0, 3861.33062d0, 4444.35421d0,
*
     :   10, 0, 848.13892d0, 843.11825d0, 846.97887d0, 867.96051d0,
     :   10,-1, 887.57773d0, 973.61835d0, 1169.52286d0, 1334.70123d0,
     :   10, 1, 950.72272d0, 1041.11393d0, 1247.81257d0, 1429.48207d0,
     :   10,-2, 1089.32731d0, 1280.70504d0, 1611.0429d0, 1891.23509d0,
     :   10, 2, 1096.74156d0, 1281.53706d0, 1600.69824d0, 1864.91668d0,
     :   10,-3, 1338.52966d0, 1619.40849d0, 2042.75489d0, 2394.83258d0,
     :   10, 3, 1338.83885d0, 1619.72358d0, 2043.09367d0, 2399.01873d0,
     :   10,-4, 1644.61399d0, 2004.77152d0, 2504.57059d0, 2937.86897d0,
     :   10, 4, 1644.61864d0, 2004.76251d0, 2504.48839d0, 2937.51165d0,
     :   10,-5, 1996.98147d0, 2425.18183d0, 2990.76616d0, 3440.50125d0,
     :   10, 5, 1996.97694d0, 2425.18243d0, 2990.76915d0, 3440.51594d0,
     :   10,-6, 2386.13345d0, 2873.18775d0, 3468.41246d0, 4056.32852d0,
     :   10, 6, 2386.13346d0, 2873.18774d0, 3468.41238d0, 4056.32684d0,
     :   10,-7, 2804.87894d0, 3343.21194d0, 4021.0854d0, 4603.88132d0,
     :   10, 7, 2804.87894d0, 3343.21194d0, 4021.08541d0, 4603.88583d0,
*
     :   11, 0, 1015.22542d0, 1010.15714d0, 1015.98885d0, 1043.03919d0,
     :   11,-1, 1051.32443d0, 1136.63654d0, 1332.61438d0, 1499.06609d0,
     :   11, 1, 1126.451d0, 1217.08686d0, 1425.98989d0, 1611.37243d0,
     :   11,-2, 1260.08037d0, 1451.6136d0, 1783.01478d0, 2065.10249d0,
     :   11, 2, 1270.45792d0, 1452.67348d0, 1768.36679d0, 2029.29418d0,
     :   11,-3, 1510.87645d0, 1791.21035d0, 2213.35765d0, 2563.53796d0,
     :   11, 3, 1511.40205d0, 1791.75613d0, 2214.81724d0, 2570.60273d0,
     :   11,-4, 1817.6864d0, 2177.59948d0, 2677.01936d0, 3111.03922d0,
     :   11, 4, 1817.69589d0, 2177.58009d0, 2676.84562d0, 3110.32986d0,
     :   11,-5, 2170.77898d0, 2598.85575d0, 3164.51121d0, 3612.86738d0,
     :   11, 5, 2170.7714d0, 2598.85732d0, 3164.51988d0, 3612.90469d0,
     :   11,-6, 2560.67361d0, 3047.66855d0, 3700.90767d0, 4230.42258d0,
     :   11, 6, 2560.67361d0, 3047.66851d0, 3700.90729d0, 4230.41678d0,
*
     :   12, 0, 1196.69378d0, 1191.8443d0, 1200.24104d0, 1234.40987d0,
     :   12,-1, 1229.57039d0, 1314.16491d0, 1510.39911d0, 1678.67604d0,
     :   12, 1, 1317.50488d0, 1408.53889d0, 1619.98985d0, 1811.088d0,
     :   12,-2, 1445.88311d0, 1637.65119d0, 1970.32376d0, 2254.6554d0,
     :   12, 2, 1459.90343d0, 1638.92638d0, 1950.31782d0, 2207.77842d0,
     :   12,-3, 1698.51669d0, 1978.11769d0, 2398.80942d0, 2746.08038d0,
     :   12, 3, 1699.36535d0, 1979.01738d0, 2401.65398d0, 2757.32942d0,
     :   12,-4, 2006.11413d0, 2365.6631d0, 2864.53293d0, 3299.15686d0,
     :   12, 4, 2006.13213d0, 2365.62419d0, 2864.19091d0, 3297.82978d0,
     :   12,-5, 2359.97324d0, 2787.82065d0, 3353.36809d0, 3800.15046d0,
     :   12, 5, 2359.96099d0, 2787.82439d0, 3353.39048d0, 3800.23657d0,
     :   12,-6, 2750.65788d0, 3237.49698d0, 3889.40513d0, 4419.55625d0,
     :   12, 6, 2750.65788d0, 3237.49686d0, 3889.40398d0, 4419.53857d0,
*
     :   13, 0, 1392.41119d0, 1388.07129d0, 1399.71302d0, 1442.05546d0,
     :   13,-1, 1422.23217d0, 1506.13449d0, 1702.85759d0, 1873.64083d0,
     :   13, 1, 1523.72741d0, 1615.34663d0, 1829.71675d0, 2024.44362d0,
     :   13,-2, 1646.62495d0, 1838.7212d0, 2172.89722d0, 2459.8503d0,
     :   13, 2, 1664.99893d0, 1840.16539d0, 2146.41734d0, 2400.45515d0,
     :   13,-3, 1901.34738d0, 2179.98846d0, 2598.70648d0, 2942.13749d0,
     :   13, 3, 1902.65933d0, 2181.41117d0, 2603.50891d0, 2959.14797d0,
     :   13,-4, 2209.80046d0, 2568.84453d0, 3066.97243d0, 3502.10005d0,
     :   13, 4, 2209.83251d0, 2568.77081d0, 3066.33811d0, 3499.74771d0,
     :   13,-5, 2564.46316d0, 2991.95233d0, 3557.18675d0, 4002.17321d0,
     :   13, 5, 2564.44407d0, 2991.96065d0, 3557.23984d0, 4002.35693d0,
     :   13,-6, 2955.98109d0, 3442.54631d0, 4093.06202d0, 4623.52521d0,
     :   13, 6, 2955.98108d0, 3442.54596d0, 4093.05884d0, 4623.47578d0,
*
     :   14, 0, 1602.25715d0, 1598.73121d0, 1614.38191d0, 1665.91049d0,
     :   14,-1, 1629.22456d0, 1712.47449d0, 1909.97465d0, 2084.07984d0,
     :   14, 1, 1744.95419d0, 1837.38088d0, 2055.06762d0, 2261.62953d0,
     :   14,-2, 1862.19078d0, 2054.7225d0, 2390.65775d0, 2680.63367d0,
     :   14, 2, 1885.64067d0, 2056.24646d0, 2356.56152d0, 2607.45684d0,
     :   14,-3, 2119.25608d0, 2396.66581d0, 2812.74498d0, 3151.45897d0,
     :   14, 3, 2121.20985d0, 2398.83722d0, 2820.289d0, 3176.01406d0,
     :   14,-4, 2428.64139d0, 2787.01885d0, 3284.19696d0, 3719.7535d0,
     :   14, 4, 2428.69547d0, 2786.88585d0, 3283.07906d0, 3715.78678d0,
     :   14,-5, 2784.14079d0, 3211.11806d0, 3775.7986d0, 4218.74345d0,
     :   14, 5, 2784.11201d0, 3211.13539d0, 3775.91856d0, 4219.11068d0,
*
*  extrapolated energies appended below
     :  15, 0,  1837.893d0,  1831.509d0,  1846.037d0,  1900.550d0,
     :  16, 0,  2082.452d0,  2075.401d0,  2092.100d0,  2154.238d0,
     :  17, 0,  2342.297d0,  2334.538d0,  2353.541d0,  2423.781d0,
     :  18, 0,  2617.426d0,  2608.917d0,  2630.362d0,  2709.180d0,
     :  19, 0,  2907.841d0,  2898.540d0,  2922.562d0,  3010.434d0,
     :  20, 0,  3213.540d0,  3203.406d0,  3230.140d0,  3327.544d0,
     :  21, 0,  3534.525d0,  3523.515d0,  3553.097d0,  3660.509d0,
     :  22, 0,  3870.794d0,  3858.868d0,  3891.434d0,  4009.329d0,
     :  23, 0,  4222.348d0,  4209.464d0,  4245.149d0,  4374.005d0,
     :  24, 0,  4589.188d0,  4575.304d0,  4614.243d0,  4754.537d0,
     :  25, 0,  4971.312d0,  4956.386d0,  4998.716d0,  5150.924d0,
*
     :  15,-1,  1855.854d0,  1937.418d0,  2133.584d0,  2306.056d0,
     :  16,-1,  2094.612d0,  2174.989d0,  2371.041d0,  2544.886d0,
     :  17,-1,  2348.292d0,  2427.408d0,  2623.339d0,  2798.644d0,
     :  18,-1,  2616.895d0,  2694.675d0,  2890.477d0,  3067.328d0,
     :  19,-1,  2900.420d0,  2976.790d0,  3172.458d0,  3350.939d0,
     :  20,-1,  3198.867d0,  3273.753d0,  3469.279d0,  3649.478d0,
     :  21,-1,  3512.237d0,  3585.565d0,  3780.941d0,  3962.943d0,
     :  22,-1,  3840.529d0,  3912.224d0,  4107.444d0,  4291.335d0,
     :  23,-1,  4183.743d0,  4253.732d0,  4448.788d0,  4634.654d0,
     :  24,-1,  4541.880d0,  4610.088d0,  4804.974d0,  4992.900d0,
     :  25,-1,  4914.939d0,  4981.292d0,  5176.000d0,  5366.073d0,
*
     :  15, 1,  1990.009d0,  2081.697d0,  2301.142d0,  2508.518d0,
     :  16, 1,  2246.507d0,  2338.381d0,  2560.808d0,  2774.156d0,
     :  17, 1,  2519.035d0,  2611.107d0,  2836.703d0,  3056.397d0,
     :  18, 1,  2807.595d0,  2899.876d0,  3128.828d0,  3355.241d0,
     :  19, 1,  3112.186d0,  3204.687d0,  3437.182d0,  3670.686d0,
     :  20, 1,  3432.808d0,  3525.542d0,  3761.765d0,  4002.735d0,
     :  21, 1,  3769.461d0,  3862.439d0,  4102.577d0,  4351.385d0,
     :  22, 1,  4122.145d0,  4215.378d0,  4459.618d0,  4716.638d0,
     :  23, 1,  4490.860d0,  4584.361d0,  4832.888d0,  5098.494d0,
     :  24, 1,  4875.607d0,  4969.386d0,  5222.388d0,  5496.951d0,
     :  25, 1,  5276.384d0,  5370.454d0,  5628.116d0,  5912.012d0,
*
     :  15,-2,  2098.697d0,  2290.929d0,  2627.458d0,  2918.587d0,
     :  16,-2,  2347.632d0,  2540.008d0,  2877.949d0,  3171.585d0,
     :  17,-2,  2612.125d0,  2804.654d0,  3144.095d0,  3440.396d0,
     :  18,-2,  2892.176d0,  3084.868d0,  3425.898d0,  3725.018d0,
     :  19,-2,  3187.786d0,  3380.649d0,  3723.356d0,  4025.453d0,
     :  20,-2,  3498.954d0,  3691.997d0,  4036.469d0,  4341.701d0,
     :  21,-2,  3825.681d0,  4018.913d0,  4365.239d0,  4673.761d0,
     :  22,-2,  4167.966d0,  4361.396d0,  4709.663d0,  5021.633d0,
     :  23,-2,  4525.810d0,  4719.447d0,  5069.744d0,  5385.317d0,
     :  24,-2,  4899.212d0,  5093.065d0,  5445.481d0,  5764.814d0,
     :  25,-2,  5288.173d0,  5482.250d0,  5836.873d0,  6160.124d0,
*
     :  15, 2,  2122.569d0,  2292.818d0,  2593.099d0,  2841.156d0,
     :  16, 2,  2375.087d0,  2542.165d0,  2838.407d0,  3082.665d0,
     :  17, 2,  2643.388d0,  2807.096d0,  3099.047d0,  3339.268d0,
     :  18, 2,  2927.471d0,  3087.612d0,  3375.019d0,  3610.965d0,
     :  19, 2,  3227.336d0,  3383.712d0,  3666.323d0,  3897.757d0,
     :  20, 2,  3542.984d0,  3695.396d0,  3972.959d0,  4199.643d0,
     :  21, 2,  3874.415d0,  4022.664d0,  4294.926d0,  4516.623d0,
     :  22, 2,  4221.627d0,  4365.516d0,  4632.226d0,  4848.697d0,
     :  23, 2,  4584.622d0,  4723.953d0,  4984.857d0,  5195.866d0,
     :  24, 2,  4963.400d0,  5097.974d0,  5352.819d0,  5558.129d0,
     :  25, 2,  5357.960d0,  5487.579d0,  5736.114d0,  5935.487d0,
*
     :  15,-3,  2356.764d0,  2634.363d0,  3051.666d0,  3391.236d0,
     :  16,-3,  2607.747d0,  2884.666d0,  3300.964d0,  3637.854d0,
     :  17,-3,  2874.416d0,  3150.613d0,  3565.843d0,  3899.886d0,
     :  18,-3,  3156.772d0,  3432.204d0,  3846.303d0,  4177.331d0,
     :  19,-3,  3454.814d0,  3729.439d0,  4142.344d0,  4470.190d0,
     :  20,-3,  3768.543d0,  4042.318d0,  4453.966d0,  4778.463d0,
     :  21,-3,  4097.958d0,  4370.841d0,  4781.169d0,  5102.149d0,
     :  22,-3,  4443.059d0,  4715.008d0,  5123.953d0,  5441.249d0,
     :  23,-3,  4803.847d0,  5074.818d0,  5482.319d0,  5795.763d0,
     :  24,-3,  5180.322d0,  5450.273d0,  5856.265d0,  6165.690d0,
     :  25,-3,  5572.483d0,  5841.372d0,  6245.793d0,  6551.030d0,
*
     :  15, 3,  2122.569d0,  2292.818d0,  2593.099d0,  2841.156d0,
     :  16, 3,  2375.087d0,  2542.165d0,  2838.407d0,  3082.665d0,
     :  17, 3,  2643.388d0,  2807.096d0,  3099.047d0,  3339.268d0,
     :  18, 3,  2927.471d0,  3087.612d0,  3375.019d0,  3610.965d0,
     :  19, 3,  3227.336d0,  3383.712d0,  3666.323d0,  3897.757d0,
     :  20, 3,  3542.984d0,  3695.396d0,  3972.959d0,  4199.643d0,
     :  21, 3,  3874.415d0,  4022.664d0,  4294.926d0,  4516.623d0,
     :  22, 3,  4221.627d0,  4365.516d0,  4632.226d0,  4848.697d0,
     :  23, 3,  4584.622d0,  4723.953d0,  4984.857d0,  5195.866d0,
     :  24, 3,  4963.400d0,  5097.974d0,  5352.819d0,  5558.129d0,
     :  25, 3,  5357.960d0,  5487.579d0,  5736.114d0,  5935.487d0,
*
     :  15,-4,  2666.868d0,  3025.495d0,  3522.994d0,  3946.037d0,
     :  16,-4,  2918.824d0,  3277.159d0,  3774.215d0,  4192.057d0,
     :  17,-4,  3186.526d0,  3544.551d0,  4041.138d0,  4453.453d0,
     :  18,-4,  3469.976d0,  3827.673d0,  4323.761d0,  4730.225d0,
     :  19,-4,  3769.173d0,  4126.524d0,  4622.087d0,  5022.373d0,
     :  20,-4,  4084.117d0,  4441.103d0,  4936.113d0,  5329.898d0,
     :  21,-4,  4414.808d0,  4771.412d0,  5265.841d0,  5652.799d0,
     :  22,-4,  4761.247d0,  5117.449d0,  5611.270d0,  5991.076d0,
     :  23,-4,  5123.432d0,  5479.216d0,  5972.400d0,  6344.729d0,
     :  24,-4,  5501.365d0,  5856.712d0,  6349.232d0,  6713.759d0,
     :  25,-4,  5895.045d0,  6249.936d0,  6741.765d0,  7098.165d0,
*
     :  15, 4,  2666.913d0,  3025.388d0,  3522.087d0,  3942.720d0,
     :  16, 4,  2918.876d0,  3277.034d0,  3773.151d0,  4188.169d0,
     :  17, 4,  3186.587d0,  3544.407d0,  4039.906d0,  4448.957d0,
     :  18, 4,  3470.046d0,  3827.508d0,  4322.353d0,  4725.087d0,
     :  19, 4,  3769.252d0,  4126.336d0,  4620.491d0,  5016.557d0,
     :  20, 4,  4084.206d0,  4440.893d0,  4934.321d0,  5323.367d0,
     :  21, 4,  4414.907d0,  4771.177d0,  5263.842d0,  5645.518d0,
     :  22, 4,  4761.356d0,  5117.189d0,  5609.055d0,  5983.009d0,
     :  23, 4,  5123.553d0,  5478.929d0,  5969.959d0,  6335.841d0,
     :  24, 4,  5501.498d0,  5856.397d0,  6346.554d0,  6704.014d0,
     :  25, 4,  5895.190d0,  6249.593d0,  6738.842d0,  7087.527d0,
*
     :  15,-5,  3023.092d0,  3450.362d0,  4016.353d0,  4457.502d0,
     :  16,-5,  3275.979d0,  3703.093d0,  4269.370d0,  4708.423d0,
     :  17,-5,  3544.672d0,  3971.620d0,  4538.201d0,  4975.027d0,
     :  18,-5,  3829.171d0,  4255.942d0,  4822.845d0,  5257.313d0,
     :  19,-5,  4129.474d0,  4556.060d0,  5123.303d0,  5555.282d0,
     :  20,-5,  4445.584d0,  4871.973d0,  5439.575d0,  5868.933d0,
     :  21,-5,  4777.498d0,  5203.682d0,  5771.660d0,  6198.267d0,
     :  22,-5,  5125.219d0,  5551.187d0,  6119.558d0,  6543.284d0,
     :  23,-5,  5488.744d0,  5914.488d0,  6483.270d0,  6903.983d0,
     :  24,-5,  5868.075d0,  6293.584d0,  6862.796d0,  7280.365d0,
     :  25,-5,  6263.212d0,  6688.476d0,  7258.135d0,  7672.429d0,
*
     :  15, 5,  3023.065d0,  3450.376d0,  4016.443d0,  4457.792d0,
     :  16, 5,  3275.947d0,  3703.109d0,  4269.477d0,  4708.766d0,
     :  17, 5,  3544.635d0,  3971.638d0,  4538.325d0,  4975.427d0,
     :  18, 5,  3829.128d0,  4255.963d0,  4822.988d0,  5257.773d0,
     :  19, 5,  4129.426d0,  4556.084d0,  5123.465d0,  5555.805d0,
     :  20, 5,  4445.529d0,  4872.001d0,  5439.758d0,  5869.523d0,
     :  21, 5,  4777.438d0,  5203.713d0,  5771.864d0,  6198.927d0,
     :  22, 5,  5125.151d0,  5551.221d0,  6119.786d0,  6544.017d0,
     :  23, 5,  5488.670d0,  5914.525d0,  6483.522d0,  6904.793d0,
     :  24, 5,  5867.994d0,  6293.625d0,  6863.072d0,  7281.255d0,
     :  25, 5,  6263.123d0,  6688.521d0,  7258.438d0,  7673.403d0  /
*
      zero = 0.d0
      two = 2.d0
*  check for consistency in the values of flaghf and csflag
      if (flaghf) then
        write (6, 5)
        write (9, 5)
5      format (' *** FLAGHF = .TRUE. FOR SINGLET SYSTEM; ABORT ***')
        stop
      end if
      if (csflag) then
        write (6, 6)
        write (9, 6)
6     format
     :   ('  *** CSFLAG = .TRUE. FOR ASYMMETRIC TOP',
     :    ' - LINEAR MOLECULE COLLISIONS NOT IMPLEMENTED.  ABORT ***')
        stop
      end if
      if (flagsu) then
        write (6, 10)
        write (9, 10)
10      format ('  *** FLAGSU = .TRUE. FOR CH2(X 3B1) (0,V2,0)',
     :     ' SURFACE COLLISIONS; NOT IMPLEMENTED.  ABORT ***')
        call exit
      end if
      xjtot = jtot
      if (nlammx .lt. nlam) then
        write (6, 40) nlam, nlammx
        write (9, 40) nlam, nlammx
40      format (' ** TOTAL NUMBER OF ANISOTROPIC TERMS=', i2,
     :          ' .GT. NLAMMX=', i2,'; ABORT')
        stop
      end if
      if (ivbend.lt.0 .or. ivbend.gt.3) then
        write (6,47) ivbend
        write (9,47) ivbend
47      format(' *** V_BEND=', i3,' OUT OF ALLOWED RANGE',
     :         '     ABORT ***')
        stop
      end if
*
*  check that ipotsy2 is either 1 or 2
      if (ipotsy2.lt.1 .or. ipotsy2.gt.2) then
        write(6,27) ipotsy2
        write(9,27) ipotsy2
27      format(' *** IPOTSY2 .NE. 1 .OR. 2:  ABORT ***')
        if (bastst) then
          return
        else
          call exit
        end if
      end if
*
*  check that iop is either 1 or -1 if ihomo=T
      if (ihomo) then
        if ( .not. (iop.eq.1 .or. iop.eq.-1) ) then
          write(6,29) 
          write(9,29)
29        format(' *** IOP .ne. 1 .or. -1 with IHOMO = true:'
     :      '   ABORT ***')  
          if (bastst) then
            return
          else
            call exit
          end if
        end if
      end if
      if (bastst) write (6, 46) nlam
      write (9, 46) nlam
46    format (' ** TOTAL NUMBER OF ANISOTROPIC TERMS IN POTENTIAL =',
     :        i4)
      if (clist) then
        if (csflag) then
          if (ihomo) then
            if (bastst) then
              write (6,65) rmu * xmconv, ivbend,
     :             ipotsy, iop, ered * econv, jtot, nu
            end if
            write (9,65) rmu * xmconv, ivbend,
     :             ipotsy, iop,ered * econv, jtot, nu
65          format(/,' **  CH2(X 3B1) (0,V_BEND,0) VIBRONIC LEVEL **',
     :        /,'     RMU=', f9.4,'    V_BEND=',i3/,'     POT-SYM=', i2,
     :       '  O/P=',i2,'  E=', f7.2,'  LBAR=', i5, 2x,' NU=', i3)
          else
            if (bastst) then
              write (6,75) rmu * xmconv,
     :              ipotsy, ered * econv, jtot, nu
            end if
            write (9,75) rmu * xmconv, ivbend,
     :              ipotsy, ered * econv, jtot, nu
75          format(/,' **  CH2(X 3B1) (0,V_BEND,0) VIBRONIC LEVEL **',
     :        /,'     RMU=', f9.4,'    V_BEND=',i3/,'     POT-SYM=', i2,
     :        '     E=', f7.2,'  LBAR=', i5, 2x,' NU=', i3)
          end if
        else
          if (ihomo) then
            if (bastst) then
              write (6,80) rmu * xmconv, ivbend, ipotsy,
     :             iop, ered * econv, jtot, jlpar
            end if
            write (9,80) rmu * xmconv, ivbend, ipotsy,
     :              iop, ered * econv, jtot, jlpar
80          format(/,' **  CH2(X 3B1) (0,V_BEND,0) VIBRONIC LEVEL **',
     :        /,'     RMU=', f9.4,'    V_BEND=',i3/
     :        '     POT-SYM=', i2,'  O/P=',i2,
     :        '  E=', f7.2, '  JTOT=', i4, 2x,' JLPAR=', i2)
          else
            if (bastst) then
              write (6,85) rmu * xmconv, ivbend,
     :             ered * econv, jtot, jlpar
            end if
            write (9,85) rmu * xmconv, ivbendt,
     :             ered * econv, jtot, jlpar
85          format(/,' **  CH2(X 3B1) (0,V_BEND,0) VIBRONIC LEVEL **',
     :        /,'     RMU=', f9.4,'    V_BEND=',i3/,
     :        '     E=', f7.2, '  JTOT=', i4, 2x,
     :        ' JLPAR=', i2)
          end if
        end if
        if (.not. flagsu) write (9,90) rcut
90      format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=',
     :          f6.2)
      end if
*
*  set up list of all CH2 j(kp,ko) states included in basis
*  for j < jmax and e < emax
      ivsub = ivbend + 3
      nlev = 0
      do i=1, num_x
        if (ch2x_e(1,i).le.jmax .and. ch2x_e(ivsub,i).le.emax) then
          nlev = nlev + 1
          etemp(nlev) = ch2x_e(ivsub,i)
          fjtemp(nlev) = ch2x_e(1,i)
          fktemp(nlev) = abs(ch2x_e(2,i))
          fistmp(nlev) = sign(1.d0,ch2x_e(2,i))
        end if
      end do
*
*  check to see if (ki/isi) corresponds to an allowed ortho/para level for
*  a symmetric molecule.  the statements below correspond to CH2 in its
*  ground 3B1 electronic state
      if (ihomo) then
        nlist = 0
        do 110  njk = 1, nlev
          ki = fktemp(njk)
          ji = fjtemp(njk)
          isi = fistmp(njk)
          if (isi .eq. 1) then
            iss = 0
          else
            iss = 1
          end if
          iph = ji + ki + iss
          if (ki .ne. 2*(ki/2)) iph = iph + 1
          iexsym = -(-1) ** iph
          if (iexsym .ne. iop) go to 110
*  here if this state is to be included
          nlist = nlist + 1
          etemp(nlist) = etemp(njk)
          fktemp(nlist) = fktemp(njk)
          fjtemp(nlist) = fjtemp(njk)
          fistmp(nlist) = fistmp(njk)
110     continue
      else
        nlist = nlev
      end if
*
*  calculate energy levels for linear molecule collision partner
*  subtract energy of j2min level
      i2num = 0
      do j2 = j2min, j2max, ipotsy2
        i2num = i2num + 1
        j2rot(i2num) = j2
        e2rot(i2num) = b2rot * (j2 * (j2 + 1) - j2min * (j2min + 1.d0))
      end do
*
*  combine asymmetric top and linear molecule levels
      nlevel = 0
      do i2 = 1, i2num
        do 160 i1 = 1, nlist
          nlevel = nlevel + 1
          ehold(nlevel) = (etemp(i1) + e2rot(i2)) / econv
          jhold(nlevel) = 10 * fjtemp(i1) + j2rot(i2)
          ishold(nlevel) = fktemp(i1) * fistmp(i1)
160     continue
      end do
*
*  now sort this list in terms of increasing energy
      if (nlevel .gt. 1) then
        do 120 i1 = 1, nlevel - 1
          esave = ehold(i1)
          do 115 i2 = i1 + 1, nlevel
            if (ehold(i2) .lt. esave) then
*  state i2 has a lower energy than state i1, switch them
              esave = ehold(i2)
              ehold(i2) = ehold(i1)
              ehold(i1) = esave
              jsave = jhold(i2)
              jhold(i2) = jhold(i1)
              jhold(i1) = jsave
              isave = ishold(i2)
              ishold(i2) = ishold(i1)
              ishold(i1) = isave
            end if
115       continue
120     continue
      end if
*
*  print this list if bastst = .true. or if clist = .true.
      if (bastst .or. clist) then
        write (6, 130)
        write (9, 130)
130     format (/2x,
     :   'LEVEL LIST SORTED BY ENERGY',/,'   #   N  ',
     :     'IS  KP  KO   S  J2    EINT(CM-1)')
        do i = 1, nlevel
          ecm = ehold(i) * econv
          j1 = jhold(i)/10
          j2 = mod(jhold(i),10)
          isi = sign(ishold(i),1)
          if (isi.eq.1) then
            slab = '+'
          else
            slab = '-'
          endif
*  compute kp and ko projection quantum numbers
          kp = abs(ishold(i))
          if (ishold(i) .ge. 0) then
            ko = j1 - kp
          else
            ko = j1 + 1 - kp
          end if
            write (6, 135) i, j1, ishold(i),
     :        kp, ko, slab, j2, ecm
            write (9, 135) i, j1, ishold(i),
     :        kp, ko, slab, j2, ecm
135         format (5i4, 3x, a1, i4, f12.3)
        end do
      end if
*
*  here for cc calculations.  first calculate range of orbital angular
*  momentum quantum numbers allowed for this state
*  determine parity of molecular state [Eq. (A21) of S. Green, J. Chem. Phys.
*  73, 2740 (1980) and Townes and Schawlow, Microwave Spectroscopy, Eq. (3-27),
*  p. 64.]  See Eq. (A3) of Green.
*  CH2(X 3B1) state requires an additional minus sign
      n = 0
      do 170 i = 1, nlevel
        j1 = jhold(i)/10
        j2 = mod(jhold(i),10)
        ki = abs(ishold(i))
        iss = sign(1, ishold(i))
        ipar = iss * (-1) ** (j1 + j2 + ki)
        j12mn = abs(j1 - j2)
        j12mx = j1 + j2
        do j12i = j12mn, j12mx
          lmin = abs(jtot - j12i)
          lmax = jtot + j12i
          do li = lmin, lmax
*  check to see if this channel has the desired parity, if so, include it
            ix = ipar * (-1) ** (li - jtot)  
            if (ix .eq. jlpar) then
              n = n + 1
              if (n .gt. nmax) then
                write (6, 150) n, nmax
                write (9, 150) n, nmax
150               format(/' *** NCHANNELS=', i5,
     :              ' .GT. MAX DIMENSION OF',i5,' ABORT ***')
                stop
              end if
              eint(n) = ehold(i) 
              j(n) = jhold(i)
              is(n) = ishold(i)
              l(n) = li
              j12(n) = j12i
              cent(n) = li * (li + 1)
            end if
          end do
        end do
170   continue
*
*  print this scattering basis list if bastst = .true. or if clist = .true.
      if (bastst .or. clist) then
        write (6, 131) jtot, jlpar
        write (9, 131) jtot, jlpar
131     format (/,6x,'JTOT = ',i4,'  JLPAR = ',i2/,2x,
     :     'CH2(X)-H2 SCATTERING BASIS LIST',/,'   #   N  ',
     :     'IS  KP  KO   S  J2  J12  L    EINT(CM-1)')
        do i = 1, n
          ecm = eint(i) * econv
          j1 = j(i)/10
          j2 = mod(j(i),10)
          isi = sign(is(i),1)
          if (isi.eq.1) then
            slab='+'
          else
            slab='-'
          endif
*  compute kp and ko projection quantum numbers
          kp = abs(is(i))
          if (is(i) .ge. 0) then
            ko = j1 - kp
          else
            ko = j1 + 1 - kp
          end if
            write (6, 136) i, j1, is(i),
     :        kp, ko, slab, j2, j12(i),l(i), ecm
            write (9, 136) i, j1, is(i),
     :        kp, ko, slab, j2, j12(i),l(i), ecm
136         format (5i4, 3x, a1, 3i4, f12.3)
        end do
      end if
*
*  also determine number of levels which are open
      nlevop = 0
      do 250  i = 1, nlevel
        if (ehold(i) .le. ered) then
          nlevop = nlevop + 1
        end if
250   continue
*  now check to see if any of the open channels are closed at r=rcut
*  this is not done for molecule-surface collisions or for rcut < 0
      if (.not.flagsu .and. rcut .gt. 0.d0 .and. .not.boundc) then
        emin = 1. e+7
        do 290  i = 1, n
          if (eint(i) .le. ered) then
*  here if channel is open asymptotically
            if ( jtot * (jtot + 1) / (two * rmu * rcut * rcut)
     :          .gt. (ered - eint(i)) ) then
*  here if channel is open asymptotically but closed at r = rcut
              if (eint(i) .lt. emin) emin = eint(i)
*  emin now contains the lowest channel energy for which this condition is met
            end if
          end if
290     continue
*  now eliminate all channels with eint .ge. emin if any of the channels
*  are open asymptotically but closed at r = rcut
        if (emin.lt.ered) then
          nn = 0
          do 300 i = 1, n
            if (eint(i) .lt. emin) then
*  here if this channel is to be included
              nn = nn + 1
              eint(nn) = eint(i)
              j(nn) = j(i)
              is(nn) = is(i)
              cent(nn) = cent(i)
              l(nn) = l(i)
            end if
300       continue
*  reset number of channels
          n = nn
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
          write (6, 303) nu, n, ntop
          write (9, 303) nu, n, ntop
303       format (' *** NCH = ',i3, ' AT NU = ',i2, ' .GT. NTOP = ',i3,
     :            '; ABORT **',/,
     :    '     CHECK RCUT')
          call exit
        end if
      end if
*
*  now calculate coupling matrix elements
*  i counts v2 elements
*  inum counts v2 elements for given lambda
*  ilam counts number of v2 matrices
*  ij is address of given v2 element in present v2 matrix
      i = 0
      if (bastst) then
        write (6, 340)
        write (9, 340)
340     format (/' ILAM  L1   M1   L2  LTOT  ICOL  IROW    I'
     :    '     IV2      VEE')
      end if
      i = 0
      do 400 ilam = 1, nlam
        ll1 = lms(ilam)%l1
        mm1 = lms(ilam)%m1
        ll2 = lms(ilam)%l2
        lltot = lms(ilam)%ltot
        inum = 0
        do 355 icol = 1, n
          j1c = j(icol)/10
          j2c = mod(j(icol),10)
          j12c = j12(icol)
          lc = l(icol)
          do 350 irow = icol, n
            ij = ntop * (icol - 1) + irow
            j1r = j(irow)/10
            j2r = mod(j(irow),10)
            j12r = j12(irow)
            lr = l(irow)
            call vch2h2(j1r,j2r,j12r,lr,j1c,j2c,j12c,lc,
     :          is(irow),is(icol),jtot,
     :          ll1,mm1,ll2,lltot,vee)
            if (abs(vee) .gt. 1.d-10) then
              i = i + 1
              if (i .le. nv2max) then
                inum = inum + 1
                v2(i) = vee
                iv2(i) = ij
                if (bastst .and. iprint.eq.2) then
                  write (6, 291) ilam, ll1, mm1, ll2, lltot,
     :              icol, irow, i, iv2(i), vee
                  write (9, 291) ilam, ll1, mm1, ll2, lltot,
     :              icol, irow, i, iv2(i), vee
291               format (i4, 4i5, 2x, 2i5, i6, i7, e18.7)
                end if
              end if
            end if
350       continue
355     continue
        if (i .le. nv2max) lamnum(ilam) = inum
        if (bastst) then
          write (6, 370) ilam, lamnum(ilam)
          write (9, 370) ilam, lamnum(ilam)
370       format ('ILAM=',i4,' LAMNUM(ILAM) = ',i7)
        end if
        lamsum = lamsum + lamnum(ilam)
400   continue
      if (i .gt. nv2max) then
        write (6, 410) i, nv2max
        write (6, 410) i, nv2max
410     format (' *** NUMBER OF NONZERO V2 ELEMENTS = ',i10,
     :      ' .GT. NV2MAX=',i10,'; ABORT ***')
      end if
      if (bastst) then
        write (6, 420) lamsum
        write (9, 420) lamsum
420     format (' *** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS',
     :      i10)
      end if
      return
      end
* --------------------------------------------------------------------
      subroutine vch2h2(j1r,j2r,j12r,lr,j1c,j2c,j12c,lc,
     :  isr,isc,jtot,l1,m1,l2,ltot,vee)
*
*  subroutine to calculate v-lambda matrices for close-coupled 
*  treatment of collisions of an asymmetric top with a diatomic molecule
*  the angular dependence of the terms in the potential is given by
*  Valiron et al. [JCP 129, 134306 (2008)].  the angular expansion
*  coefficients have been muitiplied by the normalization factor in
*  the subroutine loapot
*
*  an expression for the matrix element in the space-fixed basis is given
*  by Phillips et al. [JCP 102, 6024 (1995)].
*
*  this subroutine is based on vastp3, in hibaaztp3.f
*  current revision date:  17-jan-2021 by p.dagdigian
*  -----------------------------------------------------------------------
*  variables in call list:
*    j1c:      rotational quantum number of #1 of left side of matrix element (bra)
*    j2c:      rotational quantum number of #2 of left side of matrix element (bra)
*    j12c:     value of j12 for right side of matrix element (bra)
*    lc:       orbital angular momentum of left side of matrix element (bra)
*    j1r:      rotational quantum number of #1 of right side of matrix element (ket)
*    j2r:      rotational quantum number of #2 of right side of matrix element (ket)
*    j12r:     value of j12 for right side of matrix element (ket)
*    lr:       orbital angular momentum of right side of matrix element (ket)
*    isr:      label for bra
*    isc:      label for ket
*    jtot:     total angular momentum
*    l1,m1,l2,ltot:   parameters for the given angular expansion coefficient       
*    vee:      on return, contains desired matrix element
*
*  subroutines called:
*    xf3j, xf6j, xf9j
*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      data one, zero, two / 1.0d0, 0.0d0, 2.0d0 /
      data fourpi / 12.566370614d0 /, sqr2 / 1.414213562373095d0 /
*
      vee = zero
      xjtot = jtot
      xj1r = j1r
      xj1c = j1c
      xj2r = j2r
      xj2c = j2c
      xj12r = j12r
      xj12c = j12c
      xlr = lr
      xlc = lc
      k1r = abs(isr)
      k1c = abs(isc)
      xk1r = k1r
      xk1c = k1c
      iepr = sign(1,isr)
      if (isr .eq. 0) iepr = 1
      iepc = sign(1,isc)
      if (isc .eq. 0) iepc = 1
*  get quantum numbers for this angular term
      xl1 = l1
      xm1 = m1
      xl2 = l2
      xltot = ltot
*  check phase relations first
      iphase = iepr * iepc * (-1) ** (j1r + j1c + m1 + l2 + ltot)
      if (iphase .eq. -1) return
      x1 = xf3j(xlr, xlc, xltot, zero, zero, zero)
      if (x1 .eq. 0.d0) return
      x2 = xf3j(xj2r, xl2, xj2c, zero, zero, zero)
      if (x2 .eq. 0.d0) return
      x3 = xf6j(xlr, xlc, xltot, xj12c, xj12r, xjtot)
      if (x3 .eq.0.d0) return
      x4 = xf9j(xj12c, xj12r, xltot, xj1c, xj1r, xl1,
     +  xj2c, xj2r, xl2)
      if (x4 .eq. 0.d0) return
*
      v = zero
      vk1 = xf3j(xj1r,xl1,xj1c, -xk1r, xm1, xk1c)
      vk2 = xf3j(xj1r,xl1,xj1c, xk1r, xm1, -xk1c) 
      vk3 = xf3j(xj1r,xl1,xj1c, -xk1r, xm1, -xk1c) 
      vk4 = xf3j(xj1r,xl1,xj1c, xk1r, xm1, xk1c) 
*  factor of 2 at end takes care of [1 + eps*epsp*(-1)^xx] term
      vk = (vk1 + iepr * iepc * vk2 + iepc * vk3
     +  + iepr * vk4) * two
*  special normalization for k1r and/or k1c = 0
      if (k1r .eq. 0) vk = vk / sqr2
      if (k1c .eq. 0) vk = vk / sqr2
      v = vk * (-1.d0) ** k1r / two
*
      iph = xjtot - xj1c + xj2c - xj12r - xltot
      ph = (-1.d0) ** iph
      sqjs = sqrt((two * xj1r + one) * (two * xj1c + one)
     +  * (two * xj2r + one) * (two * xj2c + one)
     +  * (two * xj12r + one) * (two * xj12c + one)
     +  * (two * xlr + one) * (two * xlc + one)
     +  * (two * xl2 + one) * (two * xltot + one))
*
      vee = ph * sqjs * x1 * x2 * x3 * x4 * v / fourpi
      if (m1 .eq. 0) vee = vee / two
*  include normalization factor depending on l1 
      delm1 = 1.d0
      if (m1 .eq. 0) delm1 = 2.d0
      fnorm = sqrt((two * xl1 + one) * delm1 / two)
      vee = vee * fnorm
*
      return
      end
* --------------------------------------------------------------------
      subroutine raise(mesg)
      implicit none
      character*(*) mesg
      write (0, *) 'hibridon: error: ', mesg
      stop
      end
* ----------------------------------------------------------------------
      subroutine datfln(filenm, fullnm)
      character (len=*) :: filenm, fullnm
      fullnm = 'potdata/' // trim(filenm)
      return
      end subroutine datfln
* ------------------------------------eof-------------------------------
