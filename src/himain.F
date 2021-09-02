      program logair
      use mod_clseg, only: allocate_clseg
      use mod_cobuf, only: allocate_cobuf
      use mod_cdbf, only: allocate_cdbf
      
      use mod_cosout, only: allocate_cosout
      use mod_coiout, only: allocate_coiout
      use mod_cov2, only: allocate_cov2
      use mod_coiv2, only: allocate_coiv2
      use mod_cocent, only: allocate_cocent
      use mod_coeint, only: allocate_coeint
      use mod_coj12, only: allocate_coj12
      use mod_coj12p, only: allocate_coj12p
      use mod_covvl, only: allocate_covvl
      use mod_cofact, only: allocate_cofact

      use mod_coener, only: allocate_coener

      use mod_coatpi, only: allocate_coatpi
      use mod_coatpr, only: allocate_coatpr
      use mod_coatp1, only: allocate_coatp1
      use mod_coatp2, only: allocate_coatp2
      use mod_coatp3, only: allocate_coatp3
      use mod_coz, only: allocate_coz, z
      use mod_cow, only: allocate_cow, w
      use mod_cozmat, only: allocate_cozmat, zmat
      use mod_coamat, only: allocate_coamat, amat
      use mod_cobmat, only: allocate_cobmat, bmat
      use mod_cotq1, only: allocate_cotq1, tq1
      use mod_cotq2, only: allocate_cotq2, tq2
      use mod_cotq3, only: allocate_cotq3, tq3
      use mod_cojq, only: allocate_cojq, jq
      use mod_colq, only: allocate_colq, lq
      use mod_coinq, only: allocate_coinq, inq
      use mod_cojhld, only: allocate_cojhld, jhold
      use mod_coehld, only: allocate_coehld, ehold
      use mod_coinhl, only: allocate_coinhl, inhold
      use mod_coisc1, only: allocate_coisc1, isc1
      use mod_coisc2, only: allocate_coisc2, isc2
      use mod_coisc3, only: allocate_coisc3, isc3
      use mod_coisc4, only: allocate_coisc4, isc4
      use mod_coisc5, only: allocate_coisc5
      use mod_coisc6, only: allocate_coisc6
      use mod_coisc7, only: allocate_coisc7
      use mod_coisc8, only: allocate_coisc8
      use mod_coisc9, only: allocate_coisc9
      use mod_coisc10, only: allocate_coisc10
      use mod_coisc11, only: allocate_coisc11
      use mod_coisc12, only: allocate_coisc12
      use mod_colsc1, only: allocate_colsc1
      use mod_cosc1, only: allocate_cosc1, sc1
      use mod_cosc2, only: allocate_cosc2, sc2
      use mod_cosc3, only: allocate_cosc3, sc3
      use mod_cosc4, only: allocate_cosc4, sc4
      use mod_cosc5, only: allocate_cosc5, sc5
      use mod_cosc6, only: allocate_cosc6, sc6
      use mod_cosc7, only: allocate_cosc7, sc7
      use mod_cosc8, only: allocate_cosc8, sc8
      use mod_cosc9, only: allocate_cosc9, sc9
      use mod_cosc10, only: allocate_cosc10
      use mod_cosc11, only: allocate_cosc11
      use mod_coeig, only: allocate_coeig
      use mod_coeig2, only: allocate_coeig2
#if defined(HIB_UNIX)
      use mod_cokaux, only: allocate_cokaux
#endif
      use mod_cotble, only: allocate_cotble
      use mod_coqvec, only: allocate_coqvec
      use mod_coqvec2, only: allocate_coqvec2
      use mod_codim, only: allocate_codim, mairy, mmax
      use mod_comxbs, only: allocate_comxbs
      use mod_comxm, only: allocate_comxm
***********************************************************************
****   this code is not released for general public use            ****
****   all use must be by specific prior arrangement with:         ****
****     millard alexander, department of chemistry,               ****
****     university of maryland, college park, md, 20742           ****
****     tel: 1.301.405.1823; email: mha@umd.edu                   ****
****   no part of this program may be copied or used for other     ****
****   purposes without the author's permission.                   ****
***********************************************************************
*  ***  driver for log-derivative integration ***
*  author:  millard alexander
*  revised:  24-jan-2012 by p. dagdigian (added common blocks /coj12/, /co12p/)
*  revised:  replaced statements before call to propag to accommodate
*    q.ma's revised version of bound (28-jun-2013, p.dagdigian)
*  revised: 16-jun-2019 by p. dagdigian (increased kmxbas to 30)
*  added common block from ba3p2s basis routine
*  current revision date:  19-jun-2019 (p. dagdigian)
*
      implicit double precision (a-h, o-z)
#if defined(HIB_UNIX_IBM) || defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
      character *40 test
#endif
*  ----------------------------------------------------------
*  in the following parameter statements:
*     kmxbas is maximum number of included basis routines, this should
*     be updatated as basis routines are added
*     kmax is maximum number of channels set at compile time
*     klammx is maximum number of anisotropic terms in potential
*     kfact is maximum value for which log(n!) is computed
*     ken is number of total energies allowed
*     kmxpho is maximum number of different initial states allowed in
*     photodissociation calculation
*     kout is number of different values of rotational quantum number
*     for which s-matrix will be stored on disk

*  change kmax below --  supposed to be set at compile time
*  revised by p. dagdigian (13-dec-2019)
      parameter (kmax=KMAX, kairy = kmax,ktri=kmax*(kmax+1)/2,kbig=10)
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
* set size of scratch array for matrix factorization and generalized
*   eigenvalue determination
      parameter (kaux3=max(66*kmax,2*kmax*kmax+1))
#endif
#if defined(HIB_UNIX_IBM)
* set size of scratch array for matrix inversion
      parameter (kaux=100*kmax)
#endif
#if defined(HIB_UNIX) && !defined(HIB_UNIX_IBM)
* set size of scratch array for matrix inversion with lapack routines
* warning, this assumes a blocksize of 64
      parameter (kaux=128*kmax)
#endif
* modified klammx (pjd - 17-jan-2019)
* increased kfact (pjd - 15-dec-2020)
      parameter (klammx = 2000, kfact = 3500, kout = 50, ken = 25)
      parameter (kmxpho = 3, knphot = 1)
      parameter (kmxbas = 30)
      parameter (kq2 = 2 * kmax, kqmax = kmxpho * kmax)
* kv2max sets the maximum size of the v2 matrix, a reasonable size is
* kmax**2
*
* increase kv2max from 2 * kmax * kmax to 20 * kmax * kmax (pjd - 13-dec-2019)
*      parameter (kv2max = 2 * kmax * kmax)
      parameter (kv2max = 50 * kmax * kmax)
*
* krotmx set the maximum size of the asym top e.fn expansion
* increased factor multiplying kmax below and narray from 12 to 100 (pjd:  19-sep-2017)
      parameter (krotmx = 500 * kmax)
*  variable in common block /conlam/
*    nlammx:    the maximum number of angular coupling terms allowed
*    nlam:      the total number of angular coupling terms used
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx
*  variables in common block /cofil/
*     nfl:      logical unit number
*     iofbuf:   buffer offset
*     maxrec:   number of records on file associated with nfl
*     iofrec:   offset in file
*  ----------------------------------------------------------
      common /comom/  xmom(3), imom(13)
      common /cofil/  nfl,iofbuf,maxrec(60),iofrec(60),nwrec
      common /conlam/ nlam, nlammx, lamnum(klammx)
*   square matrices and vectors
      call allocate_coener(ken)
      call allocate_clseg()
      call allocate_cobuf()
      call allocate_cdbf()
      call allocate_coatpi(krotmx, 100)
      call allocate_coatpr(krotmx)
      call allocate_coatp1(krotmx)
      call allocate_coatp2(krotmx)
      call allocate_coatp3(krotmx)
      call allocate_coz(kmax)
      call allocate_cow(kmax)
      call allocate_cozmat(kmax)
      call allocate_coamat(kmax)
      call allocate_cobmat(kairy)
      call allocate_cotq1(T_MATRIX_SIZE)
      call allocate_cotq2(T_MATRIX_SIZE)
      call allocate_cotq3(T_MATRIX_SIZE)
      call allocate_cojq(kmax)
      call allocate_colq(kmax)
      call allocate_coinq(kmax)
      call allocate_cojhld(kmax)
      call allocate_coehld(kmax)
      call allocate_coinhl(kmax)
      call allocate_coisc1(kmax)
      call allocate_coisc2(kmax)
      call allocate_coisc3(kmax)
      call allocate_coisc4(kmax)
      call allocate_coisc5(kmax)
      call allocate_coisc6(kmax)
      call allocate_coisc7(kmax)
      call allocate_coisc8(kmax)
      call allocate_coisc9(kmax)
      call allocate_coisc10(kmax)
      call allocate_coisc11(kmax)
      call allocate_coisc12(kmax)
      call allocate_colsc1(kmax)
      call allocate_cosc1(kmax)
      call allocate_cosc2(kmax)
      call allocate_cosc3(kmax)
      call allocate_cosc4(kmax)
      call allocate_cosc5(kmax)
      call allocate_cosc6(kmax)
      call allocate_cosc7(kmax)
      call allocate_cosc8(kmax)
      call allocate_cosc9(kmax)
      call allocate_cosc10(kmax)
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
      call allocate_cosc11(kaux3)
#endif
#if defined(HIB_UNIX) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
      call allocate_cosc11(kaux)
#endif
*   total matrix and vector storage required is:
*     7 kmax**2 + 25 kmax + kv2max + kfact -- without airy integration
*     8 kmax**2 + 25 kmax + kv2max + kfact -- with airy integration
*   if linked with -b option, then storage requirements drop to
*     5 kmax**2 + 25 kmax + kv2max + kfact -- with airy integration
*
*  parameter below sets maximum size of asymmetric top basis fn expansion
      call allocate_cosout(kout)
      call allocate_coiout(kout)
      call allocate_cov2(kv2max)
      call allocate_coiv2(kv2max)
      call allocate_cocent(kmax)
      call allocate_coeint(kmax)
      call allocate_coj12(kmax)
      call allocate_coj12p(kmax)
      call allocate_covvl(klammx)
      call allocate_cofact(kfact)
      call allocate_coeig()
      call allocate_coeig2()
#if defined(HIB_UNIX_IBM)
      call allocate_cokaux(anaux=max(kaux,1800))
#endif
      call allocate_cotble(kfact)
      call allocate_coqvec(kqmax, kmxpho, knphot)
      call allocate_coqvec2(kq2)
      call allocate_codim(kairy, kmax, kbig)
      call allocate_comxbs(kmxbas)
      call allocate_comxm()
*
      men = ken
      nlammx = klammx
*  calculate array containing logs of factorials for use in vector
*  coupling coefficient routines:  factlg(i)=log(fact(i-1))
      call factlg(kfact)
      call finit
* start scattering calculation
      call flow (z, w, zmat, amat, bmat, jq, lq, inq, jhold, ehold,
     :           inhold, isc1, isc2, isc3, isc4, lsc1,
     :           sc2, sc1, sc3, sc4, sc5,
     :           sc6, sc7, sc8, sc9, tq1, tq2, tq3, men, mmax, mairy)
      end
*  -------------------------------------------------------------
      subroutine flow (z, w, zmat, amat, bmat, jq, lq, inq, jlev,
     :            elev, inlev, isc1, isc2, isc3, isc4, lsc1,
     :            sc2, sc1, sc3, sc4,
     :            sc5, sc6, sc7, sc8, sc9, tq1, tq2, tq3, men,
     :            nmax, nairy)
*  -------------------------------------------------------------
*  program to control the log-derivative/airy integration
*  written by:  millard alexander
*  additions by: b. follmeg, h-j werner
*  current revision date:  1-oct-2001 by mha
*  -------------------------------------------------------------
*  variables in common block /coered/
*    ered:      collision energy in atomic units (hartrees)
*    rmu:       collision reduced mass in atomic units (mass of electron = 1)
*  variables in common block /copmat/
*    rtmn,rtmx: minimum and maximum turning points
*    iflag:     variable used in determination of turning points (not used her
*           iflag = 0 if all channels are in classically forbidden region
*           iflag = 1 if some channels are open
*           iflag = 2 if all asymptotically open channels are open at r
*  variables in common block /colpar/
*    airyfl:      if .true., then airy propagation will take place
*    airypr:      if .true., then step-by-step information is printed out in
*                 airy propagation
*    bastst:      if .true., then execution terminates after the first call
*                 to basis
*    batch:       if .true., then the job is run as a batch job
*                 if .false., then the job is assumed to be interactive
*    chlist:      if .true., then the channel quantum numbers and energies are
*                 printed out at at each total-j
*                 if .false., then  this is done only at first total-j
*    csflag:      if .true., then coupled-states calculation is desired
*                 if .false., then close-coupled calculation is desired
*    flaghf:      if .true., then the system has even multiplicity
*                 (half-integer total angular momentum)
*    flagsu:      if .true., then the problem is assumed to a molecule
*                 scattering off a surface, in which case the diagonal
*                 elements of the transition probabilities are equal to the
*                 modulus squared of the s-matrix (not t-matrix elements)
*    ihomo:       if .true., then the molecule is assumed to be homonuclear
*    ipos:        if .true., then printout is suited for a 132-position printe
*                 if .false., then printout is suited for a  80 -position
*                 printer
*    logdfl:      if .true., then logd propagation will take place
*    logwr:       if .true., then the lower triangle of log-derivative matrix
*                 is printed at the end of the logd and the end of the airy
*                 integrations
*    noprin:      if .true., then most printing is suppressed
*    partw:       if .true, then the full matrix of partial cross sections
*                 (summed over m-states) is printed
*    readpt:      if .true., then potential parameters are expected
*    rsflag:      if .true., then calculation is to be restarted
*                 a check will be made to see if all required files
*                 are present:  these may include
*                    trstrt, tmp10, tmp11, xsecn (or tmpxn), smatn,
*                    psecn, tmp35, ...
*    swrit:       if .true., then the upper triangle of real and imaginary
*                 parts of s-matrix are printed
*    t2test:      if .true., then the first two columns of the square modulus
*                 of the t-matrix are printed
*    t2writ:      if .true., then the upper triangle of square modulus of
*                 t-matrix is printed
*    twomol:      .true. if molecule-molecule collision
*    writs:       if .true., and nnout is > 0, then those s-matrix elements
*                 for which both the initial and final rotational quantum
*                 numbers are in the array jout (input line 12) are written to
*                 files smat1, smat2, ...
*                 if nnout < 0, then each column of the s-matrix whose initial
*                 index is in the array jout is written to files smat1, smat2,
*    wrpart:      if .true., then input data and the matrix of partial cross
*                 sections (summed over m-states) is written to file pxsec
*    wrxsec:      if .true., then some input data and the full matrix of
*                 integral cross sections ((summed over m-states and summed
*                 from jtot1 to jtot2) is written to file xsec1, xsec2, ....
*    xsecwr:      if .true., then the full matrix of integral cross sections
*                 ((summed over m-states and summed from jtot1 to jtot2) is
*                 printed
*
*  variable in common block /cosurf/
*    surffl:      this variable is set equal to flagsu, it is held in a
*                 separate common block for compatability with subroutines
*                 smatop, soutpt, and xwrite
*                 if .true., then the problem is assumed to a molecule
*                 scattering off a surface, in which case the diagonal
*                 elements of the transition probabilities are equal to the
*                 modulus squared of the s-matrix (not t-matrix elements)
*  variable in common block /cojlpo/
*    jlpold:      parity used in xwrite subroutine to insure correct
*                 accumulation of partial waves in cases where jlpar=0
*
*  variables in common block /cophot/
*     photof        true if photodissociation calculation
*                   false if scattering calculation
*     wavefn        true if g(a,b) transformation matrices are saved
*                   to be used later in computing the wavefunction
*  variables in module constants
*    econv:       conversion factor from cm-1 to hartree
*    xmconv:      conversion factor from amu to atomic units
*
*  variable in common block /coopti/
*    optifl:      flag, signals if the calculation is an optimization
*
      use mod_cosout, only: nnout, jout
      use mod_cocent, only: cent
      use mod_coeint, only: eint
      use mod_coener, only: energ
      use mod_hiba1sg, only : basis
      use constants
      implicit double precision (a-h,o-z)
      integer :: jtotmx
      character*20 cdate
      character*10 time
      character*10 timew,cpubaw,cpuptw,cpuaiw,cpuldw,cpusmw,cpuouw,
     :             cpuphw,timew1,timew2,time1,time2
      logical logwr, swrit, t2writ, wrpart, partw, airyfl, airypr,
     :        ipos, noprin, chlist, wrxsec, xsecwr, writs, csflag,
     :        flaghf, clist, rsflag, t2test, logdfl, flagsu,
     :        batch, readpt, ihomo, bastst, twojlp, firstj,
     :        twomol, surffl, nucros, ready, photof, photfl, wavefl,
     :        wavefn, boundc, boundf, wrsmat
*  -------------------------------------------------------------
      logical optifl, first
      logical lsc1
#include "common/parpot"
#if defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IRIS)
      real secnds
      common /codec/ ttim(2)
#endif
      common /cputim/ cpuld,cpuai,cpupot,cpusmt,cpupht
      common /cosavi/ iipar, ixpar(8)
      common /cosavr/ irpar(2), rxpar(9)
      common /copmat/ rtmn, rtmx, iflag
      common /coered/ ered, rmu
      common /coipar/ jtot1, jtot2, jtotd, jlpar, nerg,numax,numin,nud,
     :                lscreen, iprint
      common /corpar/ fstfac, rincr, rcut, rendai, rendld, rstart, spac,
     :                tolhi, xmu
      common /colpar/ airyfl, airypr, bastst, batch, chlist,
     :                csflag, flaghf, flagsu, ihomo, ipos, logdfl,
     :                logwr, noprin, partw, readpt, rsflag, swrit,
     :                t2test, t2writ, twomol, writs, wrpart, wrxsec,
     :                xsecwr, nucros, photof, wavefl, boundc
      common /cophot/ photfl, wavefn, boundf, wrsmat
      common /cosurf/ surffl
      common /coselb/ ibasty
      common /cojlpo/ jlpold
      common /coopti/ optifl
      common /constp/ nsteps, isteps
*   square matrices
      dimension z(nmax,nmax), w(nmax,nmax), zmat(nmax,nmax),
     :          amat(nmax,nmax), bmat(nairy,nairy)
*  vectors
      dimension jq(10), lq(10), inq(10), jlev (1), isc1(9), isc2(1),
     :          isc3(1), isc4(1), lsc1(5), inlev(1),
     :          elev(1), sc1(2), sc2(1), sc3(1), sc4(1), nlev(25),
     :          sc5(1), sc6(1), sc7(1), sc8(1), sc9(1), tq1(1),
     :          tq2(1), tq3(1)
      data twojlp / .false. /
*   to obtain timing information calls are made to system-specific
*   time subroutine:  call mtime(cpu,wallt) where cpu is clock cpu
*   time in seconds and wallt is wall clock time in seconds
*   user will have to change this subroutine for his own installation
*  subroutine to get input data
      first=.true.
*  get default data
      call default
1     call hinput(first)
      cpupt=0
#if defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IRIS)
      ttim(1)=0.d0
      ttim(2)=secnds(0.0)
#endif
      call mtime (tcpu0, twall0)
      photfl = photof
      wavefn = wavefl
      boundf = boundc
      wrsmat = writs
      surffl = flagsu
*  subroutine to open required i/o files
      if (.not.bastst) then
            call openfi (nerg)
      end if
      call version(9)
      call acknow(9,ipos)
      call dater(cdate)
      write (9, 15) label
15    format (/' ** LABEL:     ', (a))
      write (9, 16) potnam
16    format (' ** POT NAME:  ', (a))
      write (9, 20) cdate
20    format ( ' ** DATE:      ', (a))
      if (airyfl .and. (nairy .ne. nmax)) then
        write (9, 35)
        write (6, 35)
35      format(
     :     /' *** NAIRY .NE. NMAX FOR AIRY INTEGRATION; ABORT ***')
        call exit
      end if
      if (.not. airyfl .and. (nairy. ne. 1)) write (9, 45)
45    format(
     :  /' *** WARNING:  NAIRY .NE. 1 BUT AIRYFL .EQ. .FALSE.')
*  check for proper choice of output files and dimension limits
      if (nerg .gt. men) then
        write (9, 50)
        write (6, 50)
50      format (/' *** TOO MANY TOTAL ENERGIES; ABORT ***')
        call exit
      end if
      if (nerg.gt.19) then
         if (writs.and.(xsecwr.or.wrxsec))  then
            write (6,55) nerg
55          format(
     :    /' NERG = ',i2,' > 19',
     :    ' FOR WRITS .AND. (WRXSEC .OR. WRPART) ALL TRUE; ABORT ***')
            call exit
         endif
      endif
      write (9, 60) nmax,nairy
60    format (' ** NMAX=',i4,'  NAIRY=',i4)
*  convert collision energy, rotational constant, and reduced mass
*  to atomic units
      rmu = xmu / xmconv
      jtotmx = 0
*  reduce maximum jtot if all channels are closed (not for bound states)
      if (.not.boundc) then
        jtop = sqrt (energ(1) / econv * rcut * rcut * 2.d0 * rmu )
     :        - 0.5d0
        numj = (jtot2 - jtot1) / jtotd + 1
        jtotmx = jtot1 - jtotd
        do 75  jj = 1, numj
          jtotmx = jtotmx + jtotd
          if (jtotmx .ge. jtop) then
            jtotmx = jtotmx - jtotd
            go to 76
          end if
75      continue
      endif
76    firstj=.true.
      jfirst=jtot1
      jtop =jtotmx
      cpubas=0
      cpuout=0
      cpuai=0
      cpuld=0
      cpupot=0
      cpupht=0
      cpusmt=0
      jtotd=ixpar(3)
      jlpar=ixpar(4)
      nerg =ixpar(5)
      numax=ixpar(6)
      numin=ixpar(7)
      nud=  ixpar(8)
      nufirs=numin
      nulast=numax
      nutop=numax
      jlprsv=jlpar
      fstfac=rxpar(1)
      rincr=  rxpar(2)
      rcut=  rxpar(3)
      rendai=rxpar(4)
      rendld=rxpar(5)
      rstrt0=rxpar(6)
      spac=  rxpar(7)
      tolhi= rxpar(8)
      isteps=0
      dlogd = rendld - rstart
      xmu=rxpar(9)
      rtmnla=rstrt0
      dinsid=0
      twojlp=jlpar.eq.0.and..not.csflag
      if (twojlp) jlpar = 1
      jlpold = jlpar
      jfrest=0
74    jfirst=ixpar(1)
      if(jfrest.gt.0) jfirst=jfrest
      jtotmx=jtop
      if(nucros) nulast=numin
      if (jtot1.eq.jfirst .or. jtot1.eq.jfirst+1) rstart=rxpar(6)
80    jtot1 = jfirst
      if (.not.boundc) jtot2 = jtotmx
      nchmax = 0
*  this is beginning of loop over total angular momentum
*  (partial wave index)
      
100   jtot = jtot1
*  get restart parameters if run is to be restarted
      ready=.false.
      if (rsflag) then
*  read in restart information
        call restrt (jtot,jtopo,jtotd,jlpar,nu,nutop,nud,nerg,nlev,
     :         nchmax,rtmn1,rtmx1,dinsid,writs,csflag,nucros)
        rtmnla=rtmn1
        jlpold = jlpar
*  move s-matrix files to last partial wave done
        if(writs) then
           do 88 ifile=1,nerg
              nfile=44+ifile
              ered = energ(ifile)/econv
              nlevop=nlev(ifile)
              call wrhead(nfile, cdate,
     :                  sc1(1),  sc1(2), lsc1(1), lsc1(2), lsc1(3),
     :                 lsc1(4), lsc1(5), isc1(1), isc1(2), isc1(3),
     :                 isc1(4), isc1(5), isc1(6), isc1(7), isc1(8),
     :                 isc1(9), isc2, isc3, sc2, isc4)
              call fimovs(nfile,jtot,jlpar,nu,ifile,ierr)
              if(ierr.ne.0) then
                 write(6,87) jtot,jlpar
                 write(9,87) jtot,jlpar
87               format(/'  *** JTOT=',i5,', JLPAR=',i3,' NOT FOUND',
     :                   ' IN S-MATRIX FILE, RESTART IMPOSSIBLE')
                 call exit
              end if
88         continue
        end if
        firstj=.false.
        rsflag=.false.
        if(.not. csflag) then
          jtoto=jtot
          jtot=jtot+jtotd
          if(jtot.gt.jtot2) then
            if(jlpar.lt.0.or.jlprsv.gt.0) then
              write(6,89)
              write(9,89)
89            format(/'  *** RESTART, BUT NO MORE PARTIAL WAVES',
     :                ' REQUESTED ')
              ready=.true.
              goto 105
            end if
            jlpar=-1
            jtot=ixpar(1)
          else if(jlprsv.eq.0.and.jlpold.eq.-1.and.jtoto.eq.jtopo) then
            twojlp=.true.
            jlpar=1
            jlpold=1
            jfrest=jtot
          end if
          write (6, 90) jtot,jlpar
          write (9, 90) jtot,jlpar
90        format (' ** CONTINUE CC CALCULATION AT JTOT=',i5,
     :            ', JLPAR=',i5)
        else
          if(nucros) then
            nu=nu+nud
            if(nu.gt.numax) then
              write(6,89)
              write(9,89)
              ready=.true.
              goto 105
            end if
            numin=nu
            write(6,92) nu
            write(9,92) nu
92          format(/' ** CONTINUE CS CALCULATION AT NU=',i5)
            goto 74
          else
            nu=ixpar(7)
            jtot=jtot+jtotd
            if(jtot.gt.jtot2) then
              write(6,89)
              write(9,89)
              ready=.true.
              goto 105
            end if
            write(6,93) jtot
            write(9,93) jtot
93          format (' ** CONTINUE CS CALCULATION AT LBAR=',i5)
          end if
        end if
        jtot1=jtot
      end if
*
c      write (9, 95)
c95    format (1h ,79('='))
      xjtot = jtot
      if (flaghf) xjtot = xjtot + 0.5d0
105   continue
      rtmn1 = max (rendld, rendai)
      rtmx1 = 0.
*  this is beginning of loop over coupled states projection index
*  in the case of cc calculation this loop is executed only once
      nu = numin -  1
110   nu = nu + 1
      rtmx = 0.
      rtmn = max (rendld, rendai)
      ien = 1
*  this is beginning of loop over collision energies
115   ener = energ(ien)
      ered = ener/econv
      if(.not.ready) then
      if (csflag) then
        if (.not. flaghf) then
          if(partw.and..not.noprin) then
            write (9, 120) jtot, nu, ien, ener
            write (6, 120) jtot, nu, ien, ener
120         format (/' ** LBAR=',i4,'  NU=',i2,'  IEN =',i3,
     :             '  ENERGY (CM-1) =', f11.4)
          endif
        else
          if(partw.and..not.noprin) then
            write (9, 125) jtot,nu+0.5d0,
     :              ien,ener
            write (6, 125) jtot, nu+0.5d0, ien, ener
125         format (/' ** LBAR=',i4, '  NU=', f5.1, '  IEN =',i3,
     :             '  ENERGY (CM-1) =', f11.4)
          endif
        end if
      else
        if (.not. flaghf) then
          if(partw.and..not.noprin) then
            write (9, 130) jtot,jlpar,ien,ener
            write (6, 130) jtot, jlpar, ien, ener
130         format (/' ** JTOT =',i4,'  JLPAR =', i2,'  IEN =',i3,
     :             '  ENERGY (cm-1) =', f11.4)
          endif
        else
          if(partw.and..not.noprin) then
            write (9, 135) jtot+0.5d0,
     :      jlpar, ien, ener
            write (6, 135) jtot+0.5d0, jlpar, ien, ener
135         format (/' ** JTOT =',f7.1, '  JLPAR =', i2,'  IEN =',i3,
     :             '  ENERGY (cm-1) =', f11.4)
          endif
        end if
      end if
      end if
*  ered is collision energy in hartree
      eshift = 2.d0 * rmu * (ered - energ(1) / econv)
*  if first energy (ien = 1), then set up channel indices, array of
*  internal energies, and angular coupling matrices for each vlambda(r)
      call mtime (t1,t2)
      if (ien .eq. 1) then
        clist = .false.
        if (chlist .and.(jtot. eq. jfirst)) clist = .true.
        if (noprin) clist = .false.


**        write(6,5554) nch, nmax, nchtop
**5554    format('BEFORE BASIS:',3i6)


        call basis (jq, lq, inq, jlev, elev, inlev, nlevel, nlevop,
     :              sc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu,
     :              csflag, clist, bastst, ihomo, nu, numin, jlpar,
     :              twomol, nch, nmax, nchtop)


**        write(6,5556) nch, nmax, nchtop
**5556    format('AFTER BASIS: ',3i6)


        if (ready) goto 370
        if (bastst) then
           write(6,140)
           write(9,140)
  140      format(' ** BASTST=.TRUE.; TEST OF SUBROUTINE BASIS')
           goto 420
        end if
*  on return from basis:
*  nch is number of channels
*  nchtop is the maximum row dimension of all matrices passed to the
*  subroutines propag and soutpt
*  inq is the additional quantum index of each channel
*  jq is array of rotational quantum numbers
*  lq is array of orbital angular momenta
        if (nch .gt. nchmax) nchmax = nch
* in rare cases one might come back from basis with no open channels (but
* with closed channels present).  to deal with this case set nchop to the
* number of open channels
        nchop = 0
        if (nch .ge. 1) then
          do 144 ii = 1 , nch
            if (ered - eint(ii) .gt. 0.) nchop = nchop + 1
144       continue
        end if
        if (nchop .eq. 0) then
          if (csflag) then
            if (nu .eq. numin) then
              write (6, 145)
              write (9, 145)
145           format(
     :           /' *** NCH = 0 FOR NU=NUMIN IN CS CALCULATION;',
     :            ' END CALCULATION')
*  reset jtot2 to reflect absence of all channels in cs calculations
              jtot2 = jtot - jtotd
              write (6, 150) jtot2
              write (9, 150) jtot2
150           format (/' ** RESET JTOT2=', i4)
              go to 370
            else
              if (nu - 1 .lt. nutop) then
                nutop = nu - 1
                nulast=nutop
                write (6, 155) nu, nutop
155             format
     :           (' ** NCH = 0 FOR NU=',i3,' IN CS CALCULATION;',
     :            ' SET NUTOP=', i3, ' AT NEXT JTOT'/,
     :            '    MOVE ON TO NEXT PARTIAL WAVE')
              go to 300
              end if
            end if
          else
*  here for cc calculation, if no channels and jtot .le. 2 move on to
*  next partial wave, otherwise reset jtot2 and quit calculation
            if (jtot .le. 2) then
              write (9, 160)
              write (6, 160)
160           format (/' ** NCH = 0, MOVE ON TO NEXT PARTIAL WAVE')
*--------------------------------------------------------------------
* restart with next jtot
               jfirst = jfirst + 1
               if(jtotd.gt.1) jtotmx = jtotmx + 1
               goto 80
*--------------------------------------------------------------------
            else
              write (6, 165)
              write (9, 165)
165           format (/' *** NCH = 0 AND JTOT2.GE.2; END CALCULATION')
              jtot2 = jtot2 - jtotd
              write (6, 150) jtot2
              write (9, 150) jtot2
* Claire's modification: if nch=0 then eventually go to next parity:
              if (twojlp .and. jlpar .gt. 0) then
                jlpar = -1
                goto 74
              else
                go to 370
              end if
            end if
          end if
        end if
*  store channel parameters on unit 12 if this calculation is to be
*  performed at a second energy
        if (nerg .gt. 1) then
* open file for storage of transformation matrices
          rewind (12)
          write (12, 170) nch
170       format (i4)
          if (nch .gt. 0) then
            do 180  i = 1, nch
              write (12, 175) jq(i), lq(i), inq(i), cent(i), eint(i)
175           format (3i6, 2e25.15)
180         continue
          end if
        end if
      else if(ien.gt.1) then
        rewind (12)
        read (12, 170) nch
        if (nch .gt. 0) then
          do 250  i = 1, nch
            read (12, 175) jq(i), lq(i), inq(i), cent(i), eint(i)
250       continue
        end if
      end if
      irec=(ien-1)*5+2
      nlevop=0
      do 255 i=1,nlevel
255   if(elev(i).le.ered) nlevop=nlevop+1
      nlev(ien)=nlevop
      if (ien .gt. 1) then
        if (nlev(ien) .ne. nlev(ien-1)) then
          write (6, 256) ien, nlev(ien), ien-1, nlev(ien-1)
256       format (' *** NLEV(',i2,') = ',i3,' .NE. NLEV(',i2,') = ',i3,
     :           ' ABORT ***')
          call exit
        endif
      endif
      if (wrxsec .or. xsecwr .or. partw .or. wrpart) then
        if (.not. wavefl .and. .not. photof) then
          call dres(nlevop**2+4,1,irec)
          call dres(nlevop**2+4,1,irec+1)
          call dres(nlevop**2+4,1,irec+2)
          call dres(nlevop**2+4,1,irec+3)
          call dres(nlevop**2+4,1,irec+4)
          if(nucros) then
            irec=(nerg+ien-1)*5+2
            call dres(nlevop**2+4,1,irec)
            call dres(nlevop**2+4,1,irec+1)
            call dres(nlevop**2+4,1,irec+2)
            call dres(nlevop**2+4,1,irec+3)
            call dres(nlevop**2+4,1,irec+4)
          end if
          call dsave(1)
        endif
      endif
* store header and reserve space on direct access file 2 if wavefunction
* is desired (not if bound state calculation)
      if (wavefl .and. .not.boundc) 
     :  call wavewr(jtot,jlpar,nu,nch,nchtop,rstart,rendld)
      call mtime (t11, t22)
      tb =  t11 - t1
      tbm = t22 - t2
      cpubas=cpubas + tb
*  at subsequent partial waves, adjust starting point for integration for
*  next partial wave to be dinsid inside of innermost classical turning
*  point but no larger than rendld
*  10/14/08 mha:  but don't allow an increase greater than spac/2 in rstart
      if (jtot .gt. jfirst) then
        rstart = max(rtmnla-dinsid,rstrt0)
        rstart = min(rstart,rstrt0+0.5*spac)
        rstrt0 = rstart
*  adjust starting point of airy integration to be no less than rstart + dlogd
*  but no larger than rendai
        rendld = min (rstart + dlogd, rendai)
        if (.not. logdfl) rendld = rstart
      end if

cger (next 2 lines)
      write (9,'(/" ** J =",i5," JLPAR =",i2," STARTED")') jtot,jlpar
#if defined(HIB_UNIX_IBM) || defined(HIB_UNIX_AIX)  || defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
*      call flush_(9)
      call flush(9)
#endif
#if defined(HIB_UNIX_DEC) || defined(HIB_UNIX_CONVEX) || defined(HIB_UNIX_IRIS)
      call flush(9)
#endif
* replace statement below to accommodate q.ma's revised
* version of bound (28-jun-2013, p.dagdigian)
* do not reduce maximum size of matrices if bound state calculation
*      if (boundc) then
*        ntop=nmax
*      else
*        ntop=nchtop
*      endif
      ntop = nchtop
*
#ifdef DEBUG_ISSUE35
      write(6, *) 'issue35/flow : z(1,1)=', z(1,1)
#endif
      call propag (z, w, zmat, amat, bmat,
     :             jq, lq, inq, isc1, sc1, sc2, sc3, sc4, sc5, sc6, sc7,
     :             sc8, sc9,
     :             ien, nerg, ered, eshift, rstart, rendld, spac,
     :             tolhi, rendai, rincr, fstfac, tb, tbm,
     :             ipos, logwr, noprin, airyfl, airypr,
     :             nch, nopen, nairy, ntop)
* if bound state calculation, end it now
      if (boundc) then
        endfile (9)
        close (9)
        goto 1
      endif
*  now print out s-matrix and t-matrix squared, and calculate partial
*  cross sections and print them out, if desired
      call soutpt (z, w, zmat, amat,
     :             lq, jq, inq, isc1, isc2, bmat, tq1,
     :             jlev, elev, inlev, jtot, jfirst,
     :             jtot2, jtotd, nu, numin, nulast, nud, jlpar, ien,
     :             ipos, csflag, flaghf, swrit, t2writ, t2test,
     :             writs, wrpart, partw, wrxsec, xsecwr, twomol,
     :             nucros, firstj, nlevel, nlevop, nopen, nchtop,
     :             twojlp)
      cpuout = cpuout + second() - t11
      
*  on return from soutpt:
*     if wrxsec,xsecwr, wrpart, and partw are all .false., the upper-left
*     nopen x nopen block of z contains the modulus squared of the t-matrix
*     otherwise, the upper nlevop x nlevop block of z contains the partial
*     cross sections
*     the upper-left nopen x nopen block of w contains the real part of
*     the s-matrix
*     the upper-left nopen x nopen block of zmat contains the imaginary part
*     of the s-matrix
*     the arrays eint, cent, jq, lq, inq have been packed to eliminate
*     the closed-channel components
      if (ien.eq.1 .and. nchop.gt.0) then
        rtmn1 = min(rtmn1,rtmn)
        rtmx1 = max(rtmx1,rtmx)
      end if
      ien = ien + 1
      call mtime (tcpuf, twallf)
      if (ien.eq.2) then
         tcpu1=tcpuf-tcpu0
         twall1=twallf-twall0
      endif
*  go back to start calculation at another energy
      if (ien .le. nerg) go to 115
*  first partial wave has been calculated for all energies
      if (firstj) firstj = .false.
*  go back to start calculation at another value of coupled-states projection
*  index
300   if (nu .lt. nulast. and. nu. lt. nutop) go to 110
      if( .not.nucros) nulast = nutop
*  if first partial wave, then set distance inside turning point at which
*  logd integration starts
      if (jtot .eq. jfirst)  then
        dinsid = rtmn1 - rstart
        write (9, 330) dinsid
330     format (/' ** INTEGRATION WILL START', f6.3,
     :           ' BOHR INSIDE INNER TURNING POINT')
      end if
*  save last min and max turning points for next partial wave
      rtmnla = rtmn1
      rtmxla = rtmx1
      if(partw.and..not.nucros) write (9, 350)
350   format (1h ,79('='))
      if(.not.nucros) then
c.....save restart information
        if (wrxsec .or. xsecwr .or. partw .or. wrpart) then
          if (.not. wavefl .and. .not. photof) then
           call rsave (jtot,jtop,jtotd,jlpar,nu,nutop,nud,nerg,nlev,
     :                  nchmax,rtmn1,rtmx1,dinsid,writs,csflag,
     :                  nucros)
          endif
        endif
        call mtime (tcpuf, twallf)
        tcpuf = tcpuf - tcpu0
        twallf = twallf - twall0
        if (nerg.gt.1) then
            tcpu2=(tcpuf-tcpu1)/(nerg-1)
            twall2=(twallf-twall1)/(nerg-1)
        endif
        call gettim(tcpuf,time)
        call gettim(tcpu1,time1)
        call gettim(twallf,timew)
        call gettim(twall1,timew1)
        if (nerg.gt.1) then
           call gettim(tcpu2,time2)
           call gettim(twall2,timew2)
        endif
        call gettim(cpubas,cpubaw)
        call gettim(cpuld,cpuldw)
        call gettim(cpuai,cpuaiw)
        call gettim(cpupht,cpuphw)
        call gettim(cpusmt,cpusmw)
        call gettim(cpupot,cpuptw)
        call gettim(cpuout,cpuouw)
        call dater(cdate)
        if (.not. noprin) then
          write (6, 360) jtot,jlpar,cpubaw,cpuptw,cpuldw,cpuaiw,cpuphw,
     :       cpusmw,cpuouw,time,timew,cdate
360       format (/' ** J =', i5,' JLPAR =', i2,' COMPLETED'/1x,
     :   'CPU-TIMES:',
     :   '  BASIS:',a,'  POT:',a,'  LOGD: ',a,'  AIRY: ',a/11x,
     :   '  PSI0:',a,'  SMAT:',a,'  SOUT: ',a,'  CUMULATIVE:',a,/11x,
     :   '  ELAPSED: ',a,'  CURRENT DATE:  ',a)
          write (6, 361) rtmn, rtmx
361       format (' TURNING POINTS (MIN/MAX) =', 2(f8.2) )
        else
          write (6, 362) jtot, jlpar, time,timew,cdate
          write (9, 362) jtot, jlpar, time,timew,cdate
362       format (' ** J =', i5,' JLPAR =', i2,
     :     ' FINISHED; ',
     :     ' CPU:',a,'  WALL:',a,'  DATE: ',a)
        endif
#if defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IRIS) || defined(HIB_UNIX_HP)
        call flush6
#endif
        cpubas=0
        cpuout=0
        cpuai=0
        cpupht=0
        cpuld=0
        cpupot=0
        cpusmt=0
      end if
*  go back to start calculation at the next partial wave
      jtot1 = jtot + jtotd
      if (jtot1 .le. jtot2) go to 100
      if (twojlp .and. jlpar .gt. 0) then
        jlpar = -1
        rstrt0=rxpar(6)
        rstart=rstrt0
        rtmnla=rstrt0
        dinsid=0
        goto 74
      end if
c.....next nu value if nu runs in outer loop
      if (nucros .and. .not. wavefl .and. .not. photof) then
        call nusum (z, tq1, tq2, tq3,
     :              jlev,elev, inlev, jtot, jfirst,
     :              jtop, jtotd, nu, nufirs, numax, nud, jlpar,
     :              nerg, ipos, csflag, flaghf, wrpart, partw,
     :              twomol, nucros, nlevel, nlev, nopen, nmax)
c.....save restart information
        if (wrxsec .or. xsecwr .or. partw .or. wrpart) then
          if (.not. wavefl .and. .not. photof) then
            call rsave (jtot,jtop,jtotd,jlpar,nu,nutop,nud,nerg,nlev,
     :                  nchmax,rtmn1,rtmx1,dinsid,writs,csflag,
     :                  nucros)
          endif
        endif
        call mtime (tcpuf, twallf)
        tcpuf = tcpuf - tcpu0
        twallf = twallf - twall0
        call gettim(tcpuf,time)
        call gettim(twallf,timew)
        call gettim(cpubas,cpubaw)
        call gettim(cpuld,cpuldw)
        call gettim(cpuai,cpuaiw)
        call gettim(cpupht,cpuphw)
        call gettim(cpusmt,cpusmw)
        call gettim(cpupot,cpuptw)
        call gettim(cpuout,cpuouw)
        call dater(cdate)
        if (.not. noprin) then
          write (6, 366) nu,cpubaw,cpuptw,cpuldw,cpuaiw,cpuphw,
     :   cpusmw,
     :   cpuouw,time,timew,cdate
366       format (/' ** NU =', i3,' COMPLETED'/1x,
     :   'CPU-TIMES:',
     : '  BASIS:',a,'  POT: ',a,'  LOGD: ',a,'  AIRY:    ',a/11x,
     : '  PHOTO:',a,'  SMAT: ',a,'  SOUT:',a,'  CUMULATIVE:',a,
     :   '  ELAPSED: ',a,'  CURRENT DATE:  ',a)
          write (6, 361) rtmn, rtmx
        else
          write (6, 367) nu, time,timew,cdate
          write (9, 367) nu, time,timew,cdate
367       format (' ** NU =', i2,
     :     ' COMPLETED; ',
     :      ' CPU:',a,'  ELAPSED:',a,'  CURRENT DATE: ',a)
        endif
#if defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IRIS) || defined(HIB_UNIX_HP)
        call flush6
#endif
        cpubas=0
        cpuld=0
        cpuai=0
        cpupht=0
        cpusmt=0
        cpupot=0
        cpuout=0
        if(nulast.lt.numax) then
          numin = numin + nud
          goto 74
        end if
      end if
*  calculation has now been done at all partial waves and, if
*  desired, at both values of jlpar
*  write out integral cross sections if desired
 370  if (xsecwr .or. wrxsec) then
        if(.not.bastst) 
     :  call xwrite (amat, tq3, jlev, elev, inlev, nerg, energ,
     :             jfirst, jtot2, jtotd, csflag, flaghf,
     :             wrxsec, xsecwr, ipos, twomol, nucros, nlevel,
     :             nlev, nufirs, nulast, nud, jlpar, nchtop, nmax,
     :             ihomo)
       
      endif
      if (.not. bastst .and.
     :   xsecwr .or. wrxsec .or. wrpart .or. partw) then
        do 400 ien = 1, nerg
          nfile = 70 + ien
          close (nfile)
          if (wrpart) then
            nfile = 24 + ien
            close (nfile)
            if (csflag .and.
     :          (wrpart .or. partw .or. wrxsec .or. xsecwr)) then
              nfile = 34 + ien
              close (nfile)
            end if
          end if
400     continue
      end if
      if (.not. bastst .and.
     :   xsecwr .or. wrxsec .or. wrpart .or. partw) then
        if (.not. wavefl .and. .not. photof) then
          call dclos(1)
         endif
      endif
      if (wavefl .and. .not. boundc) close(22)
      if (writs) then
        do 410 ien = 1, nerg
        nfile = 44 + ien
        call closf (nfile)
        close (nfile)
410     continue
      end if
      if (nerg .gt. 1) then
        if (airyfl) close (10)
        close (11)
      end if
420   call dater (cdate)
      if (.not. optifl) then
         write (6, 350)
         write (6,500) nchmax
         write (9,500) nchmax
500      format(' **** END OF CALCULATION ****',/,
     :    '      MAXIMUM NUMBER OF CHANNELS USED WAS:  ',i4)
         if (nerg.eq.1) then
            write (6,505) timew, time, 100d0*(tcpuf/twallf)
            write (9,505) timew, time, 100d0*(tcpuf/twallf)
505         format('      TIMING:  ELAPSED',(a),'/ CPU',(a),
     :             '/ MP RATIO',f7.1,' %')
         endif
         if (nerg.gt.1) then
            write(6,510) timew1, time1, 100d0*(tcpu1/twall1)
            write(9,510) timew1, time1, 100d0*(tcpu1/twall1)
510         format('      TIMING (FIRST ENERGY):  ELAPSED',(a),
     :            '/ CPU',(a),
     :           '/ MP RATIO',f7.1)
            write(6,515) timew2, time2, 100d0*(tcpu2/twall2)
            write(9,515) timew2, time2, 100d0*(tcpu2/twall2)
515         format('      TIMING (SECOND ENERGY): ELAPSED',(a),
     :           '/ CPU',(a),
     :           '/ MP RATIO',f7.1)
         endif
         write (6,520) cdate
         write (9,520) cdate
520      format('      CURRENT DATE:  ',(a))
         write (6, 350)
         write (9, 350)
#if defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IRIS) || defined(HIB_UNIX_HP)
         call flush6
#endif
      end if
      close (23)
      endfile (9)
      close (9)
      goto 1
      end
