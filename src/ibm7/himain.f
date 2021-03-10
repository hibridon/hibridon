cstart unix-ibm unix-aix unix-darwin
@process noopt
cend
cstart unix-hp
c;!$hp$optimize off
cend
      program logair
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
*  current revision date:  23-feb-2004 by mha
      implicit double precision (a-h, o-z)
cstart unix-ibm unix-darwin
      character *40 test
cend
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

      parameter (kmax=151, kairy = kmax,ktri=kmax*(kmax+1)/2,kbig=10)
cstart unix-darwin
c;* set size of scratch array for matrix diagonalization
c;      parameter (kaux3=3*kmax)
cend
cstart unix-ibm
* set size of scratch array for matrix inversion
      parameter (kaux=100*kmax)
cend
cstart unix .and. .not.unix-ibm
c;* set size of scratch array for matrix inversion with lapack routines
c;* warning, this assumes a blocksize of 64
c;      parameter (kaux=128*kmax)
cend
      parameter (klammx = 80, kfact = 2000, kout=21, ken = 10)
      parameter (kmxpho=3, knphot=1)
      parameter (kmxbas=14)
      parameter (kq2=2*kmax,kqmax=kmxpho*kmax)
* kv2max sets the maximum size of the v2 matrix, a reasonable size is
* kmax**2
      parameter (kv2max= 2* kmax * kmax)
*  variable in common block /cov2/
*    nv2max:    maximum core memory allocated for the v2 matrix
*    ndummy:    dummy variable for alignment
*    v2:        lower triangle of nonzero angular coupling matrix elements
*               stored in packed column form that is (1,1), (2,1), (3,1) ...
*               (n,1), (2,2), (3,2) ... (n,2), etc.
*  variable in common block /coiv2/
*   lamnum:     number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx
*   iv2:        row+column index of v2 matrix for each non-zero element
*  variable in common block /conlam/
*    nlammx:    the maximum number of angular coupling terms allowed
*    nlam:      the total number of angular coupling terms used
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx
*  variable in common block /cocent/
*    cent:      array containing centrifugal barrier of each channel
*  variable in common block /coeint/
*    eint:      array containing channel energies (in hartree)
*  variable in common block /coener/
*    energ:     array containing total energies (in cm-1) at which scattering
*               calculations are to be performed
*  variable in common block /covvl/
*    vvl        array to store r-dependence of each angular term in the
*               potential
*  variable in common block /cofact/
*    si         table of logarithms of factorials
*  variables in common block /cosout/
*    nnout:     number of different rotational levels for which s-matrix
*               elements are to be saved in file 13,
*    jout(i):   values of rotational angular momentum for these lvels
*  variables in common block /clseg/
*    lseg:      sector length on mass storage device
*    intrel:    number of integer words per real words
*    lchar:     number of characters per integer word
*  variables in common block /cobuf/
*    lbuf:      buffer length (should be multiple of lseg)
*    ibuf:      buffer used in several i/o routines
*  variables in common block /cofil/
*     nfl:      logical unit number
*     iofbuf:   buffer offset
*     maxrec:   number of records on file associated with nfl
*     iofrec:   offset in file
*  variables in common block /cotble/
*     npnt      max. number of pointer
*     jttble    array containing pointer to records in s-matrix
*  variable in common block /coqvec/
*     mxphot    maximum number of initial states (maximum columns in q)
*     nphoto    actual number of initial states used
*     q         accumulated overlap matrix with ground state
*               only calculated if photof = .true.
*  variable in common block /comxbs/
*     maxbas    maximum number of allowed basis types
*  ----------------------------------------------------------
      logical lsc1
      common /comom/  xmom(3), imom(13)
      common /cosout/ nnout, jout(21)
* should be jout(kout)
      common /cov2/ nv2max, ndummy, v2(kv2max)
      common /coiv2/ iv2(kv2max)
      common /cocent/ cent(kmax)
      common /coeint/ eint(kmax)
      common /covvl/  vvl(klammx)
      common /cofact/ si(kfact)
      common /coener/ energ(ken)
      common /clseg/  lseg,intrel,lchar
      common /cobuf/  lbuf,ibuf(1)
      common /cofil/  nfl,iofbuf,maxrec(60),iofrec(60),nwrec
      common /conlam/ nlam, nlammx, lamnum(klammx)
*   square matrices and vectors
      common /coz/ z(kmax,kmax)
      common /cow/ w(kmax,kmax)
      common /cozmat/ zmat(kmax,kmax)
      common /coamat/ amat(kmax,kmax)
      common /cobmat/ bmat(kairy,kairy)
* these matrices to store t matrices
      common /cotq1/ tq1(kmax,kmax)
      common /cotq2/ tq2(kmax,kmax)
      common /cotq3/ tq3(kmax,kmax)
      common /cojq/ jq(kmax)
      common /colq/ lq(kmax)
      common /coinq/ inq(kmax)
      common /cojhld/ jhold(kmax)
      common /coehld/ ehold(kmax)
      common /coinhl/ inhold(kmax)
      common /coisc1/ isc1(kmax)
      common /coisc2/ isc2(kmax)
      common /coisc3/ isc3(kmax)
      common /coisc4/ isc4(kmax)
      common /coisc5/ isc5(kmax)
      common /coisc6/ isc6(kmax)
      common /coisc7/ isc7(kmax)
      common /coisc8/ isc8(kmax)
      common /colsc1/ lsc1(kmax)
      common /cosc1/ sc1(kmax)
      common /cosc2/ sc2(kmax)
      common /cosc3/ sc3(kmax)
      common /cosc4/ sc4(kmax)
      common /cosc5/ sc5(kmax)
      common /cosc6/ sc6(kmax)
      common /cosc7/ sc7(kmax)
      common /cosc8/ sc8(kmax)
      common /cosc9/ sc9(kmax)
      common /cosc10/ sc10(kmax)
cstart unix-darwin
c;      common /cosc12/ sc11(kaux3)
cend
cstart unix .and. .not.unix-darwin
      common /cosc11/ sc11(kaux)
cend
cstart unix
      common /cokaux/ naux
cend
      common /cotble/ npnt, jttble(kfact)
      common /coqvec/ mxphot, nphoto, q(kqmax)
      common /coqvec2/ q2(kq2)
      common /codim/ mairy,mmax,mbig
      common /comxbs/ maxbas
      common /comxm/ ncache, mxmblk
*   total matrix and vector storage required is:
*     7 kmax**2 + 25 kmax + kv2max + kfact -- without airy integration
*     8 kmax**2 + 25 kmax + kv2max + kfact -- with airy integration
*   if linked with -b option, then storage requirements drop to
*     5 kmax**2 + 25 kmax + kv2max + kfact -- with airy integration

      mairy = kairy
      mmax = kmax
      mbig=kbig
cstart unix-ibm
      naux=max(kaux,1800)
cend
      mxphot = kmxpho
      nphoto = knphot
      maxbas = kmxbas
      men = ken
      npnt = kfact
      nv2max = kv2max
      nlammx = klammx
*  calculate array containing logs of factorials for use in vector
*  coupling coefficient routines:  factlg(i)=log(fact(i-1))
      call factlg(kfact)
      call finit
*  determine cache and block sizes for matrix multiply
c.....ncache: cache size in words
      ncache=4096
cstart unix-hp
c;      ncache=16000
cend
      mxmblk=64
* start scattering calculation
      call flow (z, w, zmat, amat, bmat, jq, lq, inq, jhold, ehold,
     :           inhold, isc1, isc2, isc3, isc4, lsc1,
     :           sc2, sc1, sc3, sc4, sc5,
     :           sc6, sc7, sc8, sc9, tq1, tq2, tq3, men, mmax, mairy)
      end
cstart unix-ibm unix-aix unix-darwin
@process noopt
cend
cstart unix-hp
c;!$hp$optimize off
cend
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
*  variable in common block /coeint/
*    eint:      array containing channel energies (in hartree)
*  variable in common block /coener/
*    energ:     array containing total energies (in cm-1) at which scattering
*               calculations are to be performed
*  variable in common block /cocent/
*    cent:      array containing centrifugal barrier of each channel
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
*  variables in common block /coconv/
*    econv:       conversion factor from cm-1 to hartree
*    xmconv:      conversion factor from amu to atomic units
*
*  variable in common block /coopti/
*    optifl:      flag, signals if the calculation is an optimization
*
      implicit double precision (a-h,o-z)
      character*20 cdate
      character*13 time,timew,cpubaw,cpuptw,cpuaiw,cpuldw,cpusmw,cpuouw,
     :             cpuphw
      logical logwr, swrit, t2writ, wrpart, partw, airyfl, airypr,
     :        ipos, noprin, chlist, wrxsec, xsecwr, writs, csflag,
     :        flaghf, clist, rsflag, t2test, logdfl, flagsu,
     :        batch, readpt, ihomo, bastst, twojlp, firstj,
     :        twomol, surffl, nucros, ready, photof, photfl, wavefl,
     :        wavefn, boundc, boundf
*  -------------------------------------------------------------
      logical optifl, first
      logical lsc1
      include "common/parpot"
cstart unix-dec unix-iris
c;      real secnds
c;      common /codec/ ttim(2)
cend
      common /cputim/ cpuld,cpuai,cpupot,cpusmt,cpupht
      common /cosavi/ iipar, ixpar(8)
      common /cosavr/ irpar(2), rxpar(9)
      common /copmat/ rtmn, rtmx, iflag
      common /cosout/ nnout, jout(21)
      common /coeint/ eint(1)
      common /cocent/ cent(1)
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
      common /cophot/ photfl, wavefn, boundf
      common /cosurf/ surffl
      common /coselb/ ibasty
      common /coener/ energ(1)
      common /cojlpo/ jlpold
      common /coconv/ econv, xmconv
      common /coopti/ optifl
      common /constp/ nsteps, isteps
*   square matrices
      dimension z(nmax,nmax), w(nmax,nmax), zmat(nmax,nmax),
     :          amat(nmax,nmax), bmat(nairy,nairy)
*  vectors
      dimension jq(10), lq(10), inq(10), jlev (1), isc1(9), isc2(1),
     :          isc3(1), isc4(1), lsc1(5), inlev(1),
     :          elev(1), sc1(2), sc2(1), sc3(1), sc4(1), nlev(10),
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
cstart unix-dec unix-iris
c;      ttim(1)=0.d0
c;      ttim(2)=secnds(0.0)
cend
      call mtime (tcpu0, twall0)
      photfl = photof
      wavefn = wavefl
      boundf = boundc
      surffl = flagsu
*  subroutine to open required i/o files
      if (.not.bastst) call openfi (nerg)
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
35      format
     :    (/' *** NAIRY .NE. NMAX FOR AIRY INTEGRATION; ABORT ***')
        call exit
      end if
      if (.not. airyfl .and. (nairy. ne. 1)) write (9, 45)
45    format
     : (/' *** WARNING:  NAIRY .NE. 1 BUT AIRYFL .EQ. .FALSE.')
*  check for proper choice of output files and dimension limits
      if (nerg .gt. men) then
        write (9, 55)
        write (6, 55)
55      format (/' *** TOO MANY TOTAL ENERGIES; ABORT ***')
        call exit
      end if
      write (9, 60) nmax,nairy
60    format (' ** NMAX=',i4,'  NAIRY=',i4)
*  convert collision energy, rotational constant, and reduced mass
*  to atomic units
      rmu = xmu / xmconv
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
      rstart=rxpar(6)
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
        call basis (jq, lq, inq, jlev, elev, inlev, nlevel, nlevop,
     :              sc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu,
     :              csflag, clist, bastst, ihomo, nu, numin, jlpar,
     :              nch, nmax, nchtop)
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
145           format
     :          (/' *** NCH = 0 FOR NU=NUMIN IN CS CALCULATION;',
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
      if (jtot .gt. jfirst) then
        rstart = max (rtmnla - dinsid, rstrt0)
*  adjust starting point of airy integration to be no less than rstart + dlogd
*  but no larger than rendai
        rendld = min (rstart + dlogd, rendai)
        if (.not. logdfl) rendld = rstart
      end if

cger (next 2 lines)
      write (9,'(/" ** J =",i5," JLPAR =",i2," STARTED")') jtot,jlpar
cstart unix-ibm unix-aix unix-darwin
*      call flush_(9)
      call flush(9)
cend
cstart unix-dec unix-convex unix-iris
c;      call flush(9)
cend
* do not reduce maximum size of matrices if bound state calculation
      if (boundc) then
        ntop=nmax
      else
        ntop=nchtop
      endif
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
      t11=second()
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
cstart unix-dec unix-iris unix-hp
c;        call flush6
cend
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
cstart unix-dec unix-iris unix-hp
c;        call flush6
cend
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
          nfile = 14 + ien
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
      if (wavefl .and. .not. boundc)  call dclos(2)
      if (writs) then
        do 410 ien = 1, nerg
        nfile = 44 + ien
        call closf (nfile)
410     continue
      end if
      if (nerg .gt. 1) then
        if (airyfl) close (10)
        close (11)
      end if
420   call dater (cdate)
      if (.not. optifl) then
         write (6, 350)
         write (6,500) nchmax, timew, time, cdate
         write (9,500) nchmax, timew, time, cdate
500      format(' **** END OF CALCULATION ****',/,
     :    '      MAXIMUM NUMBER OF CHANNELS USED WAS:  ',i3,/,
     :    '      TIMING:  ELAPSED',(a),'/ CPU',(a),/,
     :    '      CURRENT DATE:  ',(a))
         write (6, 350)
         write (9, 350)
cstart unix-dec unix-iris unix-hp
c;         call flush6
cend
      end if
      close (23)
      endfile (9)
      close (9)
      goto 1
      end
