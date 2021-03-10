*************************************************************************
*                                                                       *
*                 system dependent routines library                     *
*                                                                       *
*************************************************************************
*                          routines included:                           *
*  0.  baschk      check that basistype is allowed
*  1.  sysdat      dispatcher to select specific sysdat routine         *
*  2.  syssav      dispatcher to select specific syssav routine         *
*  3.  ptread      dispatcher to select specific ptread routine         *
*  4.  sy1sg (sav1sg/ptr1sg) defines, saves variables and reads         *
*                  potential for singlet sigma scattering               *
*  5.  sy2sg (sav2sg/ptr2sg) defines, saves variables and reads         *
*                  potential for doublet sigma scattering               *
*  6.  sy2pi (sav2pi/ptr2pi) defines, saves variables and reads         *
*                  potential for doublet pi scattering                  *
*  7.  sysgpi (savsp/ptrsp) defines, saves variables and reads          *
*                  potential for singlet/pi sigma scattering            *
*  8.  sypi (savpi/ptrpi) defines, save variables and reads             *
*                  potential for general pi scattering                  *
*  9.  systp (savstp/ptrstp) defines, save variables and reads          *
*                  potential for symmetric top/atom scattering          *
*  10.  sy13p (sav13p/ptr13p) defines, save variables and reads          *
*                  potential for 1S / 3P atom scattering          *
*  11.  sy2mol (sav2mol/ptr2mol) defines, save variables and reads          *
*                  potential for 2 singlet sigma molecule scattering
*  12.  systpln (savstpln/ptrstpln) defines, save variables and reads
*                  potential for symmetric top and singlet sigma molecule *
*                  scattering
*  13.  sy22p (sav22p/ptr22p) defines, save variables and reads         *
*                  potential for 2S / 2P atom scattering          *
*  14.  sy1del (sav1del/ptr1del) defines, saves variables and reads         *
*                  potential for singlet delta scattering                  *
*  15.  syh2p (savh2p/ptrh2p) defines, saves variables and reads         *
*                  potential for homonuclear+2P atom scattering
*  16.  syh3p (savh3p/ptrh3p) defines, saves variables and reads         *
*                  potential for homonuclear+3P atom scattering
*  17.  sy2de (sav2de/ptr2de) defines, saves variables and reads         *
*                  potential for doublet-delta scattering

**************************************************************************
* -----------------------------------------------------------------------
      subroutine baschk(ival)
* -------------------------------------------------------
* subroutine to check that allowed basis is being called for
* this should be updated as new bases are added
* allowed basis types currently are:
*  1:  singlet sigma + atom
*  2:  doublet sigma + atom
*  3:  doublet pi + atom
*  4:  sigma | pi + atom
*  5:  general pi + atom
*  6:  symmetric top + atom
*  7:  1/3 P atom + atom
*  8:  1sigma + 1sigma
*  9:  symetric top + 1sigma
*  10: 2S atom + 2P atom
*  11: singlet delta + atom
*  12: homonuclear diatomic+2P atom
*  13: homonuclear diatomic+3P atom
*  14: double delta + atom
*  99: user defined basis
* author:  millard alexander
* current revision date:  13-oct-1999
*  variable in common block /comxbs/
*     maxbas    maximum number of allowed basis types
*     THIS SHOULD BE INCREASED AS BASIS ROUTINES ARE ADDED!!
*     COMMON BLOCK COMXBS DEFINED IN HIMAIN.F
* -------------------------------------------------------
      common /comxbs/ maxbas
      common /coselb/ ibasty
      logical icheck
      icheck=.false.
      do 100 i=1, maxbas
        if (ival .eq. i) icheck=.true.
100   continue
      if (ival .eq. 99) icheck= .true.
      if (.not. icheck) then
          write(6,50) ival
50        format(' *** IBASTY =',i3,' NOT YET INCLUDED,',
     :           ' BASIS ROUTINES AVAILABLE:'/,
     :         '   1 ->  SINGLET SIGMA + ATOM',/,
     :         '   2 ->  DOUBLET SIGMA + ATOM',/,
     :         '   3 ->  DOUBLET PI + ATOM',/,
     :         '   4 ->  SIGMA | PI + ATOM',/,
     :         '   5 ->  GENERAL PI + ATOM',/,
     :         '   6 ->  SYMMETRIC TOP + ATOM   ',/,
     :         '   7 ->  1/3 P ATOM + ATOM   ',/,
     :         '   8 ->  TWO 1SIGMA MOLECULES   ',/,
     :         '   9 ->  SYM TOP + 1SIGMA MOLECULE   ',/,
     :         '   10 ->  2P ATOM + 2S ATOM   ',/,
     :         '   11 ->  SINGLET DELTA + ATOM',/,
     :         '   12 ->  HOMONUCLEAR + 2P ATOM',/,
     :         '   13 ->  HOMONUCLEAR + 3P ATOM',/,
     :         '   14 ->  DOUBLET DELTA	+ ATOM',/,
     :         '   99 -> USER DEFINED BASIS  ',/)
      endif
      return
      end
* -----------------------------------------------------------------------
      subroutine sysdat (irpot, readpt, iread)
*   dispatcher to select correct sysdat routine
*   the correct routine is selected according to value of ibasty
*   the following sysdat routines are currently available:
*
*  ibasty         sysdat routine          kind of problem
*    1              sy1sg             singlet sigma scattering
*    2              sy2sg             doublet sigma scattering
*    3              sy2pi             doublet pi scattering
*    4              sysgpi            sigma/pi scattering
*    5              sypi              general pi scattering
*    6              systp             symmetric top scattering
*    7              sy13p             1/3 P atom scattering
*    8              sy2mol            2 1sigma molecules
*    9              systpln           symetric top + 1 sigma molecule
*    10             sy22p             2/2 P atom scattering
*    11             sy1del            singlet delta atom scattering
*    12             syh2p             homonuclear + 2p atom
*    13             syh3p             homonuclear + 3p atom
*    14             syh2del           doublet delta + atom
*    99             syusr             user supplied routine
*
*
*  if iread = 1 read data from input file
*  if iread = 0 return after defining variable names
*  author: b. follmeg
*  current revision date: 13-oct-1999
*  -----------------------------------------------------------------------
      integer ibasty, irpot, iread
      logical readpt
      common /coselb/ ibasty
      include "common/parbas"
* set default for vibrational quantum numbers to zero for each term
      do 10 it=1,maxtrm
      ivrow(1,it)=0
      ivcol(1,it)=0
10    ntv(it)=1
      if (ibasty .eq. 99) then
*  user supplied routine
        call syusr(irpot, readpt, iread)
        return
      endif
      goto (100,200,300,400,500,600,700,800,900,1000,1100,1200,
     :      1300,1400)
     :     ibasty
*  singlet sigma variables
100   call sy1sg(irpot, readpt, iread)
      return
*  doublet sigma variables
200   call sy2sg(irpot, readpt, iread)
      return
*  doublet pi variables
300   call sy2pi(irpot, readpt, iread)
      return
*  sigma/pi variables
400   call sysgpi(irpot, readpt, iread)
      return
*  general pi variables
500   call sypi(irpot, readpt, iread)
      return
*  symmetric top variables
600   call systp(irpot, readpt, iread)
      return
*  1/3 P atom variables
700   call sy13p(irpot, readpt, iread)
      return
* two 1 sigma molecules
800   call sy2mol(irpot, readpt, iread)
      return
* symmetric top + 1 sigma molecule
900   call systpln(irpot, readpt, iread)
      return
* 2S atom + 2P atom
1000  call sy22p(irpot, readpt, iread)
      return
* singlet delta variables
1100  call sy1del(irpot, readpt, iread)
      return
* homonuclear + 2P atom variables
1200  call syh2p(irpot, readpt, iread)
      return
* homonuclear + 3P atom variables
1300  call syh3p(irpot, readpt, iread)
      return
* doublet delta variable
1400  call sy2del(irpot, readpt, iread)
      return

      end
* -----------------------------------------------------------------------
      subroutine syssav (readpt)
*   dispatcher to select correct syssav routine
*   the correct routine is selected according to value of ibasty
*   the following savdat routines are currently available:
*
*  ibasty         syssav routine          kind of problem
*    1              sav1sg            singlet sigma scattering
*    2              sav2sg            doublet sigma scattering
*    3              sav2pi            doublet pi scattering
*    4              savsp             sigma/pi scattering
*    5              savpi             general pi scattering
*    6              savstp            symmetric top scattering
*    7              sav13p            1/3 P atom scattering
*    8              sav2mol           1sigma+1sigma
*    9              savstpln           symetric top + 1 sigma molecule
*    10             sav22p            2/2 P atom scattering *
*    11             sav1del           singlet delta scattering
*    12             savh2p            homonuclear + 2p atom
*    13             savh3p            homonuclear + 3p atom
*    14             sav2del           doublet delta + atom
*    99             savusr            user supplied routine *
*  author: b. follmeg
*  current revision date: 13-oct-1999
*  -----------------------------------------------------------------------
      integer ibasty
      logical readpt
      common /coselb/ ibasty
      dimension readpt(1)
      if (ibasty .eq. 99) then
*  user supplied routine
         call savusr(readpt)
         return
      endif
      goto (100,200,300,400,500,600,700,800,900,1000,1100,1200,
     :      1300,1400)
     :     ibasty
*  singlet sigma variables
100   call sav1sg(readpt)
      return
*  doublet sigma variables
200   call sav2sg(readpt)
      return
*  doublet pi variables
300   call sav2pi(readpt)
      return
*  sigma/pi variables
400   call savsp(readpt)
      return
*  general pi variables
500   call savpi(readpt)
      return
*  symmetric top variables
600   call savstp(readpt)
      return
*  1/3 P atom variables
700   call sav13p(readpt)
      return
*  1sigma+1sigma variables
800   call sav2mol(readpt)
      return
* symmetric top + 1 sigma molecule
900   call savstpln(irpot, readpt, iread)
      return
*  2/2 P atom variables
1000  call sav22p(readpt)
      return
*  singlet delta variables
1100  call sav1del(readpt)
      return
* homonuclear + 2P atom variables
1200  call savh2p(irpot, readpt, iread)
      return
* homonuclear + 3P atom variables
1300  call savh3p(irpot, readpt, iread)
      return
* doublet-delta + atom variables
1400  call sav2del(irpot, readpt, iread)
      return
      end
* -----------------------------------------------------------------------
      subroutine ptread (filnam,readpt)
*   dispatcher to select correct ptread routine
*   the correct routine is selected according to value of ibasty
*   the following ptread routines are currently available:
*
*  ibasty         ptread routine          kind of problem
*    1              ptr1sg            singlet sigma scattering
*    2              ptr2sg            doublet sigma scattering
*    3              ptr2pi            doublet pi scattering
*    4              ptrsp             sigma/pi scattering
*    5              ptrpi             general pi scattering
*    6              ptrstp            symmetric top scattering
*    7              ptr13p            1/3 P atom scattering
*    8              ptr2mol           1sigma+1sigma
*    9              ptrstpln          symetric top + 1 sigma molecule
*    10             ptr22p            2/2 P atom scattering
*    11             ptr1del           singlet delta scattering
*    12             pth2p             homonuclear + 2p atom
*    13             pth3p             homonuclear + 3p atom
*    14             pt2del            doublet delta + atom
*    99             ptrusr            user supplied routine
*
*  author: b. follmeg
*  current revision date: 13-oct-1999
*  -----------------------------------------------------------------------
      integer ibasty
      logical readpt
      character*(*) filnam
      common /coselb/ ibasty
      if (ibasty .eq. 99) then
*  user supplied routine
         call ptrusr(filnam,readpt)
         return
      endif
      goto (100,200,300,400,500,600,700,800,900,1000,1100,1200,
     :      1300,1400)
     :     ibasty
*  singlet sigma potential
100   call ptr1sg(filnam,readpt)
      return
*  doublet sigma potential
200   call ptr2sg(filnam,readpt)
      return
*  doublet pi potential
300   call ptr2pi(filnam,readpt)
      return
*  sigma/pi potential
400   call ptrsp(filnam,readpt)
      return
*  general pi variables
500   call ptrpi(filnam,readpt)
      return
*  symmetric top variables
600   call ptrstp(filnam,readpt)
      return
*  1/3 P atom variables
700   call ptr13p(filnam,readpt)
      return
*  1sigma+1sigma variables
800   call ptr2mol(filnam,readpt)
      return
* symmetric top + 1 sigma molecule
900   call ptrstpln(irpot, readpt, iread)
      return
*  2/2 P atom variables
1000  call ptr22p(filnam,readpt)
      return
*  singlet delta variables
1100  call ptr1del(filnam,readpt)
      return
* homonuclear + 2P atom variables
1200  call ptrh2p(irpot, readpt, iread)
      return
* homonuclear + 3P atom variables
1300  call ptrh3p(irpot, readpt, iread)
      return
* doublet delta + atom variables
1400  call ptr2del(irpot, readpt, iread)
      return
      end
*  -----------------------------------------------------------------------
      subroutine sy1sg (irpot, readp, iread)
*  subroutine to read in system dependent parameters for singlet-sigma
*   + atom scattering using werner-follmeg potential form
*  if iread = 1 read data from input file
*  if iread = 0 return after defining variable names
*  current revision date: 10-mar-1994 by mha
*  -----------------------------------------------------------------------
*  variables in common /cosysr/
*    isrcod:  number of real parameters
*    brot, drot, hrot
*  variable in common /cosysi/
*    nscode:  total number of system dependent parameters
*             nscode = isicod + isrcod +3
*    isicod:  number of integer system dependent parameters
*    nterm:   number of different angular terms in potential
*              nb for singlet sigma molecule, nterm should be 1
*    jmin, jmax (see below)
*  variable in common /cosys/
*    scod:    character*8 array of dimension lcode, which contains names
*             of all system dependent parameters
*  line 13:
*    nlam:     total number of anistropic terms in potential
*    jmin:     minimum molecular rotational angular momenta included in
*              the channel basis
*    jmax      maximum molecular rotational angular momenta inncluded in
*              the channel basis
*  line 14:
*    brot:    rotational constant of the molecule in cm-1
*  line 16:
*    filnam:  name of file containing potential parameters
*
*  subroutines called: loapot(iunit,filnam)
*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer irpot
      logical readp, existf
      logical airyfl, airypr, bastst, batch, chlist, csflag,
     :                flaghf, flagsu, ihomo,lpar
      character*1 dot
      character*4 char
      character*8 scod
      character*(*) fname
      character*40 filnam, line, potfil
      include "common/parsys"
      common /cosys/ scod(lencod)
      common /cosysi/ nscode,isicod,iscod(maxpar)
      common /cosysr/ isrcod, junkr, rcod(maxpar)
      common /conlam/ nlam
      include "common/parbas"
      common/covib/ nvib,ivib(maxvib)
      common /coskip/ nskip,iskip
      common /colpar/ airyfl, airypr, bastst, batch, chlist, csflag,
     :                flaghf, flagsu, ihomo,lpar(18)

      common /coselb/ ibasty
      save potfil
      equivalence(iscod(1),nterm),(iscod(2),nvibmn),(iscod(3),nvibmx)
      include "common/comdot"
*     number and names of system dependent parameters
*  set default values for singlet-sigma scattering
      nterm = 1
      if (iread .eq. 0) then
        nvibmn = 0
        nvibmx = 0
        lammin(1)=0
        lammax(1)=-1
        mproj(1)=0
        niout=1
        indout(1)=0
        nvib = 1
      endif
      if (.not.readp)irpot=1
      if (iread .eq. 1) irpot=1
      if (ihomo) nskip = 2
      potfil = ' '
*  read number of vib states
      if(iread.ne.0) read (8, *, err=88) nvib, nvibmn,nvibmx
      if(nvib.gt.nvibmx-nvibmn+1) then
        write (6,40) nvib, nvibmn, nvibmx
40      format(' ** PROBABLE VIBRATIONAL LEVEL NUMBERING ERROR:',/,
     :         '   VMIN =',i2,', VMAX=',i3,', BUT NVIB =',i3)
        return
      endif
      if(nvib.gt.maxvib) stop 'nvib'
      if(nvib.le.0) nvib=1
*  read data for each vib state
*  brot, drot, hrot are bv, dv, hv (see huber herzberg, page x)
      scod(1)='NTERM'
      scod(2)='VMIN'
      scod(3)='VMAX'
      isrcod=0
      isicod=3
      iofr=2*nvib+4-1
      do  71   i = 1,nvib
      if(iread.ne.0) then
        read (8, *, err=99) ivib(i),(iscod(isicod+j),j=1,2)
        read (8, *, err=99) (rcod(isrcod+j),j=1,3)
        read (8, *, err=99) rcod(isrcod+4)
        end if
        char=' '
        if(nvib.gt.1.or.ivib(i).ne.0) then
        if(ivib(i).le.9) write(char,'(''('',i1,'')'')') ivib(i)
        if(ivib(i).gt.9) write(char,'(''('',i2,'')'')') ivib(i)
        end if
      scod(isicod+1)='JMIN'//char
      scod(isicod+2)='JMAX'//char
      scod(iofr+1)='BROT'//char
      scod(iofr+2)='DROT'//char
      scod(iofr+3)='HROT'//char
      scod(iofr+4)='EVIB'//char
      iofr=iofr+4
      isicod=isicod+2
71    isrcod=isrcod+4
      if(isicod+isrcod+3.gt.lencod) stop 'lencod'
      nscode=isicod+isrcod
      line=' '
      if(.not.readp.or.iread.eq.0) then
        call loapot(1,' ')
        close (8)
        return
      endif
      read (8, 85, end=186) line
85    format (a)
      goto 186
* here if read error occurs
88    write(6,89)
89    format(/' *** ERROR DURING READ FROM INPUT FILE ***')
      return
99    write (6, 100) i
100   format(/' ** ERROR DURING READ:',
     :  ' PROBABLY NOT ENOUGH VIBRATIONAL LEVELS SUPPLIED')
      return
* --------------------------------------------------------------
      entry ptr1sg (fname,readp)
      line = fname
      readp = .true.
186   if (readp) then
        l=1
        call parse(line,l,filnam,lc)
        if(lc.eq.0) then
          write(6,190)
190       format(' FILENAME MISSING FOR POTENTIAL INPUT')
        end if
        j=index(filnam(1:lc),dot)
        if(j.eq.0) then
          call gennam(potfil,filnam,0,'BIN',lc)
          filnam = potfil
        end if
        potfil=filnam
        inquire(file=filnam,exist=existf)
        if(.not.existf) then
          write(6,195) filnam(1:lc)
195       format(' FILE ',(a),' NOT FOUND')
          return
        end if
* now call loapot(iunit,filnam) routine to read potential parameters
        call loapot(1,filnam)
      end if
*      close (8)
      irpot=1
      return
* --------------------------------------------------------------
      entry sav1sg (readp)
*  save input parameters for singlet-sigma + atom scattering
      if (nvibmx .lt. nvibmn) then
        write (6, 210) nvibmx, nvibmn
210     format ('**  VMAX =',i3,' .LT. VMIN =',i3,' SET VMAX = VMIN')
        nvibmx=nvibmn
      endif
      write (8, 220) nvib, nvibmn,nvibmx
220   format(3i4, t34,'nvib, vmin,vmax')
      iofi=3
      iofr=0
      do 301 i=1,nvib
      write (8, 310) ivib(i),(iscod(iofi+j),j=1,2)
310   format (3i4, t50,'iv,jmin,jmax')
      write (8, 320) (rcod(iofr+j),j=1,4)
      iofi=iofi+2
      iofr=iofr+4
320   format(3g14.6,t50,'brot,drot,hrot'/f15.8,t50,'evib')
301   continue
      write (8, 85) potfil
      return
      end
*  -----------------------------------------------------------------------
      subroutine sy2sg (irpot, readpt, iread)
      implicit double precision (a-h,o-z)
*  subroutine to read in system dependent parameters for doublet-sigma
*   + atom scattering
*  if iread = 1 read data from input file
*  if iread = 0 return after defining variable names
*  current revision date: 8-apr-1997 by mha
*  -----------------------------------------------------------------------
*  variables in common bloc /cosysr/
*    isrcod:   total number of real system dependent variables
*    brot:   rotational constant for 2sigma state in cm-1
*    gsr:      2sigma state spin-rotation constant in cm-1
*  variables in common block /cosysi/
*    nscode:   total number of system dependent variables
*    isicod:   total number of integer system dependent variables
*    nterm:    number of different diabatic potentials, this is 1 here
*    nrmax:    the maximum case (b) rotational angular momenta
*    npar:     number of symmetry doublets included (npar=2 will ensure
*              both spin doublets)
*    isym:     if isym=+1, then the electronic symmetry is sigma-plus
*              if isym=-1, then the electronic symmetry is sigma-minus
*    igu:      if igu=+1, then the inversion symmetry is gerade
*              if igu=-1, then the inversion symmetry is ungerade
*    isa:      s/a symmetry index, if the molecule is homonuclear (ihomo=t)
*              then, if isa=+1 then only the s-levels are included in the
*              basis, if isa=-1, then only the a-levels are included
*    interp:   parameter to control spline interpolation in subroutine pot;
*              this should always be set equal to 1 for the first calculation
*  variable in common /cosys/
*    scod:    character*8 array of dimension nscode, which contains names
*             of all system dependent parameters.  note that the ordering
*             of the variable names in scod must correspond to the ordering
*             of the variable names in cosysi followed by the ordering of
*             variable names in cosysr followed by lammin, lammax, and mproj
* ------------------------------------------------------------------------
      logical readpt,existf
      character*8 scod
      character*(*) fname
      character*40 line,filnam,potfil
      include "common/parsys"
      include "common/parbas"
      common /cosys/ scod(lencod)
      common /cosysi/ nscode, isicod, nterm, nrmax, npar, isym, igu,
     :                isa
      common /cosysr/ isrcod, junkr, brot, gsr, drot, hrot
      common /conlam/ nlam
      logical         airyfl, airypr, bastst, batch, chlist, csflag,
     :                flaghf, flagsu, ihomo,lpar
      common /colpar/ airyfl, airypr, bastst, batch, chlist, csflag,
     :                flaghf, flagsu, ihomo,lpar(18)
      common /coskip/ nskip,iskip
      character*1 dot
      save potfil
      include "common/comdot"
      irpot = 1
*     number and names of system dependent parameters
      isicod = 6
      isrcod = 4
      nscode = isicod+isrcod
*  first all the system dependent integer variables
*  in the same order as in the common block /cosysi/
*  variable names must be in uppercase of maximum length 6 characters
      scod(1) = 'NTERM'
      scod(2) = 'NRMAX'
      scod(3) = 'NPAR'
      scod(4) = 'ISYM'
      scod(5) = 'IGU'
      scod(6) = 'ISA'
*  then all the system dependent real variables
*  in the same order as in the common block /cosysr/
      scod(7) = 'BROT'
      scod(8) = 'GSR'
      scod(9) = 'DROT'
      scod(10) = 'HROT'
*  set default values for doublet-sigma scattering
      potfil=' '
      nterm = 1
      if (iread .eq. 0) then
        nrmax = 3
        mproj(1) = 0
        lammin(1)= 1
        lammax(1) = -1
        isa=0
        isym=1
        npar=2
        niout=2
        indout(1)=-1
        indout(2)=1
      endif
      if(iread.eq.0) return
*  line 13
      read (8, *,err=888) nrmax, npar, isym, igu, isa
*  line 14
      read (8, *,err=888) brot, gsr, drot, hrot
*  line 16 name of file containing potential parameters
      line=' '
      if(.not.readpt.or.iread.eq.0) then
        call loapot(1,' ')
        return
      endif
      read (8, 285, end=286) line
      potfil=line
285   format (a)
      goto 286
* here if read error occurs
888   write(6,1000)
1000  format(/'   *** ERROR DURING READ FROM INPUT FILE ***')
      return
      entry ptr2sg (fname,readpt)
      line = fname
      readpt = .true.
286   if (readpt) then
        l=1
        call parse(line,l,filnam,lc)
        if(lc.eq.0) then
          write(6,1020)
1020      format(' FILENAME MISSING FOR POTENTIAL INPUT')
        end if
        j=index(filnam(1:lc),dot)
        if(j.eq.0) then
           call gennam(potfil,filnam,0,'BIN',lc)
           filnam = potfil
        end if
        potfil=filnam
        inquire(file=filnam,exist=existf)
        if(.not.existf) then
          write(6,1025) filnam(1:lc)
1025      format(' FILE ',(a),' NOT FOUND')
          return
        end if
* now call loapot(iunit,filnam) routine to read potential parameters
        call loapot(1,filnam)
      end if
      close (8)
      if(nskip.eq.2) ihomo=.true.
      nlam=0
      nlam=nlam+(lammax(1)-lammin(1))/nskip+1
      irpot=1
      return
      entry sav2sg (readpt)
*  save input parameters for doublet-sigma + atom scattering
*  line 13:
      write (8, 220) nrmax, npar, isym, igu, isa
220   format (5i4, t50, 'nrmax, npar, isym, igu, isa')
      write (8, 320) brot, gsr, drot, hrot
320   format (f12.6,3g12.4,t50,'brot, gsr, drot, hrot')
      write (8,330) potfil
330   format(a)
      return
      end
* ----------------------------------------------------------------
      subroutine sy2pi (irpot, readpt, iread)
*  subroutine to read in system dependent parameters for doublet-pi
*   + atom scattering using werner-follmeg potential form
*  if iread = 1 read data from input file
*  if iread = 0 return after defining variable names
*  current revision date: 19-aug-1991 by mha
*  -----------------------------------------------------------------------
*    nlam:      the total number of angular coupling terms
*  variables in common bloc /cosysr/
*    isrcod:   total number of real system dependent variables
*    brot:     rotational constant in cm-1
*    aso:      spin-orbit constant in cm-1
*    p, q:     lambda-doubling constants cm-1
*  variables in common block /cosysi/
*    nscode:  total number of system dependent parameters
*             nscode = isicod + isrcod + 3
*    isicod:   total number of integer system dependent variables
*    nterm:   number of different angular terms in potential
*              nb for singlet sigma molecule, nterm should be 1
*    jmax:     the maximum rotational angular momenta for each channel
*              in each spin-orbit manifold with convention
*              omega .le. j .le. jmax+0.5
*    igu:      permutation inversion symmetry of electronic state
*              igu=1 for gerade states, igu=-1 for ungerade states
*              for heteronuclear molecules igu should be +1
*    isa:      s/a label for molecular states. if ihomo=.true. then only
*              s states will be included if isa=1 and only a states if
*              isa=-1
*    npar:     number of symmetry doublets included (npar=2 will ensure
*              both lambda doublets)
*  variable in common /cosys/
*    scod:    character*8 array of dimension nscode, which contains names
*             of all system dependent parameters
*  variable in common block /coskip/
*   nskip  for a homonuclear molecule lamda is running in steps of nskip=2
*          for a heteronuclear molecule nskip=1
*
*   skip   same as nskip, used for consistency check
*
*  subroutines called: loapot(iunit,filnam)
      implicit double precision (a-h,o-z)
      include "common/parsys"
      logical readpt, existf
      integer irpot
      character*1 dot
      character*8 scod
      character*(*) fname
      character*40 filnam, line, potfil
      include "common/parbas"
      common /cosys/ scod(lencod)
      common /cosysi/ nscode, isicod, nterm, jmax, igu, isa, npar
      common /cosysr/ isrcod, junkr, brot, aso, p, q
      logical         airyfl, airypr, bastst, batch, chlist, csflag,
     :                flaghf, flagsu, ihomo,lpar
      common /colpar/ airyfl, airypr, bastst, batch, chlist, csflag,
     :                flaghf, flagsu, ihomo,lpar(18)
      common /conlam/ nlam
      common /coskip/ nskip,iskip
      save potfil
      include "common/comdot"
      isicod = 5
      isrcod = 4
      nscode = isicod+isrcod
      scod(1) = 'NTERM'
      scod(2) = 'JMAX'
      scod(3) = 'IGU'
      scod(4) = 'ISA'
      scod(5) = 'NPAR'
      scod(6) = 'BROT'
      scod(7) = 'ASO'
      scod(8) = 'P'
      scod(9) = 'Q'
      potfil = ' '
*  set default values for doublet pi scattering
      nterm = 2
      mproj(2) = 2
      if (iread .eq. 0) then
        mproj(1) = 0
        lammin(1)= 1
        lammax(1) = 1
        lammin(2) = 2
        lammax(2) = 2
        niout=4
        indout(1)=-1
        indout(2)=1
        indout(3)=-2
        indout(4)=2
        npar=2
        igu=1
        isa=0
        jmax=3
      endif
      if(iread.eq.0) return
*  line 13
      read (8, *, err=888) jmax, igu, isa, npar
*  line 14
      read (8, *, err=888) brot, aso, p, q
*  line 15 name of file containing potential parameters
      line=' '
      if(.not.readpt.or.iread.eq.0) then
        call loapot(1,' ')
        return
      endif
      read (8, 285, end=286) line
285   format (a)
      goto 286
* here if read error occurs
888   write(6,1000)
1000  format(/'   *** ERROR DURING READ FROM INPUT FILE ***')
      return
      entry ptr2pi (fname,readpt)
      line = fname
      readpt = .true.
286   if (readpt) then
        l=1
        call parse(line,l,filnam,lc)
        if(lc.eq.0) then
          write(6,1020)
1020      format(' FILENAME MISSING FOR POTENTIAL INPUT')
        end if
        j=index(filnam(1:lc),dot)
        if(j.eq.0) then
          call gennam(potfil,filnam,0,'BIN',lc)
          filnam = potfil
        end if
        potfil=filnam
        inquire(file=filnam,exist=existf)
        if(.not.existf) then
          write(6,1025) filnam(1:lc)
1025      format(' FILE ',(a),' NOT FOUND')
          return
        end if
* now call loapot(iunit,filnam) routine to read potential parameters
        call loapot(1,filnam)
        end if
      close (8)
      irpot=1
      return
      entry sav2pi (readpt)
*  save input parameters for doublet-pi + atom scattering
*  line 13:
      write (8, 315) jmax, igu, isa, npar
315   format(4i4,18x,'jmax, igu, isa, npar')
*  line 14
      write (8, 320) brot, aso, p, q
320   format(f10.5,f10.4,2(1pg11.3),' brot, aso, p, q')
*  line 15
      write (8, 285) potfil
      return
      end
* ----------------------------------------------------------------------
      subroutine sysgpi (irpot, readpt, iread)
* ----------------------------------------------------------------------
*
*  subroutine to read in system dependent parameters for 2pi-2sig
*   + atom scattering using werner-follmeg potential form
*  current revision date: 21-dec-1995 by pjd, mha, ade, ab
*
*  variables in common/cobspt/ must be set in loapot!!
*
      implicit double precision (a-h,o-z)
      include "common/parsys"
      logical readpt, existf
      character*8 scod,rcod
      character*8 char
      character*(*) fname
      character*1 dot
      character*40 filnam, line, potfil
      include "common/parbas"
      common /covib/ nvibs,ivibs(maxvib),nvibp,ivibp(maxvib)
      common /cosys/ scod(lencod)
      common /cosyr/ rcod(maxpar)
      common /cosysi/ nscode, isicod, ispar(maxpar)
      common /cosysr/ isrcod, junkr, rspar(maxpar)
      logical         airyfl, airypr, bastst, batch, chlist, csflag,
     :                flaghf, flagsu, ihomo,lpar
      common /colpar/ airyfl, airypr, bastst, batch, chlist, csflag,
     :                flaghf, flagsu, ihomo,lpar(18)
      common /coskip/ nskip,iskip
      common /conlam/ nlam
      save potfil
      include "common/comdot"
*  set default values for 2pi-2sigma scattering
      do 5 i=1,maxpar
      rpar=0
5     ispar(i)=0
      potfil = ' '
      isicod=0
      isrcod=0
      ispar(1)=4
      if (iread .eq. 0) then
        niout=4
        nterm=4
        indout(1)=-1
        indout(2)=1
        indout(3)=-2
        indout(4)=2
        ispar(2)=1
        ispar(3)=1
        ispar(4)=1
        ispar(5)=0
        ispar(6)=1
        ispar(7)=2
        ispar(8)=1
        ispar(12)=3
        ispar(13)=1
        ispar(14)=2
        ispar(18)=3
        lammin(1)=0
        lammin(2)=0
        lammin(3)=1
        lammin(4)=2
        lammax(1)=-1
        lammax(2)=-1
        lammax(3)=-1
        lammax(4)=-1
        mproj(1)=0
        mproj(2)=0
        mproj(3)=1
        mproj(4)=2
        nvibs=1
        nvibp=1
        nskip=1
      endif
*  read parameters which determine type
      scod(1)='NTERM'
      scod(2)='ISG'
      scod(3)='IPI'
      scod(4)='ISGPI'
      scod(5)='ISA'
      if(iread.ne.0) read (8, *, err=800) (ispar(isicod+j),j=2,5)
      nterm=ispar(1)
      isg=ispar(2)
      ipi=ispar(3)
          isgpi=ispar(4)
      if (isg .eq. 1 .and. ipi .eq. 0) then
        nterm=1
        if (isgpi .ne. 0) then
          write (6, 6) isgpi
6         format ('  ISGPI SET EQUAL TO 0 BECAUSE ONLY SIGMA STATE')
          isgpi=0
          ispar(4)=0
        endif
      endif
      if (isg .eq. 0 .and. ipi .eq. 1) then
        nterm=2
        if (isgpi .ne. 0) then
          write (6, 10) isgpi
10        format ('  ISGPI SET EQUAL TO 0 BECAUSE ONLY PI STATE')
          isgpi=0
          ispar(4)=0
        endif
      endif
      isicod=5
      if(isg.eq.0.and.ipi.eq.0) then
        isg=1
        ipi=1
      end if
      ispar(1)=nterm
      if(isg.ne.0) then
*  read symmetry parameters for sigma state
        if(iread.ne.0) read (8, *, err=800) (ispar(isicod+j),j=1,3)
        scod(isicod+1)='IGUSG'
        scod(isicod+2)='NPARSG'
        scod(isicod+3)='ISYMSG'
        isicod=isicod+3
*  read number of vib states for sigma state
        if(iread.ne.0) read (8, *, err=800) nvibs,
     :      (ispar(isicod+j),j=1,2)
        if(nvibs.gt.maxvib) stop 'nvibs'
        nvmins=ispar(isicod+1)
        nvmaxs=ispar(isicod+2)
        scod(isicod+1)='VMINSG'
        scod(isicod+2)='VMAXSG'
        isicod=isicod+2
        if(nvibs.le.0) nvibs=1
*  read data for each vib state
*  brot, drot, hrot are b, d, h
*  gs is gamma, gsd is gamma(d), gsh is gamma(h)
*  evibsg is E (see Kotlar et al., JMS 80, 86 (1980), Table II))
        do 15 i=1,nvibs
              if(isicod+2.gt.maxpar) stop 'isicod'
              if(isrcod+7.gt.maxpar) stop 'isrcod'
              if(iread.ne.0) then
            read (8,*,err=800) ivs,(ispar(isicod+j),j=1,2)
            ivibs(i)=ivs
            read (8,*,err=800) (rspar(isrcod+j),j=1,3)
            read (8,*,err=800) (rspar(isrcod+j),j=4,6)
            read (8,*,err=800) rspar(isrcod+7)
              end if
              char=' '
              if(nvibs.gt.1.or.ivs.ne.0) then
            if(ivs.le.9) write(char,'(''('',i1,'')'')') ivs
            if(ivs.gt.9) write(char,'(''('',i2,'')'')') ivs
              end if
          scod(isicod+1)='NMIN'//char
          scod(isicod+2)='NMAX'//char
          rcod(isrcod+1)='BSG'//char
          rcod(isrcod+2)='DSG'//char
          rcod(isrcod+3)='HSG'//char
          rcod(isrcod+4)='GS'//char
          rcod(isrcod+5)='GSD'//char
          rcod(isrcod+6)='GSH'//char
          rcod(isrcod+7)='ESG'//char
          isicod=isicod+2
          isrcod=isrcod+7
15      continue
      end if
      if(ipi.ne.0) then
*  read symmetry parameters for pi state
        if(iread.ne.0) read (8, *, err=800) (ispar(isicod+j),j=1,2)
        scod(isicod+1)='IGUPI'
        scod(isicod+2)='NPARPI'
        isicod=isicod+2
*  read number of vib states for pi state
        if(iread.ne.0) read (8,*,err=800) nvibp,
     :           (ispar(isicod+i),i=1,2)
        if(nvibp.gt.maxvib) stop 'nvibp'
        nvminp=ispar(isicod+1)
        nvmaxp=ispar(isicod+2)
        scod(isicod+1)='VMINPI'
        scod(isicod+2)='VMAXPI'
        isicod=isicod+2
*  read data for each vib state
        if(nvibp.le.0) nvibp=1
        do 20 i=1,nvibp
          if(isicod+3.gt.maxpar) stop 'isicod'
          if(isrcod+13.gt.maxpar) stop 'isrcod'
          if(iread.ne.0) then
            read (8,*,err=800) ivp,(ispar(isicod+j),j=1,2)
            ivibp(i)=ivp
            read (8,*,err=800) (rspar(isrcod+j),j=1,3)
            read (8,*,err=800) (rspar(isrcod+j),j=4,6)
            read (8,*,err=800) (rspar(isrcod+j),j=7,9)
            read (8,*,err=800) (rspar(isrcod+j),j=10,12)
          end if
          char=' '
          if(nvibp.gt.1.or.ivp.gt.0) then
            if(ivp.le.9) write(char,'(''('',i1,'')'')') ivp
            if(ivp.gt.9) write(char,'(''('',i2,'')'')') ivp
          end if
          scod(isicod+1)='JMIN'//char
          scod(isicod+2)='JMAX'//char
          scod(isicod+3)='IPT'//char
          rcod(isrcod+1)='BPI'//char
          rcod(isrcod+2)='DPI'//char
          rcod(isrcod+3)='HPI'//char
          rcod(isrcod+4)='ASO'//char
          rcod(isrcod+5)='P'//char
          rcod(isrcod+6)='PD'//char
          rcod(isrcod+7)='Q'//char
          rcod(isrcod+8)='QD'//char
          rcod(isrcod+9)='QH'//char
          rcod(isrcod+10)='EPI'//char
          rcod(isrcod+11)='AD'//char
          rcod(isrcod+12)='GPI'//char
* check that ad and gpi are not both defined, this would not make any
* sense, since they correspond to different hamiltonians
          if (rspar(isrcod+11) .ne. 0.d0 .and.
     :        rspar(isrcod+12) .ne. 0.d0) then
            write (6, 16) isrcod+11, rspar(isrcod+11),
     :                    isrcod+12, rspar(isrcod+12)
16          format (' *** RSPAR(',i2,') =',1pg12.5,' AND RSPAR(',i2,
     :              ') =',1pg12.5,' BOTH NOT ZERO; ABORT ***')
            call abort
          endif
          if(isgpi.ne.0) then
            if(iread.ne.0) then
              read (8,*,err=800) ispar(isicod+3),(rspar(isrcod+j),
     >                           j=13,14)
              ivs=ispar(isicod+3)
            end if
            if(ivs.le.9.and.ivp.le.9)
     :        write(char,'(''('',i1,''/'',i1,'')'')') ivs,ivp
            if(ivs.gt.9.and.ivp.le.9)
     :        write(char,'(''('',i2,''/'',i1,'')'')') ivs,ivp
            if(ivs.le.9.and.ivp.gt.9)
     :        write(char,'(''('',i1,''/'',i2,'')'')') ivs,ivp
            if(ivs.gt.9.and.ivp.gt.9)
     :        write(char,'(''('',i2,''/'',i2,'')'')') ivs,ivp
            rcod(isrcod+13)='A'//char
            rcod(isrcod+14)='B'//char
          endif
          if (isgpi.eq.0 .and.ispar(isicod+3).ne.0) then
            write (6, 19)
19          format(' ONLY PI MANIFOLD; IPERT SET TO ZERO')
            ispar(isicod+3) = 0
          endif
          isicod=isicod+3
          if (isgpi .eq. 0) then
            isrcod=isrcod+12
          else if (isgpi .eq. 1) then
            isrcod=isrcod+14
          endif
20      continue
      end if
      if(isicod+isrcod.gt.lencod) stop 'lencod'
      do 30 i=1,isrcod
30    scod(isicod+i)=rcod(i)
      nscode = isicod+isrcod
*  read file name containing potential parameters
      if(.not.readpt.or.iread.eq.0) then
        call loapot(1,' ')
        return
      endif
      read (8, 50, end=1000) line
50    format (a)
      goto 1000
* here if read error occurs
800   write(6,900)
900   format(/'   *** ERROR DURING READ FROM INPUT FILE ***')
      return
*
      entry ptrsp (fname,readpt)
*
      line = fname
      readpt = .true.
1000  if (readpt) then
        l=1
        call parse(line,l,filnam,lc)
        if(lc.eq.0) then
          write(6,1020)
1020      format(' FILENAME MISSING FOR POTENTIAL INPUT')
        end if
        j=index(filnam(1:lc),dot)
        if(j.eq.0) then
          call gennam(potfil,filnam,0,'BIN',lc)
          filnam = potfil
        end if
        potfil=filnam
        inquire(file=filnam,exist=existf)
        if(.not.existf) then
          write(6,1025) filnam(1:lc)
1025      format(' FILE ',(a),' NOT FOUND')
          return
        end if
* now call loapot(iunit,filnam) routine to read potential parameters
        call loapot(1,filnam)
      end if
* check for consistency
      close (8)
      irpot=1
      ihomo=nskip.eq.2
      return
*
      entry savsp (readpt)
*
      nterm=ispar(1)
      isg=ispar(2)
      ipi=ispar(3)
      isicod=0
      isrcod=0
      write (8, 105) (ispar(isicod+j),j=2,5),'isg, ipi, isgpi, isa'
      isicod=isicod+5
      if(isg.ne.0) then
*  save symmetry parameters for sigma state
        write (8, 103) (ispar(isicod+j),j=1,3),
     :          'igusg, nparsg, isymsg'
        isicod=isicod+3
*  save number of vib states for sigma state
        write (8, 103) nvibs,(ispar(isicod+j),j=1,2)
     :                ,'nvibs, vminsg, vmaxsg'
        nvmins=ispar(isicod+1)
        nvmaxs=ispar(isicod+2)
        isicod=isicod+2
*  save data for each vib state
        do 100 i=1,nvibs
          write (8,103) ivibs(i),(ispar(isicod+j),j=1,2),
     :                  'ivs,  nmin,  nmax'
          write (8,203) (rspar(isrcod+j),j=1,3),'bsg, dsg, hsg'
          write (8,203) (rspar(isrcod+j),j=4,6),'gs, gsd, gsh'
          write (8,201)  rspar(isrcod+7),'esg'
          isicod=isicod+2
          isrcod=isrcod+7
100     continue
      end if
      if(ipi.ne.0) then
        write (8, 102) (ispar(isicod+j),j=1,2),'igupi, nparpi'
        isicod=isicod+2
        write (8,103) nvibp,(ispar(isicod+j),j=1,2),
     :                'nvibp, vminpi, vmaxpi'
        nvminp=ispar(isicod+1)
        nvmaxp=ispar(isicod+2)
        isicod=isicod+2
        do 101 i=1,nvibp
          write (8,103) ivibp(i),(ispar(isicod+j),j=1,2),
     :           'ivp, jmin, jmax'
          write (8,203) (rspar(isrcod+j),j=1,3),'bpi, dpi, hpi'
          write (8,203) (rspar(isrcod+j),j=4,6),'aso, p, pd'
          write (8,203) (rspar(isrcod+j),j=7,9),'q, qd, qd'
          write (8,202) (rspar(isrcod+j),j=10,12),'epi, adi, gpi'
          if(ispar(4).ne.0) then
            ivs=ispar(isicod+3)
            write (8,204) ivs,(rspar(isrcod+k),k=13,14),
     >                    'ipt, a, b'
            isrcod=isrcod+14
          else if(ispar(4).eq.0) then
            isrcod=isrcod+12
          end if
          isicod=isicod+3
101     continue
      end if
102   format(2i4,t50,a)
103   format(3i4,t50,a)
105   format(4i4,t50,a)
201   format(g14.8,t50,a)
202   format(g14.8,2g14.6,t50,a)
203   format(3g14.6,t50,a)
204   format(i4,2g14.6,t50,a)
      write (8, 300) potfil
300   format(a)
      return
      end
*  ---------------------------------------------------------------------
      subroutine sypi (irpot, readpt, iread)
*  subroutine to read in system dependent parameters for singlet,
*  doublet or triplet pi molecule + atom or flat surface scattering
*  if iread = 1 read data from input file
*  if iread = 0 return after defining variable names
*  author:  didier lemoine and millard alexander
*  current revision date:  4-mar-1996 by mha
*  ---------------------------------------------------------------------
*  variable in common /cosys/
*    scod:    character*8 array of dimension nscode, which contains
*             names of all system dependent parameters. note that the
*             ordering of the variable names in scod must correspond
*             to the ordering of the variable names in cosysi, cosysr,
*             cosysl and cobspt respectively
*  variables in common block /cosysi/
*    nscode:   total number of system dependent parameters
*              nscode = isicod + isrcod + islcod + 3
*    isicod:   total number of integer system dependent variables
*    nterm:    number of different angular terms in potential
*              NB for pi molecule, nterm should be 1 or 2
*    jmax:     the maximum rotational angular momenta (-1/2 for 2pi) for
*              each channel in each spin-orbit manifold with convention
*              omega .le. j .le. jmax(+1/2 for 2pi)
*    igu:      permutation inversion symmetry of electronic state
*              igu=1 for gerade states, igu=-1 for ungerade states
*              for heteronuclear molecules igu is not used
*    isa:      s/a label for molecular states. if ihomo=.true.,
*              then only s states will be included if isa=1 and
*              only a states if isa=-1
*    npar:     number of symmetry doublets included
*              (npar=2 will ensure both lambda doublets)
*    imult:    spin multiplicity of pi state (imult = 2*S+1)
*    nman:     number of spin-orbit manifolds included
*              for a 2pi state nman=1 selects the
*              omega=1/2 manifold in a case (a) basis
*              otherwise, nman is set = imult in subroutine bapi
*  variables in common bloc /cosysr/
*    isrcod:   total number of real system dependent variables
*    brot:     rotational constant in cm-1
*    aso:      spin-orbit constant in cm-1
*    o, p, q:  lambda-doubling constants cm-1
*    dmom:     dipole moment of molecule in Debye
*    efield:   stark field in kV / cm
*  variables in common bloc /cosysl/
*    islcod:   total number of logical system dependent variables
* ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical readpt, existf
      logical         airyfl, airypr, bastst, batch, chlist, csflag,
     :                flaghf, flagsu, ihomo,lpar
      include "common/parsys"
      character*8 scod
      character*(*) fname
      character*40 filnam, line, potfil
      character*1 dot
      common /colpar/ airyfl, airypr, bastst, batch, chlist, csflag,
     :                flaghf, flagsu, ihomo,lpar(18)
      include "common/parbas"
      common /cosys/ scod(lencod)
      common /cosysi/ nscode, isicod, nterm, jmax, igu, isa,
     :                npar, imult, nman
      common /cosysr/ isrcod, junkr, brot, aso, o, p, q, dmom, efield
      common /cosysl/ islcod
      common /conlam/ nlam
      common /coskip/ nskip,iskip
      save potfil
      include "common/comdot"
      potfil = ' '
*  number and names of system dependent parameters
*  first all the system dependent integer variables
*  in the same order as in the common block /cosysi/
*  variable names must be in uppercase of maximum length 6 characters
      isicod = 7
      scod(1) = 'NTERM'
      scod(2) = 'JMAX'
      scod(3) = 'IGU'
      scod(4) = 'ISA'
      scod(5) = 'NPAR'
      scod(6) = 'IMULT'
      scod(7) = 'NMAN'
*  then all the system dependent real variables
*  in the same order as in the common block /cosysr/
      isrcod = 7
      scod(8) = 'BROT'
      scod(9) = 'ASO'
      scod(10) = 'O'
      scod(11) = 'P'
      scod(12) = 'Q'
      scod(13) = 'DMOM'
      scod(14) = 'EFIELD'
*  then all the system dependent logical variables
*  in the same order as in the common block /cosysl/
      islcod = 0
      nscode = 14
*  set default values for doublet pi scattering
      nterm = 2
      if (iread .eq. 0) then
        jmax = 3
        igu = 1
        isa = 0
        npar = 2
        nman = 0
        imult = 2
        mproj(1) = 0
        mproj(2) = 2
        lammin(1)= 1
        lammax(1) = 1
        lammin(2) = 2
        lammax(2) = 2
        niout=4
        indout(1)=-1
        indout(2)=1
        indout(3)=-2
        indout(4)=2
      endif
      if (iread .eq. 0) return
*  line 19
*  line 20
      read (8, *, err=888) jmax, igu, isa, npar, imult, nman
*  line 21
      read (8, *, err=888) brot, aso
*  line 22
      read (8, *, err=888) o, p, q
*  line 23
      read (8, *, err=888) dmom, efield
*  line 16 name of file containing potential parameters
      line=' '
      if(.not.readpt.or.iread.eq.0) then
        call loapot(1,' ')
        return
      endif
      read (8, 285, end=286) line
285   format (a)
      goto 286
* here if read error occurs
888   write(6,1000)
1000  format(/'   *** ERROR DURING READ FROM INPUT FILE ***')
      return
      entry ptrpi (fname,readpt)
      line = fname
      readpt = .true.
286   if (readpt) then
        l=1
        call parse(line,l,filnam,lc)
        if(lc.eq.0) then
          write(6,1020)
1020      format(' FILENAME MISSING FOR POTENTIAL INPUT')
        end if
        j=index(filnam(1:lc),dot)
        if(j.eq.0) then
          call gennam(potfil,filnam,0,'BIN',lc)
          filnam = potfil
        end if
        potfil=filnam
        inquire(file=filnam,exist=existf)
        if(.not.existf) then
          write(6,1025) filnam(1:lc)
1025      format(' FILE ',(a),' NOT FOUND')
          return
        end if
* now call loapot(iunit,filnam) routine to read potential parameters
        call loapot(1,filnam)
      end if
      close (8)
      irpot=1
      return
      entry savpi (readpt)
*  save input parameters for singlet, doublet or
*  triplet-pi molecule + atom scattering
*  the order of the write statements should be identical to the read
*  statements above. for consistency with the data file written by
*  gendat, format statements should reserve the first 30 spaces for
*  data, spaces 31-33 should be left blank, and the names of the
*  variables should be printed in spaces 34-80
*  line 18:
*  line 19
*  line 20
      write (8, 130) jmax, igu, isa, npar, imult, nman
130   format (6i4, 6x, '   jmax, igu, isa, npar, imult, nman')
*  line 21
      write (8, 140) brot, aso
140   format (f10.5,f10.4, 10x, '   brot, aso')
*  line 22
      write (8, 150) o, p, q
150   format (3(1pg12.4), ' o, p, q')
*  line 23
      write(8, 160) dmom, efield
160   format (2f10.5, 10x, '   dmom, efield')
*  line 16
      write (8, 285) potfil
      return
      end
*  -----------------------------------------------------------------------
      subroutine systp (irpot, readpt, iread)
*  subroutine to read in system dependent parameters for symmetric top
*   + atom scattering
*  if iread = 1 read data from input file
*  if iread = 0 return after defining variable names
*  current revision date: 3-nov-1995 by g.vs
* NOTE THAT THIS VERSION DOES NOT USE DATA FROM VFIT (FOLLMEG FORM)
*  -----------------------------------------------------------------------
*  variables in common bloc /cosysr/
*    isrcod:   total number of real system dependent variables
*    brot:     A=B rotational constant for prolate top
*    crot:     C rotational constant for prolate top
*  variables in common block /cosysi/
*    nscod:    total number of variable names which are passed to HINPUT
*              nscod must equal isrcod + isicod + 3
*    isicod:   total number of integer system dependent variables
*    nterm:    number of different associated legendre terms in
*              expansion of potential
*    numpot:   the number of the potential used, this variable is passed
*              to the pot subroutine
*    ipotsy:   cylindrical symmetry of potential.  Only terms with mu equal to
*              an integral multiple of ipotsy can be included in the potential
*              Example:  for NH3, ipotsy = 3
*    iop:      ortho/para label for molecular states. If ihomo=.true. then onl
*              para states will be included if iop=1 and only ortho states if
*              iop=-1
*    ninv:     number of inversion doublets included
*              if ninv = +1, only + inversion levels included
*              if ninv = -1, only - inversion levels included
*              if ninv = 2, both inversion levels included
*    kmax:     the maximum projection quantum number included
*    jmax0:    the maximum rotational angular momenta for the k=0 stack
*    jmax1:    the maximum rotational angular momenta for the k=1 stack
*    jmax2:    the maximum rotational angular momenta for the k=2 stack
*    jmax3:    the maximum rotational angular momenta for the k=3 stack
*    jmax4:    the maximum rotational angular momenta for the k=4 stack
*    jmax5:    the maximum rotational angular momenta for the k=5 stack
*    jmax6:    the maximum rotational angular momenta for the k=6 stack
*  variable in common /cosys/
*    scod:    character*8 array of dimension nscode, which contains names
*             of all system dependent parameters.  Note that the ordering
*             of the variable names in scod must correspond to the ordering
*             of the variable names in cosysi followed by the ordering of
*             variable names in cosysr followed by LAMMIN, LAMMAX, and MPROJ
*  -----------------------------------------------------------------------
      logical readpt, existf
      double precision brot, crot
      integer lammax, lammin, mproj, numpot
      integer icod, ircod, lencod
      integer i, iop, iread, irpot, isicod, isrcod, ipotsy, jmax0,
     :        jmax1, jmax2, jmax3, jmax4, jmax5, jmax6,
     :        jmax7, jmax8, jmax9, jmax10, jmax11, jmax12, kmax,
     :        ninv, nscode, nterm
      character*8 scod
      character*1 dot
      character*(*) fname
      character*40 line, filnam, potfil
      parameter (icod=19, ircod=2)
      parameter (lencod = icod + ircod + 3)
      include "common/parbas"
      common /coiout/ niout, indout(20)
      common /cosys/ scod(lencod)
      common /cosysi/ nscode, isicod, nterm, numpot, ipotsy, iop, ninv,
     :                kmax, jmax0, jmax1, jmax2, jmax3, jmax4, jmax5,
     :                jmax6, jmax7, jmax8, jmax9, jmax10, jmax11,
     :                jmax12
      common /cosysr/ isrcod, junkr, brot, crot
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
      scod(5)='NINV'
      scod(6)='KMAX'
      scod(7)='JMAX0'
      scod(8)='JMAX1'
      scod(9)='JMAX2'
      scod(10)='JMAX3'
      scod(11)='JMAX4'
      scod(12)='JMAX5'
      scod(13)='JMAX6'
      scod(14)='JMAX7'
      scod(15)='JMAX8'
      scod(16)='JMAX9'
      scod(17)='JMAX10'
      scod(18)='JMAX11'
      scod(19)='JMAX12'
      scod(20)='BROT'
      scod(21)='CROT'
      scod(22)='LAMMIN'
      scod(23)='LAMMAX'
      scod(24)='MPROJ'
      nscode = lencod
      isicod = icod
      isrcod = ircod
      irpot = 1
*  set default values for symmetric top scattering
        nterm = 4
        mproj(1) = 0
        mproj(2) = 3
        mproj(3) = 6
        mproj(4) = 9
      if (iread .eq. 0) then
        lammin(1)= 1
        lammin(2)= 3
        lammin(3)= 6
        lammin(4)=9
        lammax(1) = 1
        lammax(2) = 3
        lammax(3) = 6
        lammax(4) = 9
        jmax0 = 0
        jmax1 = 1
        jmax2 = 2
        jmax3 = 3
        jmax4 = 4
        jmax5 = 5
        jmax6 = 6
        jmax7 = 7
        jmax8 = 8
        jmax9 = 9
        jmax10 = 10
        jmax11 = 11
        jmax12 = 12
        ipotsy = 3
        iop = 1
        ninv = 1
        kmax = 4
        niout=5
        indout(1)=0
        indout(2)=-3
        indout(3)=+3
        indout(4)=-6
        indout(5)=6
      endif
      potfil='nopot'
      if (iread .eq. 0) return
*  line 18
      read (8, *, err=80) ipotsy, iop, ninv
*  line 19
*  line 20
      read (8, *, err=80) kmax, jmax0, jmax1, jmax2, jmax3, jmax4,
     :            jmax5, jmax6, jmax7, jmax8, jmax9, jmax10, jmax11,
     :            jmax12
*  line 21
      read (8, *, err=80) brot, crot
      if(.not.readpt.or.iread.eq.0) then
        call loapot(1,' ')
        return
      endif
      read (8, 60, end=100) line
60    format (a)
      goto 100
* here if read error occurs
80    write(6,90)
90    format(/'   *** ERROR DURING READ FROM INPUT FILE ***')
      return
      entry ptrstp (fname,readpt)
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
        inquire(file=filnam,exist=existf)
        if(.not.existf) then
          write(6,105) filnam(1:lc)
105      format(' FILE ',(a),' NOT FOUND')
         return
        end if
* now call loapot(iunit,filnam) routine to read potential parameters
        call loapot(1,filnam)
      endif
      close (8)
      return
      entry savstp (readpt)
*  save input parameters for symmetric top + atom scattering
*  the order of the write statements should be identical to the read statement
*  above. for consistency with the data file written by gendat, format
*  statements should reserve the first 30 spaces for data, spaces 31-33 should
*  be left blank, and the names of the variables should be printed in spaces
*  34-80
*  line 18:
      write (8, 220) ipotsy, iop, ninv
220   format (4i4, 14x,'   ipotsy, iop, ninv')
*  line 20
      write (8, 230) kmax, jmax0, jmax1, jmax2, jmax3, jmax4, jmax5,
     :            jmax6
230   format (i3, 7i4, '  kmax, jmax0 --- jmax6')
      write (8, 231) jmax7, jmax8, jmax9, jmax10, jmax11, jmax12
231   format (7x, 6i4, '  jmax7 --- jmax12')
*  line 21
      write (8, 250) brot, crot
250   format(2f12.6, 6x, '   brot, crot')
      write (8, 60) potfil
      return
      end
* -----------------------------------------------------------------------
      subroutine sy13p (irpot, readp, iread)
*  subroutine to read in system dependent parameters for collisions of
* atom in singlet and/or triplet P electronic state with closed shell atom
*  if iread = 1 read data from input file
*  if iread = 0 return after defining variable names
*  current revision date: 17-jun-1992 by mha
*  -----------------------------------------------------------------------
*  variables in common /cosysr/
*    isrcod:  number of real parameters
*    en(1):    asymptotic energy of j=0 fine-structure level of triplet
*    en(2):    asymptotic energy of j=1 fine-structure level of triplet
*    en(3):    asymptotic energy of j=2 fine-structure level of triplet
*    en(4):    asymptotic energy of singlet state (cm-1)
*    cmix:     mixing coefficient of J=1 singlet and triplet levels
*    de, re, be, rl, c:  msv parameters for 3Pi, 3Sig, 1Pi, 1Sig
*  variable in common /cosysi/
*    nscode:  total number of system dependent parameters
*             nscode = isicod + isrcod +3
*    isicod:  number of integer system dependent parameters
*    nterm:    number of different electronic coupling terms
*              this should be 1 here
*    nstate:   number of electronic states included
*              nstate=0:   just singlet state
*              nstate=1:   just triplet state
*              nstate=2:   both singlet and triplet state
*  variable in common /cosys/
*    scod:    character*8 array of dimension lcode, which contains names
*             of all system dependent parameters
*  line 16:
*    filnam:  name of file containing potential parameters
*
*  subroutines called: loapot(iunit,filnam)
*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include "common/parsys"
      integer irpot
      logical readp
      logical airyfl, airypr, logwr, swrit, t2writ, writs, wrpart,
     :        partw, xsecwr, wrxsec, noprin, chlist, ipos, flaghf,
     :        csflag, flagsu, rsflag, t2test, existf, logdfl, batch,
     :        readpt, ihomo, bastst, twomol, lpar

      character*1 dot
      character*8 scod
      character*(*) fname
      character*40 filnam, line, potfil
      common /cosys/ scod(lencod)
      common /cosysi/ nscode, isicod, nterm, nstate, ipol, npot
      common /cosysr/ isrcod, junkr, en(4), de(4), re(4), be(4),
     :                        rl(4), cl(4), cmix, alphg, rgaus,
     :                        agaus, demor, remor, bemor,dissmor
      common /conlam/ nlam
      include "common/parbas"
      common /coskip/ nskip,iskip
      common /colpar/ airyfl, airypr, bastst, batch, chlist, csflag,
     :                flaghf, flagsu, ihomo, ipos, logdfl, logwr,
     :                noprin, partw, readpt, rsflag, swrit,
     :                t2test, t2writ, twomol, writs, wrpart, wrxsec,
     :                xsecwr,lpar(3)
      include "common/comdot"
      save potfil
*  number and names of system dependent parameters
*  first all the system dependent integer variables
*  in the same order as in the common block /cosysi/
*  variable names must be in uppercase of maximum length 6 characters
      isicod = 4
      scod(1) = 'NTERM'
      scod(2) = 'NSTATE'
      scod(3) = 'IPOL'
      scod(4) = 'NPOT'
*  then all the system dependent real variables
*  in the same order as in the common block /cosysr/
      isrcod = 32
      scod(5) = 'E-3P0'
      scod(6) = 'E-3P1'
      scod(7) = 'E-3P2'
      scod(8) = 'E-1P1'
      scod(9) = 'DE(3PI)'
      scod(10) = 'DE(3SG)'
      scod(11) = 'DE(1PI)'
      scod(12) = 'DE(1SG)'
      scod(13) = 'RE(3PI)'
      scod(14) = 'RE(3SG)'
      scod(15) = 'RE(1PI)'
      scod(16) = 'RE(1SG)'
      scod(17) = 'BE(3PI)'
      scod(18) = 'BE(3SG)'
      scod(19) = 'BE(1PI)'
      scod(20) = 'BE(1SG)'
      scod(21) = 'RL(3PI)'
      scod(22) = 'RL(3SG)'
      scod(23) = 'RL(1PI)'
      scod(24) = 'RL(1SG)'
      scod(25) = 'C(3PI)'
      scod(26) = 'C(3SG)'
      scod(27) = 'C(1PI)'
      scod(28) = 'C(1SG)'
      scod(29) = 'CMIX'
      scod(30) = 'ALPHG'
      scod(31) = 'RGAUS'
      scod(32) = 'AGAUS'
      scod(33) = 'DEMOR'
      scod(34) = 'REMOR'
      scod(35) = 'BEMOR'
      scod(36) = 'DISSMOR'
      nscode = isicod + isrcod
*  set default values for 1/3 P atom scattering
      nterm = 1
      nstate = 2
      ipol=1
      lammin(1)= 1
      lammax(1) = 4
      mproj(1)=0
      niout=2
      indout(1)=0
      indout(2)=1
      potfil = ' '
      if(iread.eq.0) return
      read (8, *, err=888) nstate, ipol, npot
      read (8, *, err=888) (en(i), i=1,3)
      read (8, *, err=888) en(4), cmix
      read (8, *, err=888)  de(1), re(1), be(1), rl(1), cl(1)
      read (8, *, err=888)  de(2), re(2), be(2), rl(2), cl(2)
      read (8, *, err=888)  de(3), re(3), be(3), rl(3), cl(3)
      read (8, *, err=888)  de(4), re(4), be(4), rl(4), cl(4)
      read (8, *, err=888)  demor, remor, bemor, dissmor
      read (8, *, err=888)  rgaus, agaus, alphg
      rms=0.d0
      line=' '
      if(.not.readp.or.iread.eq.0) then
        call loapot(1,' ')
        return
      endif
      read (8, 285, end=286) line
285   format (a)
      goto 286
* here if read error occurs
888   write(6,1000)
1000  format(/'   *** ERROR DURING READ FROM INPUT FILE ***')
      return
* --------------------------------------------------------------
      entry ptr13p (fname,readp)
      line = fname
      readp = .true.
286   if (readp) then
        l=1
        call parse(line,l,filnam,lc)
        if(lc.eq.0) then
          write(6,1020)
1020      format(' FILENAME MISSING FOR POTENTIAL INPUT')
        end if
        j=index(filnam(1:lc),dot)
        if(j.eq.0) then
           call gennam(potfil,filnam,0,'BIN',lc)
           filnam = potfil
        end if
        potfil=filnam
        inquire(file=filnam,exist=existf)
        if(.not.existf) then
          write(6,1025) filnam(1:lc)
1025      format(' FILE ',(a),' NOT FOUND')
          return
        end if
* now call loapot(iunit,filnam) routine to read potential parameters
        call loapot(1,filnam)
      end if
      close (8)
      irpot=1
      return
      entry sav13p (readp)
*  save input parameters for singlet-sigma + atom scattering
      write (8, 300) nstate, ipol, npot
300   format(3i4, 18x,'   nstate, ipol, npot')
      write (8, 310) (en(i), i=1,3)
310   format(3(1pg12.4),4x,'   E-3P0, 3P1, 3P2')
      write (8, 315) en(4), cmix
315   format (2(1pg12.4),16x,'   E-1P1, CMIX')
      write (8, 320) de(1), re(1), be(1), rl(1), cl(1)
320   format (5(1pg12.4),' 3 Pi Parameters')
      write (8, 325) de(2), re(2), be(2), rl(2), cl(2)
325   format (5(1pg12.4),' 3 Sig Parameters')
      write (8, 330) de(3), re(3), be(3), rl(3), cl(3)
330   format (5(1pg12.4),' 1 Pi Parameters')
      write (8, 335) de(4), re(4), be(4), rl(4), cl(4)
335   format (5(1pg12.4),' 1 Sig Parameters')
      write (8, 336) demor, remor, bemor, dissmor
336   format(4(1pg12.4),12x,' 1 Sigma Morse')
      write (8, 337) rgaus, agaus, alphg
337   format (3(1pg12.4),24x,' Gaussian Coupling')
      write (8, 285) potfil
      return
      end
* --------------------------------------------------------------
      subroutine sy2mol (irpot, readpt, iread)
*  subroutine to read in system dependent parameters for hf-hf collisions
*  if iread = 1 read data from input file
*  if iread = 0 return after defining variable names
*  authors:  millard alexander and peter vohralik
*  current revision date: 19-aug-1991
*  -----------------------------------------------------------------------
*  variables in common block /cotwo/
*    numj:     number of j1-j2 values
*    nj1j2:    specific j1-j2 values (up to a maximum of 50)
*              N.B. this dimension is set in hiba2mol
*  variables in common block /cosysr/
*    isrcod:   number of real system dependent variables
*    brot:     rotational constant in cm-1
*    drot:     centrifugal distortion constant in cm-1
*    hrot:     next centrifugal distortion constant in cm-1
*  variables in common block /cosysi/
*    nscode:   total number of system dependent variables
*             nscode = isicod + isrcod
*    isicod:   numbr of integer system dependent variables
*    nterm:    number of different lambda terms included in potential
*              see below
*    nsym:     interchange symmetry of included channels
*  variable in common block /conlam/
*    nlam:     the total number of angular coupling terms
*              nlam is set equal to nterm; see above
*     nterm = 1, include lam1=1  lam2=1  lam=2   dipole-dipole
*     nterm = 2, include dipole-dipole and
*                dipole-quadrupole (lam1=1  lam2=2  lam=3)
*     nterm = 3, include dipole-dipole, dipole-quadrupole, and
*                short-range term (lam1=0  lam2=1  lam=1)
*  line 13:
*    nterm:    total number of anistropic terms in potential
*  line 14:
*    brot, drot, hrot:    rotational constants of the molecule in cm-1
*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include "common/parsys"
      logical readpt, existf
      integer irpot
      character*1 dot
      character*8 scod
      character*(*) fname
      character*40 filnam, line, potfil
      common /cosys/ scod(lencod)
      common /cosysi/ nscode, isicod, nterm, nsym
      common /cosysr/ isrcod, junkr, brot, drot, hrot
      common /conlam/ nlam
      common /cotwo/ numj,nj1j2(3)
      save potfil
      include "common/comdot"
      nscode = lencod
*     number and names of system dependent parameters
      isicod = 2
      isrcod = 3
      nscode=isicod+isrcod
      scod(1) = 'NTERM'
      scod(2) = 'NSYM'
      scod(3) = 'BROT'
      scod(4) = 'DROT'
      scod(5) = 'HROT'
*  set default values for hf-hf collisions
      nterm = 1
      nsym = 1
      numj = 3
      brot =  20.559741
      drot = 2.1204498e-3
      hrot = 0.168e-6
      nj1j2(1) = 00
      nj1j2(2) = 11
      nj1j2(3) = 02
      if(iread.eq.0) return
*  line 1
      iline=1
      read (8, *, err=888) nterm, nsym, numj
      nlam = nterm
      if (numj .le. 0) then
        write (6, 80) numj
80      format(' ** NUMJ = ',i3,'; .LE. O IN SY2MOL; ABORT')
        call exit
      endif
*  line 2
      iline=2
      read (8, *, err=888) brot, drot, hrot
      if (numj .le. 20) then
        iline=3
        read (8, *, err=888) (nj1j2(i), i=1,numj)
      else
        ilow=1
        ihigh=20
        do  100 ij = 1, numj, 20
          iline=iline+1
          itop=min(ihigh,numj)
          read (8, *, err=888) (nj1j2(i), i=ilow,itop)
          ilow=itop+1
          ihigh=ihigh+20
100     continue
      endif
*  line 15 name of file containing potential parameters
      iline=iline+1
      line=' '
      if(.not.readpt.or.iread.eq.0) then
        call loapot(1,' ')
        return
      endif
      read (8, 285, end=286) line
285   format (a)
      goto 286
* here if read error occurs
888   write(6,1000) iline
1000  format(/'   *** ERROR DURING READ OF SY2MOL INPUT; LINE: ',i2)
      return
      entry ptr2mol (fname,readpt)
      line = fname
      readpt = .true.
286   if (readpt) then
        l=1
        call parse(line,l,filnam,lc)
        if(lc.eq.0) then
          write(6,1020)
1020      format(' FILENAME MISSING FOR POTENTIAL INPUT')
        end if
        j=index(filnam(1:lc),dot)
        if(j.eq.0) then
          call gennam(potfil,filnam,0,'BIN',lc)
          filnam = potfil
        end if
        potfil=filnam
        inquire(file=filnam,exist=existf)
        if(.not.existf) then
          write(6,1025) filnam(1:lc)
1025      format(' FILE ',(a),' NOT FOUND')
          return
        end if
* now call loapot(iunit,filnam) routine to read potential parameters
        call loapot(1,filnam)
        end if
      close (8)
      irpot=1
      return
      entry sav2mol (readpt)
*  save input parameters for hf-hf scattering
*  line 13:
      write (8, 300) nterm, nsym,numj
300   format (3i4,18x,'   nterm, nsym, numj')
*  line 14
      write (8, 320) brot, drot, hrot
320   format(f8.4, 2g11.4, '   brot, drot, hrot')
*  line 15
      if (numj .le. 20) then
        iline=3
        write (8, 330) (nj1j2(i)/10,mod(nj1j2(i),10), i=1,numj)
330   format (1x,20(2i1,'  '),t65,'nj1j2')
      else
        ilow=1
        ihigh=20
        do  350 ij = 1, numj, 20
          iline=iline+1
          itop=min(ihigh,numj)
          write (8, 330) (nj1j2(i)/10,mod(nj1j2(i),10),
     :                    i=ilow,itop)
          ilow=itop+1
          ihigh=ihigh+20
350     continue
      endif
      return
      end
*  -----------------------------------------------------------------------
      subroutine systpln (irpot, readpt, iread)
*  subroutine to read in system dependent parameters for symmetric top
*   + linear molecule scattering
*  if iread = 1 read data from input file
*  if iread = 0 return after defining variable names
*  current revision date: 9-mar-1992 by cr
* NOTE THAT THIS VERSION DOES NOT USE DATA FROM VFIT (FOLLMEG FORM)
*  -----------------------------------------------------------------------
*  variables in common bloc /cosysr/
*    isrcod:   total number of real system dependent variables
*    brot:     A=B rotational constant for prolate top
*    crot:     C rotational constant for prolate top
*    drot:     rotational constant for linear molecule
*  variables in common block /cosysi/
*    nscod:    total number of variable names which are passed to HINPUT
*              nscod must equal isrcod + isicod + 3
*    isicod:   total number of integer system dependent variables
*    nterm:    number of different associated legendre terms in
*              expansion of potential
*    numpot:   the number of the potential used, this variable is passed
*              to the pot subroutine
*    ipotsy:   cylindrical symmetry of potential.  Only terms with mu equal to
*              an integral multiple of ipotsy can be included in the potential
*              Example:  for NH3, ipotsy = 3
*    iop:      ortho/para label for molecular states. If ihomo=.true. then onl
*              para states will be included if iop=1 and only ortho states if
*              iop=-1
*    ninv:     number of inversion doublets included
*              if ninv = +1, only + inversion levels included
*              if ninv = -1, only - inversion levels included
*              if ninv = 2, both inversion levels included
*    kmax:     the maximum projection quantum number included
*    jmax0:    the maximum rotational angular momenta for the k=0 stack
*    jmax1:    the maximum rotational angular momenta for the k=1 stack
*    jmax2:    the maximum rotational angular momenta for the k=2 stack
*    jmax3:    the maximum rotational angular momenta for the k=3 stack
*    jmax4:    the maximum rotational angular momenta for the k=4 stack
*    jmax5:    the maximum rotational angular momenta for the k=5 stack
*    jmax6:    the maximum rotational angular momenta for the k=6 stack
*    ipotsy2:  symmetry of potential. if linear molecule is homonuclear
*              then ipotsy=2 and only terms with lambda2  even can be
*              included in the potential,else ipotsy=1.
*    j2max:    the maximum rotational angular momentum for linear
*              molecule
*    j2min:    the minimum rotational angular momentum for linear
*              molecule
*  variable in common /cosys/
*    scod:    character*8 array of dimension nscode, which contains names
*             of all system dependent parameters.  Note that the ordering
*             of the variable names in scod must correspond to the ordering
*             of the variable names in cosysi followed by the ordering of
*             variable names in cosysr followed by LAMMIN, LAMMAX, MPROJ,
*             LAM2 and MPROJ2.
*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical readpt, existf, twomol
      double precision brot, crot
      integer lammax, lammin, mproj, numpot
      integer icod, ircod, lencod
      integer i, iop, iread, irpot, isicod, isrcod, ipotsy, jmax0,
     :        jmax1, jmax2, jmax3, jmax4, jmax5, jmax6, kmax,
     :        ninv, nscode, nterm
      character*8 scod
      character*1 dot
      character*(*) fname
      character*40 line, filnam, potfil
      parameter (icod=16, ircod=3)
      parameter (lencod = icod + ircod + 5)
      include "common/parbas"
      common /coiout/ niout, indout(20)
      common /cosys/ scod(lencod)
      common /cosysi/ nscode, isicod, nterm, numpot, ipotsy, iop, ninv,
     :                kmax, jmax0, jmax1, jmax2, jmax3, jmax4, jmax5,
     :                jmax6,ipotsy2, j2max, j2min
      common /co2mol/ twomol
      common /cosysr/ isrcod, junkr, brot, crot, drot
      include "common/comdot"
      save potfil
      twomol = .true.
*  number and names of system dependent parameters
*  first all the system dependent integer variables
*  in the same order as in the common block /cosysi/
*  variable names must be in uppercase of maximum length 6 characters
*  NTERM must be the first variable
*  followed by the system dependent real variables
*  in the same order as in the common block /cosysr/
*  then the three variable names LAMMIN, LAMMAX, MPROJ, in that order
      scod(1)='NTERM'
      scod(2)='NUMPOT'
      scod(3)='IPOTSY'
      scod(4)='IOP'
      scod(5)='NINV'
      scod(6)='KMAX'
      scod(7)='JMAX0'
      scod(8)='JMAX1'
      scod(9)='JMAX2'
      scod(10)='JMAX3'
      scod(11)='JMAX4'
      scod(12)='JMAX5'
      scod(13)='JMAX6'
      scod(14)='IPOTSY2'
      scod(15)='J2MAX'
      scod(16)='J2MIN'
      scod(17)='BROT'
      scod(18)='CROT'
      scod(19)='DROT'
      scod(20)='LAMMIN'
      scod(21)='LAMMAX'
      scod(22)='MPROJ'
      scod(23)='LAM2'
      scod(24)='M2PROJ'
      nscode = lencod
      isicod = icod
      isrcod = ircod
      irpot = 1
*  set default values for symmetric top scattering
      nterm=3
*      nterm=12
      do j=1, 4
        do i=1, 3
          k=3*(j-1)+i
          mproj(k) = 3*(i-1)
          if(j.le.2) m2proj(k) = 0
          if(j.eq.3) m2proj(k) = 1
          if(j.eq.4) m2proj(k) = 2
          if(j.eq.1) lam2(k) = 0
          if(j.ge.2) lam2(k) = 2
          lammin(k) = max(mproj(k), m2proj(k))
          lammax(k) = 6
        enddo
      enddo
      lammin(1)=1
      if (iread .eq. 0) then

        jmax0 = 0
        jmax1 = 1
        jmax2 = 2
        jmax3 = 3
        jmax4 = 4
        jmax5 = 5
        jmax6 = 6
        ipotsy = 3
        iop = 1
        ninv = 1
        kmax =6
        ipotsy2=1
        j2min=0
        j2max=2
        niout=5
        indout(1)=0
        indout(2)=-3
        indout(3)=+3
        indout(4)=-6
        indout(5)=6
      endif
      potfil='nopot'
      if (iread .eq. 0) return
*  line 18
      read (8, *, err=80) ipotsy, iop, ninv, ipotsy2
*  line 19
*  line 20
      read (8, *, err=80) kmax, jmax0, jmax1, jmax2, jmax3, jmax4,
     $     jmax5, jmax6
*  line 21
      read (8, *, err=80) j2min, j2max
*  line 22
      read (8, *, err=80) brot, crot
*  line 23
      read(8, *, err=80) drot
      if(.not.readpt.or.iread.eq.0) then
        call loapot(1,' ')
        return
      endif
      read (8, 60, end=100) line
60    format (a)
      goto 100
* here if read error occurs
80    write(6,90)
90    format(/'   *** ERROR DURING READ FROM INPUT FILE ***')
      return
*************************************************************
      entry ptrstpln (fname,readpt)
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
        inquire(file=filnam,exist=existf)
        if(.not.existf) then
          write(6,105) filnam(1:lc)
105      format(' FILE ',(a),' NOT FOUND')
         return
        end if
* now call loapot(iunit,filnam) routine to read potential parameters
        call loapot(1,filnam)
      endif
      close (8)
      return
      entry savstpln (readpt)
*  save input parameters for symmetric top + linear molecule scattering
*  the order of the write statements should be identical to the read statement
*  above. for consistency with the data file written by gendat, format
*  statements should reserve the first 30 spaces for data, spaces 31-33 should
*  be left blank, and the names of the variables should be printed in spaces
*  34-80
*  line 18:
      write (8, 220) ipotsy, iop, ninv, ipotsy2
220   format (4i4, 14x,'   ipotsy, iop, ninv, ipotsy2')
*  line 20
      write (8, 230) kmax, jmax0, jmax1, jmax2, jmax3, jmax4, jmax5,
     :            jmax6
230   format (i3, 7i4, '  kmax, jmax0 --- jmax6')
*  line 21
      write(8,231) j2min, j2max
231   format (2i4, 22x,'   j2min, j2max')
      write (8, 250) brot, crot
250   format(2f12.6, 6x, '   brot, crot' )
      write(8, 251) drot
251   format(f12.6, 18x,'   drot')
      write (8, 60) potfil
      return
      end
* -----------------------------------------------------------------------
      subroutine sy22p (irpot, readp, iread)
*  subroutine to read in system dependent parameters for collisions of
* atom in doublet S state with atom in doublet P state
*  if iread = 1 read data from input file
*  if iread = 0 return after defining variable names
*  current revision date: 14-mar-1997 by mha
*  -----------------------------------------------------------------------
*  variables in common /cosysr/
*    isrcod:  number of real parameters
*    aso:     spin-orbit constant in 2P atom (cm-1)
*  variable in common /cosysi/
*    nscode:  total number of system dependent parameters
*             nscode = isicod + isrcod +3
*    isicod:  number of integer system dependent parameters
*    nterm:    number of different electronic coupling terms
*              this should be 1 here
*    nphoto:  this should be 1 here
*    nvib:    vibrational quantum number of initial state (0-6)
*  variable in common /cosys/
*    scod:    character*8 array of dimension lcode, which contains names
*             of all system dependent parameters
*  line 16:
*    filnam:  name of file containing potential parameters
*
*  subroutines called: loapot(iunit,filnam)
*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include "common/parsys"
      integer irpot
      logical readp
      logical airyfl, airypr, logwr, swrit, t2writ, writs, wrpart,
     :        partw, xsecwr, wrxsec, noprin, chlist, ipos, flaghf,
     :        csflag, flagsu, rsflag, t2test, existf, logdfl, batch,
     :        readpt, ihomo, bastst, twomol,lpar
      character*1 dot
      character*8 scod
      character*(*) fname
      character*40 filnam, line, potfil
      common /cosys/ scod(lencod)
      common /cosysi/ nscode,isicod,iscod(maxpar)
      common /coipar/ jtot1,jtot2,jtotd,jlpar
      common /cosysr/ isrcod, junkr, rcod(maxpar)
      common /conlam/ nlam
      include "common/parbas"
      common /coskip/ nskip,iskip
      common /colpar/ airyfl, airypr, bastst, batch, chlist, csflag,
     :                flaghf, flagsu, ihomo, ipos, logdfl, logwr,
     :                noprin, partw, readpt, rsflag, swrit,
     :                t2test, t2writ, twomol, writs, wrpart, wrxsec,
     :                xsecwr,lpar(3)
      include "common/comdot"
      save potfil
*  number and names of system dependent parameters
*  first all the system dependent integer variables
*  in the same order as in the common block /cosysi/
*  variable names must be in uppercase of maximum length 6 characters
      isicod = 4
      scod(1) = 'NTERM'
      scod(2) = 'NPHOTO'
      scod(3) = 'NVIB'
      scod(4) = 'IBRAN'
*  then all the system dependent real variables
*  in the same order as in the common block /cosysr/
      isrcod = 1
      scod(5) = 'A-SO'
      nscode = isicod + isrcod
*  set default values for 2S + 2P atom scattering
      nterm = 1
      iscod(1)=1
      iscod(2)=1
      iscod(3)=0
      iscod(4)=0
      lammin(1)= 1
      lammax(1) = 13
      rcod(1)=881.
      mproj(1)=0
      niout=2
      indout(1)=0
      indout(2)=1
      potfil = ' '
      if(iread.eq.0) return
      read (8, *, err=888) nphoto, nvib, ibran, aso
      iscod(2)=nphoto
      iscod(3)=nvib
      iscod(4)=ibran
      if (nvib .gt. 6) then
        write (6, 280) nvib
280     format ('*** NVIB = ',i2, '; SHOULD BE 0-6')
        go to 888
      endif
      if (iabs(ibran) .gt. 1) then
        write (6, 281) ibran
281     format ('*** IBRAN = ',i2, '; SHOULD BE -1, 0, or 1')
        go to 888
      endif
      if (jlpar .eq. -1) then
        if (ibran .eq. 0) then
          write (6, 282) jlpar, ibran
282       format ('*** JLPAR = ',i2, '; IBRAN = ', i2,
     :            'Q-BRANCH NOT ALLOWED')
          goto 888
        endif
      else
        if (ibran .ne. 0) then
          write (6, 283) jlpar, ibran
283       format ('*** JLPAR = ',i2, '; IBRAN = ', i2,
     :            'P,R-BRANCH NOT ALLOWED')
          goto 888
        endif
      endif
      print *, ibran, jlpar
      rcod(1)=aso
      line=' '
      if(.not.readp.or.iread.eq.0) then
        call loapot(1,' ')
        return
      endif
      read (8, 285, end=286) line
285   format (a)
      goto 286
* here if read error occurs
888   write(6,1000)
1000  format(/'   *** ERROR DURING READ FROM INPUT FILE ***')
      return
* --------------------------------------------------------------
      entry ptr22p (fname,readp)
      line = fname
      readp = .true.
286   if (readp) then
        l=1
        call parse(line,l,filnam,lc)
        if(lc.eq.0) then
          write(6,1020)
1020      format(' FILENAME MISSING FOR POTENTIAL INPUT')
        end if
        j=index(filnam(1:lc),dot)
        if(j.eq.0) then
           call gennam(potfil,filnam,0,'BIN',lc)
           filnam = potfil
        end if
        potfil=filnam
        inquire(file=filnam,exist=existf)
        if(.not.existf) then
          write(6,1025) filnam(1:lc)
1025      format(' FILE ',(a),' NOT FOUND')
          return
        end if
* now call loapot(iunit,filnam) routine to read potential parameters
        call loapot(1,filnam)
      end if
      close (8)
      irpot=1
      return
      entry sav22p (readp)
*  save input parameters for singlet-sigma + atom scattering
      nphoto=iscod(2)
      nvib=iscod(3)
      ibran=iscod(4)
      aso=rcod(1)
      write (8, 300) nphoto, nvib, ibran, aso
300   format(i4,2(2x,i4),1pg12.5,3x,'  nphoto, nvib, ibran, a-so')
      write (8, 285) potfil
      return
      end
* --------------------------------------------------------------
      subroutine sy1del (irpot, readpt, iread)
*  subroutine to read in system dependent parameters for singlet-delta
*   + atom scattering using werner-follmeg potential form
*  if iread = 1 read data from input file
*  if iread = 0 return after defining variable names
*  current revision date: 20-sep-1993 by mha
*  -----------------------------------------------------------------------
*    nlam:      the total number of angular coupling terms
*  variables in common bloc /cosysr/
*    isrcod:   total number of real system dependent variables
*    brot:     rotational constant in cm-1
*    q:     lambda-doubling constant cm-1
*  variables in common block /cosysi/
*    nscode:  total number of system dependent parameters
*             nscode = isicod + isrcod + 3
*    isicod:   total number of integer system dependent variables
*    nterm:   number of different angular terms in potential
*              nb for singlet sigma molecule, nterm should be 1
*    jmax:     the maximum rotational angular momenta for each channel
*              in each spin-orbit manifold with convention
*              omega .le. j .le. jmax+0.5
*    igu:      permutation inversion symmetry of electronic state
*              igu=1 for gerade states, igu=-1 for ungerade states
*              for heteronuclear molecules igu should be +1
*    isa:      s/a label for molecular states. if ihomo=.true. then only
*              s states will be included if isa=1 and only a states if
*              isa=-1
*    npar:     number of symmetry doublets included (npar=2 will ensure
*              both lambda doublets)
*  variable in common /cosys/
*    scod:    character*8 array of dimension nscode, which contains names
*             of all system dependent parameters
*  variable in common block /coskip/
*   nskip  for a homonuclear molecule lamda is running in steps of nskip=2
*          for a heteronuclear molecule nskip=1
*
*   skip   same as nskip, used for consistency check
*
*  subroutines called: loapot(iunit,filnam)
      implicit double precision (a-h,o-z)
      include "common/parsys"
      logical readpt, existf
      integer irpot
      character*1 dot
      character*8 scod
      character*(*) fname
      character*40 filnam, line, potfil
      include "common/parbas"
      common /cosys/ scod(lencod)
      common /cosysi/ nscode, isicod, nterm, jmax, igu, isa, npar
      common /cosysr/ isrcod, junkr, brot,  q
      logical         airyfl, airypr, bastst, batch, chlist, csflag,
     :                flaghf, flagsu, ihomo,lpar
      common /colpar/ airyfl, airypr, bastst, batch, chlist, csflag,
     :                flaghf, flagsu, ihomo,lpar(18)
      common /conlam/ nlam
      common /coskip/ nskip,iskip
      save potfil
      include "common/comdot"
      isicod = 5
      isrcod = 2
      nscode = isicod+isrcod
      scod(1) = 'NTERM'
      scod(2) = 'JMAX'
      scod(3) = 'IGU'
      scod(4) = 'ISA'
      scod(5) = 'NPAR'
      scod(6) = 'BROT'
      scod(7) = 'Q'
      potfil = ' '
*  set default values for singlet delta scattering
      nterm = 2
      mproj(2) = 4
      if (iread .eq. 0) then
        mproj(1) = 0
        lammin(1)= 1
        lammax(1) = 1
        lammin(2) = 4
        lammax(2) = 4
        niout=4
        indout(1)=-1
        indout(2)=1
        indout(3)=-2
        indout(4)=2
        npar=2
        igu=1
        isa=0
        jmax=3
      endif
      if(iread.eq.0) return
*  line 13
      read (8, *, err=888) jmax, igu, isa, npar
*  line 14
      read (8, *, err=888) brot,  q
*  line 15 name of file containing potential parameters
      line=' '
      if(.not.readpt.or.iread.eq.0) then
        call loapot(1,' ')
        return
      endif
      read (8, 285, end=286) line
285   format (a)
      goto 286
* here if read error occurs
888   write(6,1000)
1000  format(/'   *** ERROR DURING READ FROM INPUT FILE ***')
      return
      entry ptr1del (fname,readpt)
      line = fname
      readpt = .true.
286   if (readpt) then
        l=1
        call parse(line,l,filnam,lc)
        if(lc.eq.0) then
          write(6,1020)
1020      format(' FILENAME MISSING FOR POTENTIAL INPUT')
        end if
        j=index(filnam(1:lc),dot)
        if(j.eq.0) then
          call gennam(potfil,filnam,0,'BIN',lc)
          filnam = potfil
        end if
        potfil=filnam
        inquire(file=filnam,exist=existf)
        if(.not.existf) then
          write(6,1025) filnam(1:lc)
1025      format(' FILE ',(a),' NOT FOUND')
          return
        end if
* now call loapot(iunit,filnam) routine to read potential parameters
        call loapot(1,filnam)
        end if
      close (8)
      irpot=1
      return
      entry sav1del (readpt)
*  save input parameters for singlet-delta + atom scattering
*  line 13:
      write (8, 315) jmax, igu, isa, npar
315   format(4i4,18x,'jmax, igu, isa, npar')
*  line 14
      write (8, 320) brot,  q
320   format(f10.5,f10.4,2(1pg11.3),' brot,  q')
*  line 15
      write (8, 285) potfil
      return
      end
*  -----------------------------------------------------------------------
      subroutine syh2p (irpot, readpt, iread)
*  subroutine to read in system dependent parameters for homonuclear
*   + 2P atom scattering
*  if iread = 1 read data from input file
*  if iread = 0 return after defining variable names
*  current revision date: 29-nov-1995 by mha
* NOTE THAT THIS VERSION DOES NOT USE DATA FROM VFIT (FOLLMEG FORM)
*  -----------------------------------------------------------------------
*  variables in common bloc /cosysr/
*    isrcod:   total number of real system dependent variables
*    brot:     rotational constant of molecule
*    aso:      spin-orbit constant of atom
*  variables in common block /cosysi/
*    nscod:    total number of variable names which are passed to HINPUT
*              nscod must equal isrcod + isicod + 3
*    isicod:   total number of integer system dependent variables
*    nterm:    number of different associated legendre terms in
*              expansion of potential
*    iop:      ortho/para label for molecular states. If ihomo=.true. then onl
*              para states will be included if iop=1 and only ortho states if
*              iop=-1
*    jmax:     the maximum rotational angular momenta for the diatomic
*  variable in common /cosys/
*    scod:    character*8 array of dimension nscode, which contains names
*             of all system dependent parameters.  Note that the ordering
*             of the variable names in scod must correspond to the ordering
*             of the variable names in cosysi followed by the ordering of
*             variable names in cosysr followed by LAMMIN, LAMMAX, and MPROJ
*  -----------------------------------------------------------------------
      logical readpt, existf
      double precision brot, aso
      character*1 dot
      character*8 scod
      character*(*) fname
      character*40 line, filnam, potfil
      parameter (icod=3, ircod=2)
      parameter (lencod = icod + ircod + 3)
      include "common/parbas"
      common /coiout/ niout, indout(20)
      common /cosys/ scod(lencod)
      common /cosysi/ nscode, isicod, nterm, iop, jmax
      common /cosysr/ isrcod, junkr, brot, aso
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
      scod(2)='IOP'
      scod(3)='JMAX'
      scod(4)='BROT'
      scod(5)='ASO'
      scod(6)='LAMMIN'
      scod(7)='LAMMAX'
      scod(8)='MPROJ'
      nscode = lencod
      isicod = icod
      isrcod = ircod
      irpot = 1
*  set default values for homonuclear+2P atom
        nterm = 1
        mproj(1) = 0
      if (iread .eq. 0) then
        lammin(1)= 1
        lammax(1) = 5
        jmax = 0
        iop = 1
      endif
      potfil='nopot'
      if (iread .eq. 0) return
*  line 18
      read (8, *, err=80) iop, jmax
*  line 21
      read (8, *, err=80) brot, aso
      if(.not.readpt.or.iread.eq.0) then
        call loapot(1,' ')
        return
      endif
      read (8, 60, end=100) line
60    format (a)
      goto 100
* here if read error occurs
80    write(6,90)
90    format(/'   *** ERROR DURING READ FROM INPUT FILE ***')
      return
      entry ptrh2p (fname,readpt)
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
        inquire(file=filnam,exist=existf)
        if(.not.existf) then
          write(6,105) filnam(1:lc)
105      format(' FILE ',(a),' NOT FOUND')
         return
        end if
* now call loapot(iunit,filnam) routine to read potential parameters
        call loapot(1,filnam)
      endif
      close (8)
      return
      entry savh2p (readpt)
*  save input parameters for symmetric top + atom scattering
*  the order of the write statements should be identical to the read statement
*  above. for consistency with the data file written by gendat, format
*  statements should reserve the first 30 spaces for data, spaces 31-33 should
*  be left blank, and the names of the variables should be printed in spaces
*  34-80
*  line 18:
      write (8, 220) iop, jmax
220   format (4i4, 14x,'   iop, jmax')
*  line 21
      write (8, 250) brot, aso
250   format(f12.4,f14.4, 6x, '   brot, aso')
      write (8, 60) potfil
      return
      end
*  -----------------------------------------------------------------------
      subroutine syh3p (irpot, readpt, iread)
*  subroutine to read in system dependent parameters for homonuclear
*   + 3p atom scattering
*  if iread = 1 read data from input file
*  if iread = 0 return after defining variable names
*  current revision date: 12-jun-1997 by mha
* NOTE THAT THIS VERSION DOES NOT USE DATA FROM VFIT (FOLLMEG FORM)
*  -----------------------------------------------------------------------
*  variables in common bloc /cosysr/
*    isrcod:   total number of real system dependent variables
*    brot:     rotational constant of molecule
*    aso1:      spin-orbit energy of j=1 state (relative to j=0)
*    aso2:      spin-orbit energy of j=2 state (relative to j=0)
*  variables in common block /cosysi/
*    nscod:    total number of variable names which are passed to HINPUT
*              nscod must equal isrcod + isicod + 3
*    isicod:   total number of integer system dependent variables
*    nterm:    number of different associated legendre terms in
*              expansion of potential
*    iop:      ortho/para label for molecular states. If ihomo=.true. then onl
*              para states will be included if iop=1 and only ortho states if
*              iop=-1
*    jmax:     the maximum rotational angular momenta for the diatomic
*  variable in common /cosys/
*    scod:    character*8 array of dimension nscode, which contains names
*             of all system dependent parameters.  Note that the ordering
*             of the variable names in scod must correspond to the ordering
*             of the variable names in cosysi followed by the ordering of
*             variable names in cosysr followed by LAMMIN, LAMMAX, and MPROJ
*  -----------------------------------------------------------------------
      logical readpt, existf
      double precision brot, aso1, aso2
      character*1 dot
      character*8 scod
      character*(*) fname
      character*40 line, filnam, potfil
      parameter (icod=3, ircod=3)
      parameter (lencod = icod + ircod + 3)
      include "common/parbas"
      common /coiout/ niout, indout(20)
      common /cosys/ scod(lencod)
      common /cosysi/ nscode, isicod, nterm, iop, jmax
      common /cosysr/ isrcod, junkr, brot, aso1, aso2
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
      scod(2)='IOP'
      scod(3)='JMAX'
      scod(4)='BROT'
      scod(5)='ASO1'
      scod(6)='ASO2'
      scod(7)='LAMMIN'
      scod(8)='LAMMAX'
      scod(9)='MPROJ'
      nscode = lencod
      isicod = icod
      isrcod = ircod
      irpot = 1
*  set default values for homonuclear+3p atom
        nterm = 1
        mproj(1) = 0
      if (iread .eq. 0) then
        lammin(1)= 1
        lammax(1) = 5
        jmax = 0
        iop = 1
      endif
      potfil='nopot'
      if (iread .eq. 0) return
*  line 18
      read (8, *, err=80) iop, jmax
*  line 21
      read (8, *, err=80) brot, aso1, aso2
      if(.not.readpt.or.iread.eq.0) then
        call loapot(1,' ')
        return
      endif
      read (8, 60, end=100) line
60    format (a)
      goto 100
* here if read error occurs
80    write(6,90)
90    format(/'   *** ERROR DURING READ FROM INPUT FILE ***')
      return
      entry ptrh3p (fname,readpt)
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
        inquire(file=filnam,exist=existf)
        if(.not.existf) then
          write(6,105) filnam(1:lc)
105      format(' FILE ',(a),' NOT FOUND')
         return
        end if
* now call loapot(iunit,filnam) routine to read potential parameters
        call loapot(1,filnam)
      endif
      close (8)
      return
      entry savh3p (readpt)
*  save input parameters for symmetric top + atom scattering
*  the order of the write statements should be identical to the read statement
*  above. for consistency with the data file written by gendat, format
*  statements should reserve the first 30 spaces for data, spaces 31-33 should
*  be left blank, and the names of the variables should be printed in spaces
*  34-80
*  line 18:
      write (8, 220) iop, jmax
220   format (4i4, 14x,'   iop, jmax')
*  line 21
      write (8, 250) brot, aso1, aso2
250   format(f12.4,2f14.4, 2x, '   brot, aso1, aso2')
      write (8, 60) potfil
      return
      end
* ----------------------------------------------------------------
      subroutine sy2del (irpot, readpt, iread)
*  subroutine to read in system dependent parameters for doublet-delta
*   + atom scattering using werner-follmeg potential form
*  if iread = 1 read data from input file
*  if iread = 0 return after defining variable names
*  current revision date: 13-oct-1999 by mha
*  NB, this subroutine is indentical to sy2pi except for name changes
*  -----------------------------------------------------------------------
*    nlam:      the total number of angular coupling terms
*  variables in common bloc /cosysr/
*    isrcod:   total number of real system dependent variables
*    brot:     rotational constant in cm-1
*    aso:      spin-orbit constant in cm-1
*    p, q:     lambda-doubling constants cm-1
*  variables in common block /cosysi/
*    nscode:  total number of system dependent parameters
*             nscode = isicod + isrcod + 3
*    isicod:   total number of integer system dependent variables
*    nterm:   number of different angular terms in potential
*              nb for singlet sigma molecule, nterm should be 1
*    jmax:     the maximum rotational angular momenta for each channel
*              in each spin-orbit manifold with convention
*              omega .le. j .le. jmax+0.5
*    igu:      permutation inversion symmetry of electronic state
*              igu=1 for gerade states, igu=-1 for ungerade states
*              for heteronuclear molecules igu should be +1
*    isa:      s/a label for molecular states. if ihomo=.true. then only
*              s states will be included if isa=1 and only a states if
*              isa=-1
*    npar:     number of symmetry doublets included (npar=2 will ensure
*              both lambda doublets)
*  variable in common /cosys/
*    scod:    character*8 array of dimension nscode, which contains names
*             of all system dependent parameters
*  variable in common block /coskip/
*   nskip  for a homonuclear molecule lamda is running in steps of nskip=2
*          for a heteronuclear molecule nskip=1
*
*   skip   same as nskip, used for consistency check
*
*  subroutines called: loapot(iunit,filnam)
      implicit double precision (a-h,o-z)
      include "common/parsys"
      logical readpt, existf
      integer irpot
      character*1 dot
      character*8 scod
      character*(*) fname
      character*40 filnam, line, potfil
      include "common/parbas"
      common /cosys/ scod(lencod)
      common /cosysi/ nscode, isicod, nterm, jmax, igu, isa, npar
      common /cosysr/ isrcod, junkr, brot, aso, p, q
      logical         airyfl, airypr, bastst, batch, chlist, csflag,
     :                flaghf, flagsu, ihomo,lpar
      common /colpar/ airyfl, airypr, bastst, batch, chlist, csflag,
     :                flaghf, flagsu, ihomo,lpar(18)
      common /conlam/ nlam
      common /coskip/ nskip,iskip
      save potfil
      include "common/comdot"
      isicod = 5
      isrcod = 4
      nscode = isicod+isrcod
      scod(1) = 'NTERM'
      scod(2) = 'JMAX'
      scod(3) = 'IGU'
      scod(4) = 'ISA'
      scod(5) = 'NPAR'
      scod(6) = 'BROT'
      scod(7) = 'ASO'
      scod(8) = 'P'
      scod(9) = 'Q'
      potfil = ' '
*  set default values for doublet delta scattering
      nterm = 2
      mproj(2) = 2
      if (iread .eq. 0) then
        mproj(1) = 0
        lammin(1)= 1
        lammax(1) = 1
        lammin(2) = 2
        lammax(2) = 2
        niout=4
        indout(1)=-1
        indout(2)=1
        indout(3)=-2
        indout(4)=2
        npar=2
        igu=1
        isa=0
        jmax=3
      endif
      if(iread.eq.0) return
*  line 13
      read (8, *, err=888) jmax, igu, isa, npar
*  line 14
      read (8, *, err=888) brot, aso, p, q
*  line 15 name of file containing potential parameters
      line=' '
      if(.not.readpt.or.iread.eq.0) then
        call loapot(1,' ')
        return
      endif
      read (8, 285, end=286) line
285   format (a)
      goto 286
* here if read error occurs
888   write(6,1000)
1000  format(/'   *** ERROR DURING READ FROM INPUT FILE ***')
      return
      entry ptr2del (fname,readpt)
      line = fname
      readpt = .true.
286   if (readpt) then
        l=1
        call parse(line,l,filnam,lc)
        if(lc.eq.0) then
          write(6,1020)
1020      format(' FILENAME MISSING FOR POTENTIAL INPUT')
        end if
        j=index(filnam(1:lc),dot)
        if(j.eq.0) then
          call gennam(potfil,filnam,0,'BIN',lc)
          filnam = potfil
        end if
        potfil=filnam
        inquire(file=filnam,exist=existf)
        if(.not.existf) then
          write(6,1025) filnam(1:lc)
1025      format(' FILE ',(a),' NOT FOUND')
          return
        end if
* now call loapot(iunit,filnam) routine to read potential parameters
        call loapot(1,filnam)
        end if
      close (8)
      irpot=1
      return
      entry sav2del (readpt)
*  save input parameters for doublet-delta + atom scattering
*  line 13:
      write (8, 315) jmax, igu, isa, npar
315   format(4i4,18x,'jmax, igu, isa, npar')
*  line 14
      write (8, 320) brot, aso, p, q
320   format(f10.5,f10.4,2(1pg11.3),' brot, aso, p, q')
*  line 15
      write (8, 285) potfil
      return
      end
