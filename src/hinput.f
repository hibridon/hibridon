cstart unix-ibm unix-aix
c;@process noopt
cend
cstart unix-hp
c;!$hp$optimize off
cend
      subroutine hinput(first)
*  subroutine to redefine system independent input parameters for
*  hibridon code
*  author:  h.-j. werner
*  revision: 13-nov-1996 by mby to include bound states
*  revised further:  24-feb-1998 by mha
*  maxbas extended to 28:  30-jul-2018
*  jobname limited to 8 characters:   04-oct-2001, format corrected 2-dec-2007
*  added command sysconf 24-dec-2007
*  added command to compute cross sections for mixed singlet-triplet states
*    e.g. CH3 X3B1 - a1A1 collisions involving perturbed levels - p. dagdigian
*  altered input list for hypxsc:  5-jan-2011 by p. dagdigian
*  added bastp1 basis subroutine (for symmetric top without inversion doubling:
*     p. dagdigian
*  added calculation of transport cross sections - p. dagdigian
*  added 2sigma-2pi (no perts.) basis routine - p. dagdigian
*  last line "quit/exit" now optional in a .com file - q. ma
*  allow calculation of differential cross sections for ibasty=9 - p.dagdigian
*  eliminate special call for ibasty=4 for printc command - p. dagdigian
*
*  revision: 15-mar-2012 by p. dagdigian
*  revision: 20-apr-2012 by q. ma (new eadiab command)
*  revision: 8-oct-2012 by q. ma (allow calculation for user-defined ibasty 100,
*     whose channels have j12's)
*  revision: 27-apr-2013 by p. dagdigian (add prsbr [pressure broadening
*     cross sections] command)
*  revision: 10-jul-2013 by q. ma (new 2pi--1sigma basis)
*  revision: 18-jul-2013 by q. ma (new symtop--1sigma basis)
*  revision: 18-sep-2014 by p. dagdigian (new 3P atom + 2S atom basis)
*  revision: 21-jul-2015 by p. dagdigian (new spherical top-atom basis)
*  revision: 10-may-2107 by p. dagdigian (hypxsc for mol-mol collisions implemented)
*  revision:  1-jun-2017 by p. dagdigian (new 1sig + 1sig basis)
*  revision:  8-jun-2017 by p. dagdigian (new 2sig + 1sig basis)
*  revision: 19-sep-2017 by p. dagdigian (new C2V asym top basis)
*  revision: 30-jul-2018 by p. dagdigian (new 3sig + 1sig basis)
*  revision: 16-jan-2019 by p. dagdigian (new chiral astop + atom basis)
*  revision: 20-jun-2019 by p. dagdigian (new C2V asym - diat molecule basis)
*
*  current revision: 16-jan-2019 by q. dagdigian
* ---------------------------------------------------------------------
      use mod_cosout, only: nnout, jout
      use mod_coiout, only: niout, indout
      use mod_codim, only: nmax => mmax
      implicit double precision (a-h,o-z)
*  iicode is the number of integer pcod's
*  ircode is the number of real pcod's
*  ncode is the number of bcod's
*  bcod stores hibridon's commands
*  fcod stores logical flags
      parameter (ncode = 39, lcode = 28, iicode = 10, ircode = 9,
     :           icode = iicode+ircode)
      character*80 line
      character*40 fnam1,fnam2,jobnam,input,output,savfil,
     :             code
      character*8 bcod(ncode)
      character*8 fcod(lcode),pcod(icode),scod,bascod
      character*9 basknd(30)
* dimension of codex, ihold, lhold, should be equal to largest number
* of identical strings of 1:nnn characters in names of all variables
* (probably 'p' is the most recurring string:  12 times in
*  pcod, fcod, and bcod)
      character*8 codex(15)
      integer ixpar
      integer ipar
      integer ibasty
      logical existf, lspar, first, openfl, is_twomol
      logical lpar, logp, opti, optifl, batch, jtrunc
      dimension a(15),ia(10), ihold(15), lhold(15),lindx(28)
      include "common/parbas"
      include "common/parpot"
      common /cosavi/ iipar, ixpar(iicode)
      common /cosavr/ irpar, junks, rxpar(ircode)
      common /coener/ energ(1)
      common /cofile/ input, output, jobnam, savfil
      common /cosys/  scod(100)
      common /cokeyl/ nncode, llcode, ijcode
      common /cobcod/ bcod
      common /cofcod/ fcod
      common /copcod/ pcod
      common /cosysi/ nscode, isicod, ispar(100)
      common /cosysr/ isrcod, idum, rspar(100)
      common /cosysl/ islcod, lspar(100)
      common /coipar/ ipar(10)
      common /corpar/ rpar(9)
      common /colpar/ lpar(28)
      common /coselb/ ibasty
      common /cobaco/ bascod
      common /coopti/ optifl
      common /coamat/ scmat(1)
      common /conlam/ nlam, nlammx, lamnum(1)
      common /cotwo/ numj,nj1j2(5)
* when adding bases, change size of array basknd and size of
* parameter kmxbas in himain.f
      data basknd /'1-SIGMA', '2-SIGMA', '2-PI', 'SIGMA|PI',
     :              'GEN-PI', 'SYM-TOP-I', '1/3-P-AT', '1SIG+1SIG',
     :              'SYMT-LIN',' 2/2-P-AT', '1-DELTA', 'HOMO+2P',
     :              'HOMO+3P', '2-DELTA', '2P-DIAT', 'ASYM-TOP',
     :              'CH2X', 'SYM-TOP-N', '2SG|2PI_N', '2PI|1SG',
     :              'SYMT|1SG','1D-3P-AT','3P-2S-AT', 'SPH-TOP',
     :              '1SG-1SG', '2SG-1SG', 'C2v-ASTP','3SG-1SG',
     :              'CASYMTOP', 'ASYM-DIAT' /
*  lindx is pointer from fcod order to order in common colpar
      data lindx/1,3,4,5,6,7,8,9,10,11,13,25,26,2,12,14,17,19,24,
     :           15,16,18,20,27,22,21,23,28/
c     common /colpar/ airyfl, airypr, bastst, batch, chlist,
c    :                csflag, flaghf, flagsu, ihomo, ipos, logdfl,
c    :                logwr, noprin, partw, readpt, rsflag, swrit,
c    :                t2test, t2writ, twomol, writs, wrpart, wrxsec,
c    :                xsecwr, nucros, photof, wavefl, boundc
      data irpot, irinp /0, 0/
      data batch /.false./
      save ipr, opti, a, a1, acc, acclas, optval, optacc, istep, inam,
     :     fnam1, fnam2, code, lc, jtot2x, irpot, irinp
      bascod='BASISTYP'
      pcod(1)='JTOT1'
      pcod(2)='JTOT2'
      pcod(3)='JTOTD'
      pcod(4)='JLPAR'
      pcod(5)='NERG'
      pcod(6)='NUMAX'
      pcod(7)='NUMIN'
      pcod(8)='NUD'
      pcod(9)='LSCREEN'
      pcod(10)='IPRINT'
      pcod(11)='FSTFAC'
      pcod(12)='RINCR'
      pcod(13)='RCUT'
      pcod(14)='RENDAI'
      pcod(15)='RENDLD'
      pcod(16)='RSTART'
      pcod(17)='SPAC'
      pcod(18)='TOLAI'
      pcod(19)='XMU'
      fcod(1)='AIRYFL'
      fcod(2)='BASTST'
      fcod(3)='BATCH'
      fcod(4)='CHLIST'
      fcod(5)='CSFLAG'
      fcod(6)='FLAGHF'
      fcod(7)='FLAGSU'
      fcod(8)='IHOMO'
      fcod(9)='IPOS'
      fcod(10)='LOGDFL'
      fcod(11)='NOPRIN'
      fcod(12)='NUCROS'
      fcod(13)='PHOTOF'
      fcod(14)='PRAIRY'
      fcod(15)='PRLOGD'
      fcod(16)='PRPART'
      fcod(17)='PRSMAT'
      fcod(18)='PRT2'
      fcod(19)='PRXSEC'
      fcod(20)='READPT'
      fcod(21)='RSFLAG'
      fcod(22)='T2TEST'
      fcod(23)='TWOMOL'
      fcod(24)='WAVEFL'
      fcod(25)='WRPART'
      fcod(26)='WRSMAT'
      fcod(27)='WRXSEC'
      fcod(28)='BOUNDC'
* addresses for commands
* check: 2700
* debrogli: 1800
* differ: 1500
* difcrs: 2000
* energy: 300
* exit: 600
* help: 75
* input: 900
* intcrs: 2200
* job: 900
* jout: 400
* label: 900
* minpot: 1700
* mrcrs: 2500
* nnout: 2400
* optimimize: 2100
* output: 900
* pot: 1000
* printc: 2600
* prints: 1900
* psi: 2800
* quit: 600
* read: 800
* run: 500
* save: 1300
* show: 700
* tenxsc: 2300
* testpot: 1200
* turn: 1600
* indout: 430
* partc: 2650
* flux: 2800
* eadiab: 2850
* j1j2:  460
* sysconf:  2900
* hypxsc: 2950
* stmix:  3000
* trnprt:  3100
* prsbr:   3200
* nb after changing the following list, check that all the variables "incode"
* that follow after address 900 are changed accordingly
      bcod(1)='CHECK'
      bcod(2)='DEBROGLI'
      bcod(3)='DIFFER'
      bcod(4)='DIFCRS'
      bcod(5)='ENERGY'
      bcod(6)='EXIT'
      bcod(7)='HELP'
      bcod(8)='INPUT'
      bcod(9)='INTCRS'
      bcod(10)='JOB'
      bcod(11)='JOUT'
      bcod(12)='LABEL'
      bcod(13)='MINPOT'
      bcod(14)='MRCRS'
      bcod(15)='NNOUT'
      bcod(16)='OPTIMIZE'
      bcod(17)='OUTPUT'
      bcod(18)='POT'
      bcod(19)='PRINTC'
      bcod(20)='PRINTS'
      bcod(21)='PSI'
      bcod(22)='QUIT'
      bcod(23)='READ'
      bcod(24)='RUN'
      bcod(25)='SAVE'
      bcod(26)='SHOW'
      bcod(27)='TENXSC'
      bcod(28)='TESTPOT'
      bcod(29)='TURN'
      bcod(30)='INDOUT'
      bcod(31)='PARTC'
      bcod(32)='FLUX'
      bcod(33)='J1J2'
      bcod(34)='EADIAB'
      bcod(35)='SYSCONF'
      bcod(36)='HYPXSC'
      bcod(37)='STMIX'
      bcod(38)='TRNPRT'
      bcod(39)='PRSBR'
*
      iipar=iicode
      irpar=ircode
      nncode=ncode
      llcode=lcode
      ijcode=icode
*   define system dependent parameter codes
      if(first) then
         islcod=0
         isrcod=0
         isicod=0
         izero=0
         call sysdat(irpot, lpar(15), izero)
         first = .false.
         call version(6)
*  in this next statement the $ sign implies no line feed
*  replace this with an equivalent formatting character if your system
*  doesn't accept this extension
         go to 1
      else
        do 3 i = 1, ircode
3       rpar(i) = rxpar(i)
        do 4 i = 1, iicode
4       ipar(i) = ixpar(i)
        if(opti) goto 2160
      end if
1     if (.not. lpar(4) .and. .not. batch) write (6, 2)
      optifl = .false.
*  in this next statement the $ sign implies no line feed
*  replace this with an equivalent formatting character if your system
*  doesn't accept this extension
cstart unix mac
2     format(' Hibridon> ',$)
cend
cstart cray
c;2      format(' Hibridon> ')
cend
      call pcoder(lpar(28),pcod,icode)
      read(5, 10, end=599) line
10    format((a))
11    if(line .eq. ' ') goto 1
      if (line(1:1) .eq. '?') then
          code='help '//line(2:)
          call vaxhlp(code)
c         call helppr(line)
          goto 1
      else if (line (1:4).eq.'help' .or. line (1:4).eq.'HELP') then
          call vaxhlp(line)
c         line = '?intro'
c         call helppr(line)
          goto 1
      else if (line(1:3) .eq.'BAT' .or. line(1:3) .eq. 'bat' .or.
     :         line(1:4) .eq.' BAT' .or. line(1:4) .eq.' bat') then
          lpar(4)=.true.
          batch = .true.
          go to 1
      end if
      call upper(line)
      l1 = 1
15    if(l1 .eq. 0) goto 1
      l = -iabs(l1)
      call parse(line,l,code,lc)
      if(lc .eq. 0) goto 1
      match = 0
      do 20 i = 1,ncode
      len = lenstr(bcod(i))
      if(bcod(i)(1:lc) .eq. code(1:lc)) then
        if (lc .eq. len) go to 40
        match = match + 1
        lhold(match) = l
        iskip = 1
        ihold(match) = i
        codex(match) = bcod(i)
      end if
20    continue
      l = iabs(l1)
      do 21 i = 1,icode
      len = index(pcod(i),' ') - 1
      if(pcod(i)(1:lc) .eq. code(1:lc)) then
        if (lc .eq. len) go to 100
        match = match + 1
        lhold(match) = l
        iskip = 2
        ihold(match) = i
        codex(match) = pcod(i)
      end if
21    continue
      do 22 i = 1,lcode
      len=lenstr(fcod(i))
      if(fcod(i)(1:lc) .eq. code(1:lc)) then
        if (lc .eq. len) go to 200
        match = match + 1
        lhold(match) = l
        iskip = 3
        ihold(match) = i
        codex(match) = fcod(i)
      end if
22    continue
      do 23 i = 1,nscode
      len = index(scod(i),' ') - 1
      if(scod(i)(1:lc) .eq. code(1:lc)) then
        if (lc .eq. len) go to 1400
        match = match + 1
        lhold(match) = l
        iskip = 4
        ihold(match) = i
        codex(match) = scod(i)
      end if
23    continue
      len = 8
      if(bascod(1:lc) .eq. code(1:lc)) then
        if (lc .eq. len) go to 50
        match = match + 1
        lhold(match) = l
        iskip = 5
        ihold(match) = i
        codex(match) = bascod
      end if
      if (match .eq. 0) then
        write(6, 27) code(1:lc),(bcod(j),j = 1,ncode)
27      format
     :   (/' *** invalid keyword "',(a),'"; valid request keys are:'//
     1   (1x,6(a8,5x)))
        write(6,31) bascod, (pcod(j),j = 1,icode)
31      format (/(1x,6(a8,5x)))
        write(6,31) (fcod(j),j = 1,lcode)
        write(6,31) (scod(j),j = 1,nscode)
        goto 1
      else if (match .gt. 1) then
        write (6, 28) code(1:lc), (codex(i), i = 1, match)
28      format (' *** ambiguity between input string ',(a),
     :           ' and request keys:',/,7(a10) )
        goto 1
      else
        i = ihold(1)
        l = lhold(1)
        goto (40, 100, 200, 1400, 50), iskip
      end if
40    goto (2700,
     :      1800,1500,2000,300,600,
     :      75,
     :      900,2200,900,400,900,1700,
     :      2500,
     :      2400,2100,900,1000,2600,
     :      1900,2800,600,
     :      800,500,1300,700,2300,
     :      1200,1600,430,2650,2800,
     :      460,2850,2900,2950,3000,
     :      3100,3200),i
* basis type and kind of calculation
50    if(l.eq.0) goto 1
      l1 = l
      call parse(line,l,code,lc)
      call getval(code(1:lc),bascod,icode,j,val)
      if(j .eq. 0) goto 15
      if(j .lt. 0) goto 1
      ibasty=int(val)
      call baschk(ibasty)
* set twomolecule true
      if (is_twomol(ibasty)) then
        lpar(20)=.true.
      else
        lpar(20)=.false.
      endif
      call sysdat(irpot, lpar(15), izero)
      l1=l
      goto 15
*  request help
c75    line = '?intro'
c      call helppr(line)
 75   call vaxhlp(line)
      goto 1
*  parameters
*  specify parameters in the form cod1=val1, cod2=val2, etc.
100   if(l .eq. 0) goto 1
      call pcoder(lpar(28),pcod,icode)
      l1 = l
      call parse(line,l,code,lc)
      call getval(code(1:lc),pcod,icode,j,val)
      if(j .eq. 0) goto 15
      if(j .lt. 0) goto 1
      if (j .eq. 5) then
        if (ipar(j) .lt. val) write (6, 101)
101     format
     :  (1x,'** NERG INCREASED; VERIFY ARRAY OF COLLISION ENERGIES')
      end if
      if(j.le.iicode) then
        ipar(j) = val
        if(ipar(8).ne.1) lpar(25)=.true.
      else
        rpar(j-iicode) = val
      end if
      call enord(energ,ipar(5))
*  numin and numax should be 0 if cc calculation, if not, then set them
*  equal to zero
* NB this is disabled for basisknd=12 (2P atom + homonuclear)
      if (.not. lpar(6)) then
        lpar(25)=.false.
        if (ipar(6) .ne. 0) then
          write (6, 105)
105       format ('  CC calculation, numax set to zero')
          ipar(6) = 0
        end if
        if (ipar(7) .ne. 0..and.ibasty.ne.12) then
          write (6, 106)
106       format ('  CC calculation, numin set to zero')
          ipar(7) = 0
        end if
      end if
      goto 100
* flags
* specify flags in the form cod1=val1,cod2=val2, etc
* where val(i) must be either "t(rue)" or "f(alse)"
200   if(l .eq. 0) goto 1
      l1 = l
      call parse(line,l,code,lc)
      call getval(code(1:lc),fcod,lcode,j,val)
      if(j .eq. 0) goto 15
      if(j .lt.0) goto 1
      logp = .false.
      if(val .eq. 1) logp = .true.
      if (j .eq. 3) then
        write (6, 201)
201     format (' ** BATCH FLAG CAN NOT BE SET INTERACTIVELY!')
        lpar(4) = batch
        go to 1
      end if
      lpar(lindx(j)) = logp
      goto 200
* energies
* specify energies in the form
* energ=e1,e2,e3...
* terminate the string with a semicolon if other parameters will follow
* on the same card, e.g. energ=e1,e2,e3;jtot1=0,jtot2=2....
300   i = 0
310   if(l .eq. 0) goto 320
      if(line(l-1:l-1) .eq. ';') goto 320
      i = i+1
      call parse(line,l,code,lc)
      call getval(code(1:lc),' ',0,j,energ(i))
      goto 310
320   if (energ(1) .gt. 0) then
         if (i .ne. nerg) then
           if (i. le. 25) then
             write (6, 321) i
321          format (' ** NERG HAS BEEN RESET TO',i3)
             nerg = i
           else
             write (6, 322)
322          format (' ** NERG RESET TO 25 (MAXIMUM VALUE)')
             nerg=25
           end if
         end if
      elseif (energ(1) .lt. 0) then
         nerg=energ(4)+0.001d0
         if (nerg.gt.25) then
            nerg=25
            write (6,322)
         endif
         e1=energ(2)
         delt_e=(energ(3)-e1)/(nerg-1)
         nde=delt_e*100d0
         delt_e=nde/100d0
         write (6,323) e1,nerg, delt_e,e1+(nerg-1)*delt_e
323      format(' ** GRID OF ENERGIES:  E-FIRST = ', f9.2,
     :     '; NERG = ',i2, '; DELTA_E = ',f7.2,
     :        '; E-LAST = ',f9.2)
         do ii=1,nerg
            energ(ii)=e1+(ii-1)*delt_e
         enddo
      endif

      call enord(energ,nerg)
      ipar(5) = nerg
      l1 = l
      goto 15
* jout values
* specify jout values in the form
* jout,nnout,jout(1),...,jout(iabs(nnout))
* terminate the string with a semicolon if other parameters will follow
* on the same card, e.g. jout,-3,0,2,4;energ=e1,e2,e3;jtot1=0,jtot2=2....
400   i = 0
      call parse(line,l,code,lc)
      call getval(code(1:lc),' ',0,j,val)
      l1 = l
      nnout=val
410   if(l .eq. 0) goto 420
      if(line(l-1:l-1) .eq. ';') goto 420
      call parse(line,l,code,lc)
      call getval(code(1:lc),' ',0,j,val)
      i = i+1
      jout(i) = val
      goto 410
420   if(nnout.ge.0) nnout = i
      if(nnout.lt.0) nnout = -i
      l1 = l
      goto 15
* indout values
* specify indout values in the form
* indout,niout,indout(1),...,indout(niout)
* terminate the string with a semicolon if other parameters will follow
* on the same card, e.g. indout,2,1,-1;energ=e1,e2,e3;jtot1=0,jtot2=2....
430   i = 0
      call parse(line,l,code,lc)
      call getval(code(1:lc),' ',0,j,val)
      l1 = l
      niout=val
      if(niout.eq.0) goto 15
440   if(l .eq. 0) goto 450
      if(line(l-1:l-1) .eq. ';') goto 450
      call parse(line,l,code,lc)
      call getval(code(1:lc),' ',0,j,val)
      i = i+1
      indout(i) = val
      goto 440
450   if(niout.ge.0) niout = i
      l1 = l
      goto 15
* j1j2 values
* specify j1j2 values in the form
* j1j2,numj,j1j2(1),...,j1j2(numj)
* terminate the string with a semicolon if other parameters will follow
* on the same card, e.g. j1j2,2,00,10;energ=e1,e2,e3;jtot1=0,jtot2=2....
460   if (.not.lpar(20)) then
        write (6, 465)
465     format(' ** NUMJ CAN ONLY BE DEFINED IF TWOMOL = .TRUE.')
        goto 15
      endif
      i = 0
      call parse(line,l,code,lc)
      call getval(code(1:lc),' ',0,j,val)
      l1 = l
      numj=val
      if(niout.eq.0) goto 15
470   if(l .eq. 0) goto 480
      if(line(l-1:l-1) .eq. ';') goto 480
      call parse(line,l,code,lc)
      call getval(code(1:lc),' ',0,j,val)
      i = i+1
      nj1j2(i) = val
      goto 470
480   if(numj.ge.0) then
        numj = i
      else
        write (6, 485)
485     format(' ** YOU MUST SPECIFY A VALUE OF NUMJ')
      endif
      l1 = l
      goto 15
* start execution,run
*  numin and numax should be 0 if cc calculation, if not, then set them
*  equal to zero
500   if (.not. lpar(6)) then
        if (ipar(6) .ne. 0) then
          write (6, 105)
          ipar(6) = 0
        end if
* NB this is disabled currently for 2P atom + homonuclear
        if (ipar(7) .ne. 0.and.ibasty.ne.12) then
          write (6, 106)
          ipar(7) = 0
        end if
      end if
      call pcoder(lpar(28),pcod,icode)
      if (lpar(6).and.ipar(8).ne.1) lpar(25)=.true.
      nerg=ipar(5)
* check to see if flags are ok if wavefunction desired or
* photodissociation calculation
      call genchk
      call enord(energ,nerg)
      do 503 i = 1,ircode
503   rxpar(i) = rpar(i)
      do 504 i = 1,iicode
504   ixpar(i) = ipar(i)
      if(irinp.eq.0) then
        write(6,505)
505     format (/,' ** SAVE DEFAULT VARIABLES OR SPECIFY INPUT',
     :          ' FILE WITH INP = filename')
        if(lpar(4)) call exit
        goto 1
      end if
      if(rpar(8).eq.0) then
        write(6,507)
507     format(/,' ** SPECIFY COLLISION REDUCED MASS WITH XMU = mass')
        goto 1
      end if
      if(irpot.ne.0.or..not.lpar(15)) then
* open output file
* first make sure it is lower case
        call lower(output)
        call upper(output(1:1))
        call openf(9,output,'sf',0)
*     write input data to file outpt
        write (9, 508) label
508     format (1x,a)
        write (9, 508) potnam
        write (9, 240)
240     format(1h ,30('='))
        write(9,710) 'Parameters:',(pcod(j),ipar(j),j = 1,iicode)
        write(9,720) (pcod(iicode+j),rpar(j),j = 1,ircode)
        if (ibasty .ne. 99) then
          length = index(basknd(ibasty),' ') - 1
          if (length .eq. -1) length=8
cstart unix cray mac
          write(9,710) basknd(ibasty)(1:length)//' system parameters:',
     :             (scod(j),ispar(j),j = 1,isicod)
cend
        else
          write(9,710) 'User defined system parameters:',
     :               (scod(j),ispar(j),j = 1,isicod)
        endif
        if(isrcod.gt.0)
     :    write(9,720) (scod(isicod+j),rspar(j),j = 1,isrcod)
        if(islcod.gt.0)
     :    write(9,735) (scod(isicod+isrcod+j),lspar(j),j = 1,islcod)
        if (.not. lpar(20) ) then
           write(9,736) 'LAMMIN: ',(lammin(j),j=1,ispar(1))
           write(9,736) 'LAMMAX: ',(lammax(j),j=1,ispar(1))
           write(9,736) 'MPROJ:  ',(mproj(j),j=1,ispar(1))
        else if (lpar(20)) then
           write (9, 738)'J1/J2: ',(nj1j2(j)/10,mod(nj1j2(j),10),
     :                    j=1,numj)
        endif
736     format(1x,(a),10i4,/,9x,10i4)
737     format(1x,(a),i2,(a),20(i3,1x))
738     format (1x,(a),1x,20(2i1,'  ') )
739     format (1x,(a),i2,(a),20(2i1,'  ') )
        write(9,730) 'Flags:',(fcod(j),lpar(lindx(j)),j = 1,lcode)
        call enord(energ,ipar(5))
        if(ipar(5).gt.0) write(9,740) (energ(j),j = 1,ipar(5))
        if(nnout.ne.0) then
          if (.not.lpar(20) ) then
            write (9,737) 'NOUT: ',nnout,
     :           '; JOUT:',(jout(j), j=1,iabs(nnout) )
          else
            write (9,739) 'NOUT: ',nnout,
     :         '; J1/J2-OUT: ',
     :        (jout(j)/10, mod (jout(j),10), j = 1,iabs(nnout))
          end if
        end if
        if(niout.ne.0) write (9,736) 'INDOUT: ',(indout(j), j=1,niout)
        write (9, 240)
        return
      else
        write(6,510)
510     format(' Potential not yet defined!')
        goto 1
      end if
* exit
599   write (6, *)
600   call exit
* show all parameters and flags
* show
700   l1 = l
      call pcoder(lpar(28),pcod,icode)
      call parse(line,l,code,lc)
      if (.not.lpar(28)) then
        write(6,710)
     :   'Parameters (scattering):',(pcod(j),ipar(j),j = 1,iicode)
      else
        write(6,710)
     :   'Parameters (bound-state):',(pcod(j),ipar(j),j = 1,iicode)
      endif
      write(6,720)   (pcod(iicode+j),rpar(j),j = 1,ircode-1)
      write(6,1720)  pcod(iicode+ircode), rpar(ircode)
      if(nnout.ne.0) then
        if (.not.lpar(20) ) then
          write (6,737) 'NOUT: ',nnout,
     :           '; JOUT:',(jout(j), j=1,iabs(nnout) )
        else
          write (6,739) 'NOUT: ',nnout,
     :         '; J1/J2-OUT: ',
     :        (jout(j)/10, mod (jout(j),10), j = 1,iabs(nnout))
        end if
      end if
      if(niout.ne.0) write (6,701) 'INDOUT:',(indout(j), j=1,niout)
701     format(1x,(a),10i5,/,8x,10i5,/,5x,10i5)
      if (ibasty .lt. 99) then
        length = index(basknd(ibasty),' ') - 1
        if (length .eq. -1) length=9
cstart unix cray mac
          write(6,710) basknd(ibasty)(1:length)//' system parameters:',
     :             (scod(j),ispar(j),j = 1,isicod)
cend
      else
        write(6,710) 'user defined system parameters:',
     :               (scod(j),ispar(j),j = 1,isicod)
      endif
      if(isrcod.gt.0)
     :  write(6,720) (scod(isicod+j),rspar(j),j = 1,isrcod)
      if(islcod.gt.0)
     :  write(6,735) (scod(isicod+isrcod+j),lspar(j),j = 1,islcod)
      if (.not. lpar(20) ) then
        write(6,736) 'LAMMIN: ',(lammin(j),j=1,ispar(1))
        write(6,736) 'LAMMAX: ',(lammax(j),j=1,ispar(1))
        write(6,736) 'MPROJ:  ',(mproj(j),j=1,ispar(1))
      else if (lpar(20)) then
        write (6, 738)'J1/J2: ',(nj1j2(j)/10,mod(nj1j2(j),10),
     :                    j=1,numj)
      end if
      write(6,730) 'Flags:',(fcod(j),lpar(lindx(j)),j = 1,lcode)
* in the next line ipar(9) is lscreen
      write(6,731)  nmax, nlammx
      if (ipar(9) .le. 24 .and. .not. batch) then
        write (6, 703)
703     format (6x,'enter <return> to continue,',
     :             ' <q> for prompt, or new data')
        read (5, 10) line
        if (line(1:1) .eq. 'q' .or. line(1:1) .eq. 'q') then
          go to 1
        else if (line(1:1) .ne. ' ') then
          go to 11
        end if
      end if
      call enord(energ,ipar(5))
      if(ipar(5).gt.0) write(6,740) (energ(j),j = 1,ipar(5))
710   format(5x,'*** ',(a)/(4(1x,a7,'=',i4,7x)))
720   format(4(1x,a7,'=',1pg11.4))
1720  format(1x,a7,'=',f10.5)
730   format(5x,'*** ',(a)/(6(1x,a6,'=',l2,3x)))
731   format(1x,'** Maximum Channels: ', i4, '; ',
     :  'Maximum Anisotropic Terms: ',i5)
735   format(3(1x,a7,'=',l2,9x))
740   format(1x,'** Energies:',(t15,5f15.6))
      lenout = index(output,' ')-1
      lenjob = index(jobnam,' ')-1
      jtrunc = .false.
      if (lenjob.gt.8) then
         jobnam=jobnam(1:8)
         jtrunc = .true.
      endif
      leninp=index(input,' ')-1
      write (6, 751) label
751   format (1x,'** Label:      ',(a))
      write (6, 752) potnam
752   format (1x,'** Pot name:      ',(a))
      if (.not. jtrunc) then
         write (6, 753) input(1:leninp),
     :     output(1:lenout), jobnam(1:lenjob)
753   format(1x,'** Input File:  ',(a),/
     :       1x,'** Output file: ',(a),/,1x,'** Jobname:     ',(a))
      else
         write (6, 754) input(1:leninp),
     :     output(1:lenout), jobnam(1:lenjob)
754   format(1x,'** Input File:  ',(a),/
     :       1x,'** Output file: ',(a),/,1x,'** Jobname:     ',(a),
     :       ' (** TRUNCATED TO 8 CHARACTERS **)')
      endif
      l1 = l
      goto 15
* read
800   call pcoder(lpar(28),pcod,icode)
      call gendat
      ione=1
      call sysdat(irpot, lpar(15),ione)
      irinp=1
      l1 = l
      if (batch) lpar(4) = .true.
      goto 15
* input, output, label and job file names
* input=infile, output=outfile, job=jobfile
* input, output, and label are now lower case:  mha 6.6.91
900   call parse(line,l,code,lc)
* incode is index of command input in list bcod
      incode=8
      if(i .eq. 8) then
        input = code(1:lc)
        call lower(input)
        call upper(input(1:1))
        inquire (file=input, exist=existf)
        if (.not. existf) then
          len = index (input,' ')
          write (6, 901) input(1:len)
901       format (' *** INPUT FILE ',(a),' DOES NOT EXIST')
          if(batch) call exit
          go to 1
        end if
        goto 800
* incode 10 is index of command job in list bcod
      else if(i .eq. 10) then
        jobnam = code(1:lc)
        call lower(jobnam)
        call upper(jobnam(1:1))
* incode 12 is index of command label in list bcod
      else if(i .eq. 12) then
        low = index (line, '=') + 1
        lend = index (line, ';') - 1
        if(lend.lt.0) then
          lend=lenstr(line)
          l=0
        else
          l=lend+2
        end if
        label = line(low:lend)
* incode 17 is index of command output in list bcod
      else if(i .eq. 17) then
        output = code(1:lc)
        call lower(output)
        call upper(output(1:1))
        inquire(9,opened=openfl)
        if(openfl) then
          endfile(9)
          close(9)
        end if
        call openf(9, output, 'sf', 0)
      end if
      l1 = l
      goto 15
* read parameters for potential
*     pot=potfile
1000  call parse(line,l,code,lc)
      call ptread(code(1:lc),lpar(15))
      l1 = l
      goto 15
* test potential
* testpot
* you will be prompted for r and theta. to exit, specify r=0
1200  if (ibasty .eq. 1 .or. ibasty .eq. 4) then
        call testptn(lpar(9))
      else if (ibasty .eq. 20) then
        call testpt20(lpar(9))
      else
        call testpt(lpar(9))
      end if
      goto 1
* save input parameters
*     save=filename
*     if filename is not specified, the inputfile is overwritten
1300  inew=0
      call pcoder(lpar(28),pcod,icode)
      if(l.ne.0) then
        call parse(line,l,code,lc)
        if(lc .eq. 0) then
          code = input
        else
          inew=1
* convert input files to lower case
          call lower(code)
          call upper(code(1:1))
        endif
      else
        code = input
      end if
      call savdat(inew,code)
      call syssav(lpar(15))
      close(8)
      l1 = l
      irinp = 1
      goto 15
* redefine system dependent parameters
* specify in the same way as other parameters, e.g.
*     jmin=0,jmax=4,brot=2.2...
1400  if(l.eq.0) goto 1
      l1=l
      call parse(line,l,code,lc)
      call getval(code(1:lc),scod,nscode,j,val)
      if(j .eq. 0) goto 15
      if(j .lt. 0) goto 1
      if(j.eq.1 .and. .not.lpar(20)) then
        write(6,1401) scod(j)
1401    format(1x,a,'CAN NOT BE MODIFIED')
        goto 1400
      end if
      if(j.le.isicod) then
        ispar(j) = val
      else if(j.le.isrcod+isicod) then
        rspar(j-isicod) = val
      else
        if(val .eq. 1) then
          lspar(j-isrcod-isicod) = .true.
        else
          lspar(j-isrcod-isicod) = .false.
        end if
      end if
      goto 1400
*.....differences of s-matrices
*  dif,jobfile1,jobfile2,iprint,ienerg,thrs
*  this compares s-matrices in the files jobfile1.smt and jobfile2.smt
*  if iprint.eq.0 largest average and absolute differences are printed
*  if iprint.eq.1 these values are given for all individual s-matrices
*  if iprint.ge.2 s-matrices are printed
*  thrs: threshold for neglect of small s-matrix elements in comparison
1500  call parse(line,l,code,lc)
      if(code.ne.' ') fnam1 = code
      call lower(fnam1)
      call upper(fnam1(1:1))
      call parse(line,l,code,lc)
      if(code.ne.' ') fnam2 = code
      call lower(fnam2)
      call upper(fnam2(1:1))
      if(fnam1 .eq. ' '.or.fnam2 .eq. ' ') goto 1
      iprint = 0
      ienerg = 1
      thrs = 1.e-5
      if(l.ne.0) then
        call parse(line,l,code,lc)
        call getval(code(1:lc),' ',0,j,val)
        iprint = val
      end if
      if(l.ne.0) then
        call parse(line,l,code,lc)
        call getval(code(1:lc),' ',0,j,val)
        ienerg = val
        ienerg = max0(1,ienerg)
      end if
      if(l.ne.0) then
        call parse(line,l,code,lc)
        call getval(code(1:lc),' ',0,j,val)
        thrs = val
      end if
      call difs(fnam1,fnam2,ienerg,iprint,acc,accmx,thrs,imx,jmx,ityp)
      code = '?'
      lc = 1
      if (thrs .lt. 0.) then
        if(ityp .eq. 1) then
          code = 'S real'
          lc = 6
        end if
        if(ityp .eq. 2) then
          code = 'S imaginary'
          lc = 11
        end if
      else
        code = 'S modulus'
        lcc = 9
      end if
      if(iprint .eq. 0) write(6,1510) acc,accmx,imx,jmx,code(1:lc),
     :  abs(thrs)
1510  format(' Average relative difference:',f10.2,'%'/
     :       ' Largest relative difference:',f10.2,'%'/
     : ' (i = ',i2,'  j = ',i2,') element of ',(a)/
     : ' Inspection threshold is ',1pg8.1)
      goto 1
*.....determine turning point from isotropic potential
*     turn
1600  if(irpot .eq. 0) then
        write(6,510)
        goto 1
      end if
      e = 0
      do 1605 i = 1,ipar(5)
1605  e = max(e,energ(i))
      if(e .eq. 0) then
        write(6,1610)
1610    format(' Total energy has not been given a value !')
        goto 1
      end if
      r = turn(e)
      write(6,1620) r
1620  format(' Turning point for isotropic potential at R = ',
     :         f5.2, ' bohr')
      l1 = l
      goto 15
*  determine minimum of isotropic potential
*  minpot
1700  if(irpot .eq. 0) then
        write(6,510)
        goto 1
      end if
      r = potmin()
      write(6,1710) r
1710  format(' Minimum of isotropic potential at r = ',f5.2, ' bohr')
      l1 = l
      goto 15
*  calculate de broglie wavelength in bohr (defined as 2pi/k)
*  debrogli
1800  e = 0
      do 1810 i = 1,ipar(5)
1810  e = max(e,energ(i))
      if(e .eq. 0) then
        write(6,1610)
        goto 1
      end if
      xmu = rpar(9)
      if(xmu .eq. 0) then
        write(6,1820)
1820    format(' Collision reduced mass has not been given a value !')
        goto 1
      end if
      r = 48.75/sqrt(xmu*e)
      write(6,1830) r,0.2*r
1830  format('   de Broglie wavelength  = ',f6.3,' bohr / 5  = ',f6.3)
      waveve = 6.283/r
      write (6, 1831) waveve, 1.8897*waveve
1831  format ('   wavevector =', g12.5,' Bohr^-1 = ',
     :   g12.5,' Angstroms^-1')
      l1 = l
      goto 15
*.....print s-matrices:
*  print,jobfile,j1,j2,jd,jlp,ienerg
*  first matrix printed is for jtot=j1
*  last matrix printed is for jtot=j2
*  increment of jtot is jd
*  default values are:
*  j1 = jtot1, j2=jtot1, jd=jtotd, jlp=jlpar,ienerg=1
*  jlp: parity; zero means both values
1900  call parse(line,l,fnam1,lc)
      if(fnam1 .eq. ' ') fnam1 = jobnam
      call lower(fnam1)
      call upper(fnam1(1:1))
      do 1910 i = 1,5
      ia(i) = 0
      if(l .eq. 0) goto 1910
      call parse(line,l,code,lc)
      call getval(code(1:lc),' ',0,j,a(i))
      ia(i)=a(i)
1910  continue
      if(ia(1).eq.0) ia(1)=ipar(1)
      if(ia(2).eq.0) ia(2)=ipar(2)
      if(ia(3).eq.0) ia(3)=ipar(3)
      if(ia(4).eq.0) ia(4)=ipar(4)
      call lower(fnam1)
      call upper(fnam1(1:1))
      call sprint(fnam1,ia)
      goto 1
*.....differential cross sections:
*  diffc,jobfile,j1,in1,j2,in2,ang1,ang2,dang,ienerg,jtotend
*
*  differential cross sections computed for atom-molecule
*  collisions or symmetric top-linear molecule collisions
2000  if (.not. lpar(20) .or. is_twomol(ibasty)) then
        call parse(line,l,fnam1,lc)
        if(fnam1 .eq. ' ') fnam1 = jobnam
        call lower(fnam1)
        call upper(fnam1(1:1))
        do 2010 i = 1,15
           a(i) = 0.d0
           if(l .eq. 0) goto 2010
           call parse(line,l,code,lc)
           call getval(code(1:lc),' ',0,j,a(i))
2010    continue
        call difcrs(fnam1,a,lpar(9),lpar(7))
      else
        write (6, 2011)
2011    format(' Sorry, differential cross sections not yet',
     :         /,'  implemented for most molecule-molecule ',
     :         'collisions')
      end if
      goto 1
*.....optimize
*  opt,code,start,end,fak,add,accav,accmx,thrs
*  code is name of variable to be optimized
*  start: start value
*  end:   end value
*  the value is updated in each step according to
*  val(n)=val(n-1)*fak+add
*  if(fak.eq.0) fak=1.0
*  optimization is stopped if the relative difference in the
*  s-matrices of subsequent runs are smaller than accav
*  (average error) and accmx (maximum error)
*  optimization is performed for jtot=jtot1 only
*  thrs:  threshold for check for s-matrix elements. all s-matrix
*  elements which are smaller than thrs are not compared
*  comparison is of s modulus if thrs .ge. 0, otherwise of s matrix
2100  if (ipar(5) .gt. 1) then
* no optimization if more than one energy requested
        write (6, 2101)
2101    format(' ** NERG SET EQUAL TO 1 FOR OPTIMIZATION')
        ipar(5)=1
      endif
      call parse(line,l,code,lc)
      do 2110 ipr = iicode+1,icode
2110  if(code(1:lc) .eq. pcod(ipr)(1:lc)) goto 2120
      write(6,2115) code(1:lc)
2115  format(' Invalid parameter: "',(a),'"')
      goto 1
2120  do 2130 i = 1,7
      a(i) = 0.d0
      if(l .eq. 0) goto 2130
      call parse(line,l,code,ld)
      call getval(code(1:ld),' ',0,j,a(i))
2130  continue
      if(a(1) .eq. 0.or.a(2) .eq. 0) then
        write(6,2140)
2140    format(' Initial or final value of parameter to optimize',
     :         ' has not been defined')
        goto 1
      end if
      if (.not.lpar(21)) then
        write (6,2141)
2141    format (' Flag WRSMAT set to .TRUE. for optimization')
        lpar(21) = .true.
      end if
      if (nnout .lt. 0) then
        write (6, 2142)
2142    format (' NNOUT set positive for optimization')
        nnout = iabs (nnout)
      else if (nnout .eq. 0) then
        write (6, 2143)
2143    format (' NNOUT=0; optimization not possible; reset NNOUT')
        go to 1
      end if
      if (.not.lpar(13)) then
         write (6,2144)
2144     format(' Flag NOPRIN set to .TRUE. for optimization')
         lpar(13)=.true.
      endif
*
      if(a(3) .eq. 0.and.a(4) .eq. 0.and.a(2).gt.a(1)) a(3) = 2.0d0
      if(a(3) .eq. 0.and.a(4) .eq. 0.and.a(2).lt.a(1)) a(3) = 0.5d0
      if(a(3) .eq. 0) a(3) = 1.0d0
      if(a(5) .eq. 0) a(5) = 1.0d0
      if(a(6) .eq. 0) a(6) = 5.0d0
      thrs = 1.d-5
      if(a(7).ne.0) thrs=a(7)
      a1 = a(1)
      rpar(ipr-iicode) = a1
      if(a(1).lt.a(2).and.(a(1)*a(3)+a(4).lt.a(1)).or.
     1   a(1).gt.a(2).and.(a(1)*a(3)+a(4).gt.a(1)).or.
     2   a(1) .eq. a(2)) then
         write(6,2150)
2150     format(' Invalid step parameters for OPTIMIZE')
         goto 1
      end if
      jtot2x = ipar(2)
      write(6,2151) pcod(ipr)(1:lc),ipar(1),(a(i),i = 1,6),
     :              abs(thrs)
2151  format(' Optimization of ',(a),' for Jtot = ',i3,/,
     :  ' Start:',f7.3,'  End:',f7.3,'  Factor:',f5.2,
     :  ' Increment:',f7.3,/,
     :  ' Average error limit:',f4.1,'%',/,
     :  ' Maximum error limit:',f4.1,'%',/,
     :  ' Threshold for S-matrix elements:',e10.2)
      ipar(2) = ipar(1)
      fnam1 = 'Joba'
      fnam2 = 'Jobb'
* avoid too long filenames on cray (cos)
      inam = 1
      istep = 1
      jobnam = fnam1
      acclas = 1.d10
      optval = a1
      opti = .true.
      optifl = .true.
      write(6,255) pcod(ipr)(1:lc),a(1)
255   format(1x,(a),' = ',f7.3)
      goto 500
2160  if(istep.ge.2) then
        call difs(fnam1,fnam2,1,0,acc,accmx,thrs,im,jm,ityp)
        code = '?'
        lcc = 1
        if(ityp .eq. 1) code = 'S real'
        if(ityp .eq. 1) lcc = 6
        if(ityp .eq. 2) code = 'S imaginary'
        if(ityp .eq. 2) lcc = 11
        if (thrs .ge. 0.) then
          code = 'S modulus'
          lcc  = 9
        end if
        write(6,2165) code(1:lcc),
     :         acc,code(1:lcc),accmx,im,jm,code(1:lcc)
2165    format(' average difference between old and new ',(a),
     :    ' = ',f10.2,'%',/,
     1         ' Largest difference between old and new ',(a),
     :    ' = ', f10.2,'%',/,
     :     ' in (i = ',i2,' j = ',i2,') element of ',(a))
        if(acc.lt.acclas) then
          optval = a(1)
          optacc = acc
          optacm = accmx
          imx = im
          jmx = jm
          itx = ityp
        end if
        acclas = acc
      end if
      a(1) = a(1)*a(3)+a(4)
      if((acclas.lt.a(5).and.accmx.lt.a(6))
     :   .or.(a(2).gt.a1.and.a(1).gt.a(2)).or.
     :   (a(2).lt.a1.and.a(1).lt.a(2))) then
         code = '?'
         lcc = 1
         if(itx .eq. 1) code = 'S real'
         if(itx .eq. 1) lcc = 6
         if(itx .eq. 2) code = 'S imaginary'
         if(itx .eq. 2) lcc = 11
          if (thrs .ge. 0.) then
            code = 'S modulus'
            lcc  = 9
          end if
         rpar(ipr-iicode) = optval
         write(6,2170) pcod(ipr)(1:lc),optval,code(1:lcc),
     :      optacc,code(1:lcc),optacm,imx,jmx,code(1:lcc)
2170     format(' optimized value for ',a,' = ',g11.4,/,
     1   ' average difference in old and new ',(a),' is',
     2     f10.2,'%',/,
     3   ' Largest difference in old and new ',(a),' is',
     4     f10.2,'%',/,
     :     ' in (i = ',i2,' j = ',i2,') element of ',(a))
         opti = .false.
         ipar(2) = jtot2x
         goto 1
      end if
      rpar(ipr-iicode) = a(1)
      if(inam .eq. 1) then
        jobnam = fnam2
        inam = 2
      else
        jobnam = fnam1
        inam = 1
      end if
      istep = istep+1
      write(6,255) pcod(ipr)(1:lc),a(1)
      goto 500
*.....integral cross sections
*  intcrs,jobfile,in1,in2,ienerg,maxjtot
2200  call parse(line,l,fnam1,lc)
      if(fnam1 .eq. ' ') fnam1 = jobnam
      call lower(fnam1)
      call upper(fnam1(1:1))
      do 2210 i = 1,4
      a(i) = 0
      if(l .eq. 0) goto 2210
      call parse(line,l,code,lc)
      call getval(code(1:lc),' ',0,j,a(i))
2210  continue
      call intcrs(fnam1,a)
      goto 1
*.....tensor cross sections
*  tenxsc,jobfile,maxn,iframe,in1,in2,ienerg,jtotend,minj,maxj
2300  call parse(line,l,fnam1,lc)
      if(fnam1 .eq. ' ') fnam1 = jobnam
      call lower(fnam1)
      call upper(fnam1(1:1))
      do 2310 i = 1,9
      a(i) = 0.d0
      if(l .eq. 0) goto 2310
      call parse(line,l,code,lc)
      call getval(code(1:lc),' ',0,j,a(i))
2310  continue
cstart mac
c;      call exit
c;      write(6,2320)
c;2320  format(' Sorry, tensor cross sections not yet implemented')
cend
cstart unix
      call tenopa(fnam1,a)
cend
      goto 1
*....nnout must be preceded by jout
2400  write (6, 2410)
2410  format(' To change NNOUT, enter the command line',/,
     :  '    jout,nnout,jout(1),...,jout(iabs(nnout))' )
      goto 1
*.....m-resolved cross sections
*  mrcrs,jobfile,ienerg
2500  call parse(line,l,fnam1,lc)
      if(fnam1 .eq. ' ') fnam1 = jobnam
      call lower(fnam1)
      call upper(fnam1(1:1))
      do 2510 i = 1,1
      a(i) = 0
      if(l .eq. 0) goto 2510
      call parse(line,l,code,lc)
      call getval(code(1:lc),' ',0,j,a(i))
2510  continue
cstart mac
c;      call exit
c;      write(6,2520)
c;2520  format(' M-resolved cross sections not yet implemented')
cend
cstart unix cray
      call mrcrs(fnam1,a)
cend
      goto 1
* print selected integral cross sections from ics file
2600  call parse(line,l,fnam1,lc)
      if(fnam1 .eq. ' ') fnam1 = jobnam
      call lower(fnam1)
      call upper(fnam1(1:1))
      do 2610 i = 1,8
      a(i) = 0
      if(l .eq. 0) goto 2610
      call parse(line,l,code,lc)
      call getval(code(1:lc),' ',0,j,a(i))
2610  continue
      if(ibasty.eq. 5) write (6, 2611)
2611  format
     : (' *** printed triplet pi cross sections should be checked')
      if(ibasty.eq. 6) write (6, 2612)
2612  format
     : (' *** printed symmetric top cross sections should be checked')
*
*  eliminate special call for ibasty = 4
      call prsg(fnam1,a)
*      if(ibasty.ne.4) call prsg(fnam1,a)
*      if(ibasty.eq.4) call prsgpi(fnam1,a)
      goto 1
* print selected partial cross sections from pcs file
2650  call parse(line,l,fnam1,lc)
      if(fnam1 .eq. ' ') fnam1 = jobnam
      call lower(fnam1)
      call upper(fnam1(1:1))
      do 2660 i = 1,8
      a(i) = 0
      if(l .eq. 0) goto 2660
      call parse(line,l,code,lc)
      call getval(code(1:lc),' ',0,j,a(i))
2660  continue
      call readpc(fnam1, a, scmat, nmax)
      goto 1
* check if inconsistencies in input parameters
2700  call genchk
      goto1
*  psi(wavefunction calculation),jobfile,mchannel
*  flux calculation,jobfile,mchannel,iflux,thresh,iprint
2800  call parse(line,l,fnam1,lc)
      call lower(fnam1)
      call upper(fnam1(1:1))
      if(fnam1 .eq. ' ') fnam1 = jobnam
      do 2810 i = 1,10
      a(i) = 0
      if(l .eq. 0) goto 2810
      call parse(line,l,code,lc)
      call getval(code(1:lc),' ',0,j,a(i))
2810  continue
      iflux=a(1)
      if (a(2) .eq. 0.d0) iflux=2
      call psi(fnam1,a)
      goto 1
* adiabatic energy calculation, jobfile
2850  call parse(line,l,fnam1,lc)
      call lower(fnam1)
      call upper(fnam1(1:1))
      if(fnam1 .eq. ' ') fnam1 = jobnam
      call parse(line,l,code,lc)
      if (code .eq. ' ') then
         l1 = 1
         l2 = 10
      else
         read (code, *, err=2860, end=2860) l2
         call parse(line,l,code,lc)
         if (code .eq. ' ') then
            l1 = 1
         else
            l1 = l2
            read (code, *, err=2860, end=2860) l2
         end if
      end if
      call eadiab1(fnam1,l1,l2)
      goto 1
2860  write (6, *) 'Parameters to EADIAB cannot be recognized'
      goto 1
*  print out system parameters
2900  call sys_conf
      goto 1
*  hyperfine xcs routine (originally written by j. klos,
*  rewritten by p.j. dagdigian
*  hypxsc,jobfile, ienerg ,nucspin, j1, j2
2950  continue
*     if (.not. lpar(20)) then
        call parse(line,l,fnam1,lc)
        if(fnam1 .eq. ' ') fnam1 = jobnam
        call lower(fnam1)
        call upper(fnam1(1:1))
        do 2013 i = 1,4
           a(i) = 0.d0
           if(l .eq. 0) goto 2013
           call parse(line,l,code,lc)
           call getval(code(1:lc),' ',0,j,a(i))
2013    continue
        call hypxsc(fnam1,a)
*      else
*        write (6, 2012)
*2012    format(' Sorry, hyperfine cross sections not yet',
*     :         /,'  implemented for molecule-molecule collisions')
*      end if
      goto 1
* singlet-triplet collisional mixing - added by p. dagdigian
3000  call parse(line,l,fnam1,lc)
      if(fnam1 .eq. ' ') fnam1 = jobnam
      call lower(fnam1)
      call upper(fnam1(1:1))
* get iener for 1st smt file
      a(1) = 0.d0
      if(l .eq. 0) goto 3005
      call parse(line,l,code,lc)
      call getval(code(1:lc),' ',0,j,a(1))
3005  call parse(line,l,fnam2,lc)
      if(fnam2 .eq. ' ') fnam2 = jobnam
      call lower(fnam2)
      call upper(fnam2(1:1))
* get iener for 2nd smt file
      a(2) = 0.d0
      if(l .eq. 0) goto 3010
      call parse(line,l,code,lc)
      call getval(code(1:lc),' ',0,j,a(2))
* get dele, emax, istata, istatx, hso
3010  do 3020 i = 3, 7
        a(i) = 0.d0
        if(l .eq. 0) goto 3020
        call parse(line,l,code,lc)
        call getval(code(1:lc),' ',0,j,a(i))
3020  continue
      call stmix(fnam1,fnam2,a)
      goto 1
* transport cross sections - added by p. dagdigian
3100  call parse(line,l,fnam1,lc)
      if(fnam1 .eq. ' ') fnam1 = jobnam
      call lower(fnam1)
      call upper(fnam1(1:1))
* get iener for 1st smt file
      a(1) = 0.d0
      if(l .eq. 0) goto 3105
      call parse(line,l,code,lc)
      call getval(code(1:lc),' ',0,j,a(1))
* get in1, in2, jtotmx, join, jmax
3105  continue
      call trnprt(fnam1,a)
      goto 1
* pressure broadening cross sections - added by p. dagdigian
3200  call parse(line,l,fnam1,lc)
      if(fnam1 .eq. ' ') fnam1 = jobnam
      call lower(fnam1)
      call upper(fnam1(1:1))
* get iener for 1st smt file
      a(1) = 0.d0
      if(l .eq. 0) goto 3205
      call parse(line,l,code,lc)
      call getval(code(1:lc),' ',0,j,a(1))
3205  call parse(line,l,fnam2,lc)
      if(fnam2 .eq. ' ') fnam2 = jobnam
      call lower(fnam2)
      call upper(fnam2(1:1))
* get iener for 2nd smt file
      a(2) = 0.d0
      if(l .eq. 0) goto 3210
      call parse(line,l,code,lc)
      call getval(code(1:lc),' ',0,j,a(2))
* get k, j1, in1, j2, in2, diag, j1p, in1p, j2p, in2p
3210  do 3220 i = 3, 12
        a(i) = 0.d0
        if(l .eq. 0) goto 3220
        call parse(line,l,code,lc)
        call getval(code(1:lc),' ',0,j,a(i))
3220  continue
      call prsbr(fnam1,fnam2,a)
      goto 1
      end
      subroutine pcoder(boundc,pcod,icode)
*  subroutine to change pcod's for bound state or scattering
      logical boundc
      character*8 pcod(icode)
      if (boundc) then
        pcod(11)='R1'
        pcod(12)='R2'
        pcod(13)='C'
        pcod(14)='SPAC'
        pcod(15)='DELR'
        pcod(16)='HSIMP'
        pcod(17)='EIGMIN'
        pcod(18)='TOLAI'
      else
        pcod(11)='FSTFAC'
        pcod(12)='RINCR'
        pcod(13)='RCUT'
        pcod(14)='RENDAI'
        pcod(15)='RENDLD'
        pcod(16)='RSTART'
        pcod(17)='SPAC'
        pcod(18)='TOLAI'
      endif
      return
      end
