* -------------------------------------------------------------------------
      subroutine prsbr(flnam1, flnam2, a)
*
*  subroutine to compute pressure broadening cross sections
*
*  for the calculation of the cross section for an isolated line,
*  only one smt file is required.  to compute the cross section
*  coupling two spectral lines, two smt files are needed, one for
*  the initial and final levels of the spectral line
*
*  author:  p. dagdigian
*  date:  18-aug-2015
*
*  revision:  31-may-3013 by q. ma (dynamically allocate memory)
*  revision:   5-jun-2013 by q. ma (no longer use codim to determine sizes)
*  revision:  27-sep-2013 by p. dagdigian (relax requirement that jfinl1 = jfinl2)
*  revision:  18-aug-2015 by p. dagdigian (skip in sums if xf6j's equal zero)
* -------------------------------------------------------------------------
      use constants
      implicit double precision (a-h, o-z)
      character*(*) flnam1, flnam2
      character*20  cdate1, cdate2
      character*40  smtfil1, smtfil2, prtfil
      character*10  elaps, cpu
      character*1 slab
      complex*8 sa, sb, term
      logical csflg1, flghf1, flgsu1, twmol1, nucrs1,
     :        csflg2, flghf2, flgsu2, twmol2, nucrs2,
     :        batch, fast, lpar2, lpar, exstfl, diagst,
     :        diagj, diagin,
     :        diagjp, daginp, diag, diagp
      include "common/parpot"
      common /coisc1/ jlev1(1)
      common /coisc2/ jlev2(1)
      common /coisc3/ inlev1(1)
      common /coisc4/ inlev2(1)
      common /coisc5/ jout1(1)
      common /coisc6/ jout2(1)
      common /cosc1/ elev1(1)
      common /cosc2/ elev2(1)
      common /coisc7/ nlevt(1)
      common /coisc8/ jlevt(1)
      common /coisc9/ inlevt(1)
      common /cosc3/ elevt(1)
      common /coisc10/ ipack(1)
      common /coisc11/ jpack(1)
      common /coisc12/ lpack(1)
      common /codim/ nairy, mmax
      dimension a(12)
*
      integer, dimension(:), allocatable :: jq, lq, inq
      double precision, dimension(:), allocatable :: sreal, simag
* storage for s-matrix elements:
*   second letter is real (r), imaginary (i)
*   third letter for 1st (a) or 2nd (b) array
*   subscripts:  jtot, jlp (= 1/2 for jlpar = +1/-1), size of channel basis array
      double precision, dimension(:, :, :), allocatable ::
     $     sra, sia, srb, sib
* arrays with values of j, in, and l:
*   subscripts:  jtot, jlp (= 1/2 for jlpar = +1/-1), length of channel basis
      integer, dimension(:, :, :), allocatable ::
     $     ja, ina, la, jb, inb, lb
* length of arrays
*   last letter in name:  1st (a) or 2nd (b) array
*   subscripts:  jtot, jlp (= 1/2 for jlpar = +1/-1)
      integer, dimension(:, :), allocatable :: lngtha, lngthb
*
* storage for partial cross sections
      double precision, dimension(:), allocatable:: partxr, partxi
*
* to compute partial cross sections, set ipartx to 1,
* than recompile and relink
      ipartx = 1
*
* initialize timer and arrays
      call mtime(cpu0,ela0)
*
* generate filename of 1st smt file and check if it is present
*
      iener1 = a(1)
      call gennam(smtfil1,flnam1,iener1,'smt',lenfs)
      inquire(file = smtfil1, exist = exstfl)
      if (.not. exstfl) then
        write(6,10) smtfil1(1:lenfs)
10      format(/' ** FILE ',(a),' NOT FOUND **'/)
        return
      end if
*
* open 1st s-matrix file
*
      call openf(1,smtfil1,'tu',0)
*
* initialize timer and arrays
      call mtime(cpu0,ela0)
*
* generate filename of 2nd smt file and check if it is present
*
      iener2 = a(2)
      call gennam(smtfil2,flnam2,iener2,'smt',lenfs)
      inquire(file = smtfil2, exist = exstfl)
      if (.not. exstfl) then
        write(6,10) smtfil2(1:lenfs)
        close(1)
        return
      end if
*
* open 2nd s-matrix file
*
      call openf(11,smtfil2,'tu',0)
*
* read header of 1st s-matrix file
*
      call sinqr(1, m1jtot, m1chmx)
      call rdhead(1,cdate1,ered1,rmu1,csflg1,flghf1,flgsu1,
     :  twmol1,nucrs1,jfrst1,jfinl1,jtotd1,numin1,numax1,nud1,
     :  nlevel1,nlvop1,nnout1,jlev1,inlev1,elev1,jout1)
*
* we need the s-matrices as lower triangles, so nnout1  m u s t  be > 0
*
      if (nnout1.lt.0) then
        write(6,11)
11      format(/' ** NNOUT < 0, ABORT **'/)
        close(1)
        return
      end if
      nout = nnout
*
* read header of 2nd s-matrix file
*
      call sinqr(11, m2jtot, m2chmx)
      call rdhead(11,cdate2,ered2,rmu2,csflg2,flghf2,flgsu2,
     :   twmol2,nucrs2,jfrst2,jfinl2,jtotd2,numin2,numax2,nud2,
     :   nlevel2,nlvop2,nnout2,jlev2,inlev2,elev2,jout2)
*
* we need the s-matrices as lower triangles, so nnout2  m u s t  be > 0
*
      if (nnout2.lt.0) then
         write(6,11)
         close(1)
         close(11)
         return
      end if
      nout = nnout
*
*  molecule-molecule cross sections not implemented
*
      if (twmol1. or. twmol2) then
         write(6,12)
12       format(/' *** PRESSURE BROADENING CROSS SECTIONS FOR',
     :     'MOLECULE - MOLECULE COLLISIONS NOT IMPLEMENTED ***'/)
         close(1)
         close(11)
         return
      end if
*
*  molecule-surface collisions not implemented
*
      if (flgsu1 .or. flgsu2) then
         write(6,14)
14       format(/' *** PRESSURE BROADENING CROSS SECTIONS',
     :     ' FOR SURFACE COLLISIONS NOT IMPLEMENTED ***'/)
         close(1)
         close(11)
         return
      end if
*
*  CS cross sections not implemented
      if (csflg1 .or. csflg2) then
         write(6,16)
16       format(/' *** CS PRESSURE BROADENING CROSS',
     :     ' SECTIONS NOT IMPLEMENTED ***'/)
         close(1)
         close(11)
         return
      end if
*
*  delta-jtot should be equal to one for both s-matrix files
      if (iabs(jtotd1).ne.1 .or. iabs(jtotd2).ne.1)
     :    then
        write(6,116)
116     format(/' *** DELTA-JTOT MUST BE EQUAL TO ONE',
     :    ' *** '/)
         close(1)
         close(11)
        return
      end if
*
*  reduced mass in the two smt files must be the same
*
      if (rmu1 .ne. rmu2) then
         write(6,16)
19       format(/' *** REDUCED MASS IN BOTH SMT FILES',
     :          ' MUST BE EQUAL ***'/)
         close(1)
         close(11)
         return
      end if
*
*  minimum jtot should be the same for both smt files
*
      if (jfrst1.ne.jfrst2) then
        write (6,117)
117     format(/' *** JTOT1 SHOULD BE THE SAME FOR BOTH ',
     :    'SMT FILES *** '/)
        close(1)
        close(11)
        return
      end if
*
*  check maximum jtot for both smt files
*
      if (jfinl1.ne.jfinl2) then
        maxjtot = min(jfinl1,jfinl2)
        write (6,17) maxjtot
17      format(/' ++ JTOT2 NOT THE SAME FOR BOTH SMT FILES.',
     :    '  SET MAXJTOT =',i4)
      else
        maxjtot = jfinl1
      end if
*
*  flaghf should be the same for both smt files
*
      if (flghf1.ne.flghf2) then
        write (6,18)
18      format(/' *** FLAGHF SHOULD BE THE SAME FOR BOTH ',
     :    'SMT FILES ***'/)
        close(1)
        close(11)
        return
      end if
      spin = 0.d0
*
*  get tensor order
      if (flghf1) spin = 0.5d0
      xk = a(3)
      if (xk.ne.0.d0 .and. xk.ne.1.d0 .and.
     :    xk.ne.2.d0) then
        write(6,1010) xk
1010    format(/' *** INVALID VALUE OF K (=',f4.1,
     :      ') ***'/)
        close(1)
        close(11)
        return
      end if
      k = xk
*
*  maximum number of channels
      mchmx = max(m1chmx, m2chmx)
*  maximum value of jtot
      jtotmx = max(jfinl1,jfinl2)
*  maximum number of the elements of (the lower triangle of) a single s-matrix
      mmax2 = mchmx * (mchmx + 1) / 2
*  allocate memory
      allocate(jq(mchmx), stat=ialloc)
      if (ialloc .ne. 0) goto 4000
      allocate(lq(mchmx), stat=ialloc)
      if (ialloc .ne. 0) goto 4001
      allocate(inq(mchmx), stat=ialloc)
      if (ialloc .ne. 0) goto 4002
      allocate(sreal(mmax2), stat=ialloc)
      if (ialloc .ne. 0) goto 4003
      allocate(simag(mmax2), stat=ialloc)
      if (ialloc .ne. 0) goto 4004
      allocate(sra(0:jtotmx, 2, mmax2), stat=ialloc)
      if (ialloc .ne. 0) goto 4005
      allocate(sia(0:jtotmx, 2, mmax2), stat=ialloc)
      if (ialloc .ne. 0) goto 4006
      allocate(srb(0:jtotmx, 2, mmax2), stat=ialloc)
      if (ialloc .ne. 0) goto 4007
      allocate(sib(0:jtotmx, 2, mmax2), stat=ialloc)
      if (ialloc .ne. 0) goto 4008
      allocate(ja(0:jtotmx, 2, mchmx), stat=ialloc)
      if (ialloc .ne. 0) goto 4009
      allocate(ina(0:jtotmx, 2, mchmx), stat=ialloc)
      if (ialloc .ne. 0) goto 4010
      allocate(la(0:jtotmx, 2, mchmx), stat=ialloc)
      if (ialloc .ne. 0) goto 4011
      allocate(jb(0:jtotmx, 2, mchmx), stat=ialloc)
      if (ialloc .ne. 0) goto 4012
      allocate(inb(0:jtotmx, 2, mchmx), stat=ialloc)
      if (ialloc .ne. 0) goto 4013
      allocate(lb(0:jtotmx, 2, mchmx), stat=ialloc)
      if (ialloc .ne. 0) goto 4014
      allocate(lngtha(0:jtotmx, 2), stat=ialloc)
      if (ialloc .ne. 0) goto 4015
      allocate(lngthb(0:jtotmx, 2), stat=ialloc)
      if (ialloc .ne. 0) goto 4016
      allocate(partxr(0:jtotmx), stat=ialloc)
      if (ialloc .ne. 0) goto 4017
      allocate(partxi(0:jtotmx), stat=ialloc)
      if (ialloc .ne. 0) goto 4018
*
*  set up initial and final state labels
      diagst = .true.
      if (a(8) .ne. 1) diagst = .false.
      j1 = a(4)
      in1 = a(5)
      j2 = a(6)
      in2 = a(7)
      if (diagst) then
        j1p = j1
        in1p = in1
        j2p = j2
        in2p = in2
      else
        j1p = a(9)
        in1p = a(10)
        j2p = a(11)
        in2p = a(12)
      end if
*
*  check that initial level(s) are present
*  in the level lists
      flg = 0
      flgp = 0
      do i=1,nlevel1
        if(j1.eq.jlev1(i) .and. in1.eq.inlev1(i)) then
          flg = 1
*  prefactor for cross section
          etrans = ered1 - elev1(i)
          prefac = ang2c * 3.1415926535897932d0 /
     :        (2.d0 * rmu1 * etrans)
        end if
        if(j1p.eq.jlev1(i) .and. in1p.eq.inlev1(i))
     :    flgp = 1
      end do
      if (flg.ne.1 .or. flgp.ne.1) then
        write(6,1001) j1,in1,j1p,in1p
1001    format(/' *** LEVELS (j=',i2,', ',i4,') AND/OR (j=',
     :    i2,', ',i4,')'/'     NOT FOUND IN LEVEL LIST IN',
     :    ' 1st SMT FILE *** '/)
        goto 1000
      end if
*  check that final level(s) are present
*  in the level lists
      flg = 0
      flgp = 0
      do i=1,nlevel2
        if (j2.ne.jlev2(i) .and. in2.ne.inlev2(i))
     :    flg = 1
        if (j2p.ne.jlev2(i) .and. in2p.ne.inlev2(i))
     :    flgp = 1
      end do
      if (flg.ne.1 .or. flgp.ne.1) then
        write(6,1002) j2,in2,j2p,in2p
1002    format('/ *** LEVELS (j=',i2,', ',i4,') AND/OR (j=',
     :    i2,', ',i4,')'/'     NOT FOUND IN LEVEL LIST IN',
     :    ' 2nd SMT FILE *** '/)
        goto 1000
      end if
*
*  open file partial cross sections, if needed
      if (ipartx.eq.1) then
        call gennam(prtfil,flnam1,iener1,'ppb',lenpb)
        call openf(2,prtfil,'sf',0)
*  initialize partial cross section arrays
        do ii=0,jtotmx
          partxr(ii) = 0.d0
          partxi(ii) = 0.d0
        end do
      end if
*
*  set values for jlpar in s-matrices read
*  this works if j=0 level has jlpar=+1
*  needs fix for bach2x (not implemented here)
      jlpmx1 = 2
      jlpmx2 = 2
*
      write (6, 24) smtfil1, cdate1, jfinl1,
     :  smtfil2, cdate2, jfinl2
      if (ipartx.eq.1) then
        write (2, 24) smtfil1, cdate1, jfinl1,
     :    smtfil2, cdate2, jfinl2
      end if
24    format(/' PRESSURE BROADENING CROSS SECTION',
     :  /,' S-MATRICES FOR INITIAL LEVEL READ FROM',
     :  ' FILE ',(a)/'   WRITTEN:   ',(a),'  JTOT2=',i4/
     :  ' S-MATRICES FOR FINAL LEVEL READ FROM FILE ',
     :  (a)/'   WRITTEN:   ',(a),'  JTOT2=',i4)
*
*  compute collision energy for initial levels
      do 50 i = 1, nlevel1
        if (j1.eq.jlev1(i) .and. in1.eq.inlev1(i)) then
          ecoll1 = ered1 - elev1(i)
        end if
50    continue
      do 60 i = 1,nlevel2
        if (j2.eq.jlev2(i) .and. in2.eq.inlev2(i)) then
          ecoll2 = ered2 - elev2(i)
        end if
60    continue
*  check that collision energies agree to within 0.1 cm^-1
      if (abs(ecoll1 - ecoll2)*econv .gt. 0.1d0) then
        write (6,21) ecoll1 * econv, ecoll2 * econv
21      format (/' *** DIFFERENCE OF COLLISION ENERGIES',
     :    ' .GT. 0.1 cm^-1 ***'/
     :    5x,'ECOLL1 =',f11.3,	4x,'ECOLL2 =',f11.3,' CM(-1)'/)
        if (ipartx.eq.1) close(2)
        goto 1000
      end if
*
* now compute the pressure broadening cross section
      sigmar = 0.d0
      sigmai = 0.d0
      write(6,1665) ecoll1 * econv
      write(2,1665) ecoll1 * econv
1665  format(' COLLISION ENERGY = ',f9.3,' CM^(-1)')
      if (maxjtot.gt.jtotmx) then
        write(6,1093) jtotmx
1093    format(/' ** MAXJTOT .GT. JTOTMAX=',i3',. ABORT **'/)
        if (ipartx.eq.1) close(2)
        goto 1000
      end if
*
* read s-matrix elements for initial level
* this assumes that jlpar=1 is stored first
* scattering calculation should be carried out with jlpar = 0,
* except for bach2x with ortho levels and only j = 0 open
*
* parameter to read lower triangle of open channel(s)
* read s-matrix for present jtot, jlpar
      iaddr = 0
20    nopen = -1
      call sread (iaddr, sreal, simag, jtot1, jlpar1,
     :   nu1, jq, lq, inq, ipack, jpack, lpack,
     :   1, mmax, nopen, lngth1, ierr)
      if (ierr .lt. -1) then
        write(6,105)
105     format(/' ** READ ERROR IN PRSBR - INITIAL ',
     :      'SMT FILE. ABORT **'/)
        if (ipartx.eq.1) close(2)
        goto 1000
        return
      end if
      jlp = 1 - (jlpar1 - 1)/2
*  copy s-matrix for this jtot1/jlpar1
      lngtha(jtot1,jlp) = lngth1
      len2 = lngth1*(lngth1 + 1)/2
      do i = 1, lngth1
        ja(jtot1,jlp,i) = jpack(i)
        ina(jtot1,jlp,i) = ipack(i)
        la(jtot1,jlp,i) = lpack(i)
      end do
      do ii = 1, len2
        sra(jtot1,jlp,ii) = sreal(ii)
        sia(jtot1,jlp,ii) = simag(ii)
      end do
* check for end of s-matrix file
      if (jtot1.eq.jfinl1) then
        if (jlpar1.eq.1 .and. jlpmx1.eq.2) goto 20
        if (jlpar1.eq.1 .and. jlpmx1.eq.1) goto 22
        if (jlpar1.eq.-1) goto 22
      endif
* reset iaddr to 0 to insure sequential read
      iaddr = 0
* loop back for next partial wave
      goto 20
*
* read s-matrix elements for final level
* this assumes that jlpar=1 is stored first
* scattering calculation should be carried out with jlpar = 0,
* except for bach2x with ortho levels and only j = 0 open
*
* parameter to read lower triangle of open channel(s)
* read s-matrix for present jtot, jlpar
22    iaddr = 0
1020  nopen = -1
      call sread (iaddr, sreal, simag, jtot2, jlpar2,
     :   nu2, jq, lq, inq, ipack, jpack, lpack,
     :   11, mmax, nopen, lngth2, ierr)
      if (ierr .lt. -1) then
        write(6,1105)
1105    format(/' ** READ ERROR IN PRSBR - FINAL ',
     :      'SMT FILE. ABORT **'/)
        goto 1000
      end if
      jlp = 1 - (jlpar2 - 1)/2
*  copy s-matrix for this jtot2/jlpar2
      lngthb(jtot2,jlp) = lngth2
      len2 = lngth2*(lngth2 + 1)/2
      do i = 1, lngth2
        jb(jtot2,jlp,i) = jpack(i)
        inb(jtot2,jlp,i) = ipack(i)
        lb(jtot2,jlp,i) = lpack(i)
      end do
      do ii = 1, len2
        srb(jtot2,jlp,ii) = sreal(ii)
        sib(jtot2,jlp,ii) = simag(ii)
      end do
* check for end of s-matrix file
      if (jtot2.eq.jfinl2) then
        if (jlpar2.eq.1 .and. jlpmx2.eq.2) goto 1020
        if (jlpar2.eq.1 .and. jlpmx2.eq.1) goto 1022
        if (jlpar2.eq.-1) goto 1022
      endif
* reset iaddr to 0 to insure sequential read
      iaddr = 0
* loop back for next partial wave
      goto 1020
1022  continue
*
*  now accummulate terms in sum for cross section
*  sum over jtot/jlpar and jtotp/jlparp
      do jtot=jfrst1,maxjtot
        xjtot = jtot + spin
        facj = 2.d0 * xjtot + 1.d0
        do jlp=1,2
          jlpar = 1 - (jlp - 1) * 2
*  boundaries of sum over jtotp
          jtpmin = max((jtot - k),0)
          if (flghf1) then
            xjtpm = max((xjtot - xk),0.5d0)
            jtpmin = xjtpm - spin
          end if
          jtpmax = min((jtot + k),maxjtot)
          do jtotp = jtpmin,jtpmax
            xjtotp = jtotp + spin
            facjpp = (2.d0 * xjtotp + 1.d0) * facj
            do jlpp=1,2
              jlparp = 1 - (jlpp - 1) * 2
*  sum over row index for jtot
              do 400 irow=1,lngtha(jtot,jlp)
                jj1 = ja(jtot,jlp,irow)
*  channel basis function not for initial level
                if (jj1.ne.j1) goto 400
                if (ina(jtot,jlp,irow).ne.in1) goto 400
*
                xj1 = jj1 + spin
                l1 = la(jtot,jlp,irow)
                xl1 = l1
*  sum over column index for jtot
                do 200 icol=1,lngtha(jtot,jlp)
                  jj1p = ja(jtot,jlp,icol)
*  channel basis function not for (other) initial level
                  if (jj1p.ne.j1p) goto 200
                  if (ina(jtot,jlp,icol).ne.in1p) goto 200
*
                  xj1p = jj1p + spin
                  sqj1j1p = sqrt((2.d0 * xj1p + 1.d0)
     :                /(2.d0 * xj1 + 1.d0))
                  l1p = la(jtot,jlp,icol)
                  xl1p = l1p
                  diagj = jj1.eq.jj1p
                  diagin = ina(jtot,jlp,irow) .eq.
     :              ina(jtot,jlp,icol)
*  sum over row index for jtotp
                  do 401 irowp=1,lngthb(jtotp,jlpp)
                    jj2 = jb(jtotp,jlpp,irowp)
*  channel basis function not for final level
                    if (jj2.ne.j2) goto 401
                    if (inb(jtotp,jlpp,irowp).ne.in2) goto 401
*
                    xj2 = jj2 + spin
                    l2 = lb(jtotp,jlpp,irowp)
                    if (l2 .ne. l1) goto 401
                    xl2 = l2
*  sum over column indx for jtotp
                    do 201 icolp=1,lngthb(jtotp,jlpp)
                      jj2p = jb(jtotp,jlpp,icolp)
*  channel basis function not for (other) final level
                      if (jj2p.ne.j2p) goto 201
                      if (inb(jtotp,jlpp,icolp).ne.in2p)
     :                    goto 201
*
                      xj2p = jj2p + spin
                      l2p = lb(jtotp,jlpp,icolp)
                      if (l2p .ne. l1p) goto 201
                      xl2p = l2p
                      diagjp = jj2.eq.jj2p
                      daginp = inb(jtotp,jlpp,irowp) .eq.
     :                    inb(jtotp,jlpp,icolp)
*  get s-matrix elements
                      if (irow.gt.icol) then
                        ii = (irow*(irow - 1))/2 + icol
                      else
                        ii = (icol*(icol - 1))/2 + irow
                      end if
                      if (irowp.gt.icolp) then
                        iip = (irowp*(irowp - 1))/2 + icolp
                      else
                        iip = (icolp*(icolp - 1))/2 + irowp
                      end if
*  phase factor
                      phase = 1.d0
                      ipower = nint(xj1 - xj1p + xl1 - xl1p)
                      if ((ipower/2)*2 .ne. ipower) phase = -1.d0
                      sa = cmplx(sra(jtot,jlp,ii),
     :                  sia(jtot,jlp,ii))
*  complex conjugate
                      sb = cmplx(srb(jtotp,jlpp,iip),
     :                  -sib(jtotp,jlpp,iip))
                      term = - sa * sb
                      diag = diagj .and. diagin
                      diagp = diagjp .and. daginp
                      if (diag .and. diagp .and. (l1.eq.l1p))
     :                    term = cmplx(1.d0, 0.d0) + term
                      factor = phase * facjpp * sqj1j1p
                      x61 = xf6j(xj1p,xk,xj2p,xjtotp,xl1p,xjtot)
                      if (x61 .eq. 0.d0) goto 201
                      x62 = xf6j(xj1,xk,xj2,xjtotp,xl1,xjtot)
                      if (x62 .eq. 0.d0) goto 201
                      factor = factor * x61 * x62
                      sigmar = sigmar + real(factor*term)
                      sigmai = sigmai + aimag(factor*term)
                      if (ipartx.eq.1) then
                        partxr(jtot) = partxr(jtot)
     :                    + real(factor*term)
                        partxi(jtot) = partxi(jtot)
     :                    + aimag(factor*term)
                      end if
201                 continue
401               continue
200             continue
400           continue
            end do
          end do
        end do
      end do
      sigmar = sigmar * prefac
      sigmai = sigmai * prefac
*
*  print real and imag part of cross section
      if (.not. flghf1) then
        if (diagst) then
          write(6,2020) k,j1,in1,j2,in2,sigmar,sigmai
          if (ipartx.eq.1)
     :      write(2,2020) k,j1,in1,j2,in2,sigmar,sigmai
2020      format(/' DIAGONAL PRESSURE BROADENING ',
     :      'CROSS SECTION (ANG^2) - K =',i2/4x,'J1=',i3,
     :      ' IN1=',i3,' -> J2=',i3,' IN2=',i3/4x,
     :      'REAL PART =',1pe12.4,'  IMAG PART =',1pe12.4)
        else
          write(6,1030) k,j1,in1,j2,in2,j1p,in1p,j2p,in2p,
     :       sigmar,sigmai
          if (ipartx.eq.1)
     :      write(2,1030) k,j1,in1,j2,in2,j1p,in1p,j2p,in2p,
     :         sigmar,sigmai
1030      format(/' PRESSURE BROADENING ',
     :      'CROSS SECTION (ANG^2) - K =',i2/6x,'J1= ',i3,
     :      ' IN1= ',i3,' -> J2= ',i3,' IN2= ',i3/4x,'J1P=',i3,
     :      ' IN1P=',i3,' -> J2P=',i3,' IN2P=',i3/4x,
     :      'REAL PART =,'1pe12.4,'  IMAG PART =',1pe12.4)
        end if
      else
        xj1 = j1 + spin
        xj2 = j2 + spin
        if (diagst) then
          write(6,2120) k,xj1,in1,xj2,in2, sigmar,sigmai
          if (ipartx.eq.1)
     :        write(2,2120) k,xj1,in1,xj2,in2, sigmar,sigmai
2120      format(/' DIAGONAL PRESSURE BROADENING ',
     :      'CROSS SECTION (ANG^2) - K =',i2/6x,'J1= ',f4.1,
     :      ' IN1= ',i3,' -> J2=',f4.1,' IN2=',i3/
     :      ' REAL PART =',1pe12.4,'  IMAG PART =',1pe12.4)
        else
        xj1p = j1p + spin
        xj2p = j2p + spin
          write(6,1130) k,xj1,in1,xj2,in2,xj1p,in1p,xj2p,in2p,
     :       sigmar, sigmai
          if (ipartx.eq.1)
     :       write(2,1130) k,xj1,in1,xj2,in2,xj1p,in1p,
     :       xj2p,in2p,sigmar, sigmai
1130      format(/' PRESSURE BROADENING ',
     :      'CROSS SECTION (ANG^2) - K =',i2/4x,'J1= ',f4.1,
     :      ' IN1= ',i3,' -> J2= ',f4.1,' IN2= ',i3/
     :      6x,'J1P=',f4.1,' IN1P=',i3,' -> J2P=',f4.1,
     :      ' IN2P=',i3/4x,'REAL PART =,'1pe12.4,
     :      '  IMAG PART =',1pe12.4)
        end if
      end if
*
      if (ipartx.eq.1) then
        write (2,1140)
1140    format(/' PARTIAL CROSS SECTIONS (ANG^2)'/
     :    6x,'JTOT',8x,'REAL PART',10x,'IMAG PART')
        do jtot=jfrst1,maxjtot
          realp = partxr(jtot) * prefac
          aimagp = partxi(jtot) * prefac
          if (.not. flghf1) then
            write(2,1150) jtot, realp, aimagp
1150        format(i10,6x,1pe12.4,7x,1pe12.4)
          else
            xjtot = jtot + spin
            write(2,1152) xjtot, realp, aimagp
1152        format(f10.1,6x,1pe12.4,7x,1pe12.4)
          end if
        end do
      end if
*
      call mtime(cpu1,ela1)
      ela1 = ela1 - ela0
      cpu1 = cpu1 - cpu0
      call gettim(ela1,elaps)
      call gettim(cpu1,cpu)
      write(6,900) elaps, cpu
900   format(/,' ** PRSBR FINAL TIMING, ELAPSED: ',a,'  CPU: ',a,' **'/)
      ialloc = 0
 1000 continue
 4019 deallocate(partxi)
 4018 deallocate(partxr)
 4017 deallocate(lngthb)
 4016 deallocate(lngtha)
 4015 deallocate(lb)
 4014 deallocate(inb)
 4013 deallocate(jb)
 4012 deallocate(la)
 4011 deallocate(ina)
 4010 deallocate(ja)
 4009 deallocate(sib)
 4008 deallocate(srb)
 4007 deallocate(sia)
 4006 deallocate(sra)
 4005 deallocate(simag)
 4004 deallocate(sreal)
 4003 deallocate(inq)
 4002 deallocate(lq)
 4001 deallocate(jq)
 4000 close(1)
      close(11)
      if (ialloc .ne. 0) write (6, 4100)
 4100 format (' *** INSUFFICIENT MEMORY OR SMT FILE CORRUPTED. ***')
      return
      end
* -------------------------------eof---------------------------------------
