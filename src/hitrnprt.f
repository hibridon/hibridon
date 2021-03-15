* ------------------------------------------------------------------
      subroutine trnprt(filnam,a)
*
* subroutine to calculate effective cross sections for transport
* from s-matrix elements
*
* see G. C. Maitland, M. Mustafa, W. A. Wakeham, and F. R. McCourt,
* Mol. Phys. 61, 359 (1987); L. Monchick and S. Green, J. Chem. Phys.
* 63, 2000 (1975); G. A. Parker and R. T Pack, J. Chem. Phys. 68,
* 1585 (1978)
*
* this subroutine returns the effective cross sections for diffusion and
* viscosity, Q(n)(j1->j2)(E), where n=1 and 2.
*
* this subroutine requires CC s-matrix for both parities
*
* author: p.j. dagdigian
*
* changed open command to open smt files for streaming reads
*   revised:  9-jan-2012 by q. ma
* increase LENMX to 800 (9-sep-2012 by pjd)
* simplify input to TRNPRT.  calculate full table of cross sections
*   (similar to INTCRS command) (17-oct-2012 by pjd)
* fix for low-energy calculation where only j=0 is open (e.g. ch2a-he
*   for para levels and ch2x-he for ortho levels) (7-jan-2013)
* allocate memory dynamically; clean up unused variable (4-jun-2013
*   by q. ma)
*
* wrote second subroutine for computing transport cross sections for
* sysfems for which the j12 array is required (5-jun-2015 by pjd)
*
* current revision date:  24-jul-2015 by pjd
* ------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      character*(*) filnam
      character*40  trnfil, smtfil
      character*20  cdate
      character*10  elaps, cpu
      logical csflag, flaghf, flagsu, twomol, exstfl, lpar,
     :        batch, nucros,lpar2
*
      include "common/parpot"
      common /colpar/ lpar(3), batch, lpar2(23)
      common /cospbf/ lnbufs, lnbufl, nbuf, maxlsp, maxllb, ihibuf,
     :                igjtp
      common /cotrn/  spin, maxjt,
     :     nwaves, jfsts, jlparf, jlpars, njmax, jpmax
      common /cosout/ nnout, jout(25)
      common /coisc2/ inlev(1)
      common /coisc3/ jlev(1)
      common /coisc4/ jpack(1)
      common /coisc5/ lpack(1)
      common /coisc6/ ipack(1)
      common /coisc8/ jlist(1)
      common /coisc9/ jslist(1)
      common /coisc10/ inlist(1)
      common /cojq/   jq(1)
      common /colq/   lq(1)
      common /coinq/  inq(1)
      common /cosc1/  elev(1)
      common /cosc2/  prefac(1)
      common /cosc3/  etrans(1)
      common /cozmat/ jtotpa(1)
      common /coamat/ labadr(1)
      common /cotble/ npnt, jttble(1)
      common /coconv/ econv, xmconv, ang2c
      common /codim/  nairy, mmax
      common /coselb/ ibasty
      save nout
      dimension a(6)
      data  tol,   zero,   nstep
     :    /1.d-7,  0.d0,     2/
*
* initialize timer
      call mtime(cpu0,ela0)
* input
      iener = a(1)
      if (iener.le. 0) iener = 1
*
* generate filename and check if it is present
*
      call gennam(smtfil,filnam,iener,'smt',lenfs)
      inquire(file = smtfil, exist = exstfl)
      if (.not. exstfl) then
          write(6,10) smtfil(1:lenfs)
10        format(' ** FILE ',(a),' NOT FOUND **')
          return
      end if
*
* open smatrix-file
*
      call openf(1,smtfil,'tu',0)
*
* open file for transport cross sections
      call gennam(trnfil,filnam,iener,'trn',lenft)
      call openf(2,trnfil,'sf',0)
*
* read header of smatrix-file
*
      call sinqr(1, m1jtot, m1chmx)
      call rdhead(1,cdate,ered,rmu,csflag,flaghf,flagsu,twomol,
     :   nucros,jfirst,jfinal,jtotd,numin,numax,nud,nlevel,nlevop,nnout,
     :   jlev,inlev,elev,jout)
*
* we need the s-matrices as lower triangles, so nnout  m u s t  be > 0
*
      if (nnout.lt.0) then
         write(2,11)
         if(.not.batch) write(6,11)
11       format(' ** NNOUT < 0, ABORT **')
         return
      end if
      nout = nnout
*

*  molecule-molecule cross sections not yet implemented
*
*      if (twomol) then
*         write(2,12)
*         if (.not. batch) write(6,12)
*12       format(' *** TRANSPORT CROSS SECTIONS FOR MOLECULE -',
*     :          ' MOLECULE COLLISIONS NOT YET IMPLEMENTED ***')
*         goto 300
*      end if

      if (flagsu) then
         write(2,14)
         if(.not. batch) write(6,14)
14       format(' *** TRANSPORT CROSS SECTIONS FOR SURFACE -',
     :          ' COLLISIONS NOT IMPLEMENTED ***')
         goto 300
      end if
      if (csflag) then
         write(2,16)
         if(.not. batch) write(6,16)
16       format(' *** CS TRANSPORT CROSS SECTIONS',
     :          ' NOT IMPLEMENTED ***')
         goto 300
      end if
*
      spin = 0.d0
      if(flaghf) then
         spin = 0.5d0
      end if
*
      write (2, 20) smtfil, cdate, label, potnam
      if(.not. batch) write (6, 20) smtfil, cdate, label,potnam
20    format(/' CLOSE COUPLED TRANSPORT CROSS SECTIONS',/,
     :        ' S-MATRICES READ FROM FILE ',(a),/,
     : '      WRITTEN:   ',(a),/,
     : '      LABEL:     ',(a)/,
     : '      POT NAME:  ',(a))
*
* obtain new date
*
      call dater(cdate)
      write(2, 22) cdate
      if(.not. batch) write(6, 22) cdate
22    format(' DATE:    ',(a))
* set maxjt to jfinal
      maxjt = jfinal
*
* set up list of open levels and calculate prefactors
* save pointer in array jlist
* also save in values in array inlist (pjd)
      nj = 0
      jpmax = -1
      do 40 i=1, iabs(nout)
        jo = jout(i)
        do 30 j=1, nlevop
          j1 = jlev(j)
          in = inlev(j)
          if (j1.ne.jo) goto 30
          nj = nj + 1
          jlist(nj) = j
          jslist(nj) = j1
          inlist(nj) = in
          etrans(nj) = ered - elev(j)
          prefac(nj) = ang2c * 3.1415926535897932d0 /
     :        (2.d0* rmu * etrans(nj))
30      continue
40    continue
      njmax = nj
* set values for jlpar in s-matrix read
* this works if j=0 level has jlpar=+1
* needs fix for bach2x
      jlpmin = 1
      jlpmax = 2
*
      if (nlevop.eq.1 .and. jlev(1).eq.0) then
        jlpmax = 1
* special fix for bach2x since jlpar=-1 for j=0 here
        if (ibasty.eq.17) then
          jlpmin = 2
          jlpmax = 2
        end if
* special fix for bah3p/O(3PJ)
        if (ibasty.eq.13) jlpmax = 2
      end if
* write header
      write(2,60) ered*econv,rmu*xmconv,jfirst,maxjt
      if(.not.batch) write(6,60)
     :  ered*econv,rmu*xmconv,jfirst,maxjt
60    format(/,' ENERGY: ',f11.3,' cm(-1)    MASS: ',f11.3,/,
     :     ' SUMMING PARTIAL WAVES FROM JTOT=',i3,' TO JTOT=',i3)
*
* write level list
      write(2,64)
      if(.not. batch) write(6,64)
64    format(/,' COLUMNS ARE INITIAL STATES; ROWS ARE FINAL STATES')
      write(2,65)
      if(.not. batch) write(6,65)
65    format (/,' LEVEL LIST FOR TRANSPORT CROSS SECTIONS',
     :        ' (OPEN CHANNELS)',
     :        /,'   N     J   INDEX  EINT(cm-1)')
      do 68 i = 1, nj
         ii=jlist(i)
         if (.not. flaghf) then
           write(2,66) i,jlev(ii),inlev(ii),elev(ii)*econv
           if(.not.batch) write(6,66) i,jlev(ii),inlev(ii),
     :                                elev(ii)*econv
66         format (i4, 1x, i5, i6, f11.3)
         else
           write(2,67) i,jlev(ii)+spin,inlev(ii),elev(ii)*econv
           if(.not.batch) write(6,67) i,jlev(ii)+spin,inlev(ii),
     :                                elev(ii)*econv
67         format (i4, 1x, f5.1, i6, f11.3)
         endif

68    continue
*
* now compute transport cross sections
      if (is_twomol(ibasty) .or. is_j12(ibasty)) then
        call sigkp1(1,m1jtot,m1chmx,
     :   nnout,jfirst,jfinal,jtotd,nj,mmax,jpack,
     :   lpack,ipack,jttble,prefac,
     :   etrans,labadr,
     :   jtotpa,jlpmin,jlpmax,flaghf,ierr)
      else
        call sigkap(1,m1jtot,m1chmx,
     :   nnout,jfirst,jfinal,jtotd,nj,mmax,jpack,
     :   lpack,ipack,jttble,prefac,
     :   etrans,labadr,
     :   jtotpa,jlpmin,jlpmax,ierr)
      endif
      goto 300
300   close (1)
      close (2)
      close (3)
      close (4)
      return
999   write(2,1000)
      if(.not.batch) write(6,1000)
1000  format(' *** I/O ERROR IN TRNPRT.  ABORT')
      close (1)
      close (2)
      close (3)
      rewind (4)
      close (4)
      return
      end
* ------------------------------------------------------------------
      subroutine sigkap(iunit,mjtot,mchmx,
     :                nnout,jfirst,jfinal,jtotd,nj,mmax,jpack,
     :                lpack,ipack,jttble,prefac,
     :                etrans,labadr,
     :                jtotpa,jlpmin,jlpmax,ierr)
*
* subroutine to calculate Q(n) effective cross sections
*
* see G. C. Maitland, M. Mustafa, W. A. Wakeham, and F. R. McCourt,
* Mol. Phys. 61, 359 (1987); L. Monchick and S. Green, J. Chem. Phys.
* 63, 2000 (1975); G. A. Parker and R. T Pack, J. Chem. Phys. 68,
* 1585 (1978)
*
* subr first computes integral of product of differential cross section
* P_k(cos_theta).  the Q(n) cross sections are weighted sums of these
* integrals
*
* author: p. dagdigian
*
* revision: 7-jan-2013 by pj dagdigian
* allocate memory dynamically; clean up unused variable (4-jun-2013
*   by q. ma)
*
* this subroutine computes transport cross sections for systems in which
* the j12 array is not needed
*
* current revision date:  15-jun-2015 by pjd
*------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      complex*8 t, tp
      logical diag, diagj, diagin, diagp, diagjp, diagnp,
     :  lpar1, lpar2, batch, ipos
      character*10 elaps, cpu
      common /colpar/ lpar1(3), batch, lpar2(5), ipos, slpar3(17)
      common /cospbf/ lnbufs, lnbufl, nbuf, maxlsp, maxllb, ihibuf,
     :                igjtp
      common /cotrn/  spin, maxjt,
     :     nwaves, jfsts, jlparf, jlpars, njmax, jpmax
      common /coipar/ ipar(9),iprnt
* common blocks for levels for which xs's to be computed
      common /coisc9/ jslist(1)
      common /coisc10/ inlist(1)
      common /cojq/ jq(1)
      common /colq/ lq(1)
      common /coinq/ inq(1)
      common /coz/ sreal(1)
      common /cow/ simag(1)
      common /coselb/ ibasty
*
      dimension etrans(1)
      dimension jpack(1),ipack(1),lpack(1)
      dimension prefac(1), labadr(1), jtotpa(1)
      dimension jttble(1)
      dimension f36j(0:2)
*
* sr, si: s-matrix elements
*   second letter is real (r), imaginary (i)
*   subscripts:  jtot, jlp (=1/2 for jlpar = +1/-1), size of channel basis array
* j, in, l: arrays with values of j ,in, and l:
*   subscripts:  jtot, jlp (=1/2 for jlpar = +1/-1), length of channel basis
* sigma: array to hold transport cross sections
      double precision, dimension(:, :, :), allocatable ::
     $     sr, si, j, in, l, sigma
* length of arrays
*   subscripts:  jtot, jlp (=1/2 for jlpar = +1/-1)
      double precision, dimension(:, :), allocatable :: length
*
      onesix = 1.d0/6.d0
      twothr = 2.d0/3.d0
*
* initialize timer
      call mtime(cpu0,ela0)
      ierr = 0
* allocate and clear sigma array
      allocate(sigma(njmax, njmax, 2), stat=ialloc)
      if (ialloc .ne. 0) goto 4000
      sigma = 0d0
* allocate storage for s-matrices
      mmax2 = mchmx * (mchmx + 1) / 2
      allocate(sr(0:mjtot, 2, mmax2), stat=ialloc)
      if (ialloc .ne. 0) goto 4001
      allocate(si(0:mjtot, 2, mmax2), stat=ialloc)
      if (ialloc .ne. 0) goto 4002
      allocate(j(0:mjtot, 2, mchmx), stat=ialloc)
      if (ialloc .ne. 0) goto 4003
      allocate(in(0:mjtot, 2, mchmx), stat=ialloc)
      if (ialloc .ne. 0) goto 4004
      allocate(l(0:mjtot, 2, mchmx), stat=ialloc)
      if (ialloc .ne. 0) goto 4005
*
      allocate(length(0:mjtot, 2), stat=ialloc)
      if (ialloc .ne. 0) goto 4006
*
* read S-matrix elements
* this assumes that jlpar=1 is stored first
* scattering calculation should be carried out with jlpar = 0,
* except for bach2x with ortho levels and only j = 0 open
*
* parameter to read lower triangle of open channel(s)
* read s-matrix for present jtot, jlpar
      iaddr = 0
20    nopen = -1
      call sread (iaddr, sreal, simag, jtot, jlpar,
     :   nu, jq, lq, inq, ipack, jpack, lpack,
     :   iunit, mmax, nopen, lngth, ierr)
      if (ierr .lt. -1) then
        write(6,105)
105     format(/' ** READ ERROR IN TRNPRT. ABORT **'/)
        goto 1000
      end if
      jlp = 1 - (jlpar - 1)/2
*  copy s-matrix for this jtot1/jlpar1
      length(jtot,jlp) = lngth
      len2 = lngth*(lngth + 1)/2
      do i = 1, lngth
        j(jtot,jlp,i) = jpack(i)
        in(jtot,jlp,i) = ipack(i)
        l(jtot,jlp,i) = lpack(i)
      end do
      do ii = 1, len2
        sr(jtot,jlp,ii) = sreal(ii)
        si(jtot,jlp,ii) = simag(ii)
      end do
* check for end of s-matrix file
      if (jtot.eq.maxjt) then
          if (jlpar.eq.1 .and. jlpmax.eq.2) goto 20
          if (jlpar.eq.1 .and. jlpmax.eq.1) goto 22
          if (jlpar.eq.-1) goto 22
      endif
* reset iaddr to 0 to insure sequential read
      iaddr = 0
* loop back for next partial wave
      goto 20
*
* now sum over the jtot/jlpar pairs
22    do jtot = jfirst, maxjt
*
        if (jtot.eq.0) write(6,*) '  '
        if (jtot.eq.10*(jtot/10) .and. jtot.ne.0) write(6,210) jtot
210     format('++ STARTING JTOT =',i4)
*
        xjtot = jtot + spin
        facj = 2.d0 * xjtot + 1.d0
        do jlp = jlpmin, jlpmax
          jlpar = 1 - (jlp - 1)*2
*
* boundaries of sum over jtotp
          jtpmin = max((jtot - 2), jfirst)
          jtpmax = min((jtot + 2), maxjt)
          do jtotp = jtpmin, jtpmax
            xjtotp = jtotp + spin
            facjjp = (2.d0 * xjtotp + 1.d0) * facj
            do jlpp = jlpmin, jlpmax
              jlparp = 1 - (jlpp - 1)*2
*
* sum over row index for jtot
              do 400 irow = 1, length(jtot,jlp)
                j1 = j(jtot,jlp,irow)
                xj1 = j1 + spin
                l1 = l(jtot,jlp,irow)
                xl1 = l1
                sl1 = sqrt(2.d0 * xl1 + 1.d0)
                do 74 ij = 1, njmax
                  if (j1 .ne. jslist(ij)) goto 74
                  if (in(jtot,jlp,irow) .ne. inlist(ij)) goto 74
                  ir = ij
                  denrow = prefac(ij)
                  etrow = etrans(ij)
                  goto 75
74              continue
* level not in list
                goto 400
75              continue
*
* sum over column index for jtot
                do 200 icol = 1, length(jtot,jlp)
                  j2 = j(jtot,jlp,icol)
                  xj2 = j2 + spin
                  l2 = l(jtot,jlp,icol)
                  xl2 = l2
                  sl2 = sqrt(2.d0 * xl2 + 1.d0)
* replace check for diagonal s matrix element - this allows for checking
* matching j, l, and in
                  diag = irow .eq. icol
*                  diagj = j1 .eq. j2
*                  diagin = in(jtot,j1p,irow) .eq. in(jtot,j1p,icol)
*                  diag = diagj .and. diagin
                  do 84 ij = 1, njmax
                    if (j2 .ne. jslist(ij)) goto 84
                    if (in(jtot,jlp,icol) .ne. inlist(ij)) goto 84
                    ic = ij
                    etcol = etrans(ij)
                    goto 85
84                continue
* level not in list
                  goto 200
85                continue
*
* sum over row index for jtotp
                  do 401 irowp = 1, length(jtotp,jlpp)
                    j1p = j(jtotp,jlpp,irowp)
                    xj1p = j1p + spin
                    if (j1p .ne. j1) goto 401
                    if (in(jtotp,jlpp,irowp) .ne. in(jtot,jlp,irow))
     :                  goto 401
                    l1p = l(jtotp,jlpp,irowp)
                    xl1p = l1p
                    sl1p = sqrt(2.d0 * xl1p + 1.d0)
*
* sum over column index for jtotp
                    do 201 icolp = 1, length(jtotp,jlpp)
                      j2p = j(jtotp,jlpp,icolp)
                      xj2p = j2p + spin
                      if (j2p .ne. j2) goto 201
                      if (in(jtotp,jlpp,icolp) .ne. in(jtot,jlp,icol))
     :                    goto 201
                      l2p = l(jtotp,jlpp,icolp)
                      xl2p = l2p
                      sl2p = sqrt(2.d0 * xl2p + 1.d0)
* replace check for diagonal s matrix element - this allows for checking
* matching j, l, and in
                       diagp = irowp .eq.icolp
*                      diagjp = j1p .eq. j2p
*                      diagnp = in(jtotp,j1pp,irowp)
*     :                    .eq. in(jtotp,j1pp,icolp)
*                      diagp = diagjp .and. diagnp
* get s-matrix elements
                      if (irow .gt. icol) then
                        ii = (irow*(irow - 1))/2 + icol
                      else
                        ii = (icol*(icol - 1))/2 + irow
                      end if
                      if (irowp .gt. icolp) then
                        iip = (irowp*(irowp - 1))/2 + icolp
                      else
                        iip = (icolp*(icolp - 1))/2 + irowp
                      end if
* phase factor
                      phase = 1.d0
                      ipower = nint(xj1 - xj2
     :                    - 0.5d0*(xl1 - xl1p + xl2 - xl2p))
                      if ((ipower/2)*2 .ne. ipower) phase = -1.d0
* get T-matrix elements
                      t = - cmplx(sr(jtot,jlp,ii), si(jtot,jlp,ii))
*                      if (diag .and. (l1.eq.l2))
                      if (diag)
     :                    t = (1.d0, 0.d0) + t
                      tp=-cmplx(sr(jtotp,jlpp,iip),si(jtotp,jlpp,iip))
*                      if (diagp .and. (l1p.eq.l2p))
                      if (diagp)
     :                     tp = (1.d0, 0.d0) + tp
                      factor = phase * facjjp * sl1 * sl2 * sl1p
     :                    * sl2p * real(t*conjg(tp))/(2.d0*xj1 + 1.d0)
*
* loop over k
                      do 100 k = 0, 2
                        xk = k
                        f36j(k) = 0.d0
                        x1 = xf3j(xl1,xl1p,xk,0.d0,0.d0,0.d0)
                        if (x1 .eq. 0.d0) goto 100
                        x2 = xf3j(xl2,xl2p,xk,0.d0,0.d0,0.d0)
                        if (x2 .eq. 0.d0) goto 100
                        x3 = xf6j(xl1,xjtot,xj1,xjtotp,xl1p,xk)
                        if (x3 .eq. 0.d0) goto 100
                        x4 = xf6j(xl2,xjtot,xj2,xjtotp,xl2p,xk)
                        if (x4 .eq. 0.d0) goto 100
                        f36j(k) = factor * x1*x2*x3*x4*(-1.d0)**k
100                   continue
                      epove = etcol/etrow
*
* N = 1
                      fac1 = f36j(0) - sqrt(epove)*f36j(1)
                      sigma(ic,ir,1) = sigma(ic,ir,1) + fac1*denrow
*
* N = 2
                      fac2 = onesix*(5.d0 - epove**2)*f36j(0)
     :                    - twothr*epove*f36j(2)
                      sigma(ic,ir,2) = sigma(ic,ir,2) + fac2*denrow
201                 continue
401               continue
200             continue
400           continue
            end do
          end do
        end do
      end do
*
* here when calculation has been done
999   do 1100 k = 1, 2
        write(2,800) k
        if(.not.batch) write(6,800) k
800     format(/' N =',i3)
        call mxoutd (2, sigma(1,1,k), njmax, njmax, 0, ipos)
* use diag as scratch variable (ipos = .false. for screen output)
        diag = .false.
        if(.not. batch)
     :    call mxoutd (6, sigma(1,1,k), njmax, njmax, 0, diag)
1100  continue
*
      call mtime(cpu1,ela1)
      cpu1 = cpu1 - cpu0
      ela1 = ela1 - ela0
      call gettim(ela1,elaps)
      call gettim(cpu1,cpu)
      write(2,1200) elaps, cpu
      if(.not. batch) write(6,1200) elaps, cpu
1200  format(/' ** TIMING ELAPSED: ',a,
     :        ' CPU: ',a,/)
*
 1000 continue
 4007 deallocate(length)
 4006 deallocate(l)
 4005 deallocate(in)
 4004 deallocate(j)
 4003 deallocate(si)
 4002 deallocate(sr)
 4001 deallocate(sigma)
 4000 continue
      if (ialloc .ne. 0) write (6, 4100)
 4100 format (' *** INSUFFICIENT MEMORY. ***')
      return
      end
* ------------------------------------------------------------------
      subroutine sigkp1(iunit,mjtot,mchmx,
     :                nnout,jfirst,jfinal,jtotd,nj,mmax,jpack,
     :                lpack,ipack,jttble,prefac,
     :                etrans,labadr,
     :                jtotpa,jlpmin,jlpmax,flaghf,ierr)
*
* subroutine to calculate Q(n) effective cross sections
*
* see G. C. Maitland, M. Mustafa, W. A. Wakeham, and F. R. McCourt,
* Mol. Phys. 61, 359 (1987); L. Monchick and S. Green, J. Chem. Phys.
* 63, 2000 (1975); G. A. Parker and R. T Pack, J. Chem. Phys. 68,
* 1585 (1978)
*
* subr first computes integral of product of differential cross section
* P_k(cos_theta).  the Q(n) cross sections are weighted sums of these
* integrals
*
* author: p. dagdigian
*
* revision: 7-jan-2013 by pj dagdigian
* allocate memory dynamically; clean up unused variable (4-jun-2013
*   by q. ma)
*
* this version computes transport cross sections for systems in which
* the j12 array is needed, i.e. those for which is_twomol(ibasty) or
* is_j12(ibasty) is true.
*
* current revision date:  23-jun-2015 by pjd
*------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      complex*8 t, tp
      logical diag, diagj, diagin, diagp, diagjp, diagnp,
     :  lpar1, lpar2, batch, ipos, flaghf
      character*10 elaps, cpu
      common /colpar/ lpar1(3), batch, lpar2(5), ipos, slpar3(17)
      common /cospbf/ lnbufs, lnbufl, nbuf, maxlsp, maxllb, ihibuf,
     :                igjtp
      common /cotrn/  spin, maxjt,
     :     nwaves, jfsts, jlparf, jlpars, njmax, jpmax
      common /coipar/ ipar(9),iprnt
* common blocks for levels for which xs's to be computed
      common /coisc9/ jslist(1)
      common /coisc10/ inlist(1)
      common /cojq/ jq(1)
      common /colq/ lq(1)
      common /coinq/ inq(1)
      common /coj12/ j12(1)
      common /coz/ sreal(1)
      common /cow/ simag(1)
      common /coselb/ ibasty
*
      dimension etrans(1)
      dimension jpack(1),ipack(1),lpack(1)
      dimension prefac(1), labadr(1), jtotpa(1)
      dimension jttble(1)
      dimension f36j(0:2)
*
* sr, si: s-matrix elements
*   second letter is real (r), imaginary (i)
*   subscripts:  jtot, jlp (=1/2 for jlpar = +1/-1), size of channel basis array
* j, in, l, jj12: arrays with values of j ,in, l, and j12:
*   subscripts:  jtot, jlp (=1/2 for jlpar = +1/-1), length of channel basis
* sigma: array to hold transport cross sections
      double precision, dimension(:, :, :), allocatable ::
     $     sr, si, j, in, l, jj12, sigma
     :  , plam
*
* length of arrays
*   subscripts:  jtot, jlp (=1/2 for jlpar = +1/-1)
      double precision, dimension(:, :), allocatable :: length
*
      onesix = 1.d0/6.d0
      twothr = 2.d0/3.d0
*
      spnj2 = 0.d0
      spnj12 = 0.d0
      spntot = 0.d0
      if (flaghf) then
        spnj12 = 0.5d0
        spntot = 0.5d0
        if (ibasty.eq.12 .or. ibasty.eq.15) then
          spin = 0.d0
          spnj2 = 0.5d0
        endif
      endif
      if (ibasty.eq.23) then
        spnj2 = 0.5d0
        spnj12 = 0.5d0
        spntot = 0.5d0
      endif
*
* initialize timer
      call mtime(cpu0,ela0)
      ierr = 0
* allocate and clear sigma array
      allocate(sigma(njmax, njmax, 2), stat=ialloc)
      if (ialloc .ne. 0) goto 4000
      sigma = 0d0
*
*
      iplmt = 1
*
*
      allocate(plam(njmax, njmax, 3), stat=ialloc)
      if (ialloc .ne. 0) goto 4400
      plam = 0d0
*
*
* allocate storage for s-matrices
      mmax2 = mchmx * (mchmx + 1) / 2
      allocate(sr(0:mjtot, 2, mmax2), stat=ialloc)
      if (ialloc .ne. 0) goto 4001
      allocate(si(0:mjtot, 2, mmax2), stat=ialloc)
      if (ialloc .ne. 0) goto 4002
      allocate(j(0:mjtot, 2, mchmx), stat=ialloc)
      if (ialloc .ne. 0) goto 4003
      allocate(in(0:mjtot, 2, mchmx), stat=ialloc)
      if (ialloc .ne. 0) goto 4004
      allocate(l(0:mjtot, 2, mchmx), stat=ialloc)
      if (ialloc .ne. 0) goto 4005
      allocate(jj12(0:mjtot, 2, mchmx), stat=ialloc)
      if (ialloc .ne. 0) goto 4005
*
      allocate(length(0:mjtot, 2), stat=ialloc)
      if (ialloc .ne. 0) goto 4006
*
* read S-matrix elements
* this assumes that jlpar=1 is stored first
* scattering calculation should be carried out with jlpar = 0,
* except for bach2x with ortho levels and only j = 0 open
*
* parameter to read lower triangle of open channel(s)
* read s-matrix for present jtot, jlpar
      iaddr = 0
20    nopen = -1
      call sread (iaddr, sreal, simag, jtot, jlpar,
     :   nu, jq, lq, inq, ipack, jpack, lpack,
     :   iunit, mmax, nopen, lngth, ierr)
      if (ierr .lt. -1) then
        write(6,105)
105     format(/' ** READ ERROR IN TRNPRT. ABORT **'/)
        goto 1000
      end if
      jlp = 1 - (jlpar - 1)/2
*  copy s-matrix for this jtot1/jlpar1
      length(jtot,jlp) = lngth
      len2 = lngth*(lngth + 1)/2
      do i = 1, lngth
        j(jtot,jlp,i) = jpack(i)
        in(jtot,jlp,i) = ipack(i)
        l(jtot,jlp,i) = lpack(i)
        jj12(jtot,jlp,i) = j12(i)
      end do
      do ii = 1, len2
        sr(jtot,jlp,ii) = sreal(ii)
        si(jtot,jlp,ii) = simag(ii)
      end do
* check for end of s-matrix file
      if (jtot.eq.maxjt) then
          if (jlpar.eq.1 .and. jlpmax.eq.2) goto 20
          if (jlpar.eq.1 .and. jlpmax.eq.1) goto 22
          if (jlpar.eq.-1) goto 22
      endif
* reset iaddr to 0 to insure sequential read
      iaddr = 0
* loop back for next partial wave
      goto 20
*
* now sum over the jtot/jlpar pairs
22    do jtot = jfirst, maxjt
*
        if (jtot.eq.0) write(6,*) '  '
        if (jtot.eq.10*(jtot/10) .and. jtot.ne.0) write(6,210) jtot
210     format('++ STARTING JTOT =',i4)
*
        xjtot = jtot + spntot
        facj = 2.d0 * xjtot + 1.d0
        do jlp = jlpmin, jlpmax
          jlpar = 1 - (jlp - 1)*2
*
* boundaries of sum over jtotp
          jtpmin = max((jtot - 2), jfirst)
          jtpmax = min((jtot + 2), maxjt)
          do jtotp = jtpmin, jtpmax
            xjtotp = jtotp + spntot
            facjjp = (2.d0 * xjtotp + 1.d0) * facj
            do jlpp = jlpmin, jlpmax
              jlparp = 1 - (jlpp - 1)*2
*
* sum over row index for jtot
              do 400 irow = 1, length(jtot,jlp)
                if (is_twomol(ibasty)) then
                  j1_i = j(jtot,jlp,irow)/10
                  j2_i = mod(j(jtot,jlp,irow),10.d0)
                end if
                if (.not. is_twomol(ibasty) .and. is_j12(ibasty)) then
                  j1_i = j(jtot,jlp,irow)
                  if (ibasty.eq.23) then
                    j2_i = 0
                  else
                    j2_i = in(jtot,jlp,irow)
                  end if
                end if
                j12_i = jj12(jtot,jlp,irow)
                xj1_i = j1_i + spin
                xj2_i = j2_i + spnj2
                xj12_i = j12_i + spnj12
                l1 = l(jtot,jlp,irow)
                xl1 = l1
                sl1 = sqrt(2.d0 * xl1 + 1.d0)
                do 74 ij = 1, njmax
                  if (int(j(jtot,jlp,irow)) .ne. jslist(ij)) goto 74
                  if (int(in(jtot,jlp,irow)) .ne. inlist(ij)) goto 74
                  ir = ij
                  denrow = prefac(ij)
                  etrow = etrans(ij)
                  goto 75
74              continue
* level not in list
                goto 400
75              continue
*
* sum over column index for jtot
                do 200 icol = 1, length(jtot,jlp)
                  j12_f = jj12(jtot,jlp,icol)
                  xj12_f = j12_f + spnj12
                  l2 = l(jtot,jlp,icol)
                  xl2 = l2
                  sl2 = sqrt(2.d0 * xl2 + 1.d0)
* check for diagonal S matrix element
                  diag = irow .eq. icol
                  do 84 ij = 1, njmax
                    if (int(j(jtot,jlp,icol)) .ne. jslist(ij)) goto 84
                    if (int(in(jtot,jlp,icol)) .ne. inlist(ij)) goto 84
                    ic = ij
                    etcol = etrans(ij)
                    goto 85
84                continue
* level not in list
                  goto 200
85                continue
*
* sum over row index for jtotp
                  do 401 irowp = 1, length(jtotp,jlpp)
                    j12_ip = jj12(jtotp,jlpp,irowp)
                    if (j(jtotp,jlpp,irowp) .ne. j(jtot,jlp,irow))
     :                  goto 401
                    if (in(jtotp,jlpp,irowp) .ne. in(jtot,jlp,irow))
     :                  goto 401
                    if (j12_ip .ne. j12_i) goto 401
                    l1p = l(jtotp,jlpp,irowp)
                    xl1p = l1p
                    sl1p = sqrt(2.d0 * xl1p + 1.d0)
*
* sum over column index for jtotp
                    do 201 icolp = 1, length(jtotp,jlpp)
                      j12_fp = jj12(jtotp,jlpp,icolp)
                      if (j(jtotp,jlpp,icolp) .ne. j(jtot,jlp,icol))
     :                    goto 201
                      if (in(jtotp,jlpp,icolp) .ne. in(jtot,jlp,icol))
     :                    goto 201
                      if (j12_fp .ne. j12_f) goto 201
                      l2p = l(jtotp,jlpp,icolp)
                      xl2p = l2p
                      sl2p = sqrt(2.d0 * xl2p + 1.d0)
* check for diagonal S matrix element
                      diagp = irowp.eq.icolp
* get s-matrix elements
                      if (irow .gt. icol) then
                        ii = (irow*(irow - 1))/2 + icol
                      else
                        ii = (icol*(icol - 1))/2 + irow
                      end if
                      if (irowp .gt. icolp) then
                        iip = (irowp*(irowp - 1))/2 + icolp
                      else
                        iip = (icolp*(icolp - 1))/2 + irowp
                      end if
* phase factor
                      phase = 1.d0
                      ipower = nint(xj12_i - xj12_f
     :                    - 0.5d0*(xl1 - xl1p + xl2 - xl2p))
                      if ((ipower/2)*2 .ne. ipower) phase = -1.d0
* get T-matrix elements
                      t = - cmplx(sr(jtot,jlp,ii), si(jtot,jlp,ii))
                      if (diag)
     :                    t = (1.d0, 0.d0) + t
                      tp = -cmplx(sr(jtotp,jlpp,iip),si(jtotp,jlpp,iip))
                      if (diagp)
     :                     tp = (1.d0, 0.d0) + tp
                      factor = phase * facjjp * sl1 * sl2 * sl1p
     :                    * sl2p * real(t*conjg(tp))
     :                    / ((2.d0*xj1_i + 1.d0) * (2.d0*xj2_i + 1.d0))
*
* loop over k
                      do 100 k = 0, 2
                        xk = k
                        f36j(k) = 0.d0
                        x1 = xf3j(xl1,xl1p,xk,0.d0,0.d0,0.d0)
                        if (x1 .eq. 0.d0) goto 100
                        x2 = xf3j(xl2,xl2p,xk,0.d0,0.d0,0.d0)
                        if (x2 .eq. 0.d0) goto 100
                        x3 = xf6j(xl1,xjtot,xj1,xjtotp,xl1p,xk)
                        if (x3 .eq. 0.d0) goto 100
                        x4 = xf6j(xl2,xjtot,xj2,xjtotp,xl2p,xk)
                        if (x4 .eq. 0.d0) goto 100
                        f36j(k) = factor * x1*x2*x3*x4*(-1.d0)**k
                        if (iplmt.eq.1) plam(ic,ir,k+1)
     :                    = plam(ic,ir,k+1) + f36j(k)
100                   continue
                      epove = etcol/etrow
*
* N = 1
                      fac1 = f36j(0) - sqrt(epove)*f36j(1)
                      sigma(ic,ir,1) = sigma(ic,ir,1) + fac1*denrow
*
* N = 2
                      fac2 = onesix*(5.d0 - epove**2)*f36j(0)
     :                    - twothr*epove*f36j(2)
                      sigma(ic,ir,2) = sigma(ic,ir,2) + fac2*denrow
201                 continue
401               continue
200             continue
400           continue
            end do
          end do
        end do
      end do
*
* here when calculation has been done
999   do 1100 k = 1, 2
        write(2,800) k
        if(.not.batch) write(6,800) k
800     format(/' N =',i3)
        call mxoutd (2, sigma(1,1,k), njmax, njmax, 0, ipos)
* use diag as scratch variable (ipos = .false. for screen output)
        diag = .false.
        if(.not. batch)
     :    call mxoutd (6, sigma(1,1,k), njmax, njmax, 0, diag)
1100  continue
*
*


      if (iplmt.eq.1) then
      write(6,4411)
4411  format(/'PLAM:')
      do k = 0,2
      write(6,800)k
      diag=.false.
      call mxoutd (6, plam(1,1,k+1), njmax, njmax, 0, diag)
      enddo
      endif


*
*
*
      call mtime(cpu1,ela1)
      cpu1 = cpu1 - cpu0
      ela1 = ela1 - ela0
      call gettim(ela1,elaps)
      call gettim(cpu1,cpu)
      write(2,1200) elaps, cpu
      if(.not. batch) write(6,1200) elaps, cpu
1200  format(/' ** TIMING ELAPSED: ',a,
     :        ' CPU: ',a,/)
*
 1000 continue
 4007 deallocate(length)
 4014 deallocate(jj12)
 4006 deallocate(l)
 4005 deallocate(in)
 4004 deallocate(j)
 4003 deallocate(si)
 4002 deallocate(sr)
*
*
 4400 deallocate(plam)
*
*
 4001 deallocate(sigma)
 4000 continue
      if (ialloc .ne. 0) write (6, 4100)
 4100 format (' *** INSUFFICIENT MEMORY. ***')
      return
      end
c=====eof=====
