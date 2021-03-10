* ------------------------------------------------------------------
      subroutine tenopa(filnam,a)
*
* subroutine to calculate tensor opacities from s-matrix elements

* all references to B. Follmeg, P. Rosmus, and H.-J. Werner,
* J. Chem. Phys. 93, 4687 (1990)

* this subroutine returns the sigma(lambda,ki,kf,ji,jf) of Eq. 30
* of this paper, or, for lambda=0, sigma(ji,jf,K) which is Eq. 31
*
* author: b. follmeg
* current revision date:  4-may-1997 by mha
* ------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      character*(*) filnam
      character*40  tcsfil, smtfil, tcbfil
      character*20  cdate
      character*12  elaps, cpu
      logical csflag, flaghf, flagsu, twomol, exstfl, lpar,
     :        batch, fast, nucros,lpar2
*
      include "common/partens"
      include "common/parpot"
      common /colpar/ lpar(3), batch, lpar2(23)
      common /cospbf/ lnbufs, lnbufl, nbuf, maxlsp, maxllb, ihibuf,
     :                igjtp
      common /comom/  spin, xj1,xj2, j1, in1, j2, in2, maxjt, maxjot,
     :             nwaves, jfsts, jlparf, jlpars, njmax, j1min, j2max
      common /cosout/ nnout, jout(21)
      common /coisc2/ inlev(1)
      common /coisc3/ jlev(1)
      common /coisc4/ jpack(1)
      common /coisc5/ lpack(1)
      common /coisc6/ ipack(1)
      common /coisc7/ matel(1)
      common /coisc8/ jlist(1)
      common /cojq/ jq(1)
      common /colq/ lq(1)
      common /coinq/ inq(1)
      common /cosc1/  elev(1)
      common /cosc2/  prefac(1)
      common /cozmat/ jtotpa(1)
      common /coamat/ labadr(1)
      common /colnlb/ lenlab(1000)
      common /coz/    sreal(20)
      common /cow/    simag(20)
      common /cobmat/ sigma(kkmx*jmx*jmx)
      common /cotble/ npnt, jttble(1)
      common /ckli/  kplist(0:kmx)
      common /cf9a/ f9pha(kkmx)
      common /coconv/ econv, xmconv, ang2c
      common /codim/  nairy, mmax
      save nout
      dimension a(8)
      data  tol,   zero,   nstep
     :    /1.d-7,  0.d0,     2/
*
* initialyze timer
      srealp(lbufs)=zero
      simagp(lbufs)=zero
      ipackp(lbuflb)=0
      jpackp(lbuflb)=0
      lpackp(lbuflb)=0
      call mtime(cpu0,ela0)
      lnbufs = lbufs
      lnbufl = lbuflb
* input
      n = a(2)
* NB this variable is called lambda in Follmeg's paper and thesis
* it is the moment of the velocity distribution
      if(n.lt.0) then
         minn = iabs(n)
         maxn = minn
      else
         minn = 0
         maxn = n
      end if
      n = minn
      maxk=a(3)
      in1 = a(4)
*      if(in1.eq.0) in1=1
* previous statement deleted, mha 10/27/95
      in2 = a(5)
*      if(in2.eq.0) in2=1
* previous statement deleted, mha 10/27/95
      iener = a(1)
      maxjot = a(6)
      j1min = a(7)
      j2max = a(8)
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
      call openf(1,smtfil,'du',0)
*
* open file for tensor opacities
*
      call gennam(tcsfil,filnam,iener,'tcs',lenft)
      call openf(2,tcsfil,'sf',0)
      call gennam(tcbfil,filnam,iener,'tcb',lenfb)
      call openf(4,tcbfil,'su',0)
* rewind .tcb file
      rewind (4)
*
* read header of smatrix-file
*
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
      if (twomol) then
         write(2,12)
         if (.not. batch) write(6,12)
12       format(' *** TENSOR OPACITIES FOR MOLECULE -',
     :          ' MOLECULE COLLISIONS NOT YET IMPLEMENTED ***')
         goto 300
      end if
      if (flagsu) then
         write(2,14)
         if(.not. batch) write(6,14)
14       format(' *** TENSOR OPACITIES FOR SURFACE -',
     :          ' COLLISIONS NOT YET IMPLEMENTED ***')
         goto 300
      end if
      if (csflag) then
         write(2,16)
         if(.not. batch) write(6,16)
16       format(' *** CS TENSOR OPACITIES',
     :          ' NOT YET IMPLEMENTED ***')
         goto 300
      end if
*
      fast = .true.
      spin = 0.
      if(flaghf) then
         spin = 0.5d0
         fast = .false.
      end if
*
      write (2, 20) smtfil, cdate, label, potnam
      if(.not. batch) write (6, 20) smtfil, cdate, label,potnam
20    format(/' CLOSE COUPLED TENSOR OPACITIES',/,
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
* reset maxjt to jfinal if necessary
      if (maxjot.eq.0) maxjot = jfinal
      maxjt = min0(jfinal,maxjot)
      if (maxjt .lt. matjot) then
        write (2, 23) maxjot, jfinal
        if (.not. batch) write (6, 23) maxjot, jfinal
23      format (' MAX(JTOT) RESET TO JFINAL = ',i4, ' IN TENOPA')
      endif
      maxjot = maxjt
* nwaves is the number of partial waves
      nwaves = jfinal - jfirst + 1
* calculate prefactors and check if j' s are in jout list
* save pointer in array jlist
      nj = 0
      jmax = -1
      if (j1min.eq.0) j1min = jout(1)
      if (j2max.eq.0) j2max = jout(iabs(nnout))
      do 40 i=1, iabs(nout)
         jo = jout(i)
         if (jo .gt. j2max) goto 45
         if (jo .lt. j1min) goto 40
         do 30 j=1, nlevel
            j1 = jlev(j)
            in = inlev(j)
            if(j1.ne.jo .or. in.ne.in1 .and. in.ne.in2) goto 30
            if(j1 .gt. jmax) jmax = j1
* include only open channels
            if (ered .gt. elev(j)) then
              nj = nj + 1
              jlist(nj) = j
              prefac(j1+1) = ang2c * 3.1415926535897932d0 /
     :                  (2.d0* rmu * (ered - elev(j)))
              matel(j1+1) = nj
            endif
30      continue
40    continue
* check if there had been any match
45    if(nj.eq.0) then
         write(2,50)
         if(.not. batch) write(6,50)
50       format(' *** NO TRANSITIONS FOUND, ABORT ***')
         goto 300
      end if
*      maxk = 2 * jmax
*      maxk=0
      njmax = nj
      if(nj.gt.jmx) then
        write(6,51) jmx,nj
51      format(/' TENOPA: NJ =',i4, ' .GT. JMX = ',i4)
        return
      end if
* build up pointer table for direct access i/o
      call saddr(1,jfirst,jfinal,numin,numax,csflag,jttble,npnt,
     :           jfsts,jlparf,jlpars,ierr)
      if (ierr .ne. 0) goto 300
* read s-matrix for jfinal to find out the maximum buffer length
      jlp = 1
      if (jlparf .eq. 1) jlp = 0
      jj = jlp * nwaves + jfinal + 1
      iaddr = jttble(jj)
      nopen = -1
      call sread ( iaddr, sreal, simag, jtot, jlpar, nu,
     :            jq, lq, inq, ipack, jpack, lpack,
     :             1, mmax, nopen, length, ierr)
      maxlsp = (length*(length+1))/2
      maxllb = length
      nbuf1 = lnbufs/maxlsp
      nbuf2 = lnbufl/maxllb
      nbuf = min(nbuf1,nbuf2)
* store information on job.tcb file for later transformation into
* m-resolved cross sections
      write(4, err=999) label, cdate, ered, rmu, flaghf, nlevel,
     :                  nlevop, njmax, minn, maxn, nstep, maxk
      write(4 ,err=999) (jlev(i),i=1, nlevel)
      write(4 ,err=999) (inlev(i),i=1, nlevel)
      write(4 ,err=999) (elev(i),i=1, nlevel)
      write(4, err=999) (jlist(i),i=1, njmax)
* write header
      write(2,60) ered*econv,rmu*xmconv,jfirst,maxjt,maxk,minn,maxn,nbuf
      if(.not.batch) write(6,60)
     :  ered*econv,rmu*xmconv,jfirst,maxjt,maxk,minn,maxn,nbuf
60    format(/,' ENERGY: ',f11.3,' cm(-1)    MASS: ',f11.3,/,
     :     ' SUMMING PARTIAL WAVES FROM JTOT=',i3,' TO JTOT=',i3,/
     :     ' MAX(KI,KF) = ',i3,';',i3,' .LE. LAMDA .LE.',i3,
     :     '; MAXIMUM NUMBER OF BUFFERS =',i3)

* write level list
      write(2,65)
      if(.not. batch) write(6,65)
65    format (/,' ROWS ARE INITIAL STATES; COLUMNS ARE FINAL STATES',/,
     :        /,' LEVEL LIST FOR TENSOR OPACITIES (OPEN CHANNELS)',
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
* loop over n
      n = minn
* fast algorithm if n = 0
70    if (n.eq.0) then
         call sigk(maxk,nnout,jfirst,jfinal,jtotd,nj,mmax,jpack,
     :             lpack,ipack,prefac,sigma,
     :             sreal,simag,matel,lenlab,labadr,
     :             jtotpa,jttble,fast,ierr)
         if(ierr.ne.0) return
      else
* here for n > 0
         nk = 0
         xn = n
         do 75 k = 0, maxk
         xk = k
         do 75 kp = 0, maxk
         xkp = kp
         x = xf3jm0(xkp,xk,xn)
         if (abs(x) .lt. tol) goto 75
         nk = nk + 1
75       continue
         if(nk.gt.kkmx) then
           write(6,76) kkmx,nk
76         format(/' TENOPA: kkmx too small:',2i5)
           return
         end if
         call sigkkp(n,maxk,nk,nnout,jfirst,jfinal,jtotd,nj,mmax,
     :        jpack,lpack,ipack,prefac,sigma,
     :        sreal,simag,matel,lenlab,labadr,
     :        jtotpa,jttble,kplist,f9pha,fast,ierr)
         if(ierr.ne.0) return
      end if
* next n
      n = n + nstep
      if (n .le. maxn) goto 70
300   call mtime(cpu1,ela1)
      ela1 = ela1 - ela0
      cpu1 = cpu1 - cpu0
      call gettim(ela1,elaps)
      call gettim(cpu1,cpu)
      write(2,400) elaps, cpu
      if(.not. batch) write(6,400) elaps, cpu
400   format(/,' ** TENXSC FINAL TIMING, ELAPSED: ',a,'  CPU: ',a,' **')
      call closf(1)
      close (2)
      close (3)
      close (4)
      return
999   write(2,1000)
      if(.not.batch) write(6,1000)
1000  format(' *** I/O ERROR IN TENOPA; ABORT')
      call closf(1)
      close (2)
      close (3)
      close (4)
      return
      end
* ------------------------------------------------------------------
      subroutine addsp(jtmin,jtmax,jlp,
     :                 labadr,lenlab,jtotpa,jttble)
*
      implicit double precision (a-h,o-z)
      logical lpar, batch,lpar2
      include "common/partens"
      common /colpar/ lpar(3), batch,lpar2(23)
      common /comom/  spin, xj1,xj2, j1, in1, j2, in2, maxjt, maxjot,
     :             nwaves, jfsts, jlparf, jlpars, njmax
      common /cospbf/ lnbufs, lnbufl, nbuf, maxlsp, maxllb, ihibuf,
     :                igjtp
      common /cojq/ jqp(1)
      common /colq/ lqp(1)
      common /coinq/ inqp(1)
      dimension labadr(1),lenlab(1),jtotpa(1),jttble(1)
*
      ierr = 0
      ibuf = ihibuf
      jtpmin = max((igjtp+1),jtmin)
      jtpmax = jtmax
      if(jtmax-jtmin+1.gt.nbuf) then
        write(6,10) jtmax-jtmin+1,nbuf
10      format(/' NOT ENOUGH BUFFERS IN ADDSP:',2i5)
        stop
      end if
      do 100 jtp=jtpmin,jtpmax
* calculate address of s-matrix in file
         jjp = jlp * nwaves + jtp + 1
         iaddrp = jttble(jjp)
         if(iaddrp .lt. 0) goto 100
* calculate address of s-matrix in buffer
         ibuf = ibuf + 1
         if (ibuf .gt. nbuf) ibuf = 1
         ioffs = (ibuf-1)*maxlsp + 1
         if (ioffs .gt. 35525) write (6,11) jtp, ioffs
11       format (' *** JTP = ',i3, '; IOFFS = ', i8)

         ioff = (ibuf-1)*maxllb + 1
         if (ioff .gt. 1421) write (6,12) jtp, ioff
12       format (' *** JTP = ',i3, '; IOFF = ', i5)
* read s-matrix for jtot' into buffer
         nopenp = -1
         call sread ( iaddrp, srealp(ioffs), simagp(ioffs), jtotp,
     :                jlparp, nu, jqp, lqp, inqp, ipackp(ioff),
     :                jpackp(ioff), lpackp(ioff),
     :                1, maxlsp, nopenp, lengtp, ierr)
         if(ierr.eq.-1) goto 999
         if(ierr.lt.-1) then
            write(2,20)
            if(.not.batch) write(6,20)
20          format(' ** READ ERROR IN ADDSP, ABORT **')
            return
         end if
* save pointer for s-matrix and for label arrays
c        print*,'read jtp=',jtp,'  jtotp=',jtotp,
c    1          ' jlp=',jlp,'  jlpar=',jlparp
         jtotpa(jtotp+1) = ioffs - 1
         labadr(jtotp+1) = ioff - 1
         lenlab(jtotp+1) = lengtp
100   continue
      ihibuf = ibuf
      igjtp = jtotp
*
* here if end of file detected
*
999   return
      end
* ------------------------------------------------------------------
      subroutine mrcrs(filnam,a)
*
*
* author: b. follmeg
* current revision date: 18-apr-1997 by mha
*
* subroutine to compute integral fully m-resolved cross sections
* sigma(j1,j2,m1,m2) from tensorial cross sections sigma(n,k,k',j1,j2),
*
*                                    j1+j2-m1-m2                1/2
* sigma(lam,j1,j2,m1,m2) =  sum  (-1)            [(2k+1)(2k'+1)]
*                           k,k'
*
*                   ! j1  j1  k ! ! j2  j2  k'!
*                   !           ! !           ! sigma(lam,k,k',j1,j2)  ,
*                   ! m1 -m1  0 ! ! m2 -m2  0 !
*
*
* partially degeneracy averaged cross sections sigma(lam,j1,j2,m1),
*
* sigma(lam,j1,j2,m1) = sum sigma(lam,j1,j2,m1,m2) ,
*                      m2
*
*
* and fully degeneracy averaged cross sections sigma(lam,j1,j2),
*
*                          -1
* sigma(lam,j1,j2) = (2j1+1)  sum  sigma(lam,j1,j2,m1,m2) .
*                            m1,m2
*
*
* ------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      character*(*) filnam
      character*40  tcbfil, mcsfil
      character*20  cdate
      character*12  elaps, cpu
      logical flaghf, exstfl, lpar,
     :        batch,lpar2
*
      include "common/parpot"
      common /colpar/ lpar(3), batch,lpar2(23)
      common /comom/  spin, xj1,xj2, j1, in1, j2, in2, maxjt, maxjot,
     :             nwaves, jfsts, jlparf, jlpars, njmax
      common /coinhl/ jlev(1)
      common /coisc2/ inlev(1)
      common /coisc3/ jlist(1)
      common /cosc1/  elev(1)
      common /coz/    xm1lab(1)
      common /cow/    xm2lab(1)
      common /cozmat/ sigma(1)
      common /coamat/ sigmam(1)
      common /cobmat/ sigmak(1)
      common /coconv/ econv, xmconv, ang2c
      common /codim/  nairy, mmax
      data zero /0.d0/
*
      dimension a(1)
* initialyze timer
      call mtime(cpu0,ela0)
      zero =  zero
* input
      iener = a(1)
      if (iener .le. 0) iener = 1
*
* generate filename and check if it is present
*
      call gennam(tcbfil,filnam,iener,'tcb',lenft)
      inquire(file = tcbfil, exist = exstfl)
      if (.not. exstfl) then
          write(6,10) tcbfil(1:lenft)
10        format(' ** FILE ',(a),' NOT FOUND **')
          return
      end if
      call openf(1,tcbfil,'su',0)
      call gennam(mcsfil,filnam,iener,'mcs',lenm)
      call openf(2,mcsfil,'sf',0)
      read(1, err=999) label, cdate, ered, rmu, flaghf, nlevel,
     :                 nlevop, njmax, minn, maxn, nstep, maxk
      read(1 ,err=999) (jlev(i),i=1, nlevel)
      read(1 ,err=999) (inlev(i),i=1, nlevel)
      read(1 ,err=999) (elev(i),i=1, nlevel)
      read(1, err=999) (jlist(i),i=1, njmax)
*
      spin = 0.
      if(flaghf) spin = 0.5
*
      write (2, 20) tcbfil, cdate, label, maxk, maxk
      if(.not. batch) write (6, 20) tcbfil, cdate, label, maxk, maxk
20    format(/' CLOSE COUPLED M-RESOLVED CROSS SECTIONS',/,
     :        ' K K''-MATRICES READ FROM FILE ',(a),/,
     :        ' WRITTEN: ',(a),/,' LABEL:   ',(a),/,
     :        ' MAX(K, K'') IN PREVIOUS CALCULATION = ',i3,/,/,
     :        ' WARNING: M-DEPENDENCE CORRECT ONLY FOR J+J'' .LE.',i3,/,
     :        /,' ROWS ARE INITIAL STATES, COLUMNS ARE FINAL STATES',/,
     :        ' LAST COLUMN IS SUM OF THE ROW')
* NB the variable n, (minn .le. n .le. maxn in even steps)
* is called lambda in Follmeg's paper and thesis
* it is the moment of the velocity distribution

* loop over transitions
         do 400 i = 1, njmax
            ii = jlist(i)
            j1 = jlev(ii)
            xj1 = j1 + spin
            m1comp = nint(2.d0*xj1 + 1.d0)
            do 300 j = 1, njmax
               jj = jlist(j)
               j2 = jlev(jj)
               xj2 = j2 + spin
               m2comp = nint(2.d0*xj2 + 1.d0)
               call sigms(numk,i,j,m1comp,m2comp,minn,maxn,
     :                nstep,xm1lab,xm2lab,sigmak,sigmam,mmax,ierr)
               if (ierr .ne. 0) return
300         continue
400       continue
* obtain timing information
      call mtime(cpu1,ela1)
      cpu1 = cpu1 - cpu0
      ela1 = ela1 - ela0
      call gettim(ela1,elaps)
      call gettim(cpu1,cpu)
      write(2,600) elaps, cpu
      if(.not. batch) write(6,600) elaps, cpu
600   format(/,' ** MRCRS FINAL TIMING ELAPSED: ',
     :         (a),' CPU: ',(a),/)
      close(1)
      close(2)
      return
999   write(2,1000)
      write(6,1000)
1000  format(' ** READ ERROR IN MRCRS, ABORT')
      close(1)
      close(2)
      return
      end
* -------------------------------------------------------------------
      subroutine sigms(numk,ii,jj,m1comp,m2comp,minn,maxn,nstep,
     :                 xm1lab,xm2lab,sigmak,sigmam,mmax,ierr)
* ------------------
* current revision date: 18-apr-1997 by mha
* ------------------
      implicit double precision(a-h,o-z)
      logical lpar, batch, flaghf, lpar2
      character*20  cdate
*
      include "common/parpot"
      common /colpar/ lpar(3), batch,lpar2(23)
      common /comom/  spin, xj1,xj2, j1, in1, j2, in2, maxjt, maxjot,
     :             nwaves, jfsts, jlparf, jlpars, njmax
      common /coisc4/ isc1(1)
      common /coisc5/ isc2(1)
      common /coisc6/ isc3(1)
      common /coisc7/ isc4(1)
      common /cosc2/  sc1(1)
      dimension xm1lab(1),xm2lab(1),sigmak(mmax,1),sigmam(mmax,1)
      data tol /1.0d-10/
      data zero /0.0d0/
      ierr = 0
      rewind (1)
* dummy read
      read(1, err=999) label, cdate, ered, rmu, flaghf, nlevel,
     :                 id2, id3, id4, id5, id6, id7
      read(1 ,err=999) (isc1(i),i=1, nlevel)
      read(1 ,err=999) (isc2(i),i=1, nlevel)
      read(1 ,err=999) (sc1(i),i=1, nlevel)
      read(1, err=999) (isc3(i),i=1, id3)
* loop over n
      do 600 n = minn, maxn, nstep
* clear array
      do 30 i = 1, m1comp
      do 20 j = 1, m2comp+1
         sigmam(i,j) = zero
20    continue
30    continue
      sigmat = 0.
* loop over k and k'
      read(1, err=999, end=400) numk
      do 300 kk = 1, numk
         read(1, err = 999, end = 400) k,kp
         xk = k
         xkp = kp
         fack = sqrt(2.d0*xk+1)*sqrt(2.d0*xkp+1.d0)
* read matrix for k,k' values
         do 50 i = 1, njmax
50       read(1, err=999, end = 400) (sigmak(i,j), j = 1, njmax)
         trans = sigmak(ii,jj)
         if (abs(trans) .lt. tol) goto 300
* loop over components
         xm1 = -xj1
         do 200 m1 = 1, m1comp
            xm1lab(m1) = xm1
            xm2 = -xj2
            do 100 m2 = 1, m2comp
                xm2lab(m2) = xm2
* evaluate phase
                 phase = 1.d0
                 ipower = nint(xj1+xj2-xm1-xm2)
                 if((ipower/2)*2 .ne. ipower) phase = -1.d0
                 f3ja = xf3j(xj1,xj1,xk,xm1,-xm1,zero)
                 f3jb = xf3j(xj2,xj2,xkp,xm2,-xm2,zero)
                 s = phase * fack * f3ja * f3jb * trans
                 sigmam(m1,m2) = sigmam(m1,m2) + s
                 sigmam(m1,m2comp+1) = sigmam(m1,m2comp+1) + s
                 sigmat = sigmat + s
100         xm2 = xm2 + 1.d0
200      xm1 = xm1 + 1.d0
300   continue
* print out results
400   continue
      if (flaghf) then
         write(2, 410) xj1, xj2,n, sigmat/(2.d0*xj1+1.d0)
         if (.not. batch)
     :    write(6, 410) xj1, xj2, n, sigmat/(2.d0*xj1+1.d0)
410       format
     :      (/' TRANSITION J1 = ',f4.1,' -> J2 = ',f4.1,' LAM = ',i2,
     :       ' , TOTAL CROSS SECTION = ',1pe10.3,/)
      else
         write(2, 415) nint(xj1), nint(xj2),n, sigmat/(2.d0*xj1+1.d0)
         if (.not. batch)
     :   write(6, 415) nint(xj1), nint(xj2),n, sigmat/(2.d0*xj1+1.d0)
415      format
     :   (/' TRANSITION J1 = ',i2,' -> J2 = ',i2,', LAM = ',i2,
     :       ' , TOTAL CROSS SECTION = ',1pe10.3,/)
      endif
      lmax = 0
420   lmin = lmax + 1
      lmax = lmax + 9
      lmax = min0(lmax,m2comp)
      if (flaghf) then
        write(2,430) (xm2lab(l),l=lmin,lmax)
        if(.not. batch) write(6,430) (xm2lab(l),l=lmin,lmax)
430     format(10x,9(f5.1,6x))
      else
        write(2,435) (nint(xm2lab(l)),l=lmin,lmax)
        if(.not. batch) write(6,435) (nint(xm2lab(l)),l=lmin,lmax)
435     format(8x,9(i5,6x))
      endif
      do 500 m1 = 1, m1comp
          xm1 = xm1lab(m1)
          if (lmax .eq. m2comp) then
             if (flaghf) then
               write(2,440) xm1,(sigmam(m1,l),l=lmin,lmax),
     :                         sigmam(m1,m2comp+1)
               if (.not. batch) write(6,440)
     :                        xm1,(sigmam(m1,l),l=lmin,lmax),
     :                         sigmam(m1,m2comp+1)
440            format(1x,f4.1,12(1pe11.3))
             else
               write(2,445) nint(xm1),(sigmam(m1,l),l=lmin,lmax),
     :                         sigmam(m1,m2comp+1)
               if (.not. batch) write(6,445)
     :                        nint(xm1),(sigmam(m1,l),l=lmin,lmax),
     :                         sigmam(m1,m2comp+1)
445            format(1x,i4,12(1pe11.3))
             endif
          else
             if (flaghf) then
               write(2,440) xm1,(sigmam(m1,l),l=lmin,lmax)
               if (.not. batch) write(6,440)
     :                    xm1,(sigmam(m1,l),l=lmin,lmax)
             else
               write(2,445) nint(xm1),(sigmam(m1,l),l=lmin,lmax)
               if (.not. batch) write(6,445)
     :                    nint(xm1),(sigmam(m1,l),l=lmin,lmax)
             endif
          end if
500    continue
       write(2,'(a)') ' '
       if(.not. batch) write(6,'(a)') ' '
       if((m2comp-lmax)) 600,600,420
600    continue
       return
999   write(2,1000)
      write(6,1000)
1000  format(' ** READ ERROR IN SIGMS, ABORT')
      ierr = 1
      return
      end
*------- -----------------------------------------------------------------
      subroutine sigk(maxk,nnout,jfirst,jfinal,jtotd,nj,mmax,jpack,
     :                lpack,ipack,prefac,sigma,
     :                sreal,simag,matel,lenlab,labadr,
     :                jtotpa,jttble,fast,ierr)
*
* subroutine to calculate sigma(k,j1,j2) cross sections:
* ( see also " m.h. alexander and s.l. davis, jcp 78(11),6748(1983)"
*   equation 23 )
*
*         k          pi                                l1+l2-j1-j2+2j
*   sigma       =  -----      sum   (2j+1)(2j'+1)  (-1)
*         j1,j2    k(j1)^2  j,j',l1,l2
*
*             ! j1 j1 k ! ! j2 j2 k !    j             * j'
*             {         } {         }  t              t
*             ! j  j' l1! ! J  J' l2!    j1,l1,j2,l2    j1,l1,j2,l2
*
* author: b. follmeg
* current revision date: 18-apr-1997 by fdw
*
****** present version not correct for half-integer j (iadr)*************
*
*------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      complex*8 t, tp
      logical diag, diagj, diagin, lpar1, lpar2, batch, ipos,
     :        twopar, fast,lpar3
      character*12 elaps, cpu
      common /colpar/ lpar1(3), batch, lpar2(5), ipos,lpar3(17)
      common /comom/  spin, xj1,xj2, j1, in1, j2, in2, maxjt, maxjot,
     :             nwaves, jfsts, jlparf, jlpars, njmax, j1min, j2max
      common /cospbf/ lnbufs, lnbufl, nbuf, maxlsp, maxllb, ihibuf,
     :                igjtp
      include "common/partens"
      common/cadr/ iadr(0:jmx,lmx)
      common/c6jt/ f6a(kmx,0:jmx,lmx),f6p(kmx)
      common /cojq/ jq(1)
      common /colq/ lq(1)
      common /coinq/ inq(1)
      dimension jpack(1), lpack(1), ipack(1)
      dimension sreal(1), simag(1)
      dimension prefac(1), matel(1), labadr(1), jtotpa(1)
      dimension lenlab(1), jttble(1)
      dimension sigma(nj,nj,maxk+1)
* initialize timer
      call mtime(cpu0,ela0)
      ierr = 0
* clear sigma array
      do 6 k = 0, maxk
      do 5 j = 1, njmax
      do 4 i = 1, njmax
      sigma(i,j,k+1) = 0.d0
4     continue
5     continue
6     continue
      if(maxk.gt.kmx.or.j2max.gt.jmx) then
        write(6,7) maxk,kmx,j2max,jmx
7       format(/' SIGK: jmx or kmx too small:',4i5)
        ierr=-1
        return
      end if
      diagin = .false.
      if (in1. eq. in2) diagin = .true.
      twopar = .false.
      if(jlpars.ne.0) twopar = .true.
      jt = jfirst
      jstart = jfirst
      jlp = 1
      if (jlparf .eq. 1) jlp = 0
      jjoff = jlp * nwaves + 1
      ihibuf = 0
      igjtp = -1
* for integer jtot we need only even/even or odd/odd combinations of
* jtot/jtotp
      jtpstp = 1
      if (fast) jtpstp = 2
* sum over jtot
10    continue
* find address of jtot in smatrix file
      jj = jjoff + jt
      iaddr = jttble(jj)
      if(iaddr .lt. 0) goto 700
* read s-matrix for jtot
      nopen = -1
      call sread ( iaddr, sreal, simag, jtot, jlpar, nu,
     :            jq, lq, inq, ipack, jpack, lpack,
     :             1, mmax, nopen, length, ierr)
      if(ierr.eq.-1) goto 999
      if(ierr.lt.-1) then
        write(2,20)
        if(.not.batch) write(6,20)
20      format(' *** READ ERROR, ABORT')
        return
      end if
      xjtot = jtot + spin
      facj = (2.d0 * xjtot + 1.d0)
* bounderies for sum over jtot'
      jtpmin = jtot
      jtpmax = min((jtot + maxk), jfinal)
*
* fill buffer with required s' matrices
*
      call addsp(jtpmin,jtpmax,jlp,
     :           labadr,lenlab,jtotpa,jttble)
* sum over jtot'
      if (srealp(35525) .ne. 0.d0) print *, 'srealp error in sigk'
      if (simagp(35525) .ne. 0.d0) print *, 'simagp error in sigk'
      if (ipackp(1421) .ne. 0) print *, 'ipackp error in sigk'
      if (lpackp(1421) .ne. 0) print *, 'lpackp error in sigk'
      if (jpackp(1421) .ne. 0) print *, 'jpackp error in sigk'
      jtotp = jtpmin
60    continue
      switch = 1.d0
      if (jtotp .eq. jtot .or. jtotp .gt. maxjt) switch = 0.
* find address of jtot' in s-matrix buffer
      jj = jtotp + 1
      ioffs = jtotpa(jj)
      ioff  = labadr(jj)
      lengtp = lenlab(jj)
      xjtotp = jtotp + spin
      facjjp = (2.d0 * xjtotp + 1.d0) * facj
      jminjp = iabs(jtot-jtotp)
      jplujp = nint(xjtot + xjtotp)
      kmin=jminjp
      jpmax=jpackp(ioff+lengtp)
* now loop over all transitions
      lmax=jtotp+jpmax+1
      do 70 irowp = 1, lengtp
      irp = ioff + irowp
      if(ipackp(irp).ne.in1) goto 70
      j1p = jpackp(irp)
      l1p = lpackp(irp)
      iadr(j1p,lmax-l1p)=irowp
70    continue
      jmx1=0
      do 400 irow = 1, length
         if(ipack(irow).ne.in1) goto 400
         j1 = jpack(irow)
         if (j1 .gt. j2max ) goto 500
         if (j1 .gt. jpmax ) goto 500
         if (j1 .lt. jminjp) goto 400
         if (j1 .lt. j1min ) goto 400
         l1 = lpack(irow)
         if (l1.lt.abs(jtotp-j1).or.l1.gt.jtotp+j1) goto 400
         xj1 = j1 + spin
         j1t2 = nint(2.d0*xj1)
         xl1 = l1
         ir = matel(j1+1)
         denrow = prefac(j1+1)
         jmx1=max(j1,jmx1)
         kmax=min(jplujp,2*jmx1)
         kk=0
         l1m=lmax-l1
         do 80 k=kmin,kmax
         xk=k
         kk=kk+1
         if (kk.gt.kmx) print *, 'kk error in sigk'
         if (j1.gt.jmx) print *, 'j1 error in sigk'
         if (l1m.gt.lmx) print *, 'l1m error in sigk'
80       f6a(kk,j1,l1m) = xf6j(xj1,xj1,xk,xjtot,xjtotp,xl1)
         do 200 icol = 1, irow
            if(ipack(icol).ne.in2) goto 200
            j2 = jpack(icol)
            if (j2 .gt. j2max) goto 400
            if (j2 .gt. jpmax) goto 400
            if (j2 .lt. j1min) goto 200
            l2 = lpack(icol)
            if (l2.lt.abs(jtotp-j2).or.l2.gt.jtotp+j2) goto 200
            diagj = j1.eq.j2
            xj2 = j2 + spin
            j2t2 = nint(2.d0*xj2)
            min2j = min(j1t2,j2t2)
            if (min2j .lt. jminjp) goto 200
            xl2 = l2
            ic = matel(j2+1)
            dencol = prefac(j2+1)
            ii = (irow*(irow-1))/2 + icol
            l2m=lmax-l2
            irowp=iadr(j1,l1m)
            icolp=iadr(j2,l2m)
            if(irowp.ge.icolp) then
              iip = ioffs + (irowp*(irowp-1))/2 + icolp
            else
              iip = ioffs + (icolp*(icolp-1))/2 + irowp
            end if
            if (diagj) dencol = 0.
            diag = diagj .and. diagin
            wd = 1.d0
            if (diag .and. l1.ne.l2) wd = 2.d0
* phase factor
             phasea = 1.d0
             ipower= nint(xl1+xl2-xj1-xj2+2.d0*xjtot)
             if ((ipower/2)*2 .ne. ipower) phasea = -1.d0
             phaseb = 1.d0
             ipower= ipower-2*nint(xjtotp-xjtot)
             if ((ipower/2)*2 .ne. ipower) phaseb = -1.d0
*
* convert s-matrix to t-matrix: t(j,j1,in1,l1 ; j2,in2,l2) =
*     delta(j1,in1,l1 ; j2,in2,l2) - s(j,j1,in1,l1 ; j2,in2,l2)
*
             t  = -cmplx(sreal(ii),simag(ii))
             tp = -cmplx(srealp(iip),simagp(iip))
         if (iip .gt. 35525) write (6,11)iip
11       format (' *** IIP = ', i8)
             if (diag .and. l1.eq.l2) then
                 t = t + (1.d0,0)
                 tp = tp + (1.d0,0)
             end if
             factor = wd * facjjp * real(t*conjg(tp))
             factor = factor * (phasea + phaseb*switch)
* loop over k
            kmax = min(min2j,jplujp)
            kk=0
            do 90 k = kmin, kmax
               kk=kk+1
90             f6p(kk) = f6a(kk,j1,l1m)*f6a(kk,j2,l2m)*factor
            kk=0
            do 100 k = kmin, kmax
               kk=kk+1
               sigma(ir,ic,k+1) = sigma(ir,ic,k+1) + f6p(kk)*denrow
               sigma(ic,ir,k+1) = sigma(ic,ir,k+1) + f6p(kk)*dencol
100         continue
200      continue
400   continue
* next jtot'
500   jtotp = jtotp + jtpstp
      if(jtotp.le.jtpmax) goto 60
* next jtot
700   jt = jt + 1
      if (jt.le.maxjt) goto 10
* next parity
      if (twopar) then
          twopar = .false.
          jt = jfsts
          jstart = jfsts
          jlp = 1
          jjoff = jlp * nwaves + 1
          ihibuf = 0
          igjtp = -1
          goto 10
      end if
* here if calculation has been done
999   write(4) maxk + 1
      do 1100 k = 0, maxk
      write(4) k,k
      write(2,800) k
      if(.not.batch) write(6,800) k
800   format(/' TENSOR RANK K =',i3)
      do 1000 i = 1, njmax
1000  write(4) (sigma(i,j,k+1),j=1,njmax)
      call mxoutd (2, sigma(1,1,k+1), njmax, njmax, 0, ipos)
* use diag as scratch variable (ipos = .false. for screen output)
      diag = .false.
      if(.not. batch)
     :   call mxoutd (6, sigma(1,1,k+1), njmax, njmax, 0, diag)
1100  continue
      call mtime(cpu1,ela1)
      cpu1 = cpu1 - cpu0
      ela1 = ela1 - ela0
      call gettim(ela1,elaps)
      call gettim(cpu1,cpu)
      write(2,1200) elaps, cpu
      if(.not. batch) write(6,1200) elaps, cpu
1200  format(/' ** N = 0 COMPLETED, TIMING ELAPSED: ',a,
     :        ' CPU: ',a,/)
      return
      end
*------------------------------------------------------------------------
      subroutine sigkkp(n,maxk,nk,nnout,jfirst,jfinal,jtotd,nj,mmax,
     :         jpack,lpack,ipack,prefac,sigma,
     :         sreal,simag,matel,lenlab,labadr,
     :         jtotpa,jttble,kplist,f9pha,fast,ierr)
*
* author: b. follmeg
* current revision date: 18-apr-1997 by mha
*
*------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      complex*8 t, tp, ai, cphase
      logical diag, diagj, diagin, lpar1, lpar2, batch, ipos,
     :        twopar, fast,lpar3
      character*12 elaps, cpu
      common /colpar/ lpar1(3), batch, lpar2(5), ipos,lpar3(17)
      common /comom/  spin, xj1,xj2, j1, in1, j2, in2, maxjt, maxjot,
     :             nwaves, jfsts, jlparf, jlpars, njmax, j1min, j2max
      common /cospbf/ lnbufs, lnbufl, nbuf, maxlsp, maxllb, ihibuf,
     :                igjtp
      include "common/partens"
      common/cccpu/ tchi,tlog,tsetup,tdelt,lenk,max1,max2,max3,maxkk
      common/cadr/ iadr(0:jmx,lmx)
      common/c6jt/ f6a(kmx,0:jmx,lmx),f9a(3*kmx,lmx)
      common /cojq/ jq(1)
      common /colq/ lq(1)
      common /coinq/ inq(1)
      dimension jpack(1), lpack(1), ipack(1)
      dimension sreal(1), simag(1)
      dimension prefac(1), matel(1), labadr(1), jtotpa(1)
      dimension lenlab(1), jttble(1), kplist(0:maxk)
      dimension sigma(nk,nj,nj), f9pha(nk)
      data tol /1.d-10/
      data zero /0.d0/
cstart unix
      data ai / (0.d0,1.d0)/
cend
cstart cray
c;      data ai / (0.d0,1.d0)/
cend

* initialize timer
      call mtime(cpu0,ela0)
      t6j=0
      t9j=0
      ierr = 0
      icount = 0
      xn = n
* clear sigma array
      if(maxk.gt.kmx-1) then
        write(6,4) kmx,maxk
4       format(/' SIGKKP: KMX too small:',i4,'  need:',i4)
        ierr=-1
        return
      end if
      if(j2max.gt.jmx) then
        write(6,5) jmx,j2max
5       format(/' SIGKKP: JMX too small:',i4,'  need:',i4)
        ierr=-1
        return
      end if
      nkk=0
      iof=0
      do 6 kp=0,maxk
      kplist(kp)=nkk
      do 6 k=abs(kp-n),min(maxk,kp+n),n
      nkk=nkk+1
      if (nkk .gt. nk) print *, 'nkk error',kp, k, nk

      do 6 i=1,nj
      do 6 j=1,nj
6     sigma(nkk,j,i)=zero

*      if(nkk.ne.nk) then
*        print *, maxk, n,k,kp,nkk
*        write(6,7) nk,nkk
*7       format(/' nk ne nkk:',2i4)
*        ierr=-1
*        return
*      end if
      if(nkk.gt.kmx*3) then
        write(6,8) kmx**2,nkk
8       format(/' SIGKKP: KMX**2 too small:',i4,'  need:',i4)
        ierr=-1
        return
      end if
      diagin = .false.
      if (in1. eq. in2) diagin = .true.
      twopar = .false.
      if(jlpars.ne.0) twopar = .true.
      jt = jfirst
      jstart = jfirst
      jlp = 1
      if (jlparf .eq. 1) jlp = 0
      jjoff = jlp * nwaves + 1
      ihibuf = 0
      igjtp = -1
* for integer jtot we need only even/even or odd/odd combinations of
* jtot/jtotp
      jtpstp = 1
      if (fast) jtpstp = 2
* sum over jtot
10    continue
* find address of jtot in smatrix file
      jj = jjoff + jt
      iaddr = jttble(jj)
      if(iaddr .lt. 0) goto 700
* read s-matrix for jtot
      nopen = -1
      call sread ( iaddr, sreal, simag, jtot, jlpar, nu,
     :            jq, lq, inq, ipack, jpack, lpack,
     :             1, mmax, nopen, length, ierr)
      if(ierr.eq.-1) goto 999
      if(ierr.lt.-1) then
        write(2,20)
        if(.not.batch) write(6,20)
20      format(' *** READ ERROR IN SIGKKP, ABORT')
        return
      end if
      xjtot = jtot + spin
      facj = (2.d0* xjtot + 1.d0)
* bounderies for sum over jtot'
      jtpmin = jtot
      jtpmax = min((jtot + maxk), jfinal)
* fill buffer with required s' matrices
      call addsp(jtpmin,jtpmax,jlp,
     :           labadr,lenlab,jtotpa,jttble)
* sum over jtot'
      jtotp = jtpmin
60    continue
      switch = 1.d0
      if (jtot .eq. jtotp .or. jtotp .gt. maxjt) switch = 0.
* find address of jtot' in s-matrix buffer
      jj = jtotp + 1
      ioffs = jtotpa(jj)
      ioff  = labadr(jj)
      lengtp = lenlab(jj)
      xjtotp = jtotp + spin
      facjjp = (2.d0* xjtotp + 1.d0) * facj
      jminjp = iabs(jtot-jtotp)
      jplujp = nint(xjtot + xjtotp)
* precompute 6j symbols
      j2mxp=jpackp(ioff+lengtp)
      j2mx =jpack(length)
      lmax =jtot+j2mx+1
      lmaxp=jtotp+j2mxp+1
      kpmin=jminjp
      kmin=abs(kpmin-n)
      t1=second()
      do 70 icol=1,length
      if(ipack(icol).ne.in2) goto 70
      j2=jpack(icol)
      l2=lpack(icol)
      xj2=j2+spin
      xl2=l2
      kpmax=min(jplujp,nint(2*xj2))
      kkp=0
      do 65 kp=kpmin,kpmax
      kkp=kkp+1
      xkp=kp
65    f6a(kkp,j2,lmax-l2) = xf6j(xkp,xjtotp,xjtot,xl2,xj2,xj2)
70    continue
      t6j=t6j+second()-t1
      kpmx=kpmax
      do 80 icolp=1,lengtp
      irp = ioff + icolp
      if(ipackp(irp).ne.in1) goto 80
      j1=jpackp(irp)
      l1=lpackp(irp)
      iadr(j1,lmaxp-l1)=icolp
80    continue
* now loop over all transitions
      do 400 irow = 1, length
         if(ipack(irow).ne.in1) goto 400
         j1 = jpack(irow)
         if (j1 .gt. j2max) goto 500
         if (j1 .gt. j2mxp) goto 500
         if (j1 .lt. j1min) goto 400
         l1 = lpack(irow)
         l1pmin=max(abs(l1-n),abs(j1-jtotp))
         if((-1)**(l1pmin+j1-jtotp).ne.jlpar) l1pmin=l1pmin+1
         l1pmax=min(l1+n,j1+jtotp)
         if(l1pmax.lt.l1pmin) goto 400
         xj1 = j1 + spin
         xl1 = l1
         j1t2 = nint(2.d0*xj1)
         ir = matel(j1+1)
         xxj=xjtot+xj1+xjtotp+xj1
         denrow = prefac(j1+1)
         kmax=nint(2*xj1)
         factor=denrow*facjjp
         ll=0
         t1=second()
      do 100 l1p=l1pmin,l1pmax,2
         ll=ll+1
         xl1p = l1p
         cphase = 1.d0
         ipower = l1 + l1p
         if (ipower .ne. 0) cphase = ai ** ipower
         f3a = real(cphase)*xf3jm0(xl1,xl1p,xn)
     1      *sqrt((2.d0*xl1+1.d0)*(2.d0*xl1p+1.d0))*factor
* precalculate 9j symbols and phase factors
         do 90 kp=kpmin,kpmx
         kmin1=abs(kp-n)
         kmax1=min(kp+n,kmax)
         ik=kplist(kp)
         xkp = kp
         do 90 k=kmin1,kmax1,n
         xk = k
         ik=ik+1
         f9a(ik,ll) =  f3a*(-1)**(k+kp)
     :            * sqrt((2.d0*xk+1.d0)*(2.d0*xkp+1.d0))
     :            * xf3jm0(xkp,xk,xn)
c     :            * xf9j(xkp,xk,xn,xjtot,xj1,xl1,xjtotp,xj1,xl1p)
     :             * xf9j(xn,xkp,xk,xl1,xjtot,xj1,xl1p,xjtotp,xj1)
         f9pha(ik) = switch
         ipower = nint(xk+xkp+xxj)
         if(mod(ipower,2).ne.0) f9pha(ik) = -switch
90    continue
100   continue
         t9j=t9j+second()-t1
      do 300 icol = 1, length
         if(ipack(icol).ne.in2) goto 300
         j2 = jpack(icol)
         if (j2 .gt. j2max) goto 400
         if (j2. gt. j2mxp) goto 400
         if (j2. lt. j1min) goto 300
         l2 = lpack(icol)
         if(l2.gt.jtotp+j2.or.l2.lt.abs(jtotp-j2)) goto 300
         xj2 = j2 + spin
         xl2 = l2
         kpmax=min(jplujp,nint(2*xj2))
         isrow = max(irow,icol)
         iscol = min(irow,icol)
         ii = (isrow * (isrow-1)) / 2 + iscol
         t  = -cmplx(sreal(ii),simag(ii))
         diagj = j1 .eq. j2
         diag = diagj .and. diagin
         if (diag .and. l1.eq.l2) t = t + (1.d0,0.d0)
         l2m=lmax-l2
         icolp=iadr(j2,lmaxp-l2)
         if(icolp.eq.0) goto 300
         potenz = ((2.d0* xj1) - xj2 + xl2)
         ipower = nint(xjtotp+potenz)
         ifaka=1
         if(mod(ipower,2).ne.0) ifaka=-1
         ipower = nint(xjtot+potenz)
         ifakb=1
         if(mod(ipower,2).ne.0) ifakb=-1
         ic = matel(j2+1)
         ll=0
      do 200 l1p=l1pmin,l1pmax,2
         ll=ll+1
         irowp=iadr(j1,lmaxp-l1p)
         if(irowp.eq.0) goto 200
         isrowp = max(irowp,icolp)
         iscolp = min(irowp,icolp)
         iip = ioffs + (isrowp * (isrowp-1)) / 2 + iscolp
*
* convert s-matrix to t-matrix: t(j,j1,in1,l1 ; j2,in2,l2) =
*     delta(j1,in1,l1 ; j2,in2,l2) - s(j,j1,in1,l1 ; j2,in2,l2)
*
         if (iip .gt. 35525) write (6,11)iip
11       format (' *** IIP = ', i8)
         tp = -cmplx(srealp(iip),simagp(iip))
         if (diag .and. l1p.eq.l2) tp = tp + (1.d0,0d0)
         fac = real(t*conjg(tp))
         factra = fac*ifaka
         factrb = fac*ifakb
* loop over k
         kkp=0
c      print*,'jtot,jtotp,j1,l1,j2,l2,l1p,icolp,irowp',
c    1         jtot,jtotp,j1,l1,j2,l2,l1p,icolp,irowp
       do 160 kp=kpmin,kpmax
         kmin1=abs(kp-n)
         kmax1=min(kp+n,kmax)
         ikk=kplist(kp)
         kkp=kkp+1
         fac=f6a(kkp,j2,l2m)
         faca=factra*fac
         facb=factrb*fac
       do 150 k=kmin1,kmax1,n
         ikk=ikk+1
         sigma(ikk,ir,ic)=sigma(ikk,ir,ic)
     1       +f9a(ikk,ll)*(faca+f9pha(ikk)*facb)
150   continue
160   continue
200   continue
300   continue
400   continue
* next jtot'
500   jtotp = jtotp + jtpstp
      if(jtotp.le.jtpmax) goto 60
* next jtot
700   jt = jt + 1
      if (jt.le.maxjt) goto 10
* next parity
      if (twopar) then
          twopar = .false.
          jt = jfsts
          jstart = jfsts
          jlp = 1
          jjoff = jlp * nwaves + 1
          ihibuf = 0
          igjtp = -1
          goto 10
      end if
* here if calculation has been done
999   write(4) nk
      ikk=0
      do 1000 kp= 0, maxk
      do 1000 k = abs(kp-n),min(maxk,kp+n),n
      ikk=ikk+1
      write(4) k,kp
      write(2,800) n,k,kp
      if(.not.batch) write(6,800) n,k,kp
800   format(/' LAMBDA =',i2,'; TENSOR RANK KI =',i3,' KF =',i3/)
      do 1000 i = 1, njmax
      write(2,1010) (sigma(ikk,i,j),j=1,njmax)
      if(.not. batch) write(6,1010) (sigma(ikk,i,j),j=1,njmax)
1000  write(4) (sigma(ikk,i,j),j=1,njmax)
1010  format(1x,10(1pd12.4))
* use diag as scratch variable (ipos = .false. for screen output)
      diag = .false.
      call mtime(cpu1,ela1)
      cpu1 = cpu1 - cpu0
      ela1 = ela1 - ela0
      call gettim(ela1,elaps)
      call gettim(cpu1,cpu)
*	      write(2,1200) n,elaps, cpu
*      if(.not. batch) write(6,1200) n,elaps, cpu
*1200  format(/' ** N =',i2,' COMPLETED, TIMING ELAPSED: ',a,
*     :        ' CPU: ',a)
*	      write(2,1200) n,elaps, cpu,t6j,t9j,
*     :      lenk,max1,max2,max3,maxkk,tsetup,tdelt,tchi,tlog
*      if(.not. batch) write(6,1200) n,elaps, cpu,t6j,t9j,
*     :      lenk,max1,max2,max3,maxkk,tsetup,tdelt,tchi,tlog
*1200  format(/' ** N =',i2,' COMPLETED, TIMING ELAPSED: ',a,
*     :        ' CPU: ',a/
*     :        ' 6J Symbols:',f10.2,'  9J Symbols:',f10.2/
*     : ' LK=',i4,'  L1=',i4,'  L2=',i4,'  L3=',i4,'  LKK',i4/
*     : '  SETUP:',f10.2,'  TDELT:',f10.2,'  TCHI:',f10.2,
*     : '  TLOG:',f10.2)
      return
      end
