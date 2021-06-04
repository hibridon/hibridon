* -------------------------------------------------------------------------
      subroutine stmix(flnam1, flnam2, a)
*
*  subroutine to compute cross sections for mixed singlet-triplet states
*    e.g. CH3 X3B1 - a1A1 collisions involving perturbed levels
*
*  this subroutine requires two smt files, one for the singlet level and
*     one for the triplet level
*
*  author:  p. dagdigian
*  revision:  9-may-2011
*  changed open command to open smt files for streaming reads (jan-2012,
*     q. ma)
*  fixed typo in iflag name (one line had flag instead of iflag)
*     (jun-2013, p. dagdigian)
*  allocate memory dynamically (jun-2013, q. ma)
*  fixed problem with calculation of endoergic transitions involving
*     unperturbed levels (21-jun-2013)
*
*  current revision:  21-jun-2013 by p. dagdigian
* -------------------------------------------------------------------------
      implicit double precision (a-h, o-z)
      character*(*) flnam1, flnam2
      character*20  cdate1, cdate2
      character*40  smtfil1, smtfil2
      character*10  elaps, cpu
      character*1 slab
      logical csflg1, flghf1, flgsu1, twmol1, nucrs1,
     :     csflg2, flghf2, flgsu2, twmol2, nucrs2,
     :     batch, fast, lpar2, lpar, exstfl, sngsmt, trpsmt
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
      common /coconv/ econv, xmconv, ang2
      common /codim/ nairy, mmax
      dimension a(8)
      dimension epert(2), cpert(2,2)
*     scratch arrays for diagonalization subroutine
      parameter (narray = 2)
      dimension e(narray,narray), eig(narray), vec(narray,narray),
     :     sc1(narray), sc2(narray), work(144)
*     cross section array  -- first 2 columns will contain cross sections
*     for transitions out of the 2 perturbed levels,the second 2 columns
*     will contain cross sections for transitions into these levels
      double precision, dimension(:, :), allocatable :: sigma
      integer, dimension(:), allocatable :: iptsng, ipttrp, jq, lq, inq
      double precision, dimension(:), allocatable :: sreal, simag
*     storage for s-matrix elements:
*     second letter is real (r), imaginary (i)
*     third letter is singlet (s), triplet (t)
*     subscripts:  jtot, jlp (= 1/2 for jlpar = +1/-1), size of channel basis array
      double precision, dimension(:, :, :), allocatable ::
     $     srs, sis, srt, sit
*     arrays with values of j, in, and l:
*     subscripts:  jtot, jlp (= 1/2 for jlpar = +1/-1), length of channel basis
      integer, dimension(:, :, :), allocatable ::
     $     js, ins, ls, jt, intt, lt
*     length of arrays
*     last letter in name:  singlet (s), triplet (t)
*     subscripts:  jtot, jlp (= 1/2 for jlpar = +1/-1)
      integer, dimension(:, :), allocatable :: lngths, lngtht
*     singlet or triplet (2nd index) - for positive (p) and (n) negative parity
      logical, dimension(:, :), allocatable :: exsmtp, exsmtn
*
*     initialize timer and arrays
      call mtime(cpu0,ela0)
*
*     generate filename of 1st smt file and check if it is present
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
*     open 1st s-matrix file
*
      call openf(1,smtfil1,'tu',0)
*
*     initialize timer and arays
      call mtime(cpu0,ela0)
*
*     generate filename of 2nd smt file and check if it is present
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
*     open 2nd s-matrix file
*
      call openf(11,smtfil2,'tu',0)
*
*     read header of 1st s-matrix file
*
      call sinqr(1, mjtots, mmaxs)
      call rdhead(1,cdate1,ered1,rmu1,csflg1,flghf1,flgsu1,
     :     twmol1,nucrs1,jfrst1,jfinl1,jtotd1,numin1,numax1,nud1,
     :     nlevel1,nlvop1,nnout1,jlev1,inlev1,elev1,jout1)
      mmax2s = mmaxs * (mmaxs + 1) / 2
*
*     we need the s-matrices as lower triangles, so nnout1  m u s t  be > 0
*
      if (nnout1.lt.0) then
         write(6,11)
 11      format(/' ** NNOUT < 0, ABORT **'/)
         goto 4000
      end if
      nout = nnout
*
*     read header of 2nd s-matrix file
*
      call sinqr(11, mjtott, mmaxt)
      call rdhead(11,cdate2,ered2,rmu2,csflg2,flghf2,flgsu2,
     :     twmol2,nucrs2,jfrst2,jfinl2,jtotd2,numin2,numax2,nud2,
     :     nlevel2,nlvop2,nnout2,jlev2,inlev2,elev2,jout2)
      mmax2t = mmaxt * (mmaxt + 1) / 2
*
      mmaxst = max(mmaxs, mmaxt)
*
*     we need the s-matrices as lower triangles, so nnout1  m u s t  be > 0
*
      if (nnout2.lt.0) then
         write(6,11)
         goto 4000
      end if
      nout = nnout
*
*     molecule-molecule cross sections not implemented
*
      if (twmol1. or. twmol2) then
         write(6,12)
 12      format(/' *** SINGLET-TRIPLET MIXING FOR MOLECULE -',
     :        ' MOLECULE COLLISIONS NOT IMPLEMENTED ***'/)
         goto 4000
      end if
*
*     molecule-surface collisions not implemented
*
      if (flgsu1 .or. flgsu2) then
         write(6,14)
 14      format(/' *** SINGLET-TRIPLET MIXING FOR SURFACE -',
     :        ' COLLISIONS NOT IMPLEMENTED ***'/)
         goto 4000
      end if
*
*     CS cross sections not implemented
      if (csflg1 .or. csflg2) then
         write(6,16)
 16      format(/' *** CS SINGLET-TRIPLET CROSS SECTIONS',
     :        ' NOT IMPLEMENTED ***'/)
         goto 4000
      end if
*
*     delta-jtot should be equal to one for both s-matrix files
      if (iabs(jtotd1).ne.1 .or. iabs(jtotd2).ne.1)
     :     then
         write(6,116)
 116     format(/' *** DELTA-JTOT MUST BE EQUAL TO ONE',
     :        ' *** '/)
         goto 4000
      end if
*
*     electron spin cannot be half-integral
*
      if (flghf2 .or. flghf2) then
         write(6,16)
 18      format(/' *** ELECTRON SPIN CANNOT BE',
     :        ' HALF-INTEGRAL ***'/)
         goto 4000
      end if
*
*     reduced mass in the two smt files must be the same
*
      if (rmu1 .ne. rmu2) then
         write(6,16)
 19      format(/' *** REDUCED MASS IN BOTH SMT FILES',
     :        ' MUST BE EQUAL ***'/)
         goto 4000
      end if
*
*     difference of energies in the two smt files must equal DELE
*     (to within 0.1 cm^-1)
*
      dele = a(3)
      emax = a(4)
      istats = a(5)
      istatt = a(6)
      wso = a(7)
      diffe = (ered1 - ered2) * econv
      if (abs(dele - diffe) .gt. 0.1d0) then
         write(2,16)
         if(.not. batch) write(6,21) diffe, dele
 21      format(/' *** ENERGY DIFFERENCE OF SMT FILES',
     :        ' NOT EQUAL TO DELE ***'/
     :        5x,'DIFFE = ',f11.3,4x,'DELE ',f11.3/)
         goto 4000
      end if
      write (6, 22) label, potnam
 22   format(/' SINGLET-TRIPLET MIXING CROSS SECTIONS'/
     :     '      LABEL:     ',(a)/
     :     '      POT NAME:  ',(a))
      write (6, 24) smtfil1, cdate1, smtfil2, cdate2
 24   format(/' S-MATRICES FOR SINGLET READ FROM FILE ',(a),/,
     :     '      WRITTEN:   ',(a),//
     :     ' S-MATRICES FOR TRIPLET READ FROM FILE ',(a),/,
     :     '      WRITTEN:   ',(a))
      write(6,26) dele
 26   format(/' ENERGY OF TRIPLET MANIFOLD RELATIVE TO',
     :     ' SINGLET:',f10.3)

      write (6,38) emax
 38   format (/' LEVELS WITH ENERGY .LT. ',f9.2,
     :     ' CM(-1):')
*
*     print list of levels in singlet state
      write (6, 40)
 40   format (/,5x,'SINGLET LEVELS'/2x,
     :     'LEVEL LIST SORTED BY ENERGY',/,'   #   J  ',
     :     'IS  KP  KO  S   EINT(CM-1)')
      nlevels = 0
*
      allocate(iptsng(3 * nlevel1), stat=ialloc)
      if (ialloc .ne. 0) goto 4000
*
      do i=1, nlevel1
         if (elev1(i)*econv.le.emax) then
            kp = iabs(inlev1(i))
            if (inlev1(i) .ge. 0) then
               ko = jlev1(i) - kp
            else
               ko = jlev1(i) + 1 - kp
            end if
            is = sign(1, inlev1(i))
            if (is.eq.1) then
               slab='+'
            else
               slab='-'
            end if
            ecm = elev1(i) * econv
            nlevels = nlevels + 1
*     iptsng points to levels in nlevel1 list
            iptsng(nlevels) = i
            write (6, 1351) nlevels, jlev1(i), inlev1(i),
     :           kp, ko, slab, ecm
 1351       format (5i4, a3, f10.3)
         end if
      end do
*
*     print list of levels in triplet state
      write (6, 55)
 55   format (/,5x,'TRIPLET LEVELS'/2x,
     :     'LEVEL LIST SORTED BY ENERGY',/,'   #   N',
     :     '   J  IS  KP  KO  S   EINT(CM-1)')
      nlevelt = 0
*
      allocate(ipttrp(3 * nlevel2), stat=ialloc)
      if (ialloc .ne. 0) goto 4001
*
      do i=1, nlevel2
         ecm = elev2(i) * econv + dele
         if (ecm.le.emax) then
            kp = iabs(inlev2(i))
            if (inlev2(i) .ge. 0) then
               ko = jlev2(i) - kp
            else
               ko = jlev2(i) + 1 - kp
            end if
            is = sign(1, inlev2(i))
            if (is.eq.1) then
               slab='+'
            else
               slab='-'
            end if
*     shift energies of triplet levels
            jmin = abs(jlev2(i) - 1)
            jmax = jlev2(i) + 1
            do jtrip = jmin, jmax
               nlevelt = nlevelt + 1
*     ipttrp points to levels in nlevel2 list
               ipttrp(nlevelt) = i
               jlevt(nlevelt) = jtrip
               nlevt(nlevelt) = jlev2(i)
               inlevt(nlevelt) = inlev2(i)
               elevt(nlevelt) = ecm / econv
               write (6, 1361) nlevelt, jlev2(i), jtrip, inlev2(i),
     :              kp, ko, slab, ecm
 1361          format (6i4, a3, f10.3)
            end do
         end if
      end do
*
      if (istats.gt.nlevel1 .or. istats.le.0) then
         write(6,62) istats
 62      format(/' PERTURBED SINGLET STATE INDEX OUT OF BOUNDS'/
     :        5x,'ISTATS =',i4/)
         goto 4002
      end if
      if (istatt.gt.nlevelt .or. istatt.le.0) then
         write(6,64) istatt
 64      format(/' PERTURBED TRIPLET STATE INDEX OUT OF BOUNDS'/
     :        5x,'ISTATT = ',i4/)
         goto 4002
      end if
*
      if (wso .eq. 0.d0) then
         write(6,1369)
 1369    format(/' WSO .eq. 0.  DO NOT COMPUTE CROSS SECTIONS'/)
         goto 4002
      end if
*
* make sure J for the singlet and triplet levels are equal
      if (jlev1(istats) .ne. jlevt(istatt)) then
         write(6,61) istats, jlev1(istats), inlev1(istats),
     :        istatt, jlevt(istatt), nlevt(istatt), inlevt(istatt)
 61      format(/,' *** TOTAL ANGULAR MOMENTUM J OF',
     :        ' LEVELS TO BE MIXED NOT EQUAL ***'/
     :        5x,'SINGLET STATE #',i3,':  J =',i3,' IS =',i3/
     :        5x,'TRIPLET STATE #',i3,':  J =',i3,' N  =',i3,
     :        ' IS =',i3/)
         goto 4002
      end if
*
* set up 2 X 2 secular determinant for mixing
      isize = 2
      e(1,1) = elev1(istats) * econv
      e(2,2) = elevt(istatt) * econv
      e(1,2) = wso
      e(2,1) = wso
      lwork = 144
      call dsyev('V','L',isize,e,narray,eig,work,
     :     lwork,ierr)
      epert(1) = eig(1)
      epert(2) = eig(2)
      if (e(1,1) .ge. 0.d0) then
         cpert(1,1) = e(1,1)
         cpert(2,1) = e(2,1)
      else
         cpert(1,1) = -e(1,1)
         cpert(2,1) = -e(2,1)
      end if
      if (e(2,2) .ge. 0.d0) then
         cpert(1,2) = e(1,2)
         cpert(2,2) = e(2,2)
      else
         cpert(1,2) = -e(1,2)
         cpert(2,2) = -e(2,2)
      end if
      es = elev1(istats) * econv
      et = elevt(istatt) * econv
      write(6,1661) istats, es, istatt, et,
     :     wso, epert(1), cpert(1,1), cpert(2,1),
     :     epert(2), cpert(1,2), cpert(2,2)
 1661 format(/' SPECTROSCOPIC PERTURBATION - ZERO-ORDER LEVELS'/
     :     10x,'SINGLET:  STATE #',i3,'  ENERGY =',f10.3/
     :     10x,'TRIPLET:  STATE #',i3,'  ENERGY =',f10.3/
     :     10x,'SPIN-ORBIT MATRIX ELEMENT =',f10.5/
     :     5x,'MIXED LEVELS:',
     :     /10x,'E(1)=',f10.3,4x, 'COEFFS =',2f8.4
     :     /10x,'E(2)=',f10.3,4x, 'COEFFS =',2f8.4)
*
*     now compute cross sections for transitions out of perturbed levels
      write(6,1665) jfinl1, jfinl2
 1665 format(
     :     /' INTEGRAL CROSS SECTIONS CALCULATED FROM LOADED S-MATRICES'
     :     /5x,'JTOT-MAX =',i4,' (SINGLET)  ',i4,' (TRIPLET)'/)
*     clear sigma array
      allocate(sigma(nlevels + nlevelt, 4), stat=ialloc)
      if (ialloc .ne. 0) goto 4002
      sigma = 0d0
*
      maxjtot = max(jfinl1, jfinl2)
*
      jtot = 0
      jlpar = 1
      irdsng = 1
      irdtrp = 1
*
      allocate(jq(mmaxst), stat=ialloc)
      if (ialloc .ne. 0) goto 4003
      allocate(lq(mmaxst), stat=ialloc)
      if (ialloc .ne. 0) goto 4004
      allocate(inq(mmaxst), stat=ialloc)
      if (ialloc .ne. 0) goto 4005
      allocate(sreal(mmaxst * (mmaxst + 1) / 2), stat=ialloc)
      if (ialloc .ne. 0) goto 4006
      allocate(simag(mmaxst * (mmaxst + 1) / 2), stat=ialloc)
      if (ialloc .ne. 0) goto 4007
      allocate(srs(0:jfinl1, 2, mmax2s), stat=ialloc)
      if (ialloc .ne. 0) goto 4008
      allocate(sis(0:jfinl1, 2, mmax2s), stat=ialloc)
      if (ialloc .ne. 0) goto 4009
      allocate(srt(0:jfinl2, 2, mmax2t), stat=ialloc)
      if (ialloc .ne. 0) goto 4010
      allocate(sit(0:jfinl2, 2, mmax2t), stat=ialloc)
      if (ialloc .ne. 0) goto 4011
      allocate(js(0:jfinl1, 2, mmaxs), stat=ialloc)
      if (ialloc .ne. 0) goto 4012
      allocate(ins(0:jfinl1, 2, mmaxs), stat=ialloc)
      if (ialloc .ne. 0) goto 4013
      allocate(ls(0:jfinl1, 2, mmaxs), stat=ialloc)
      if (ialloc .ne. 0) goto 4014
      allocate(jt(0:jfinl2, 2, mmaxt), stat=ialloc)
      if (ialloc .ne. 0) goto 4015
      allocate(intt(0:jfinl2, 2, mmaxt), stat=ialloc)
      if (ialloc .ne. 0) goto 4016
      allocate(lt(0:jfinl2, 2, mmaxt), stat=ialloc)
      if (ialloc .ne. 0) goto 4017
      allocate(lngths(0:jfinl1, 2), stat=ialloc)
      if (ialloc .ne. 0) goto 4018
      allocate(lngtht(0:jfinl2, 2), stat=ialloc)
      if (ialloc .ne. 0) goto 4019
      allocate(exsmtp(0:maxjtot, 2), stat=ialloc)
      if (ialloc .ne. 0) goto 4020
      allocate(exsmtn(0:maxjtot, 2), stat=ialloc)
      if (ialloc .ne. 0) goto 4021
*
      lngtht = 0
*
*     read s-matrices.  Note that jtot need not be the same for the singlet
*     and triplet lebels
 100  if (irdsng.eq.1) then
*
*     read next s-matrix for singlet
*     this assumes that jlpar=1 is stored first, and also for triplet below
*
*     parameter to read lower triangle of s-matrix in sread
*     (overwritten with number of open channels)
         nopen = -1
         call sread (0, sreal, simag, jtot1, jlpar1,
     :        nu1, jq, lq, inq, ipack, jpack, lpack,
     :        1, mmax, nopen, lngth1, ierr)
         if (ierr .lt. -1) then
            write(6,102)
 102        format(/' ** READ ERROR IN STMIX. ABORT **'/)
            goto 1000
         end if
         sngsmt = .true.
         jlp = 1 - (jlpar1 - 1)/2
*     copy s-matrix for this jtot1/jlpar1
         lngths(jtot1,jlp) = lngth1
         len2 = lngth1*(lngth1+1)/2
         if (jlpar1.eq.1) then
            exsmtp(jtot1,1) = .true.
         else
            exsmtn(jtot1,1) = .true.
         end if
         do i = 1, lngth1
            js(jtot1,jlp,i) = jpack(i)
            ins(jtot1,jlp,i) = ipack(i)
            ls(jtot1,jlp,i) =lpack(i)
         end do
         do ii = 1, len2
            srs(jtot1,jlp,ii) = sreal(ii)
            sis(jtot1,jlp,ii) = simag(ii)
         end do
*
         if (jtot1.eq.jfinl1 .and. jlpar1.eq.1) then
            irdsng = 2
         end if
         if (jtot1.eq.jfinl1 .and. jlpar1.eq.-1) then
            irdsng = 0
         end if
      else
         if (irdsng.eq.2) sngsmt = .false.
      end if
*
*     read next s-matrix for triplet
      if (irdtrp.eq.1) then
*     parameter to read lower triangle of s-matrix in sread
*     (overwritten with number of open channels)
         nopen = -1
         call sread (0, sreal, simag, jtot2, jlpar2,
     :        nu2, jq, lq, inq, ipack, jpack, lpack,
     :        11, mmax, nopen, lngth2, ierr)
         if (ierr .lt. -1) then
            write(6,102)
            goto 1000
         end if
         trpsmt = .true.
         jlp = 1 - (jlpar2 - 1)/2
*     copy s-matrix for this jtot1/jlpar1
         lngtht(jtot2,jlp) = lngth2
         len2 = lngth2*(lngth2+1)/2
         jlp = 1 -(jlpar2 - 1)/2
         if (jlpar2.eq.1) then
            exsmtp(jtot2,2) = .true.
         else
            exsmtn(jtot2,2) = .true.
         end if
         do i = 1, lngth2
            jt(jtot2,jlp,i) = jpack(i)
            intt(jtot2,jlp,i) = ipack(i)
            lt(jtot2,jlp,i) = lpack(i)
         end do
         do ii = 1, len2
            srt(jtot2,jlp,ii) = sreal(ii)
            sit(jtot2,jlp,ii) = simag(ii)
         end do
         if (jtot2.eq.jfinl2 .and. jlpar2.eq.1) then
            irdtrp = 2
         end if
         if (jtot2.eq.jfinl2 .and. jlpar2.eq.-1) then
            irdtrp = 0
         end if
      else
         if (irstrp.eq.2) trpsmt = .false.
      end if
*
      if (sngsmt) then
         jtot = jtot1
         jlpar = jlpar1
      end if
      if (trpsmt) then
         jtot = jtot2
         jlpar = jlpar2
      end if
*
*     make sure that the first jtot for both the singlet and triplet is jtot=0
      if (jtot.eq.0 .and. jlpar.eq.1) then
         if (jtot1.ne.jtot .or. jlpar1.ne.jlpar1) goto 70
         if (jtot2.ne.jtot .or. jlpar2.ne.jlpar2) goto 70
      end if
      goto 80
 70   write(6,71) jtot1, jlpar1, jtot2, jlpar2
 71   format(' *** FIRST VALUES OF JTOT/JLPAR SHOULD BE',
     :     ' 0/1'/5x,'JTOT1=',i4,'  JLPAR1=',i3/
     :     5x,'JTOT2=',i4,'  JLPAR2=',i3/)
      goto 1000
*
*     loop back to the next jtot/jlpar, if either available for singlet or triplet
 80   if (irdsng.eq.2  .and. irdtrp.eq.2) then
         irdsng = 1
         irdtrp = 1
         if (jlpar.eq.-1) goto 1000
      end if
      sngsmt = .false.
      trpsmt = .false.
      if (jtot1.lt.jfinl1 .or. jlpar1.eq.1) goto 100
      if (jtot2.lt.jfinl2 .or. jlpar2.eq.1) goto 100
*
*     now compute squares of t-matrix elements for each jtot/jlpar
      do jtot = 0, maxjtot
         do jlp = 1, 2
            jlpar = 1 - (jlp - 1)*2
            if (exsmtp(jtot,1) .or. exsmtp(jtot,2)) then
               sngsmt = exsmtp(jtot,1)
               trpsmt = exsmtp(jtot,2)
            else
               sngsmt = exsmtn(jtot,1)
               trpsmt = exsmtn(jtot,2)
            end if
            lngth1 = lngths(jtot,jlp)
            call mixint(mmaxs, mmax2s, mmaxt, mmax2t, sigma,
     $           jtot, jlpar, istats, istatt, srs, sis, srt, sit,
     $           js, ins, ls, jt, intt, lt, lngth1, lngtht,
     $           jfinl1, jfinl2, nlevel1, nlevels, nlevel2, nlevelt,
     $           cpert, iptsng, ipttrp, sngsmt, trpsmt)
         end do
      end do
*
*     compute cross sections for sums of squares of t-matrix elements
*
*     transitions to/from unperturbed singlet levels
      fak = acos(-1.d0) * ang2 / (2.d0 * rmu1)
      do is = 1,nlevels
         isub = iptsng(is)
         xj = jlev1(isub)
         xjsng = jlev1(istats)
         densng = (2.d0 * xjsng + 1.d0)
     :        * (ered1 - elev1(istats))
         denrow = (2.d0 * xj + 1.d0)
     :        * (ered1 - elev1(isub))
*     transition into perturbed levels
         sigma(is,3) = sigma(is,1) * fak / denrow
         sigma(is,4) = sigma(is,2) * fak / denrow
*     transition from perturbed levels
         sigma(is,1) = sigma(is,1) * fak / densng
         sigma(is,2) = sigma(is,2) * fak / densng
      end do
*
*     transitions to/from unperturbed triplet levels
      isubt = ipttrp(istatt)
      do it = 1,nlevelt
         isub = ipttrp(it)
         xj = jlevt(it)
         xjtrp = jlevt(istatt)
         dentrp = (2.d0 * xjtrp + 1.d0)
     :        * (ered2 - elev2(isubt))
         denrow = (2.d0 * xj + 1.d0)
     :        * (ered2 - elev2(isub))
*     transition into perturbed levels
         sigma(it+nlevels,3) = sigma(it+nlevels,1)
     :        * fak / denrow
         sigma(it+nlevels,4) = sigma(it+nlevels,2)
     :        * fak / denrow
*     transition from perturbed levels
         sigma(it+nlevels,1) = sigma(it+nlevels,1)
     :        * fak / dentrp
         sigma(it+nlevels,2) = sigma(it+nlevels,2)
     :        * fak / dentrp
      end do
*
*     print cross sections for transitions to and from perturbed levels
      write (6,290)
 290  format(' CROSS SECTIONS TO AND FROM PERTURBED LEVELS'
     :     //4x,'#',2x,'S/T',2x,'N',3x,'J',3x,'IS',12x,
     :     'FROM',20x,'TO'/29x,'1',11x,'2',11x,'1',11x,'2')
      do is = 1,nlevels
         isub = iptsng(is)
         write(6,121) is,jlev1(isub),inlev1(isub),
     :        (sigma(is,j),j=1,4)
 121     format(1x,i4,3x,'S',4x,2i4,2x,4(1pe12.4))
      end do
      do it=1,nlevelt
         isub = ipttrp(it)
         write(6,122) it,nlevt(it),jlevt(it),inlevt(it),
     :        (sigma(it+nlevels,j),j=1,4)
 122     format(1x,i4,3x,'T',3i4,2x,4(1pe12.4))
      end do
*
      call mtime(cpu1,ela1)
      ela1 = ela1 - ela0
      cpu1 = cpu1 - cpu0
      call gettim(ela1,elaps)
      call gettim(cpu1,cpu)
      write(6,400) elaps, cpu
 400  format(/,' ** STMIX FINAL TIMING, ELAPSED: ',a,'  CPU: ',a,' **'/)
*
 1000 continue
 4022 deallocate(exsmtn)
 4021 deallocate(exsmtp)
 4020 deallocate(lngtht)
 4019 deallocate(lngths)
 4018 deallocate(lt)
 4017 deallocate(intt)
 4016 deallocate(jt)
 4015 deallocate(ls)
 4014 deallocate(ins)
 4013 deallocate(js)
 4012 deallocate(sit)
 4011 deallocate(srt)
 4010 deallocate(sis)
 4009 deallocate(srs)
 4008 deallocate(simag)
 4007 deallocate(sreal)
 4006 deallocate(inq)
 4005 deallocate(lq)
 4004 deallocate(jq)
 4003 deallocate(sigma)
 4002 deallocate(ipttrp)
 4001 deallocate(iptsng)
 4000 close(1)
      close(11)
      if (ialloc .ne. 0) write (6, 4100)
 4100 format (' *** INSUFFICIENT MEMORY OR SMT FILE CORRUPTED. ***')
      return
      end
* -------------------------------------------------------------------------
      subroutine mixint(mmaxs, mmax2s, mmaxt, mmax2t, sigma,
     $     jtot, jlpar, istats, istatt, srs, sis, srt, sit,
     $     js, ins, ls, jt, intt, lt, lngth1, lngtht,
     $     jfinl1, jfinl2, nlevel1, nlevels, nlevel2, nlevelt,
     $     cpert, iptsng, ipttrp, sngsmt, trpsmt)
*
* subroutine to compute the contribution of a given jtot/jlpar to the integral
* cross sections for collisional transitions out of a mixed singlet/triplet
* pair of rotational levels
*
*  author:  p. dagdigian
*  current revision date:  1-jun-2010
* -------------------------------------------------------------------------
      implicit double precision (a-h, o-z)
      logical sngsmt, trpsmt
      complex*8 t
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
      common /coconv/ econv, xmconv, ang2
      double precision sigma(nlevels + nlevelt, 4)
      dimension cpert(2,2), iptsng(1), ipttrp(1)
* storage for s-matrix elements:
*   second letter is real (r), imaginary (i)
*   third letter is singlet (s), triplet (t)
*   subscripts:  jtot, jlp (= 1/2 for jlpar = +1/-1), size of channel basis array
      dimension srs(0:jfinl1,2,mmax2s), sis(0:jfinl1,2,mmax2s),
     :   srt(0:jfinl2,2,mmax2t), sit(0:jfinl2,2,mmax2t)
* arrays with values of j, in, and l:
*   subscripts:  jtot, jlp (= 1/2 for jlpar = +1/-1), length of channel basis
      dimension js(0:jfinl1,2,mmaxs), ins(0:jfinl1,2,mmaxs),
     :  ls(0:jfinl1,2,mmaxs), jt(0:jfinl2,2,mmaxt),
     :  intt(0:jfinl2,2,mmaxt), lt(0:jfinl2,2,mmaxt)
* size of s-matrices for triplet levels
*   subscripts:  jtot, jlp (= 1/2 for jlpar = +1/-1)
      dimension lngtht(0:jfinl2,2)
* storage for singlet and triplet t-matrix elements
      double precision, dimension(:, :), allocatable ::
     $     tmatrs, tmatis, tmatrt, tmatit
*
* electron spin in the triplet state
      ispin = 1
*
      xjtot = jtot
      jlp = 1 - (jlpar - 1)/2
*
* compute t-matrix elements connecting unperturbed singlet levels
* with a perturbed level
      if (sngsmt) then
        isubs = iptsng(istats)
        do i=1,nlevels
          if (i.ne.isubs) then
            isub = iptsng(i)
            do irow=1,lngth1
              do icol=1,irow
                if ((js(jtot,jlp,irow).eq.jlev1(isubs) .and.
     :              ins(jtot,jlp,irow).eq.inlev1(isubs) .and.
     :              js(jtot,jlp,icol).eq.jlev1(isub) .and.
     :              ins(jtot,jlp,icol).eq.inlev1(isub)) .or.
     :              (js(jtot,jlp,icol).eq.jlev1(isubs) .and.
     :              ins(jtot,jlp,icol).eq.inlev1(isubs) .and.
     :              js(jtot,jlp,irow).eq.jlev1(isub) .and.
     :              ins(jtot,jlp,irow).eq.inlev1(isub))) then
* include contribution of this t-matrix element
* does not include diagonal term for singlet levels since perturbed singlet level
* is not in the list of singlet levels scanned through here
* hence, t-mat = - s-mat
                  ii = (irow*(irow-1))/2 + icol
                  t2 = srs(jtot,jlp,ii)**2
     :              + sis(jtot,jlp,ii)**2
* scale cross sections by fractional singlet character
                  do ipert=1,2
                    sigma(i,ipert) = sigma(i,ipert)
     :                + t2 * (2.d0 * xjtot + 1.d0)
     :                * cpert(1,ipert)**2
                  end do
                end if
              end do
            end do
          end if
        end do
      end if
*
* compute t-matrix elements connecting unperturbed triplet levels
* with a perturbed level
      if (trpsmt) then
        do i=1,nlevelt
          if (i.ne.istatt) then
            xn = nlevt(istatt)
            xnp = nlevt(i)
            xj = jlevt(istatt)
            xjp = jlevt(i)
            fjjp = sqrt((2.d0*xj+1.d0)*(2.d0*xjp+1.d0))
            jtmin = max(0,jtot-ispin)
            jtmax = min(jfinl2,jtot+ispin)
            xspin = ispin
* clear t-matrix array
            lmax2 = jtot + max(jlevt(istatt),jlevt(i)) + 2
            allocate(tmatrt(lmax2 + 1, lmax2 + 1), stat=ialloc)
            if (ialloc .ne. 0) goto 4000
            allocate(tmatit(lmax2 + 1, lmax2 + 1), stat=ialloc)
            if (ialloc .ne. 0) goto 4001
            tmatrt = 0d0
            tmatit = 0d0
*
            do jc=jtmin,jtmax
              xjc = jc
              fjc = 2.d0 * xjc + 1.d0
* total parity must be the same for all T-matrix elements included in sum over jc
              jlpt = jlp
              if (abs(jc-jtot).eq.1) then
                if (jlp.eq.1) then
                  jlpt = 2
                else
                  jlpt = 1
                end if
              end if
*
              lngth2 = lngtht(jc,jlpt)
* include check below if min jtot in triplet smt file > 0
              if (lngth2.gt.0) then
                do irow=1,lngth2
* flag to make sure perturbed level is the bra, and unperturbed level is theket
                  iflag = 0
                  if (jt(jc,jlpt,irow).eq.nlevt(istatt) .and.
     :              intt(jc,jlpt,irow).eq.inlevt(istatt)) then
                    iflag = 1
                  end if
                  do icol=1,irow
                    if ((jt(jc,jlpt,irow).eq.nlevt(istatt) .and.
     :                  intt(jc,jlpt,irow).eq.inlevt(istatt) .and.
     :                  jt(jc,jlpt,icol).eq.nlevt(i) .and.
     :                  intt(jc,jlpt,icol).eq.inlevt(i)) .or.
     :                  (jt(jc,jlpt,icol).eq.nlevt(istatt) .and.
     :                  intt(jc,jlpt,icol).eq.inlevt(istatt) .and.
     :                  jt(jc,jlpt,irow).eq.nlevt(i) .and.
     :                  intt(jc,jlpt,irow).eq.inlevt(i))) then
                      ii = (irow*(irow-1))/2 + icol
                      if (iflag.ne.1) then
                        xl = lt(jc,jlpt,icol)
                        xlp = lt(jc,jlpt,irow)
                      else
                        xl = lt(jc,jlpt,irow)
                        xlp = lt(jc,jlpt,icol)
                      end if
* convert s-matrix element to t-matrix element
                      t = -cmplx(srt(jc,jlpt,ii),
     :                  sit(jc,jlpt,ii))
* below for diagonal T-matrix element
                      if (irow.eq.icol) then
                        t = t + cmplx(1.d0, 0.d0)
                      end if
                      iph = xjp - xlp - xj + xl
                      phase = 1.d0
                      if (iph.ne.2*(iph/2)) phase = -1.d0
                      t = t * phase * fjjp * fjc
     :                  * xf6j(xspin,xn,xj,xl,xjtot,xjc)
     :                  * xf6j(xspin,xnp,xjp,xlp,xjtot,xjc)
                      l = xl
                      lp = xlp
                      if (l .gt. lmax2) print *, l, lmax2
                      if (lp .gt. lmax2) print *, lp, lmax2
                      tmatrt(l+1,lp+1) =
     :                  tmatrt(l+1,lp+1) + real(t)
                      tmatit(l+1,lp+1) =
     :                  tmatit(l+1,lp+1) + aimag(t)
* for initial level = final level but l.ne.lp, need to include both T(l,lp)
* and T(lp,l)
                      if (nlevt(i).eq.nlevt(istatt) .and.
     :                    inlevt(i).eq.inlevt(istatt) .and.
     :                    xl.ne.xlp) then
* note that T(lp,l) = T(l,lp)* (hermitean matrix)
                        tmatrt(lp+1,l+1) =
     :                    tmatrt(lp+1,l+1) + real(t)
                        tmatit(lp+1,l+1) =
     :                    tmatit(lp+1,l+1) - aimag(t)
                      end if
                    end if
                  end do
                end do
* if statement below for check on nonzero value of lngth2
              end if
            end do
            t2sum = 0.d0
            do i2=1,lmax2+1
              do j2=1,lmax2+1
                t2sum = t2sum + tmatrt(i2,j2)**2
     :               + tmatit(i2,j2)**2
              end do
            end do
* scale cross sections by fractional triplet character
            do ipert=1,2
              sigma(i+nlevels,ipert) =
     :          sigma(i+nlevels,ipert)
     :          + t2sum * (2.d0 * xjtot + 1.d0)
     :          * cpert(2,ipert)**2
            end do
            deallocate(tmatrt)
            deallocate(tmatit)
          end if
        end do
      end if
*
* finally, compute t-matrix elements between the two perturbing levels
*
* determine maximum value of orbital angular momentum
* first, find max j in singlet and n in triplet perturbed level
      if (sngsmt) then
        isubs = iptsng(istats)
        jsgmax = jlev1(isubs)
      end if
      if (trpsmt) then
        jtpmax = jlevt(istatt)
      end if
      if (sngsmt.and.trpsmt) then
        jmolmx = max(jsgmax,ntpmax)
      else
        if (sngsmt) jmolmx = jsgmax
        if (trpsmt) jmolmx = jtpmax
      end if
* now compute max l
      lmax = jmolmx + jtot
*
* next, clear arrays with real and imag parts of the t-matrices
* for the perturbed singlet and triplet level
      allocate(tmatrt(lmax + 1, lmax + 1), stat=ialloc)
      if (ialloc .ne. 0) goto 4000
      allocate(tmatit(lmax + 1, lmax + 1), stat=ialloc)
      if (ialloc .ne. 0) goto 4001
      allocate(tmatrs(lmax + 1, lmax + 1), stat=ialloc)
      if (ialloc .ne. 0) goto 4002
      allocate(tmatis(lmax + 1, lmax + 1), stat=ialloc)
      if (ialloc .ne. 0) goto 4003
      tmatrs = 0.d0
      tmatis = 0.d0
      tmatrt = 0.d0
      tmatit = 0.d0
*
* compute t-matrix for perturbed singlet level
      if (sngsmt) then
        isubs = iptsng(istats)
        do irow=1,lngth1
          do icol=1,irow
            if (js(jtot,jlp,irow).eq.jlev1(isubs) .and.
     :          ins(jtot,jlp,irow).eq.inlev1(isubs) .and.
     :          js(jtot,jlp,icol).eq.jlev1(isubs) .and.
     :          ins(jtot,jlp,icol).eq.inlev1(isubs)) then
* include contribution of this t-matrix element
              ii = (irow*(irow-1))/2 + icol
              lsub = ls(jtot,jlp,irow) + 1
              lpsub = ls(jtot,jlp,icol) + 1
              tmatrs(lsub,lpsub) = tmatrs(lsub,lpsub)
     :            - srs(jtot,jlp,ii)
* conversion from s-matrix element to t-matrix element for
* initial channel = final channel
              if (irow.eq.icol) tmatrs(lsub,lpsub)
     :             = tmatrs(lsub,lpsub) + 1.d0
              tmatis(lsub,lpsub) = tmatis(lsub,lpsub)
     :            + sis(jtot,jlp,ii)
* need to include both T(l,lp) and T(lp,l)
* note that T(lp,l) = T(l,lp)* (hermitean matrix)
              if (lsub.ne.lpsub) then
                tmatrs(lpsub,lsub) = tmatrs(lsub,lpsub)
                tmatis(lpsub,lsub) = -tmatis(lsub,lpsub)
              end if
            end if
          end do
        end do
      end if
*
* compute t-matrix for perturbed triplet level
      if (trpsmt) then
        xn = nlevt(istatt)
        xj = jlevt(istatt)
        fjjp = 2.d0 * xj + 1.d0
        jtmin = max(0,jtot-ispin)
        jtmax = min(jfinl2,jtot+ispin)
        xspin = ispin
        jlp = 1 - (jlpar - 1)/2
        do jc=jtmin,jtmax
          xjc = jc
          fjc = 2.d0 * xjc + 1.d0
* total parity must be the same for all T-matrix elements included in sum over jc
          jlpt = jlp
          if (abs(jc-jtot).eq.1) then
            if (jlp.eq.1) then
              jlpt = 2
            else
              jlpt = 1
            end if
          end if
*
          lngth2 = lngtht(jc,jlpt)
* include check below if min jtot in triplet smt file > 0
          if (lngth2.gt.0) then
            do irow=1,lngth2
              do icol=1,irow
                if (jt(jc,jlpt,irow).eq.nlevt(istatt) .and.
     :              intt(jc,jlpt,irow).eq.inlevt(istatt) .and.
     :              jt(jc,jlpt,icol).eq.nlevt(istatt) .and.
     :              intt(jc,jlpt,icol).eq.inlevt(istatt)) then
                  ii = (irow*(irow-1))/2 + icol
                  lsub = lt(jc,jlpt,icol) + 1
                  lpsub = lt(jc,jlpt,irow) + 1
                  xl = lsub - 1
                  xlp = lpsub - 1
* convert s-matrix element to t-matrix element
                  t = -cmplx(srt(jc,jlpt,ii),
     :                sit(jc,jlpt,ii))
                  if (irow.eq.icol) then
                    t = t + cmplx(1.d0, 0.d0)
                  end if
                  iph = xjp - xlp - xj + xl
                  phase = 1.d0
                  if (iph.ne.2*(iph/2)) phase = -1.d0
                  t = t * phase * fjjp * fjc
     :                * xf6j(xspin,xn,xj,xl,xjtot,xjc)
     :                * xf6j(xspin,xn,xj,xlp,xjtot,xjc)
                  tmatrt(lsub,lpsub) = tmatrt(lsub,lpsub)
     :                + real(t)
                  tmatit(lsub,lpsub) = tmatit(lsub,lpsub)
     :                + aimag(t)
* need to include both T(l,lp) and T(lp,l)
* note that T(lp,l) = T(l,lp)* (hermitean matrix)
                  if (lsub.ne.lpsub) then
                    tmatrt(lpsub,lsub) = tmatrt(lpsub,lsub)
     :                + real(t)
                    tmatit(lpsub,lsub) = tmatit(lpsub,lsub)
     :                - aimag(t)
                  end if
                end if
              end do
            end do
* if statement below for check on nonzero value of lngth2
          end if
        end do
      end if
*
* finally, compute sum of squares of t-matrix elements
* involving the two perturbed levels
      if (cpert(1,1)**2 .gt. 0.5d0) then
* singlet lower in energy than triplet
        do i=1,lmax+1
          do j=1,lmax+1
            t = cpert(1,1)**2
     :          * cmplx(tmatrs(i,j),tmatis(i,j))
     :          + cpert(2,1)**2
     :          * cmplx(tmatrt(i,j), tmatit(i,j))
            t2 = real(t * conjg(t))
            sigma(istats,1) = sigma(istats,1)
     :          + (2.d0 * xjtot + 1.d0) * t2
*
            t = cpert(1,1) * cpert(1,2)
     :          * cmplx(tmatrs(i,j),tmatis(i,j))
     :          + cpert(2,1) * cpert(2,2)
     :          * cmplx(tmatrt(i,j),tmatit(i,j))
            t2 = real(t * conjg(t))
            sigma(istats,2) = sigma(istats,2)
     :          + (2.d0 * xjtot + 1.d0) * t2
*
            sigma(istatt+nlevels,1) = sigma(istats,2)
*
            t = cpert(1,2)**2
     :          * cmplx(tmatrs(i,j),tmatis(i,j))
     :          + cpert(2,2)**2
     :          * cmplx(tmatrt(i,j),tmatit(i,j))
            t2 = real(t * conjg(t))
            sigma(istatt+nlevels,2) =
     :          sigma(istatt+nlevels,2)
     :          + (2.d0 * xjtot + 1.d0) * t2
          end do
        end do
      else
* triplet lower in energy than singlet
        do i=1,lmax+1
          do j=1,lmax+1
            t = cpert(1,2)**2
     :          * cmplx(tmatrs(i,j),tmatis(i,j))
     :          + cpert(2,2)**2
     :          * cmplx(tmatrt(i,j),tmatit(i,j))
            t2 = real(t * conjg(t))
            sigma(istats,2) = sigma(istats,2)
     :          + (2.d0 * xjtot + 1.d0) * t2
*
            t = cpert(1,1) * cpert(1,2)
     :          * cmplx(tmatrs(i,j),tmatis(i,j))
     :          + cpert(2,1) * cpert(2,2)
     :          * cmplx(tmatrt(i,j),tmatit(i,j))
            t2 = real(t * conjg(t))
            sigma(istats,1) = sigma(istats,1)
     :          + (2.d0 * xjtot + 1.d0) * t2
*
            sigma(istatt+nlevels,2) = sigma(istats,1)
*
            t = cpert(1,1)**2
     :          * cmplx(tmatrs(i,j),tmatis(i,j))
     :          + cpert(2,1)**2
     :          * cmplx(tmatrt(i,j),tmatit(i,j))
            t2 = real(t * conjg(t))
            sigma(istatt+nlevels,1) =
     :          sigma(istatt+nlevels,1)
     :          + (2.d0 * xjtot + 1.d0) * t2
          end do
        end do
      end if
 4004 deallocate(tmatis)
 4003 deallocate(tmatrs)
 4002 deallocate(tmatit)
 4001 deallocate(tmatrt)
 4000 continue
      if (ialloc .ne. 0) write (6, 4100)
 4100 format ('  *** INSUFFICIENT MEMORY. ***')
      return
      end
* -------------------------------------------------------------------------
