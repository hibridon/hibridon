*     -------------------------------------------------------------
      subroutine hypxsc(flname, a)
*
*     subroutine to compute hyperfine-resolved integral
*     cross sections
*
*     reference:  alexander and dagdigian, jcp 83, 2191 (1985)
*     see also corey and mccourt, jpca 87, 2723 (1983) for
*     expression for spin-resolved T-matrix elements
*
*     this subroutine requires close-coupled s-matrix for
*     both parities
*
*     cross sections written to both terminal output and
*     {jobname}n.hfx file
*
*     author:  p.j. dagdigian
*
*     this subroutine is a complete rewrite of an earlier
*     subroutine written by j. klos and f. lique
*
*     this module has been moved out of hibrid1.f (5-jun-2013 by q. ma)
*
*     modified to treat molecule-molecule collisions (p. dagdigian)
*     used formalism presented by offer et al., jcp 100, 362 (1994)
*
*     modified to deal with cases for which jtot=0 has only one parity
*
*     current revision:  8-jan-2018 by p. dagdigian
*     -------------------------------------------------------------
      implicit double precision (a-h,o-z)
      complex(8) t, tf
      character*(*) flname
      character*20 cdate
      character*40 smtfil, hfxfil
      character*10 elaps, cpu
      logical csflg, flaghf, flgsu, twmol, nucrs,
     :     batch, fast, lpar2, lpar, exstfl
      include "common/parpot"
      common /coisc1/ jlev(1)
      common /coisc3/ inlev(1)
      common /coisc5/ jout(1)
      common /cosc1/ elev(1)
      common /coisc7/ jlevh(1)
      common /coisc8/ iflevh(1)
      common /coisc9/ inlevh(1)
      common /cosc3/ elevh(1)
      common /coisc10/ ipack(1)
      common /coisc11/ jpack(1)
      common /coisc12/ lpack(1)
      common /cojq/ jq(1)
      common /colq/ lq(1)
      common /coinq/ inq(1)
      common /coj12/ j12q(1)
      common /coconv/ econv, xmconv, ang2
      common /codim/ nairy, mmax
      common /coselb/ ibasty
      dimension a(4)
*
*     storage for S-matrix elements red from .smt file
      real(8), dimension(:), allocatable :: sreal, simag
*     storage for cross sections
      real(8), dimension(:, :), allocatable :: sigma
*     storage for T-matrix elements
      real(8), dimension(:, :), allocatable :: tmatr, tmati
*     storage for s-matrix elements:
*     second letter is real (r), imaginary (i)
*     subscripts: jtot, jlp (= 1/2 for jlpar = +1/-1), size of channel
*     basis array
      real(8), dimension(:, :, :), allocatable :: sr, si
*     arrays with values of j, in, and l:
*     subscripts: jtot, jlp (= 1/2 for jlpar = +1/-1), length of channel
*     basis
      integer, dimension(:, :, :), allocatable :: j, in, l, j12
*     length of arrays
*     subscripts: jtot, jlp (= 1/2 for jlpar = +1/-1)
      integer, dimension(:, :), allocatable :: length
*     for positive (p) and (n) negative parity
      logical, dimension(:), allocatable :: exsmtp, exsmtn
*
*     initialize timer and arrays
      call mtime(cpu0, ela0)
*
*     nucspin is nuclear spin I times 2
*     j1min and j2max are the min and max values of the rotational
*     angular momentum for which hyperfine-resolved cross sections
*     computed
      nucspin = a(2)
      finuc = nucspin/2.d0
      j1min = a(3)
      j2max = a(4)
*
*     generate filename of smt file and check if it is present
*
      iener = a(1)
      call gennam(smtfil,flname,iener,'smt',lends)
      inquire(file = smtfil(1:lends), exist = exstfl)
      if (.not. exstfl) then
         write(6,10) smtfil(1:lends)
 10      format(/' ** FILE ',(a),' NOT FOUND **'/)
         return
      end if
*
*     open s-matrix file
*
      call openf(1, smtfil, 'tu', 0)
*
*     initialize timer and arrays
      call mtime(cpu0, ela0)
*
*     read header of s-matrix file
*
      call sinqr(1, mjtot, mchmx)
      mchmx2 = mchmx * (mchmx + 1) / 2
      call rdhead(1,cdate,ered,rmu,csflg,flaghf,flgsu,
     :     twmol,nucrs,jfrst,jfinl,jtotd,numin,numax,nud,
     :     nlevel,nlvop,nnout,jlev,inlev,elev,jout)
*
*     we need the s-matrices as lower triangles, so nnout  m u s t  be > 0
*
      if (nnout.lt.0) then
         write(6,11)
 11      format(/' ** NNOUT < 0, ABORT **'/)
         close(1)
         return
      end if
      nout = nnout
*
      if (flaghf) then
         fspin = 0.5d0
      else
         fspin = 0.d0
      end if
*
*     unsupported basis types
      if (ibasty.eq.10 .or. ibasty.eq.12 .or. ibasty.eq.13 .or.
     :  ibasty.eq.15 .or. ibasty.eq.22 .or. ibasty.eq.23) then
         write(6,14) ibasty
 14      format(/' *** HYPERFINE CROSS SECTIONS FOR BASIS TYPE =',
     :        i3,' NOT IMPLEMENTED ***'/)
         close(1)
         return
      end if
*
*     molecule-surface collisions not implemented
      if (flgsu) then
         write(6,15)
 15      format(/' *** HYPERFINE CROSS SECTIONS FOR SURFACE -',
     :        ' COLLISIONS NOT IMPLEMENTED ***'/)
         close(1)
         return
      end if
*
*     CS cross sections not implemented
      if (csflg) then
         write(6,16)
 16      format(/' *** CS HYPERFINE CROSS SECTIONS CROSS ',
     :        'SECTIONS NOT IMPLEMENTED ***'/)
         close(1)
         return
      end if
*
*     delta-jtot should be equal to one
      if (iabs(jtotd).ne.1)
     :     then
         write(6,116)
 116     format(/' *** DELTA-JTOT MUST BE EQUAL TO ONE',
     :        ' *** '/)
         close(1)
         return
      end if
*
*     generate output file name
      call gennam(hfxfil,flname,iener,'hfx',lendh)
      call openf(11, hfxfil(1:lendh),'sf',0)
*
      ee = ered * econv
      write (6, 22) label, potnam, smtfil, cdate
      write (11, 22) label, potnam, smtfil, cdate
 22   format(/'% HYPERFINE-RESOLVED CROSS SECTIONS'/
     :     '%      LABEL:     ',(a)/
     :     '%      POT NAME:  ',(a)/
     :     '%      S-MATRICES READ FROM FILE ',(a),/,
     :     '%      WRITTEN:   ',(a))
      write (6,24) jfinl, ee, finuc
      write (11,24) jfinl, ee, finuc
 24   format('%      JTOT2:     ',i10/
     :     '%      ETOT:      ',f10.3,' CM(-1)'
     :     //'% NUCLEAR SPIN = ',f4.1/)
*
*     set up list of hyperfine levels for which cross sections
*     to be calculated
      nlevelh = 0
      do 60 i=1, nlevel
         if (twmol) then
           ij1 = jlev(i)/10
         else
           ij1 = jlev(i)
         endif
         if (ij1.lt.j1min .or. ij1.gt.j2max)
     :        goto 60
*     check that level energetically allowed
         if (elev(i).gt.ered) goto 60
         fj = ij1 + fspin
         ffmin = abs(fj - finuc)
         ffmax = fj + finuc
         nhyp = ffmax - ffmin + 1
         do if=1,nhyp
            nlevelh = nlevelh + 1
            jlevh(nlevelh) = jlev(i)
            inlevh(nlevelh) = inlev(i)
*     iflevh is an integer which equals F if F is integral,
*     equals F - 0.5 if F is half-integral
            ff = ffmin + (if - 1)
            iflevh(nlevelh) = ff
            elevh(nlevelh) = elev(i)
         end do
 60   continue
*
*     clear sigma array
      allocate(sigma(nlevelh, nlevelh), stat=ialloc)
      if (ialloc .ne. 0) goto 4000
      sigma = 0d0
*
*     clear length array, in case minimum jtot > 0
      allocate(length(0:jfinl, 2), stat=ialloc)
      if (ialloc .ne. 0) goto 4001
      length = 0d0
*
      jtot = 0
      jlpar = 1
*
*     read s-matrix elements
*     this assumes that jlpar=1 is stored first
*
      allocate(sreal(mchmx2), stat=ialloc)
      if (ialloc .ne. 0) goto 4002
      allocate(simag(mchmx2), stat=ialloc)
      if (ialloc .ne. 0) goto 4003
      allocate(sr(0:jfinl, 2, mchmx2), stat=ialloc)
      if (ialloc .ne. 0) goto 4004
      allocate(si(0:jfinl, 2, mchmx2), stat=ialloc)
      if (ialloc .ne. 0) goto 4005
      allocate(j(0:jfinl, 2, mchmx), stat=ialloc)
      if (ialloc .ne. 0) goto 4006
      allocate(in(0:jfinl, 2, mchmx), stat=ialloc)
      if (ialloc .ne. 0) goto 4007
      allocate(l(0:jfinl, 2, mchmx), stat=ialloc)
      if (ialloc .ne. 0) goto 4008
      allocate(j12(0:jfinl, 2, mchmx), stat=ialloc)
      if (ialloc .ne. 0) goto 4009
      allocate(exsmtp(0:jfinl), stat=ialloc)
      if (ialloc .ne. 0) goto 4015
      allocate(exsmtn(0:jfinl), stat=ialloc)
      if (ialloc .ne. 0) goto 4010
*     fix for case where jtot=0 has only one parity
      jfrst = 2000
      jfin = 0
      do ij = 1, jfinl+1
        do ip = 1, 2
          length(ij,ip) = 0
        end do
      end do
*     parameter to read lower triangle of open channel(s)
 100  nopen = -1
      call sread (0, sreal, simag, jtot, jlpar,
     :     nu, jq, lq, inq, ipack, jpack, lpack,
     :     1, mmax, nopen, lngth, ierr)
      if (ierr .lt. -1) then
         write(6,102)
 102     format(/' ** READ ERROR IN HYPXSC. ABORT **'/)
         goto 1000
      end if
      jfrst = min(jfrst, jtot)
      jfin = max(jfin, jtot)
      jlp = 1 - (jlpar - 1)/2
*     copy s-matrix for this jtot1/jlpar1
      length(jtot,jlp) = lngth
      len2 = lngth*(lngth + 1)/2
      if (jlpar.eq.1) then
         exsmtp(jtot) = .true.
      else
         exsmtn(jtot) = .true.
      end if
      do i = 1, lngth
         j(jtot,jlp,i) = jpack(i)
         in(jtot,jlp,i) = ipack(i)
         l(jtot,jlp,i) = lpack(i)
         j12(jtot,jlp,i) = j12q(i)
      end do
      do ii = 1, len2
         sr(jtot, jlp, ii) = sreal(ii)
         si(jtot, jlp, ii) = simag(ii)
      end do
*
*     loop back to the next jtot/jlpar
      if (jtot.lt.jfinl .or. jlpar.eq.1) goto 100
      jfinl =jfin
*
*     now compute squares of T-matrix elements for each tot (vector sum
*     of jtot + nucspin)/parity pair
*
      if (.not. twmol) then
*
*     this option for atom - molecule collisions
*
      fhspin = 0.d0
      if (flaghf .and. nucspin.eq.2*(nucspin/2) .or.
     :     .not.flaghf .and. nucspin.ne.2*(nucspin/2))
     :     fhspin = 0.5d0
      iftmn = jfrst + fspin - finuc - fhspin
      iftmn = max(0,iftmn)
      iftmx = jfinl + fspin + finuc - fhspin

      jmx = 0
      do i = i, nlevel
        jmx = max(jmx, jlev(i))
      end do
      lmax = iftmx + jmx + 1
      jrmx = lmax + jmx + 1
      idimr = jrmx + 1
      idim = idimr * (lmax + 1)
      allocate(tmatr(idim, idim), stat=ialloc)
      if (ialloc .ne. 0) goto 4011
      allocate(tmati(idim, idim), stat=ialloc)
      if (ialloc .ne. 0) goto 4011
      do iftot = iftmn, iftmx
        xftot = iftot + fhspin
        do jlp = 1, 2
          jlpar = 1 - (jlp -1)*2
          write(6,334) xftot,jlpar
334       format(' Computing partial wave J_tot =',f5.1,
     :        ', jlpar=',i4)
          do i=1,nlevelh
            do ii=i,nlevelh
              xj = jlevh(i) + fspin
              xjp = jlevh(ii) + fspin
              xf = iflevh(i) + fhspin
              xfp = iflevh(ii) + fhspin
              fffp = sqrt((2.d0*xf+1.d0)*(2.d0*xfp+1.d0))
              xjttmn = max(fhspin, xftot - finuc)
              jttmin = xjttmn - fspin
              jttmin = max(jfrst, jttmin)
              jttmax = xftot + finuc
              jttmax = min(jfinl, jttmax)
*     clear T-matrix array
              tmatr = 0.d0
              tmati = 0.d0
*     sum over jtot consistent with vector addition
*     ftot = jtot + nucspin
              do jtot=jttmin,jttmax
                xjtot = jtot + fspin
                fjtot = 2.d0 * xjtot + 1.d0
*     total parity must be the same for all T-matrix elements
*     in sum over jtot (for the same parity, jlparf
*     changes sign for each increase in jtot by 1)
                jlparf = jlpar*(-1)**(iftot - jtot)
                jlpt = 1 + (1 - jlparf)/2
                if (length(jtot,jlpt) .gt. 0) then
                  do irow=1,length(jtot,jlpt)
*     flag to make sure initial level is the bra, final level the ket
                    iflag = 0
                    if (j(jtot,jlpt,irow).eq.jlevh(i) .and.
     :                  in(jtot,jlpt,irow).eq.inlevh(i)) then
                      iflag = 1
                    end if
                    do icol=1,irow
                      if (j(jtot,jlpt,irow).eq.jlevh(i) .and.
     :                    in(jtot,jlpt,irow).eq.inlevh(i) .and.
     :                    j(jtot,jlpt,icol).eq.jlevh(ii) .and.
     :                    in(jtot,jlpt,icol).eq.inlevh(ii) .or.
     :                    j(jtot,jlpt,icol).eq.jlevh(i) .and.
     :                    in(jtot,jlpt,icol).eq.inlevh(i) .and.
     :                    j(jtot,jlpt,irow).eq.jlevh(ii) .and.
     :                    in(jtot,jlpt,irow).eq.inlevh(ii)) then
                        is = (irow*(irow - 1))/2 + icol
                        if (iflag.ne.1) then
                          xl = l(jtot,jlpt,icol)
                          xlp = l(jtot,jlpt,irow)
                        else
                          xl = l(jtot,jlpt,irow)
                          xlp = l(jtot,jlpt,icol)
                        end if
*     convert S-matrix element to T-matrix element
                        t = -cmplx(sr(jtot,jlpt,is),
     :                      si(jtot,jlpt,is))
*     next statement for diagonal T-matrix element
                        if (irow.eq.icol) then
                          t = t + cmplx(1.d0, 0.d0)
                        end if
                        iph = xfp - xlp - xf + xl
                        phase = 1.d0
                        if (iph.ne.2*(iph/2)) phase = -1.d0
                        t = t * phase * fffp * fjtot
     :                      * xf6j(finuc,xj,xf,xl,xftot,xjtot)
     :                    * xf6j(finuc,xjp,xfp,xlp,xftot,xjtot)
                        ll = xl
                        lp = xlp
                        tmatr(ll+1,lp+1) =
     :                      tmatr(ll+1,lp+1) + real(t)
                        tmati(ll+1,lp+1) =
     :                      tmati(ll+1,lp+1) + aimag(t)
*     for initial level = final level, but l.ne.lp, need to include
*     both T(l,lp) and T(lp,l)
                        if (jlevh(i).eq.jlevh(ii) .and.
     :                      inlevh(i).eq.inlevh(ii) .and.
     :                      xl.ne.xlp) then
*     note that T(l,lp) = T(lp,l)* (Hermitean matrix)
                          tmatr(lp+1,ll+1) =
     :                        tmatr(lp+1,ll+1) + real(t)
                          tmati(lp+1,ll+1) =
     :                        tmati(lp+1,ll+1) - aimag(t)
                        end if
*     if statement below is end of tests for triangle relations
                      end if
                    end do
                  end do
*     end of if statement checking in length
                end if
              end do
              t2sum = 0.d0
              do i1=1,idim
                do i2=1,idim
                  t2sum = t2sum + tmatr(i1,i2)**2
     :                + tmati(i1,i2)**2
                end do
              end do
              sigma(i,ii) = sigma(i,ii) + t2sum
     :            * (2.d0 * xftot + 1.d0)
              if (i.ne.ii) then
                sigma(ii,i) = sigma(ii,i) + t2sum
     :              * (2.d0 * xftot + 1.d0)
              end if
            end do
          end do
        end do
      end do
      deallocate(tmatr)
      deallocate(tmati)
*
*  end of atom-molecule section
*
*     this option for molecule - molecule collisions
*
      else
*
      fhspin = 0.d0
      if (flaghf .and. nucspin.eq.2*(nucspin/2) .or.
     :     .not.flaghf .and. nucspin.ne.2*(nucspin/2))
     :     fhspin = 0.5d0
      iftmn = jfrst + fspin - finuc - fhspin
      iftmn = max(0,iftmn)
      iftmx = jfinl + fspin + finuc - fhspin
*     determine maximum values of j1 and j2
      j1mx = 0
      j2mx = 0
      do i = 1, nlevel
        j1 = jlev(i)/10
        j2 = mod(jlev(i),10)
        j1mx = max(j1mx,j1)
        j2mx = max(j2mx,j2)
      end do
*     determine limits on jR and dimensions of tmat arrays
      lmax = iftmx + j1mx + j2mx + 1
      jrmx = lmax + j2mx + 1
      idimr = jrmx + 1
      idim = idimr * (lmax + 1)
*     next allocate scratch arrays
      allocate(tmatr(idim, idim), stat=ialloc)
      if (ialloc .ne. 0) goto 4011
      allocate(tmati(idim, idim), stat=ialloc)
      if (ialloc .ne. 0) goto 4011
      do iftot = iftmn, iftmx
        xftot = iftot + fhspin
        do jlp = 1, 2
          jlpar = 1 - (jlp -1) * 2
          write(6,334) xftot, jlpar
          do i=1,nlevelh
            do ii=i,nlevelh
              xj = (jlevh(i)/10) + fspin
              xjp = (jlevh(ii)/10) + fspin
              xj2 = mod(jlevh(i),10)
              xj2p = mod(jlevh(ii),10)
              xf = iflevh(i) + fhspin
              xfp = iflevh(ii) + fhspin
              fffp = sqrt((2.d0*xf+1.d0)*(2.d0*xfp+1.d0))
              xjttmn = max(fhspin, xftot - finuc)
              jttmin = xjttmn - fspin
              jttmin = max(jfrst, jttmin)
              jttmax = xftot + finuc
              jttmax = min(jfinl, jttmax)
*     clear tmat arrays
              tmatr = 0.d0
              tmati = 0.d0
*     sum over jtot consistent with vector addition
*     ftot = jtot + nucspin
              do jtot=jttmin,jttmax
                xjtot = jtot + fspin
                fjtot = 2.d0 * xjtot + 1.d0
*     total parity must be the same for all T-matrix elements
*     in sum over jtot (for the same parity, jlparf
*     changes sign for each increase in jtot by 1)
                jlparf = jlpar*(-1)**(iftot - jtot)
                jlpt = 1 + (1 - jlparf)/2
                do irow=1,length(jtot,jlpt)
*     flag to make sure initial level is the bra, final level the ket
                  iflag = 0
                  if (j(jtot,jlpt,irow).eq.jlevh(i) .and.
     :                in(jtot,jlpt,irow).eq.inlevh(i)) then
                    iflag = 1
                  end if
                  do icol=1,irow
                    if (j(jtot,jlpt,irow).eq.jlevh(i) .and.
     :                  in(jtot,jlpt,irow).eq.inlevh(i) .and.
     :                  j(jtot,jlpt,icol).eq.jlevh(ii) .and.
     :                  in(jtot,jlpt,icol).eq.inlevh(ii) .or.
     :                  j(jtot,jlpt,icol).eq.jlevh(i) .and.
     :                  in(jtot,jlpt,icol).eq.inlevh(i) .and.
     :                  j(jtot,jlpt,irow).eq.jlevh(ii) .and.
     :                  in(jtot,jlpt,irow).eq.inlevh(ii)) then
                      is = (irow*(irow - 1))/2 + icol
                      if (iflag.ne.1) then
                        xl = l(jtot,jlpt,icol)
                        xlp = l(jtot,jlpt,irow)
                        xj12 = j12(jtot,jlpt,icol) + fspin
                        xj12p = j12(jtot,jlpt,irow) + fspin
                      else
                        xl = l(jtot,jlpt,irow)
                        xlp = l(jtot,jlpt,icol)
                        xj12 = j12(jtot,jlpt,irow) + fspin
                        xj12p = j12(jtot,jlpt,icol) + fspin
                      end if
                      fj12 = sqrt(2.d0 * xj12 + 1.d0)
                      fj12p = sqrt(2.d0 * xj12p + 1.d0)
*     convert S-matrix element to T-matrix element
                      t = -cmplx(sr(jtot,jlpt,is),
     :                    si(jtot,jlpt,is))
*     next statement for diagonal T-matrix element
                      if (irow.eq.icol) then
                        t = t + cmplx(1.d0, 0.d0)
                      end if
*     sums over jR and jRp
                      jrmin = abs(xj2 - xl)
                      jrmax = xj2 + xl
                      jrpmin = abs(xj2p - xlp)
                      jrpmax = xj2p + xlp
                      do jr = jrmin, jrmax
                        xjr = jr
                        fjr = sqrt(2.d0 * xjr + 1.d0)
                        do jrp = jrpmin, jrpmax
                          xjrp = jrp
                          fjrp = sqrt(2.d0 * xjrp + 1.d0)
                          iph = xjr + xjrp + xl + xlp
     :                        + xj2 + xj2p
                          ph = 1.d0
                          if (iph.ne.2*(iph/2)) ph = -1.d0
                          tf = t * ph * fjtot * fffp
     :                        * fjr * fjrp * fj12 * fj12p
     :                        * xf6j(xj,xj2,xj12,xl,xjtot,xjr)
     :                        * xf6j(xjp,xj2p,xj12p,xlp,xjtot,xjrp)
     :                        * xf6j(xjr,xj,xjtot,finuc,xftot,xf)
     :                        * xf6j(xjrp,xjp,xjtot,finuc,xftot,xfp)
                          ll = xl
                          lp = xlp
                          is = (ll + 1) * (idimr - 1) + (jr + 1)
                          isp = (lp + 1) * (idimr - 1) + (jrp + 1)
                          tmatr(is,isp) =
     :                        tmatr(is,isp) + real(tf)
                          tmati(is,isp) =
     :                        tmati(is,isp) + aimag(tf)
*     for initial level = final level, but l.ne.lp or j12.ne.j12p, need to include
*     both T(l,lp) and T(lp,l)
                          if (jlevh(i).eq.jlevh(ii) .and.
     :                        inlevh(i).eq.inlevh(ii) .and.
     :                        irow.ne.icol) then
*     note that T(l,lp) = T(lp,l)* (Hermitean matrix)
                            tmatr(isp,is) =
     :                          tmatr(isp,is) + real(tf)
                            tmati(isp,is) =
     :                          tmati(isp,is) - aimag(tf)
                          end if
                        end do
                      end do
                    end if
                  end do
                end do
              end do
              t2sum = 0.d0
              do i1 = 1, idim
                do i2 = 1, idim
                  t2sum = t2sum + tmatr(i1,i2)**2
     :                + tmati(i1,i2)**2
                end do
              end do
              sigma(i,ii) = sigma(i,ii) + t2sum
     :               * (2.d0 * xftot + 1.d0)
              if (i.ne.ii) then
                 sigma(ii,i) = sigma(ii,i) + t2sum
     :               * (2.d0 * xftot + 1.d0)
              end if
            end do
          end do
        end do
      end do
      deallocate(tmatr)
      deallocate(tmati)
*
*  end of molecule-molecule section
*
      end if
*
*     compute cross sections for sums of squares of T-matrix elements
      fak = acos(-1.d0) * ang2 / (2.0d0 * rmu)
      do i=1,nlevelh
         if (twmol) then
           ij2 = mod(jlevh(i),10)
         else
           ij2 = 0.d0
         endif
         ffi = iflevh(i) + fhspin
         denrow = (2.d0 * ffi + 1.d0) * (2.d0 * ij2 + 1.d0)
     :        * (ered - elevh(i))
         do ii=i,nlevelh
            if (twmol) then
              ij2p = mod(jlevh(ii),10)
            else
              ij2p = 0.d0
            endif
            fff = iflevh(ii) + fhspin
            dencol = (2.d0 * fff + 1.d0) * (2.d0 * ij2p + 1.d0)
     :           * (ered - elevh(ii))
            sigma(i,ii) = sigma(i,ii) * fak / denrow
            if (i.ne.ii) then
               sigma(ii,i) = sigma(ii,i) * fak / dencol
            end if
         end do
      end do
*
*     print out cross sections
*
      if (.not. twmol) then
        write(6,1001)
        write(11,1001)
 1001   format('%',5x,'E(CM-1)',5x,'JI',5x,'INI',3x,'FI',6x,'JF',5x,
     :       'INF',3x,'FF',6x,'CROSS SECTION (ANG^2)')
      else
        write(6,2001)
        write(11,2001)
 2001   format(/'%',4x,'E(CM-1)',5x,'JI',5x,'INI',3x,'FI',4x,'J2',6x,
     :       'JF',5x,'INF',
     :       3x,'FF',4x,'J2P',6x,'CROSS SECTION (ANG^2)')
      endif
      ee = ered * econv
      do 90 i=1,nlevelh
        do 85 ii=1,nlevelh
          if (.not. twmol) then
            xj = jlevh(i) + fspin
            xf = iflevh(i) + fhspin
            xjp = jlevh(ii) + fspin
            xfp = iflevh(ii) + fhspin
            if (sigma(i,ii).ne.0.d0) then
              write(6,1002) ee,xj,inlevh(i),xf,
     :             xjp,inlevh(ii),xfp,sigma(i,ii)
                write(11,1002) ee,xj,inlevh(i),xf,
     :           xjp,inlevh(ii),xfp,sigma(i,ii)
 1002         format(f12.3,f8.1,i6,f6.1,3x,f6.1,i6,f6.1,
     :             5x,1pe15.4)
            end if
          else
            xj = (jlevh(i)/10) + fspin
            ij2 = mod(jlevh(i),10)
            xf  = iflevh(i) + fhspin
            xjp = (jlevh(ii)/10) + fspin
            ij2p = mod(jlevh(ii),10)
            xfp  = iflevh(ii) + fhspin
            if (sigma(i,ii).ne.0.d0) then
              write(6,2002) ee,xj,inlevh(i),xf,ij2,
     :             xjp,inlevh(ii),xfp,ij2p,sigma(i,ii)
              write(11,2002) ee,xj,inlevh(i),xf,ij2,
     :             xjp,inlevh(ii),xfp,ij2p,sigma(i,ii)
 2002         format(f12.3,f8.1,i6,f6.1,i6,3x,f6.1,i6,f6.1,
     :             i6,5x,1pe15.4)
            endif
          endif
 85     continue
 90   continue
*
 1000 goto 4021

 4011 write(6,4031)
 4031 format('Problem with Tmatr/i de-allocaton')

 4021 deallocate(exsmtn)
 4010 deallocate(exsmtp)
 4015 deallocate(j12)
 4009 deallocate(l)
 4008 deallocate(in)
 4007 deallocate(j)
 4006 deallocate(si)
 4005 deallocate(sr)
 4004 deallocate(simag)
 4003 deallocate(sreal)
 4002 deallocate(length)
 4001 deallocate(sigma)
 4000 close(1)
      close(11)
      if (ialloc .ne. 0) write (6, 4100)
 4100 format (' *** INSUFFICIENT MEMORY OR SMT FILE CORRUPTED. ***')
      return
      end
*     -------------------------------------------------------------
