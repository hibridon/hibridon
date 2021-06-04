* --------------------------------------------------------------------
*  subroutine to determine angular coupling potential for collision
*  of a 3sigma- diatomic with a closed-shell (1sigma+) diatomic.
*  the second molecule can be homonuclear.
*
*  this subroutine is an extension of ba2sg1sg, which treats 2sigma -
*  1sigma molecules
*
*  pot matrix elements are computed assuming the angular dependence defined
*  in the appendix of green, jcp 62, 2271 (1975)
*
*  author:  paul dagdigian
*  current revision date:  1-nov-2018 by p.dagdigian
* --------------------------------------------------------------------
c     This module contains (explictly) the number of terms and their
c     indices in the expansion of the PES.  Its contents should be
c     filled in the pot routine.
c
c     This module replaces lammin, lammax, mproj in hibridon and
c     is compiled with ba1sg1sg
c
c      module mod_1sg1sg
c      implicit none
c
c      type lm_type
c      integer :: l1, l2, ltot
c      end type lm_type
c
c      type(lm_type), dimension(:), allocatable :: lms
c      end module mod_1sg1sg
* --------------------------------------------------------------------
      subroutine ba3s1sg (j, l, is, jhold, ehold, ishold, nlevel,
     :     nlevop, e1, sc2, sc3, sc4, rcut, jtot,
     :     flaghf, flagsu, csflag, clist, bastst, ihomo,
     :     nu, numin, jlpar, n, nmax, ntop)
* --------------------------------------------------------------------
*  variables in call list:
*    j:        on return contains combined rotational quantum numbers for each
*              channel.  in terms of the rotational quantum numbers of each
*              molecule we have:  j = 10 * j1 + j2
*    l:        on return contains orbital angular momentum for each
*              channel
*    is:       additional quantum index of each channel:  on return this is
*              set equal to the Fi fine-structure label
*    jhold:    on return contains rotational quantum numbers for each
*              rotational level, irrespective of projection degeneracy and the
*              degeneracy associated with different values of j12
*              note that jhold = 10 * j1hold + j2hold
*    ehold:    on return contains energy in hartrees of each rotational
*              level
*    ishold:   on return this is set the Fi label
*    nlevel:   on return contains number of energetically distinct
*              rotational levels used in channel basis
*    nlevop:   on return contains number of energetically distinct
*              rotational levels used in channel basis which are open
*              asymptotically
*    e1        scratch vector for energy (denoted sc1 in calling statement)
*    sc2:      scratch vector (not used here)
*    sc3:      scratch vector (not used here)
*    sc4:      scratch vector (not used here)
*    rcut:     cut-off point for keeping higher energy channels
*              NOTE:  selection with rcut not implemented in this basis routine
*    jtot:     total angular momentum
*              in cc calculation jtot+1/2 is the total angular momentum
*              in cs calculation jtot is the l-bar quantum number
*    flaghf:   if .true., then system has half-integer spin
*              if .false., then system has integer spin (this is the case
*              here)
*    flagsu:   if .true., then molecule-surface collisons (this variable
*              should always be false in the present case)
*    csflag:   if .true., then coupled-states calculation
*              if .false., the close-coupled calculation
*              NOTE:  CS calculation not implemented in this basis routine
*    clist:    if .true., then quantum numbers and energies listed for
*              each channel
*    bastst:   if .true., then execution terminates after the first call
*              to basis
*    ihomo:    if .true. , then homonuclear molecule
*              if .false., then heteronuclear molecule
*              the homonuclear option is not specifically implemented here
*    nu:       coupled-states projection index
*    numin:    minimum coupled states projection index
*              for cc calculations nu and numin are both set = 0 by calling
*              program
*    jlpar:    total parity of included channels in cc calculation
*              parity=(-1)**jlpar (by definition parity=eps*(-1)**(j1+j2+l)
*    n:        on return equals number of channels
*    nmax:     maximum dimension of arrays
*    ntop:     maximum row dimension of all matrices passed to subroutines
*              propag and soutpt.  ntop is set in basis only if nu = numin
*              otherwise it is unchanged from the value supplied by the
*              calling program
*  variables in common block /cosysr/
*    isrcod:   number of real system dependent variables
*    b1rot:    rotational constant in cm-1 of molecuole 1
*    d1rot:    centrifugal distortion constant in cm-1 of molecule 1
*    gamma:    spin-rotation constant of molecule 1
*    flmbda    spin-spin constant of molecule1
*    b2rot:    rotational constant in cm-1 of molecuole 2
*  variables in common block /cosysi/
*    nscode:   total number of system dependent variables
*    isicod:   numbr of integer system dependent variables
*    j1max:    maximum rotational angular momentum of molecule 1
*    j2min:    minimum rotational angular momentum of molecule 2
*    j2max:    maximum rotational angular momentum of molecule 2
*    ipotsy2:  symmetry of potential.  set to 2 for homonuclear
*              molecule 2, set to 1 for heteronuclear molecule 2
*  variable in common block /cocent/
*    cent:      array containing centrifugal barrier of each channel
*  variable in common block /coeint/
*    eint:      array containing channel energies (in hartree)
*  variable in common block /coj12/
*    j12:       vector resultant of j1 + j2
*  variables in common block /coered/
*    ered:      collision energy in atomic units (hartrees)
*    rmu:       collision reduced mass in atomic units
*               (mass of electron = 1)
*  variable in common block /conlam/
*    nlam:      the number of angular coupling terms actually used
*               nlam is set equal to nterm; see above
*    nlammx:    the maximum number of angular coupling terms
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx
*  variable in common block /cov2/
*    nv2max:    maximum core memory allocated for the v2 matrix
*    v2:        lower triangle of nonzero angular coupling matrix elements
*               only nonzero elements are stored
*  variable in common block /coiv2/
*   iv2:        matrix address of v2 matrix for each non-zero element
*  variables in common block /coconv/
*   econv:      conversion factor from cm-1 to hartree
*   xmconv:     conversion factor from amu to atomic units
*  variable in common block /coatpr/
*   c:          expansion coefficienfs for the rotational/fine-structure
*               wave functions in the symmetrized basis:
*               (1) sigma=1, eps=+1; (2) sigma=1, eps=-1, (3) sigma = 0, eps=+1
*  variable in common block /coatp1/
*   ctemp:      temporary storage for rot. e.fn coeffs.
*  variable in common block /coatp2/
*   chold:      ditto
*
*  subroutines called:
*   v3sgsg:    returns molecule-molecule angular coupling coefficient for
*              particular choice of channel index
* --------------------------------------------------------------------
      use mod_1sg1sg
      implicit double precision (a-h,o-z)
      logical ihomo, flaghf, csflag, clist, flagsu, bastst
      include "common/parbas"
      include "common/parbasl"
      common /coipar/ iiipar(9), iprint
      common /cosysi/ nscode, isicod, j1max, j2min, j2max,
     :  ipotsy2
      common /cosysr/ isrcod, junkr, b1rot, d1rot, flmbda, gamma, b2rot
      common /coselb/ ibasty
      common /cov2/ nv2max, junkv, v2(1)
      common /coiv2/ iv2(1)
      common /conlam/ nlam, nlammx, lamnum(1)
      common /cocent/ cent(1)
      common /coeint/ eint(1)
      common /coj12/  j12(1)
      common /coered/ ered, rmu
      common /coconv/ econv, xmconv
      common /coatpr/ c(1)
      common /coatp1/ ctemp(1)
      common /coatp2/ chold(1)
      dimension j(1), l(1), is(1), jhold(1), ehold(1),
     :    e1(1), sc2(1), sc3(1), sc4(1), ishold(1)
      dimension j1(2000),is1(2000), fnn1(3,3), fnn2(3,3),
     :  hss(3,3), hsr(3,3), ham(3,3), u(3,3), uu(3,3), wsig(3,3),
     :  eig(3), eigkp(3)
      data zero, one, two, three, twthr, sq2 / 0.d0, 1.d0, 2.d0, 3.d0,
     :  0.666666666666667d0,  1.414213562373095/
      data spin, ss1 / 1.d0, 2.d0 /
      data uu / 1, 1, 0, 1, -1, 0, 0, 0, 1.414213562373095 /
*  check for consistency in the values of flaghf and csflag
      if (flaghf) then
        write (9, 5)
5       format (' *** FLAGHF = .TRUE. FOR TRIPLET SYSTEM; ABORT ***' )
        if (bastst) then
          return
        else
          call exit
        end if
      end if
      if (flagsu) then
        write (6, 10)
        write (9, 10)
10      format ('  *** FLAGSU = .TRUE. FOR',
     :         ' MOLECULE-MOLECULE COLLISION; ABORT ***')
        if (bastst) then
          return
        else
          call exit
        end if
      end if
      if (csflag) then
        write (6, 25)
        write (9, 25)
25      format (' *** CSFLAG=.TRUE.; NOT ALLOWED IN 3SIG-1SIG BASIS;',
     :          ' ABORT ***')
        if (bastst) then
          return
        else
          call exit
        end if
      end if
*  check that ipotsy2 is either 1 or 2
      if (ipotsy2.lt.1 .or. ipotsy2.gt.2) then
        write (6, 27) ipotsy2
        write (9, 27) ipotsy2
27      format (' *** IPOTSY2 .NE. 1 .OR.2;',
     :          ' ABORT ***')
        if (bastst) then
          return
        else
          call exit
        end if
      end if
      if (bastst) then
        write (6,  90) rmu * xmconv, ered * econv, jtot, jlpar
        write (9,  90) rmu * xmconv, ered * econv, jtot, jlpar
90      format(/2x, ' **  CC  3SIGMA-1SIGMA ** RMU=',f10.6,
     :     '  E=', f9. 3, 2x, 'JTOT=', i4,
     :     '  JLPAR=', i2)
      end if
*
*  assign quantum numbers and energies for rotational levels for molecule 1
*  one level for j = 0
      n1 = 1
      j1(n1) = 0
      is1(n1) = 1
      fj = j1(n1)
      e1(n1) = b1rot * (fj * (fj + one) + ss1)
     :  - twthr * flmbda * ss1 - gamma * ss1
      ctemp(1) = zero
      ctemp(2) = zero
      ctemp(3) = one
      do jj1 = 1, j1max
        fj = jj1
        off = -sqrt(two * fj * (fj + one))
* set up matrices in signed sigma basis:  1+, 1-, 0
        fnn1(1,1) = fj * (fj + one)
        fnn1(1,2) = zero
        fnn1(1,3) = off
        fnn1(2,1) = zero
        fnn1(2,2) = fnn1(1,1)
        fnn1(2,3) = off
        fnn1(3,1) = off
        fnn1(3,2) = off
        fnn1(3,3) = fnn1(1,1) + two
*  evaluate nn1 * nn1 (matrix multiply)
        call dgemm ('n','n',3,3,3,1.d0,fnn1,3,fnn1,3,0.d0,fnn2,3)
*  spin-spin matrix in signed sigma basis
        do ii = 1, 3
          do jj = 1, 3
            hss(ii,jj) = zero
            if (ii .eq. jj) then
              if (ii .lt. 3) then
                hss(ii,jj) = twthr * (three - ss1)
              else
                hss(ii,jj) = - twthr * ss1
              end if
            end if
          end do
        end do
*  spin-rotation matrix in signed sigma basis
        off = 0.5d0 * sqrt( fj * (fj + one) * ss1)
        hsr(1,1) = one - ss1
        hsr(1,2) = zero
        hsr(1,3) = off
        hsr(2,1) = zero
        hsr(2,2) = hsr(1,1)
        hsr(2,3) = off
        hsr(3,1) = off
        hsr(3,2) = off
        hsr(3,3) = -ss1
*  set up Hamiltonian in signed sigma basis
        do ii = 1, 3
          do jj = 1, 3
            ham(ii,jj) = b1rot * fnn1(ii,jj) - d1rot * fnn2(ii,jj)
     :        + flmbda * hss(ii,jj) + gamma * hsr(ii,jj)
          end do
        end do
*  transform to symmetrized basis:  (1) sigma=1, eps=+1; (2) sigma=1, eps=-1;
*  sigma = 0, eps=+1
        u = uu / sq2
        call dgemm ('n','n',3,3,3,1.d0,u,3,ham,3,0.d0,wsig,3)
        call dgemm ('n','n',3,3,3,1.d0,wsig,3,u,3,0.d0,ham,3)
*  diagonalize hamiltonian
        lwork = 144
        call dsyev('V','L',3,ham,3,eig,work,lwork,ierr)
        n1 = n1 + 1
        e1(n1) = eig(1)
        j1(n1) = jj1
        is1(n1) = 1
        n1 = n1 + 1
        e1(n1) = eig(2)
        j1(n1) = jj1
        is1(n1) = 2
        n1 = n1 + 1
        e1(n1) = eig(3)
        j1(n1) = jj1
        is1(n1) = 3
*  make sure sign of (2,2) element of eignvector matrix is positive (i.e. +1)
        if (ham(2,2) .lt. zero) then
          do kk = 1, 3
            do ll = 1,3
              ham(kk,ll) = - ham(kk,ll)
            end do
          end do
        end if
*  save eigenvectors
        do jj = 1, 3
          do kk = 1, 3
           isub = (n1 + jj - 3 - 1) * 3 + kk
           ctemp(isub) = ham(kk,jj)
         end do
        end do
      end do
*
*  now add rotational levels of molecule 2 and shift zero
*  of energy to that of the (j1=1,F1;j2=j2min) pair
      nlevel = 0
      ezero = e1(2) + b2rot * j2min * (j2min + one)
      do jj1 = 1, n1
        do j2 = j2min, j2max, ipotsy2
          nlevel = nlevel + 1
          jhold(nlevel) = 10 * j1(jj1) + j2
          ishold(nlevel) = is1(jj1)
          ehold(nlevel) = (e1(jj1) + b2rot * j2 * (j2 + one)
     :      - ezero) / econv
          do kk = 1, 3
            isub = (nlevel - 1) * 3 + kk
            isub2 = (jj1 - 1) * 3 + kk
            chold(isub) = ctemp(isub2)
          end do
          end do
      end do
*
*  order rotational levels in order of increasing energy
*  sort algorithm from Numerical Recipes
      if (nlevel.gt. 1) then
        do 212 jj = 2, nlevel
          ekeep = ehold(jj)
          jkeep = jhold(jj)
          ikeep = ishold(jj)
          do kk = 1, 3
            isub = (jj - 1) * 3 + kk
            eigkp(kk) = chold(isub)
          end do
          do 211 ii = jj-1, 1, -1
            if (ehold(ii).le.ekeep) goto 210
            ehold(ii+1) = ehold(ii)
            jhold(ii+1) = jhold(ii)
            ishold(ii+1) = ishold(ii)
            do kk = 1, 3
              isub = (ii + 1 - 1) * 3 + kk
              isub2 = (ii - 1) * 3 + kk
              chold(isub) = chold(isub2)
            end do
211       continue
          ii = 0
210       ehold(ii+1) = ekeep
          jhold(ii+1) = jkeep
          ishold(ii+1) = ikeep
          do kk = 1, 3
            isub = (ii + 1 - 1) * 3 + kk
            chold(isub) = eigkp(kk)
          end do
212     continue
      end if
*
*  form list of all energetically open rotational levels included in the
*  calculations and their energies
      nlevop = 0
      do i = 1, nlevel
        if (ehold(i) .le. ered) then
          nlevop = nlevop + 1
        end if
      end do
*
*  list levels if requested
      if (bastst) then
        write(6,221)
        write(9,221)
221     format(/'  LEVEL LIST SORTED BY ENERGY'/
     :    '   N    J1    Fi    J2     EINT(CM-1)       COEFFS')
        do i = 1, nlevel
          jj1 = jhold(i)/10
          j2 = mod(jhold(i),10)
          ecm = ehold(i) * econv
          isub = (i - 1) * 3
          write(6,222) i,jj1,ishold(i),j2,ecm,
     :      (chold(isub + kk), kk=1,3)
          write(9,222) i,jj1,ishold(i),j2,ecm,
     :      (chold(isub + kk), kk=1,3)
222       format(i4, 3i6, f14.3, 4x, 3f7.3)
        end do
      end if
*
*  set up CC scattering channel basis
      n = 0
      xjtot = jtot
      do i = 1, nlevel
        jj1 = jhold(i)/10
        j2 = mod(jhold(i),10)
        j12mn = abs(jj1 - j2)
        j12mx = jj1 + j2
        j12min = j12mn
        j12max = j12mx
        do j12i = j12min, j12max
          lmin = abs(jtot - j12i)
          lmax = jtot + j12i
          do li = lmin, lmax
c  parity of 3sigma molecule is ieps1 * (-1) ^ (j1 - 1/2)
            ix = (- 1) ** (jj1 + j2 - jtot + li)
            if (ishold(i) .eq. 2) ix = -ix
            if (ix .eq. jlpar) then
              n = n + 1
              if (n .le. nmax) then
                j(n) = jhold(i)
                is(n) = ishold(i)
                l(n) = li
                j12(n) = j12i
                eint(n) = ehold(i)
                cent(n) = li * (li + 1.d0)
                do kk = 1, 3
                  isub = (n - 1) * 3 + kk
                  isub2 = (i - 1) * 3 + kk
                  c(isub) = chold(isub2)
                end do
              else
                write (9, 230) n, nmax
                write (6, 230) n, nmax
230             format(/' *** NCHANNELS=', i5,' .GT. ',
     :              'MAX DIMENSION OF ',i5,'; ABORT')
                if (bastst) then
                  return
                else
                  call exit
                end if
              end if
            end if
          end do
        end do
      end do
      if (n .eq. 0) return
***      if (nn .eq. numin) then
        ntop = max(n, nlevop)
*  ntop is the maximum row dimension of all matrices passed in the
*  call list of subroutines propag and soutpt.
*  for fps make sure that this is an odd number, for faster bank access.
*  this has no effect on vax or cray
        if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax) then
          ntop = ntop + 1
        else
          if (n .gt. ntop) then
            write(6,160) nn, n, ntop
            write(9,160) nn, n, ntop
160         format (' *** NCH = ',i3,' AT NU = ',i2,' .GT. NTOP = ',i3,
     :          '; ABORT ***',/,
     :      '      CHECK RCUT')
            call exit
          end if
***        end if
      end if
*
*  now list channels if requested
      if (bastst) then
        write (6, 200)
        write(9, 200)
200     format(/'  CHANNEL LIST SORTED BY ENERGY'/
     :      '   N    J1     Fi    J2    J12      L   ',
     :      ' EINT(CM-1)     COEFFS')
        do i = 1, n
          ecm = eint(i) * econv
          jj1 = j(i)/10
          fj1 = jj1
          j2 = mod(j(i),10)
          fj12 = j12(i)
          isub = (i - 1) * 3
          write (6, 220) i, fj1, is(i), j2, fj12, l(i), ecm,
     :      (c(isub + kk), kk=1,3)
          write (9, 220) i, fj1, is(i), j2, fj12, l(i), ecm,
     :      (c(isub + kk), kk=1,3)
220       format (i4, f7.1, 2i6, f7.1, i7, f13.3, 2x, 3f7.3)
        end do
      end if
c
c     Calculate coupling matrix elements
      if (bastst .and. iprint.eq.2) then
        write (6, 285)
        write (9, 285)
285     format (/' ILAM  L1   L2  LTOT   ICOL IROW      I',
     :    '      IV2         VEE')
      end if
      i = 0
      lamsum = 0
      do 400 ilam = 1, nlam
c     ilam denotes a particular L1,L2,L term
        ll1 = lms(ilam)%l1
        ll2 = lms(ilam)%l2
        lltot = lms(ilam)%ltot
        inum = 0
        do 355 icol = 1, n
          j1c = j(icol)/10
          j2c = mod(j(icol),10)
          ifc = is(icol)
          j12c = j12(icol)
          lc = l(icol)
          do 350 irow = icol, n
            ij = ntop * (icol - 1) + irow
            j1r = j(irow)/10
            j2r = mod(j(irow),10)
            ifr = is(irow)
            j12r = j12(irow)
            lr = l(irow)
c     initialize potential to zero
            vee = 0.d0
            call v3sgsg(irow,icol,j1r,ifr,j2r,j12r,lr,
     :          j1c,ifc,j2c,j12c,lc,jtot,ll1,ll2,lltot,vee)
            if (abs(vee) .gt. 1.d-8) then
              i = i + 1
              if (i .le. nv2max) then
                inum = inum + 1
                v2(i) = vee
                iv2(i) = ij
                if (bastst .and. iprint.ge.2) then
                  write (6, 290) ilam, ll1, ll2, lltot,
     :                icol, irow, i, iv2(i), vee
                  write (9, 290) ilam, ll1, ll2, lltot,
     :                icol, irow, i, iv2(i), vee
290               format (i4, 3i5, 2x, 2i5, 2i8, e20.7)
                end if
              end if
            end if
350       continue
355     continue
        if (i .le. nv2max) lamnum(ilam) = inum
        if (bastst) then
          write (6, 370) ilam, ll1, ll2, lltot, lamnum(ilam)
          write (9, 370) ilam, ll1, ll2, lltot, lamnum(ilam)
370       format ('ILAM=',i4,'  L1=',i3,' L2=',i3,' LTOT=',i3,
     :      '  LAMNUM(ILAM) = ',i7)
        end if
        lamsum = lamsum + lamnum(ilam)
400   continue
      if (i .gt. nv2max) then
        write (6, 410) i, nv2max
        write (6, 410) i, nv2max
410     format (' *** NUMBER OF NONZERO V2 ELEMENTS = ',i6,
     :      ' .GT. NV2MAX=',i6,'; ABORT ***')
      end if
      if (bastst) then
        write (6, 420) lamsum
        write (9, 420) lamsum
420     format (' *** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS',
     :      i8)
      end if
      return
      end
* --------------------------------------------------------------------
      subroutine v3sgsg(irow,icol,j1r,ifr,j2r,j12r,lr,
     :  j1c,ifc,j2c,j12c,lc,jtot,l1,l2,ltot,vee)
*  subroutine to calculate v-lambda matrices for close-coupled
*  treatment of collisions of a 3sigma and a 1sigma molecule
*  electronic state with a (different) molecule in a 1sigma state
*
*  matrix elements adapted from 1sigma-1sigma pot matrix elements
*  originally derived by s. green, jcp 62, 2271 (1973)
*
*  subroutines called:
*     xf3j, xf6j, xf9j
*
*  written  by p. dagdigian
*  current revision date:  1-aug-2018 (p.dagdigian)
* --------------------------------------------------------------------
      use mod_1sg1sg
      implicit double precision (a-h,o-z)
      common /coatpr/ c(1)
      data one, two/ 1.d0, 2.d0/,
     :  sq4pi3 / 44.546623974653656d0 /
      vee = 0.d0
      xj1r = j1r
      xj2r = j2r
      xj12r = j12r
      xlr = lr
      xj1c = j1c
      xj2c = j2c
      xj12c = j12c
      xlc = lc
      xjtot = jtot
      xl1 = l1
      xl2 = l2
      xltot = ltot
      isubr = (irow - 1) * 3
      isubc = (icol - 1) * 3
*
      x1 = xf3j(xltot,xlr,xlc,0.d0,0.d0,0.d0)
      if (x1 .eq. 0.d0) return
*
      x2 = xf3j(xl2,xj2r,xj2c,0.d0,0.d0,0.d0)
      if (x2 .eq. 0.d0) return
*
      x6 = xf6j(xlr,xlc,xltot,xj12c,xj12r,xjtot)
      if (x6 .eq. 0.d0) return
*
      x9 = xf9j(xj12r,xj2r,xj1r,xj12c,xj2c,xj1c,
     :  xltot,xl2,xl1)
      if (x9 .eq. 0.d0) return
*
      iph1 = xl1 + xj1r + xj1c
      x3j1 = xf3j(xl1,xj1r,xj1c,0.d0,1.d0,-1.d0)
      x3j0 = xf3j(xl1,xj1r,xj1c,0.d0,0.d0,0.d0)
*
*  basis:  (1) sigma = 1+, (2) sigma = 1-, (3) sigma = 0
      c1r = c(isubr + 1)
      c0r = c(isubr + 3)
      c1c = c(isubc + 1)
      c0c = c(isubc + 3)
*
      if (ifr .eq. 2) then
*  here for F2 bra
        if (ifc .eq. 2) then
*    here for F2 ket
          if ( (-1.d0) ** iph1 .eq. 1.d0) then
            vee = x3j1
          end if
        else
*    here for F1 or F3 ket
          if ( (-1.d0) ** iph1 .eq. -1.d0) then
            vee = c1c * x3j1
          end if
          vee = vee - c0c * x3j0
        end if
      else
*  here for F1 or F3 bra
        if (ifr .eq. 2) then
*    here for F2 ket
          if ( (-1.d0) ** iph1 .eq. -1.d0) then
            vee = c1r * x3j1
          end if
          vee = vee - c0r * x3j0
        else
*    here for F1 or F3 ket
          if ( (-1.d0) ** iph1 .eq. 1.d0) then
            vee = c1r *c1c * x3j1
          end if
          vee = vee - c0r * c0c * x3j0
        end if
      end if
*
      facj = (two*xltot + one) * sqrt((two*xl1 + one) * (two*xl2 + one)
     :    * (two*xj1r + one) * (two*xj2r + one) * (two*xj12r + one)
     :    * (two*xlr + one)
     :    * (two*xj1c + one) * (two*xj2c + one) * (two*xj12c + one)
     :    * (two*xlc + one) )
      iph = xjtot + xj1c + xj2c + xj12r - 1.d0
      phase = (-1.d0) ** iph
      vee = vee * phase * facj * x1 * x2 * x6 * x9 / sq4pi3
*
      return
      end
* --------------------------------eof---------------------------------

