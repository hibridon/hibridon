* ----------------------------------------------------------------------
*  subroutine to determine angular coupling potential for collision
*  of an asymmetric top molecule having no symmetry elements, e.g. a chiral
*  molecule, with a structureless atom or with an uncorrugated surface
*
*  author:  paul dagdigian
*  revision:  6-jun-2013 by q. ma (fix a bug in counting anisotropic
*     terms)
*  revision:  30-jan-2018 by pjd -- switch row and column in assembling
*     eigenvectors
*  revision:  16-jan-2019 by  pjd -- extend baastp to treat collisions
*     of chiral molecules
*
*  author:  paul dagdigian
*  current revision date:  20-feb-2019 by p.dagdigian
* --------------------------------------------------------------------
c     This module contains (explictly) the number of terms and their
c     indices in the expansion of the PES.  Its contents should be
c     filled in the pot routine.
c
c     This module replaces lammin, lammax, mproj in hibridon
c
      module mod_chiral
      implicit none
c
      type lm_type
      integer :: l1, m1
      end type lm_type
c
      type(lm_type), dimension(:), allocatable :: lms
      end module mod_chiral
* --------------------------------------------------------------------
      subroutine baastp2(j, l, is, jhold, ehold, ishold, nlevel, nlevop,
     :                  etemp, fjtemp, fktemp, fistmp,
     :                  rcut, jtot, flaghf, flagsu,
     :                  csflag, clist, bastst, ihomo, nu, numin, jlpar,
     :                  n, nmax, ntop)
* --------------------------------------------------------------------
*  variables in call list:
*    j:        on return contains rotational quantum number for each
*              channel
*    l:        on return contains orbital angular momentum for each
*              channel
*    is:       on return contains symmetry index (ieps * kp) for each
*              channel
*  Note that we have adopted the following convention for the symmetry
*  index "is" so that on return is = 0/1, with the symmetric top basis
*  functions written as [|j,kp> + (-1)^is*|j,-kp)\, respectively.
*    jhold:    on return contains rotational quantum numbers for each
*              rotational level
*    ehold:    on return contains energy in hartrees of each rotational
*              level
*    ishold:   on return contains symmetry index of each energetically
*              distinct level
*    nlevel:   on return contains number of energetically distinct
*              rotational levels used in channel basis
*    nlevop:   on return contains number of energetically distinct
*              rotational levels used in channel basis which are open
*              asymptotically
*    etemp:    scratch array used to create channel list
*    fjtemp:   scratch array used to create channel list
*    fktemp:   scratch array used to create channel list
*    fistmp:   scratch array used to create channel list
*    rcut:     cut-off point for keeping higher energy channels
*              NOTE:  selection with rcut not implemented in this basis routine
*    jtot:     total angular momentum
*              in cc calculation jtot is the total angular momentum
*              in cs calculation jtot is the l-bar quantum number
*    flaghf:   if .true., then system has half-integer spin
*              if .false., then system has integer spin
*              NOTE:  flaghf should be f for this basis routine
*    flagsu:   if .true., then molecule-surface collisons
*              NOTE:  flagsu should be f for this basis routine
*    csflag:   if .true., then coupled-states calculation
*              if .false., the close-coupled calculation
*              NOTE:  csflag should be f for this routine.
*              CS approximation not implemented here
*    clist:    if .true., then quantum numbers and energies listed for
*              each channel
*    bastst:   if .true., then execution terminates after the first call
*              to basis
*    ihomo:    if .true., then the molecule posseses interchange symmetry
*              NOTE:  ihomo not used here
*    nu:       coupled-states projection index
*    numin:    minimum coupled states projection index
*              for cc calculations nu and numin are both set = 0 by calling
*              program
*    jlpar:    total parity of included channels in cc calculation
*              only those channels are included for which
*                  (-1)**(is+j+kp+l-jtot)=jlpar
*              where parity designates the parity of the molecular state
*              in cs calculation jlpar is set equal to 1 in calling program
*    n:        on return equals number of channels
*    nmax:     maximum dimension of arrays
*    ntop:     maximum row dimension of all matrices passed to subroutines
*              propag and soutpt.  ntop is set in basis only if nu = numin
*              otherwise it is unchanged from the value supplied by the
*              calling program
*  variables in common bloc /cosysr/
*    isrcod:   total number of real system dependent variables
*    arot:     A rotational constant for the asymmetric top
*    brot:     B rotational constant for the asymmetric top
*    crot:     C rotational constant for the asymmetric top
*    emax:     the maximum rotational energy (in cm^-1) for a channel to be
*              included in the basis
*  variables in common block /cosysi/
*    nscode:   total number of variable names which are passed to HINPUT
*              nscod must equal isrcod + isicod + 3
*    isicod:   total number of integer system dependent variables
*    nterm:    number of different associated legendre terms in
*              expansion of potential
*    numpot:   the number of the potential used, this variable is passed
*              to the pot subroutine
*    jmax:     the maximum rotational angular momentum for the asymmetric top
*  variable in common block /cocent/
*    cent:      array containing centrifugal barrier of each channel
*  variable in common block /coeint/
*    eint:      array containing channel energies (in hartree)
*               the zero of energy is assumed to be the 0(0,0) level
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
*  variables in common block /cov2/
*    nv2max:    maximum core memory allocated for the v2 matrix
*    v2:        lower triangle of nonzero angular coupling matrix elements
*               stored in packed column form that is (1,1), (2,1), (3,1) ...
*               (n,1), (2,2), (3,2) ... (n,2), etc.
*  variable in common block /coiv2/
*   iv2:        row+column index of v2 matrix for each non-zero element
*  variables in common block /coconv/
*   econv:      conversion factor from cm-1 to hartrees
*   xmconv:     converson factor from amu to atomic units
*  variables in common block /coatpi/
*   narray:     maximum size of asymmetric top basis fn expansion
*               set to 12 in himain (suitable for j <= 40)
*   isiz:       length of eigenfunction expansion for each rot. level
*  variable in common block /coatp3/
*   isizh:      temporary storage for length (related to isiz)
*  variable in common block /coatpr/
*   c:          expansion coefficients for asymmetric top rotor wave fns.
*               in a symmetrized symmetric top basis.  these basis functions
*               are given by green, jcp 64, 3463 (1976) as:
*                  |j,k,m,s> = [2*(1+delta(k,0)]^-1 *
*                      (|j,k,m> + (-1)*s * |j,-k,m>)
*               (note that green uses eps = (-1)*s for the symmetry index.)
*               in this basis, the asymmetric top hamiltonian block diagonalizes
*               into 4 groups:  (1) k even, s = +1, (2) k even, s= -1,
*               (3) k odd, s = =1, and (4) k odd, s = _1.
*               the expansion coefficients stored in the array c are:
*               c(k=0), c(k=2), ... c(k=j) for even k, and
*               c(k=1), c(k=3), ... c(k=j) for odd k.
*               the levels are denoted in the channel list by j and the
*               symmetry index is.  the prolate-limit projection quantum number kp
*               equals abs(is), and s is the sign of is.  the oblate limit projection
*               quantum number ko can be obtained as follows:
*                  if is >= 0, then ko = j - kp
*                  if is < 0, then ko = j + 1 - kp
*               if bastst=t, the values of j, is, kp, and ko are printed out, as
*               well as the values of the expansion coefficients
*  variable in common block /coatp1/
*   ctemp:      temporary storage for rot. e.fn coeffs.
*  variable in common block /coatp2/
*   chold:      ditto
*  subroutines called:
*   vlmtp2:     returns angular coupling coefficient for particular
*               choice of channel index
*   prmtp2:     computes primitive cc and cs v-lambda matrix elements
*               between symmetric top basis fns.
*   rotham:     computes matrix elements of asymmetric top hamiltonian
* --------------------------------------------------------------------
      use mod_chiral
      implicit double precision (a-h,o-z)
      logical flaghf, csflag, clist, flagsu, ihomo, bastst
      character*1 slab
      include "common/parbas"
      include "common/parbasl"
      common /cosysi/ nscode, isicod, nterm, numpot, jmax
      common /coipar/ iiipar(9), iprint
      common /cosysr/ isrcod, junkr, arot, brot, crot, emax
      common /coselb/ ibasty
      common /cov2/ nv2max, junkv, v2(1)
      common /coiv2/ iv2(1)
      common /conlam/ nlam, nlammx, lamnum(1)
      common /cocent/ cent(1)
      common /coeint/ eint(1)
      common /coered/ ered, rmu
      common /coconv/ econv, xmconv
      common /coatpi/ narray, isiz(1)
      common /coatp3/ isizh(1)
      common /coatpr/ c(1)
      common /coatp1/ ctemp(1)
      common /coatp2/ chold(1)
      dimension j(1), l(1), is(1), jhold(1), ehold(1),
     :          ishold(1), etemp(1), fjtemp(1), fktemp(1),
     :          fistmp(1)
*  scratch arrays for computing asymmetric top energies and wave fns.
      dimension e(narray,narray), eig(narray), vec(narray,narray),
     :  sc1(narray), sc2(narray), work(288)
*
      zero = 0.d0
      two = 2.d0
*  check for consistency in the values of flaghf and csflag
      if (flaghf) then
        write (6, 5)
        write (9, 5)
5      format (' *** FLAGHF = .TRUE. FOR SINGLET SYSTEM; ABORT ***')
        stop
      end if
      if (csflag) then
        write (6, 6)
        write (9, 6)
6     format
     :   ('  *** CSFLAG = .TRUE. FOR ASYMMETRIC TOP',
     :     ' COLLISIONS; NOT IMPLEMENTED.  ABORT ***')
        stop
      end if
      if (flagsu) then
        write (6, 10)
        write (9, 10)
10      format ('  *** FLAGSU = .TRUE. FOR ASYMMETRIC TOP',
     :     ' COLLISIONS; NOT IMPLEMENTED.  ABORT ***')
        call exit
      end if
      xjtot = jtot

      if (bastst) then

        write (6,85) rmu * xmconv, arot, brot, crot,
     :    ered * econv, jtot, jlpar
        write (9,85) rmu * xmconv, arot, brot, crot,
     :    ered * econv, jtot, jlpar
85      format(/,' **  CC CHIRAL ASYMMETRIC TOP **',
     :    /,'     RMU=', f9.4,'  AROT=', f7.3, '  BROT=',f7.3,
     :    '  CROT=',f7.3,/,  '     E=', f7.2, '  JTOT=', i4, 2x,
     :    ' JLPAR=', i2)

      end if

*
*  first set up list of all j(kp,ko) states included in basis
      nlist = 0
      do ji = 0, jmax
* set up and diagonalize rotational hamiltonian for each value of ji
* hamiltonian block diagonalize into 1-4 subblocks
*
* E+ submatrix
        isize = (ji + 2)/2
        do mm = 1, isize
          do nn = 1, isize
            e(mm,nn) = 0.d0
          end do
        end do
        e(1,1) = rotham(ji,0,ji,0)
        if (isize .gt. 1) then
          do mm = 2, isize
            kk = 2*mm - 2
            e(mm,mm) = rotham(ji, kk, ji, kk)
            e(mm,mm-1) = rotham(ji,kk,ji,kk-2)
            e(mm-1,mm) = e(mm,mm-1)
          end do
          e(2,1) = sqrt(2.d0)*e(2,1)
          e(1,2) = e(2,1)
        end if
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
        lwork = 144
        call dsyev('V','L',isize,e,narray,eig,work,lwork,ierr)
#endif
#if !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
      write (6,1225)
1225    format('  *** DSYEV NOT ACCESSIBLE IN THIS VERSION OF CODE.')
        call exit
#endif
        do mm = 1, isize
          nlist = nlist + 1
*  eigenvalues stored in order of increasing energy (lowest KP first)
          etemp(nlist) = eig(mm)
          fjtemp(nlist) = ji
          fktemp(nlist) = 2*mm - 2
          fistmp(nlist) = 1
*  eigenfunctions expressed over basis KP = 0, 2, 4, ...
          nbas = int(fktemp(nlist)/2) + 1
          do nn = 1, isize
            isub = (nlist - 1)*narray + nn
*  make sure coeff for prolate k is positive
            if (e(mm, nbas) .gt. zero) then
*              ctemp(isub) = e(mm, nn)
              ctemp(isub) = e(nn, mm)
            else
*              ctemp(isub) = -e(mm, nn)
              ctemp(isub) = -e(nn, mm)
            end if
          end do
        end do
*
* E- submatrix
        isize = ji/2
        if (isize .gt. 0) then
          do mm = 1, isize
            do nn = 1, isize
              e(mm,nn) = 0.d0
            end do
          end do
          e(1,1) = rotham(ji,2,ji,2)
          if (isize .gt. 1) then
          do mm = 2, isize
            kk = 2*mm
            e(mm,mm) = rotham(ji, kk, ji, kk)
            e(mm,mm-1) = rotham(ji,kk,ji,kk-2)
            e(mm-1,mm) = e(mm,mm-1)
          end do
          end if
          call dsyev('V','L',isize,e,narray,eig,work,lwork,ierr)
          do mm = 1, isize
            nlist = nlist + 1
*  eigenvalues stored in order of increasing energy (lowest KP first)
            etemp(nlist) = eig(mm)
            fjtemp(nlist) = ji
            fktemp(nlist) =  2*mm
            fistmp(nlist) = -1
*  eigenfunctions expressed over basis KP = 0, 2, 4, ...
            nbas = int(fktemp(nlist)/2)
            do nn = 1, isize
              isub = (nlist - 1)*narray + nn
*  make sure coeff for prolate k is positive
              if (e(mm, nbas) .gt. zero) then
*                ctemp(isub) = e(mm, nn)
                ctemp(isub) = e(nn, mm)
              else
*                ctemp(isub) = -e(mm, nn)
                ctemp(isub) = -e(nn, mm)
              end if
            end do
          end do
        end if
*
* O+ submatrix
        isize = (ji + 1)/2
        if (isize .gt. 0) then
          do mm = 1, isize
            do nn = 1, isize
              e(mm,nn) = 0.d0
            end do
          end do
          e(1,1) = rotham(ji,1,ji,1) + rotham(ji,1,ji,-1)
          if (isize .gt. 1) then
          do mm = 2, isize
            kk = 2*mm - 1
            e(mm,mm) = rotham(ji, kk, ji, kk)
            e(mm,mm-1) = rotham(ji,kk,ji,kk-2)
            e(mm-1,mm) = e(mm,mm-1)
          end do
          end if
          call dsyev('V','L',isize,e,narray,eig,work,lwork,ierr)
          do mm = 1, isize
            nlist = nlist + 1
*  eigenvalues stored in order of increasing energy (lowest KP first)
            etemp(nlist) = eig(mm)
            fjtemp(nlist) = ji
            fktemp(nlist) =  2*mm - 1
            fistmp(nlist) = 1
*  eigenfunctions expressed over basis KP = 1, 3, 5, ...
            nbas = int((fktemp(nlist) - 1)/2) + 1
            do nn = 1, isize
              isub = (nlist - 1)*narray + nn
*  make sure coeff for prolate k is positive
              if (e(mm, nbas) .gt. zero) then
*                ctemp(isub) = e(mm, nn)
                ctemp(isub) = e(nn, mm)
              else
*                ctemp(isub) = -e(mm, nn)
                ctemp(isub) = -e(nn, mm)
              end if
            end do
          end do
        end if
*
* O- submatrix
        isize = (ji + 1)/2
        if (isize .gt. 0) then
          do mm = 1, isize
            do nn = 1, isize
              e(mm,nn) = 0.d0
            end do
          end do
          e(1,1) = rotham(ji,1,ji,1) - rotham(ji,1,ji,-1)
          if (isize .gt. 1) then
          do mm = 2, isize
            kk = 2*mm - 1
            e(mm,mm) = rotham(ji, kk, ji, kk)
            e(mm,mm-1) = rotham(ji,kk,ji,kk-2)
            e(mm-1,mm) = e(mm,mm-1)
          end do
          end if
          call dsyev('V','L',isize,e,narray,eig,work,lwork,ierr)
          do mm = 1, isize
            nlist = nlist + 1
*  eigenvalues stored in order of increasing energy (lowest KP first)
            etemp(nlist) = eig(mm)
            fjtemp(nlist) = ji
            fktemp(nlist) =  2*mm - 1
            fistmp(nlist) = -1
*  eigenfunctions expressed over basis KP = 1, 3, 5, ...
            nbas = int((fktemp(nlist) - 1)/2) + 1
            do nn = 1, isize
              isub = (nlist - 1)*narray + nn
*  make sure coeff for prolate k is positive
              if (e(mm, nbas) .gt. zero) then
*                ctemp(isub) = e(mm, nn)
                ctemp(isub) = e(nn, mm)
              else
*                ctemp(isub) = -e(mm, nn)
                ctemp(isub) = -e(nn, mm)
              end if
            end do
          end do
        end if
      end do
*
*  now sort this list in terms of increasing energy
      if (nlist .gt. 1) then
        do 120 i1 = 1, nlist - 1
          esave = etemp(i1)
          do 115 i2 = i1 + 1, nlist
            if (etemp(i2) .lt. esave) then
*  state i2 has a lower energy than state i1, switch them
              esave = etemp(i2)
              etemp(i2) = etemp(i1)
              etemp(i1) = esave
              fjsave = fjtemp(i2)
              fjtemp(i2) = fjtemp(i1)
              fjtemp(i1) = fjsave
              fksave = fktemp(i2)
              fktemp(i2) = fktemp(i1)
              fktemp(i1) = fksave
              fissav = fistmp(i2)
              fistmp(i2) = fistmp(i1)
              fistmp(i1) = fissav
*  also move e.fn coeffs (don't worry about size of e.fn)
              do mm = 1, narray
                isub1 = (i1 - 1)*narray + mm
                isub2 = (i2 - 1)*narray + mm
                sc1(mm) = ctemp(isub2)
                ctemp(isub2) = ctemp(isub1)
                ctemp(isub1) = sc1(mm)
              end do
            end if
115       continue
120     continue
      end if
*
*  now set up channel and level list
*  print this list if bastst = .true. or if clist = .true.
      if (bastst .or. clist) then
        write (6, 130)
130     format (/,2x,
     :   'LEVEL LIST SORTED BY ENERGY',/,'   N   J  ',
     :     'IS  KP  KO  S   EINT(CM-1)  COEFFS')
      end if
      n = 0
      nlevel = 0
      do 170  njk = 1, nlist
        ki = fktemp(njk)
        ji = fjtemp(njk)
        isi = fistmp(njk)
*
*  delete state if energy > emax
        if (etemp(njk) .gt. emax) go to 170
*  here if this state is to be included
        nlevel = nlevel + 1
        ehold(nlevel) = etemp(njk) / econv
        jhold(nlevel) = ji
        ishold(nlevel) = isi * ki
*  also move e.fn coeffs
        do mm = 1, narray
          isub = (nlevel - 1)*narray + mm
          isub1 = (njk - 1)*narray + mm
          chold(isub) = ctemp(isub1)
        end do
*
*  print this level if bastst = .true.
*  determine size of wave function expansion for printout
        if (ki .eq. 2*(ki/2)) then
          if (isi .eq. 1) then
*  E+ block
            isize = (ji + 2)/2
          else
*  E- block
            isize = ji/2
         end if
        else
*  O+, O- blocks
            isize = (ji + 1)/2
        end if
        isizh(nlevel) = isize
        if (bastst .or. clist) then
          ecm = ehold(nlevel) * econv
          isub = (nlevel - 1)*narray
          if (isi.eq.1) then
            slab='+'
          else
            slab='-'
          endif
*  compute kp and ko projection quantum numbers
          kp = abs(ishold(nlevel))
          if (ishold(nlevel) .ge. 0) then
            ko = jhold(nlevel) - kp
          else
            ko = jhold(nlevel) + 1 - kp
          end if
          if (isize .gt. 12) then
            write (6, 135) nlevel, jhold(nlevel), ishold(nlevel),
     :        kp, ko, slab, ecm, (chold(isub + mm), mm=1,isize)
135         format (5i4, a3, f10.3, 3x, 6f8.4/36x, 6f8.4/
     :        36x,6f8.4)
          else
            if (isize .gt. 6) then
              write (6, 1351) nlevel, jhold(nlevel), ishold(nlevel),
     :          kp, ko, slab, ecm, (chold(isub + mm), mm=1,isize)
1351           format (5i4, a3, f10.3, 3x, 6f8.4/36x, 6f8.4)
            else
              write (6, 1352) nlevel, jhold(nlevel), ishold(nlevel),
     :          kp, ko, slab, ecm, (chold(isub + mm), mm=1,isize)
1352            format (5i4, a3, f10.3, 3x, 6f8.4)
            end if
          end if
        end if
*
*  here for CC calculations.  first calculate range of orbital angular
*  momentum quantum numbers allowed for this state
*  determine parity of molecular state [Eq. (A21) of S. Green, J. Chem. Phys.
*  73, 2740 (1980) and Townes and Schawlow, Microwave Spectroscopy, Eq. (3-27),
*  p. 64.]  See Eq. (A3) of Green.
          iss = sign(1, ishold(nlevel))
          ipar = (-1) ** (ji + ki) * iss
          lmax = jtot + ji
          lmin = iabs (jtot - ji)
          do 155  li = lmin, lmax
            ix = ipar * (-1) ** (li - jtot)
            if (ix .eq. jlpar) then
              n = n + 1
              if (n .gt. nmax) then
                write (6, 150) n, nmax
                write (9, 150) n, nmax
150             format(/' *** NCHANNELS=', i5,
     :             ' .GT. MAX DIMENSION OF',i5,' ABORT ***')
                stop
              end if
              is(n) = ishold(nlevel)
              j(n) = ji
              eint(n) = ehold(nlevel)
              l(n) = li
              cent(n) = li * (li + 1)
*  move e.fn also
              isiz(n) = isizh(nlevel)
              do mm = 1, narray
                isub = (n - 1)*narray + mm
                isub1 = (nlevel - 1)*narray + mm
                c(isub) = chold(isub1)
              end do
            end if
155       continue
*
170   continue
*
*  also determine number of levels which are open
      nlevop = 0
      do 250  i = 1, nlevel
        if (ehold(i) .le. ered) then
          nlevop = nlevop + 1
        end if
250   continue
*  return if no channels
      if (n .eq. 0) return
      if (nu .eq. numin) then
        ntop = max(n, nlevop)
*  ntop is the maximum row dimension of all matrices passed in the
*  call list of subroutines propag and soutpt.
*  for fps make sure this is an odd number, for faster bank access.
*  this has no effect on vax or cray
        if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax)
     :       ntop = ntop + 1
      else
        if (n.gt.ntop) then
          write (6, 303) nu, n, ntop
          write (9, 303) nu, n, ntop
303       format (' *** NCH = ',i3, ' AT NU = ',i2, ' .GT. NTOP = ',i3,
     :            '; ABORT **',/,
     :    '     CHECK RCUT')
          call exit
        endif
      end if
*  now list channels if requested
      if (bastst) then
        write (6,310)
        write (9,310)
310     format
     :    (/,2x,'CC CHANNEL LIST SORTED BY ENERGY',/,
     :    '   N   J  IS   L   EINT(CM-1)   COEFFS')
        do 330  i = 1, n
          ecm = eint(i) * econv
          if (bastst .or. clist) then
            isize = isiz(i)
            isub = (i - 1)*narray
            if (isize .gt. 12) then
              write (6, 320) i, j(i), is(i), l(i), ecm,
     :          (c(isub + mm), mm=1,isize)
              write (9, 320) i, j(i), is(i), l(i), ecm,
     :          (c(isub + mm), mm=1,isize)
320           format (4i4, f10.3, 3x, 6f8.4/29x, 6f8.4/
     :          29x, 6f8.4)
            else
              if (isize .gt. 6) then
                write (6, 3201) i, j(i), is(i), l(i), ecm,
     :            (c(isub + mm), mm=1,isize)
                write (9, 3201) i, j(i), is(i), l(i), ecm,
     :            (c(isub + mm), mm=1,isize)
3201             format (4i4, f10.3, 3x, 6f8.4/29x, 6f8.4)
              else
                write (6, 3202) i, j(i), is(i), l(i), ecm,
     :            (c(isub + mm), mm=1,isize)
                write (9, 3202) i, j(i), is(i), l(i), ecm,
     :            (c(isub + mm), mm=1,isize)
3202            format (4i4, f10.3, 3x, 6f8.4)
              end if
            end if
          end if
330     continue
      end if
*
*  Calculate coupling matrix elements
*  i counts v2 elements
*  inum counts v2 elements for given lambda
*  ilam counts number of v2 matrices
*  ij is address of given v2 element in present v2 matrix
      if (bastst) then
        write (6, 340)
        write (9, 340)
340     format (/'  ILAM LAMBDA MU    ICOL IROW      I      IV2',
     :    '        VEE')
       end if
      i = 0
      lamsum = 0
*  ilam denotes a particular lambda,mu term
      do 400 ilam = 1, nlam
c     ilam denotes a particular LAMBDA,MU term
        lam = lms(ilam)%l1
        mu = lms(ilam)%m1
c
c  IMPORTANT:  THE MU < 0 TERMS DO NOT CONTRIBUTE TO THE
c  COUPLING MATRIX ELEMENTS
        if (mu .lt. 0) goto 400
c
        inum = 0
        do 355 icol = 1, n
          do 350 irow = icol, n
            ij = ntop * (icol - 1) + irow
c     initialize potential to zero
            vee = 0.d0
            call vchirl(j(irow), l(irow), j(icol), l(icol),
     :        is(irow), is(icol), jtot, lam, mu, vee,
     :        irow, icol)
c     check for nonzro matrix element
            if (abs(vee) .gt. 1.d-10) then
              i = i + 1
              if (i .le. nv2max) then
                inum = inum + 1
                v2(i) = vee
                iv2(i) = ij
                if (bastst .and. iprint.ge.2) then
                  write (6, 390) ilam, lam, mu,
     :                icol, irow, i, iv2(i), vee
                  write (9, 390) ilam, lam,mu,
     :                icol, irow, i, iv2(i), vee
390               format (i6, 2i5, 2x, 2i5, 2i8, e20.7)
                end if
              else
                write (6, 410) i, nv2max
                write (9, 410) i, nv2max
410             format (' *** NUMBER OF NONZERO V2 ELEMENTS = ',i16,
     :            ' .GT. NV2MAX=',i16,'; ABORT ***')
                stop
              end if
            end if
350       continue
355     continue
        if (i .le. nv2max) lamnum(ilam) = inum
        if (bastst) then
          write (6, 370) ilam, lam, mu, lamnum(ilam)
          write (9, 370) ilam, lam, mu, lamnum(ilam)
370       format ('ILAM=',i4,'  L1=',i3,' M1=',i3,
     :      '  LAMNUM(ILAM) = ',i7)
        end if
        lamsum = lamsum + lamnum(ilam)
400   continue
      if (bastst) then
        write (6, 420) lamsum
        write (9, 420) lamsum,nv2max
420     format (' *** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS ',
     :      i16)
        write (6, 441) nv2max
        write (9, 441) nv2max
441     format('     MAXIMUM ALLOWED NUMBER IS',i16)
      end if
      return
      end
* --------------------------------------------------------------------
      subroutine vchirl(jp, lp, j, l, isp, is, jtot, lam, mu, vee,
     :        indp, ind)
*  subroutine to calculate v-lambda matrices for close-coupled
*  treatment of collisions of a chiral asymmetric top with an atom
*
*  the angular dependence is characterized by sums and differences
*  for Y(lambda,+mu) and Y(lambda,-mu) terms.  see notes by Claire
*  Rist and Paul Dagdigian.
*  matrix element of Y(lambda,mu) between symmetric top basis functions
*  given in eq. (26) of s. green, j. chem. phys. 64, 3463 (1976).
*
*  revised from vlmstp subr for asymmetric top levels by paul dagdigian
*
*  current revision date:  21-jan-2019 by pjd
*  -----------------------------------------------------------------------
*  variables in call list:
*    jp:       rotational quantum number of left side of matrix element (bra)
*    lp:       orbital angular momentum of left side of matrix element (bra)
*    j:        rotational quantum number of right side of matrix element (ket)
*    l:        orbital angular momentum of right side of matrix element (ket)
*    isp:      eigenvalue index of bra
*    is:       eigenvalue index of ket
*    jtot:     total angular momentum
*    lam:      value of lambda in expansion of potential
*    mu:       value of mu in expansion of potential
*    vee:      on return, contains desired matrix element
*    indp:     number of left-side level in channel list
*    ind:      number of right-side level in channel list
*  variable in common block /coatpi/
*    narray:   maximum size of asymmetric top basis fn expansion
*              set to 12 in himain (suitable for j <= 40)
*  variable in common block /coatpr/
*    c:        expansion coefficients for asymmetric top rotor wave fns.
*  subroutines called:
*    xf3j, xf6j, prmtp2
*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      common /coatpi/ narray
      common /coatpr/ c(1)
*
      data one, zero, two / 1.0d0, 0.0d0, 2.0d0 /
      data onsq4p /0.282094791773878d0 /
      vee = 0.d0
*   compute and check that 3j and 6j symbols are nonzero
      xlam = lam
      xlp = lp
      xl = l
      x3j = xf3j(xlp, xl, xlam, zero, zero, zero)
      if (x3j .eq. zero) goto 1000
      xjtot = jtot
      xjp = jp
      xj = j
      x6j = xf6j(xj, xl, xjtot, xlp, xjp, xlam)
      if (x6j .eq. zero) goto 1000
*  determine k_prolate and symmetry index of the asymmetric top wave functions
      kp = abs(isp)
      iepsp = sign(1, isp)
      k = abs(is)
      ieps = sign(1, is)
*  first determine size and symmetry block of wave function expansions
*  ket
      if (k .eq. 2*(k/2)) then
        if (ieps .eq. 1) then
*  E+ block
          isize = (j + 2)/2
          kmin = 0
        else
*  E- block
          isize = j/2
          kmin = 2
        end if
      else
*  O+, O- blocks
          isize = (j + 1)/2
          kmin = 1
      end if
*  bra
      if (kp .eq. 2*(kp/2)) then
        if (iepsp .eq. 1) then
*  E+ block
          isizp = (jp + 2)/2
          kminp = 0
        else
*  E- block
          isizp = jp/2
          kminp = 2
        end if
      else
*  O+,O- blocks
          isizp = (jp + 1)/2
          kminp = 1
      end if
*
      do 100 mm = 1,isize
      do 100 nn = 1,isizp
        isubp = (indp - 1)*narray + nn
        isub = (ind - 1)*narray + mm
        kbas = kmin + (mm - 1)*2
        kbasp = kminp + (nn - 1)*2
        call prmtp2 (jp, lp, j, l, jtot, kbasp, kbas,
     :    iepsp, ieps, lam, mu, vprm)
        xnorm = 0.5d0
        if (kbas .eq. 0) xnorm = xnorm / sqrt(two)
        if (kbasp .eq. 0) xnorm = xnorm / sqrt(two)
        iph = jp + j + kbasp - jtot
        vee = vee + (-one)**iph * vprm * xnorm * c(isubp) * c(isub)
100   continue
c
      xlam2 = two*xlam + one
      xjp2 = two*xjp + one
      xj2 = two*xj + one
      xlp2 = two*xlp + one
      xl2 = two*xl + one
      vee = vee * sqrt(xlam2 * xjp2 * xj2 * xlp2 * xl2) * x3j * x6j
     :  * onsq4p
c
1000  return
      end
* ----------------------------------------------------------------------
      subroutine prmtp2 (jp, lp, j, l, jtot, kp, k, iepsp, ieps,
     :  lam, mu, vprm)
*  subroutine to calculate primitive v-lambda matrix elements for close-coupled
*  treatment of collisions of a chiral symmetric top with an atom
*  the primitive cc matrix element are given in eqs. (26)
*  of s. green, j. chem. phys. 64, 3463 (1976)
*  note, that in this article the bra indices are unprimed and the ket indices
*  primed, while in the conventions of the present subroutine the bra indices
*  are primed and the ket indices, unprimed.
*
*  extensive modification of prmstp (in hibastp.f by m.h.alexander) and
*  prmatp (in hibaastp.f by p.dagdigian)
*
*  current revision date:  22-jan-2019 (p.dagdigian)
*  -----------------------------------------------------------------------
*  variables in call list:
*    jp:       rotational quantum number of left side of matrix element (bra)
*    lp:       orbital angular momentum of left side of matrix element (bra)
*    j:        rotational quantum number of right side of matrix element (ket)
*    l:        orbital angular momentum of right side of matrix element (ket)
*    jtot:     total angular momentum
*    kp:       k quantum number of bra
*    k:        k quantum number of ket
*    iepsp:    symmtry index of bra
*    ieps:     symmmefry index of ket
*    lam:      order of legendre term in expansion of potential
*    mu:       index of legendre term in expansion of potential
*    vrpm:     on return, contains primitive matrix element
*  subroutine called:
*     xf3j
*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      data one, two, zero, onsq2 / 1.d0, 2.d0, 0.d0,
     :  0.707106781186547d0 /
      vprm = zero
      xjp = jp
      xj = j
      xkp = kp
      xk = k
      xjtot = jtot
      xlp = float(lp)
      xl = float(l)
      xlam = float(lam)
      xmu = float(mu)
      if (mu .eq. 0) then
        iph = ieps * iepsp * (-1) ** (jp + j + lam)
        if (iph .eq. -1) return
        f3j1 = xf3j(xjp, xj, xlam, xkp, -xk, zero)
        f3j2 = xf3j(xjp, xj, xlam, -xkp, -xk, zero)
        vprm = two * (f3j1 + iepsp * f3j2)
      else
        if (mu .gt. 0) then
          iph = ieps * iepsp * (-1) ** (jp + j + lam + mu)
          if (iph .eq. -1) return
          f3j1 = xf3j(xjp, xj, xlam, xkp, -xk, xmu)
          f3j2 = xf3j(xjp, xj, xlam, -xkp, -xk, xmu)
          f3j3 = xf3j(xjp, xj, xlam, xkp, xk, xmu)
          f3j4 = xf3j(xjp, xj, xlam, -xkp, xk, xmu)
          vprm = onsq2 * two * (f3j1 + iepsp * f3j2
     :      + ieps * f3j3 + iepsp * ieps * f3j4)
        endif
      endif
c
      return
      end
* ----------------------------------------------------------------------
*      double precision function rotham(ji, ki, jf, kf)
*
*  subroutine to compute matrix elements of the asymmmetric top hamiltionian
*  in a prolate (case Ia) basis between unsymmetrized basis functions
*  (ji,ki) and (jf,kf)
*
*  THIS ROUTINE IS IN THE FILE hibaastp.f
*
*  author:  paul dagdigian
*  current revision date:  16-aug-2009
* ---------------------------------eof----------------------------------
