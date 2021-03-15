c   pot_ch3he_vib.f
c   authors: Qianli Ma
c
c   Pot routine for vibrational relaxation of CH3-He.
c
c   PES calculated by p.dagdigian and q.ma at RCCSD(T)/aug-cc-pVQZ level
c
c   Q. Ma, P. J. Dagdigian, and M. H. Alexander, J. Chem. Phys. 138, 104317 (2013).
c
*   Note:  this pot routine requires the following data files to be 
*   in hibxx/bin/progs/potdata:
c       pot_ch3he_vib_ylmsym
c       pot_ch3he_vib_ylmasym
c       pot_ch3he_vib_data
c   Included source file  (placed in hibxx/src/pot):
c       pot_ch3he_vib_common.f
c
c   These dummy subroutines are not used in this pot file
      include "common/ground"
c
c
c   Routine for compatibility with hib44
      subroutine datfln(filenm, fullnm)
      implicit none
      character (len=255) :: filenm, fullnm
      fullnm = "potdata/" // filenm
      return
      end subroutine datfln
c
c   The `regular' pot routine
c -------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h, o-z)
c
      include "pot_ch3he_vib_common.f"
c   Main function for `makepot', no arguments, print vvl's interactively
c
      double precision r, vv0
      integer i, j, k, iblock
c   Temporary for print
      double precision vvltmp(NVLM)
c
      print *, 'CH3-He vibrational relaxation'
   1  print *, 'Please input the maximum value of v2 (0--3):'
      read (5, *, end=99) vmax
      if (vmax .gt. 3 .or. vmax .lt. 0) goto 1
  10  print *, 'R (bohr), Ctrl+D to exit:'
      read (5, *, end=99) r
      call pot(vv0, r)
      iblock = 0
      do i = 0, vmax
        do j = i, vmax
          iblock = iblock + 1
          if (mod((i + j), 2) .eq. 0) then
            print 51, iblock, i, j
          else
            print 52, iblock, i, j
          endif
          call dcopy(NVLM, vvl((iblock - 1) * NVLM + 1), 1, vvltmp, 1)
          print 55, vvltmp * ECONV
        enddo
      enddo
      goto 10
  51  format (1x, 'vibrational block #', i2, ', from ', i1, ' to ',
     +        i1, ', (lambda + mu) is even')
  52  format (1x, 'vibrational block #', i2, ', from ', i1, ' to ',
     +        i1, ', (lambda + mu) is odd')
  55  format(' v(lam,0):', 5(1pe13.5), /, ' v(lam,3):', 4(1pe13.5), /,
     +       ' v(lam,6):', 2(1pe13.5), /, ' v(lam,9):', 1(1pe13.5))
  99  end
c   end subroutine driver
c
c
c -------------------------------------------------------------------
      subroutine loapot(iunit, filnam)
      implicit double precision (a-h, o-z)
c
      include "pot_ch3he_vib_common.f"
      include "common/parpot"
c   Initialize parameters for the potential
c
c   Arguments:
c       Arguments are not refered to in this basis routine
      character*(*) filnam
      integer iunit
c
c
c   Hidden returned value:
c       /conlam/ nlam
c       /cosysi/ ipotsy
c       /coptnm/ potnam
c
c
      potnam = 'CH3-He vibrational relaxation'
c   vvl now represent <v_2'|v_{\lambda\mu}(Q_2, R)|v_2> in the
c   following order:
c       <0|v_00|0>, <0|v_20|0>, ..., <0|v_99|0>, <0|v_10|1>, ...
c   where v2 >= v2'
c
c   vv0 will not be used (set to zero).
c
c   The dimension of vvl array:
      nlam = NVLM * (vmax + 1) * (vmax + 2) / 2
c
c   Set symmetry (forced)
      ipotsy = 3
c
      return
      end
c   end subroutine laopot
c
c
c -------------------------------------------------------------------
      subroutine pot(vv0, r)
      implicit double precision (a-h, o-z)
c   Subroutine to calculate <v_2'|v_{\lambda\mu}(Q_2, R)|v_2> in the
c   following order:
c       <0|v_00|0>, <0|v_20|0>, ..., <0|v_99|0>, <0|v_10|1>, ...
c   where v2 >= v2'
c
c   vv0 will not be used (set to zero).
c
      include "pot_ch3he_vib_common.f"
c
c
c   Arguments:
c       r: intermolecular distance
      double precision r
c
c
c   Returned value:
c       vv0: zero
      double precision vv0
c
c
c   Hidden returned value:
c       /covvl/ vvl
c
c
c   Function called:
      integer gblkid
c
c
c   vsp: potential for all theta/phi tuples, to be obtained from splch3
      double precision vsp(NANGLE)
c   ivi, ivf, iblock: indeces in loop
      integer ivi, ivf, iblock
c   The following are variables required by dgelsd
      double precision RCOND
      parameter (RCOND=1d-6)
      integer LWORK
      parameter (LWORK=NANGLE*NVLM)
      integer info, irank, iwork(LWORK/3)
      double precision swork(NVLM), work(LWORK)
c
c   The following are used in treating long range potential
c     RC: at which R the switching function is centered
c     C6: long term C6 parameter (V = C6/R^6, no angular dependence)
c     C6 is fitted from ab initio potential at R = 20 Bohrs
      double precision S, RC, C6
      parameter (S=1d0, RC=1.5d1, C6=8.752841d6)
      double precision stepfc, lrpot
      integer nvvlp, i, j
c
c   Y_lm terms, ylm coefficients for calculation (will be modified in
c   dgelsd), and YLMC(S/A) is a stationary copy of the array.
c   s/a stands for symmetric/anti-symmetric (of <v2'|V|v2>)
      double precision ylm(NANGLE, NVLM)
      double precision YLMCS(NANGLE, NVLM), YLMCA(NANGLE, NVLM)
c
c   Load Y_lm terms in the first call
      logical isfst
      data isfst /.true./
      character*255 datfl
      if (isfst) then
        call datfln('pot_ch3he_vib_ylmsym', datfl)
        open (unit=10, file=datfl)
        read (10, *) YLMCS
        close(10)
        call datfln('pot_ch3he_vib_ylmasym', datfl)
        open (unit=10, file=datfl)
        read (10, *) YLMCA
        close(10)
        isfst = .false.
      endif
c   Calculate for each (v2, v2') block
      iblock = 0
      do ivi = 0, vmax
        do ivf = ivi, vmax
          iblock = iblock + 1
c   Get potential values for the current block at distance R
          call splch3(vsp, r, ivi, ivf)
c   Check the symmetry of the (v2, v2') block and obtain ylm coefficients
          if (mod((ivi + ivf), 2) .eq. 0) then
            ylm = YLMCS
          else
            ylm = YLMCA
          endif
c   Linear least-square fit
          call dgelsd(NANGLE, NVLM, 1, ylm, NANGLE, vsp, NANGLE,
     +                swork, RCOND, irank, work, LWORK, iwork, info)
          call dcopy(NVLM, vsp, 1, vvl((iblock - 1) * NVLM + 1), 1)
        enddo
      enddo
c   Coefficients for long-range potential
      stepfc = 1d0 - 5d-1 * (dtanh(S * (r - RC)) + 1d0)
      lrpot = (stepfc - 1d0) * C6 / r ** 6
c   Total number of vvl terms for the given vmax
      nvvlp = (vmax + 1) * (vmax + 2) / 2 * NVLM
c   Damp the potential at long range
      call dscal(nvvlp, stepfc, vvl, 1)
      do i = 0, vmax
        iblock = gblkid(i, i, vmax)
        j = (iblock - 1) * NVLM + 1
        vvl(j) = vvl(j) + lrpot
      enddo
c   Convert potential to hartree
      call dscal(nvvlp, 1d0 / ECONV, vvl, 1)
c   vv0 is not used here
      vv0 = 0
      return
      end
c   end subroutine pot
c
c
c -------------------------------------------------------------------
      subroutine splch3(vsp, r, v2_i, v2_f)
      implicit double precision (a-h, o-z)
      include "pot_ch3he_vib_common.f"
c   Spline the interaction potential (integreted over vibrational
c   coordinate Q2) to the given R for a particular (v2, v2') tuple.
c
c   Arguments:
c       r: intermolecular distance
c       v2_i, v2_f: initial and final vibrational level
      double precision r
      integer v2_i, v2_f
c
c
c   Returned value:
c       vsp: <v2'|V|v2> at R for all (theta, phi) tuples
      double precision vsp(*)
c
c
c   Function called:
      integer gblkid
c
c
c   Potential as calculated
      double precision v(NDIST, NANGLE, V2TMAX)
c   Parameters from linear fit
      double precision b(NDIST, NANGLE, V2TMAX)
      double precision c(NDIST, NANGLE, V2TMAX)
      double precision d(NDIST, NANGLE, V2TMAX)
c   Loop indeces
      integer i, k
c   Block index
      integer iblock
c   Intermolecular distances used in the ab initio calculation
      double precision RR(NDIST)
      data RR /3.5d0, 4d0, 4.5d0, 5d0, 5.5d0, 6d0, 6.5d0, 7d0, 7.5d0,
     +         8d0, 8.5d0, 9d0, 9.5d0, 1d1, 1.1d1, 1.2d1, 1.3d1,
     +         1.5d1, 2d1/
c   Function to evaluate splined potential
      double precision seval
c
      character*255 datfl
      logical isfst
      data isfst /.true./
c   In the first call, read data and determine fitted coefficients
      if (isfst) then
c   Read data file
        call datfln('pot_ch3he_vib_data', datfl)
        open (unit=10, file=datfl)
        read (10, *) v
        close(10)
c   Calculate coefficients for each block
        do k = 1, V2TMAX
          do i = 1, NANGLE
c   calculate coefficients
            call spline(NDIST, RR, v(1, i, k), b(1, i, k), c(1, i, k),
     +                  d(1, i, k))
          enddo
        enddo
        isfst = .false.
      endif
c
c   On a regular run, calculate potential from previously determined
c   spline coefficients.
c
c   Calculate potential from spline coefficients
      iblock = gblkid(v2_i, v2_f, V2MAX)
      do i = 1, NANGLE
        vsp(i) = seval(NDIST, r, RR, v(1, i, iblock), b(1, i, iblock),
     +                 c(1, i, iblock), d(1, i, iblock))
      enddo
      return
      end
c   end subroutine splch3
c
c
c
c   User defined basis
c -------------------------------------------------------------------
      subroutine syusr(irpot, readpt, iread)
      implicit double precision (a-h, o-z)
c
      include "pot_ch3he_vib_common.f"
c   Subroutine to read parameters for CH3 v2 vibrational relaxation
c
c   Parameters:
c       iread: 1 to read data from input file, 0 to set default
c       irpot, readpt: not refered to in this basis
      integer irpot, iread
      logical readpt
c
c
c   Hidden returned value:
c       /cosys/ scod
c       /cosysi/ nscode, isicod, ipotsy, iop, jmax, vmax
c       /cosysr/ isrcod, emax0, emax1, emax2, emax3
c       /coiout/ niout, indout
c
c
c   Constants:
c       CNIOUT, CIOUT: defaults to indout
      integer CNIOUT, CIOUT(18)
      parameter (CNIOUT=18)
      data CIOUT /0, 1, -1, 2, -2, 3, -3, 4, -4,
     +            100, 101, -101, 102, -102, 103, -103, 104, -104/
c
c   Loop variable:
      integer i
c
      character*(*) fname
c
c
c   Set system dependent parameters
      scod(1) = 'NTERM'
      scod(2) = 'IPOTSY'
      scod(3) = 'IOP'
      scod(4) = 'JMAX'
      scod(5) = 'VMAX'
      scod(6) = 'EMAX0'
      scod(7) = 'EMAX1'
      scod(8) = 'EMAX2'
      scod(9) = 'EMAX3'
      scod(10) = 'LAMMIN'
      scod(11) = 'LAMMAX'
      scod(12) = 'MPROJ'
c   Set lengths of scod array
      nscode = LENCOD
      isicod = ICOD
      isrcod = IRCOD
c   Set default indout data
c     This should not be done if an input file is read!
      if (niout .le. 0) then
         niout = CNIOUT
         do i = 1, CNIOUT
            indout(i) = CIOUT(i)
         enddo
      endif
c   Read data if required
      if (iread .eq. 1) then
        read (8, *, err=80) ipotsy, iop
        read (8, *, err=80) jmax, vmax
        read (8, *, err=80) emax0, emax1, emax2, emax3
        call loapot(1, ' ')
        close(8)
      endif
      return
c   On read error
   80 write(6,90)
   90 format(/'   *** ERROR DURING READ FROM INPUT FILE ***')
      return
c
c -------------------------------------------------------------------
      entry ptrusr(fname, readpt)
c   This subroutine will not be used.
      return
c
c -------------------------------------------------------------------
      entry savusr(readpt)
c   Save parameters
      write (8, 201) ipotsy, iop
      write (8, 202) jmax, vmax
      write (8, 203) emax0, emax1, emax2, emax3
  201 format (i4, 4x, i4, 21x, 'ipotsy, iop')
  202 format (i4, 4x, i4, 21x, 'jmax, vmax')
  203 format (4(f7.2, 1x), 1x, 'emax0, emax1, emax2, emax3')
      return
      end
c
c
c
c
c   The `regular' basis routine
c -------------------------------------------------------------------
      subroutine bausr(j, l, is, jhold, ehold, ishold, nlevel,
     +                  nlevop, k, ieps, jtemp, ktemp, rcut, jtot,
     +                  flaghf, flagsu, csflag, clist, bastst,
     +                  ihomo, nu, numin, jlpar, n, nmax, ntop)
      implicit double precision (a-h, o-z)
c
      include "pot_ch3he_vib_common.f"
c
c   Arguments (13):
c       rcut: cut-off point for keeping higher energy channels
      double precision rcut
c       jtot: total angular momentum
c       nu: not refered to
c       numin: not refered to
c       jlpar: total parity of included channels in cc calculation
c       nmax: maximum number of channels
      integer jtot, nu, numin, jlpar, nmax
c       flaghf: should be true as the system has integer spin
c       flagsu: should be false for atom-molecule collision
c       csflag: should be false for CC calculations
c       clist: if true,quantum numbers and energies listed for each channel (unable to set in hibridon?)
c       bastst: if true, execution terminates after the first call to basis and show channel list
c       ihomo: not used
      logical flaghf, flagsu, csflag, clist, bastst, ihomo
c
c
c   Returned value (11):
c       j: rotational quantum numbers for each channel
c       k: projection quantum number for each channel
c       l: orbital angular momentum for each channel
c       is: index [eps * (100 * v2 + k)] for each channel
c       ieps: epsilon label for each channel
      integer j(*), k(*), l(*), is(*), ieps(*)
c       jhold: rotational quantum numbers for each rotational level
c       ishold: symmetry index of each rotational level
      integer jhold(*), ishold(*)
c       ehold: energy in hartrees of each rotational level
      double precision ehold(*)
c       nlevel: number of rotational levels used in channel basis
c       nlevop: number of rotational levels used in channel basis which are open asymptotically
c       n: number of channels
c       ntop: n or n+1 (must be odd)
      integer nlevel, nlevop, n, ntop
c
c
c   Workspace as arguments (2):
c       jtemp, ktemp
      integer jtemp(*), ktemp(*)
c
c
c   Hidden arguments:
c       /conlam/ nlam, nlammx
c       /ch3he/ brot, crot, evib, lamsym, lamasy, musym, muasy
c       /cosysi/ nscode, isicod, nterm, ipotsy, iop, jmax, vmax
c       /cosysr/ isrcod, emax0, emax1, emax2, emax3
c       /coered/ ered, rmu
c       /cov2/ nv2max
c       /coipar/ iprint
c
c
c   Hidden returned value:
c       /conlam/ lamnum
c       /cov2/ v2
c       /coiv2/ iv2
c       /cocent/ cent
c       /coeint/ eint
c
c
c   Variables used only in the subroutine:
c
c   Epsilon of double precision number, used to check if a
c   number is zero
      double precision EPS
      parameter (EPS=1.12d-16)
c
      integer nlist
      integer i, vi, ki, ji, iep, i1, j1, njk, nn
      integer ilmmin, ilmmax, vibblk, ij
      double precision lvleng
      integer vtemp(KMAX), ietemp(KMAX)
      double precision esave, etemp
      integer vsave, jsave, ksave, iesave
      integer lamsum
      double precision vee
      double precision emax(V2MAX+1), emin
      integer ipar, lpar, lmax, lmin, li, ilm, ilms, lambda, mu
      integer v(KMAX)
c
c
c   Function called:
      integer gblkid
c
c -------------------------------------------------------------------
c
c   Construct emax array from input parameters
      emax(1) = emax0
      emax(2) = emax1
      emax(3) = emax2
      emax(4) = emax3
c
c   Pre-run check
c
c   CS calculation not implemented yet
      if (csflag) then
        write (6, 5)
        write (9, 5)
    5   format (' *** CS CALCULATION NOT IMPLEMENTED; ABORT ***')
        stop
      endif
c   Electron spin not considered here
      if (flaghf) then
        write (6, 6)
        write (9, 6)
    6   format (' *** FLAGHF = .TRUE. FOR SINGLET SYSTEM; ABORT ***')
        stop
      endif
c   Check if numbers of V_lm terms exceeds limit
      if (nlammx .lt. nlam) then
        write (6, 7) nlam, nlammx
        write (9, 7) nlam, nlammx
    7   format (' *** NLAMMX NOT BIG ENOUGH; ABORT ***', i3, i3)
        stop
      endif
c
c   Write info when clist is true
      if (clist) then
        write (6, 10) rmu*XMCONV, ered*econv, jtot, jlpar
        write (9, 10) rmu*XMCONV, ered*econv, jtot, jlpar
   10   format (/, ' **  CC SYMMETRIC TOP VIBRATION **',
     +          /, '    RMU=', f9.4, '  E=', f7.2, '  JTOT=', i4, 2x,
     +          '  JLPAR=', i2)
      endif
c
c
c   Set up ro-vibrational basis
c
c   Set up the lisk of levels
      nlist = 0
      do 15 vi = 0, vmax
        do 15 ki = 0, jmax
          do 15 ji = ki, jmax
            if (abs(iop) .eq. 1) then
              iep = (-1) ** (vi + ji)
            else
              iep = - (-1) ** (vi + ji)
            endif
            if (iop .le. 0) then
              if (mod(ki, ipotsy) .ne. 0) goto 15
              if (ki .eq. 0) then
                if (mod(vi + ji, 2) .ne. 0) goto 15
                if (iep .eq. -1) goto 15
              endif
            else
              if (mod(ki, ipotsy) .eq. 0) goto 15
            endif
c   Calculate level energy and make sure it does not exceeds emax
            lvleng = evib(vi+1) + brot(vi+1) * ji * (ji+1)
     +               + (crot(vi+1) - brot(vi+1)) * ki ** 2
            if (lvleng .gt. emax(vi+1)) goto 15
c   Levels survived here should be included in the channel basis
            nlist = nlist + 1
            ehold(nlist) = lvleng / econv
            vtemp(nlist) = vi
            jtemp(nlist) = ji
            ktemp(nlist) = ki
            ietemp(nlist) = iep
   15 continue
c   Check if any channel exist
      if (nlist .eq. 0) then
        write (6, 20)
        write (9, 20)
   20   format (' *** NO LEVEL EXISTS; ABORT ***')
        stop
      endif
c   Sort the list of levels in increasing energy - bubble sort
      do i1 = 1, nlist - 1
        do j1 = nlist, i1 + 1, -1
          if (ehold(j1) .lt. ehold(i1)) then
            call iswap(vtemp(i1), vtemp(j1))
            call iswap(jtemp(i1), jtemp(j1))
            call iswap(ktemp(i1), ktemp(j1))
            call iswap(ietemp(i1), ietemp(j1))
            etemp = ehold(i1)
            ehold(i1) = ehold(j1)
            ehold(j1) = etemp
          endif
        enddo
      enddo
c     Determint the number of open levels
      nlevop = 0
      do while (nlevop .lt. nlist .and. ehold(nlevop + 1) .le. ered)
         nlevop = nlevop + 1
      end do
c
c   Print list for bastst run
      if (bastst) then
        print 130
        print 131
  130   format (/, 9('-'), ' SORTED LEVEL LIST ', 9('-'))
  131   format (3x, 'N', 4x, 'V', 4x, 'J', 4x, 'K', 2x, 'EPS',
     +          3x, 'EINT(CM-1)')
        do i1 = 1, nlist
          print 135, i1, vtemp(i1), jtemp(i1), ktemp(i1), ietemp(i1),
     +                   ehold(i1)*econv
        enddo
  135   format (4(i4, 1x), 1x, i2, 2x, f12.4)
        print *
      endif
c
c
c   Set up channel and level list for scattering calculation
      n = 0
      nlevel = 0
      do njk = 1, nlist
        vi = vtemp(njk)
        ki = ktemp(njk)
        ji = jtemp(njk)
        iep = ietemp(njk)
        lvleng = ehold(njk)
        nlevel = nlevel + 1
        jhold(nlevel) = ji
        ishold(nlevel) = (100 * vi + ki) * iep
c   JLPAR check
        ipar = iep * (-1) ** (ji + ki)
        lmax = jtot + ji
        lmin = iabs(jtot - ji)
        do li = lmin, lmax
          lpar = (-1) ** (li - jtot)
          if (ipar * lpar .eq. jlpar) then
            n = n + 1
            if (n .gt. nmax) then
              write (6, 150)
              write (9, 150)
  150         format (' *** # CHANNELS EXCEEDS LIMIT; ABORT ***')
              stop
            endif
            is(n) = ishold(nlevel)
            eint(n) = lvleng
            v(n) = vi
            k(n) = ki
            ieps(n) = iep
            j(n) = ji
            l(n) = li
            cent(n) = li * (li + 1)
          endif
        enddo
      enddo
c
      if (bastst) then
        print 140, nlevel, iop
  140   format (' ** ', i3, ' LEVELS FOR IOP=', i2)
        print *
      endif
c
c   rcut mechanism
      if (rcut .gt. 0d0) then
c   First determine the lowest channel energy for which:
c       1. open asymptotically
c       2. closed at r = rcut
        emin = 1d7
        do i = 1, n
          if (eint(i) .le. ered) then
            if (jtot * (jtot + 1) / (2d0 * rmu * rcut ** 2)
     +          .gt. (ered - eint(i))) then
              if (eint(i) .lt. emin) emin = eint(i)
            endif
          endif
        enddo
c   Keep channels whose energy is less than emin only
        if (emin .lt. ered) then
          nn = 0
          do i = 1, n
            if (eint(i) .lt. emin) then
              nn = nn + 1
              eint(nn) = eint(i)
              v(nn) = v(n)
              j(nn) = j(i)
              ieps(nn) = ieps(i)
              is(nn) = is(i)
              cent(nn) = cent(i)
              k(nn) = k(i)
              l(nn) = l(i)
            endif
          enddo
          n = nn
        endif
      endif
c
c
c   Return if no channel
      if (n .eq. 0) return
c
c   Setting ntop - mechnism not fully understood
      if (nu .eq. numin) then
        ntop = max(n, nlevop)
c   ntop is the maximum row dimension of all matrices passed in the
c   call list of subroutines propag and soutpt.
c   for fps make sure this is an odd number, for faster bank access.
c   this has no effect on vax or cray
        if (mod(ntop, 2) .eq. 0) ntop = ntop + 1
      else
        if (n .gt. ntop) then
          write (6, 303) nu, n, ntop
          write (9, 303) nu, n, ntop
          call exit
        endif
      endif
  303 format (' *** NCH = ', i3, ' AT NU = ', i2, ' .GT. NTOP = ', i3,
     +        '; ABORT **', /, '     CHECK RCUT')
c
c
c   List channels if requested (iprint >= 1)
      if (bastst .and. iprint .ge. 1) then
        print 310
        print 311
  310   format (19('-'), ' CHANNEL LIST ', 18('-'))
  311   format (3x, 'N', 4x, 'V', 4x, 'J', 4x, 'K', 3x, 'EPS', 4x,
     +          'L', 4x, 'IND', 4x, 'EINT(CM-1)')
        do i = 1, n
          print 315, i, v(i), j(i), k(i), ieps(i), l(i), is(i),
     +               eint(i)*ECONV
        enddo
        print *
  315   format (4(i4, 1x), 2x, i2, 2x, i4, 1x, i6, 2x, f12.4)
        print 320, nlevop, n, jtot, jlpar
  320   format (' ** ', i4, '/', i4, ' CHANNELS FOR JTOT=', i3,
     +          ', JLPAR=', i2)
        print *
      endif
c
c   Calculate coupling matrix
c
c   Reset the lamnum array
      do i = 1, nlam
        lamnum(i) = 0
      enddo
      lamsum = 0
c   v2 matrix is to be expanded as
c       (lam1, ij1), (lam1, ij2), ..., (lam1, ijn), (lam2, ij1), ...
c   non-zero elements are not stored at all
      i = 0
      do 160 ilm = 1, nlam
        do 160 i1 = 1, n
          do 160 j1 = i1, n
            vibblk = gblkid(v(i1), v(j1), vmax)
            ilmmin = NVLM * (vibblk - 1) + 1
            ilmmax = NVLM * vibblk
            if (ilm .gt. ilmmax .or. ilm .lt. ilmmin) goto 160
            ilms = ilm - ilmmin + 1
            ij = ntop * (i1 - 1) + j1
c   Calculate and write non-zero v2 elements
            if (mod(v(i1) + v(j1), 2) .eq. 0) then
c   Symmetric vibrational coupling potential
              lambda = lamsym(ilms)
              mu = musym(ilms)
            else
c   Anti-symmetric vibrational coupling potential
              lambda = lamasy(ilms)
              mu = muasy(ilms)
            endif
            call vlmstp(j(i1), l(i1), j(j1), l(j1), jtot,
     +                  k(i1), k(j1), lambda, mu,
     +                  ieps(i1), ieps(j1), vee, .false.)
            if (dabs(vee) .gt. EPS) then
              i = i + 1
c   Protect from segmentation fault
              if (i .ge. nv2max) goto 450
c   Save the non-zero v2 element
              v2(i) = vee
              iv2(i) = ij
              lamnum(ilm) = lamnum(ilm) + 1
              lamsum = lamsum + 1
c   Print non-zero v2 element if requested
                if (bastst .and. iprint .ge. 2) then
                  print 431, lambda, mu, i1, j1, vee
                endif
  431   format (' V2 FOR LAMBDA=', i2, ' MU=', i2,
     +          ', FOR CHANNELS ', i4, ' AND', i4,
     +          ' IS ', f13.6)
            endif
  160 continue
c
c
      if (bastst) then
        write (6, 430) lamsum
        write (9, 430) lamsum
  430   format (' ** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS',
     +          1x, i8)
c   Print lamnum array (non-zero elements only)
        if (iprint .ge. 1) then
          do i = 1, nlam
            if (lamnum(i) .gt. 0) then
              ilm = mod((i - 1) / 12, 2)
              if (ilm .eq. 0) then
                lambda = lamsym(mod(i - 1, 12) + 1)
                mu = musym(mod(i - 1, 12) + 1)
              else
                lambda = lamasy(mod(i - 1, 12) + 1)
                mu = muasy(mod(i - 1, 12) + 1)
              endif
              print 170, i, lambda, mu, lamnum(i)
            endif
          enddo
          print *
        endif
  170   format (' ILAM=', i3, ' LAM=', i3, ' MU=', i3,
     +          ' LAMNUM(ILAM) = ', i6)
      endif
c
      return
c
c   Error: the number of non-zero v2 elements exceeds limit
  450 write (6, 451)
      write (9, 451)
  451 format (' *** TO MANY NON-ZERO V2 ELEMENTS; ABORT ***')
      stop
c
      end
c   end subroutine bastpv
c
c
c -------------------------------------------------------------------
      function gblkid(v2, v2p, vmax)
      implicit none
      integer gblkid, v2, v2p, vmax
c
c   Subroutine to determine the ID of a (v2, v2') block in the
c   compact form - 00, 01, 02, ..., 11, 12, ...
c
c   Arguments
c       v2, v2p: v_2 and v_2'
c       vmax: maximum value of v2
c   Return
c       The corresponding block ID for the v2, v2p combination
c
c   Example
c       gblkid(1, 0, 4) = 2
c           when vmax = 4, (10)=(01) is the second block in the list
c       gblkid(1, 2, 3) = 6
c           when vmax = 3, (12)=(21) is the sixth block in the list
c
      integer vl, vg
c   The lesser/greater v2p
      if (v2 .lt. v2p) then
        vl = v2
        vg = v2p
      else
        vl = v2p
        vg = v2
      endif
c   In the compact form, the lesser v_2 (v_2l) is indexed first.
c   The number of blocks for v_2l < vl is
c       (vmax + 1 - 0) + (vmax + 1 - 1) + ... + (vmax + 1 - (vl - 1))
c   For v_2l = vl, v_2g = vg is the
c       vg - vl + 1
c   th block.
      gblkid = ((2 * vmax + 3 - vl) * vl) / 2
     +         + (vg - vl + 1)
      return
      end
c   end function gblkid
