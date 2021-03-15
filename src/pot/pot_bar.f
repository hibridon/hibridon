*  BAr(X2Pi) potential
*  E. Hwang, Y.-L. Huang, P.J.Dagdigian, and M.H.Alexander, 
*  J. Chem. Phys. 98, 8484 (1993).
*
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parpot"
      include "common/parbas"
      common /coselb/ ibasty
      potnam='B-AR X-B'
      ibasty=99
      lammin(1)=1
      lammax(1)=1
      mproj(1)=0
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* --------------------------------------------------------------------------

      subroutine pot (vv0, r)
*  -----------------------------------------------------------------------

*  subroutine to calculate the r-dependent coefficients in the
*  model B-Ar B state potential curve

*  on return:
*  vv0 contains the isotropic term (n=0) in the potential
*  vvl(1) contains the (1,1) potential (this is zero here)

*  variable in common block /conlam/
*    nlam:      the total number of angular coupling terms
*    nlammx:    the maximum number of anisotropic terms
*  variable in common block /covvl/
*    vvl        array to store r-dependence of each angular term in the
*               potential

      implicit double precision (a-h, o-z)
      common /conlam/ nlam, nlammx

      common /covvl/ vvl(1)

      dimension rpaul(54),vpaul(54),tt(60),bb(60),work(350)
      data rpaul/3.4770942,3.4959914,3.5054400,3.5148887,
     : 3.5243373,3.5337859,3.5432345,3.5526832,3.5621318,
     : 3.5715804,3.5810290,3.5904777,3.5999263,3.6093749,
     : 3.6131544,3.6207133,3.6301619,3.6433900,3.6660667,
     : 3.7000817,3.7567735,3.8682673,4.0855856,4.3955006,
     : 4.6997463,4.9643078,5.2213104,5.4877617,5.7500555,
     : 5.7712204,6.0000661,6.0792456,6.4231756,6.5000874,
     : 7,7.5,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25/
      data vpaul/
     : 2.3081000e+3,1.9208000e+3,1.7372000e+3,1.5598000e+3,
     : 1.3886000e+3,1.2232000e+3,1.0635000e+3,9.0917000e+2,7.6012000e+2,
     : 6.1612000e+2,4.7699000e+2,3.4255000e+2,2.1264000e+2,8.7076000e+1,
     : 3.8000000e+1,-5.8100000e+1,-1.7020000e+2,-2.9850000e+2,
     : -4.4280000e+2,
     : -6.0320000e+2,-7.7970000e+2,-9.7220000e+2,-1.0745000e+3,
     : -9.7220000e+2,
     : -7.7970000e+2,-6.0320000e+2,-4.4280000e+2,-2.9850000e+2,
     : -1.7910000e+2,
     : -1.7020000e+2,-8.3453000e+1,-5.8100000e+1,3.8000000e+1,
     : 5.4525000e+1,
     : 1.2345000e+2,1.4047000e+2,1.2631000e+2,6.7020000e+1,1.7780000e+1,
     : -7.6048000,-1.6132000e+1,-1.6122000e+1,-1.3110000e+1,-9.7064000,
     : -6.8738000,-4.9409000,-3.6036000,-2.6646000,-1.9960000,
     : -1.5134000,-1.1605000,-8.9933000e-1,-7.0379000e-1,-5.5582000e-1/
      data interp / 1/
      data izero /0/
      data zero,half,one,two /0.d0,0.5d0,1.d0,2.d0/
      data inb  /1/
*
      if (interp .eq. 1) then
*
*  here for first pass through pot subroutine
*  define spline initial conditions
      ibcl = 1
      ibcr = 2
      fbcl = zero
      fbcr = zero
      kntopt = 2
      inbv=1
      dr=rpaul(2)-rpaul(1)
*  compute initial derivatives
      dd=(vpaul(2)-vpaul(1))/dr


*  interpolate the  coefficients for the energies
          call bint4 (rpaul, vpaul, 54, ibcl,ibcr,dd,fbcr,kntopt,
     :                tt, bb,nc,kord, work)
        interp = 0
      endif
* here is  the entry point for subsequent entries into pot subroutine
* here to determine potential
      kord=4
      izero=0
      vvpaul= bvalu (tt, bb, nc, kord,izero, r, inb, work)
      vv0=vvpaul/219474.6
      vvl(1)=0.d0
      return
      end
c  -----------------------------------------------------------------
      subroutine ground(wf, r, nch, nphoto, mxphot)
c  -----------------------------------------------------------------
*  driven equation source vector for BAr
*  author:  millard alexander
*  current revision date:  10-dec-1992
c  -----------------------------------------------------------------
*     variables in call list:

*     wf        array of dimension nch*nphoto, containing, on return,
*               ground state wavefunction in each of nch components
*               nphoto is number of difference ground state wavefunctions
*     r         value of separation coordinate
*     nch       total number of channels (row dimension of q)
*     nphoto    number of different wavefunctions calculated
*               column index of q vector
*     mxphot    maximum size of q vector (mxphot .ge. nch*nphoto)
*  -------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (nl =  25)
      dimension wf(6), sc(6)
      common /cotrans/ t(6,6)
      common /cojtot/ j,jlpar
      common /coered/ ered, rmu
      data two /2.d0/
      dimension rl(nl), dip1(nl), ddip1(nl), dip3(nl)
* these are the morse parameters for bar
* note that g0 is more commonly denoted beta
      data de, r0, g0 /0.000463379, 7.01094, 0.75517/
* determine morse function for ground vibrational state
      rmlmda = two*rmu
      d = dsqrt(rmlmda*de)/g0
      twod = two*d
      twodm1 = twod-one
      dln2d = log(twod)
      gln2d = gammln(twodm1)
      z = g0*(r-r0)
      arg = twodm1*(dln2d-z)-gln2d-twod*exp(-z)
      p = sqrt(g0*dexp(arg))
* determine transition moment
        dipfunc=5.5d0
* initialize source vector to zero
* since 1Pi state is channel 6, projection of ground state is only onto
* this channel
      wf(1)=rmlmda*p*dipfunc
      entry wfintern(wf,yymin,nnvib,nny)
      return
      end
* -----------------------------------------
      function gammln (xx)
      implicit double precision (a-h,o-z)
c
c     -----------------------------------------------------------------
c     natural log of the gamma function (numerical recipes)
c     -----------------------------------------------------------------
c
      dimension cof(6)
      data cof,stp /76.18009173d0,-86.50532033d0,24.01409822d0,
     + -1.231739516d0,0.120858003d-2,-0.536382d-5,2.50662827465d0/
c
      x = xx-1.0d0
      tmp = x+5.5d0
      tmp = (x+0.5d0)*dlog(tmp)-tmp
      ser = 1.0d0
      do 11 j = 1,6
         x = x+1.0d0
         ser = ser+cof(j)/x
  11  continue
      gammln = tmp+dlog(stp*ser)
      return
      end
* -----------------------------------------
      subroutine syusr (irpot, readp, iread)
*  dummy syusr subroutine
      logical readp
      character*(*)fname
      if(.not.readp.or.iread.eq.0) then
        call loapot(1,' ')
        return
      endif
      entry ptrusr (fname)
      entry savusr (readp)
      return
      end
*  -----------------------------------------------------------------------
      subroutine bausr (j, l, is, jhold, ehold, ishold, nlevel, nlevop,
     :                  sc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu,
     :                  csflag, clist, bastst, ihomo, nu, numin, jlpar,
     :                  n, nmax, ntop)
* --------------------------------------------------------------------
*  subroutine to determine angular coupling potential
*  for B-AR B state
*  authors:  millard alexander
*  current revision date:  10-dec-92
* --------------------------------------------------------------------
*  variables in call list:
*    j:        on return contains rotational quantum numbers for each
*              channel (dummy)
*    l:        on return contains orbital angular momentum for each
*              channel (zero)
*    is:       on return contains symmetry index for each channel (dummy)
*    jhold:    on return contains rotational quantum numbers for each
*              rotational level (dummy)
*    ehold:    on return contains energy in hartrees of each level
*              
*    ishold:   on return contains symmetry index of each rotational level
*    nlevel:   on return contains number of energetically distinct
*              rotational levels used in channel basis
*    nlevop:   on return contains number of energetically distinct
*              rotational levels used in channel basis which are open
*              asymptotically
*    sc1,sc2:  scratch vectors of length at least nmax
*    sc3,sc4:  scratch vectors of length at least nmax
*              these scratch vectors are not used here
*    rcut:     cut-off point for keeping higher energy channels
*              if any open channel is still closed at r=rcut, then
*              all closed channels as well any open channels which are
*              still closed at r=rcut are dropped from basis
*              note that this is ignored in molecule-surface collisions!!!
*    jtot:     total angular momentum
*    flaghf:   if .true., then system with half-integer spin
*              if .false., then system with integer spin (this is the case
*              here)
*    flagsu:   if .true., then molecule-surface collisons
*    csflag:   if .true., then coupled-states calculation
*              if .false., the close-coupled calculation
*    clist:    if .true., then quantum numbers and energies listed for
*              each channel
*    bastst:   if .true., then execution terminates after the first call
*              to basis
*    ihomo:    if .true. , then homonuclear molecule
*              if .false., then heteronuclear molecule
*              if the molecule is homonuclear (ihomo = .true.), the
*              rotational levels included go from jmin to jmax in steps
*              of 2 and only even lambda terms in the anisotropy are
*              included
*    nu:       coupled-states projection index
*    numin:    minimum coupled states projection index
*              for cc calculations nu and numin are both set = 0 by calling
*              program
*    jlpar:    total parity of included channels in cc calculation
*              only those channels are included for which
*                  (-1)**(j+l-jtot)=jlpar
*    n:        on return equals number of channels
*    nmax:     maximum dimension of arrays
*    ntop:     maximum row dimension of all matrices passed to subroutines
*              propag and soutpt.  ntop is set in basis only if nu = numin
*              otherwise it is unchanged from the value supplied by the
*              calling program
*    note!!!   if flaghf = .true., then the true values of the rotational
*    quantum numbers, the total angular momentum, and the coupled-states
*    projection index are equal to the values stored in j, jtot, and nu
*    plus 1/2
*  variables in common block /cobspt/
*    lammin:   array containing minimum value of lambda for each term
*    lammax:   array containing maximum value of lambda for each term
*    mproj:    array containing the order of the reduced rotation matrix
*              elements for each term.  lammin can not be less than mproj.
*              for homonuclear molecules, the allowed values of lambda for
*              each term range from lammin to lammax in steps of 2
*  variable in common block /cocent/
*    cent:      array containing centrifugal barrier of each channel
*  variable in common block /coeint/
*    eint:      array containing channel energies (in hartree)
*  variables in common block /coered/
*    ered:      collision energy in atomic units (hartrees)
*    rmu:       collision reduced mass in atomic units
*               (mass of electron = 1)
*  variable in common block /conlam/
*    nlam:      the number of angular coupling terms actually used
*    nlammx:    the maximum number of angular coupling terms
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx
*  variable in common block /cov2/
*    nv2max:    maximum core memory allocated for the v2 matrix
*    v2:        lower triangle of nonzero angular coupling matrix elements
*               only nonzero elements are stored
*  variable in common block /coiv2/
*   iv2:        matrix address of v2 matrix for each non-zero element
* --------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical ihomo, flaghf, csflag, clist, flagsu, bastst
      include "common/parbas"
      common /cov2/ nv2max, junkv, v2(1)
      common /coiv2/ iv2(1)
      common /conlam/  nlam, nlammx,lamnum(1)
      common /cocent/ cent(1)
      common /coeint/ eint(1)
      common /coered/ ered, rmu
      common /coconv/ econv, xmconv
      dimension j(1), l(1), jhold(1), ehold(1), sc1(1), sc2(1), sc3(1),
     :          sc4(1), ishold(1), is(1)
*   econv is conversion factor from cm-1 to hartrees
*   xmconv is converson factor from amu to atomic units
      zero = 0.d0
      one = 1.d0
      two = 2.d0
        if (.not. flagsu) write (9,30) rcut
30      format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=',
     :          f6.2)
      n=1
      nlevop=1
      nlevel=1
      eint(1)=0
      do 40 i=1, 1
        cent(i)=zero
        j(i)=0
        jhold(i)=0
        l(i)=0
        
        is(i)=(-1)**i
        ishold(i)=is(i)
        ehold(i)=eint(i)
40    continue
      if (nu .eq. numin) then
        ntop = max(n, nlevel)
*  ntop is the maximum row dimension of all matrices passed in the
*  call list of subroutines propag and soutpt.
*  for fps make sure this is an odd number, for faster bank access.
*  this has no effect on vax or cray
        if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax)
     :       ntop = ntop + 1
      end if
*  now list channels if requested
      if (clist) then
        write (6, 255)
        write (9, 255)
255     format(/'   N   V   J    L      EINT(CM-1)')
        do 265  i = 1, n
          write (6, 260) i, is(i), j(i), l(i), eint(i) * econv
          write (9, 260) i, is(i), j(i), l(i), eint(i) * econv
260       format (3i4, i5, f13.3)
265     continue
        write (6, 256)
256     format(/' Open channels:'//
     1          '   N   V   J      EINT(CM-1)')
        do 266  i = 1, nlevop
          write (6, 261) i, ishold(i), jhold(i),  ehold(i) * econv
261       format (3i4, f13.3)
266     continue
      end if
*  now calculate coupling matrix elements
      if (bastst) then
        write (6, 280)
        write (9, 280)
280     format (/' ILAM  LAMBDA  ICOL  IROW   I    IV2    VEE')
      end if
* i counts v2 elements
* inum counts v2 elements for given lambda
* ilam counts numver of v2 matrices
* ij is address of given v2 element in present v2 matrix
      i = 0
      ilam=0
      do 320 il = lammin(1), lammax(1), 1
        inum = 0
        ilam=ilam+1
        lb=il
        if(ilam.gt.nlammx) then
          write(6,311) ilam
311       format(/' ILAM.GT.NLAMMX IN BAUSR')
          call exit
        end if
        do 310  icol= 1, n
          do 300  irow = icol, n
            ij = ntop * (icol - 1) +irow
            vee=zero            
            if (il .eq. 1) then
              if (icol .eq. 1 .and. irow .eq. 1) vee=1.d0
            endif
            if (il .eq. 2) then
              if (icol .eq. 2 .and. irow .eq. 2) vee=1.d0
            endif
            if (il .eq. 3) then
              if (irow .eq. 2 .and. icol .eq. 1) vee=1.d0
            endif
            if (vee .eq. zero) goto 300
              i = i + 1
              inum = inum + 1
              if (i .gt. nv2max) goto 325
                v2(i) = vee
                iv2(i) = ij
                if (bastst) then
                  write (6, 290) ilam, lb, icol, irow, i, iv2(i),
     :                           vee
                  write (9, 290) ilam, lb, icol, irow, i, iv2(i),
     :                           vee
290               format (i4, 2i7, 2i6, i6, g17.8)
                endif
300       continue
310     continue
        lamnum(ilam) = inum
        if (bastst) then
          write (6, 315) ilam, lamnum(ilam)
          write (9, 315) ilam, lamnum(ilam)
315       format ('ILAM=',i3,' LAMNUM(ILAM) = ',i3)
        end if
320   continue
325   if ( i.gt. nv2max) then
        write (6, 350) i, nv2max
        write (9, 350) i, nv2max
350     format (' *** NUMBER OF NONZERO V2 ELEMENTS = ',i6,
     :           ' .GT. NV2MAX=',i6,'; ABORT ***')
        if (bastst) then
          return
        else
          call exit
        end if
      end if
      if (clist) then
        write (6, 360) i
        write (9, 360) i
360     format (' ** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS',
     :           i6)
      end if
      return
      end
