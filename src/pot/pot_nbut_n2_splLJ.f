*  System:  n-butane-N2 spherically averaged splined potential
*
*   from A. Jasper (jun-2015)
*   computed vdW Rmin/eps and spherically averaged repulsivze potential
*   splined repulsive and vdW potentials
*  
*   written by p. dagdigian
*   current revision date:  23-jun-2015
      include "common/syusr"
      include "common/bausr"
      include "common/ground"
* ------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(15)
      include "common/parpot"
      econv=219474.6d0
      potnam='nButane-N2 splined pot (data from Jasper)'
      print *, potnam
1     print *, 'R (bohr):'
      read (5, *, end=99) r
      call pot(vv0,r)
      if (r.le.0.d0) goto 99
*     vlm coefficient is returned in atomic units (hartree)
*     convert from atomic units for printout
      write (6, 100) vv0*econv
100   format(1pe16.8)
      goto 1
99    rr=5.d0
      dr=0.2d0
      open (unit=12,file='nbut_n2_splined_pec.txt')
      write(12,109)
109   format(' %R/bohr V0')
      do i=1,250
        call pot(vv0,rr)
        write (12,110) rr,vv0*econv
110     format(f7.2,1pe16.8)
        rr = rr + dr
      enddo
      close(12)
      end
* --------------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      character*(*) filnam
      include "common/parbas"
      include "common/parpot" 
      common /conlam/ nlam, nlammx, lamnum(2)
      common /cosysi/ nscode, isicod, nterm
      potnam='nButane-N2 splined pot (data from Jasper)'
      npot=1
      nterm=1
      lammin(1)=2
      lammax(1)=2
      mproj(1)=0
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* --------------------------------------------------------------------------
      subroutine pot (vv0, r)

*  subroutine to calculate the r-dependent coefficients in the
*  collision of a homonuclear diatomic with a structureless target
*  in units of hartree for energy and bohr for distance

*  on return:
*  vv0 contains the isotropic term (n=0) in the potential
*  the coefficients for each angular term in the coupling potential
*  [ vvl(i) for i = 1, nlam ] are returned in common block /covvl/
*  vvl(1) contains the anisotropic (n=2) term in the potential

*  variable in common block /conlam/
*    nlammx:    the maximum number of anisotropic terms
*    nlam:      the total number of angular coupling terms
*  variable in common block /covvl/
*    vvl        array to store r-dependence of each angular term in the
*               potential

* author:  paul dagdigian
* latest revision date:  23-jun-2015
* ----------------------------------------------------------------------
* CH4-N2 spherically averaged potential - exp6 fit
*
      implicit double precision (a-h,o-z)
      common /conlam/ nlam,nlammx
      include "common/parbas"
      common /covvl/ vvl(1)
      dimension rr(99),vl(99)
      dimension csplin(99),vec(99)
      data ifirst /0/
*
*  99 values of R
      data rr /
     + 4.8000000e+00, 4.9000000e+00, 5.0000000e+00, 5.1000000e+00, 
     + 5.2000000e+00, 5.3000000e+00, 5.4000000e+00, 5.5000000e+00, 
     + 5.6000000e+00, 5.7000000e+00, 5.8000000e+00, 5.9000000e+00, 
     + 6.0000000e+00, 6.1000000e+00, 6.2000000e+00, 6.3000000e+00, 
     + 6.4000000e+00, 6.5000000e+00, 6.6000000e+00, 6.7000000e+00, 
     + 6.7028600e+00, 6.8028600e+00, 6.9028600e+00, 7.0028600e+00, 
     + 7.1028600e+00, 7.2028600e+00, 7.3028600e+00, 7.4028600e+00, 
     + 7.5028600e+00, 7.6028600e+00, 7.7028600e+00, 7.8028600e+00, 
     + 7.9028600e+00, 8.0028600e+00, 8.1028600e+00, 8.2028600e+00, 
     + 8.3028600e+00, 8.4028600e+00, 8.5028600e+00, 8.6028600e+00, 
     + 8.7028600e+00, 8.8028600e+00, 8.9028600e+00, 9.0028600e+00, 
     + 9.1028600e+00, 9.2028600e+00, 9.3028600e+00, 9.4028600e+00, 
     + 9.5028600e+00, 9.6028600e+00, 9.7028600e+00, 9.8028600e+00, 
     + 9.9028600e+00, 1.0002860e+01, 1.0102860e+01, 1.0202860e+01, 
     + 1.0302860e+01, 1.0402860e+01, 1.0502860e+01, 1.0602860e+01, 
     + 1.0702860e+01, 1.0802860e+01, 1.0902860e+01, 1.1002860e+01, 
     + 1.1102860e+01, 1.1202860e+01, 1.1302860e+01, 1.1402860e+01, 
     + 1.1502860e+01, 1.1602860e+01, 1.1702860e+01, 1.1802860e+01, 
     + 1.1902860e+01, 1.2002860e+01, 1.2102860e+01, 1.2202860e+01, 
     + 1.2302860e+01, 1.2402860e+01, 1.2502860e+01, 1.2602860e+01, 
     + 1.2702860e+01, 1.2802860e+01, 1.2902860e+01, 1.3002860e+01, 
     + 1.3102860e+01, 1.3202860e+01, 1.3302860e+01, 1.3402860e+01, 
     + 1.3502860e+01, 1.3602860e+01, 1.3702860e+01, 1.3802860e+01, 
     + 1.3902860e+01, 1.5000000e+01, 2.0000000e+01, 2.5000000e+01, 
     + 3.0000000e+01, 3.5000000e+01, 4.0000000e+01 /
*  splined pot (data from Jasper)
      data vl /
     + 2.1078891e+05, 1.7689416e+05, 1.4844967e+05, 1.2457903e+05, 
     + 1.0454678e+05, 8.7735712e+04, 7.3627852e+04, 6.1788529e+04, 
     + 5.1852964e+04, 4.3515032e+04, 3.6517836e+04, 3.0645786e+04, 
     + 2.5717959e+04, 2.1582523e+04, 1.8112064e+04, 1.5199653e+04, 
     + 1.2755556e+04, 1.0704469e+04, 8.9831956e+03, 7.5387021e+03, 
     + 7.5010000e+03, 6.1860324e+03, 5.1754301e+03, 4.2753020e+03, 
     + 3.5100219e+03, 2.8551939e+03, 2.3135991e+03, 1.8654724e+03, 
     + 1.4913653e+03, 1.1729747e+03, 9.2423353e+02, 7.1085295e+02, 
     + 5.3218236e+02, 3.8449956e+02, 2.6414435e+02, 1.6797382e+02, 
     + 9.2847568e+01, 3.5625184e+01, -8.3508972e+00, -4.7595926e+01, 
     + -7.6857001e+01, -9.8344194e+01, -1.1371983e+02, -1.2430109e+02, 
     + -1.3113265e+02, -1.3504364e+02, -1.3669234e+02, -1.3660144e+02, 
     + -1.3518586e+02, -1.3277471e+02, -1.2962869e+02, -1.2595400e+02, 
     + -1.2191326e+02, -1.1763434e+02, -1.1321729e+02, -1.0873995e+02, 
     + -1.0426246e+02, -9.9830763e+01, -9.5479520e+01, -9.1234406e+01, 
     + -8.7113952e+01, -8.3131029e+01, -7.9294042e+01, -7.5607889e+01, 
     + -7.2074727e+01, -6.8694593e+01, -6.5465896e+01, -6.2385818e+01, 
     + -5.9450625e+01, -5.6655927e+01, -5.3996871e+01, -5.1468308e+01, 
     + -4.9064911e+01, -4.6781280e+01, -4.4612014e+01, -4.2551769e+01, 
     + -4.0595306e+01, -3.8737519e+01, -3.6973464e+01, -3.5298371e+01, 
     + -3.3707660e+01, -3.2196939e+01, -3.0762016e+01, -2.9398893e+01, 
     + -2.8103765e+01, -2.6873017e+01, -2.5703219e+01, -2.4591117e+01, 
     + -2.3533633e+01, -2.2527850e+01, -2.1571011e+01, -2.0660509e+01, 
     + -1.9793881e+01, -1.2589515e+01, -2.2463591e+00, -5.8894070e-01, 
     + -1.9723725e-01, -7.8218394e-02, -3.5104060e-02 /  
*
* spline fit
      if (ifirst .eq. 0) then
* spline fit of the vl coefficients
         call dcopy(99,vl,1,vec,1)
*    evaluate derivative at first point
         der1=(vec(2)-vec(1))/(rr(2)-rr(1))
         call dspline(rr,vec,99,der1,0d0,csplin)
         ifirst = 1
       end if
* r^-6 fit to at R = 30 bohr 
       c6sum = vl(97)*rr(97)**6
* switching function for long-range
       switch_lr=0.5*(tanh(0.5*(r - 25.d0)) + 1.d0)
* determine splined coefficients at R=r
       call dcopy(99,vl,1,vec,1)
       call dsplint(rr,vec,csplin,99,r,vvx)
* merge with asymptotic form
       if (r.gt.25.d0) then
         vvx = (1.d0 - switch_lr)*vvx 
     +     + switch_lr*c6sum/(r**6)
       endif
* convert to hartree
       econv=1.d0/219474.6d0
       call dscal(1,econv,vvl,1)
* isotropic term 
       vv0 = vvx*econv
       vvl(1) = 0.d0
*
       return
       end
*===========================eof===============================