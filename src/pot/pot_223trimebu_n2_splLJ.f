*  System:  2,2,3-trimethylbuttane-N2 spherically averaged splined potential
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
      potnam='233TriMeButane-N2 splined pot (data from Jasper)'
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
99    rr=6.5d0
      dr=0.2d0
      open (unit=12,file='223trimebu_n2_splined_pec.txt')
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
      potnam='233TriMeButane-N2 splined pot (data from Jasper)'
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
     + 6.8000000e+00, 6.9000000e+00, 7.0000000e+00, 7.1000000e+00, 
     + 7.2000000e+00, 7.3000000e+00, 7.4000000e+00, 7.5000000e+00, 
     + 7.6000000e+00, 7.6482901e+00, 7.7482901e+00, 7.8482901e+00, 
     + 7.9482901e+00, 8.0482901e+00, 8.1482901e+00, 8.2482901e+00, 
     + 8.3482901e+00, 8.4482901e+00, 8.5482901e+00, 8.6482901e+00, 
     + 8.7482901e+00, 8.8482901e+00, 8.9482901e+00, 9.0482901e+00, 
     + 9.1482901e+00, 9.2482901e+00, 9.3482901e+00, 9.4482901e+00, 
     + 9.5482901e+00, 9.6482901e+00, 9.7482901e+00, 9.8482901e+00, 
     + 9.9482901e+00, 1.0048290e+01, 1.0148290e+01, 1.0248290e+01, 
     + 1.0348290e+01, 1.0448290e+01, 1.0548290e+01, 1.0648290e+01, 
     + 1.0748290e+01, 1.0848290e+01, 1.0948290e+01, 1.1048290e+01, 
     + 1.1148290e+01, 1.1248290e+01, 1.1348290e+01, 1.1448290e+01, 
     + 1.1548290e+01, 1.1648290e+01, 1.1748290e+01, 1.1848290e+01, 
     + 1.1948290e+01, 1.2048290e+01, 1.2148290e+01, 1.2248290e+01, 
     + 1.2348290e+01, 1.2448290e+01, 1.2548290e+01, 1.2648290e+01, 
     + 1.2748290e+01, 1.2848290e+01, 1.2948290e+01, 1.3048290e+01, 
     + 1.3148290e+01, 1.3248290e+01, 1.3348290e+01, 1.3448290e+01, 
     + 1.3548290e+01, 1.3648290e+01, 1.3748290e+01, 1.3848290e+01, 
     + 1.3948290e+01, 1.5000000e+01, 2.0000000e+01, 2.5000000e+01, 
     + 3.0000000e+01, 3.5000000e+01, 4.0000000e+01 /
*  splined pot (data from Jasper)
      data vl /
     + 9.6484615e+05, 8.1358277e+05, 6.8603366e+05, 5.7848100e+05, 
     + 4.8778987e+05, 4.1131681e+05, 3.4683278e+05, 2.9245821e+05, 
     + 2.4660819e+05, 2.0794629e+05, 1.7534559e+05, 1.4785586e+05, 
     + 1.2467581e+05, 1.0512981e+05, 8.8648128e+04, 7.4750353e+04, 
     + 6.3031396e+04, 5.3149674e+04, 4.4817155e+04, 3.7790963e+04, 
     + 3.1866300e+04, 2.6870474e+04, 2.2657866e+04, 1.9105688e+04, 
     + 1.6110401e+04, 1.3584699e+04, 1.1454964e+04, 9.6591163e+03, 
     + 8.1448121e+03, 7.5010000e+03, 6.2219221e+03, 5.1592206e+03, 
     + 4.2594732e+03, 3.4698256e+03, 2.8299395e+03, 2.2873635e+03, 
     + 1.8274977e+03, 1.4518402e+03, 1.1262475e+03, 8.6911188e+02, 
     + 6.5456517e+02, 4.8207956e+02, 3.3998894e+02, 2.2379592e+02, 
     + 1.3041193e+02, 5.6748397e+01, -6.8748299e-01, -4.9114383e+01, 
     + -8.5833005e+01, -1.1302876e+02, -1.3272448e+02, -1.4652365e+02, 
     + -1.5569865e+02, -1.6125989e+02, -1.6401020e+02, -1.6458761e+02, 
     + -1.6349917e+02, -1.6114764e+02, -1.5785265e+02, -1.5386754e+02, 
     + -1.4939269e+02, -1.4458614e+02, -1.3957213e+02, -1.3444783e+02, 
     + -1.2928875e+02, -1.2415314e+02, -1.1908543e+02, -1.1411900e+02, 
     + -1.0927848e+02, -1.0458149e+02, -1.0004014e+02, -9.5662161e+01, 
     + -9.1451850e+01, -8.7410827e+01, -8.3538633e+01, -7.9833218e+01, 
     + -7.6291331e+01, -7.2908829e+01, -6.9680924e+01, -6.6602383e+01, 
     + -6.3667685e+01, -6.0871142e+01, -5.8206999e+01, -5.5669505e+01, 
     + -5.3252976e+01, -5.0951833e+01, -4.8760642e+01, -4.6674130e+01, 
     + -4.4687208e+01, -4.2794978e+01, -4.0992742e+01, -3.9276001e+01, 
     + -3.7640459e+01, -2.4433368e+01, -4.3600996e+00, -1.1430536e+00, 
     + -3.8280792e-01, -1.5180999e-01, -6.8131617e-02 / 
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