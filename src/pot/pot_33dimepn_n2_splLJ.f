*  System:  3,3-dimethylpentane-N2 spherically averaged splined potential
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
      potnam='33DiMePentane-N2 splined pot (data from Jasper)'
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
      open (unit=12,file='33dimepn_n2_splined_pec.txt')
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
      potnam='33DiMePentane-N2 splined pot (data from Jasper)'
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
     + 7.6000000e+00, 7.6235347e+00, 7.7235347e+00, 7.8235347e+00, 
     + 7.9235347e+00, 8.0235347e+00, 8.1235347e+00, 8.2235347e+00, 
     + 8.3235347e+00, 8.4235347e+00, 8.5235347e+00, 8.6235347e+00, 
     + 8.7235347e+00, 8.8235347e+00, 8.9235347e+00, 9.0235347e+00, 
     + 9.1235347e+00, 9.2235347e+00, 9.3235347e+00, 9.4235347e+00, 
     + 9.5235347e+00, 9.6235347e+00, 9.7235347e+00, 9.8235347e+00, 
     + 9.9235347e+00, 1.0023535e+01, 1.0123535e+01, 1.0223535e+01, 
     + 1.0323535e+01, 1.0423535e+01, 1.0523535e+01, 1.0623535e+01, 
     + 1.0723535e+01, 1.0823535e+01, 1.0923535e+01, 1.1023535e+01, 
     + 1.1123535e+01, 1.1223535e+01, 1.1323535e+01, 1.1423535e+01, 
     + 1.1523535e+01, 1.1623535e+01, 1.1723535e+01, 1.1823535e+01, 
     + 1.1923535e+01, 1.2023535e+01, 1.2123535e+01, 1.2223535e+01, 
     + 1.2323535e+01, 1.2423535e+01, 1.2523535e+01, 1.2623535e+01, 
     + 1.2723535e+01, 1.2823535e+01, 1.2923535e+01, 1.3023535e+01, 
     + 1.3123535e+01, 1.3223535e+01, 1.3323535e+01, 1.3423535e+01, 
     + 1.3523535e+01, 1.3623535e+01, 1.3723535e+01, 1.3823535e+01, 
     + 1.3923535e+01, 1.5000000e+01, 2.0000000e+01, 2.5000000e+01, 
     + 3.0000000e+01, 3.5000000e+01, 4.0000000e+01 /
*  splined pot (data from Jasper)
      data vl /
     + 6.5487419e+05, 5.5900139e+05, 4.7716424e+05, 4.0730796e+05, 
     + 3.4767856e+05, 2.9677883e+05, 2.5333076e+05, 2.1624344e+05, 
     + 1.8458565e+05, 1.5756253e+05, 1.3449556e+05, 1.1480557e+05, 
     + 9.7998173e+04, 8.3651356e+04, 7.1404897e+04, 6.0951304e+04, 
     + 5.2028105e+04, 4.4411252e+04, 3.7909498e+04, 3.2359592e+04, 
     + 2.7622186e+04, 2.3578331e+04, 2.0126491e+04, 1.7179996e+04, 
     + 1.4664865e+04, 1.2517946e+04, 1.0685334e+04, 9.1210136e+03, 
     + 7.7857080e+03, 7.5010000e+03, 6.3136578e+03, 5.2203182e+03, 
     + 4.3014138e+03, 3.5349140e+03, 2.8923085e+03, 2.3375559e+03, 
     + 1.8889510e+03, 1.4901984e+03, 1.1787625e+03, 9.1189456e+02, 
     + 7.0031233e+02, 5.1245960e+02, 3.6170604e+02, 2.4586201e+02, 
     + 1.5583678e+02, 8.2537443e+01, 1.6871070e+01, -3.6389945e+01, 
     + -7.5159913e+01, -1.0429760e+02, -1.2578500e+02, -1.4120527e+02, 
     + -1.5182430e+02, -1.5865511e+02, -1.6250880e+02, -1.6403506e+02, 
     + -1.6375426e+02, -1.6208311e+02, -1.5935510e+02, -1.5583681e+02, 
     + -1.5174100e+02, -1.4723711e+02, -1.4245969e+02, -1.3751518e+02, 
     + -1.3248735e+02, -1.2744170e+02, -1.2242906e+02, -1.1748839e+02, 
     + -1.1264918e+02, -1.0793326e+02, -1.0335635e+02, -9.8929310e+01, 
     + -9.4659105e+01, -9.0549634e+01, -8.6602379e+01, -8.2816934e+01, 
     + -7.9191434e+01, -7.5722904e+01, -7.2407532e+01, -6.9240898e+01, 
     + -6.6218153e+01, -6.3334161e+01, -6.0583615e+01, -5.7961128e+01, 
     + -5.5461300e+01, -5.3078781e+01, -5.0808305e+01, -4.8644727e+01, 
     + -4.6583048e+01, -4.4618426e+01, -4.2746194e+01, -4.0961865e+01, 
     + -3.9261132e+01, -2.5253420e+01, -4.5118271e+00, -1.1828984e+00, 
     + -3.9615378e-01, -1.5710264e-01, -7.0506942e-02 / 
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