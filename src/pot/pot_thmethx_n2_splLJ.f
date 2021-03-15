*  System:  3-methylhexane-N2 spherically averaged splined potential
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
      potnam='3-MeHexane-N2 splined pot (data from Jasper)'
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
      open (unit=12,file='thmehx_n2_splined_pec.txt')
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
      potnam='3-MeHexane-N2 splined pot (data from Jasper)'
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
     + 7.5360404e+00, 7.6360404e+00, 7.7360404e+00, 7.8360404e+00, 
     + 7.9360404e+00, 8.0360404e+00, 8.1360404e+00, 8.2360404e+00, 
     + 8.3360404e+00, 8.4360404e+00, 8.5360404e+00, 8.6360404e+00, 
     + 8.7360404e+00, 8.8360404e+00, 8.9360404e+00, 9.0360404e+00, 
     + 9.1360404e+00, 9.2360404e+00, 9.3360404e+00, 9.4360404e+00, 
     + 9.5360404e+00, 9.6360404e+00, 9.7360404e+00, 9.8360404e+00, 
     + 9.9360404e+00, 1.0036040e+01, 1.0136040e+01, 1.0236040e+01, 
     + 1.0336040e+01, 1.0436040e+01, 1.0536040e+01, 1.0636040e+01, 
     + 1.0736040e+01, 1.0836040e+01, 1.0936040e+01, 1.1036040e+01, 
     + 1.1136040e+01, 1.1236040e+01, 1.1336040e+01, 1.1436040e+01, 
     + 1.1536040e+01, 1.1636040e+01, 1.1736040e+01, 1.1836040e+01, 
     + 1.1936040e+01, 1.2036040e+01, 1.2136040e+01, 1.2236040e+01, 
     + 1.2336040e+01, 1.2436040e+01, 1.2536040e+01, 1.2636040e+01, 
     + 1.2736040e+01, 1.2836040e+01, 1.2936040e+01, 1.3036040e+01, 
     + 1.3136040e+01, 1.3236040e+01, 1.3336040e+01, 1.3436040e+01, 
     + 1.3536040e+01, 1.3636040e+01, 1.3736040e+01, 1.3836040e+01, 
     + 1.3936040e+01, 1.5000000e+01, 2.0000000e+01, 2.5000000e+01, 
     + 3.0000000e+01, 3.5000000e+01, 4.0000000e+01 /
*  splined pot (data from Jasper)
      data vl /
     + 6.4783767e+05, 5.5042015e+05, 4.6765163e+05, 3.9732929e+05, 
     + 3.3758155e+05, 2.8681829e+05, 2.4368846e+05, 2.0704421e+05, 
     + 1.7591028e+05, 1.4945806e+05, 1.2698355e+05, 1.0788861e+05, 
     + 9.1665038e+04, 7.7881059e+04, 6.6169824e+04, 5.6219646e+04, 
     + 4.7765709e+04, 4.0583020e+04, 3.4480415e+04, 2.9295480e+04, 
     + 2.4890220e+04, 2.1147394e+04, 1.7967390e+04, 1.5265573e+04, 
     + 1.2970038e+04, 1.1019690e+04, 9.3626222e+03, 7.9547333e+03, 
     + 7.5010000e+03, 6.2786501e+03, 5.2805213e+03, 4.4260234e+03, 
     + 3.6636346e+03, 2.9578735e+03, 2.4966259e+03, 2.0133719e+03, 
     + 1.6481769e+03, 1.3172978e+03, 1.0363024e+03, 8.1053073e+02, 
     + 6.1385372e+02, 4.5271966e+02, 3.2066471e+02, 2.1375679e+02, 
     + 1.2831284e+02, 6.0649789e+01, 7.0845564e+00, -4.0018869e+01, 
     + -7.6867738e+01, -1.0472125e+02, -1.2538804e+02, -1.4032155e+02, 
     + -1.5069095e+02, -1.5743741e+02, -1.6131901e+02, -1.6294649e+02, 
     + -1.6281195e+02, -1.6131176e+02, -1.5876501e+02, -1.5542831e+02, 
     + -1.5150770e+02, -1.4716830e+02, -1.4254210e+02, -1.3773418e+02, 
     + -1.3282788e+02, -1.2788888e+02, -1.2296856e+02, -1.1810672e+02, 
     + -1.1333383e+02, -1.0867278e+02, -1.0414038e+02, -9.9748527e+01, 
     + -9.5505222e+01, -9.1415320e+01, -8.7481192e+01, -8.3703255e+01, 
     + -8.0080392e+01, -7.6610300e+01, -7.3289777e+01, -7.0114946e+01, 
     + -6.7081440e+01, -6.4184550e+01, -6.1419346e+01, -5.8780771e+01, 
     + -5.6263718e+01, -5.3863086e+01, -5.1573829e+01, -4.9390994e+01, 
     + -4.7309742e+01, -4.5325373e+01, -4.3433340e+01, -4.1629256e+01, 
     + -3.9908903e+01, -2.5838967e+01, -4.6221551e+00, -1.2119104e+00, 
     + -4.0587274e-01, -1.6095705e-01, -7.2236795e-02  /  
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