*  System:  i-butane-N2 spherically averaged splined potential
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
      potnam='iButane-N2 splined pot (data from Jasper)'
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
      open (unit=12,file='ibut_n2_splined_pec.txt')
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
      potnam='iButane-N2 splined pot (data from Jasper)'
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
* latest revision date:  22-jun-2015
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
     + 6.7115527e+00, 6.8115527e+00, 6.9115527e+00, 7.0115527e+00, 
     + 7.1115527e+00, 7.2115527e+00, 7.3115527e+00, 7.4115527e+00, 
     + 7.5115527e+00, 7.6115527e+00, 7.7115527e+00, 7.8115527e+00, 
     + 7.9115527e+00, 8.0115527e+00, 8.1115527e+00, 8.2115527e+00, 
     + 8.3115527e+00, 8.4115527e+00, 8.5115527e+00, 8.6115527e+00, 
     + 8.7115527e+00, 8.8115527e+00, 8.9115527e+00, 9.0115527e+00, 
     + 9.1115527e+00, 9.2115527e+00, 9.3115527e+00, 9.4115527e+00, 
     + 9.5115527e+00, 9.6115527e+00, 9.7115527e+00, 9.8115527e+00, 
     + 9.9115527e+00, 1.0011553e+01, 1.0111553e+01, 1.0211553e+01, 
     + 1.0311553e+01, 1.0411553e+01, 1.0511553e+01, 1.0611553e+01, 
     + 1.0711553e+01, 1.0811553e+01, 1.0911553e+01, 1.1011553e+01, 
     + 1.1111553e+01, 1.1211553e+01, 1.1311553e+01, 1.1411553e+01, 
     + 1.1511553e+01, 1.1611553e+01, 1.1711553e+01, 1.1811553e+01, 
     + 1.1911553e+01, 1.2011553e+01, 1.2111553e+01, 1.2211553e+01, 
     + 1.2311553e+01, 1.2411553e+01, 1.2511553e+01, 1.2611553e+01, 
     + 1.2711553e+01, 1.2811553e+01, 1.2911553e+01, 1.3011553e+01, 
     + 1.3111553e+01, 1.3211553e+01, 1.3311553e+01, 1.3411553e+01, 
     + 1.3511553e+01, 1.3611553e+01, 1.3711553e+01, 1.3811553e+01, 
     + 1.3911553e+01, 1.5000000e+01, 2.0000000e+01, 2.5000000e+01, 
     + 3.0000000e+01, 3.5000000e+01, 4.0000000e+01 /
*  splined pot (data from Jasper)
      data vl /
     + 1.8156813e+05, 1.5368869e+05, 1.3009009e+05, 1.1011500e+05, 
     + 9.3207052e+04, 7.8895287e+04, 6.6781065e+04, 5.6526960e+04, 
     + 4.7847353e+04, 4.0500483e+04, 3.4281711e+04, 2.9017819e+04, 
     + 2.4562189e+04, 2.0790711e+04, 1.7598337e+04, 1.4896146e+04, 
     + 1.2608871e+04, 1.0672803e+04, 9.0340144e+03, 7.6468587e+03, 
     + 7.5010000e+03, 6.2505691e+03, 5.1550004e+03, 4.2477362e+03, 
     + 3.4934439e+03, 2.8310740e+03, 2.2886004e+03, 1.8387832e+03, 
     + 1.4680149e+03, 1.1627550e+03, 8.8980640e+02, 6.8928101e+02, 
     + 5.1386522e+02, 3.6863728e+02, 2.5262360e+02, 1.6105539e+02, 
     + 8.9161887e+01, 3.2172358e+01, -1.4961843e+01, -5.2630744e+01, 
     + -8.0702286e+01, -1.0124413e+02, -1.1586900e+02, -1.2585480e+02, 
     + -1.3221496e+02, -1.3575377e+02, -1.3710970e+02, -1.3678962e+02, 
     + -1.3519578e+02, -1.3264714e+02, -1.2939627e+02, -1.2564283e+02, 
     + -1.2154418e+02, -1.1722395e+02, -1.1277879e+02, -1.0828380e+02, 
     + -1.0379690e+02, -9.9362237e+01, -9.5013053e+01, -9.0773861e+01, 
     + -8.6662259e+01, -8.2690372e+01, -7.8866008e+01, -7.5193583e+01, 
     + -7.1674877e+01, -6.8309621e+01, -6.5095989e+01, -6.2030974e+01, 
     + -5.9110701e+01, -5.6330664e+01, -5.3685932e+01, -5.1171291e+01, 
     + -4.8781373e+01, -4.6510748e+01, -4.4353994e+01, -4.2305761e+01, 
     + -4.0360806e+01, -3.8514027e+01, -3.6760487e+01, -3.5095427e+01, 
     + -3.3514279e+01, -3.2012667e+01, -3.0586414e+01, -2.9231539e+01, 
     + -2.7944254e+01, -2.6720962e+01, -2.5558248e+01, -2.4452879e+01, 
     + -2.3401789e+01, -2.2402081e+01, -2.1451012e+01, -2.0545992e+01, 
     + -1.9684570e+01, -1.2566255e+01, -2.2421663e+00, -5.8784093e-01, 
     + -1.9686892e-01, -7.8072324e-02, -3.5038504e-02 /  
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