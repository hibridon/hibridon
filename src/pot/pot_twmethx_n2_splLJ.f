*  System:  2-methylhexane-N2 spherically averaged splined potential
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
      potnam='2-MeHexane-N2 splined pot (data from Jasper)'
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
      open (unit=12,file='twmehx_n2_splined_pec.txt')
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
      potnam='2-MeHexane-N2 splined pot (data from Jasper)'
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
     + 7.5148755e+00, 7.6148755e+00, 7.7148755e+00, 7.8148755e+00, 
     + 7.9148755e+00, 8.0148755e+00, 8.1148755e+00, 8.2148755e+00, 
     + 8.3148755e+00, 8.4148755e+00, 8.5148755e+00, 8.6148755e+00, 
     + 8.7148755e+00, 8.8148755e+00, 8.9148755e+00, 9.0148755e+00, 
     + 9.1148755e+00, 9.2148755e+00, 9.3148755e+00, 9.4148755e+00, 
     + 9.5148755e+00, 9.6148755e+00, 9.7148755e+00, 9.8148755e+00, 
     + 9.9148755e+00, 1.0014875e+01, 1.0114875e+01, 1.0214875e+01, 
     + 1.0314875e+01, 1.0414875e+01, 1.0514875e+01, 1.0614875e+01, 
     + 1.0714875e+01, 1.0814875e+01, 1.0914875e+01, 1.1014875e+01, 
     + 1.1114875e+01, 1.1214875e+01, 1.1314875e+01, 1.1414875e+01, 
     + 1.1514875e+01, 1.1614875e+01, 1.1714875e+01, 1.1814875e+01, 
     + 1.1914875e+01, 1.2014875e+01, 1.2114875e+01, 1.2214875e+01, 
     + 1.2314875e+01, 1.2414875e+01, 1.2514875e+01, 1.2614875e+01, 
     + 1.2714875e+01, 1.2814875e+01, 1.2914875e+01, 1.3014875e+01, 
     + 1.3114875e+01, 1.3214875e+01, 1.3314875e+01, 1.3414875e+01, 
     + 1.3514875e+01, 1.3614875e+01, 1.3714875e+01, 1.3814875e+01, 
     + 1.3914875e+01, 1.5000000e+01, 2.0000000e+01, 2.5000000e+01, 
     + 3.0000000e+01, 3.5000000e+01, 4.0000000e+01 /
*  splined pot (data from Jasper)
      data vl /
     + 7.8673202e+05, 6.6282041e+05, 5.5842509e+05, 4.7047220e+05, 
     + 3.9637204e+05, 3.3394278e+05, 2.8134624e+05, 2.3703373e+05, 
     + 1.9970052e+05, 1.6824736e+05, 1.4174812e+05, 1.1942255e+05, 
     + 1.0061330e+05, 8.4766537e+04, 7.1415666e+04, 6.0167579e+04, 
     + 5.0691085e+04, 4.2707154e+04, 3.5980706e+04, 3.0313685e+04, 
     + 2.5539229e+04, 2.1516757e+04, 1.8127832e+04, 1.5272668e+04, 
     + 1.2867197e+04, 1.0840592e+04, 9.1331809e+03, 7.6946897e+03, 
     + 7.5010000e+03, 6.2154516e+03, 5.2835388e+03, 4.3854443e+03, 
     + 3.6307459e+03, 2.9977120e+03, 2.4935870e+03, 2.0557249e+03, 
     + 1.6421045e+03, 1.3309488e+03, 1.0463825e+03, 8.1316568e+02, 
     + 6.3849780e+02, 4.8302314e+02, 3.4811098e+02, 2.3392821e+02, 
     + 1.4029423e+02, 6.7028478e+01, 1.3950362e+01, -3.1256068e+01, 
     + -7.0149197e+01, -9.9636790e+01, -1.2160903e+02, -1.3758390e+02, 
     + -1.4878222e+02, -1.5618660e+02, -1.6058838e+02, -1.6262510e+02, 
     + -1.6281041e+02, -1.6155810e+02, -1.5920129e+02, -1.5600795e+02, 
     + -1.5219331e+02, -1.4792994e+02, -1.4335582e+02, -1.3858095e+02, 
     + -1.3369261e+02, -1.2875967e+02, -1.2383612e+02, -1.1896386e+02, 
     + -1.1417505e+02, -1.0949395e+02, -1.0493847e+02, -1.0052141e+02, 
     + -9.6251473e+01, -9.2134071e+01, -8.8172031e+01, -8.4366121e+01, 
     + -8.0715500e+01, -7.7218080e+01, -7.3870819e+01, -7.0669960e+01, 
     + -6.7611223e+01, -6.4689962e+01, -6.1901283e+01, -5.9240153e+01, 
     + -5.6701471e+01, -5.4280135e+01, -5.1971089e+01, -4.9769359e+01, 
     + -4.7670087e+01, -4.5668548e+01, -4.3760164e+01, -4.1940521e+01, 
     + -4.0205367e+01, -2.5799706e+01, -4.6151059e+00, -1.2100619e+00, 
     + -4.0525365e-01, -1.6071154e-01, -7.2126609e-02 /  
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