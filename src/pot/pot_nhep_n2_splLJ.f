*  System:  n-C7H16-N2 spherically averaged splined potential
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
      potnam='nC7H16-N2 splined pot (data from Jasper)'
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
      open (unit=12,file='nhep_n2_splined_pec.txt')
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
      potnam='nC7H16-N2 splined pot (data from Jasper)'
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
     + 7.2000000e+00, 7.3000000e+00, 7.3914763e+00, 7.4914763e+00, 
     + 7.5914763e+00, 7.6914763e+00, 7.7914763e+00, 7.8914763e+00, 
     + 7.9914763e+00, 8.0914763e+00, 8.1914763e+00, 8.2914763e+00, 
     + 8.3914763e+00, 8.4914763e+00, 8.5914763e+00, 8.6914763e+00, 
     + 8.7914763e+00, 8.8914763e+00, 8.9914763e+00, 9.0914763e+00, 
     + 9.1914763e+00, 9.2914763e+00, 9.3914763e+00, 9.4914763e+00, 
     + 9.5914763e+00, 9.6914763e+00, 9.7914763e+00, 9.8914763e+00, 
     + 9.9914763e+00, 1.0091476e+01, 1.0191476e+01, 1.0291476e+01, 
     + 1.0391476e+01, 1.0491476e+01, 1.0591476e+01, 1.0691476e+01, 
     + 1.0791476e+01, 1.0891476e+01, 1.0991476e+01, 1.1091476e+01, 
     + 1.1191476e+01, 1.1291476e+01, 1.1391476e+01, 1.1491476e+01, 
     + 1.1591476e+01, 1.1691476e+01, 1.1791476e+01, 1.1891476e+01, 
     + 1.1991476e+01, 1.2091476e+01, 1.2191476e+01, 1.2291476e+01, 
     + 1.2391476e+01, 1.2491476e+01, 1.2591476e+01, 1.2691476e+01, 
     + 1.2791476e+01, 1.2891476e+01, 1.2991476e+01, 1.3091476e+01, 
     + 1.3191476e+01, 1.3291476e+01, 1.3391476e+01, 1.3491476e+01, 
     + 1.3591476e+01, 1.3691476e+01, 1.3791476e+01, 1.3891476e+01, 
     + 1.3991476e+01, 1.5000000e+01, 2.0000000e+01, 2.5000000e+01, 
     + 3.0000000e+01, 3.5000000e+01, 4.0000000e+01 /
*  splined pot (data from Jasper)
      data vl /
     + 3.3327231e+05, 2.8788463e+05, 2.4867821e+05, 2.1481123e+05, 
     + 1.8555653e+05, 1.6028596e+05, 1.3845694e+05, 1.1960077e+05, 
     + 1.0331259e+05, 8.9242655e+04, 7.7088880e+04, 6.6590302e+04, 
     + 5.7521504e+04, 4.9687768e+04, 4.2920892e+04, 3.7075583e+04, 
     + 3.2026335e+04, 2.7664733e+04, 2.3897130e+04, 2.0642628e+04, 
     + 1.7831351e+04, 1.5402935e+04, 1.3305241e+04, 1.1493227e+04, 
     + 9.9279879e+03, 8.5759155e+03, 7.5010000e+03, 6.4028530e+03, 
     + 5.3831416e+03, 4.4753269e+03, 3.8077701e+03, 3.1297849e+03, 
     + 2.6169759e+03, 2.1446252e+03, 1.7760038e+03, 1.4203151e+03, 
     + 1.1429459e+03, 9.0923311e+02, 6.9464544e+02, 5.3005182e+02, 
     + 3.9716803e+02, 2.8612274e+02, 1.9367054e+02, 1.1660186e+02, 
     + 5.1707132e+01, -3.9716248e+00, -4.8076427e+01, -8.2010106e+01, 
     + -1.0780046e+02, -1.2703851e+02, -1.4101251e+02, -1.5076647e+02, 
     + -1.5714698e+02, -1.6084093e+02, -1.6240571e+02, -1.6229366e+02, 
     + -1.6087172e+02, -1.5843743e+02, -1.5523176e+02, -1.5144963e+02, 
     + -1.4724841e+02, -1.4275478e+02, -1.3807046e+02, -1.3327673e+02, 
     + -1.2843822e+02, -1.2360598e+02, -1.1881997e+02, -1.1411111e+02, 
     + -1.0950302e+02, -1.0501331e+02, -1.0065478e+02, -9.6436341e+01, 
     + -9.2363734e+01, -8.8440198e+01, -8.4666962e+01, -8.1043666e+01, 
     + -7.7568710e+01, -7.4239525e+01, -7.1052812e+01, -6.8004721e+01, 
     + -6.5091007e+01, -6.2307150e+01, -5.9648457e+01, -5.7110138e+01, 
     + -5.4687373e+01, -5.2375360e+01, -5.0169354e+01, -4.8064700e+01, 
     + -4.6056852e+01, -4.4141397e+01, -4.2314058e+01, -4.0570712e+01, 
     + -3.8907389e+01, -2.5826627e+01, -4.6294319e+00, -1.2140116e+00, 
     + -4.0658408e-01, -1.6123964e-01, -7.2363667e-02 /  
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