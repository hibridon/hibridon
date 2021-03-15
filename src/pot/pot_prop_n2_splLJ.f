*  System:  n-C3H8-N2 spherically averaged splined potential
*
*   from A. Jasper (jun-2015)
*   computed vdW Rmin/eps and spherically averaged repulsivze potential
*   splined repulsive and vdW potentials
*  
*   written by p. dagdigian
*   current revision date:  22-jun-2015
      include "common/syusr"
      include "common/bausr"
      include "common/ground"
* ------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(15)
      include "common/parpot"
      econv=219474.6d0
      potnam='propane-N2 splined pot (data from Jasper)'
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
      open (unit=12,file='prop_n2_splined_pec.txt')
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
      potnam='propane-N2 splined pot (data from Jasper)'
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
     + 6.0000000e+00, 6.1000000e+00, 6.2000000e+00, 6.2886319e+00, 
     + 6.3886319e+00, 6.4886319e+00, 6.5886319e+00, 6.6886319e+00, 
     + 6.7886319e+00, 6.8886319e+00, 6.9886319e+00, 7.0886319e+00, 
     + 7.1886319e+00, 7.2886319e+00, 7.3886319e+00, 7.4886319e+00, 
     + 7.5886319e+00, 7.6886319e+00, 7.7886319e+00, 7.8886319e+00, 
     + 7.9886319e+00, 8.0886319e+00, 8.1886319e+00, 8.2886319e+00, 
     + 8.3886319e+00, 8.4886319e+00, 8.5886319e+00, 8.6886319e+00, 
     + 8.7886319e+00, 8.8886319e+00, 8.9886319e+00, 9.0886319e+00, 
     + 9.1886319e+00, 9.2886319e+00, 9.3886319e+00, 9.4886319e+00, 
     + 9.5886319e+00, 9.6886319e+00, 9.7886319e+00, 9.8886319e+00, 
     + 9.9886319e+00, 1.0088632e+01, 1.0188632e+01, 1.0288632e+01, 
     + 1.0388632e+01, 1.0488632e+01, 1.0588632e+01, 1.0688632e+01, 
     + 1.0788632e+01, 1.0888632e+01, 1.0988632e+01, 1.1088632e+01, 
     + 1.1188632e+01, 1.1288632e+01, 1.1388632e+01, 1.1488632e+01, 
     + 1.1588632e+01, 1.1688632e+01, 1.1788632e+01, 1.1888632e+01, 
     + 1.1988632e+01, 1.2088632e+01, 1.2188632e+01, 1.2288632e+01, 
     + 1.2388632e+01, 1.2488632e+01, 1.2588632e+01, 1.2688632e+01, 
     + 1.2788632e+01, 1.2888632e+01, 1.2988632e+01, 1.3088632e+01, 
     + 1.3188632e+01, 1.3288632e+01, 1.3388632e+01, 1.3488632e+01, 
     + 1.3588632e+01, 1.3688632e+01, 1.3788632e+01, 1.3888632e+01, 
     + 1.3988632e+01, 1.5000000e+01, 2.0000000e+01, 2.5000000e+01, 
     + 3.0000000e+01, 3.5000000e+01, 4.0000000e+01 /
*  splined pot (data from Jasper)
      data vl /
     + 9.5040535e+04, 8.0136085e+04, 6.7568981e+04, 5.6972676e+04, 
     + 4.8038105e+04, 4.0504672e+04, 3.4152647e+04, 2.8796759e+04, 
     + 2.4280793e+04, 2.0473029e+04, 1.7262407e+04, 1.4555281e+04, 
     + 1.2272693e+04, 1.0348064e+04, 8.7252595e+03, 7.5010000e+03, 
     + 6.2215002e+03, 5.1868066e+03, 4.1979548e+03, 3.4518438e+03, 
     + 2.8124038e+03, 2.2860492e+03, 1.8365368e+03, 1.4649221e+03, 
     + 1.1436325e+03, 8.9353318e+02, 6.9520295e+02, 5.2247968e+02, 
     + 3.7860677e+02, 2.6245913e+02, 1.7002465e+02, 9.7285505e+01, 
     + 4.0223846e+01, -5.4753967e+00, -4.2824693e+01, -7.0650501e+01, 
     + -9.0993526e+01, -1.0547182e+02, -1.1536403e+02, -1.2168154e+02, 
     + -1.2522488e+02, -1.2662792e+02, -1.2639255e+02, -1.2491603e+02, 
     + -1.2251252e+02, -1.1943013e+02, -1.1586441e+02, -1.1196904e+02, 
     + -1.0786435e+02, -1.0364413e+02, -9.9380935e+01, -9.5130497e+01, 
     + -9.0935102e+01, -8.6826367e+01, -8.2827444e+01, -7.8954782e+01, 
     + -7.5219542e+01, -7.1628734e+01, -6.8186121e+01, -6.4892949e+01, 
     + -6.1748530e+01, -5.8750704e+01, -5.5896212e+01, -5.3180994e+01, 
     + -5.0600419e+01, -4.8149478e+01, -4.5822920e+01, -4.3615377e+01, 
     + -4.1521446e+01, -3.9535760e+01, -3.7653038e+01, -3.5868128e+01, 
     + -3.4176032e+01, -3.2571923e+01, -3.1051166e+01, -2.9609318e+01, 
     + -2.8242139e+01, -2.6945585e+01, -2.5715813e+01, -2.4549175e+01, 
     + -2.3442210e+01, -2.2391645e+01, -2.1394379e+01, -2.0447486e+01, 
     + -1.9548198e+01, -1.8693904e+01, -1.7882141e+01, -1.7110584e+01, 
     + -1.6377042e+01, -1.5679449e+01, -1.5015856e+01, -1.4384429e+01, 
     + -1.3783437e+01, -9.0897646e+00, -1.6217676e+00, -4.2519881e-01, 
     + -1.4240051e-01, -5.6471841e-02, -2.5344312e-02 /  
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