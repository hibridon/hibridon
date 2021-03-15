*  System:  CH4-N2 spherically averaged - LJ fit
*
*   from A. Jasper (may-2015).  see jasper et al., JCP 141, 124313 (2014)
*  
*   written by p. dagdigian
*   current revision date:  12-jun-2015
      include "common/syusr"
      include "common/bausr"
      include "common/ground"
* ------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(15)
      include "common/parpot"
      econv=219474.6d0
      potnam='CH4-N2 sph avg - exp6 fit'
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
99    rr=3.4d0
      dr=0.2d0
      open (unit=12,file='ch4n2_LJ_pec.txt')
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
      potnam='CH4-N2 sph avg - exp6 fit'
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
* latest revision date:  12-jun-2015
* ----------------------------------------------------------------------
* CH4-N2 spherically averaged potential - exp6 fit
*
      implicit double precision (a-h,o-z)
      common /conlam/ nlam,nlammx
      include "common/parbas"
      common /covvl/ vvl(1)
      dimension rr(60),vl(60)
      dimension csplin(60),vec(60)
      data ifirst /0/
      data rconv / 0.5291771d0 /
*
*  60 values or R in Ang
      data rr /
     + 1.7000000e+00, 1.8000000e+00, 1.9000000e+00, 2.0000000e+00, 
     + 2.1000000e+00, 2.2000000e+00, 2.3000000e+00, 2.4000000e+00, 
     + 2.5000000e+00, 2.6000000e+00, 2.7000000e+00, 2.8000000e+00, 
     + 2.9000000e+00, 3.0000000e+00, 3.1000000e+00, 3.2000000e+00, 
     + 3.3000000e+00, 3.4000000e+00, 3.5000000e+00, 3.6000000e+00, 
     + 3.7000000e+00, 3.8000000e+00, 3.9000000e+00, 4.0000000e+00, 
     + 4.1000000e+00, 4.2000000e+00, 4.3000000e+00, 4.4000000e+00, 
     + 4.5000000e+00, 4.6000000e+00, 4.7000000e+00, 4.8000000e+00, 
     + 4.9000000e+00, 5.0000000e+00, 5.1000000e+00, 5.2000000e+00, 
     + 5.3000000e+00, 5.4000000e+00, 5.5000000e+00, 5.6000000e+00, 
     + 5.7000000e+00, 5.8000000e+00, 5.9000000e+00, 6.0000000e+00, 
     + 6.2500000e+00, 6.5000000e+00, 6.7500000e+00, 7.0000000e+00, 
     + 7.5000000e+00, 8.0000000e+00, 8.5000000e+00, 9.0000000e+00, 
     + 1.0000000e+01, 1.1000000e+01, 1.2000000e+01, 1.3000000e+01, 
     + 1.5000000e+01, 2.0000000e+01, 2.5000000e+01, 3.0000000e+01 /
*  exp6 pot
      data vl /  
     + 4.0411093e+06, 2.0269489e+06, 1.0536975e+06, 5.6535116e+05, 
     + 3.1193135e+05, 1.7640307e+05, 1.0193936e+05, 6.0023658e+04, 
     + 3.5911825e+04, 2.1770383e+04, 1.3332841e+04, 8.2221199e+03, 
     + 5.0860272e+03, 3.1405926e+03, 1.9233008e+03, 1.1569014e+03, 
     + 6.7276758e+02, 3.6697654e+02, 1.7471543e+02, 5.5131067e+01, 
     + -1.7769068e+01, -6.0655870e+01, -8.4292830e+01, -9.5665107e+01, 
     + -9.9321422e+01, -9.8228179e+01, -9.4318007e+01, -8.8844923e+01, 
     + -8.2615965e+01, -7.6143224e+01, -6.9744177e+01, -6.3608193e+01, 
     + -5.7840762e+01, -5.2492973e+01, -4.7581168e+01, -4.3100035e+01, 
     + -3.9031293e+01, -3.5349424e+01, -3.2025408e+01, -2.9029126e+01, 
     + -2.6330861e+01, -2.3902195e+01, -2.1716520e+01, -1.9749274e+01, 
     + -1.5644658e+01, -1.2475119e+01, -1.0014778e+01, -8.0933719e+00, 
     + -5.3884501e+00, -3.6747666e+00, -2.5615413e+00, -1.8213020e+00, 
     + -9.7001111e-01, -5.4813058e-01, -3.2538442e-01, -2.0135381e-01, 
     + -8.5349159e-02, -1.5192998e-02, -3.9828660e-03, -1.3338620e-03 / 
*
      rang = r * rconv
* spline fit
      if (ifirst .eq. 0) then
* spline fit of the vl coefficients
         call dcopy(60,vl,1,vec,1)
*    evaluate derivative at first point
         der1=(vec(2)-vec(1))/(rr(2)-rr(1))
         call dspline(rr,vec,60,der1,0d0,csplin)
         ifirst = 1
       end if
* r^-6 fit to at R = 30 bohr 
       c6sum = vl(60)*rr(60)**6
* switching function for long-range
       switch_lr=0.5*(tanh(0.5*(rang - 27.d0)) + 1.d0)
* determine splined coefficients at R=r
       call dcopy(60,vl,1,vec,1)
       call dsplint(rr,vec,csplin,60,rang,vvx)
* merge with asymptotic form
       if (r.gt.25.d0) then
         vvx = (1.d0 - switch_lr)*vvx 
     +     + switch_lr*c6sum/(rang**6)
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