*  System:  CH4-N2 spherically averaged splined potential
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
      potnam='methane-N2 splined pot (data from Jasper)'
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
      open (unit=12,file='meth_n2_splined_pec.txt')
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
      potnam='methane-N2 splined pot (data from Jasper)'
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
     + 5.1028285e+00, 5.2028285e+00, 5.3028285e+00, 5.4028285e+00, 
     + 5.5028285e+00, 5.6028285e+00, 5.7028285e+00, 5.8028285e+00, 
     + 5.9028285e+00, 6.0028285e+00, 6.1028285e+00, 6.2028285e+00, 
     + 6.3028285e+00, 6.4028285e+00, 6.5028285e+00, 6.6028285e+00, 
     + 6.7028285e+00, 6.8028285e+00, 6.9028285e+00, 7.0028285e+00, 
     + 7.1028285e+00, 7.2028285e+00, 7.3028285e+00, 7.4028285e+00, 
     + 7.5028285e+00, 7.6028285e+00, 7.7028285e+00, 7.8028285e+00, 
     + 7.9028285e+00, 8.0028285e+00, 8.1028285e+00, 8.2028285e+00, 
     + 8.3028285e+00, 8.4028285e+00, 8.5028285e+00, 8.6028285e+00, 
     + 8.7028285e+00, 8.8028285e+00, 8.9028285e+00, 9.0028285e+00, 
     + 9.1028285e+00, 9.2028285e+00, 9.3028285e+00, 9.4028285e+00, 
     + 9.5028285e+00, 9.6028285e+00, 9.7028285e+00, 9.8028285e+00, 
     + 9.9028285e+00, 1.0002829e+01, 1.0102829e+01, 1.0202829e+01, 
     + 1.0302829e+01, 1.0402829e+01, 1.0502829e+01, 1.0602829e+01, 
     + 1.0702829e+01, 1.0802829e+01, 1.0902829e+01, 1.1002829e+01, 
     + 1.1102829e+01, 1.1202829e+01, 1.1302829e+01, 1.1402829e+01, 
     + 1.1502829e+01, 1.1602829e+01, 1.1702829e+01, 1.1802829e+01, 
     + 1.1902829e+01, 1.2002829e+01, 1.2102829e+01, 1.2202829e+01, 
     + 1.2302829e+01, 1.2402829e+01, 1.2502829e+01, 1.2602829e+01, 
     + 1.2702829e+01, 1.2802829e+01, 1.2902829e+01, 1.3002829e+01, 
     + 1.3102829e+01, 1.3202829e+01, 1.3302829e+01, 1.3402829e+01, 
     + 1.3502829e+01, 1.3602829e+01, 1.3702829e+01, 1.3802829e+01, 
     + 1.3902829e+01, 1.5000000e+01, 2.0000000e+01, 2.5000000e+01, 
     + 3.0000000e+01, 3.5000000e+01, 4.0000000e+01 /
*  splined pot (data from Jasper)
      data vl /
     + 1.2675767e+04, 1.0659381e+04, 8.9637491e+03, 7.5378486e+03, 
     + 7.5010000e+03, 6.2014386e+03, 5.0838957e+03, 4.1524986e+03, 
     + 3.3677841e+03, 2.7352041e+03, 2.1818828e+03, 1.7528710e+03, 
     + 1.3940111e+03, 1.0971261e+03, 8.5246391e+02, 6.5314057e+02, 
     + 4.9311926e+02, 3.6226285e+02, 2.5598191e+02, 1.7094190e+02, 
     + 1.0380829e+02, 5.1246540e+01, 9.9221220e+00, -2.3369526e+01, 
     + -4.8493816e+01, -6.6893232e+01, -8.0030846e+01, -8.9060867e+01, 
     + -9.4894754e+01, -9.8252784e+01, -9.9704367e+01, -9.9699641e+01, 
     + -9.8594297e+01, -9.6669137e+01, -9.4145500e+01, -9.1197476e+01, 
     + -8.7961571e+01, -8.4544393e+01, -8.1028755e+01, -7.7478537e+01, 
     + -7.3942567e+01, -7.0457711e+01, -6.7051351e+01, -6.3743358e+01, 
     + -6.0547680e+01, -5.7473600e+01, -5.4526757e+01, -5.1709952e+01, 
     + -4.9023795e+01, -4.6467227e+01, -4.4037930e+01, -4.1732655e+01, 
     + -3.9547482e+01, -3.7478032e+01, -3.5519621e+01, -3.3667395e+01, 
     + -3.1916422e+01, -3.0261775e+01, -2.8698586e+01, -2.7222091e+01, 
     + -2.5827661e+01, -2.4510825e+01, -2.3267288e+01, -2.2092937e+01, 
     + -2.0983847e+01, -1.9936285e+01, -1.8946708e+01, -1.8011757e+01, 
     + -1.7128256e+01, -1.6293207e+01, -1.5503779e+01, -1.4757304e+01, 
     + -1.4051270e+01, -1.3383314e+01, -1.2751211e+01, -1.2152871e+01, 
     + -1.1586329e+01, -1.1049739e+01, -1.0541366e+01, -1.0059581e+01, 
     + -9.6028563e+00, -9.1697536e+00, -8.7589245e+00, -8.3691022e+00, 
     + -7.9990970e+00, -7.6477915e+00, -7.3141361e+00, -6.9971447e+00, 
     + -6.6958907e+00, -6.4095033e+00, -6.1371641e+00, -5.8781035e+00, 
     + -5.6315980e+00, -3.5838373e+00, -6.4095242e-01, -1.6812941e-01, 
     + -5.6313542e-02, -2.2333015e-02, -1.0023067e-02 /  
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