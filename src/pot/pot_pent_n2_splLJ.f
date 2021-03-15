*  System:  nC5H12-N2 spherically averaged splined potential
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
      potnam='npentane-N2 splined pot (data from Jasper)'
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
      open (unit=12,file='npent_n2_splined_pec.txt')
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
      potnam='npentane-N2 splined pot (data from Jasper)'
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
     + 6.8000000e+00, 6.9000000e+00, 6.9672327e+00, 7.0672327e+00, 
     + 7.1672327e+00, 7.2672327e+00, 7.3672327e+00, 7.4672327e+00, 
     + 7.5672327e+00, 7.6672327e+00, 7.7672327e+00, 7.8672327e+00, 
     + 7.9672327e+00, 8.0672327e+00, 8.1672327e+00, 8.2672327e+00, 
     + 8.3672327e+00, 8.4672327e+00, 8.5672327e+00, 8.6672327e+00, 
     + 8.7672327e+00, 8.8672327e+00, 8.9672327e+00, 9.0672327e+00, 
     + 9.1672327e+00, 9.2672327e+00, 9.3672327e+00, 9.4672327e+00, 
     + 9.5672327e+00, 9.6672327e+00, 9.7672327e+00, 9.8672327e+00, 
     + 9.9672327e+00, 1.0067233e+01, 1.0167233e+01, 1.0267233e+01, 
     + 1.0367233e+01, 1.0467233e+01, 1.0567233e+01, 1.0667233e+01, 
     + 1.0767233e+01, 1.0867233e+01, 1.0967233e+01, 1.1067233e+01, 
     + 1.1167233e+01, 1.1267233e+01, 1.1367233e+01, 1.1467233e+01, 
     + 1.1567233e+01, 1.1667233e+01, 1.1767233e+01, 1.1867233e+01, 
     + 1.1967233e+01, 1.2067233e+01, 1.2167233e+01, 1.2267233e+01, 
     + 1.2367233e+01, 1.2467233e+01, 1.2567233e+01, 1.2667233e+01, 
     + 1.2767233e+01, 1.2867233e+01, 1.2967233e+01, 1.3067233e+01, 
     + 1.3167233e+01, 1.3267233e+01, 1.3367233e+01, 1.3467233e+01, 
     + 1.3567233e+01, 1.3667233e+01, 1.3767233e+01, 1.3867233e+01, 
     + 1.3967233e+01, 1.5000000e+01, 2.0000000e+01, 2.5000000e+01, 
     + 3.0000000e+01, 3.5000000e+01, 4.0000000e+01 /
*  splined pot (data from Jasper)
      data vl /
     + 2.6414291e+05, 2.2411480e+05, 1.9015253e+05, 1.6133690e+05, 
     + 1.3688797e+05, 1.1614402e+05, 9.8543601e+04, 8.3610343e+04, 
     + 7.0940065e+04, 6.0189836e+04, 5.1068692e+04, 4.3329763e+04, 
     + 3.6763588e+04, 3.1192449e+04, 2.6465558e+04, 2.2454978e+04, 
     + 1.9052159e+04, 1.6165003e+04, 1.3715365e+04, 1.1636944e+04, 
     + 9.8734859e+03, 8.3772618e+03, 7.5010000e+03, 6.2683467e+03, 
     + 5.2256044e+03, 4.3718061e+03, 3.6450476e+03, 2.9725212e+03, 
     + 2.4330826e+03, 1.9697053e+03, 1.6213457e+03, 1.2925426e+03, 
     + 1.0126548e+03, 7.8779071e+02, 6.0760983e+02, 4.5542337e+02, 
     + 3.2457535e+02, 2.1461546e+02, 1.2567398e+02, 5.7881215e+01, 
     + 1.1367472e+01, -2.8910225e+01, -6.4746507e+01, -9.1653751e+01, 
     + -1.1147426e+02, -1.2567830e+02, -1.3544114e+02, -1.4170292e+02, 
     + -1.4521617e+02, -1.4658341e+02, -1.4628703e+02, -1.4471314e+02, 
     + -1.4217053e+02, -1.3890588e+02, -1.3511592e+02, -1.3095715e+02, 
     + -1.2655369e+02, -1.2200357e+02, -1.1738384e+02, -1.1275458e+02, 
     + -1.0816229e+02, -1.0364249e+02, -9.9221908e+01, -9.4920230e+01, 
     + -9.0751496e+01, -8.6725252e+01, -8.2847477e+01, -7.9121332e+01, 
     + -7.5547765e+01, -7.2126006e+01, -6.8853959e+01, -6.5728526e+01, 
     + -6.2745864e+01, -5.9901588e+01, -5.7190943e+01, -5.4608929e+01, 
     + -5.2150412e+01, -4.9810203e+01, -4.7583124e+01, -4.5464057e+01, 
     + -4.3447983e+01, -4.1530010e+01, -3.9705394e+01, -3.7969554e+01, 
     + -3.6318081e+01, -3.4746746e+01, -3.3251498e+01, -3.1828469e+01, 
     + -3.0473971e+01, -2.9184491e+01, -2.7956686e+01, -2.6787385e+01, 
     + -2.5673572e+01, -1.6806689e+01, -3.0023951e+00, -7.8719933e-01, 
     + -2.6363579e-01, -1.0455014e-01, -4.6921635e-02 /  
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