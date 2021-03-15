*  System:  CH4-N2 spherically averaged - exp6 fit
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
      open (unit=12,file='ch4n2_exp6_pec.txt')
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
     + 2.7127698e+05, 1.9220762e+05, 1.3585417e+05, 9.5799065e+04, 
     + 6.7392488e+04, 4.7286928e+04, 3.3083255e+04, 2.3067687e+04, 
     + 1.6019094e+04, 1.1069117e+04, 7.6012907e+03, 5.1786096e+03, 
     + 3.4917016e+03, 2.3218295e+03, 1.5145371e+03, 9.6090361e+02, 
     + 5.8422997e+02, 3.3059230e+02, 1.6214732e+02, 5.2392137e+01, 
     + -1.7188798e+01, -5.9493291e+01, -8.3474353e+01, -9.5320060e+01, 
     + -9.9284718e+01, -9.8273624e+01, -9.4253882e+01, -8.8542541e+01, 
     + -8.2008279e+01, -7.5212271e+01, -6.8506300e+01, -6.2100864e+01, 
     + -5.6112274e+01, -5.0595057e+01, -4.5564121e+01, -4.1009787e+01, 
     + -3.6907898e+01, -3.3226517e+01, -2.9930289e+01, -2.6983217e+01, 
     + -2.4350354e+01, -2.1998781e+01, -1.9898112e+01, -1.8020689e+01, 
     + -1.4146697e+01, -1.1197977e+01, -8.9364815e+00, -7.1878422e+00, 
     + -4.7531225e+00, -3.2273813e+00, -2.2433045e+00, -1.5920244e+00, 
     + -8.4606890e-01, -4.7758390e-01, -2.8334680e-01, -1.7528537e-01, 
     + -7.4277664e-02, -1.3219828e-02, -3.4654990e-03, -1.1605890e-03 / 
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
       switch_lr=0.5*(tanh(0.5*(rang - 25.d0)) + 1.d0)
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