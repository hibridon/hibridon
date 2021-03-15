*  System:  2,3-dimethylpentane-N2 spherically averaged splined potential
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
      potnam='23DiMePentane-N2 splined pot (data from Jasper)'
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
      open (unit=12,file='23dimepn_n2_splined_pec.txt')
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
      potnam='23DiMePentane-N2 splined pot (data from Jasper)'
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
     + 6.8000000e+00, 6.9000000e+00, 7.0000000e+00, 7.1000000e+00, 
     + 7.2000000e+00, 7.3000000e+00, 7.4000000e+00, 7.5000000e+00, 
     + 7.6000000e+00, 7.6038816e+00, 7.7038816e+00, 7.8038816e+00, 
     + 7.9038816e+00, 8.0038816e+00, 8.1038816e+00, 8.2038816e+00, 
     + 8.3038816e+00, 8.4038816e+00, 8.5038816e+00, 8.6038816e+00, 
     + 8.7038816e+00, 8.8038816e+00, 8.9038816e+00, 9.0038816e+00, 
     + 9.1038816e+00, 9.2038816e+00, 9.3038816e+00, 9.4038816e+00, 
     + 9.5038816e+00, 9.6038816e+00, 9.7038816e+00, 9.8038816e+00, 
     + 9.9038816e+00, 1.0003882e+01, 1.0103882e+01, 1.0203882e+01, 
     + 1.0303882e+01, 1.0403882e+01, 1.0503882e+01, 1.0603882e+01, 
     + 1.0703882e+01, 1.0803882e+01, 1.0903882e+01, 1.1003882e+01, 
     + 1.1103882e+01, 1.1203882e+01, 1.1303882e+01, 1.1403882e+01, 
     + 1.1503882e+01, 1.1603882e+01, 1.1703882e+01, 1.1803882e+01, 
     + 1.1903882e+01, 1.2003882e+01, 1.2103882e+01, 1.2203882e+01, 
     + 1.2303882e+01, 1.2403882e+01, 1.2503882e+01, 1.2603882e+01, 
     + 1.2703882e+01, 1.2803882e+01, 1.2903882e+01, 1.3003882e+01, 
     + 1.3103882e+01, 1.3203882e+01, 1.3303882e+01, 1.3403882e+01, 
     + 1.3503882e+01, 1.3603882e+01, 1.3703882e+01, 1.3803882e+01, 
     + 1.3903882e+01, 1.5000000e+01, 2.0000000e+01, 2.5000000e+01, 
     + 3.0000000e+01, 3.5000000e+01, 4.0000000e+01 /
*  splined pot (data from Jasper)
      data vl /
     + 6.9660769e+05, 5.9265851e+05, 5.0422082e+05, 4.2897999e+05, 
     + 3.6496674e+05, 3.1050567e+05, 2.6417140e+05, 2.2475122e+05, 
     + 1.9121339e+05, 1.6268015e+05, 1.3840469e+05, 1.1775167e+05, 
     + 1.0018053e+05, 8.5231391e+04, 7.2512993e+04, 6.1692460e+04, 
     + 5.2486588e+04, 4.4654435e+04, 3.7991011e+04, 3.2321917e+04, 
     + 2.7498776e+04, 2.3395354e+04, 1.9904252e+04, 1.6934100e+04, 
     + 1.4407160e+04, 1.2257295e+04, 1.0428237e+04, 8.8721146e+03, 
     + 7.5482000e+03, 7.5010000e+03, 6.2888081e+03, 5.2102823e+03, 
     + 4.3299268e+03, 3.5389384e+03, 2.9212055e+03, 2.3604640e+03, 
     + 1.9063243e+03, 1.5167456e+03, 1.2001585e+03, 9.4318034e+02, 
     + 7.2401810e+02, 5.4782024e+02, 3.9826339e+02, 2.6871453e+02, 
     + 1.6107936e+02, 7.7379268e+01, 1.9635653e+01, -2.2898616e+01, 
     + -6.3970263e+01, -9.5184793e+01, -1.1853576e+02, -1.3560844e+02, 
     + -1.4767779e+02, -1.5577006e+02, -1.6071176e+02, -1.6316888e+02, 
     + -1.6367809e+02, -1.6267178e+02, -1.6049822e+02, -1.5743763e+02, 
     + -1.5371520e+02, -1.4951159e+02, -1.4497137e+02, -1.4020988e+02, 
     + -1.3531874e+02, -1.3037037e+02, -1.2542162e+02, -1.2051672e+02, 
     + -1.1568971e+02, -1.1096637e+02, -1.0636585e+02, -1.0190192e+02, 
     + -9.7584084e+01, -9.3418379e+01, -8.9408127e+01, -8.5554486e+01, 
     + -8.1856916e+01, -7.8313557e+01, -7.4921541e+01, -7.1677235e+01, 
     + -6.8576447e+01, -6.5614583e+01, -6.2786785e+01, -6.0088028e+01, 
     + -5.7513209e+01, -5.5057208e+01, -5.2714946e+01, -5.0481416e+01, 
     + -4.8351720e+01, -4.6321093e+01, -4.4384911e+01, -4.2538714e+01, 
     + -4.0778204e+01, -2.6047978e+01, -4.6597293e+00, -1.2217641e+00, 
     + -4.0917283e-01, -1.6226577e-01, -7.2824142e-02 / 
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