*  System:  C2H6-N2 spherically averaged splined potential
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
      potnam='ethane-N2 splined pot (data from Jasper)'
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
      open (unit=12,file='eth_n2_splined_pec.txt')
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
      potnam='ethane-N2 splined pot (data from Jasper)'
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
     + 5.6000000e+00, 5.7000000e+00, 5.7870985e+00, 5.8870985e+00, 
     + 5.9870985e+00, 6.0870985e+00, 6.1870985e+00, 6.2870985e+00, 
     + 6.3870985e+00, 6.4870985e+00, 6.5870985e+00, 6.6870985e+00, 
     + 6.7870985e+00, 6.8870985e+00, 6.9870985e+00, 7.0870985e+00, 
     + 7.1870985e+00, 7.2870985e+00, 7.3870985e+00, 7.4870985e+00, 
     + 7.5870985e+00, 7.6870985e+00, 7.7870985e+00, 7.8870985e+00, 
     + 7.9870985e+00, 8.0870985e+00, 8.1870985e+00, 8.2870985e+00, 
     + 8.3870985e+00, 8.4870985e+00, 8.5870985e+00, 8.6870985e+00, 
     + 8.7870985e+00, 8.8870985e+00, 8.9870985e+00, 9.0870985e+00, 
     + 9.1870985e+00, 9.2870985e+00, 9.3870985e+00, 9.4870985e+00, 
     + 9.5870985e+00, 9.6870985e+00, 9.7870985e+00, 9.8870985e+00, 
     + 9.9870985e+00, 1.0087098e+01, 1.0187098e+01, 1.0287098e+01, 
     + 1.0387098e+01, 1.0487098e+01, 1.0587098e+01, 1.0687098e+01, 
     + 1.0787098e+01, 1.0887098e+01, 1.0987098e+01, 1.1087098e+01, 
     + 1.1187098e+01, 1.1287098e+01, 1.1387098e+01, 1.1487098e+01, 
     + 1.1587098e+01, 1.1687098e+01, 1.1787098e+01, 1.1887098e+01, 
     + 1.1987098e+01, 1.2087098e+01, 1.2187098e+01, 1.2287098e+01, 
     + 1.2387098e+01, 1.2487098e+01, 1.2587098e+01, 1.2687098e+01, 
     + 1.2787098e+01, 1.2887098e+01, 1.2987098e+01, 1.3087098e+01, 
     + 1.3187098e+01, 1.3287098e+01, 1.3387098e+01, 1.3487098e+01, 
     + 1.3587098e+01, 1.3687098e+01, 1.3787098e+01, 1.3887098e+01, 
     + 1.3987098e+01, 1.5000000e+01, 2.0000000e+01, 2.5000000e+01, 
     + 3.0000000e+01, 3.5000000e+01, 4.0000000e+01 /
*  splined pot (data from Jasper)
      data vl /
     + 4.3215441e+04, 3.6190365e+04, 3.0307281e+04, 2.5380547e+04, 
     + 2.1254701e+04, 1.7799549e+04, 1.4906065e+04, 1.2482944e+04, 
     + 1.0453724e+04, 8.7543733e+03, 7.5010000e+03, 6.1702855e+03, 
     + 5.1782025e+03, 4.2033234e+03, 3.3922495e+03, 2.7995843e+03, 
     + 2.2453794e+03, 1.7968155e+03, 1.4476029e+03, 1.1290187e+03, 
     + 8.8009108e+02, 6.7972159e+02, 5.0527545e+02, 3.6799428e+02, 
     + 2.6408936e+02, 1.8479041e+02, 1.2132697e+02, 6.4928578e+01, 
     + 6.8247669e+00, -3.1591549e+01, -5.8511730e+01, -7.8265583e+01, 
     + -9.2398957e+01, -1.0213402e+02, -1.0843797e+02, -1.1207675e+02, 
     + -1.1365719e+02, -1.1366004e+02, -1.1246602e+02, -1.1037643e+02, 
     + -1.0762932e+02, -1.0441244e+02, -1.0087344e+02, -9.7128005e+01, 
     + -9.3266313e+01, -8.9358235e+01, -8.5457439e+01, -8.1604692e+01, 
     + -7.7830508e+01, -7.4157254e+01, -7.0600852e+01, -6.7172134e+01, 
     + -6.3877934e+01, -6.0721958e+01, -5.7705491e+01, -5.4827948e+01, 
     + -5.2087331e+01, -4.9480583e+01, -4.7003869e+01, -4.4652812e+01, 
     + -4.2422664e+01, -4.0308447e+01, -3.8305071e+01, -3.6407412e+01, 
     + -3.4610379e+01, -3.2908969e+01, -3.1298299e+01, -2.9773635e+01, 
     + -2.8330412e+01, -2.6964244e+01, -2.5670936e+01, -2.4446482e+01, 
     + -2.3287070e+01, -2.2189079e+01, -2.1149073e+01, -2.0163800e+01, 
     + -1.9230183e+01, -1.8345316e+01, -1.7506453e+01, -1.6711005e+01, 
     + -1.5956532e+01, -1.5240732e+01, -1.4561439e+01, -1.3916612e+01, 
     + -1.3304328e+01, -1.2722778e+01, -1.2170260e+01, -1.1645171e+01, 
     + -1.1146003e+01, -1.0671335e+01, -1.0219831e+01, -9.7902336e+00, 
     + -9.3813579e+00, -6.1843933e+00, -1.1042342e+00, -2.8954761e-01, 
     + -9.6972759e-02, -3.8456746e-02, -1.7259243e-02 /  
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