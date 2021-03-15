*  System:  H + CO2(X) RCCSD(T)
*
*   calculation of potential energy curves by p.dagdigian 
*  
*   written by p. dagdigian
*   current revision date:  11-mar-2015
      include "common/syusr"
      include "common/bausr"
      include "common/ground"
* ------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(33)
      include "common/parpot"
      econv=219474.6d0
      potnam='H-CO2 RCCSD(T)'
      print *, potnam
1     print *, 'R (bohr):'
      read (5, *, end=99) r
      call pot(vv0,r)
      if (r.le.0.d0) goto 99
*     vlm coefficient is returned in atomic units (hartree)
*     convert from atomic units for printout
      write (6, 100) vv0*econv, (econv*vvl(i), i=1,6)
100   format(7(1pe16.8))
      goto 1
99    rr=3.0d0
      dr=0.1d0
      open (unit=12,file='hco2_vlms.dat')
      write(12,109)
109   format(' %R/bohr V00  ...')
      do i=1,250
        call pot(vv0,rr)
        write (12,110) rr,vv0*econv, (econv*vvl(j),j=1,6)
110     format(f7.2,7(1pe16.8))
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
      potnam='H-CO2 RCCSD(T)'
      npot=1
      nterm=1
      lammin(1)=2
      lammax(1)=12
      mproj(1)=0
      ipotsy = 2
*
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
*  count number of anisotropic terms
      nlam = 0
      do il=lammin(1),lammax(1),ipotsy
        nlam = nlam + 1
      enddo
      nlammx = nlam
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

*  variable in common block /conlam/
*    nlammx:    the maximum number of anisotropic terms
*    nlam:      the total number of angular coupling terms
*  variable in common block /covvl/
*    vvl        array to store r-dependence of each angular term in the
*               potential
* 
* author:  paul dagdigian
* latest revision date:  11-mar-2015
* ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical ljunk, ihomo, csflag, lljunk
      include "common/parbas"
      dimension v(8)
      dimension csplin(37,8)
      dimension rr(35), vl(245),vec(35)
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(6)
*
*  35 values or R
      data rr /
     +  2.500,  2.750,  3.000,  3.250,  3.500,
     +  3.750,  4.000,  4.250,  4.500,  4.750,
     +  5.000,  5.250,  5.500,  5.750,  6.000,
     +  6.250,  6.500,  6.750,  7.000,  7.250,
     +  7.500,  7.750,  8.000,  8.500,  9.000,
     +  9.500, 10.000, 10.500, 11.000, 11.500,
     + 12.000, 13.000, 14.000, 15.000, 20.000 /
*  values of the vlam coefficients
      data vl /  
     +  4.3406982e+04,  3.0078349e+04,  2.3142415e+04,  1.4847399e+04,
     +  1.0086316e+04,  7.2972329e+03,  5.4795115e+03,  4.0646877e+03,
     +  2.9113884e+03,  2.0164591e+03,  1.3543587e+03,  8.8111309e+02,
     +  5.5199059e+02,  3.2890598e+02,  1.8163998e+02,  8.7351270e+01,
     +  2.9187097e+01, -4.9091012e+00, -2.3415029e+01, -3.2110628e+01,
     + -3.4870961e+01, -3.4239054e+01, -3.1814598e+01, -2.5114878e+01,
     + -1.8701847e+01, -1.3630449e+01, -9.8740336e+00, -7.1584460e+00,
     + -5.2472939e+00, -3.8883872e+00, -2.9016193e+00, -1.6757062e+00,
     + -9.8478212e-01, -6.2023096e-01, -3.5471463e-02,
     +  1.1254623e+05,  7.1628736e+04,  5.5703273e+04,  3.1015153e+04,
     +  1.9352967e+04,  1.4061080e+04,  1.1421225e+04,  9.1950193e+03,
     +  7.0211468e+03,  5.1398734e+03,  3.6522870e+03,  2.5346829e+03,
     +  1.7215539e+03,  1.1444082e+03,  7.4339286e+02,  4.7038346e+02,
     +  2.8818394e+02,  1.6917666e+02,  9.3208262e+01,  4.6013383e+01,
     +  1.7674736e+01,  1.4163259e+00, -7.2812761e+00, -1.2696284e+01,
     + -1.1608361e+01, -9.0312239e+00, -6.5753736e+00, -4.6896795e+00,
     + -3.3008974e+00, -2.4374385e+00, -1.6989056e+00, -9.1982562e-01,
     + -5.4440354e-01, -3.0279045e-01, -3.9082200e-02,
     +  1.0927333e+05,  6.4546119e+04,  5.0073721e+04,  2.2030164e+04,
     +  1.0129395e+04,  5.7740604e+03,  4.9348348e+03,  4.3970286e+03,
     +  3.5223405e+03,  2.6216094e+03,  1.8727897e+03,  1.3043316e+03,
     +  8.9119555e+02,  5.9899059e+02,  3.9639521e+02,  2.5826143e+02,
     +  1.6549136e+02,  1.0418170e+02,  6.4254622e+01,  3.8682636e+01,
     +  2.2582152e+01,  1.2658976e+01,  6.7160427e+00,  1.0874705e+00,
     + -5.9497791e-01, -8.6983055e-01, -7.4767641e-01, -5.4350965e-01,
     + -3.9238235e-01, -2.2426636e-01, -1.4838764e-01, -7.1476666e-02,
     + -2.9186607e-03,  1.5200519e-02,  1.1981719e-02,
     +  8.5770676e+04,  5.3422302e+04,  4.3581418e+04,  1.9207919e+04,
     +  8.2309122e+03,  3.1079997e+03,  2.0863143e+03,  1.8506902e+03,
     +  1.4566547e+03,  1.0428596e+03,  7.1343482e+02,  4.7749074e+02,
     +  3.1532428e+02,  2.0623767e+02,  1.3379210e+02,  8.6046470e+01,
     +  5.4844540e+01,  3.4669007e+01,  2.1724762e+01,  1.3492002e+01,
     +  8.3046946e+00,  5.0102162e+00,  2.9647209e+00,  1.0273213e+00,
     +  2.4066399e-01,  8.5494690e-02,  1.8651668e-02, -4.4943549e-02,
     + -6.4681462e-03,  3.3558598e-03, -2.9757964e-02, -1.8418127e-02,
     + -3.5142910e-02,  6.7244259e-03,  4.4544615e-03,
     +  4.7244358e+04,  3.2726576e+04,  2.9517003e+04,  1.4660669e+04,
     +  7.4411957e+03,  2.6039654e+03,  1.2470469e+03,  8.9231275e+02,
     +  6.1515581e+02,  3.9300532e+02,  2.4314605e+02,  1.4888702e+02,
     +  9.1045447e+01,  5.5735106e+01,  3.4168086e+01,  2.0940801e+01,
     +  1.2803681e+01,  7.8002692e+00,  4.7364633e+00,  2.8876779e+00,
     +  1.7699023e+00,  1.1086649e+00,  6.4599941e-01,  2.3842353e-01,
     +  5.5855702e-02,  2.9580922e-02,  1.7173202e-02,  1.3520661e-02,
     + -8.4744019e-03,  4.5045222e-02, -1.0116614e-02,  1.2856473e-02,
     +  1.8616834e-02,  8.0824187e-03,  8.4913629e-04,
     +  1.5454389e+04,  1.2118253e+04,  1.2883245e+04,  7.3257636e+03,
     +  4.9445058e+03,  1.8129463e+03,  7.7521082e+02,  4.6892428e+02,
     +  2.7529091e+02,  1.5086544e+02,  8.1657740e+01,  4.4324551e+01,
     +  2.4276459e+01,  1.3468884e+01,  7.5305886e+00,  4.3237965e+00,
     +  2.4956561e+00,  1.4046679e+00,  8.1800803e-01,  4.5777289e-01,
     +  2.6212637e-01,  1.9249224e-01,  1.7489314e-01,  3.1780541e-02,
     +  2.2704190e-02,  1.5367503e-02, -1.0902946e-02,  8.9215315e-03,
     +  3.0670005e-02,  7.3943462e-03,  5.9241603e-04,  5.3661639e-03,
     + -6.4142472e-04, -9.4216291e-04, -1.0327569e-02,
     +  2.3422056e+03,  2.1332882e+03,  2.9447791e+03,  1.9607254e+03,
     +  2.4149294e+03,  9.2316554e+02,  4.2306929e+02,  2.6227830e+02,
     +  1.3824917e+02,  6.4376098e+01,  2.9936110e+01,  1.4044621e+01,
     +  6.6532693e+00,  3.2233595e+00,  1.6661549e+00,  8.6566614e-01,
     +  4.2881829e-01,  1.8366578e-01,  9.8427586e-02,  1.8183072e-02,
     + -3.2279459e-02, -4.6614110e-02, -6.6513656e-02, -1.8672761e-02,
     +  1.3128248e-02, -9.3205479e-03,  1.2426859e-02, -6.2588209e-03,
     + -9.0098172e-03, -3.1499962e-02, -5.4168048e-05, -1.1119610e-02,
     + -1.0427306e-02, -3.3079862e-03, -3.1799011e-03 /
*
* spline fit
      if (ifirst .eq. 0) then
* spline fit of the vl coefficients
         ind=1
         do ilam=1,7
           call dcopy(35,vl(ind),1,vec,1)
*    evaluate derivative at first point
           der1=(vec(2)-vec(1))/(rr(2)-rr(1))
           call dspline(rr,vec,35,der1,0d0,csplin(1,ilam))
           ind = ind + 35
         enddo
         ifirst = 1
       end if
* r^-6 fit to at R = 15 bohr for isotropic part of potential
       c6sum = -7.0648e+06
* switching function for long-range
       switch_lr=0.5*(tanh(0.5*(r - 18.d0)) + 1.d0)
* determine splined coefficients at R=r
       ind=1
       do ilam=1,7
         call dcopy(35,vl(ind),1,vec,1)
         call dsplint(rr,vec,csplin(1,ilam),35,r,vvx)
* kill anisotropic terms at large R
         vvx = (1.d0 - switch_lr)*vvx
         if (ilam.eq.1) then
* merge with asymptotic form
            vvx = vvx + switch_lr*c6sum/(r**6)
         endif
         v(ilam)=vvx
         call dcopy(35,vl(ind),1,vec,1)
         ind = ind + 35
       enddo
       call dcopy(6,v(2),1,vvl(1),1)
* convert to hartree
       econv=1.d0/219474.6d0
       call dscal(6,econv,vvl,1)
* isotropic term 
       vv0 = v(1)*econv
*
       return
       end
*===========================eof===============================