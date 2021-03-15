*  System:  H + O2(X) - 4App PES, r = re
*
*   calculation of potential energy curves by p.dagdigian 
*  
*   written by p. dagdigian
*   current revision date:  4-feb-2014
      include "common/syusr"
      include "common/bausr"
      include "common/ground"
* ------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(33)
      include "common/parpot"
      econv=219474.6d0
      potnam='H+O2(X) a4App'
      print *, potnam
1     print *, 'R (bohr):'
      read (5, *, end=99) r
      call pot(vv0,r)
      if (r.le.0.d0) goto 99
*     vlm coefficient is returned in atomic units (hartree)
*     convert from atomic units for printout
      write (6, 100) vv0*econv, (econv*vvl(i), i=1,5)
100   format(6(1pe16.8))
      goto 1
99    rr=3.0d0
      dr=0.1d0
      open (unit=12,file='ho2_a4app_vlms.dat')
      write(12,109)
109   format(' %R/bohr V00  ...')
      do i=1,250
        call pot(vv0,rr)
        write (12,110) rr,vv0*econv, (econv*vvl(j),j=1,5)
110     format(f7.2,6(1pe16.8))
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
      potnam='H+O2(X) a4App'
      npot=1
      nterm=1
      lammin(1)=2
      lammax(1)=10
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
* latest revision date:  4-feb-2014
* ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical ljunk, ihomo, csflag, lljunk
      include "common/parbas"
      dimension v(6)
      dimension csplin(37,6)
      dimension rr(37), vl(222),vec(37)
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(5)
*
*  37 values or R
      data rr /
     +  2.500,  2.600,  2.700,  2.800,  2.900,
     +  3.000,  3.100,  3.200,  3.300,  3.400,
     +  3.500,  3.600,  3.700,  3.800,  3.900,
     +  4.000,  4.100,  4.200,  4.300,  4.400,
     +  4.500,  5.000,  5.500,  6.000,  6.500,
     +  7.000,  7.500,  8.000,  8.500,  9.000,
     +  9.500, 10.000, 11.000, 13.000, 15.000,
     + 20.000, 30.000 /
*  values of the vlam coefficients
      data vl /  
     +  4.3182854e+04,  3.7333987e+04,  3.2760879e+04,  2.8671620e+04,
     +  2.5272258e+04,  2.2254274e+04,  1.9535200e+04,  1.7118867e+04,
     +  1.4953047e+04,  1.3029503e+04,  1.1362799e+04,  9.7137073e+03,
     +  8.3947335e+03,  7.2022656e+03,  6.1505212e+03,  5.2609484e+03,
     +  4.4809267e+03,  3.8176182e+03,  3.2247595e+03,  2.7175925e+03,
     +  2.2885369e+03,  9.0202657e+02,  3.0973877e+02,  7.4710299e+01,
     + -6.9105826e+00, -2.7730490e+01, -2.7546899e+01, -2.1729835e+01,
     + -1.6265447e+01, -1.1731526e+01, -8.2970742e+00, -5.8534268e+00,
     + -2.9245934e+00, -6.1281657e-01,  6.4829804e-02,  4.2212040e-01,
     +  4.8186712e-01,
     + -3.3296100e+04, -2.6942824e+04, -2.2683284e+04, -1.8897384e+04,
     + -1.6254796e+04, -1.4079188e+04, -1.2197521e+04, -1.0719881e+04,
     + -9.3811625e+03, -8.2042996e+03, -7.1791747e+03, -6.3454141e+03,
     + -5.4866453e+03, -4.6838658e+03, -4.0074735e+03, -3.4527489e+03,
     + -2.9558004e+03, -2.4981254e+03, -2.1417391e+03, -1.8185650e+03,
     + -1.5387110e+03, -6.4131759e+02, -2.3992287e+02, -7.5847073e+01,
     + -1.4927755e+01,  4.1122890e+00,  7.7016864e+00,  7.6511553e+00,
     +  5.4800274e+00,  3.6781164e+00,  2.5176190e+00,  1.7111909e+00,
     +  8.4289274e-01,  2.6660300e-01,  9.7011605e-02,  3.1667089e-02,
     +  1.8329084e-02,
     +  8.8946983e+03,  5.1297294e+03,  2.8210365e+03,  3.6738887e+02,
     + -8.9561982e+02, -1.6951132e+03, -2.1730311e+03, -2.3683077e+03,
     + -2.2731993e+03, -2.0488678e+03, -2.0791266e+03, -1.5829358e+03,
     + -1.4663698e+03, -1.3107937e+03, -1.1128792e+03, -9.9986939e+02,
     + -9.0486713e+02, -8.4744322e+02, -7.3633703e+02, -6.3615459e+02,
     + -5.7337467e+02, -2.7022216e+02, -1.3307747e+02, -6.3602296e+01,
     + -2.9371142e+01, -1.3273132e+01, -6.2632491e+00, -3.3968957e+00,
     + -1.2992577e+00, -4.0631149e-01, -1.1132504e-01, -4.6842213e-02,
     + -2.0908885e-02, -4.9047375e-02, -2.0764906e-02, -3.6010813e-02,
     + -2.9798032e-02,
     + -2.2399284e+03, -2.1277802e+03, -2.1016200e+03, -8.2455618e+02,
     + -1.8327091e+02,  3.5319403e+02,  7.4301901e+02,  8.8349239e+02,
     +  1.0412877e+03,  6.6666802e+02,  6.9818936e+02,  7.2734306e+02,
     +  6.5962651e+02,  6.1798184e+02,  4.2085286e+02,  3.9905607e+02,
     +  3.3503294e+02,  2.3208859e+02,  2.6156848e+02,  2.0209137e+02,
     +  1.8104496e+02,  8.0111377e+01,  3.4312977e+01,  1.4313104e+01,
     +  5.5391056e+00,  2.2691260e+00,  1.3400029e+00, -3.1572414e-01,
     + -2.1648841e-01,  1.4885449e-02,  4.1046932e-02,  1.4630358e-02,
     + -2.0872527e-02,  6.0710015e-02,  2.5191819e-02,  2.4132646e-02,
     +  2.5589046e-02,
     + -6.0922915e+03, -2.8131054e+03, -3.0905112e+02, -2.1475206e+02,
     +  5.8793653e+01,  1.8072723e+01, -8.4157482e+01, -1.0613691e+02,
     + -2.0444612e+02, -1.8384376e+02, -6.6588894e+01, -2.8986013e+02,
     + -2.7375637e+02, -1.5070308e+02, -1.4588696e+02, -1.4301938e+02,
     + -9.8824050e+01, -4.6106024e+01, -3.0229952e+01, -3.3279432e+01,
     + -7.9568769e+00, -1.0026680e+01, -2.8509342e+00,  2.0333062e-01,
     +  1.2107011e+00,  9.6847568e-01,  3.7734313e-01,  3.4433586e-01,
     +  1.3270275e-01,  2.9517143e-02, -2.8343852e-02, -1.4143855e-02,
     +  1.5724583e-02, -2.7558604e-02, -6.5540959e-02, -5.6976194e-02,
     + -5.0082920e-02,
     +  5.6954359e+03,  2.9662256e+03,  8.4565716e+02,  4.9620446e+02,
     +  1.0252568e+02, -1.0189035e+01, -5.5588866e+01,  8.2587525e+01,
     + -5.1195939e+01,  2.2332605e+02,  1.3205850e+02,  1.6412302e+02,
     +  1.4971101e+02, -2.7994293e+01,  4.5545491e+01,  2.7019510e+01,
     +  2.9799830e+01,  6.4447703e+01, -1.2219828e+01,  1.7485764e+01,
     +  1.4599684e+00,  4.2932934e+00,  2.1225935e+00,  4.5729196e-01,
     + -4.9594646e-01, -6.1218259e-01, -4.1741397e-01,  8.6864202e-02,
     +  1.5499782e-01,  1.0487241e-01,  3.3686939e-02,  7.4201661e-02,
     +  6.9919075e-02, -4.4123055e-02,  5.1658351e-03,  5.5029152e-03,
     + -2.1595399e-03 /
*
* spline fit
      if (ifirst .eq. 0) then
* spline fit of the vl coefficients
         ind=1
         do ilam=1,6
           call dcopy(37,vl(ind),1,vec,1)
*    evaluate derivative at first point
           der1=(vec(2)-vec(1))/(rr(2)-rr(1))
           call dspline(rr,vec,37,der1,0d0,csplin(1,ilam))
           ind = ind + 37
         enddo
         ifirst = 1
       end if
* r^-6 fit to at R = 10 bohr for isotropic part of potential
       c6sum = -5.853427e+06
* switching function for long-range
       switch_lr=0.5*(tanh(0.5*(r - 25.d0)) + 1.d0)
* determine splined coefficients at R=r
       ind=1
       do ilam=1,6
         call dcopy(37,vl(ind),1,vec,1)
         call dsplint(rr,vec,csplin(1,ilam),37,r,vvx)
* kill anisotropic terms at large R
         vvx = (1.d0 - switch_lr)*vvx
         if (ilam.eq.1) then
* merge with asymptotic form
            vvx = vvx + switch_lr*c6sum/(r**6)
         endif
         v(ilam)=vvx
         call dcopy(37,vl(ind),1,vec,1)
         ind = ind + 37
       enddo
       call dcopy(5,v(2),1,vvl(1),1)
* convert to hartree
       econv=1.d0/219474.6d0
       call dscal(5,econv,vvl,1)
* isotropic term 
       vv0 = v(1)*econv
*
       return
       end
*===========================eof===============================