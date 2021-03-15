*  system:  Ar-CH4 PES
*  written by P. Dagdigian, Jul 2015
*
*  PES calculated by Heijmen et al. (JCP 107, 902 (1997).  Fortran code
*  to generate the tetrahedral expansion coefficients provided by
*  A. van der Avoird (July 2015)
*
*  expansion coefficients determined by least squares fit to (theta, phi)
*  points on R grid
*
*  the PES is fitted with 5 angular terms
*
      include "common/syusr"
      include "common/ground"
      include "common/bausr"
* ------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(4)
      include "common/parpot"
      s4pi = sqrt ( 4.d0 * acos(-1.d0) )
      econv=219474.6d0
      potnam='Ar-CH4 Nijmegen 1997 LsqFit'
      print *, potnam
1     print *, 'R (bohr):'
      read (5, *, end=99) r
      call pot(vv0,r)
      if (r.le.0.d0) goto 99
*     vlm coefficient is returned in atomic units (hartree)
*     convert from atomic units for printout
      write (6, 100) vv0*econv, (econv*vvl(i), i=1,4)
100   format(6(1pe16.8))
      goto 1
99    rr=3.5d0
      dr=0.5d0
      open (unit=12,file='arch4_vlms_lsqfit.txt')
      write(12,109)
109   format(' %R/bohr V0   V3    V4     V6    V7')
      do i=1,40
        call pot(vv0,rr)
        write (12,110) rr,vv0*econv, (econv*vvl(j),j=1,4)
110     format(f7.2,5(1pe16.8))
        rr = rr + dr
      enddo
      close(12)
      end
* ------------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
      implicit double precision (a-h,o-z)
      character*(*) filnam
      include "common/parbas"
      include "common/parpot" 
      common /conlam/ nlam, nlammx, lamnum(4)
      common /cosysi/ nscode, isicod, nterm
      potnam='Ar-CH4 Nijmegen 1997 LsqFit'
      ibasty = 24
*
      nterm = 4
      lammin(1) = 3
      lammin(2) = 4
      lammin(3) = 6
      lammin(4) = 7
      do i=1,nterm
        mproj(i) = 0
        lammax(i) = lammin(i)
        lamnum(i) = 1
      enddo
*  total number of anisotropic terms
      nlammx = nterm
      nlam = nterm
*
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* ----------------------------------------------------------------------
      subroutine pot (vv0, r)
*  subroutine to calculate the r-dependent coefficients
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:          interparticle distance
*  on return:
*    vv0         contains isotropic term (Y00)
*  variable in common block /covvl/
*    vvl:        vector of length 4 to store r-dependence of each term
*                in potential expansion
*  variable in common block /coloapot/
*    s4pi:       normalization factor for isotropic potential
*
*  uses linear least squares routines from lapack
*
* author:  paul dagdigian
* latest revision date:  29-jul-2015
* ----------------------------------------------------------------------
* this pes is a fit to 50 values of R
*
      implicit double precision (a-h,o-z)
      logical ljunk, ihomo, csflag, lljunk
      include "common/parbas"
      dimension v(5)
      dimension csplin(50,5)
      dimension rr(50), vl(250),vec(50)
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(4)
* here are the 50 values of R
      data rr/
     +  3.000,  3.250,  3.500,  3.750,  4.000,
     +  4.250,  4.500,  4.750,  5.000,  5.250,
     +  5.500,  5.750,  6.000,  6.250,  6.500,
     +  6.750,  7.000,  7.250,  7.500,  7.750,
     +  8.000,  8.250,  8.500,  8.750,  9.000,
     +  9.250,  9.500,  9.750, 10.000, 10.250,
     + 10.500, 10.750, 11.000, 11.250, 11.500,
     + 11.750, 12.000, 13.000, 14.000, 15.000,
     + 16.000, 17.000, 18.000, 19.000, 20.000,
     + 21.000, 22.000, 23.000, 24.000, 25.000 /
* here are column ordered angular expansion coefficients (5 total) for the
* potential at each of 50 values of R (250 values)
      data vl/
     +  1.5495619e+04,  4.3645109e+04,  4.6160481e+04,  3.9195183e+04,
     +  3.0044545e+04,  2.1657144e+04,  1.4951808e+04,  9.9746966e+03,
     +  6.4541967e+03,  4.0505779e+03,  2.4560758e+03,  1.4252867e+03,
     +  7.7574101e+02,  3.7774612e+02,  1.4206322e+02,  8.8296060e+00,
     + -6.1256868e+01, -9.3497884e+01, -1.0386969e+02, -1.0225385e+02,
     + -9.4617825e+01, -8.4458314e+01, -7.3740744e+01, -6.3501131e+01,
     + -5.4225747e+01, -4.6086688e+01, -3.9085344e+01, -3.3137826e+01,
     + -2.8124357e+01, -2.3916671e+01, -2.0392241e+01, -1.7440836e+01,
     + -1.4966761e+01, -1.2888789e+01, -1.1138964e+01, -9.6609692e+00,
     + -8.4083980e+00, -4.9854983e+00, -3.0955933e+00, -1.9970356e+00,
     + -1.3301303e+00, -9.1031933e-01, -6.3781692e-01, -4.5619315e-01,
     + -3.3230425e-01, -2.4604488e-01, -1.8487276e-01, -1.4076822e-01,
     + -1.0848919e-01, -8.4540431e-02,
     + -1.6041534e+04,  1.1671725e+04,  1.9980926e+04,  1.9867522e+04,
     +  1.6681597e+04,  1.2869098e+04,  9.4267855e+03,  6.6629890e+03,
     +  4.5846928e+03,  3.0869132e+03,  2.0398575e+03,  1.3248856e+03,
     +  8.4602275e+02,  5.3067017e+02,  3.2622306e+02,  1.9570448e+02,
     +  1.1371475e+02,  6.3127380e+01,  3.2574786e+01,  1.4617172e+01,
     +  4.4493127e+00, -9.9208330e-01, -3.6333881e+00, -4.6679715e+00,
     + -4.8223325e+00, -4.5293563e+00, -4.0398289e+00, -3.4934009e+00,
     + -2.9632113e+00, -2.4835866e+00, -2.0669776e+00, -1.7141324e+00,
     + -1.4200648e+00, -1.1774517e+00, -9.7847799e-01, -8.1577033e-01,
     + -6.8280448e-01, -3.5077937e-01, -1.9365112e-01, -1.1337123e-01,
     + -6.9380694e-02, -4.3936175e-02, -2.8612199e-02, -1.9088319e-02,
     + -1.3012210e-02, -9.0456065e-03, -6.4016569e-03, -4.6053033e-03,
     + -3.3631038e-03, -2.4899993e-03,
     +  3.1169888e+04,  7.0046121e+03, -2.6024351e+03, -5.5139044e+03,
     + -5.5981378e+03, -4.6806074e+03, -3.5673317e+03, -2.5730124e+03,
     + -1.7878270e+03, -1.2082769e+03, -7.9863876e+02, -5.1786736e+02,
     + -3.2991447e+02, -2.0651427e+02, -1.2686318e+02, -7.6263100e+01,
     + -4.4626213e+01, -2.5180080e+01, -1.3458465e+01, -6.5606457e+00,
     + -2.6285030e+00, -4.8739732e-01,  5.9534735e-01,  1.0701917e+00,
     +  1.2095647e+00,  1.1749690e+00,  1.0597130e+00,  9.1581733e-01,
     +  7.7076275e-01,  6.3778406e-01,  5.2211005e-01,  4.2467760e-01,
     +  3.4429930e-01,  2.7888797e-01,  2.2612157e-01,  1.8377907e-01,
     +  1.4988945e-01,  6.9603254e-02,  3.5255864e-02,  1.9289345e-02,
     +  1.1183218e-02,  6.7599085e-03,  4.2161226e-03,  2.6971544e-03,
     +  1.7638834e-03,  1.1768254e-03,  7.9979589e-04,  5.5298552e-04,
     +  3.8850445e-04,  2.7703324e-04,
     + -1.0525680e+04, -4.3542655e+03, -1.5195086e+03, -3.0329494e+02,
     +  1.5614493e+02,  2.8124143e+02,  2.7253965e+02,  2.2080330e+02,
     +  1.6398628e+02,  1.1565399e+02,  7.8814587e+01,  5.2415980e+01,
     +  3.4238429e+01,  2.2068511e+01,  1.4090111e+01,  8.9443916e+00,
     +  5.6683525e+00,  3.6037572e+00,  2.3123267e+00,  1.5080877e+00,
     +  1.0075099e+00,  6.9446865e-01,  4.9640693e-01,  3.6851327e-01,
     +  2.8340988e-01,  2.2453662e-01,  1.8196540e-01,  1.4978396e-01,
     +  1.2448221e-01,  1.0396748e-01,  8.6973842e-02,  7.2710881e-02,
     +  6.0659214e-02,  5.0454170e-02,  4.1820449e-02,  3.4535666e-02,
     +  2.8412270e-02,  1.2550060e-02,  5.2611963e-03,  2.1065138e-03,
     +  8.0593143e-04,  2.9184248e-04,  9.7017703e-05,  2.6898336e-05,
     +  3.6156396e-06, -2.9159064e-06, -3.9033580e-06, -3.3254066e-06,
     + -2.4951526e-06, -1.7844325e-06,
     +  9.8495728e+03,  5.2791156e+03,  2.7884128e+03,  1.4450862e+03,
     +  7.2993677e+02,  3.5551632e+02,  1.6378052e+02,  6.8561394e+01,
     +  2.3361559e+01,  3.4161200e+00, -4.2457215e+00, -6.2666907e+00,
     + -5.9443149e+00, -4.8413962e+00, -3.6446609e+00, -2.6147979e+00,
     + -1.8171533e+00, -1.2363432e+00, -8.3057275e-01, -5.5548954e-01,
     + -3.7318605e-01, -2.5437991e-01, -1.7779516e-01, -1.2861332e-01,
     + -9.6849220e-02, -7.5957949e-02, -6.1757651e-02, -5.1639814e-02,
     + -4.4017321e-02, -3.7948317e-02, -3.2888084e-02, -2.8529785e-02,
     + -2.4703090e-02, -2.1313890e-02, -1.8307243e-02, -1.5647323e-02,
     + -1.3306795e-02, -6.6368678e-03, -3.1061478e-03, -1.3862262e-03,
     + -5.9715586e-04, -2.5044238e-04, -1.0287983e-04, -4.1577394e-05,
     + -1.6582779e-05, -6.5444776e-06, -2.5581007e-06, -9.9403312e-07,
     + -3.8314834e-07, -1.4717113e-07 /
*
* spline fit
      if (ifirst .eq. 0) then
* spline fit of the vl coefficients
         ind=1
         do ilam=1,5
           call dcopy(50,vl(ind),1,vec,1)
*    evaluate derivative at first point
           der1=(vec(2)-vec(1))/(rr(2)-rr(1))
           call dspline(rr,vec,50,der1,0d0,csplin(1,ilam))
           ind = ind + 50
         enddo
         ifirst = 1
       end if
* r^-6 fit to isotropic part of potential
       c6 = vl(47)*rr(47)**6
* switching function for long-range
       switch_lr=0.5*(tanh(0.5*(r - 22.d0)) + 1.d0)
* determine splined coefficients at R=r
       ind=1
       do ilam=1,5
         call dcopy(50,vl(ind),1,vec,1)
         call dsplint(rr,vec,csplin(1,ilam),50,r,vvx)
* kill anisotropic terms at large R
         vvx = (1.d0 - switch_lr)*vvx
         if (ilam.eq.1) then
* merge with asymptotic form
            vvx = vvx + switch_lr*c6/(r**6)
         endif
         v(ilam)=vvx
         call dcopy(50,vl(ind),1,vec,1)
         ind = ind + 50
       enddo
       call dcopy(4,v(2),1,vvl(1),1)
* convert to hartree
       econv=1.d0/219474.6d0
       call dscal(4,econv,vvl,1)
       vv0 = v(1) * econv
*
       return
       end
