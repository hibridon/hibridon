*  System:  H2-He - re(H2) - Patkowski et al. PES
*
*   Brandon W. Bakr, Daniel G. A. Smith, and Konrad Patkowski,
*   JCP 139, 144305 (2013)
*  
*   written by p. dagdigian
*   current revision date:  17-oct-2014
      include "common/syusr"
      include "common/bausr"
      include "common/ground"
* ------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(2)
      include "common/parpot"
      econv=219474.6d0
      potnam='H2He PES-Patkowski'
      print *, potnam
1     print *, 'R (bohr):'
      read (5, *, end=99) r
      call pot(vv0,r)
      if (r.le.0.d0) goto 99
*     vlm coefficient is returned in atomic units (hartree)
*     convert from atomic units for printout
      write (6, 100) vv0*econv, (econv*vvl(i), i=1,2)
100   format(3(1pe16.8))
      goto 1
99    rr=3.0d0
      dr=0.2d0
      open (unit=12,file='h2he_vlms.txt')
      write(12,109)
109   format(' %R/bohr V0  V2  V4')
      do i=1,250
        call pot(vv0,rr)
        write (12,110) rr,vv0*econv, (econv*vvl(j),j=1,2)
110     format(f7.2,3(1pe16.8))
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
      potnam='H2He PES-Patkowski'
      npot=1
      nterm=1
      lammin(1)=2
      lammax(1)=4
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
* latest revision date:  17-oct-2014
* ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical ljunk, ihomo, csflag, lljunk
      include "common/parbas"
      dimension v(3)
      dimension csplin(36,3)
      dimension rr(36), vl(108),vec(36)
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(2)
*
*  36 values or R
      data rr /
     +  2.500,  2.750,  3.000,  3.250, 
     +  3.500,  3.750,  4.000,  4.250, 
     +  4.500,  4.750,  5.000,  5.250, 
     +  5.500,  5.750,  6.000,  6.250, 
     +  7.000,  7.250,  7.500,  7.750, 
     +  8.000,  8.500,  9.000,  9.500, 
     + 10.000, 11.000, 12.000, 13.000, 
     + 15.000, 20.000, 25.000, 30.000, 
     + 35.000, 40.000, 50.000, 100.000 / 
*  values of the vlam coefficients
      data vl /  
     +  2.9552141e+04,  1.8849590e+04,  1.1838093e+04,  7.3223120e+03,
     +  4.4586961e+03,  2.6689956e+03,  1.5662269e+03,  8.9651800e+02,
     +  4.9613603e+02,  2.6104752e+02,  1.2604064e+02,  5.0753168e+01,
     +  1.0511631e+01, -9.5760069e+00, -1.8375876e+01, -2.1084725e+01,
     + -1.6543708e+01, -1.4169255e+01, -1.1972892e+01, -1.0041558e+01,
     + -8.3919768e+00, -5.8563473e+00, -4.1204590e+00, -2.9390206e+00,
     + -2.1291839e+00, -1.1701846e+00, -6.7952147e-01, -4.1333010e-01,
     + -1.7093593e-01, -2.9447940e-02, -7.6025971e-03, -2.5249221e-03,
     + -9.9625910e-04, -4.4564924e-04, -1.1637271e-04, -1.8089302e-06,
*
     +  6.4532511e+03,  4.4336922e+03,  2.9492343e+03,  1.9116069e+03,
     +  1.2120406e+03,  7.5333038e+02,  4.5931560e+02,  2.7452295e+02,
     +  1.6045048e+02,  9.1253129e+01,  5.0034475e+01,  2.5976351e+01,
     +  1.2276929e+01,  4.7249348e+00,  7.5308767e-01, -1.1807151e+00,
     + -2.1140736e+00, -1.9000589e+00, -1.6454938e+00, -1.3946160e+00,
     + -1.1673959e+00, -8.0480698e-01, -5.5368250e-01, -3.8476939e-01,
     + -2.7151572e-01, -1.4225178e-01, -7.9416907e-02, -4.6799202e-02,
     + -1.8472448e-02, -2.9872962e-03, -7.4856999e-04, -2.4457760e-04,
     + -9.5552349e-05, -4.2468348e-05, -1.1005993e-05, -1.6935928e-07,
*
     +  6.6711478e+02,  3.8229451e+02,  2.2347974e+02,  1.3279904e+02,
     +  7.9891016e+01,  4.8393339e+01,  2.9359352e+01,  1.7749290e+01,
     +  1.0630412e+01,  6.2720037e+00,  3.6200653e+00,  2.0250849e+00,
     +  1.0813651e+00,  5.3604986e-01,  2.3118796e-01,  6.8731005e-02,
     + -5.5417107e-02, -5.3248886e-02, -4.6002425e-02, -3.7405904e-02,
     + -2.9215995e-02, -1.6445093e-02, -8.6292852e-03, -4.3112209e-03,
     + -2.0845644e-03, -4.5307389e-04, -9.0637142e-05, -1.7391758e-05,
     + -3.0476190e-07, -3.2970097e-08,  3.0476191e-09, -2.5047792e-09,
     + -1.0723230e-09, -3.5206523e-10,  4.7303321e-11,  0.0000000e+00 /
*
* spline fit
      if (ifirst .eq. 0) then
* spline fit of the vl coefficients
         ind=1
         do ilam=1,3
           call dcopy(36,vl(ind),1,vec,1)
*    evaluate derivative at first point
           der1=(vec(2)-vec(1))/(rr(2)-rr(1))
           call dspline(rr,vec,36,der1,0d0,csplin(1,ilam))
           ind = ind + 36
         enddo
         ifirst = 1
       end if
* r^-6 fit to at R = 30 bohr for isotropic part of potential
       c6sum = -1.840668e+06
* switching function for long-range
       switch_lr=0.5*(tanh(0.5*(r - 28.d0)) + 1.d0)
* determine splined coefficients at R=r
       ind=1
       do ilam=1,3
         call dcopy(36,vl(ind),1,vec,1)
         call dsplint(rr,vec,csplin(1,ilam),36,r,vvx)
         if (ilam.eq.1) then
* merge with asymptotic form for V0
           vvx = vvx + switch_lr*c6sum/(r**6)
         else
* kill anisotropic terms at large R
           vvx = (1.d0 - switch_lr)*vvx
         endif
         v(ilam)=vvx
         call dcopy(36,vl(ind),1,vec,1)
         ind = ind + 36
       enddo
       call dcopy(2,v(2),1,vvl(1),1)
* convert to hartree
       econv=1.d0/219474.6d0
       call dscal(2,econv,vvl,1)
* isotropic term 
       vv0 = v(1)*econv
*
       return
       end
*===========================eof===============================