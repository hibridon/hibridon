*  system:  H2O-H, PES
*  written by P. Dagdigian, Feb 2013
*
*  H2O defined to lie on the xz plane, with the origin at center of mass
*  x axis is the C2 symmetry axis of the CH2 molecule
*  z axis is the a inertial axis of the molecule (perpendicular to C2 axis)
*  when theta = 90, phi = 0, He is on C side of molecule
*  phi = 0 has all 4 atoms coplanar
*
*  the PES is fitted with 12 angular terms

      include "common/syusr"
      include "common/ground"
      include "common/bausr"
* ------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(11)
      include "common/parpot"
      s4pi = sqrt ( 4.d0 * acos(-1.d0) )
      potnam='H2O-H CCSD(T) PES-12 coeffs'
      print *, potnam
1     print *, 'R (bohr):'
      read (5, *, end=99) r
      call pot(vv0,r)
*     vlm coefficient is returned in atomic units (hartree)
*     convert from atomic units for printout
      econv=219474.6d0
      write (6, 100) vv0*econv*s4pi, (econv*vvl(i), i=1,11)
100   format(' v(lam,0):',4(1pe16.8)/' v(lam,1):',3(1pe16.8)/
     :    ' v(lam,2):',3(1pe16.8)/' v(lam,3):',2(1pe16.8))
99    goto 1
      end
* ------------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
      implicit double precision (a-h,o-z)
      character*(*) filnam
      include "common/parbas"
      include "common/parpot" 
      common /conlam/ nlam, nlammx, lamnum(2)
      common /cosysi/ nscode, isicod, nterm
      potnam='H2O-H CCSD(T) PES-12 coeffs'
      ibasty = 16
*
      nterm = 4
      do i=1,nterm
        mproj(i)=i-1
      enddo
      lammin(1) = 2
      do i=2,nterm
        lammin(i)=i-1
      enddo
      lammax(1) = 6
      lammax(2) = 5
      lammax(3) = 6
      lammax(4) = 5
* 
      ipotsy = 2
*
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
*  calculate total number of anisotropic terms
      nlam = 0
      do i=1,nterm
        lmin = lammin(i)
        lmax = lammax(i)
        lamnum(i)=(lammax(i)-lammin(i))/2+1
        do lb = lmin, lmax, ipotsy
          nlam = nlam + 1
        end do
      end do
      nlammx = nlam
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
*    vvl:        vector of length 11 to store r-dependence of each term
*                in potential expansion
*    vvl(1-3):   expansion coefficients in Yl0 (l=2:6:2) of v(lam,0)
*    vvl(4-6):   expansion coefficients in [Yl1 - Y(l,-1)] (l=1:6:2) of v(lam,1)
*    vvl(7-9):   expansion coefficients in [Yl2 + Y(l,-2)] (l=2:6:2) of v(lam,2)
*    vvl(10-11): expansion coefficients in [Yl3 - Y(l,-3)] (l=3:5:2) of v(lam,3)
*  variable in common block /coloapot/
*    s4pi:       normalization factor for isotropic potential
*
*  uses linear least squares routines from lapack
*
* author:  paul dagdigian
* latest revision date:  feb-5-2013
* ----------------------------------------------------------------------
* this pes is a fit to 19 values of R and 190 orientations, namely
* R = [3:0.5:10 11 12 13 15 20]
* theta=[0:10:90] and phi=[0:10:180]
*
* in this version, 12 terms are employed in the angular fit
      implicit double precision (a-h,o-z)
      logical ljunk, ihomo, csflag, lljunk
      include "common/parbas"
      dimension v(12)
      dimension csplin(20,12)
      dimension rr(20), vl(240),vec(12)
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(11)
* hyperbolic tangent scaling factor
      data alph /1.2d0/
      data rmax /13d0/
* sqrt(4*pi)
      s4pi = 3.544907701811032d0
* here are the values of R
      data rr/3d0,3.5d0,4d0,4.5d0,5d0,5.5d0,6d0,6.5d0,7d0,7.5d0,
     +  8d0,8.5d0,9d0,9.5d0,10d0,11d0,12d0,13d0,15d0,20d0/
* here are column ordered legendre expansion coefficients (12 total) for the
* potential at each of 19 values of R (228 values)
      data vl/
     + 2.3909061e+04, 1.0488210e+04, 4.3975353e+03, 1.6153655e+03,
     + 4.2192764e+02, -2.8526061e+01, -1.5814277e+02, -1.6549200e+02,
     + -1.3582592e+02, -1.0162217e+02, -7.3052842e+01, -5.1743632e+01,
     + -3.6611919e+01, -2.6050142e+01, -1.8713251e+01, -9.9908336e+00,
     + -5.5750548e+00, -3.2366736e+00, -1.1444278e+00, 1.3789457e-01,
     + 6.6362121e+03, 2.5008613e+03, 9.7030643e+02, 3.4356879e+02,
     + 9.8354611e+01, 1.4081773e+01, -8.5328913e+00, -1.0792934e+01,
     + -7.9838369e+00, -4.9947125e+00, -2.8849723e+00, -1.5882531e+00,
     + -8.4306108e-01, -4.4207222e-01, -2.3391583e-01, -5.6546254e-02,
     + -1.8629871e-03, 1.6235961e-02, -3.2560430e-03, 5.3806311e-03,
     + -2.1449004e+01, 1.3281706e+02, 5.0059415e+01, 8.1449024e+00,
     + -2.4013939e+00, -2.8694668e+00, -1.5816423e+00, -5.8149798e-01,
     + -7.2256579e-02, 1.1226450e-01, 1.5094426e-01, 1.2712790e-01,
     + 9.9540794e-02, 7.0647165e-02, 5.7524474e-02, 2.2123814e-02,
     + 1.1573667e-02, 8.0778095e-03, 1.1936292e-02, 5.9034256e-04,
     + -7.4535083e+02, -1.0732530e+02, -1.6328494e+01, -2.6562118e+00,
     + -1.6076867e-01, 2.9501996e-01, 3.1678476e-01, 2.5841423e-01,
     + 1.7690564e-01, 1.1088461e-01, 6.7206653e-02, 4.1009302e-02,
     + 2.3938059e-02, 1.6203924e-02, 6.2312047e-03, 4.4413610e-03,
     + 2.2931754e-03, 1.6791146e-03, -4.7371759e-03, 7.2828836e-03,
     + 4.1825573e+03, 1.6149275e+03, 7.3891197e+02, 3.3427163e+02,
     + 1.3546777e+02, 4.4564507e+01, 7.7411501e+00, -4.5775536e+00,
     + -7.0934235e+00, -6.2738120e+00, -4.7337970e+00, -3.3361328e+00,
     + -2.2740091e+00, -1.5483283e+00, -1.0663421e+00, -5.1299253e-01,
     + -2.5216959e-01, -1.3337413e-01, -4.3281305e-02, 4.2043798e-03,
     + 3.6711810e+03, 1.2624845e+03, 5.3688736e+02, 2.2818093e+02,
     + 8.5650002e+01, 2.4612323e+01, 1.8084755e+00, -4.8501486e+00,
     + -5.5303066e+00, -4.4749537e+00, -3.2190975e+00, -2.2043201e+00,
     + -1.4802686e+00, -9.8575392e-01, -6.5502265e-01, -3.1746457e-01,
     + -1.5942549e-01, -8.1059124e-02, -4.1636902e-02, -1.0401855e-03,
     + 6.3407633e+02, 1.7597303e+02, 5.5413827e+01, 1.6858131e+01,
     + 4.1820810e+00, 3.7555475e-01, -5.1980856e-01, -5.6038075e-01,
     + -4.0335380e-01, -2.5731082e-01, -1.5454466e-01, -9.1385869e-02,
     + -5.3121414e-02, -3.0585453e-02, -2.3133471e-02, -6.0910454e-03,
     + -3.6736580e-03, 2.3664151e-03, 5.2395207e-03, -7.1115836e-03,
     + 2.3996891e+03, 9.2796509e+02, 3.8361391e+02, 1.4769119e+02,
     + 4.6817379e+01, 8.0929291e+00, -4.1563335e+00, -6.3733466e+00,
     + -5.4629800e+00, -4.0165508e+00, -2.7905299e+00, -1.9041486e+00,
     + -1.3066512e+00, -9.0433933e-01, -6.3078774e-01, -3.4746836e-01,
     + -2.0443162e-01, -1.1823403e-01, -3.8162715e-02, -1.1193742e-02,
     + 1.8158201e+03, 4.1161052e+02, 1.3149440e+02, 4.8490443e+01,
     + 1.5551756e+01, 2.6117924e+00, -1.5974791e+00, -2.3589378e+00,
     + -2.0063042e+00, -1.4488087e+00, -9.7383677e-01, -6.4017433e-01,
     + -4.1161810e-01, -2.6220177e-01, -1.7469611e-01, -7.1379072e-02,
     + -3.4210220e-02, -1.5480219e-02, -7.0288466e-03, -4.2105372e-03,
     + 8.0719994e+02, 1.4808499e+02, 3.2301119e+01, 7.7819771e+00,
     + 1.5546963e+00, -2.3310697e-02, -3.4464288e-01, -2.9924940e-01,
     + -2.0712705e-01, -1.3230129e-01, -7.8317196e-02, -4.2411065e-02,
     + -2.6276767e-02, -1.5163462e-02, -4.4247122e-03, -3.3399530e-03,
     + -2.0941946e-03, 2.7463038e-03, 2.1380541e-03, -5.2700257e-03,
     + 7.6628269e+02, 2.5447246e+02, 1.0573224e+02, 4.3908966e+01,
     + 1.5639054e+01, 3.7502811e+00, -5.6016642e-01, -1.7523604e+00,
     + -1.6785974e+00, -1.3040443e+00, -9.3164477e-01, -6.4732787e-01,
     + -4.4198730e-01, -3.0129786e-01, -2.1629647e-01, -1.0052449e-01,
     + -5.4491122e-02, -3.2099542e-02, -7.7142701e-03, -1.1508951e-03,
     + 1.0195569e+03, 1.8567257e+02, 4.3024372e+01, 1.2167001e+01,
     + 3.0140835e+00, 7.5768360e-02, -6.6962185e-01, -7.1672289e-01,
     + -5.4560859e-01, -3.6932183e-01, -2.3432134e-01, -1.4072135e-01,
     + -9.0008437e-02, -5.5624135e-02, -3.0880934e-02, -1.3953582e-02,
     + -7.1605299e-03, -3.6538400e-03, 1.3037250e-03, -1.8180218e-02 /
* spline fit
      if (ifirst .eq. 0) then
* spline fit of the vl coefficients
         ind=1
         do ilam=1,12
           call dcopy(20,vl(ind),1,vec,1)
*    evaluate derivative at first point
           der1=(vec(2)-vec(1))/(rr(2)-rr(1))
           call dspline(rr,vec,20,der1,0d0,csplin(1,ilam))
           ind = ind + 20
         enddo
         ifirst = 1
       end if
* r^-6 fit to isotropic part of potential
       c6sum = -1.9457075e+07
* switching function for long-range
       switch_lr=0.5*(tanh(0.5*(r - 25.d0)) + 1.d0)
* determine splined coefficients at R=r
       ind=1
       do ilam=1,12
         call dcopy(20,vl(ind),1,vec,1)
         call dsplint(rr,vec,csplin(1,ilam),20,r,vvx)
* kill anisotropic terms at large R
         vvx = (1.d0 - switch_lr)*vvx
         if (ilam.eq.1) then
* merge with asymptotic form
            vvx = vvx + switch_lr*c6sum/(r**6)
         endif
         v(ilam)=vvx
         call dcopy(20,vl(ind),1,vec,1)
         ind = ind + 20
       enddo
       call dcopy(11,v(2),1,vvl(1),1)
* convert to hartree
       econv=1.d0/219474.6d0
       call dscal(11,econv,vvl,1)
* isotropic term - divide by sqrt(4*pi)
       vv0 = v(1)*econv/s4pi
*
       return
       end
