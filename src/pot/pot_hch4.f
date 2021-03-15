*  system:  H-CH4 PES
*  written by P. Dagdigian, Jul 2015
*
*  coordinate definition follows that of Hutson, van der Avoird et al.
*  RCCSD(T)-F12a calculations with AVTZ basis
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
      potnam='H-CH4 RCCSD(T)-F12a'
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
      open (unit=12,file='hch4_vlms.dat')
      write(12,109)
109   format(' %R/bohr V00   V3    V4     V6    V7')
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
      potnam='H-CH4 RCCSD(T)-F12a'
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
* latest revision date:  15-jul-2015
* ----------------------------------------------------------------------
* this pes is a fit to 19 values of R and 190 orientations, namely
* R = [3:0.5:10 11 12 13 15 20]
*
      implicit double precision (a-h,o-z)
      logical ljunk, ihomo, csflag, lljunk
      include "common/parbas"
      dimension v(5)
      dimension csplin(35,5)
      dimension rr(35), vl(175),vec(35)
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(4)
* here are the 35 values of R
      data rr/
     +  3.000,  3.250,  3.500,  3.750,  4.000,
     +  4.250,  4.500,  4.750,  5.000,  5.250,
     +  5.500,  5.750,  6.000,  6.250,  6.500,
     +  6.750,  7.000,  7.250,  7.500,  7.750,
     +  8.000,  8.500,  9.000,  9.500, 10.000,
     + 10.500, 11.000, 11.500, 12.000, 13.000,
     + 14.000, 15.000, 20.000, 25.000, 30.000 /
* here are column ordered angular expansion coefficients (5 total) for the
* potential at each of 35 values of R (175 values)
      data vl/
     +  1.4481944e+04,  9.6893876e+03,  6.5969883e+03,  4.5379829e+03,
     +  3.1261246e+03,  2.1371803e+03,  1.4379350e+03,  9.4459534e+02,
     +  6.0038961e+02,  3.6444050e+02,  2.0637886e+02,  1.0343851e+02,
     +  3.8704190e+01, -1.2910664e-01, -2.1855578e+01, -3.2611340e+01,
     + -3.6585841e+01, -3.6569740e+01, -3.4367150e+01, -3.1113426e+01,
     + -2.7495742e+01, -2.0585795e+01, -1.5009962e+01, -1.0873983e+01,
     + -7.9039661e+00, -5.7949091e+00, -4.2963553e+00, -3.2247515e+00,
     + -2.4491139e+00, -1.4661471e+00, -9.2279815e-01, -6.0521493e-01,
     + -1.0651339e-01, -2.0000000e-02,  0.0000000e+00,
     +  1.0197217e+04,  5.5019767e+03,  3.2554959e+03,  2.1321261e+03,
     +  1.5109034e+03,  1.1147227e+03,  8.2930210e+02,  6.1041824e+02,
     +  4.4057439e+02,  3.1070922e+02,  2.1385098e+02,  1.4354689e+02,
     +  9.3836541e+01,  5.9573440e+01,  3.6510241e+01,  2.1343037e+01,
     +  1.1613281e+01,  5.5370123e+00,  1.8672027e+00, -2.3867824e-01,
     + -1.3505563e+00, -1.9799010e+00, -1.7639126e+00, -1.3611720e+00,
     + -9.8269284e-01, -6.9675798e-01, -5.0081964e-01, -3.5729230e-01,
     + -2.4949618e-01, -1.1988938e-01, -6.5654218e-02, -3.3435635e-02,
     + -1.0989161e-02,  7.9302322e-19,  0.0000000e+00,
     + -5.1458252e+03, -2.1985285e+03, -9.9548010e+02, -5.2306292e+02,
     + -3.3704941e+02, -2.5195516e+02, -1.9857791e+02, -1.5552024e+02,
     + -1.1811960e+02, -8.6659441e+01, -6.1465084e+01, -4.2254217e+01,
     + -2.8172619e+01, -1.8231332e+01, -1.1410807e+01, -6.8535399e+00,
     + -3.8893220e+00, -2.0143675e+00, -8.6813639e-01, -1.9947416e-01,
     +  1.5800252e-01,  3.7428511e-01,  3.4439362e-01,  2.7718094e-01,
     +  2.0808549e-01,  1.5296648e-01,  1.1665077e-01,  8.9328941e-02,
     +  6.2675951e-02,  2.0732359e-02,  5.6867261e-03,  2.5925071e-03,
     +  6.0599393e-03,  1.9949320e-17,  0.0000000e+00,
     +  2.4069580e+03,  9.0207205e+02,  3.1311336e+02,  9.6728979e+01,
     +  2.6076403e+01,  7.4076164e+00,  4.5348779e+00,  4.9967922e+00,
     +  5.2321267e+00,  4.7490287e+00,  3.9014523e+00,  2.9952338e+00,
     +  2.1601028e+00,  1.4806564e+00,  9.6840730e-01,  6.0503042e-01,
     +  3.6680011e-01,  2.2006671e-01,  1.2830983e-01,  6.8503158e-02,
     +  2.7578457e-02, -3.5557929e-03, -4.0324711e-04,  2.3518757e-03,
     +  1.7298284e-03,  9.0057491e-04, -1.8841267e-05, -2.1762591e-03,
     + -3.3390192e-03, -1.9664656e-03, -1.1454656e-03,  1.1228692e-03,
     +  5.0172523e-04, -4.3368087e-18,  0.0000000e+00,
     + -2.7707297e+03, -9.9536748e+02, -3.2706183e+02, -8.9225503e+01,
     + -1.4477131e+01,  3.5113845e+00,  4.6039492e+00,  2.2569635e+00,
     +  2.3062971e-01, -7.3824060e-01, -1.0076018e+00, -9.6954537e-01,
     + -7.8189330e-01, -5.6228303e-01, -3.7324043e-01, -2.2131689e-01,
     + -1.1847090e-01, -6.1065190e-02, -2.5700136e-02, -5.6516364e-04,
     +  2.5364421e-02,  4.9333552e-02,  3.6675543e-02,  2.1141356e-02,
     +  9.0306313e-03, -7.2056066e-04, -4.7925516e-03, -5.1083260e-03,
     + -1.5927120e-03,  2.2408680e-03,  2.2772771e-03,  5.2273280e-04,
     + -8.5232098e-04, -5.3247025e-19,  0.0000000e+00 /
*
* spline fit
      if (ifirst .eq. 0) then
* spline fit of the vl coefficients
         ind=1
         do ilam=1,5
           call dcopy(35,vl(ind),1,vec,1)
*    evaluate derivative at first point
           der1=(vec(2)-vec(1))/(rr(2)-rr(1))
           call dspline(rr,vec,35,der1,0d0,csplin(1,ilam))
           ind = ind + 35
         enddo
         ifirst = 1
       end if
* r^-6 fit to isotropic part of potential
       c6 = vl(32)*rr(32)**6
* switching function for long-range
       switch_lr=0.5*(tanh(0.5*(r - 18.d0)) + 1.d0)
* determine splined coefficients at R=r
       ind=1
       do ilam=1,5
         call dcopy(35,vl(ind),1,vec,1)
         call dsplint(rr,vec,csplin(1,ilam),35,r,vvx)
* kill anisotropic terms at large R
         vvx = (1.d0 - switch_lr)*vvx
         if (ilam.eq.1) then
* merge with asymptotic form
            vvx = vvx + switch_lr*c6/(r**6)
         endif
         v(ilam)=vvx
         call dcopy(35,vl(ind),1,vec,1)
         ind = ind + 35
       enddo
       call dcopy(4,v(2),1,vvl(1),1)
* convert to hartree
       econv=1.d0/219474.6d0
       call dscal(4,econv,vvl,1)
       vv0 = v(1) * econv
*
       return
       end
