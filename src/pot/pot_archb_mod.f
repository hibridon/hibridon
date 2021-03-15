* System:  CH(B 2Sigma-)+Ar,  modified ab initio MRCI+Q PES's
* Reference:  M. H. Alexander, S. Gregurick, P. J. Dagdigian, G. W. Lemire, 
* M. J. McQuaid, and R. C. Sausa, J. Chem. Phys. 101, 4547 (1994).

      include "common/syusr"
      include "common/bausr"
      include "common/ground"
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(6)
1      print *, ' r (bohr)'
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vsum',/,7(1pe16.8))
      goto 1
99    end
      subroutine loapot(iunit,filnam)
* ------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      potnam='ALEXANDER Ar-CH(B) MRCI(D) mod'
      lammin(1)=1
      lammax(1)=6
      mproj(1)=0
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* ----------------------------------------------------------------------
      subroutine pot (vv0, rr)
*  subroutine to calculate the r-dependent coefficients in the
*  modified Ar-CH(B) potential of Alexander
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*  vv0        contains isotropic term (d00 term in vsum)
*  variable in common block /covvl/
*    vvl:     vector of length 6 to store r-dependence of each term
*             in potential expansion
*    vvl(1-6) expansion coefficients in dl0 (l=1:6) of vsum

* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  19-jan-1993
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension xlam1(7),xlam2(7),r0(7),c1(7),c2(7),c3(7),
     :          clr(7),vsum(7),xsum(7),
     :          vap(7), d0(49),aa(64)
      dimension kpvt(8),qraux(7),work(55),rsd(7),re(7)

      common /covvl/ vvl(6)
      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
* re's for ab initio potentials
      data re /
     :  7.3327, 7.2928, 7.4059, 7.6015, 7.3922, 7.1451, 7.1272/
* coefficicients for d0 rotation matrices
* stored (by column) for each of 7 angles and for l=0:6
* angles are 0 30 60 90 120 150 180
      data d0/
     : 1d0,  1d0,  1d0,  1d0,  1d0,  1d0,  1d0,
     : 1d0,  8.6602541d-1,  5.0000001d-1,  1.3349125d-8, -4.9999998d-1,
     : -8.6602539d-1, -1d0,
     : 1d0,  6.2500001d-1, -1.2499999d-1, -5.0000000d-1, -1.2500002d-1,
     : 6.2499997d-1,  1d0,
     : 1d0,  3.2475954d-1, -4.3750000d-1, -2.0023687d-8,  4.3750001d-1,
     : -3.2475948d-1, -1d0,
     : 1d0,  2.3437511d-2, -2.8906251d-1,  3.7500000d-1, -2.8906248d-1,
     : 2.3437446d-2,  1d0,
     : 1d0, -2.2327216d-1,  8.9843733d-2,  2.5029609d-8, -8.9843784d-2,
     : 2.2327222d-1, -1d0,
     : 1d0, -3.7402343d-1,  3.2324218d-1, -3.1250000d-1,  3.2324220d-1,
     : -3.7402346d-1,  1d0/

* coefficients for expansion of vap
      data xlam1/
     : 6.4953405d-1, 6.3347935d-1, 6.4288911d-1, 6.2457456d-1,
     : 5.6814151d-1, 7.0037076d-1, 6.2418782d-1/
       data xlam2/
     : 2.2644513, 2.3738468, 2.0359748, 2.2032478, 2.1282801, 1.7973112,
     : 2.0914083/
      data r0 /
     : 8.4233257, 8.4862349, 7.3262907, 7.1690938, 8.3482171, 7.3940358,
     : 9.2886423/
      data c1 /
     : -2.3995201d+4, -1.6447725d+4, -1.3977187d+4, -6.2142578d+3,
     : -7.1355996d+3, -2.5122153d+4, -1.1358545d+4/
      data c2/
     : -7.5157914d+8, -1.6919849d+9, -7.3484587d+7, -2.3980916d+8,
     : -1.3385108d+8, -7.4074586d+6, -8.3423887d+7/
      data c3/
     : 2.5752582d+8, 4.6371148d+8, 3.1794198d+7, 8.0564858d+7,
     : 4.6937646d+7, 4.9938662d+6, 3.0799479d+7/
      data clr/
     : -6.6386126d+6, -4.1506859d+6, -2.6123837d+6, 6.3007597d+6,
     : -1.5735142d+6, -2.2048239d+6, -2.9746218d+6/

* determine A' potentials at angles
* rshift is factor for shift of atan midpoint in potential correction
* rrshift is inward shift of ab initio potential
* xfact is factor for additional potential correction
* default 0.5,0.4,0.6667
      rshift=0.5
      rrshift=0.4
      xfact=0.6667
      rshift=0.37
      rrshift=0.32
      xfact=0.4
* first shift radius out
      r=rr+rrshift

      rm3=one/r**3
      rm6=rm3*rm3
      do 100 i=1,7
        vap(i)=c1(i)*exp(-xlam1(i)*r) +
     :        (c2(i)+c3(i)*r)*exp(-xlam2(i)*r)-
     :         half*(tanh(alph*(r-r0(i)))+1)*clr(i)*rm6
* determine switching function
        fact=xfact*0.5*(tanh(alph*(r-re(i)+rshift))+1)
        vap(i)=vap(i)+fact*vap(i)

100   continue
* solve simultaneous equations for solutions
* first for vsigma
      tol=1.e-10
      call dcopy(49,d0,1,aa,1)
      call dqrank(aa,7,7,7,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,7,7,7,kr,vap,xsum,rsd,kpvt,qraux)
* convert to hartree
      conv=1.d0/219474.6
      call dscal(7,conv,xsum,1)
      vv0=xsum(1)
      call dcopy(6,xsum(2),1,vvl,1)
      end
