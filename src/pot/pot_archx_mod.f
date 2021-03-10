* System:  CH(X 2Pi)+Ar, modified ab initio CEPA PES's
* Reference: M. H. Alexander, S. Gregurick, P. J. Dagdigian, G. W. Lemire, M. J. McQuaid, and R. C. Sausa, J. Chem. Phys. 101, 4547 (1994).

      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(11)
1      print *, ' r (bohr)'
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vsum',/,7(1pe16.8),/,
     :    '  vdif',/,5e16.8)
      goto 1
99    end
      include "common/syusr"
      include "common/bausr"
      include "common/ground"

      subroutine loapot(iunit,filnam)
* ------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      potnam='CEPA Ar-CH(X) MOD'
      lammin(1)=1
      lammax(1)=6
      lammin(2)=2
      lammax(2)=6
      mproj(1)=0
      mproj(2)=2
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* ----------------------------------------------------------------------
      subroutine pot (vv0, rr)
*  subroutine to calculate the r-dependent coefficients in the
*  Ar-CH(X) potentials of alexander deepened and shifted
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
*    vvl(7-11) expansion coefficients in dl2 (l=2:6) of vdif

* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  10-dec-1993
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension xlam1(14),xlam2(14),r0(14),c1(14),c2(14),c3(14),
     :          clr(14),vsum(7),xsum(7),vdif(7),xdif(7),
     :          ddif(7),vap(7),va2p(7),
     :          d0(49),d2(25),aa(64)
      dimension kpvt(8),qraux(7),work(55),rsd(7), re(14)

      common /covvl/ vvl(11)
      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
* for distances beyond rmax difference potential is damped
      data rmax /20d0/
* re's for ab initio potentials
      data re /
     : 8.1258, 7.9320, 7.5123, 7.5725, 7.7333, 7.7505, 7.7548,
     : 8.1258, 7.8395, 6.5726, 5.8221, 7.1128, 7.6277, 7.7548/

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
* coefficicients for d2 rotation matrices
* stored (by column) for each of 5 angles and for l=2:6
* angles are 30 60 90 120 150
      data d2/
     : 1.5309311d-1,  4.5927932d-1,  6.1237244d-1,  4.5927934d-1,
     : 1.5309312d-1,
     : 2.9646353d-1,  5.1348990d-1,  1.8279042d-8, -5.1348989d-1,
     : -2.9646355d-1,
     : 4.1999000d-1,  2.2234766d-1, -3.9528471d-1,  2.2234762d-1,
     : 4.1999002d-1,
     : 4.9023048d-1, -1.6982081d-1, -2.4180900d-8,  1.6982085d-1,
     : -4.9023049d-1,
     : 4.8532921d-1, -3.4523418d-1,  3.2021721d-1, -3.4523418d-1,
     : 4.8532920d-1/

* coefficients for expansion of vap(1st 7 entries) and
* for va2p (entries 8:14)
      data xlam1/
     : 6.5305582d-1, 6.3486156d-1, 6.1971171d-1, 6.1599334d-1,
     : 6.2306111d-1, 6.2321511d-1, 6.2503203d-1, 6.1582705d-1,
     : 6.2861111d-1, 6.4423571d-1, 6.4712894d-1, 6.2922170d-1,
     : 6.2906088d-1, 6.2503203d-1/
      data xlam2/
     : 2.1193482, 2.1565991, 2.1196403, 2.0223474,
     : 2.0390668, 2.1210251, 2.1244199, 2.2240155,
     : 2.2203704, 2.5209673, 2.9092415, 2.1534760,
     : 2.0986297, 2.1244199/
      data r0 /
     : 8.3669084, 8.4279840, 8.1839418, 8.2026676,
     : 8.2430145, 8.3527381, 8.4232322, 8.9159275,
     : 8.5629630, 7.8113484, 7.6480887, 8.0069730,
     : 8.1803877, 8.4232322/
      data c1 /
     : -2.1170875d+4, -1.4770608d+4, -9.9753262d+3, -9.1138991d+3,
     : -1.0676544d+4, -1.1369648d+4, -1.1895944d+4, -1.4191307d+4,
     : -1.3652552d+4, -1.2492420d+4, -1.1642870d+4, -1.0092237d+4,
     : -1.1770926d+4, -1.1895958d+4/
      data c2/
     : -7.1836205d+8, -5.6387151d+8, -1.3028309d+8, -5.7650414d+7,
     : -1.0850939d+8, -2.5825162d+8, -2.5935498d+8, -1.8852582d+9,
     : -8.4779037d+8, -8.4606365d+8, -2.4010801d+9, -1.4551906d+8,
     : -1.8620290d+8, -2.5935552d+8/
      data c3/
     : 2.4496656d+8, 1.8912105d+8, 5.3778569d+7, 2.6006001d+7,
     : 4.1825791d+7, 8.8327696d+7, 9.2080766d+7, 5.4222063d+8,
     : 2.6461100d+8, 2.6495857d+8, 6.8641284d+8, 4.6197052d+7,
     : 6.5508436d+7, 9.2080878d+7/
       data clr/
     : -4.5483260d+6, -3.2510755d+6, -2.4384527d+6, -2.2434577d+6,
     : -2.5953016d+6, -2.5612100d+6, -2.6453840d+6, -4.6305257d+6,
     : -3.4706581d+6, -3.5920008d+6, -3.6103524d+6, -2.6343987d+6,
     : -2.8206101d+6, -2.6455254d+6/

* rshift is factor for shift of atan midpoint in potential correction
* rrshift is inward shift of ab initio potential
* xfact is factor for additional potential correction
* default 0.5,0.4,0.6667
      rshift=0.37
      rrshift=0.32
      xfact=0.4

* determine A' and A" potentials at angles
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
        j=i+7
        va2p(i)=c1(j)*exp(-xlam1(j)*r) +
     :        (c2(j)+c3(j)*r)*exp(-xlam2(j)*r)-
     :         half*(tanh(alph*(r-r0(j)))+1)*clr(j)*rm6
* determine switching function
        fact=xfact*0.5*(tanh(alph*(r-re(j)+rshift))+1)
        va2p(i)=va2p(i)+fact*va2p(i)
        vsum(i)=half*(va2p(i)+vap(i))
* don't compute vdif for colinear geometries
        if (i.ne.1 .and. i.ne.7) then
          vdif(i-1)=half*(-vap(i)+va2p(i))
* for long range damp out difference potential

          if (r .gt. rmax) then
            damp=-half*(tanh(1.d0*(r-rmax))-one)
            vdif(i-1)=vdif(i-1)*damp
          endif
        endif
100   continue
* solve simultaneous equations for solutions
* first for vsigma
      tol=1.e-10
      call dcopy(49,d0,1,aa,1)
      call dqrank(aa,7,7,7,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,7,7,7,kr,vsum,xsum,rsd,kpvt,qraux)
* convert to hartree
      conv=1.d0/219474.6
      call dscal(7,conv,xsum,1)
      vv0=xsum(1)
      call dcopy(6,xsum(2),1,vvl,1)
      call dcopy(25,d2,1,aa,1)
      call dqrank(aa,5,5,5,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,5,5,5,kr,vdif,xdif,rsd,kpvt,qraux)
      call dscal(5,conv,xdif,1)
      call dcopy(5,xdif,1,vvl(7),1)
      end
