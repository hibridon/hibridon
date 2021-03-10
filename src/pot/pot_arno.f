* System:  NO(X 2Pi)+Ar, original ab initio CEPA PES's
* Reference: M. H. Alexander, J. Chem. Phys. 99, 7725 (1993).


      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(19)
      common /coconv/ econv
      include "common/parpot"
      potnam='ALEXANDER Ar-NO CEPA'
      econv=219474.6d0
1      print *, ' r (bohr)'
      rshift=0.5
      xfact=0.8
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vsum',/,11(1pe16.8),/,
     :    '  vdif',/,9e16.8)
      goto 1
99    end
      include "common/syusr"
      include "common/bausr"
      include "common/ground"
* --------------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      potnam='ALEXANDER Ar-NO CEPA'
      lammin(1)=1
      lammax(1)=10
      lammin(2)=2
      lammax(2)=10
      mproj(1)=0
      mproj(2)=2
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
*  -----------------------------------------------------------------------
      subroutine pot (vv0, r)
*  subroutine to calculate the r-dependent coefficients in the
*  Ar-NO potentials of alexander
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*  vv0      spherically symmetric term in V0
*  variable in common block /covvl/
*    vvl:     vector of length 19 to store r-dependence of each term
*             in potential expansion
*    vvl(1-10) expansion coefficients in dl0 (l=1:10) of vsum
*    vvl(11-19) expansion coefficients in dl2 (l=2:10) of vdif

* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  26-jan-1993
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension xlam1(22),xlam2(22),r0(22),c1(22),c2(22),c3(22),
     :          clr(22),vsum(11),xsum(11),vdif(11),xdif(11),
     :          ddif(11),vap(11),va2p(11),
     :          d0(121),d2(81),aa(121)
      dimension kpvt(11),qraux(11),work(55),rsd(11)

      common /covvl/ vvl(19)
      common /coconv/ econv

      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
* for distances beyond rmax difference potential is damped
      data rmax /9d0/
* coefficicients for d0 rotation matrices
* stored (by column) for each of 11 angles and for l=0:10
      data d0/
     : 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0,
     : 1d0, 9.3969262d-1, 7.6604445d-1, 5.0000001d-1, 1.7364819d-1,
     : 1.3349125d-8, -1.7364816d-1, -4.9999998d-1, -7.6604443d-1,
     : -9.3969261d-1, -1d0,
     : 1d0, 8.2453334d-1, 3.8023614d-1, -1.2499999d-1, -4.5476946d-1,
     : -5.0000000d-1, -4.5476947d-1, -1.2500002d-1, 3.8023610d-1,
     : 8.2453331d-1, 1d0,
     : 1d0, 6.6488474d-1, -2.5233323d-2, -4.3750000d-1, -2.4738195d-1,
     : -2.0023687d-8, 2.4738191d-1, 4.3750001d-1, 2.5233373d-2,
     : -6.6488469d-1, -1d0,
     : 1d0, 4.7497774d-1, -3.1900434d-1, -2.8906251d-1, 2.6590160d-1,
     : 3.7500000d-1, 2.6590163d-1, -2.8906248d-1, -3.1900437d-1,
     : 4.7497767d-1, 1d0,
     : 1d0, 2.7149176d-1, -4.1968205d-1, 8.9843733d-2, 2.8101755d-1,
     : 2.5029609d-8, -2.8101752d-1, -8.9843784d-2, 4.1968205d-1,
     : -2.7149167d-1, -1d0,
     : 1d0, 7.1903012d-2, -3.2357074d-1, 3.2324218d-1, -1.3212132d-1,
     : -3.1250000d-1, -1.3212137d-1, 3.2324220d-1, -3.2357069d-1,
     : 7.1902917d-2, 1d0,
     : 1d0, -1.0722615d-1, -1.0060172d-1, 2.2314455d-1, -2.8347993d-1,
     : -2.9201210d-8, 2.8347991d-1, -2.2314450d-1, 1.0060165d-1,
     : 1.0722624d-1, -1d0,
     : 1d0, -2.5183942d-1, 1.3862678d-1, -7.3638895d-2, 2.3307822d-2,
     : 2.7343750d-1, 2.3307885d-2, -7.3638959d-2, 1.3862685d-1,
     : -2.5183950d-1, 1d0,
     : 1d0, -3.5169654d-1, 2.9001294d-1, -2.6789855d-1, 2.5962717d-1,
     : 3.2851362d-8, -2.5962718d-1, 2.6789857d-1, -2.9001298d-1,
     : 3.5169659d-1, -1d0,
     : 1d0, -4.0126914d-1, 2.9734523d-1, -1.8822863d-1, 6.4682158d-2,
     : -2.4609375d-1, 6.4682090d-2, -1.8822857d-1, 2.9734520d-1,
     : -4.0126916d-1, 1d0/
* coefficicients for d2 rotation matrices
* stored (by column) for each of 9 angles and for l=2:10
      data d2/
     : 7.1633966d-2, 2.5301754d-1, 4.5927932d-1, 5.9390714d-1,
     : 6.1237244d-1, 5.9390715d-1, 4.5927934d-1, 2.5301755d-1,
     : 7.1633976d-2,
     : 1.5051848d-1, 4.3340068d-1, 5.1348990d-1, 2.3060770d-1,
     : 1.8279042d-8, -2.3060767d-1, -5.1348989d-1, -4.3340070d-1,
     : -1.5051850d-1,
     : 2.3957418d-1, 5.0756736d-1, 2.2234766d-1, -3.0244623d-1,
     : -3.9528471d-1, -3.0244626d-1, 2.2234762d-1, 5.0756736d-1,
     : 2.3957421d-1,
     : 3.2835759d-1, 4.3600554d-1, -1.6982081d-1, -2.7746878d-1,
     : -2.4180900d-8, 2.7746875d-1, 1.6982085d-1, -4.3600551d-1,
     : -3.2835762d-1,
     : 4.0592179d-1, 2.3830029d-1, -3.4523418d-1, 1.5131754d-1,
     : 3.2021721d-1, 1.5131759d-1, -3.4523418d-1, 2.3830023d-1,
     : 4.0592182d-1,
     : 4.6231022d-1, -1.3906524d-2, -1.9131360d-1, 2.8490319d-1,
     : 2.8675019d-8, -2.8490317d-1, 1.9131355d-1, 1.3906594d-2,
     : -4.6231024d-1,
     : 4.8973053d-1, -2.2700358d-1, 1.1374297d-1, -3.5240934d-2,
     : -2.7731624d-1, -3.5240995d-2, 1.1374303d-1, -2.2700364d-1,
     : 4.8973054d-1,
     : 4.8345442d-1, -3.2461587d-1, 2.7905800d-1, -2.6334950d-1,
     : -3.2484296d-8, 2.6334950d-1, -2.7905801d-1, 3.2461588d-1,
     : -4.8345440d-1,
     : 4.4236809d-1, -2.7891372d-1, 1.6870458d-1, -5.7117526d-2,
     : 2.4836194d-1, -5.7117459d-2, 1.6870452d-1, -2.7891367d-1,
     : 4.4236805d-1/

* coefficients for expansion of vap(1st ll entries) and
* for va2p (entries 11:22)
      data xlam1/
     : 6.0159279d-1, 6.2856352d-1, 6.5406250d-1, 7.0687500d-1,
     : 6.7834662d-1, 6.5389916d-1, 5.5684845d-1, 5.9871302d-1,
     : 6.4898866d-1, 6.1208333d-1, 6.3173970d-1, 7.3068470d-1,
     : 7.1833124d-1, 6.5282381d-1, 6.0362863d-1, 6.0043267d-1,
     : 6.9676612d-1, 6.2956116d-1, 6.4287569d-1, 7.0687909d-1,
     : 6.9982578d-1, 6.3217800d-1/
      data xlam2/
     : 1.5398153d0, 1.5004927d0, 1.4343750d0, 1.4906250d0,
     : 1.7574606d0, 1.6852388d0, 1.4858927d0, 1.6090658d0,
     : 1.6771081d0, 1.6203125d0, 1.6219209d0, 1.5925521d0,
     : 1.5308013d0, 1.5202854d0, 1.5111278d0, 1.5954466d0,
     : 1.7737654d0, 1.6785045d0, 1.6661772d0, 1.7681239d0,
     : 1.6893319d0, 1.6175197d0/
      data r0 /
     : 7.7586501d0, 7.4530646d0, 7.2187500d0, 7.2046875d0,
     : 7.8355119d0, 7.7892111d0, 6.8886196d0, 7.6953928d0,
     : 7.7282604d0, 7.0145833d0, 6.8757716d0, 7.6438629d0,
     : 7.4636526d0, 7.7898239d0, 7.3688745d0, 7.4120364d0,
     : 7.1274074d0, 7.3670173d0, 7.4408467d0, 7.6121919d0,
     : 7.4384372d0, 6.6168347d0/
      data c1 /
     : -1.1037338d+4, -9.8602740d+3, -1.1625749d+4, -1.7798472d+4,
     : -1.4514572d+4, -9.9102550d+3, -3.3318844d+3, -5.8508550d+3,
     : -1.3123780d+4, -6.7151420d+3, -8.6103775d+3, -4.5712868d+4,
     : -3.4407936d+4, -1.3500519d+4, -6.2852500d+3, -5.8015148d+3,
     : -1.7000840d+4, -8.3306955d+3, -8.7477676d+3, -2.1917320d+4,
     : -2.0452743d+4, -8.5686448d+3/
      data c2/
     : 1.0483913d+8, 7.6337096d+7, 3.5723022d+7, 2.3383558d+7,
     : 3.6913855d+7, 3.1800782d+7, 1.5433575d+7, 3.6784210d+7,
     : 7.1655165d+7, 9.4020778d+7, 1.0896266d+8, 1.1829204d+8,
     : 7.5494915d+7, 4.2486344d+7, 2.3126487d+7, 2.3427037d+7,
     : 3.6559614d+7, 3.1954439d+7, 4.4196415d+7, 9.1816104d+7,
     : 1.1220657d+8, 1.0726562d+8/
      data c3/
     : -1.2515702d+7, -9.3782630d+6, -4.4314740d+6, -2.7473276d+6,
     : -3.4393310d+6, -4.0993982d+6, -2.3040482d+6, -4.8408921d+6,
     : -7.9724458d+6, -1.2064178d+7, -1.4040614d+7, -1.2276643d+7,
     : -8.6107099d+6, -5.2632281d+6, -3.1874715d+6, -3.2842322d+6,
     : -3.7019839d+6, -4.1854779d+6, -5.8830472d+6, -9.5379384d+6,
     : -1.3137186d+7, -1.3876586d+7/
      data clr /
     : -6.6402317d+6, -2.0732930d+6, -3.2933224d+6, -7.2318127d+5,
     : -2.2620458d+6, -1.4571115d+6, -2.6133349d+6, -1.4965277d+6,
     : -3.9534923d+6, -5.0564638d+4, -1.3824359d+6, -7.7687365d+6,
     : -6.8707815d+6, -3.0742673d+6, -1.3156279d+6, -1.0717414d+6,
     : -2.8216755d+6, -2.2791322d+6, -1.1552067d+6, -3.2215939d+6,
     : -2.7258252d+6, -1.3967760d+6/
      data ddif /2.6179097d-3, 1.6518007d-2, 2.0609152d-2,
     :  8.8988168d-3, 2.1167721d-3, 4.4815416d-3, 1.2467578d-2,
     : -1.7149729d-3, 2.8956750d-3, 3.5467902d-3, 1.9690444d-3/
* determine A' and A" potentials at angles
      rm3=one/r**3
      rm6=rm3*rm3
      do 100 i=1,11
        vap(i)=c1(i)*exp(-xlam1(i)*r) +
     :        (c2(i)+c3(i)*r)*exp(-xlam2(i)*r)-
     :         half*(tanh(alph*(r-r0(i)))+1)*clr(i)*rm6
        j=i+11
        va2p(i)=c1(j)*exp(-xlam1(j)*r) +
     :        (c2(j)+c3(j)*r)*exp(-xlam2(j)*r)-
     :         half*(tanh(alph*(r-r0(j)))+1)*clr(j)*rm6
        vsum(i)=half*(va2p(i)+vap(i))
* don't compute vdif for colinear geometries
        if (i.ne.1 .and. i.ne.11) then
          vdif(i-1)=half*(va2p(i)-vap(i))
* for long range damp out difference potential
          if (r .gt. rmax) then
            damp=-half*(tanh(3.d0*(r-9.5d0))-one)
            vdif(i-1)=damp*vdif(i-1)
          endif
        endif
100   continue
* solve simultaneous equations for solutions
      lsum=11
      ldif=9
*      lsum=9
*      ldif=7
      tol=1.d-10
      call dscal(19,zero,vvl,1)
      call dcopy(121,d0,1,aa,1)
      call dqrank(aa,11,11,lsum,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,11,11,lsum,kr,vsum,xsum,rsd,kpvt,qraux)
* remove terms less than 0.2 cm-1 in size
      do 110 i=1,11
        if (abs(xsum(i)) .lt. 0.2d0) xsum(i)=zero
110   continue
      vv0=xsum(1)/econv
      call dcopy(lsum-1,xsum(2),1,vvl,1)
      call dcopy(81,d2,1,aa,1)
      call dqrank(aa,9,9,ldif,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,9,9,ldif,kr,vdif,xdif,rsd,kpvt,qraux)
      do 120 i=1,9
        if (abs(xdif(i)) .lt. 0.2d0) xdif(i)=zero
120   continue
      call dcopy(ldif,xdif,1,vvl(11),1)
* convert potential to hartree
      call dscal(19,one/econv,vvl,1)
      end
