* System:  BH(A 1Pi)+Ar, original ab initio MRCI+Q PES's
* Reference: M. H. Alexander, S. Gregurick, and P. J.  Dagdigian, J. Chem. Phys. 101, 2887 (1994).

      include "common/syusr"
      include "common/bausr"
      include "common/ground"

      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(11)
1      print *, ' r (bohr)'
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vsum',/,7(1pe16.8),/,
     :    '  vdif',/,5e16.8)
      goto 1
99    end
      subroutine loapot(iunit,filnam)
* ------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      potnam='ALEXANDER Ar-BH(A) MRCI(D)'
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
      subroutine pot (vv0, r)
*  subroutine to calculate the r-dependent coefficients in the
*  Ar-BH(a) potentials of alexander, gregurick, and dagdigian
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
* latest revision date:  8-oct-1993
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension xlam1(14),xlam2(14),r0(14),c1(14),c2(14),c3(14),
     :          clr(14),vsum(7),xsum(7),vdif(7),xdif(7),
     :          ddif(7),vap(7),va2p(7),
     :          d0(49),d2(25),aa(64)
      dimension kpvt(8),qraux(7),work(55),rsd(7)

      common /covvl/ vvl(11)
      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
* for distances beyond rmax difference potential is damped
      data rmax /12.5d0/
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
* coefficicients for d4 rotation matrices
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
     : 6.6079628d-1, 6.1339020d-1, 4.9331378d-1, 6.0330636d-1,
     : 5.4841370d-1, 6.3310444d-1, 6.5920197d-1, 6.5710215d-1,
     : 6.5020685d-1, 6.6040433d-1, 6.9876427d-1, 6.1075198d-1,
     : 6.5121729d-1, 6.5565742d-1/
      data xlam2/
     : 2.1378828, 2.0831808, 1.2183893, 1.5571739, 1.7410957, 1.7715806,
     : 1.6946412, 2.1469134, 2.1610867, 2.2102108, 1.8358129, 1.7795761,
     : 1.7687715, 1.6998204/
      data r0 /
     : 8.3536798, 8.5394774, 7.4104360, 7.5517302, 8.9606031, 8.1646702,
     : 8.2750767, 8.3980686, 8.1695398, 7.6278149, 7.0703293, 8.0080266,
     : 8.0937206, 8.3108120/
      data c1 /
     : -2.9152420d+4, -1.4801623d+4, -3.1303784d+3, -1.2879022d+4,
     : -8.6575067d+3, -1.8932316d+4, -2.4132743d+4, -2.7517653d+4,
     : -2.2474258d+4, -1.9666274d+4, -2.6803662d+4, -1.3540910d+4,
     : -2.0351188d+4, -2.2895274d+4/
      data c2/
     : -9.9238882d+8, -4.9095533d+8, 5.3604750d+6, -1.4612413d+6,
     : -1.6923298d+7, -3.6690024d+7, -1.8711557d+7, -1.0828427d+9,
     : -6.7100461d+8, -2.1888289d+8, -4.3982442d+6, -1.8970448d+7,
     : -3.8680766d+7, -1.9938947d+7/
      data c3/
     : 3.0713289d+8, 1.5626409d+8, -6.0121820d+5, 2.5856655d+6,
     : 7.9446601d+6, 1.2388285d+7, 7.7468654d+6, 3.3019130d+8,
     : 2.0896274d+8, 7.5410352d+7, 2.9777318d+6, 5.7588152d+6,
     : 1.1974545d+7, 8.1396652d+6/
       data clr/
     : -6.2696519d+6, -2.6961564d+6, -2.2810611d+5, 1.0231110d+6,
     : -5.3427602d+6, -3.3704757d+6, -3.7217315d+6, -5.9330264d+6,
     : -5.2533655d+6, -4.2234358d+6, -4.0029865d+6, -5.5887863d+6,
     : -3.1555615d+6, -3.7919772d+6/

* determine A' and A" potentials at angles
      rm3=one/r**3
      rm6=rm3*rm3
      do 100 i=1,7
        vap(i)=c1(i)*exp(-xlam1(i)*r) +
     :        (c2(i)+c3(i)*r)*exp(-xlam2(i)*r)-
     :         half*(tanh(alph*(r-r0(i)))+1)*clr(i)*rm6
        j=i+7
        va2p(i)=c1(j)*exp(-xlam1(j)*r) +
     :        (c2(j)+c3(j)*r)*exp(-xlam2(j)*r)-
     :         half*(tanh(alph*(r-r0(j)))+1)*clr(j)*rm6
        vsum(i)=half*(va2p(i)+vap(i))
* don't compute vdif for colinear geometries
        if (i.ne.1 .and. i.ne.7) then
          vdif(i-1)=half*(-vap(i)+va2p(i))
* for long range damp out difference potential

          if (r .gt. rmax) then
            damp=-half*(tanh(3.d0*(r-rmax))-one)
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
