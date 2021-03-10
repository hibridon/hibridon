* System:  NO(X 2Pi)+He, original ab initio CEPA PES's
* Reference: M. Yang and M. H. Alexander, J. Chem. Phys. 103, 6973 (1995).


      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(11)
      include "common/parpot"
      potnam='YANG-ALEXANDER He-NO(X) CEPA PES'
      print *, potnam
1      print *, ' r (bohr)'
      rshift=0.5
      xfact=0.8
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
      potnam='YANG-ALEXANDER He-NO(X) CEPA PES'
** 1-18-95 Cepa potential using avqz-f basis  for 7 angles **
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
*  subroutine to calculate the r-dependent coefficients
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
* revised for He-NO(X) : 1-20-95 by Moonbong Yang
* 0 degree for He-NO and 180 for He-ON
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension xlam1(14),xlam2(14),r0(14),c1(14),c2(14),c3(14),
     :          clr(14),vsum(7),xsum(7),vdif(7),xdif(7),
     :          ddif(7),vap(7),va2p(7),
     :          d0(49),d2(25),aa(64)
      dimension kpvt(8),qraux(7),work(55),rsd(7),re(14)

      common /covvl/ vvl(11)
      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
* for distances beyond rmax difference potential is damped
      data rmax /13d0/
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
     : 5.044212d-1, 4.756620d-1, 4.140096d-1, 5.928979d-1,
     : 5.273224d-1, 5.319981d-1, 5.506655d-1,
     : 5.044212d-1, 4.646972d-1, 5.158683d-1, 3.119052d-1,
     : 5.116694d-1, 5.114204d-1, 5.506655d-1/


      data xlam2/
     : 1.6465543, 1.5760182, 1.4568913, 1.7132306,
     : 1.6604188, 1.7149652, 1.7704004,
     : 1.6465543, 1.5442198, 1.4678931, 1.3906384,
     : 1.6215820, 1.6905983, 1.7704004/

      data r0 /
     : 6.7251418, 6.5387895, 6.0356751, 5.2609200,
     : 5.7461632, 6.2598366, 6.5188687,
     : 6.7251418, 6.6651158, 6.5354033, 5.2953720,
     : 5.8893572, 6.2760350, 6.5188687/

      data c1 /
     : -8.22429853d+2, -5.73432822d+2, -2.38335791d+2, -6.85028847d+2,
     : -5.48746882d+2, -7.38896951d+2, -9.61949294d+2,
     : -8.22429853d+2, -5.79003281d+2, -9.59133863d+2, -1.02989852d+2,
     : -5.56042223d+2, -6.44684250d+2, -9.61949294d+2/

      data c2/
     : 3.1271372d+7, 1.6217298d+7, 3.9975835d+6, 4.1273603d+6,
     : 6.8194520d+6, 1.8866965d+7, 3.2367286d+7,
     : 3.1271372d+7, 1.2548594d+7, 3.2884061d+6, 2.0199058d+6,
     : 5.6787621d+6, 1.6120788d+7, 3.2367286d+7/

      data c3/
     : -4.55749683d+6, -2.34581367d+6, -6.25702337d+5, -6.72813315d+5,
     : -1.10734903d+6, -2.90399889d+6, -5.03243538d+6,
     : -4.55749683d+6, -1.94049777d+6, -5.52395616d+5, -3.98089552d+5,
     : -9.83566584d+5, -2.61897134d+6, -5.03243538d+6/

       data clr/
     : -2.72349026d+6, -2.73168996d+6, -2.35814548d+6, 2.06640628d+5,
     : -8.38835271d+5, -1.36158327d+6, -1.44934246d+6,
     : -2.72349026d+6, -3.72309666d+6, -3.95883700d+6, -4.88126505d+6,
     : -1.55406297d+6, -1.72042269d+6, -1.44934246d+6/

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
          vdif(i-1)=half*(va2p(i)-vap(i))
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
