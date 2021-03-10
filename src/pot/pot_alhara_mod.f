* System:  AlH(A 1Pi)+Ar, modified MRCI+Q PES's
* Reference: M. Yang, M. H. Alexander, S. Gregurick, and P.  J. Dagdigian, J. Chem. Phys. 102, 2413 (1994).

      include "common/syusr"
      include "common/bausr"
      include "common/ground"

      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(11)
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
      subroutine loapot(iunit,filnam)
* ------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      potnam='Ar-AlH(A) MRCI(D)'         
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
*  Ar-AlH(A) potentials of alexander deepened and shifted
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*  vv0        contains isotropic term (d00 term in vsum)
*  variable in common block /covvl/
*    vvl:     vector of length 11 to store r-dependence of each term
*             in potential expansion
*    vvl(1-6) expansion coefficients in dl0 (l=1:6) of vsum
*    vvl(7-11) expansion coefficients in dl2 (l=2:6) of vdif

* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  14-dec-1993
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
cmb      dimension xlam1(14),xlam2(14),r0(14),c1(14),c2(14),c3(14),
cmb     :          clr(14),vsum(7),xsum(7),vdif(7),xdif(7),
      dimension xlam(14),r0(14),c1(14),c2(14),clr6(14),clr8(14),
     :          clr10(14),vsum(7),xsum(7),vdif(7),xdif(7),
     :          ddif(7),vap(7),va2p(7),
     :          d0(49),d2(25),aa(64)
      dimension kpvt(8),qraux(7),work(55),rsd(7), re(14)

* rshift is factor for shift of atan midpoint in potential correction
* rrshift is inward shift of ab initio potential
* xfact is factor for additional potential correction
      common /cosysr/ xjunk(2),rshift,xfact,rrshift
      common /covvl/ vvl(11)
      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
* for distances beyond rmax difference potential is damped
      data rmax /15d0/
* re's for ab initio potentials
      data re /
     : 9.211d0, 9.200d0, 10.046d0, 10.071d0, 9.223d0,
     : 8.487d0, 8.469d0,
     : 9.211d0, 8.642d0,  7.378d0,  6.922d0, 7.741d0,
     : 8.355d0, 8.469d0/

      data reav /8.673d0/

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
      data xlam/
     : 1.2084119,  1.2171049,  0.9618411,  0.8899939,
     : 0.9854853,  1.1634865,  0.8882226,
     : 1.2084119,  1.2468105,  1.4175226,  1.0921686,
     : 1.15782108, 1.1359800,  0.8882226/
      data r0 /
     : 6.9279653,  6.3683601, 6.4496028, 4.4346027,
     : 5.1059595,  4.0001789, 7.7937500,
     : 6.9279653,  6.1715856, 3.8852187, 6.5498749,
     : 4.1875611,  4.0999413, 7.7937500/

      data c1 /
     :  3.4333203d+7,  2.2372074d+7,  2.4352647d+6,  1.3892381d+6,
     :  1.8470357d+6,  7.4125888d+6,  1.3660039d+6,
     :  3.4333203d+7,  2.6942644d+7,  2.5451916d+7,  1.6803620d+6,
     :  4.7481305d+6,  6.0768235d+6,  1.3660039d+6/
      data c2/
     : -3.3619811d+6, -2.4606848d+6, -2.5137586d+5, -1.4297978d+5,
     : -2.0262505d+5, -7.0897320d+5, -1.8154911d+5,
     : -3.3619811d+6, -3.3765624d+6, -2.4978916d+6, -2.7494124d+5,
     : -5.2262086d+5, -6.4452948d+5, -1.8154911d+5/

       data clr6/
     : -8.5829810d+6,  6.1076712d+7, 7.1526857d+7,  5.6931404d+7,
     :  6.8332013d+7,  4.1111992d+7, 1.4448104d+7,
     : -8.5829810d+6,  5.5509095d+7, 3.4965457d+7,  5.3434510d+7,
     :  3.6489357d+7,  3.4093573d+7, 1.4448104d+7/

       data clr8/
     :  1.0415885d+10, -1.6287410d+9, -5.5648128d+9, -3.1245946d+9,
     :  -3.7900310d+9,  6.0639852d+8, -6.3919742d+9,
     :  1.0415885d+10, -2.0199823d+9,  4.5880136d+8, -4.2480783d+9,
     :  -1.8705610d+8,  1.9531044d+8, -6.3919742d+9/

       data clr10/
     : -2.5401762d+11, -3.8105803d+10,  1.3636159d+11, 4.4722047d+10,
     :  5.1958223d+10,  3.1131737d+10,  3.3070758d+11,
     : -2.5401762d+11, -5.2560961d+10,  4.2776892d+10, 1.0396183d+11,
     :  3.1957749d+10,  3.5838137d+10,  3.3070758d+11/

* determine A' and A" potentials at angles

* first shift radius out
      r=rr+rrshift
      rm2=one/r**2
      rm3=one/r**3
      rm6=rm3*rm3
      rm8=rm6*rm2
      rm10=rm8*rm2

      do 100 i=1,7
        vap(i)=
     :        (c1(i)+c2(i)*r)*exp(-xlam(i)*r)-
     :         half*(tanh(alph*(r-r0(i)))+1)*
     :        (clr6(i)*rm6+clr8(i)*rm8+clr10(i)*rm10)

* determine switching function
        fact=xfact*0.5*(tanh(alph*(r-re(i)+rshift))+1)
* determine amplification factor
        vap(i)=vap(i)+fact*vap(i)
        j=i+7
        va2p(i)=
     :        (c1(j)+c2(j)*r)*exp(-xlam(j)*r)-
     :         half*(tanh(alph*(r-r0(j)))+1)*
     :        (clr6(j)*rm6+clr8(j)*rm8+clr10(j)*rm10)

* determine switching function
        fact=xfact*0.5*(tanh(alph*(r-re(j)+rshift))+1)
        va2p(i)=va2p(i)+fact*va2p(i)
        vsum(i)=half*(va2p(i)+vap(i))
* don't compute vdif for colinear geometries
        if (i.ne.1 .and. i.ne.7) then
          vdif(i-1)=half*(-vap(i)+va2p(i))
* for long range damp out difference potential

          if (r .gt. rmax) then
            damp=-half*(tanh(3.d0*(r-rmax-5))-one)
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
