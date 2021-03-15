* System:  NH(A 3Pi)+Ne
* Reference:  G. Kerenskaya, U. Schnupf, M. C. Heaven, A. van der Avoird,
* and G. C. Groenenboom, Phys. Chem. Chem. Phys. 7, 846 (2005)

      include "common/syusr"
      include "common/bausr"
      include "common/ground"
*      subroutine driver
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
      potnam='Ne-NH(A)'
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
*  Ne-NH(A) potentials
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
* Ne(3Pi)-Ne coefficients supplied by Galina Kerenskaya 08/20/02
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

     : 6.5 , 6.5, 6.5 , 6.5 , 6.5 , 6.0 , 6.0,
     : 6.5 , 6.5, 7.5 , 7.5 , 7.5 , 6.5 , 6.0/

c     : 6.258 , 6.372 , 6.732 , 6.983 , 6.709 , 6.262 , 6.124 ,
c     : 6.258 , 6.573 , 7.215 , 7.468 , 7.2   , 6.578 , 6.124 /

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
     : 7.7640034d-1, 7.4814750d-1, 7.1437963d-1, 7.1399004d-1,
     : 7.1179086d-1, 7.3849813d-1, 7.8129395d-1, 7.6372532d-1,
     : 7.1498944d-1, 6.0133836d-1, 5.8340090d-1, 5.9279401d-1,
     : 7.0028940d-1, 7.7665018d-1/
      data xlam2/
     : 2.6500153, 2.7012185, 2.7408298, 2.5885109,
     : 2.7482043, 2.8540335, 2.7938430, 2.6790185,
     : 2.3988975, 2.3207935, 2.2709357, 2.3678556,
     : 2.4224086, 2.8109752/
      data r0 /
     : 10.61016954, 10.34727831, 10.29669670, 10.27874888,
     : 10.27588763, 10.20073722, 9.458221310, 10.60789245,
     : 10.69276918, 6.656171890, 6.640556930, 6.53960961,
     : 10.32096312, 9.45967514/
      data c1 /
     : -1.7913620d+4, -1.1591044d+4, -5.4884446d+3, -4.6542580d+3,
     : -4.5190769d+3, -5.7282950d+3, -8.1129478d+3, -1.7851152d+4,
     : -8.7427099d+3, -1.7379968d+3, -1.2833637d+3, -1.4859452d+3,
     : -4.3425383d+3, -7.7846793d+3/
      data c2/
     : -3.7563729d+8, -6.0942897d+8, -5.5066560d+8, -1.5568718d+8,
     : -5.4169305d+8, -8.4416572d+8, -4.3725266d+8, -4.7640989d+8,
     : -4.9222288d+7, -9.6619967d+7, -7.8487097d+7, -1.1839805d+8,
     : -5.6580447d+7, -4.6876253d+8/
      data c3/
     : 1.9164186d+8, 2.2818227d+8, 1.8559002d+8, 6.2042472d+7,
     : 1.7289255d+8, 2.5950424d+8, 1.4765004d+8, 2.2846990d+8,
     : 4.0015914d+7, 3.9741200d+7, 3.2436021d+7, 4.5643685d+7,
     : 2.7496610d+7, 1.5685680d+8/
       data clr/
     : -1.0681784d+6, -1.0214664d+6, -8.1001927d+5, -7.8815822d+5,
     : -7.8551515d+5, -9.3095723d+5, -7.7577423d+5, -1.0739091d+6,
     : -8.1008649d+5, -4.2825085d+5, -5.0719813d+5, -3.8771557d+5,
     : -8.2031915d+5, -7.7732198d+5/

* rshift is factor for shift of atan midpoint in potential correction
* rrshift is inward shift of ab initio potential
* xfact is factor for additional potential correction
      rshift=0.5
      rrshift=0.4
      xfact=0.6667
* unmodified potential
      rshift=0
      rrshift=0
      xfact=0
* determine A' and A" potentials at angles
* first shift radius out
      r=rr+rrshift
      rm3=one/r**3
      rm6=rm3*rm3
      do 100 i=1,7

* A' POTENTIAL
        vap(i)=c1(i)*exp(-xlam1(i)*r) +
     :        (c2(i)+c3(i)*r)*exp(-xlam2(i)*r)-
     :         half*(tanh(alph*(r-r0(i)))+1)*clr(i)*rm6
* determine switching function
        fact=xfact*0.5*(tanh(alph*(r-re(i)+rshift))+1)
        vap(i)=vap(i)+fact*vap(i)
        j=i+7

* A" POTENTIAL
        va2p(i)=c1(j)*exp(-xlam1(j)*r) +
     :        (c2(j)+c3(j)*r)*exp(-xlam2(j)*r)-
     :         half*(tanh(alph*(r-r0(j)))+1)*clr(j)*rm6
* determine switching function
        fact=xfact*0.5*(tanh(alph*(r-re(j)+rshift))+1)
        va2p(i)=va2p(i)+fact*va2p(i)

* SUM POTENTIAL
        vsum(i)=half*(va2p(i)+vap(i))
* don't compute vdif for colinear geometries
* DIF POTENTIAL

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
      vv0=1.7*xsum(1)
      call dcopy(6,xsum(2),1,vvl,1)
      call dcopy(25,d2,1,aa,1)
      call dqrank(aa,5,5,5,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,5,5,5,kr,vdif,xdif,rsd,kpvt,qraux)
      call dscal(5,conv,xdif,1)
      call dcopy(5,xdif,1,vvl(7),1)
      end

