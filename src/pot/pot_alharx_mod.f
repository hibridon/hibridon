* System:  AlH(X 1Sig+)+Ar, modified MRCI+Q PES's
* Reference: M. Yang, M. H. Alexander, S. Gregurick, and P.  J.  Dagdigian, J. Chem. Phys. 102, 2413 (1994).

      include "common/syusr"
      include "common/bausr"
      include "common/ground"

      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(6)
1      print *, ' r (bohr)'
      rshift=0.5
      xfact=0.8
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vvl',/,7(1pe16.8))
      goto 1
99    end
      subroutine loapot(iunit,filnam)
* ------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      potnam='Ar-AlH(X) mod MRCI(D)'
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
*  Ar-AlH(X) potential of Yang and Alexander deepened and shifted
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
* latest revision date:  3-dec-1993
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension xlam(7),r0(7),c1(7),c2(7),clr6(7),
     :          clr8(7),clr10(7),vsum(7),xsum(7),
     :          vap(7), d0(49),aa(64)
      dimension kpvt(8),qraux(7),work(55),rsd(7),re(7)

      common /covvl/ vvl(6)
* rshift is factor for shift of atan midpoint in potential correction
* rrshift is inward shift of ab initio potential
* xfact is factor for additional potential correction
      data rshift, rrshift, xfact /0.6d0, 0.58d0, 0.54d0/
      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
* re's for ab initio potentials
      data re /
     :  9.806, 9.220, 7.505, 7.617, 8.711, 8.967, 8.994/
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
      data xlam/
     : 1.2394876,  1.2061023,  1.3793313,  1.2575876,
     : 1.1437756,  1.1876108,  1.2137655/
      data r0 /
     : 6.1835688,  6.0372674, 3.8955235, 3.9095714,
     : 4.0660949,  4.2047630, 4.1888993/
      data c1 /
     :  6.1826435d+7,  2.9897078d+7,  2.2114368d+7,  7.1129735d+6,
     :  6.9895353d+6,  1.2913208d+7,  1.6453041d+7/
      data c2/
     : -7.2109466d+6, -3.5534298d+6, -2.3133198d+6, -6.9821519d+5,
     : -7.1767842d+5, -1.2336281d+6, -1.4895783d+6/
       data clr6/
     :  2.2618187d+7,  5.7280574d+7, 2.8753793d+7,  2.7264556d+7,
     :  3.0635240d+7,  3.3560220d+7, 3.1319811d+7/
       data clr8/
     :  1.0151952d+10, -3.5347815d+9,  5.4294310d+8,  8.0807568d+8,
     :  2.3611635d+8,   7.9237172d+8,  1.6320728d+9/
       data clr10/
     : -1.0781053d+12, -1.3395149d+10,  3.8863790d+10, 8.2234809d+9,
     :  3.7773219d+10,  6.0188786d+10,  6.0973192d+10/

* determine A' potentials at angles
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
        vap(i)=vap(i)+fact*vap(i)
100   continue
* solve simultaneous equations for solutions
* first for vsigma
      eol=1.e-10
      call dcopy(49,d0,1,aa,1)
      call dqrank(aa,7,7,7,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,7,7,7,kr,vap,xsum,rsd,kpvt,qraux)
* convert to hartree
      conv=1.d0/219474.6
      call dscal(7,conv,xsum,1)
      vv0=xsum(1)
      call dcopy(6,xsum(2),1,vvl,1)
      end
