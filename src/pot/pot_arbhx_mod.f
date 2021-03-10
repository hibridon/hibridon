* System:  BH(X 1Sigma+)+Ar, modified ab initio MRCI+Q PES's
* Reference: M. H. Alexander, S. Gregurick, and P. J.  Dagdigian, J. Chem. Phys. 101, 2887 (1994).

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
      potnam='ALEXANDER Ar-BH(X) avtz mod MRCI(D)'
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
*  Ar-BH(X) potential of Alexander deepened and shifted
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
      dimension xlam1(7),xlam2(7),r0(7),c1(7),c2(7),c3(7),
     :          clr(7),vsum(7),xsum(7),
     :          vap(7), d0(49),aa(64)
      dimension kpvt(8),qraux(7),work(55),rsd(7),re(7)

      common /covvl/ vvl(6)
* rshift is factor for shift of atan midpoint in potential correction
* rrshift is inward shift of ab initio potential
* xfact is factor for additional potential correction
      common /cosysr/ xjunk(2),rshift,xfact,rrshift

      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
* re's for ab initio potentials
      data re /
     :  9.3535, 8.8943544, 7.2427491, 7.2227887, 8.4517169, 8.9341612, 
     :  8.9573738/
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
     : 6.3485064d-1, 6.5829143d-1, 6.7460519d-1,
     : 6.4808585d-1, 6.5243750d-1, 6.2750309d-1, 6.4130420d-1/
       data xlam2/
     : 2.0693477d+0, 2.0707777d+0, 2.2009408d+0,
     : 2.0005953d+0, 1.8136442d+0, 1.8589631d+0, 1.8431337d+0/
      data r0 /
     : 8.9862670d+0, 8.5031660d+0, 7.8453976d+0,
     : 7.9956398d+0, 8.1732049d+0, 8.6532697d+0, 8.4861143d+0/
      data c1 /
     : -1.5005962d+4, -1.6303604d+4, -1.6382683d+4,
     : -1.1387341d+4, -1.3985392d+4, -1.3153971d+4,
     : -1.6667068d+4/
      data c2/
     : -1.2993133d+9, -6.2429475d+8, -2.4391523d+8,
     : -4.2656866d+7, -5.0720574d+7, -1.0343820d+8,
     : -9.0233352d+7/
      data c3/
     : 3.7947750d+8, 1.9297270d+8, 8.3020018d+7,
     : 1.5799605d+7, 1.6801537d+7, 3.5236561d+7, 3.4074165d+7/
      data clr/
     : -1.9144423d+6, -1.4559694d+6, -2.7908688d+6,
     : -1.8366535d+6, -1.6232663d+6, -1.6842528d+6,
     : -1.7503310d+6/

* determine A' potentials at angles
* first shift radius out
      r=rr+rrshift
* xfact is factor for additional potential correction
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
