* System:  C2H2(X 1Sigma+)+Ar,  ab initio CCSDT(1.1) Surfaces
* Reference: M. Yang, M. H. Alexander, H.-J. Werner, and R. Bemish,
* J. Chem. Phys. 105, 10462 (1996)
      include "common/syusr"
      include "common/bausr"
      include "common/ground"
      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(3)
1      print *, ' r (bohr)'
      rshift=0.5
      xfact=0.8
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
      common /coselb/ ibasty
      potnam='YANG-ALEXANDER-WERNER Ar-C2H2(X) avtz CCSDT(1.1)'
      ibasty=1
      lammin(1)=2
      lammax(1)=6
      mproj(1)=0
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* ----------------------------------------------------------------------
      subroutine pot (vv0, r)
*  subroutine to calculate the r-dependent coefficients in the
*  Ar-C2H2(X) potential of Yang-Alexander-Werner
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*  vv0        contains isotropic term (d00 term in vsum)
*  variable in common block /covvl/
*    vvl:     vector of length 6 to store r-dependence of each term
*             in potential expansion
*    vvl(1-3) expansion coefficients in dl0 (l=2:2:6) of vsum

* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  4-mar-1998
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension xlam1(4),xlam2(4),r0(4),c1(4),c2(4),c3(4),
     :          clr(4),vsum(4),xsum(4),
     :          vap(4), d0(16),aa(25)
      dimension kpvt(8),qraux(7),work(25),rsd(7)

      common /covvl/ vvl(6)
      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
* coefficicients for d0 rotation matrices
* stored (by column) for each of 4 angles and for l=0:2:6
* angles are 0 30 60 90 120 150 180
* angles are 0 30 60 90
      data d0/
     : 1d0, 1d0, 1d0, 1d0, 1d0, 6.25d-1, -1.25d-1, -5.d-1, 1d0,
     : 2.3438d-2, -2.8906d-1, 3.75d-1, 1d0, -3.7402d-1, 3.2324d-1,
     : -3.125d-1 /

* coefficients for expansion of vap
      data xlam1/
     : 1.8421086d0, 1.7281821d0, 1.7410293d0, 1.7921001d0/
       data xlam2/
     : 1.0042094d0, 9.1408595d-1, 7.2900589d-1, 6.8840575d-1/
      data r0 /
     : 7.1923912d0, 6.7583033d0, 6.3473221d0, 6.4016057d0/
      data c1 /
     : 1.0214178d9, 1.9765905d8, 4.3638903d7, 3.5146118d7/
      data c2/
     :  -1.0029934d6, -3.3017587d5, -4.4972923d4, -2.7434825d4/
      data c3/
     : 1.2214575d2, 7.0633665d1, -6.3854740d-1, -1.1203465d0/
      data clr/
     : 4.0225977d7, 2.8540575d7, 9.0216847d6, 5.2012133d6/

* determine A' potentials at angles
      rm3=one/r**3
      rm6=rm3*rm3
      do 100 i=1,4
        vap(i)=c1(i)*exp(-xlam1(i)*r) +
     :        (c2(i)+c3(i)*r)*exp(-xlam2(i)*r)-
     :         half*(tanh(alph*(r-r0(i)))+1)*clr(i)*rm6
100   continue
* solve simultaneous equations for solutions
* first for vsigma
      tol=1.e-10
      call dcopy(16,d0,1,aa,1)
      call dqrank(aa,4,4,4,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,4,4,4,kr,vap,xsum,rsd,kpvt,qraux)
* convert to hartree
      conv=1.d0/219474.6
      call dscal(4,conv,xsum,1)
      vv0=xsum(1)
      call dcopy(3,xsum(2),1,vvl,1)
      end
