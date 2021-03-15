* System:  OH(X 2Pi)+Ar, original ab initio RCCSD(T) PES's
* Calculations by R. Tobola and J. Klos
* PES CORRECTED BY J. KLOS FOR CENTER OF MASS OF OH
* BY J. KLOS 2009.06.12
* ORIGINALL CCSD(T) PES WAS CALCULATED WITH RESPECT TO
* ORIGIN LYING 0.86192 Angstrom instead 0.9122 ANgstrom
* FROM HYDROGEN ATOM

*Reference:  Grant Paterson , Sarantos Marinakis , 
* Matthew L. Costen , Kenneth G. McKendrick , 
* Jacek Klos , and Robert Tobola,J. Chem. Phys. 129, 074304 (2008)
*[ERRATUM Grant Paterson et al. J. Chem. Phys. 131, 159901 (2009)]



      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(19)
      include "common/parpot"
      potnam='R. Tobola and J.Klos Ar-OH(X) RCCSDT PES'
      print *, potnam
1      print *, ' r (bohr)'
      rshift=0.5
      xfact=0.8
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vsum',/,11(1pe16.8),/,
     :    '  vdif',/,9(1pe16.8))
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
*  variable in common block /conlam/
*    nlam:      the number of angular coupling terms actually used
*    nlammx:    the maximum number of angular coupling terms
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx
      common /conlam/ nlam, nlammx, lamnum(2)
      potnam='R. Tobola and J.Klos Ar-OH(X) RCCSDT PES'
** ccsdt potential using avqz+bf basis  for 11 angles **
      lammin(1)=1
      lammax(1)=10
      lammin(2)=2
      lammax(2)=10
      lamnum(1)=lammax(1)-lammin(1)+1
      lamnum(2)=lammax(2)-lammin(2)+1
      nlam=lamnum(1)+lamnum(2)
      nlammx=nlam
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
*    vvl:     vector of length 19 to store r-dependence of each term
*             in potential expansion
*    vvl(1-10) expansion coefficients in dl0 (l=1:10) of vsum
*    vvl(11-19) expansion coefficients in dl2 (l=2:10) of vdif

* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  8-oct-1993
* revised for Ar-OH  Theta=0 for Ar-HO
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension xlam1(22),xlam2(22),r0(22),c1(22),c2(22),c3(22),
     :          clr(22),vsum(11),xsum(11),vdif(11),xdif(11),
     :          ddif(11),
     :          d0(121),d2(81),aa(121),thta(11),cthta(11)
      dimension kpvt(11),qraux(11),work(200),rsd(11),re(22)

      common /covvl/ vvl(19)
      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
* for distances beyond rmax difference potential is damped
      data rmax /20d0/
* coefficicients for d0 rotation matrices
* stored (by column) for each of 11 angles and for l=0:10
* angles are  0 20 40 60 80 90 100 120 140 160 180
      data d0/
     :  1.d0,!l=0
     :  1.d0,
     :  1.d0,
     :  1.d0,
     :  1.d0,
     :  1.d0,
     :  1.d0,
     :  1.d0,
     :  1.d0,
     :  1.d0,
     :  1.d0,!l=0
     :                    1.d0, !l=1
     :    0.939692620785908d0,
     :    0.766044443118978d0,
     :                  0.5d0,
     :     0.17364817766693d0,
     :                  1d-16,
     :    -0.17364817766693d0,
     :                 -0.5d0,
     :   -0.766044443118978d0,
     :   -0.939692620785908d0,
     :                   -1d0, !l=1
     :                    1d0, !l=2
     :    0.824533332339233d0,
     :    0.380236133250197d0,
     :               -0.125d0,
     :   -0.454769465589431d0,
     :                 -0.5d0,
     :   -0.454769465589431d0,
     :               -0.125d0,
     :    0.380236133250197d0,
     :    0.824533332339233d0,
     :                    1d0, !l=2
     :                    1d0, !l=3
     :    0.664884732794716d0,
     :  -0.0252333338303835d0,
     :              -0.4375d0,
     :   -0.247381933374901d0,
     :                  1d-16,
     :    0.247381933374901d0,
     :               0.4375d0,
     :   0.0252333338303835d0,
     :   -0.664884732794716d0,
     :                   -1d0, !l=3
     : 1d0,                   !l=4
     :    0.474977735636283d0,
     :   -0.319004346471378d0,
     :           -0.2890625d0,
     :    0.265901610835095d0,
     :                0.375d0,
     :    0.265901610835095d0,
     :           -0.2890625d0,
     :   -0.319004346471378d0,
     :    0.474977735636283d0,
     :                    1d0,!l=4
     :  1.d0,                 !l=5
     :    0.271491745551255d0,
     :   -0.419682045437054d0,
     :           0.08984375d0,
     :    0.281017540988309d0,
     :                  1d-16,
     :   -0.281017540988309d0,
     :          -0.08984375d0,
     :    0.419682045437054d0,
     :   -0.271491745551255d0,
     :                   -1d0,!l=5
     :                   1.d0,!l=6
     :   0.0719030017842305d0,
     :   -0.323570725710931d0,
     :         0.3232421875d0,
     :   -0.132121338573299d0,
     :              -0.3125d0,
     :   -0.132121338573299d0,
     :         0.3232421875d0,
     :   -0.323570725710931d0,
     :   0.0719030017842305d0,
     :                    1d0, !l=6
     :                    1d0, !l=7
     :   -0.107226158692938d0,
     :   -0.100601708629502d0,
     :        0.22314453125d0,
     :   -0.283479918813435d0,
     :                  1d-16,
     :    0.283479918813435d0,
     :       -0.22314453125d0,
     :    0.100601708629502d0,
     :    0.107226158692938d0,
     :                   -1d0,!l=7
     :                    1d0,!l=8
     :   -0.251839432959275d0,
     :    0.138626797752243d0,
     :   -0.073638916015625d0,
     :   0.0233078500507821d0,
     :            0.2734375d0,
     :   0.0233078500507821d0,
     :   -0.073638916015625d0,
     :    0.138626797752243d0,
     :   -0.251839432959275d0,
     :                    1d0,!l=8
     :                    1d0, !l=9
     :   -0.351696543958561d0,
     :    0.290012951832139d0,
     :   -0.267898559570312d0,
     :    0.259627174131175d0,
     :                  1d-16,
     :   -0.259627174131175d0,
     :    0.267898559570312d0,
     :   -0.290012951832139d0,
     :    0.351696543958561d0,
     :                   -1d0, !l=9
     :                    1d0, !l=10
     :   -0.401269139852809d0,
     :    0.297345221371711d0,
     :   -0.188228607177734d0,
     :   0.0646821277096134d0,
     :          -0.24609375d0,
     :   0.0646821277096134d0,
     :   -0.188228607177734d0,
     :    0.297345221371711d0,
     :   -0.401269139852809d0,
     :                    1d0/!l=10
     
* coefficicients for d2 rotation matrices
* stored (by column) for each of 9 angles and for l=2:10
* angles are 20 40 60 80 90 100 120 140 160  
      data d2/
     :   0.0716339671058987d0,!l=2
     :    0.25301753909188d0,
     :    0.459279326771846d0,
     :    0.593907147345913d0,
     :    0.612372435695794d0,
     :    0.593907147345913d0,
     :    0.459279326771846d0,
     :    0.25301753909188d0,
     :    0.0716339671058987d0, !l=2
     :   0.150518479233129d0, !l=3
     :    0.433400687687707d0,
     :    0.513489897661093d0,
     :    0.230607689206516d0,
     :                  1d-16,
     :   -0.230607689206516d0,
     :   -0.513489897661093d0,
     :   -0.433400687687707d0,
     :   -0.150518479233129d0, !l=3
     :    0.239574181663041d0, !l=4
     :    0.507567357301874d0,
     :    0.222347647980589d0,
     :   -0.302446243003737d0,
     :   -0.395284707521047d0,
     :   -0.302446243003737d0,
     :    0.222347647980589d0,
     :    0.507567357301874d0,
     :    0.239574181663041d0, !l=4      
     :    0.328357589465018d0, !l=5
     :    0.436005533330449d0,
     :   -0.169820821244407d0,
     :   -0.277468765109838d0,
     :                  1d-16,
     :    0.277468765109838d0,
     :    0.169820821244407d0,
     :   -0.436005533330449d0,
     :   -0.328357589465018d0, !l=5  
     :    0.405921793012549d0, !l=6
     :    0.238300278063999d0,
     :   -0.345234181079693d0,
     :    0.151317563585355d0,
     :    0.320217211436238d0,
     :    0.151317563585355d0,
     :   -0.345234181079693d0,
     :    0.238300278063999d0,
     :    0.405921793012549d0, !l=6
     :     0.46231022215767d0, !l=7
     :  -0.0139065397386551d0,
     :   -0.191313584921252d0,
     :    0.284903176975073d0,
     :                  1d-16,
     :   -0.284903176975073d0,
     :    0.191313584921252d0,
     :   0.0139065397386551d0,
     :    -0.46231022215767d0, !l=7
     :    0.489730532881397d0, !l=8
     :   -0.227003593559681d0,
     :     0.11374298899392d0,
     :  -0.0352409613338776d0,
     :   -0.277316239832795d0,
     :  -0.0352409613338776d0,
     :     0.11374298899392d0,
     :   -0.227003593559681d0,
     :    0.489730532881397d0, !l=8
     :     0.48345441474116d0, !l=9
     :    -0.32461587132088d0,
     :    0.279058005816313d0,
     :   -0.263349502692564d0,
     :                  1d-16,
     :    0.263349502692564d0,
     :   -0.279058005816313d0,
     :     0.32461587132088d0,
     :    -0.48345441474116d0, !l=9
     :    0.442368089279743d0, !l=10
     :   -0.278913706604953d0,
     :    0.168704562081535d0,
     :   -0.057117496039293d0,
     :     0.24836194310956d0,
     :   -0.057117496039293d0,
     :    0.168704562081535d0,
     :   -0.278913706604953d0,
     :    0.442368089279743d0/!l=10


      thta(1)=0.D0
      thta(2)=20.D0
      thta(3)=40.D0
      thta(4)=60.D0
      thta(5)=80.D0
      thta(6)=90.D0
      thta(7)=100.D0
      thta(8)=120.D0
      thta(9)=140.D0
      thta(10)=160.D0
      thta(11)=180.D0

* determine A' and A" potentials at angles
      do 100 i=1,11
        vsum(i)=VsumPES(r,thta(i))
* don't compute vdif for colinear geometries
        if (i.ne.1 .and. i.ne.11) then
          vdif(i-1)=VdifPES(r,thta(i))
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
       lsum=11
       ldif=9
*      lsum=9
*      ldif=7
      call dscal(19,zero,vvl,1)
      call dcopy(121,d0,1,aa,1)
      call dqrank(aa,11,11,lsum,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,11,11,lsum,kr,vsum,xsum,rsd,kpvt,qraux)
* convert to hartree
      conv=1.d0/219474.6
      call dscal(11,conv,xsum,1)
      vv0=xsum(1)
      call dcopy(lsum-1,xsum(2),1,vvl,1)
      call dcopy(81,d2,1,aa,1)
      call dqrank(aa,9,9,ldif,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,9,9,ldif,kr,vdif,xdif,rsd,kpvt,qraux)
      call dscal(9,conv,xdif,1)
      call dcopy(ldif,xdif,1,vvl(11),1)
      end


      FUNCTION VsumPES(R_GOOD, Theta_GOOD)
      IMPLICIT DOUBLE PRECISION  ( A-H, O-Z )
C   DISTANCE OF H FROM CENTER OF MASS IN ANGSTROM
C   THIS FUNCTION IS TO CORRECT FOR PROPER CENTER 
C   OF MASS OF OH FOR THE FIT THAT WAS CALCULATED
C   IN WRONG CENTER OF MASS ORIGINALLY
C   REVISION 2009.06.12 
      RHX_GOOD=0.912221737264442D0
      RHX_BAD=0.861920D0
      TOAU=0.529177249D0
      RHXDIFF=(RHX_GOOD-RHX_BAD)/TOAU
      R_BAD=DSQRT(R_GOOD**2+RHXDIFF**2-
     .      2.D0*R_GOOD*RHXDIFF*DCOSD(Theta_GOOD))
      IF (Theta_GOOD.EQ.0.D0) THEN
        Theta_BAD=0.D0 
      ELSE
      IF (Theta_GOOD.EQ.180.D0) THEN
        Theta_BAD=180.D0 
      ELSE
      Theta_BAD=DACOSD((R_GOOD**2-R_BAD**2-RHXDIFF**2)/
     .                (2.D0*R_BAD*RHXDIFF))
      ENDIF
      ENDIF
      VsumPES=VsumPES_BadCOM(R_BAD,Theta_BAD)
      RETURN
      END

      FUNCTION VdifPES(R_GOOD, Theta_GOOD)
      IMPLICIT DOUBLE PRECISION  ( A-H, O-Z )
C   DISTANCE OF H FROM CENTER OF MASS IN ANGSTROM
C  THIS FUNCTION IS TO CORRECT FOR PROPER CENTER 
C   OF MASS OF OH FOR THE FIT THAT WAS CALCULATED
C   IN WRONG CENTER OF MASS ORIGINALLY
      RHX_GOOD=0.912221737264442D0
      RHX_BAD=0.861920D0
      TOAU=0.529177249D0
      RHXDIFF=(RHX_GOOD-RHX_BAD)/TOAU
      R_BAD=DSQRT(R_GOOD**2+RHXDIFF**2-
     .      2.D0*R_GOOD*RHXDIFF*DCOSD(Theta_GOOD))
      IF (Theta_GOOD.EQ.0.D0) THEN
        Theta_BAD=0.D0 
      ELSE
      IF (Theta_GOOD.EQ.180.D0) THEN
        Theta_BAD=180.D0
      ELSE
      Theta_BAD=DACOSD((R_GOOD**2-R_BAD**2-RHXDIFF**2)/
     .                (2.D0*R_BAD*RHXDIFF))
      ENDIF
      ENDIF
      VdifPES=VdifPES_BadCOM(R_BAD,Theta_BAD)
      RETURN
      END

   
      FUNCTION VsumPES_BadCOM(R,theta)
C*********************************
C THIS PES WAS CALCULATED IN WRONG CENTER OF MASS!!!
C System: Ar-OH(X2Pi)
C Method:RCCSD(T)
C Basis:aug-cc-pvqz+3s3p2d2f1g
C CP and BSSE corrected
C PES: Vsum=1/2(A''+A')
C Variables:
C R in Bohr
C Theta n Degrees
C Output in cm-1
C Theta=0 for Ar---HO geom.
C********************************
C Authors:
C Dr Jacek Klos
C Department of Chemistry and Biochemistry
C University of Maryland
C College Park, MD 20742
C email:jklos@mail.umd.edu
C 
C Dr Robert Tobola
C Department of Chemistry
C Oakland University
C Rochester, MI 48309
C email:
C**********************************
C Needs link with LAPACK
C*********************************
      implicit double precision(a-h, o-z)
      dimension V0(11)
      dimension T(11)
      pi=dacos(-1.d0)
      call Vlsum(R,V0)
       do j=1,11
       T(j)=PLGNDR((j-1),0,DCOS(theta*pi/180.D0))
       enddo
       s=0.D0
       do i=1,11
       s=s+V0(i)*T(i)
       enddo
       VsumPES_BadCOM=s
       return
       end

      FUNCTION VdifPES_BadCOM(R,theta)
C*********************************
C THIS PES WAS CALCULATED IN WRONG CENTER OF MASS!!!
C System: Ar-OH(X2Pi)
C Method:RCCSD(T)
C Basis:aug-cc-pvqz+3s3p2d2f1g
C CP and BSSE corrected
C PES: Vdif=1/2(A''-A')
C Variables:
C R in Bohr
C Theta n Degrees
C Output in cm-1
C Theta=0 for Ar---HO geom.
C********************************
C Authors:
C Dr Jacek Klos
C Department of Chemistry and Biochemistry
C University of Maryland
C College Park, MD 20742
C email:jklos@mail.umd.edu
C 
C Dr Robert Tobola
C Department of Chemistry
C Oakland University
C Rochester, MI 48309
C email:
C**********************************
C Needs link with LAPACK
C*********************************
      implicit double precision(a-h, o-z)
      dimension V0(9)
      dimension T(9)
      dimension FN(400)
      common/FACTOR/FN
      call setup
      pi=dacos(-1.d0)
      call Vldif(R,V0)
       do j=1,9
       TERM1=FN(j+1-2+1)
       TERM2=FN(j+1+2+1)
       FACT = DEXP(TERM1-TERM2)
       FACT = DSQRT(FACT)
       T(j)=FACT*PLGNDR(j+1,2,DCOS(theta*pi/180.D0))
       enddo
       s=0.D0
       do i=1,9
       s=s+V0(i)*T(i)
       enddo
       VdifPES_BadCOM=s
       return
       end

      subroutine Vlsum(R,V0)
      implicit double precision(a-h, o-z)
      dimension T(11,11),ipvt(11),work(100)
      dimension V0(11),theta(11)
      LWORK=100
      pi=dacos(-1.d0)
      theta(1)=0.D0
      theta(2)=20.D0
      theta(3)=40.D0
      theta(4)=60.D0
      theta(5)=80.D0
      theta(6)=90.D0
      theta(7)=100.D0
      theta(8)=120.D0
      theta(9)=140.D0
      theta(10)=160.D0
      theta(11)=180.D0
      do i=1,11
       do j=1,11
       T(i,j)=PLGNDR((j-1),0,DCOS(theta(i)*pi/180.D0)) 
       enddo
      enddo
      do k=1,11
      V0(k)=0.5D0*(VAPRIMI(R,k)+VABISI(R,k))
      enddo
      call dgesv(11,1,T,11,ipvt,V0,11,info)
c      call dgels('N',11,9,1,T,11,V0,11,WORK,LWORK,info)
      return
      end
      subroutine Vldif(R,V0)
      implicit double precision(a-h, o-z)
      dimension T(9,9),ipvt(9)
      dimension V0(9),theta(9)
      dimension FN(400)
      common/FACTOR/FN
      call setup
      pi=dacos(-1.d0)
      theta(1)=20.D0
      theta(2)=40.D0
      theta(3)=60.D0
      theta(4)=80.D0
      theta(5)=90.D0
      theta(6)=100.D0
      theta(7)=120.D0
      theta(8)=140.D0
      theta(9)=160.D0
      do i=1,9
       do j=1,9
       TERM1=FN(j+1-2+1)
       TERM2=FN(j+1+2+1)
       FACT = DEXP(TERM1-TERM2)
       FACT = DSQRT(FACT)
       T(i,j)=FACT*PLGNDR(j+1,2,DCOS(theta(i)*pi/180.D0))
       enddo
      enddo
      do k=1,9
      V0(k)=0.5D0*(VABISI(R,k+1)-VAPRIMI(R,k+1))
      enddo
      call dgesv(9,1,T,9,ipvt,V0,9,info)

      return
      end
C*******************************************************
      SUBROUTINE SETUP
C ***CALCULATES AND STORES FACTORIALS FOR 3NJ PROGRAMS
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FACL(400)
      COMMON/FACTOR/FACL
C     ***N CAN BE INCREASED OR DECREASED AS DESIRED
C
      FACL(1) = 0.D0
      DO 100 I=2,400
             I1 = I-1
             FACL(I) = FACL(I-1) + DLOG(DFLOAT(I1))
100   CONTINUE
      RETURN
      END
C******************************************************

      FUNCTION PLGNDR(L,M,x)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.)pause
     *'bad arguments in plgndr'
      pmm=1.
      if(m.gt.0) then
        somx2=dsqrt((1.-x)*(1.+x))
        fact=1.
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.
11      continue
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
         do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndr=pll
        endif
      endif
      return
      END


      FUNCTION VAPRIMI(R,N)
CSYSTEM: OH-Ar Aprim state T
CLEVEL:RHF/RCCSD(T)/aug-cc-pvtz+33221
C units R=au E=cm-1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(17),XX1(17),XX2(17),XX3(17),
     * XX4(17),XX5(17),XX6(17),XX7(17),XX8(17),
     * XX9(17),XX10(17),XX11(17)
C theta=0
      DATA XX1/   
     *     -0.5835284219401391D+00,        
     *       38.04952789322221D+00,        
     *      -985.9019652754446D+00,  
     *       11518.96508535197D+00,  
     *      -21506.60006303519D+00,  
     *      -1007794.615872079D+00,  
     *       12184854.96419048D+00,  
     *      -59734675.55243188D+00,  
     *       112230765.7760808D+00,  
     *      0.9298477185884633D+00,  
     *      0.3179101021369883D+00,  
     *      -1.217028123076279D+00,  
     *      -33195530529219.63D+00,  
     *       240524473649379.3D+00,  
     *      -453384045468199.6D+00,  
     *       22784339.37813496D+00,
     *       673550266.2988292D+00/
C theta=20
      DATA XX2/
     *     -0.2467609222528962D+00,              
     *       17.26753563883851D+00,              
     *      -525.4799329248770D+00,  
     *       8974.015625243632D+00,  
     *      -92780.78637041834D+00,  
     *       579553.7916646295D+00,  
     *      -2013450.798753579D+00,  
     *       2891015.672154362D+00,  
     *       523861.2216936092D+00,  
     *       1.047845764413967D+00,  
     *      -2.505191160849684D+00,  
     *      -1.443465811465117D+00,  
     *      -213134237815293.5D+00,  
     *       1573589709217317.D+00,  
     *      -2961141900801557.D+00,  
     *       26015682.89534215D+00,
     *       580544634.5541881D+00/
C theta=40
      DATA XX3/
     *     -0.2241288995485578D+00,          
     *       15.32075097919870D+00,          
     *      -451.2869984552427D+00,  
     *       7339.398703454794D+00,  
     *      -70099.07106983257D+00,  
     *       377099.8265538289D+00,  
     *      -881359.6529619953D+00,  
     *      -726688.3184882962D+00,  
     *       5539392.221177019D+00,  
     *       1.014101448099209D+00,  
     *      -1.768854473835698D+00,  
     *      -1.062112799312645D+00,  
     *      -201765878.6255219D+00,  
     *       1092954390159.944D+00,  
     *      -3195510843870.502D+00,  
     *       30653317.28656752D+00,
     *       402679067.9708959D+00/
C theta=60
      DATA XX4/
     *    -0.2460320887341368D+00,         
     *      15.54610016819009D+00,         
     *     -410.0095744463382D+00, 
     *      5550.346939390212D+00, 
     *     -35597.68624180995D+00, 
     *      428.7837750242383D+00, 
     *      1534790.540634057D+00, 
     *     -9276931.321156578D+00, 
     *      18383761.67349069D+00, 
     *     0.9854644314568246D+00, 
     *    -0.6766692056683580D+00, 
     *    -0.8329984500493898D+00, 
     *     -8718759637.581518D+00, 
     *      176303656271.3300D+00, 
     *     -431637373554.6143D+00, 
     *      32532339.69886902D+00,
     *      265041085.2395844D+00/
C theta=80
      DATA XX5/
     *     -0.1816786054107331D+00,            
     *       11.16319999924390D+00,            
     *      -287.6081870872177D+00, 
     *       3814.947395240347D+00, 
     *      -24029.51864146546D+00, 
     *      -191.2997031438089D+00, 
     *       1013775.173412235D+00, 
     *      -6059977.173475049D+00, 
     *       11878381.56127396D+00, 
     *      0.9967232084828510D+00, 
     *     -0.8524794188392986D+00, 
     *      -1.342123228725824D+00, 
     *      -980212630023.9100D+00, 
     *       6958024088525.354D+00, 
     *      -12160209509583.75D+00, 
     *       32337795.33680038D+00,
     *       204186225.6738213D+00/
C theta=90
      DATA XX6/
     *     -0.1703823198957445D+00,     
     *       10.63619049567595D+00,     
     *      -278.1849758104884D+00,  
     *       3741.762623253233D+00,  
     *      -23863.47000810176D+00,  
     *      -291.0841380767020D+00,  
     *       1028765.068266686D+00,  
     *      -6194827.972336392D+00,  
     *       12210728.90458321D+00,  
     *      0.9833079607952939D+00,  
     *     -0.7157542798286749D+00,  
     *      -1.258812152565030D+00,  
     *      -699887474769.4435D+00,  
     *       4677455693272.880D+00,  
     *      -7876534745188.580D+00,  
     *       32246018.53578549D+00,
     *       195306248.4096138D+00/
C theta=100
      DATA XX7/
     *     -0.1728006534180419D+00,              
     *       10.72274703832234D+00,              
     *      -279.1825813640263D+00, 
     *       3744.798126721926D+00, 
     *      -23881.84308870698D+00, 
     *       476.1716308762578D+00, 
     *       1023339.224037222D+00, 
     *      -6199317.157442045D+00, 
     *       12354719.93357890D+00, 
     *      0.9940778707477032D+00, 
     *     -0.8516895773886548D+00, 
     *     -0.1618676362211799D+00, 
     *      -4028049211.139320D+00, 
     *       20537646895.05984D+00, 
     *      -30982419237.70815D+00, 
     *       32564168.05940399D+00,
     *       193690194.4198301D+00/
C theta=120
      DATA XX8/
     *     -0.1691016268642641D+00,            
     *       10.31418509253226D+00,            
     *      -263.2477801682790D+00, 
     *       3428.306888128504D+00, 
     *      -20292.82991756392D+00, 
     *      -23289.46849971248D+00, 
     *       1111174.439560930D+00, 
     *      -6359656.342728067D+00, 
     *       12460588.57757912D+00, 
     *       1.011227410876825D+00, 
     *      -1.158333328464265D+00, 
     *     -0.7601444704311453D+00, 
     *       36142579346.71783D+00, 
     *       61111755383.42400D+00, 
     *      -491377392110.9773D+00, 
     *       34444322.94723273D+00,
     *       200515052.0984532D+00/
C theta=140
      DATA XX9/
     *     -0.1241153142193939D+00,                       
     *       8.625773665852453D+00,                       
     *      -264.4262405893857D+00, 
     *       4653.355940988027D+00, 
     *      -51376.21459413337D+00, 
     *       364375.3706912036D+00, 
     *      -1622040.388552431D+00, 
     *       4145894.450466237D+00, 
     *      -4669985.813322528D+00, 
     *       1.168375058764390D+00, 
     *      -5.017720269989782D+00, 
     *     -0.5163018938498151D+00, 
     *      -586043854446.1162D+00, 
     *       3137090073156.085D+00, 
     *      -4434704096123.796D+00, 
     *       37265189.03008173D+00,
     *       209166778.2669497D+00/
C theta=160
      DATA XX10/
     *     -0.4670287750285096D+00,           
     *       30.15344063538946D+00,           
     *      -843.4983188541162D+00,  
     *       13072.85997994218D+00,  
     *      -119154.5474006862D+00,  
     *       608642.3152853674D+00,  
     *      -1289424.002214431D+00,  
     *      -1699839.658747914D+00,  
     *       9623907.851426231D+00,  
     *       1.066638520669773D+00,  
     *      -1.473845749728721D+00,  
     *      -1.658879660246379D+00,  
     *       82101608770995.61D+00,  
     *      -436180651761959.3D+00,  
     *       662739262453405.3D+00,  
     *       39656681.61530334D+00,
     *       215036700.4001355D+00/
C theta=180
      DATA XX11/
     *      -1.265675249002982D+00,               
     *       80.70694360877671D+00,               
     *      -2209.628545387617D+00, 
     *       32897.44970382784D+00, 
     *      -275837.7534370786D+00, 
     *       1119962.034941880D+00, 
     *       104065.2025598504D+00, 
     *      -18309799.86185131D+00, 
     *       48548363.21235426D+00, 
     *       1.049525066930224D+00, 
     *     -0.1485877045164195D+00, 
     *     -0.9569735416759715D+00, 
     *      -690239600787.7011D+00, 
     *       4679771181527.622D+00, 
     *      -7850908473831.484D+00, 
     *       40559439.26241566D+00,
     *       217292197.4162477D+00/
C
      IF(N.EQ.1) THEN
      DO I=1,17
      XX(I)=XX1(I)
      ENDDO
      ENDIF
      IF(N.EQ.2) THEN
      DO I=1,17
      XX(I)=XX2(I)
      ENDDO
      ENDIF
      IF(N.EQ.3) THEN
      DO I=1,17
      XX(I)=XX3(I)
      ENDDO
      ENDIF
      IF(N.EQ.4) THEN
      DO I=1,17
      XX(I)=XX4(I)
      ENDDO
      ENDIF
      IF(N.EQ.5) THEN
      DO I=1,17
      XX(I)=XX5(I)
      ENDDO
      ENDIF
      IF(N.EQ.6) THEN
      DO I=1,17
      XX(I)=XX6(I)
      ENDDO
      ENDIF
      IF(N.EQ.7) THEN
      DO I=1,17
      XX(I)=XX7(I)
      ENDDO
      ENDIF
      IF(N.EQ.8) THEN
      DO I=1,17
      XX(I)=XX8(I)
      ENDDO
      ENDIF
      IF(N.EQ.9) THEN
      DO I=1,17
      XX(I)=XX9(I)
      ENDDO
      ENDIF
      IF(N.EQ.10) THEN
      DO I=1,17
      XX(I)=XX10(I)
      ENDDO
      ENDIF
      IF(N.EQ.11) THEN
      DO I=1,17
      XX(I)=XX11(I)
      ENDDO
      ENDIF
C
      TMP1 = DEXP(-XX(10)*R-XX(11))
      TMP2 = 0.5d+00*(1.0d+00 + dtanh(1.0d+00 + XX(12)*R))
C
      VAPRIMI = (XX(1)*R**8 + XX(2)*R**7 + XX(3)*R**6 + XX(4)*R**5 +
     &           XX(5)*R**4 + XX(6)*R**3 + XX(7)*R**2 + XX(8)*R + 
     &           XX(9)) * TMP1 
     &           - TMP2 * (XX(16)/R**6 + XX(17)/R**7 +                    
     &                     XX(13)/R**8 + XX(14)/R**9 +XX(15)/R**10)
      TOCM = 0.219474643545745D0
      VAPRIMI = VAPRIMI*TOCM
C
      RETURN
      END
      FUNCTION VABISI(R,N)
CSYSTEM: OH-Ar Abis state Theta=0 for Ar--OH
CLEVEL:RHF/RCCSD(T)/aug-cc-pvqz+bf
C units R=au E=cm-1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(17),XX1(17),XX2(17),XX3(17),
     * XX4(17),XX5(17),XX6(17),XX7(17),XX8(17),
     * XX9(17),XX10(17),XX11(17)
C theta=0
      DATA XX1/   
     *     -0.5835284219401391D+00,        
     *       38.04952789322221D+00,        
     *      -985.9019652754446D+00,  
     *       11518.96508535197D+00,  
     *      -21506.60006303519D+00,  
     *      -1007794.615872079D+00,  
     *       12184854.96419048D+00,  
     *      -59734675.55243188D+00,  
     *       112230765.7760808D+00,  
     *      0.9298477185884633D+00,  
     *      0.3179101021369883D+00,  
     *      -1.217028123076279D+00,  
     *      -33195530529219.63D+00,  
     *       240524473649379.3D+00,  
     *      -453384045468199.6D+00,  
     *       22784339.37813496D+00,
     *       673550266.2988292D+00/
C theta=20
      DATA XX2/
     *     -0.3269836892734609D+00,        
     *       23.27164321351734D+00,        
     *      -720.1816275482132D+00,    
     *       12511.70461733529D+00,    
     *      -131746.2516154350D+00,    
     *       840789.5280734375D+00,    
     *      -3012019.429192393D+00,    
     *       4651583.990301123D+00,    
     *       20683.57602622089D+00,    
     *       1.032779970175058D+00,    
     *      -2.053136978342261D+00,    
     *     -0.8035968671578722D+00,    
     *      -826268544890.2833D+00,    
     *       5116877951189.167D+00,    
     *      -8395055720906.515D+00,    
     *       26106126.63138232D+00,
     *       583173073.2842053D+00/
C theta=40
      DATA XX3/
     *     -0.1739299410142206D+00,           
     *       11.81554386501734D+00,           
     *      -347.9344798507028D+00, 
     *       5700.196685941467D+00, 
     *      -55500.30163193998D+00, 
     *       311906.6960722521D+00, 
     *      -834774.3076728931D+00, 
     *      -8014.753722898944D+00, 
     *       3563905.138441560D+00, 
     *       1.028003938319646D+00, 
     *      -2.169309853474742D+00, 
     *      -1.451850145199881D+00, 
     *      -43376931175914.34D+00, 
     *       313446666655857.2D+00, 
     *      -571386219832440.1D+00, 
     *       30564549.86610075D+00,
     *       421732359.1335933D+00/
C theta=60
      DATA XX4/
     *    -0.2266416477567826D+00,        
     *      14.56059035794553D+00,        
     *     -402.4456292530374D+00, 
     *      6065.645739741473D+00, 
     *     -51695.04816430719D+00, 
     *      214441.7830612679D+00, 
     *      19329.15091081641D+00, 
     *     -3717093.223449220D+00, 
     *      10229579.25738579D+00, 
     *      1.020946360415995D+00, 
     *     -1.330639371814747D+00, 
     *     -1.611671740292514D+00, 
     *     -141336220954002.1D+00, 
     *      1033972436285405.D+00, 
     *     -1903076617500239.D+00, 
     *      31569450.45284114D+00,
     *      319896898.8036590D+00/
C theta=80
      DATA XX5/
     *     -0.2472596758875793D+00,           
     *       15.18982773709356D+00,           
     *      -394.3331150656819D+00,  
     *       5312.843808070574D+00,  
     *      -34272.26422073379D+00,  
     *      -1441.054535451151D+00,  
     *       1588160.806436515D+00,  
     *      -10094178.39025672D+00,  
     *       21440670.54848527D+00,  
     *       1.005013549903144D+00,  
     *     -0.6518161108745924D+00,  
     *      -1.262070519745715D+00,  
     *      -9103535780124.301D+00,  
     *       63414157753372.48D+00,  
     *      -111578556277434.7D+00,  
     *       30601142.78600957D+00,
     *       294712494.7967941D+00/
C theta=90
      DATA XX6/
     *     -0.2125507754136597D+00,         
     *       13.52541749465000D+00,         
     *      -364.0253759292552D+00, 
     *       5091.695039300854D+00, 
     *      -34366.49789146336D+00, 
     *       7546.784122432019D+00, 
     *       1573084.849199305D+00, 
     *      -10335205.87144623D+00, 
     *       22330342.26970545D+00, 
     *      0.9809511701695504D+00, 
     *     -0.4618685855351585D+00, 
     *     -0.6987866404894252D+00, 
     *      -140419739301.0930D+00, 
     *       817264623803.9906D+00, 
     *      -1218808549890.561D+00, 
     *       30402323.21417024D+00,
     *       295817496.1387828D+00/
C theta=100
      DATA XX7/
     *     -0.2197263282412382D+00,        
     *       13.77552177949496D+00,        
     *      -365.5395754202103D+00,   
     *       5034.671001542191D+00,   
     *      -33133.73926978195D+00,   
     *      -2706.307648953639D+00,   
     *       1610241.431812328D+00,   
     *      -10362667.70073688D+00,   
     *       22237615.26718034D+00,   
     *      0.9931456312428649D+00,   
     *     -0.5769689524755101D+00,   
     *     -0.7460304096054683D+00,   
     *      -181443967975.5101D+00,   
     *       1086471330536.616D+00,   
     *      -1653871868817.420D+00,   
     *       30776363.72900048D+00,
     *       296947897.1959993D+00/
C theta=120
      DATA XX8/
     *      -1.071702356938537D+00,          
     *       68.62042755134195D+00,          
     *      -1916.642352189165D+00, 
     *       29957.08803666413D+00, 
     *      -280079.1870313794D+00, 
     *       1521774.969163422D+00, 
     *      -3956914.235229969D+00, 
     *      -217534.9857797138D+00, 
     *       17244877.74783726D+00, 
     *       1.073265376184860D+00, 
     *     -0.6110497681519810D+00, 
     *      -1.428964427858465D+00, 
     *      -46024223834937.15D+00, 
     *       333208283889870.8D+00, 
     *      -606633602357302.8D+00, 
     *       33112959.99239098D+00,
     *       287174839.4862120D+00/
C theta=140
      DATA XX9/
     *     -0.9685731944485279D+00,            
     *       62.60229745685347D+00,            
     *      -1762.764777013899D+00,    
     *       27727.50987180000D+00,    
     *      -260211.3133289935D+00,    
     *       1414029.818794777D+00,    
     *      -3655636.065699264D+00,    
     *      -260217.1026901831D+00,    
     *       15793836.24302723D+00,    
     *       1.070840991358049D+00,    
     *     -0.7633104489080703D+00,    
     *     -0.9461812462741014D+00,    
     *      -830202803702.4766D+00,    
     *       5359810538923.443D+00,    
     *      -8683666715388.859D+00,    
     *       36595352.19192874D+00,
     *       259101301.2057308D+00/
C theta=160
      DATA XX10/
     *      -1.192453674950103D+00,           
     *       74.83241190524375D+00,           
     *      -2019.789077699435D+00,  
     *       29688.33236055743D+00,  
     *      -245931.1619116676D+00,  
     *       981866.9813181594D+00,  
     *       206113.9433240369D+00,  
     *      -16725510.68550027D+00,  
     *       44085686.38402973D+00,  
     *       1.059836525806998D+00,  
     *     -0.2858906663225153D+00,  
     *      -2.044295686649257D+00,  
     *       180671376196946.4D+00,  
     *       468862781411974.9D+00,  
     *      -3590390962861951.D+00,  
     *       39526663.85093726D+00,
     *       229026546.9722313D+00/
C theta=180
      DATA XX11/
     *      -1.265675249002982D+00,           
     *       80.70694360877671D+00,           
     *      -2209.628545387617D+00, 
     *       32897.44970382784D+00, 
     *      -275837.7534370786D+00, 
     *       1119962.034941880D+00, 
     *       104065.2025598504D+00, 
     *      -18309799.86185131D+00, 
     *       48548363.21235426D+00, 
     *       1.049525066930224D+00, 
     *     -0.1485877045164195D+00, 
     *     -0.9569735416759715D+00, 
     *      -690239600787.7011D+00, 
     *       4679771181527.622D+00, 
     *      -7850908473831.484D+00, 
     *       40559439.26241566D+00,
     *       217292197.4162477D+00/
C
      IF(N.EQ.1) THEN
      DO I=1,17
      XX(I)=XX1(I)
      ENDDO
      ENDIF
      IF(N.EQ.2) THEN
      DO I=1,17
      XX(I)=XX2(I)
      ENDDO
      ENDIF
      IF(N.EQ.3) THEN
      DO I=1,17
      XX(I)=XX3(I)
      ENDDO
      ENDIF
      IF(N.EQ.4) THEN
      DO I=1,17
      XX(I)=XX4(I)
      ENDDO
      ENDIF
      IF(N.EQ.5) THEN
      DO I=1,17
      XX(I)=XX5(I)
      ENDDO
      ENDIF
      IF(N.EQ.6) THEN
      DO I=1,17
      XX(I)=XX6(I)
      ENDDO
      ENDIF
      IF(N.EQ.7) THEN
      DO I=1,17
      XX(I)=XX7(I)
      ENDDO
      ENDIF
      IF(N.EQ.8) THEN
      DO I=1,17
      XX(I)=XX8(I)
      ENDDO
      ENDIF
      IF(N.EQ.9) THEN
      DO I=1,17
      XX(I)=XX9(I)
      ENDDO
      ENDIF
      IF(N.EQ.10) THEN
      DO I=1,17
      XX(I)=XX10(I)
      ENDDO
      ENDIF
      IF(N.EQ.11) THEN
      DO I=1,17
      XX(I)=XX11(I)
      ENDDO
      ENDIF
C
      TMP1 = DEXP(-XX(10)*R-XX(11))
      TMP2 = 0.5d+00*(1.0d+00 + dtanh(1.0d+00 + XX(12)*R))
C
      VABISI = (XX(1)*R**8 + XX(2)*R**7 + XX(3)*R**6 + XX(4)*R**5 +
     &          XX(5)*R**4 + XX(6)*R**3 + XX(7)*R**2 + XX(8)*R + 
     &          XX(9)) * TMP1 
     &       - TMP2 * (XX(16)/R**6 + XX(17)/R**7 +                    
     &                 XX(13)/R**8 + XX(14)/R**9 + XX(15)/R**10)
      TOCM=0.219474643545745D0
      VABISI = VABISI*TOCM
C
      RETURN
      END
      

