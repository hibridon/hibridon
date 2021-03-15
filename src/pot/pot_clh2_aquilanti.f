*system:  Cl(2P)+H2, Dubernet-Hutson expansion of Aquilanti et al PES's
*references:  
*   V. Aquilanti, D. Cappelletti, S. Cavalli, F. Pirani, and A. Volpi,
*        J. Phys. Chem. A 105, 2401 (2001).
*   M.-L. Dubernet and J. M. Hutson, J. Chem. Phys. 101, 1939 (1994).
*   M. H. Alexander, J. Chem. Phys. 99, 6014 (1993).
*   M. H. Alexander and M. Yang, J. Chem. Phys. 103, 7956 (1995).


      include "common/syusr"
      include "common/ground"
      include "common/bausr"
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      common /coselb/ ibasty
      potnam='Cl(2P)-H2 Aquilanti PES expan.  DUBERNET-HUTSON'
      ibasty=12
      lammin(1)=1
      lammax(1)=9
      mproj(1)=0
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* --------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      character *48 potnam
      character *2 frame
      logical csflag, ljunk, ihomo, lljunk
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(9)
      potnam='Cl(2P)-H2 Aquilanti  PES expan.  DUBERNET-HUTSON'
      print *, potnam
      print *
1      print *, ' r (bohr) frame (sf/bf)'
      read (5, *, end=93) r, frame
      if (frame.eq.'sf') then
        csflag=.false.
        ihomo=.true.
      else if (frame .eq. 'bf') then
        csflag=.true.
        ihomo=.false.
      else
        print *, 'frame must be either "sf" or "bf"'
        go to 1
      endif
      if (r .lt. 0.d0) go to 93
      call pot(vv0,r)
      if (.not. csflag .or. (csflag .and. ihomo)) write (6, 100) vv0,vvl
100   format(' V000, V220, V022, V202:  ',4(1pe16.8),/,
     :       ' V222, V224, V404:  ',3(1pe16.8),/,
     :       ' V422, V424, V426:  ',3(1pe16.8))
      if (csflag .and. .not.ihomo) write (6, 110) vv0,vvl
110   format(' v000, v220, v020, v200:  ',4(1pe16.8),/,
     :       ' v222, v221, v400:  ',3(1pe16.8),/,
     :       ' v420, v422, v421:  ',3(1pe16.8))
      goto 1
93    r=4
      do i=1,100
       call pot(vv0,r)
       write(2,101) r,vv0,vvl
101    format(f8.4,10(1pe16.8))
       r=r+0.2
      enddo

99    end
      subroutine pot (vv0, r)
* ----------------------------------------------------------------------
*  subroutine to calculate the r-dependent coefficients in the
*  al-h2 potentials of williams and alexander using body-frame
*  expansion of dubernet and flower
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*  vv0       totally symmetric potential (v000)
*  variable in common block /covvl/
*    vvl:     vector of length 5 to store r-dependence of each term
*             in potential expansion
*  CS calculations case 1A (csflag=.true. and ihomo=.false)
*  body frame expansion coefficients
*    vvl(1)  v220(R)
*    vvl(2)  v020(R)
*    vvl(3)  v200(R)
*    vvl(4)  v222(R)
*    vvl(5)  v221(R)
*    vvl(6)  v400(R)
*    vvl(7)  v420(R)
*    vvl(8)  v422(R)
*    vvl(9)  v421(R)
*  CC calculations (csflag=.false.)
*    vvl(1)  V220(R)
*    vvl(2)  V022(R)
*    vvl(3)  V202(R)
*    vvl(4)  V222(R)
*    vvl(5)  V224(R)
*    vvl(6)  V404(R)
*    vvl(7)  V422(R)
*    vvl(8)  V424(R)
*    vvl(9)  V426(R)

* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  1-feb-1998
* ----------------------------------------------------------------------


      implicit double precision (a-h,o-z)
      logical ljunk, ihomo, csflag, lljunk
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(9)
*      dimension beta(5)
      dimension d0(25),d1(12),d2(16),aa(25)
      dimension fthet(3,5)
      dimension kpvt(5),
     :          qraux(5), work(24),rsd(5)
      dimension vvll(18)
*    vvll(1-5) expansion coefficients in dl0 (l=0,2,4,6,8) of vsigma
*    vvll(6-10) expansion coefficients in dl0 (l=0,2,4,6,8) of vsum
*    vvll(11-14) expansion coefficients in dl2 (l=2,4,6,8) of vdif
*    vvll(15-18) expansion coefficients in dl1 (l=2,4,6,8) of v1
* angles (not needed here, but included for completeness)
*      data beta / 0d0, 22.5d0, 45d0, 67.5d0, 90d0/
      data zero, one, half /0.d0,1.d0,0.5d0/
      data two, three, four, five /2.d0,3.d0,4.d0,5.d0/
      data seven,eleven /7.d0,11.d0/
      data sq2 /1.414213562d0/
* coefficients for d0 rotation matrices
* stored (by column) for each of 5 angles and for l=0,2,4,6,8
      data d0/
     : 1.0000000d0, 1.0000000d0, 1.0000000d0, 1.0000000d0, 1.0000000d0,
     : 1.0000000d0, 7.8033009d-01, 2.5000001d-01, -2.8033008d-01,
     : -5.0000000d-01, 1.0000000d0, 3.6159588d-01, -4.0625000d-01,
     : -8.0345887d-02, 3.7500000d-01, 1.0000000d0, -7.6358299d-02,
     : -1.4843752d-01, 2.7167082d-01, -3.1250000d-01, 1.0000000d0,
     : -3.5735359d-01, 2.9833984d-01, -2.7863272d-01, 2.7343750d-01/
* coefficients for d1 rotation matrices
* stored (by column) for each of 3 angles and for l=2,4,6,8
      data d1/
     : -4.3301270d-01, -6.1237244d-01, -4.3301270d-01,
     : -5.8796105d-01, -1.3975425d-01,  3.9031869d-01,
     : -4.9200540d-01,  3.5441551d-01, -1.8822068d-01,
     : -2.1012777d-01,  1.1186650d-01, -5.0614421d-02/
* coefficients for d2 rotation matrices
* stored (by column) for each of 4 angles and for l=2,4,6,8
      data d2/
     : 8.9679865d-02,  3.0618621d-01,  5.2269256d-01,  6.1237244d-01,
     : 2.8798601d-01,  4.9410589d-01,  8.4775334d-03, -3.9528471d-01,
     : 4.5386125d-01,  4.0027167d-02, -2.5372551d-01,  3.2021721d-01,
     : 4.8368899d-01, -3.2931303d-01,  2.8759680d-01, -2.7731624d-01/
* hyperbolic tangent scaling factor
      data alph /0.6d0/
      cm=2.194746d05
      pi=acos(-1.0d0)
* first for vzz
      tol=1.e-10
c      call dcopy(25,d0,1,aa,1)
c      call dqrank(aa,5,5,5,tol,kr,kpvt,qraux,work)
c      call dqrlss(aa,5,5,5,kr,vzz,xzz,rsd,kpvt,qraux)
c      call dcopy(5,xzz,1,vvll,1)
c      call dqrlss(aa,5,5,5,kr,vs,xs,rsd,kpvt,qraux)
c      call dcopy(5,xs,1,vvll(6),1)
*
c      call dcopy(16,d2,1,aa,1)
c      call dqrank(aa,4,4,4,tol,kr,kpvt,qraux,work)
c      call dqrlss(aa,4,4,4,kr,vd(2),xd,rsd,kpvt,qraux)
c      call dcopy(4,xd,1,vvll(11),1)
*
c      call dcopy(12,d1,1,aa,1)
c      call dqrank(aa,3,3,4,tol,kr,kpvt,qraux,work)
c      call dqrlss(aa,3,3,4,kr,vxz(2),x1,rsd,kpvt,qraux)
c      call dcopy(4,x1,1,vvll(15),1)
      rl=r*0.529177249d0
      rh2=1.40112d0*0.529177249d0
c  Vsigma
      CALL diabats(rl,0.D0,VSG,VPI,VPI2,V12)
      CALL diabats(rl,90.D0,VA1,VB2,VB1,V122)
      vvll(1)=1.D0/3.D0*(2.D0*VA1+VSG)
      vvll(2)=2.D0/3.D0*(VSG-VA1)
      vvll(3)=0.D0
      vvll(4)=0.D0
      vvll(5)=0.D0
c Vsum
      vvll(6)=1.D0/3.D0*(VB1+VB2+
     >  VPI)
      vvll(7)=2.D0/3.D0*(VPI-
     >       1.D0/2.D0*(VB1+VB2))
      vvll(8)=0.D0
      vvll(9)=0.D0
      vvll(10)=0.D0
cVdif minus sign before Vdif needed in Hibridon for
c correct result
      vvll(11)=-(2.D0/DSQRT(6.D0))*(VB1-VB2)
      vvll(12)=0.d0
      vvll(13)=0.D0
      vvll(14)=0.D0
c V12
c      vvll(15)=0.D0
c      vvll(16)=0.D0
c      vvll(17)=0.D0
c      vvll(18)=0.D0
c      vvll(15)=POTV12VLM(rl,rh2,2,1)/DSQRT(2.D0)
c      vvll(16)=POTV12VLM(rl,rh2,4,1)/DSQRT(2.D0)
c      vvll(17)=POTV12VLM(rl,rh2,6,1)/DSQRT(2.D0)
c      vvll(18)=POTV12VLM(rl,rh2,8,1)/DSQRT(2.D0)
       vvll(15)=POTV12VLM(rl,rh2,2,1)
       vvll(16)=POTV12VLM(rl,rh2,4,1)
       vvll(17)=POTV12VLM(rl,rh2,6,1)
       vvll(18)=POTV12VLM(rl,rh2,8,1)
* determine body frame expansion coefficients in dubernet and flower
* expansion
*----------------Jason's testing
      junk=3.0d0
*---------------------
* here is totally symmetric term
      vv0=(two*vvll(6)+vvll(1))/three
      sq6=1.d0/sqrt(6.d0)
* v200
      vvl(1)=(two*vvll(7)+vvll(2))/three
* v400
      vvl(6)=(two*vvll(8)+vvll(3))/three
* v020
      vvl(2)=(vvll(1)-vvll(6))*five/three
* v220
      vvl(3)=(vvll(2)-vvll(7))*five/three
* v420
      vvl(7)=(vvll(3)-vvll(8))*five/three
* v222
      vvl(4)=five*vvll(11)*sq6
* v422
      vvl(9)=five*vvll(12)*sq6
* v221
      vvl(5)=five*vvll(15)*sq6
* v421
      vvl(8)=five*vvll(16)*sq6
* transform to space frame terms if CC or CS case 1C
      if (.not. csflag .or. (csflag .and. ihomo)) then
        call dcopy(9,vvl,1,vvll,1)
        vvl(1)=(two*(vvll(4)-vvll(5))+vvll(3))/sqrt(5.d0)
        vvl(2)=vvll(2)
        vvl(3)=vvll(1)
        vvl(4)=(two*vvll(4)+vvll(5)-vvll(3))*sqrt(two/7.d0)
        vvl(5)=(vvll(4)+four*vvll(5)+3*vvll(3))*sqrt(two/35.d0)
        vvl(6)=vvll(6)
        vvl(7)=sqrt(two/seven)*vvll(7)-sqrt(20.d0/21.d0)*vvl(8)
     :        +sqrt(10d0/21d0)*vvl(9)
        vvl(8)=sqrt(one/77d0)*(-sqrt(20d0)*vvll(7)+sqrt(6d0)*vvl(8)
     :        +sqrt(108d0)*vvl(9))
        vvl(9)=sqrt(five/eleven)*vvll(7)+sqrt(32.d0/33.d0)*vvl(8)
     :        +sqrt(four/33d0)*vvl(9)
      else
* reorder body-frame terms
        vvll(1)=vvl(1)
        vvl(1)=vvl(3)
        vvl(3)=vvll(1)
        vvll(8)=vvl(9)
        vvl(9)=vvl(8)
        vvl(8)=vvll(8)
      endif
* convert to hartree
      econv=3.674930886769D-5
      do i=1,18
      vvl(i)=vvl(i)*econv
      enddo
c      econv=1./219474.6
c      call dscal(9,econv,vvl,1)
      vv0=vv0*econv

      end
C FUNCTIONS FOR A1,B1,B2,SIGMA AND PI STATES
C CC 2D MODEL OF THOSE STATES
C BY J. KLOS COPYRIGHT HE,HE :)
C*************START V12*******************************
      FUNCTION POTV12VLM(R,rsmall,L,M)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     THIS PROGRAM CALCULATES THE RADIAL V_L1 COEFFICIENTS
C     FOR V12=H12 MRCISD F-H2 DIABATIC POTENTIAL
C     REFERENCE TO AB INITIO CALCULATIONS:
C     Int.J.Quant.Chem., in process
      CALL SETUP
      CALL QGAUS12(R,rsmall,-1.d0,1.d0,L,M,SS)
      POTV12VLM=SS
      end

      SUBROUTINE QGAUS12(R,rsmall,A,B,L,M,SS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C        RETURNS AS SS THE INTEGRAL OF THE FUNCTION FUNC
C         BETWEEN A AND B,BY TEN-POINT GAUSS-LEGENDRE INTEGRATION:
C        THE FUNCTION IS EVALUATED EXACTLY TEN TIMES AT INTERIOR
C        POINTS IN THE RANGE OF INTEGRATION.
C      INTEGER J
       DIMENSION W(50),X(50)
       DIMENSION FN(400)
       COMMON/FACTOR/FN
       SAVE W,X
c      DATA W/0.2955242247,0.2692667193,0.2190863625,
c     1 0.1494513491,0.0666713443/
c      DATA X/0.1488743389,0.4333953941,0.6794095682,
c     1 0.8650633666,0.9739065285/
       CALL gauleg(-1.D0,1.D0,X,W,50)
      XM = 0.5D0*(B + A)
      XR = 0.5D0*(B - A)
      SS=0.D0
      DO  J=1,50
      TERM1=FN(L-M+1)
      TERM2=FN(L+M+1)
      FACT = DEXP(TERM1-TERM2)
      FACT = DSQRT(FACT)
       DX = XR*X(J)
      CALL diabats(R,dacos(XM+DX)*180.d0/dacos(-1.D0),V1,V2,V3,V12)
       SS = SS + W(J)*(V12*
     1 plgndr(L,M,(XM+DX))*(dfloat(2*L+1)/2.d0))*FACT
      ENDDO
       SS = XR*SS
      RETURN
      END
C************END V12***************************
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
      SUBROUTINE gauleg(x1,x2,x,w,n)
      INTEGER n
      DOUBLE PRECISION x1,x2,x(n),w(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END





CFrom volpinho@impact.dyn.unipg.it Mon Oct 14 12:46:28 2002
CDate: Mon, 14 Oct 2002 10:56:02 +0200 (MET DST)
CFrom: Volpi Alessandro <volpinho@impact.dyn.unipg.it>
CTo: jakl@theochem.kun.nl
CSubject: Cl-H2 routine

CComputation of the three diabatic surfaces and
Ccouplings for the Cl + H_2 system.
CThe routine works exactly the same as the one 
Cyou got for F + H_2, so no spin-orbit is included.
CAs before, you will find inside the routine 
Ccomments about call, units, etc. 
C
CBest regards

CAlessandro




C===================================================
       subroutine diabats(R,theta,V1,V2,V3,V12)
C===================================================
c  The routine supplies the value of the diabatic
c  electronic states for the Cl + H2 system as a function
c  of the Jacobi vectors R and of the theta
c  angle between R and the internuclear distance
c  r (H-H). The vibrational coordinate r(H-H) is kept fixed
c  at the equilibrium value of 1.40 bohr (~0.74 Angstrom).
C===================================================
C   IMPORTANT !!!!!!!!!!!
c The surface computed does NOT include spin orbit
c coupling and is widely described in
c
c   V. Aquilanti, D. Cappelletti, S. Cavalli,
c   F. Pirani, and A. Volpi,
c   J. Phys. Chem. A 105, 2401 (2001).
c   Please see that reference for all details.
c
C===================================================
C===================================================
c  R         --->  Angstrom
c  theta     --->  degrees
c  V1        --->  1A'  (millieV)
c  V2        --->  2A'  (millieV)
c  V3        --->  A"   (millieV)
c  V12       --->  diabatic coupling between
c                  1A' and 2A'   (millieV)
C===================================================
      IMPLICIT REAL*8 (A-H,O-Z)                                

	parameter(pigreco = 3.1415926535897932)
      COMMON/BLK1/EPSI0,RM0,BETA,C000,X3,X4,B1,B2,B3,B4    
      COMMON/BLK2/X1,X2,A1,A2,A3,A4,AEXP,ALFA1       
      COMMON/BLK3/A202,ALFA202,C202
      COMMON/BLK4/C022,A022,ALFA022,GAMMA022           
      COMMON/BLK5/A220,ALFA220      
      COMMON/BLK6/A222,ALFA222      
      COMMON/BLK7/C224      

c parametri termini radiali

c V000 (MSV)

	RM0 = 3.7000d0                                           
	EPSI0= 5.700d0                                
	BETA = 6.300D0                                  
	C000 = 26000.0d0                               
 	AEXP= 0.0000D0                               
 	ALFA1=0.00000D0                                    
c	AEXP= 0.3800D0                               
c	ALFA1=17.00000D0                                    
c	X1 = 0.6900D0                                  
c	X2 = 0.8600D0                                  
 	X1 = 0.0000D0                                  
 	X2 = 0.0000D0                                  
 	X3 = 1.1000D0                                            
 	X4 = 1.6000D0   

c V202
c	A202=2.0d5
c	ALFA202=4.000d0

 	A202=4.26d5
 	ALFA202=3.500d0

	C202=C000*0.0928d0

c V022
	C022=4680d0 
	A022=7.0000d5
	ALFA022=3.0d0
	GAMMA022=0.0000D0                                                              
c V220
 	A220=A202
c	A220=0.0d0
	ALFA220=ALFA202

c V222
 	A222=-A202
c	A222=0.0d0
	ALFA222=ALFA202

c V224
	C224=0.49d0*4.059d0*1129.0d0


      CALL SPLINE_PES                    
c     write(*,11) b1,b2,b3,b4
11    format(4f11.4)

	c1=5.d0*dsqrt(5.0d0)
	c2=2.0d0*dsqrt(2.0d0)/(5.0d0*dsqrt(7.0d0))
	c3=6.0d0*dsqrt(2.0d0)/(5.0d0*dsqrt(35.d0))
	c4=9.0d0/(10.0d0*dsqrt(70.0d0))
C
C=======chiamate alle routines che calcolano i coeff. Vabc===
 	call v000(r,p000)
	call v202(r,p202)
 	call v022(r,p022)
 	call v220(r,p220)
 	call v222(r,p222)
 	call v224(r,p224)
C===============================
          thetarad = theta*pigreco/180.0d0
c========espressioni delle armoniche===================
          c20 = 0.5d0*(3.0d0*dcos(thetarad)*dcos(thetarad) - 1.0d0)
          c21 = - dsqrt(1.5d0)*dcos(thetarad)*dsin(thetarad)
          c22 = dsqrt(3.0d0/8.0d0)*dsin(thetarad)*dsin(thetarad)
c-----------------------------------------------------------
          var1 = p220/dsqrt(5.0d0) - p222*dsqrt(2.0d0/7.0d0)
     >           + p224*3.0d0*dsqrt(2.0d0/35.0d0)
          var2 = p220/dsqrt(5.0d0) + dsqrt(2.0d0/7.0d0)*p222       
     >           + p224/dsqrt(70.0d0)
          var3 = -p220/dsqrt(5.0d0) + p222/dsqrt(14.0d0) +
     >           p224*2.0d0*dsqrt(2.0d0/35.0d0)

          var4 = p000 - p022/5.0d0 + (p202 - var1/5.0d0)*c20
          var5 = dsqrt(6.0d0)*var2*c22/5.0d0

c==========espressioni degli stati diabatici================
        V1 = p000 + 2.0d0/5.0d0*p022 + 2.0d0*(var1 + p202
     >       *5.0d0/2.0d0)*c20/5.0d0

        V2 = var4 + var5

        V3 = var4 -var5
  
        V12 = dsqrt(6.0d0)*var3*c21/5.0d0

C=======write sulle superfici=============================
c         write(6,12) R,V1,V2,V3,V12
C=========================================================

c12	format(5(f11.5))
	return
	end

C=======SUBROUTINES UTILIZZATE=============================


       REAL*8 FUNCTION P2(x)
       IMPLICIT REAL*8(A-H,O-Z)
       P2=0.5d0*(3.d0*(DCOS(x))**2-1.d0)
       RETURN
       END


      SUBROUTINE V000(R,V)                        
      IMPLICIT REAL*8 (A-H,O-Z)                             
      COMMON/BLK1/EPSI0,RM0,BETA1,C6,X3,X4,B1,B2,B3,B4
      COMMON/BLK2/X1,X2,A1,A2,A3,A4,AEXP,ALFA1              
      X=R/RM0                                                          
      IF(X.LE.X1) GOTO160                   
      IF(X.LT.X2) GOTO161                                      
      IF(X.LE.X3) GOTO162                                     
      IF(X.LT.X4) GOTO163                               
C.......................VAN DER WAALS...........          
        R2=R/100.d0                                                    
        V=-(1.D-12*C6/R2**6)                                  
      GOTO164                                               
C....................EXPONENTIAL........................           
160     SP=-ALFA1*(X-1.d0)                                            
        V= EPSI0*AEXP*DEXP(SP)                               
      GOTO164                                                 
C.......................SPLINE1...................          
161     SP=A3+(X-X1)*A4                                              
        SP=A2+(X-X2)*SP                                         
        SP=A1+(X-X1)*SP                                        
        V=EPSI0*DEXP(SP)                                      
      GOTO164                                              
C.........................MORSE...............................           
162     SP=-BETA1*(X-1.d0)                         
        V=EPSI0*(DEXP(2.d0*SP)-2.d0*DEXP(SP))                 
      GOTO164                                                   
C..........................SPLINE2.............              
163     SP=B3+(X-X3)*B4                                          
        SP=B2+(X-X4)*SP                                            
        SP=B1+(X-X3)*SP                                           
        V=EPSI0*SP                                          
164   RETURN                                               
      END                                                       

      SUBROUTINE V022(R,V)
      IMPLICIT REAL*8(A-H,O-Z)                                    
      COMMON/BLK4/C022,A022,ALFA022,GAMMA022           
        V=-A022*DEXP(-ALFA022*R-GAMMA022*R*R)+C022/R**6.d0           
      RETURN 
      END                                                                     

      SUBROUTINE V202(R,V)
      IMPLICIT REAL*8(A-H,O-Z)                                    
      COMMON/BLK3/A202,ALFA202,C202
        V=A202*DEXP(-ALFA202*R)-C202/R**6.d0           
      RETURN 
      END                                                                     

      SUBROUTINE V220(R,V)
      IMPLICIT REAL*8(A-H,O-Z)                                    
      COMMON/BLK5/A220,ALFA220      
        V=A220*DEXP(-ALFA220*R)          
      RETURN 
      END                                                                     

      SUBROUTINE V222(R,V)
      IMPLICIT REAL*8(A-H,O-Z)                                    
      COMMON/BLK6/A222,ALFA222      
        V=A222*DEXP(-ALFA222*R)          
      RETURN 
      END                                                                     

      SUBROUTINE V224(R,V)
      IMPLICIT REAL*8(A-H,O-Z)                                    
      COMMON/BLK7/C224      
	R70=8.3666003d0
	V=R70*C224/R**5
      RETURN
      END


      SUBROUTINE V2BUCK(R,V)
      IMPLICIT REAL*8(A-H,O-Z)                                    

	RD=2.063d0*2.29d0
	A2=360.d3
	B2=4.0d0
	C62=475.6d0
	C82=6.0312d3
	C102=31.842d3

	IF(R.LT.RD) THEN
	AEXP=-1.0d0*(RD/R-1.0d0)**2
	F=DEXP(AEXP)
	ELSE
	F=1.0d0
	END IF

	VA=-C62/R**6-C82/R**8-C102/R**10
	VA=VA*F
	VR=A2*DEXP(-B2*R)

	V=VR+VA

      RETURN
      END

      SUBROUTINE SPLINE_PES                                              
      IMPLICIT REAL*8 (A-H,O-Z)                                   
      COMMON/BLK1/EPSI0,RM0,BETA1,C6,X3,X4,B1,B2,B3,B4       
      COMMON/BLK2/X1,X2,A1,A2,A3,A4,AEXP,ALFA1            
        SP=-BETA1*(X3-1.d0)                                  
        V1 =EPSI0*(DEXP(2.d0*SP)-2.d0*DEXP(SP))                     
        SP=DEXP(-BETA1*(X3-1.d0))                                
        DV1 = EPSI0*(-2.d0*BETA1*SP*(SP-1.d0))/RM0            
        R1=X4*RM0/100.d0                               
        V2=-(1.0D-12*C6/R1**6)                                   
        DV2 = 6.D-14*C6/R1**7                                  
        B1 = V1/EPSI0                                         
        B2 = (V2/EPSI0-B1)/(X4-X3)                             
        B3 = (DV1*RM0/EPSI0-B2)/(X3-X4)                          
        B4 = (DV2*RM0/EPSI0-B2-B3*(X4-X3))/(X4-X3)**2      
        IF(AEXP.LT.0.01d0) GOTO 17                        
        SP=-ALFA1*(X1-1.d0)                                        
        A1=DLOG(AEXP)+SP                                
        SP=-BETA1*(X2-1.d0)                                 
        V =(DEXP(2.d0*SP)-2.d0*DEXP(SP))                      
        A2=(DLOG(V)-A1)/(X2-X1)                                
        A3= (ALFA1+A2)/(X2-X1)                                     
        SP=DEXP(-BETA1*(X2-1.d0))                                  
        DV=(-2.d0*BETA1*(SP-1.d0))/(SP-2.d0)                           
        A4= (DV-A2-A3*(X2-X1))/((X2-X1)**2)                    
17    RETURN                                                 
      END                                                      



C******************************************************
C       Dr. Alessandro Volpi                           
C       Dipartimento di Chimica                        
C       Universita' degli Studi di Perugia             
C       Via Elce di Sotto, 8                           
C       06123 Perugia (PG) Italia                      
C                                                      
C                              Tel.:+(39) 075 5855510  
C                              Fax: +(39) 075 5855606  
C                e-mail: volpinho@hermes.dyn.unipg.it  
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
