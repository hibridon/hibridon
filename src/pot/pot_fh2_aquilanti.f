*system:  F(2P)+H2, Dubernet-Hutson expansion of Aquilanti et al PES's
*references: .
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
      potnam='F(2P)-H2 Aquilanti et al PES expan.  DUBERNET-HUTSON'
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
      potnam='F(2P)-H2 Aquilanti et al. PES  expan.  DUBERNET-HUTSON'
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




c***H12 COUPLING 3-D MRCI****************
      FUNCTION POT3DH12(RR,rsmall,T)
C F-H2 system STATE:H12 MRCI 
C units R=A E=mikrohartree
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(35)
        DATA XX/
     *    0.146785045039485578D+02,
     *    0.639375370693218148D+01,
     *    0.151309417676178299D+03,
     *    0.274355254615448914D+02,
     *    0.182771692231812133D+01,
     *    0.106138270858910547D+00,
     *   -0.257751121801726955D-01,
     *    0.178071479524053472D+00,
     *   -0.199256700458829187D+03,
     *   -0.348623967196420921D+02,
     *   -0.189030316773231233D+01,
     *    0.141927268991115241D+00,
     *    0.396013774912161512D-01,
     *   -0.201682930487279699D+00,
     *    0.880168246801674172D+02,
     *    0.146941474850966749D+02,
     *    0.625184715211922315D+00,
     *   -0.166144804878691321D+00,
     *   -0.867688162199622276D-02,
     *    0.763923456040812299D-01,
     *   -0.132195819894863451D+02,
     *   -0.206769189296769618D+01,
     *   -0.628160580822506864D-01,
     *    0.358512231466525388D-01,
     *   -0.105914455880216072D-02,
     *   -0.970636190221498517D-02,
     *    0.799138571452743122D+03,
     *    0.102809848102600995D+03,
     *    0.733001295294727928D+01,
     *    0.399438638172067906D+01,
     *    0.500726073851736242D+01,
     *    0.111641442903222288D+02,
     *    0.263797551918497923D+02,
     *   -0.147978375673764972D+01,
     *   -0.192970042731921265D+01/
      PI = DACOS(-1.D0)
      CT=DCOS(T*PI/180.d0)
      DTERM = 0.D0
      LLEGD=1
      LLEGB=1
      LLEGG=6
      DO 80 L=1,LLEGD
      PLNORMD = DSQRT((2.D0*DFLOAT(L-1)+1.D0)/2.D0)
       DTERM = DTERM + XX(L)*PLNORMD*PLGNDR(L-1,0,CT)
80    CONTINUE
      BTERM = 0.D0
      DO 90 L=1,LLEGB
      PLNORMB = DSQRT((2.D0*DFLOAT(L-1)+1.D0)/2.D0)
       BTERM = BTERM + XX(LLEGD+L)*PLNORMB*PLGNDR(L-1,0,CT)
90    CONTINUE
      GTERM = 0.D0
      DO 100 L=1,LLEGG
      PLNORMG = DSQRT((2.D0*DFLOAT(2*L)+1.D0)/2.D0)
       GTERM = GTERM +
     * (XX(LLEGD+LLEGB+L)+
     * RR*XX(LLEGD+LLEGB+LLEGG+L)+
     * RR*RR*XX(LLEGD+LLEGB+LLEGG+LLEGG+L)+
     * RR*RR*RR*XX(LLEGD+LLEGB+LLEGG+LLEGG+LLEGG+L))*
     * PLNORMG*PLGNDR(2*L,1,CT)
100   CONTINUE
      LLEGC = LLEGD+LLEGB+LLEGG+LLEGG+LLEGG+LLEGG
       RAD3 = RR**3
       RAD4 = RR*RAD3
       RAD5 = RR*RAD4
       RAD6 = RR*RAD5
       RAD7 = RR*RAD6
       RAD8 = RR*RAD7
       RAD9 = RR*RAD8
       RAD10 = RR*RAD9
       RAD11 = RR*RAD10
       RAD12 = RR*RAD11
       RAD13 = RR*RAD12
C VAN DER WAALS PART
       re=0.7408D0
       CTERM5 = XX(LLEGC+1)/(RAD5)*D(5,BTERM,RR)*
     * DSQRT((2.D0*2.D0+1.D0)/2.D0)*PLGNDR(2,1,CT)
       CTERM6 = XX(LLEGC+2)/(RAD6)*D(6,BTERM,RR)*
     * DSQRT((2.D0*4.D0+1.D0)/2.D0)*PLGNDR(4,1,CT)
       CTERM7 = XX(LLEGC+3)/(RAD7)*D(7,BTERM,RR)*
     * DSQRT((2.D0*6.D0+1.D0)/2.D0)*PLGNDR(6,1,CT)
C INTERNAL COORDINATE
       ZVAR = (rsmall-re)/re
       VTERM = XX(LLEGC+4)+XX(LLEGC+5)*ZVAR+
     *         XX(LLEGC+6)*ZVAR**2+XX(LLEGC+7)*ZVAR**3
       VTERMEXP=DEXP(-XX(LLEGC+8)*ZVAR-XX(LLEGC+9)*ZVAR**2)
         CTERM = CTERM5+CTERM6+CTERM7
C FUNCTION VALUE

       POT3DH12 =(GTERM*DEXP(DTERM -RR*BTERM)*VTERMEXP
     *          + CTERM)*VTERM*4.556335d0

       RETURN
C LAST CARD IN POT3DH12 FUNCTION
       END

       FUNCTION D(nn,beta,r)
c
c     calculate the damping factor (small R correct)
c
      implicit real*8 (a-h,o-z)
      br=beta*r
      sum=1.0d0
      term=1.0d0
      ncn=nn
      do i=1,ncn
        term=term*br/i
        sum=sum+term
      enddo
      d=1.0d0 - dexp(-br)*sum
c     in case of d --> 0 use
c     d=1.0d0 - dexp(-br)*sum = sum_m=ncn+1^\infty br^m/m!
      if(dabs(d).lt.1.0d-8) then
        d=0.0d0
        do i=ncn+1,1000
          term=term*br/i
          d=d+term
          if(term/d .lt. 1.0d-8) go to 111
        enddo
        write(6,*) 'No convergence in d'
  111 continue
      d=d*dexp(-br)
      endif
      return
      end

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

      FUNCTION POT2DA1(R,rh2)
C FUNCTION FOR (R,r) DEPENDENCE OF
C F-H2 A1 STATE
C R,rh2 IN Angstrems
C OUTPUT IN MIKROHARTREE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension V0(9)
       CALL VlambdaA1(R,V0)
       re=0.74085D0
       z=(rh2-re)/re
       term1=V0(1)*z+VA1I(R,7)
       term2=V0(2)*z**2
       term3=V0(3)*z**3
       term4=V0(4)*z**4
       term5=V0(5)*z**5
       term6=V0(6)*z**6
       term7=V0(7)*z**7
       term8=V0(8)*z**8
       term9=V0(9)*z**9
       POT2DA1=term1+term2+term3+term4+
     * term5+term6+term7+term8+term9
       RETURN
      END
C
C
C
C
      subroutine VlambdaA1(R,V0)
      implicit double precision(a-h, o-z)
      dimension T(9,9),ipvt(9)
      dimension V0(9)
      pi=dacos(-1.d0)
      re=0.74085D0
      r1=0.42334D0
      r2=0.47626D0
      r3=0.52918D0
      r4=0.58210D0
      r5=0.63501D0
      r6=0.68793D0
      r7=re
      r8=0.79377D0
      r9=0.84668D0
      r10=0.89960D0
      t1=(r1-re)/re
      t2=(r2-re)/re
      t3=(r3-re)/re
      t4=(r4-re)/re
      t5=(r5-re)/re
      t6=(r6-re)/re
      t7=(r8-re)/re
      t8=(r9-re)/re
      t9=(r10-re)/re

      T(1,1)=t1
      T(1,2)=t1**2
      T(1,3)=t1**3
      T(1,4)=t1**4
      T(1,5)=t1**5
      T(1,6)=t1**6
      T(1,7)=t1**7
      T(1,8)=t1**8
      T(1,9)=t1**9
      T(2,1)=t2
      T(2,2)=t2**2
      T(2,3)=t2**3
      T(2,4)=t2**4
      T(2,5)=t2**5
      T(2,6)=t2**6
      T(2,7)=t2**7
      T(2,8)=t2**8
      T(2,9)=t2**9
      T(3,1)=t3
      T(3,2)=t3**2
      T(3,3)=t3**3
      T(3,4)=t3**4
      T(3,5)=t3**5
      T(3,6)=t3**6
      T(3,7)=t3**7
      T(3,8)=t3**8
      T(3,9)=t3**9
      T(4,1)=t4
      T(4,2)=t4**2
      T(4,3)=t4**3
      T(4,4)=t4**4
      T(4,5)=t4**5
      T(4,6)=t4**6
      T(4,7)=t4**7
      T(4,8)=t4**8
      T(4,9)=t4**9
      T(5,1)=t5
      T(5,2)=t5**2
      T(5,3)=t5**3
      T(5,4)=t5**4
      T(5,5)=t5**5
      T(5,6)=t5**6
      T(5,7)=t5**7
      T(5,8)=t5**8
      T(5,9)=t5**9
      T(6,1)=t6
      T(6,2)=t6**2
      T(6,3)=t6**3
      T(6,4)=t6**4
      T(6,5)=t6**5
      T(6,6)=t6**6
      T(6,7)=t6**7
      T(6,8)=t6**8
      T(6,9)=t6**9
      T(7,1)=t7
      T(7,2)=t7**2
      T(7,3)=t7**3
      T(7,4)=t7**4
      T(7,5)=t7**5
      T(7,6)=t7**6
      T(7,7)=t7**7
      T(7,8)=t7**8
      T(7,9)=t7**9
      T(8,1)=t8
      T(8,2)=t8**2
      T(8,3)=t8**3
      T(8,4)=t8**4
      T(8,5)=t8**5
      T(8,6)=t8**6
      T(8,7)=t8**7
      T(8,8)=t8**8
      T(8,9)=t8**9
      T(9,1)=t9
      T(9,2)=t9**2
      T(9,3)=t9**3
      T(9,4)=t9**4
      T(9,5)=t9**5
      T(9,6)=t9**6
      T(9,7)=t9**7
      T(9,8)=t9**8
      T(9,9)=t9**9
      V0(1)=VA1I(R,1)-VA1I(R,7)
      V0(2)=VA1I(R,2)-VA1I(R,7)
      V0(3)=VA1I(R,3)-VA1I(R,7)
      V0(4)=VA1I(R,4)-VA1I(R,7)
      V0(5)=VA1I(R,5)-VA1I(R,7)
      V0(6)=VA1I(R,6)-VA1I(R,7)
      V0(7)=VA1I(R,8)-VA1I(R,7)
      V0(8)=VA1I(R,9)-VA1I(R,7)
      V0(9)=VA1I(R,10)-VA1I(R,7)
      call dgesv(9,1,T,9,ipvt,V0,9,info)
C
      return
      end
C

      FUNCTION VA1I(R,N)
CSYSTEM: F-H2 A1
CLEVEL:RHF/UCCSD(T)/aug-cc-pvqz+332
C units R=A E=microEh
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(14),XX1(14),XX2(14),XX3(14),
     * XX4(14),XX5(14),XX6(14),XX7(14),XX8(14),
     * XX9(14),XX10(14)
C r=r1=0.42334 A
      DATA XX1/
     *  0.301418858448982929D+01,
     * -0.603705273584179025D+01,
     *  0.682915960700290307D+05,
     * -0.187666032260949782D+06,
     *  0.238301308334246220D+06,
     * -0.164244683880907105D+06,
     *  0.654163507393592809D+05,
     * -0.152247928229055105D+05,
     *  0.193569781326728071D+04,
     * -0.105288765445411940D+03,
     * -0.267482472856100245D-01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r2=0.47626 A
      DATA XX2/
     *  0.302063810440341030D+01,
     * -0.638833011929894568D+01,
     *  0.628740498834297978D+05,
     * -0.176096661305668618D+06,
     *  0.219680223613722919D+06,
     * -0.148219855928161065D+06,
     *  0.580248577947132362D+05,
     * -0.133253428905331384D+05,
     *  0.167675037673261022D+04,
     * -0.904316540418584083D+02,
     * -0.218200084608609485D-01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r3=0.52918 A
      DATA XX3/
     *  0.301685686644561724D+01,
     * -0.551775847247352313D+01,
     *  0.185856955473403039D+06,
     * -0.523721464309566189D+06,
     *  0.642625085935761104D+06,
     * -0.425782884896181815D+06,
     *  0.164142937796302431D+06,
     * -0.372112660932159197D+05,
     *  0.463005236787568083D+04,
     * -0.247007765851240379D+03,
     * -0.825900761957491514D-01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r4=0.58210 A
      DATA XX4/
     *  0.302848416668694886D+01,
     * -0.528421029366862793D+01,
     *  0.299444914273685019D+06,
     * -0.846202614471198991D+06,
     *  0.102471454895413003D+07,
     * -0.669979653675300884D+06,
     *  0.255722534873673751D+06,
     * -0.575729189852614654D+05,
     *  0.713228182899861167D+04,
     * -0.379708891873872460D+03,
     * -0.810470636018972429D-01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r5=0.63501 A
      DATA XX5/
     *  0.301903904174249416D+01,
     * -0.580166995255255191D+01,
     *  0.208814868503257632D+06,
     * -0.588130130126829026D+06,
     *  0.702723792203680612D+06,
     * -0.453320331730610400D+06,
     *  0.171024947143409197D+06,
     * -0.381150275527972772D+05,
     *  0.467918343905950132D+04,
     * -0.247067072552540992D+03,
     * -0.403090685500953963D-01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r6=0.68793 A
      DATA XX6/
     *  0.302492933202342806D+01,
     * -0.693082706820432204D+01,
     *  0.808806175727103546D+05,
     * -0.227090708047667606D+06,
     *  0.268582342694093648D+06,
     * -0.171656644661493512D+06,
     *  0.643261826310891774D+05,
     * -0.142737683146155960D+05,
     *  0.174882220895869500D+04,
     * -0.924922544120467336D+02,
     *  0.170244251061168403D-01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=re=0.74085 A
      DATA XX7/
     *  0.301980277505242789D+01,
     * -0.617615636080535335D+01,
     *  0.195263148647961119D+06,
     * -0.545196221413553925D+06,
     *  0.637964371158752940D+06,
     * -0.403625246676922659D+06,
     *  0.149959957821118820D+06,
     * -0.330330960155646899D+05,
     *  0.402097847767468920D+04,
     * -0.211193237785435201D+03,
     *  0.124209082604196828D-02,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r7=0.79377 A
      DATA XX8/
     *  0.301021338673093108D+01,
     * -0.624001964617968508D+01,
     *  0.202324783012565807D+06,
     * -0.561443567235070281D+06,
     *  0.650422295567868161D+06,
     * -0.407610947389787179D+06,
     *  0.150199530809161515D+06,
     * -0.328513252261865491D+05,
     *  0.397414152611111467D+04,
     * -0.207595239703559344D+03,
     *  0.282830888754052878D-02,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r8=0.84668 A
      DATA XX9/
     *  0.299847282996920583D+01,   
     * -0.655756498688726186D+01,
     *  0.159401107390797115D+06,
     * -0.439483762258756673D+06,
     *  0.504276070480560011D+06,
     * -0.313130611333813518D+06,
     *  0.114440752744809768D+06,
     * -0.248454294757396201D+05,
     *  0.298491272761102528D+04,
     * -0.154777916253962246D+03,
     * -0.203724640779095814D-01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r9=0.89960 A
      DATA XX10/
     *  0.380511452384088811D+01,
     * -0.908785189377794111D+01,
     * -0.122544933518701597D+06,
     *  0.458995191575567995D+06,
     * -0.732561192073940299D+06,
     *  0.652036493840829120D+06,
     * -0.355327266345501412D+06,
     *  0.121969156890437516D+06,
     * -0.258844174728603793D+05,
     *  0.311940518867429682D+04,
     * -0.164759837360351469D+03,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/

      IF(N.EQ.1) THEN
      DO I=1,14
      XX(I)=XX1(I)
      ENDDO
      ENDIF
      IF(N.EQ.2) THEN
      DO I=1,14
      XX(I)=XX2(I)
      ENDDO
      ENDIF
      IF(N.EQ.3) THEN
      DO I=1,14
      XX(I)=XX3(I)
      ENDDO
      ENDIF
      IF(N.EQ.4) THEN
      DO I=1,14
      XX(I)=XX4(I)
      ENDDO
      ENDIF
      IF(N.EQ.5) THEN
      DO I=1,14
      XX(I)=XX5(I)
      ENDDO
      ENDIF
      IF(N.EQ.6) THEN
      DO I=1,14
      XX(I)=XX6(I)
      ENDDO
      ENDIF
      IF(N.EQ.7) THEN
      DO I=1,14
      XX(I)=XX7(I)
      ENDDO
      ENDIF
      IF(N.EQ.8) THEN
      DO I=1,14
      XX(I)=XX8(I)
      ENDDO
      ENDIF
      IF(N.EQ.9) THEN
      DO I=1,14
      XX(I)=XX9(I)
      ENDDO
      ENDIF
      IF(N.EQ.10) THEN
      DO I=1,14
      XX(I)=XX10(I)
      ENDDO
      ENDIF


      TERM1 = DEXP(-XX(1)*R-XX(2))
      TERM2 = XX(3) + XX(4)*R+XX(5)*R**2+
     *        XX(6)*R**3+XX(7)*R**4+XX(8)*R**5
     *        +XX(9)*R**6+XX(10)*R**7+XX(11)*R**8
      TERM3=0.5D0*
     *   (1.D0+ DTANH(1.D0+XX(12)*R))
      TERM4 = XX(13)/(R**5)+XX(14)/(R**6)

       VA1I = (TERM1*TERM2-TERM3*TERM4)
       RETURN
       END
C
      FUNCTION POT2DB1(R,rh2)
C FUNCTION FOR (R,r) DEPENDENCE OF
C F-H2 B1 STATE
C R,rh2 IN Angstrems
C OUTPUT IN MIKROHARTREE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension V0(6)
       CALL VlambdaB1(R,V0)
       re=0.74085D0
       z=(rh2-re)/re
       term1=V0(1)*z+VB1I(R,4)
       term2=V0(2)*z**2
       term3=V0(3)*z**3
       term4=V0(4)*z**4
       term5=V0(5)*z**5
       term6=V0(6)*z**6
       POT2DB1=term1+term2+term3+term4
     * +term5+term6
       RETURN
      END
C
C
C
C
      subroutine VlambdaB1(R,V0)
      implicit double precision(a-h, o-z)
      dimension T(6,6),ipvt(6)
      dimension V0(6)
      pi=dacos(-1.d0)
      re=0.74085D0
      r1=0.42334D0
      r2=0.52918D0
      r3=0.63501D0
      r5=0.84668D0
      r6=0.95252D0
      r7=1.05835D0
      t1=(r1-re)/re
      t2=(r2-re)/re
      t3=(r3-re)/re
      t5=(r5-re)/re
      t6=(r6-re)/re
      t7=(r7-re)/re
      T(1,1)=t1
      T(1,2)=t1**2
      T(1,3)=t1**3
      T(1,4)=t1**4
      T(1,5)=t1**5
      T(1,6)=t1**6
      T(2,1)=t2
      T(2,2)=t2**2
      T(2,3)=t2**3
      T(2,4)=t2**4
      T(2,5)=t2**5
      T(2,6)=t2**6
      T(3,1)=t3
      T(3,2)=t3**2
      T(3,3)=t3**3
      T(3,4)=t3**4
      T(3,5)=t3**5
      T(3,6)=t3**6
      T(4,1)=t5
      T(4,2)=t5**2
      T(4,3)=t5**3
      T(4,4)=t5**4
      T(4,5)=t5**5
      T(4,6)=t5**6
      T(5,1)=t6
      T(5,2)=t6**2
      T(5,3)=t6**3
      T(5,4)=t6**4
      T(5,5)=t6**5
      T(5,6)=t6**6
      T(6,1)=t7
      T(6,2)=t7**2
      T(6,3)=t7**3
      T(6,4)=t7**4
      T(6,5)=t7**5
      T(6,6)=t7**6
      V0(1)=VB1I(R,1)-VB1I(R,4)
      V0(2)=VB1I(R,2)-VB1I(R,4)
      V0(3)=VB1I(R,3)-VB1I(R,4)
      V0(4)=VB1I(R,5)-VB1I(R,4)
      V0(5)=VB1I(R,6)-VB1I(R,4)
      V0(6)=VB1I(R,7)-VB1I(R,4)
      call dgesv(6,1,T,6,ipvt,V0,6,info)
C
      return
      end
C

      FUNCTION VB1I(R,N)
CSYSTEM: F-H2 B1 
CLEVEL:RHF/UCCSD(T)/aug-cc-pvqz+332
C units R=A E=microEh
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(14),XX1(14),XX2(14),XX3(14),
     * XX4(14),XX5(14),XX6(14),XX7(14)
C r=r1=0.42334 A
      DATA XX1/
     *  0.283357195276282159D+01,
     * -0.447996395278260007D+01,
     *  0.120474830125722481D+06,
     *  0.115015940999464419D+05,
     * -0.718971929144628084D+05,
     *  0.322034585517986015D+05,
     * -0.491576304350523515D+04,
     * -0.215086627651632853D+03,
     *  0.112713490885495432D+03,
     *  0.146634989562883871D-02,
     * -0.150112789126050239D+01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r2=0.52918 A
      DATA XX2/
     *  0.287565694372276637D+01,
     * -0.366724084316262822D+01,
     *  0.168661718634924589D+06,
     *  0.214034530736243818D+06,
     * -0.274929454983082716D+06,
     *  0.105687925716860424D+06,
     * -0.159634578179779619D+05,
     * -0.370851298572992278D+03,
     *  0.305974904741581042D+03,
     * -0.516782167836614750D-03,
     * -0.443307822670590834D+01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r3=0.63501 A
      DATA XX3/
     *  0.272357501448750794D+01,
     * -0.595235034416498454D+01,
     * -0.318407904728104040D+04,
     *  0.707250451042203640D+05,
     * -0.912179444956080115D+05,
     *  0.576819887385092807D+05,
     * -0.227460512564833298D+05,
     *  0.586382841600056418D+04,
     * -0.964449377738941280D+03,
     *  0.919700063136916839D+02,
     * -0.391895466702581707D+01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r4=0.74085 A
      DATA XX4/
     *  0.303641286457554482D+01,
     * -0.339987361360767126D+01,
     *  0.134927154486689280D+06,
     *  0.176307448709671007D+06,
     *  0.134958616753002862D+06,
     * -0.312351622227754153D+06,
     *  0.183894987760278251D+06,
     * -0.545020965378945184D+05,
     *  0.840349416828280118D+04,
     * -0.551336151607739453D+03,
     * -0.340617927152318768D-03,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r5=0.84668 A
      DATA XX5/
     *  0.314670355154628645D+01,
     * -0.421978869926994360D+01,
     *  0.380897949482160030D+02,
     *  0.203825493417480378D+06,
     * -0.799566701096517354D+05,
     *  0.396907951227411848D+00,
     * -0.691896869880620670D+04,
     *  0.997650545223350673D+04,
     * -0.425004450481778258D+04,
     *  0.803118691836126004D+03,
     * -0.606881967303535959D+02,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r6=0.95252 A
      DATA XX6/
     *  0.276653030787949072D+01,
     * -0.601197717883762461D+01,
     * -0.211874523650474330D+05,
     *  0.914183935555512144D+05,
     * -0.987666230306430953D+05,
     *  0.598198122497683580D+05,
     * -0.240408687476286286D+05,
     *  0.646133974122356358D+04,
     * -0.111220215202421628D+04,
     *  0.110682612636269241D+03,
     * -0.491308095624267693D+01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r7=1.05835 A
      DATA XX7/
     *  0.317856973781594565D+01,
     * -0.475694080460601221D+01,
     * -0.741846712837660471D+04,
     *  0.824400855359893321D+05,
     * -0.910671521711559184D+00,
     * -0.594514654130900908D+04,
     * -0.163809303451785072D+05,
     *  0.139203513135412231D+05,
     * -0.485591071112742793D+04,
     *  0.823138186948582870D+03,
     * -0.576795108934895779D+02,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/

        IF(N.EQ.1) THEN
        DO I=1,14
        XX(I)=XX1(I)
        ENDDO
        ENDIF
        IF(N.EQ.2) THEN
        DO I=1,14
        XX(I)=XX2(I)
        ENDDO
        ENDIF
        IF(N.EQ.3) THEN
        DO I=1,14
        XX(I)=XX3(I)
        ENDDO
        ENDIF
        IF(N.EQ.4) THEN
        DO I=1,14
        XX(I)=XX4(I)
        ENDDO
        ENDIF
        IF(N.EQ.5) THEN
        DO I=1,14
        XX(I)=XX5(I)
        ENDDO
        ENDIF
        IF(N.EQ.6) THEN
        DO I=1,14
        XX(I)=XX6(I)
        ENDDO
        ENDIF
        IF(N.EQ.7) THEN
        DO I=1,14
        XX(I)=XX7(I)
        ENDDO
        ENDIF

      TERM1 = DEXP(-XX(1)*R-XX(2))
      TERM2 = XX(3) + XX(4)*R+XX(5)*R**2+
     *        XX(6)*R**3+XX(7)*R**4+XX(8)*R**5
     *        +XX(9)*R**6+XX(10)*R**7+XX(11)*R**8
      TERM3=0.5D0*
     *   (1.D0+ DTANH(1.D0+XX(12)*R))
      TERM4 = XX(13)/(R**5)+XX(14)/(R**6)

       VB1I = (TERM1*TERM2-TERM3*TERM4)
       RETURN
       END
C
      FUNCTION POT2DB2(R,rh2)
C FUNCTION FOR (R,r) DEPENDENCE OF
C F-H2 B2 STATE
C R,rh2 IN Angstrems
C OUTPUT IN MIKROHARTREE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension V0(6)
       CALL VlambdaB2(R,V0)
       re=0.74085D0
       z=(rh2-re)/re
       term1=V0(1)*z+VB2I(R,4)
       term2=V0(2)*z**2
       term3=V0(3)*z**3
       term4=V0(4)*z**4
       term5=V0(5)*z**5
       term6=V0(6)*z**6
       POT2DB2=term1+term2+term3+term4
     * +term5+term6
       RETURN
      END
C
C
C
C
      subroutine VlambdaB2(R,V0)
      implicit double precision(a-h, o-z)
      dimension T(6,6),ipvt(6)
      dimension V0(6)
      pi=dacos(-1.d0)
      re=0.74085D0
      r1=0.42334D0
      r2=0.52918D0
      r3=0.63501D0
      r5=0.84668D0
      r6=0.95252D0
      r7=1.05835D0
      t1=(r1-re)/re
      t2=(r2-re)/re
      t3=(r3-re)/re
      t5=(r5-re)/re
      t6=(r6-re)/re
      t7=(r7-re)/re
      T(1,1)=t1
      T(1,2)=t1**2
      T(1,3)=t1**3
      T(1,4)=t1**4
      T(1,5)=t1**5
      T(1,6)=t1**6
      T(2,1)=t2
      T(2,2)=t2**2
      T(2,3)=t2**3
      T(2,4)=t2**4
      T(2,5)=t2**5
      T(2,6)=t2**6
      T(3,1)=t3
      T(3,2)=t3**2
      T(3,3)=t3**3
      T(3,4)=t3**4
      T(3,5)=t3**5
      T(3,6)=t3**6
      T(4,1)=t5
      T(4,2)=t5**2
      T(4,3)=t5**3
      T(4,4)=t5**4
      T(4,5)=t5**5
      T(4,6)=t5**6
      T(5,1)=t6
      T(5,2)=t6**2
      T(5,3)=t6**3
      T(5,4)=t6**4
      T(5,5)=t6**5
      T(5,6)=t6**6
      T(6,1)=t7
      T(6,2)=t7**2
      T(6,3)=t7**3
      T(6,4)=t7**4
      T(6,5)=t7**5
      T(6,6)=t7**6
      V0(1)=VB2I(R,1)-VB2I(R,4)
      V0(2)=VB2I(R,2)-VB2I(R,4)
      V0(3)=VB2I(R,3)-VB2I(R,4)
      V0(4)=VB2I(R,5)-VB2I(R,4)
      V0(5)=VB2I(R,6)-VB2I(R,4)
      V0(6)=VB2I(R,7)-VB2I(R,4)
      call dgesv(6,1,T,6,ipvt,V0,6,info)
C
      return
      end
C

      FUNCTION VB2I(R,N)
CSYSTEM: F-H2 B2
CLEVEL:RHF/UCCSD(T)/aug-cc-pvqz+332
C units R=A E=microEh
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(14),XX1(14),XX2(14),XX3(14),
     * XX4(14),XX5(14),XX6(14),XX7(14)
C r=r1=0.42334 A
      DATA XX1/
     *  0.295390872587323061D+01,
     * -0.453947093032389226D+01,
     *  0.118156228584792130D+06,
     *  0.224167464154384543D+05,
     * -0.690554894091893366D+05,
     *  0.266362877211812229D+05,
     * -0.355625183723824512D+04,
     * -0.612217862939555546D-02,
     * -0.479362034306542313D+02,
     *  0.293949486242627032D+02,
     * -0.351786738641760000D+01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/

C r=r2=0.52918 A
      DATA XX2/
     *  0.275155005023126487D+01,
     * -0.615709853193104717D+01,
     *  0.945635569747667978D+04,
     *  0.341928142287404407D+05,
     * -0.520340841086358196D+05,
     *  0.333593620514415379D+05,
     * -0.128791165468607669D+05,
     *  0.321851674394995507D+04,
     * -0.512455089089091530D+03,
     *  0.473214849876440127D+02,
     * -0.195438603790498422D+01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/

C r=r3=0.63501 A
      DATA XX3/
     *  0.333725803136194044D+01,
     * -0.267880823655760159D+01,
     *  0.596950710992069566D+06,
     *  0.163943950042255544D+02,
     *  0.867914629123724997D+06,
     * -0.981230853226367501D+06,
     *  0.450003677472476382D+06,
     * -0.103125263823762230D+06,
     *  0.814011750871661661D+04,
     *  0.934650894937899579D+03,
     * -0.174645165794070209D+03,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/

C r=r4=0.74085 A
      DATA XX4/
     *  0.287466199339493267D+01,
     * -0.740940376009246826D+01,
     * -0.219249015717386010D+04,
     *  0.197722690470349771D+05,
     * -0.226850242782013411D+05,
     *  0.136138783780608283D+05,
     * -0.526133850926898958D+04,
     *  0.133656314205604372D+04,
     * -0.214818339537281190D+03,
     *  0.197390619170379082D+02,
     * -0.801984425698802683D+00,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/

C r=r5=0.84668 A
      DATA XX5/
     *  0.344635486563916960D+01,
     * -0.226638385192963199D+01,
     *  0.124263179003213765D+06,
     *  0.107743554458265635D+07,
     *  0.867351491690326133D+06,
     * -0.133065270364417532D+07,
     *  0.697818993751274538D+06,
     * -0.195084434963865147D+06,
     *  0.243627862149680877D+05,
     *  0.111579050816088432D+01,
     * -0.233502575631271128D+03,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/

C r=r6=0.95252 A
      DATA XX6/
     *  0.334637810745046593D+01,
     * -0.687657540987928062D+01,
     * -0.817659920031785987D+04,
     *  0.352538520695417465D+05,
     * -0.240665573384614945D+05,
     *  0.959739155079603006D+04,
     * -0.220114290909493639D+04,
     *  0.866968678060150069D+02,
     *  0.575741998160614372D+02,
     * -0.889339266248669880D+01,
     * -0.655469724938011173D-04,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/

C r=r7=1.05835 A
      DATA XX7/
     *  0.349595516987215493D+01,
     * -0.371815832436441784D+01,
     *  0.301121766462376783D+02,
     * -0.115853613700391384D+04,
     *  0.778787608245793148D+06,
     * -0.872544104222255759D+06,
     *  0.506782723350841203D+06,
     * -0.176643484336668451D+06,
     *  0.355901030709365878D+05,
     * -0.375907062809635636D+04,
     *  0.143508214783077364D+03,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/

        IF(N.EQ.1) THEN
        DO I=1,14
        XX(I)=XX1(I)
        ENDDO
        ENDIF
        IF(N.EQ.2) THEN
        DO I=1,14
        XX(I)=XX2(I)
        ENDDO
        ENDIF
        IF(N.EQ.3) THEN
        DO I=1,14
        XX(I)=XX3(I)
        ENDDO
        ENDIF
        IF(N.EQ.4) THEN
        DO I=1,14
        XX(I)=XX4(I)
        ENDDO
        ENDIF
        IF(N.EQ.5) THEN
        DO I=1,14
        XX(I)=XX5(I)
        ENDDO
        ENDIF
        IF(N.EQ.6) THEN
        DO I=1,14
        XX(I)=XX6(I)
        ENDDO
        ENDIF
        IF(N.EQ.7) THEN
        DO I=1,14
        XX(I)=XX7(I)
        ENDDO
        ENDIF

      TERM1 = DEXP(-XX(1)*R-XX(2))
      TERM2 = XX(3) + XX(4)*R+XX(5)*R**2+
     *        XX(6)*R**3+XX(7)*R**4+XX(8)*R**5
     *        +XX(9)*R**6+XX(10)*R**7+XX(11)*R**8
      TERM3=0.5D0*
     *   (1.D0+ DTANH(1.D0+XX(12)*R))
      TERM4 = XX(13)/(R**5)+XX(14)/(R**6)

       VB2I = (TERM1*TERM2-TERM3*TERM4)
       RETURN
       END
C
      FUNCTION POT2DPI(R,rh2)
C FUNCTION FOR (R,r) DEPENDENCE OF
C F-H2 PI STATE
C R,rh2 IN Angstrems
C OUTPUT IN MIKROHARTREE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension V0(5)
       CALL VlambdaPI(R,V0)
       re=0.74085D0
       z=(rh2-re)/re
       term1=V0(1)*z+VPII(R,3)
       term2=V0(2)*z**2
       term3=V0(3)*z**3
       term4=V0(4)*z**4
       term5=V0(5)*z**5
       POT2DPI=term1+term2+term3
     *         +term4+term5 
       RETURN
      END
C
C
C
C
      subroutine VlambdaPI(R,V0)
      implicit double precision(a-h, o-z)
      dimension T(5,5),ipvt(5)
      dimension V0(5)
      pi=dacos(-1.d0)
      re=0.74085D0
      r1=0.52918D0
      r2=0.63501D0
      r3=0.84668D0
      r4=0.95252D0
      r5=1.05835D0
      t1=(r1-re)/re
      t2=(r2-re)/re
      t3=(r3-re)/re
      t4=(r4-re)/re
      t5=(r5-re)/re
      T(1,1)=t1
      T(1,2)=t1**2
      T(1,3)=t1**3
      T(1,4)=t1**4
      T(1,5)=t1**5
      T(2,1)=t2
      T(2,2)=t2**2
      T(2,3)=t2**3
      T(2,4)=t2**4
      T(2,5)=t2**5
      T(3,1)=t3
      T(3,2)=t3**2
      T(3,3)=t3**3
      T(3,4)=t3**4
      T(3,5)=t3**5
      T(4,1)=t4
      T(4,2)=t4**2
      T(4,3)=t4**3
      T(4,4)=t4**4
      T(4,5)=t4**5
      T(5,1)=t5
      T(5,2)=t5**2
      T(5,3)=t5**3
      T(5,4)=t5**4
      T(5,5)=t5**5



      V0(1)=VPII(R,1)-VPII(R,3)
      V0(2)=VPII(R,2)-VPII(R,3)
      V0(3)=VPII(R,4)-VPII(R,3)
      V0(4)=VPII(R,5)-VPII(R,3)
      V0(5)=VPII(R,6)-VPII(R,3)
      call dgesv(5,1,T,5,ipvt,V0,5,info)
C
      return
      end
C

      FUNCTION VPII(R,N)
CSYSTEM: F-H2 PI
CLEVEL:RHF/UCCSD(T)/aug-cc-pvqz+332
C units R=A E=microEh
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(14),XX1(14),XX2(14),XX3(14),
     * XX4(14),XX5(14),XX6(14)
C r=r1=0.52918 A
      DATA XX1/
     *   0.288101387672307796D+01,
     *  -0.354829457167905948D+01,
     *   0.291136354246576608D+06,
     *   0.369464696956159969D+05,
     *   0.237963337296140143D+01,
     *  -0.141753258082832588D+06,
     *   0.104161230265923543D+06,
     *  -0.333158711217429227D+05,
     *   0.519580324533040766D+04,
     *  -0.321386881706345662D+03,
     *  -0.241476449572259400D+01,
     *  -0.289573768192804302D+02,
     *   0.876525683507813606D+06,
     *   0.387887406511697778D+07/

C r=r2=0.63501 A
      DATA XX2/
     *   0.242584637231939126D+01,
     *  -0.445927465567822878D+01,
     *   0.441012349872669001D+05,
     *   0.167905792207280698D+06,
     *  -0.261928162331646425D+06,
     *   0.163923933653356362D+06,
     *  -0.593695738761498578D+05,
     *   0.136318252623030912D+05,
     *  -0.198993785386692002D+04,
     *   0.170439090478737825D+03,
     *  -0.665517431122651004D+01,
     *  -0.289573768192804302D+02,
     *   0.876525683507813606D+06,
     *   0.387887406511697778D+07/

C r=r3=0.74085 A
      DATA XX3/
     *   0.270604825982857111D+01,
     *  -0.379856573477150539D+01,
     *   0.219345198441248707D+06,
     *  -0.942320080548216938D+05,
     *   0.159497867944249418D+06,
     *  -0.178063998875649500D+06,
     *   0.837444404448538407D+05,
     *  -0.192982070065093685D+05,
     *   0.192860335836597301D+04,
     *   0.111114209799324945D+00,
     *  -0.107682709499718499D+02,
     *  -0.289573768192804302D+02,
     *   0.876525683507813606D+06,
     *   0.387887406511697778D+07/

C r=r4=0.84668 A
      DATA XX4/
     *   0.295616710734156696D+01,
     *  -0.582731330195263908D+01,
     *   0.119370392994385213D+06,
     *  -0.257451496362622100D+06,
     *   0.304612783285240177D+06,
     *  -0.186315242876922130D+06,
     *   0.599656510847001628D+05,
     *  -0.902007764541772121D+04,
     *  -0.399783258697925872D-01,
     *   0.175809304071151104D+03,
     *  -0.166037450690838355D+02,
     *  -0.289573768192804302D+02,
     *   0.876525683507813606D+06,
     *   0.387887406511697778D+07/

C r=r5=0.95252 A
      DATA XX5/
     *   0.249851276811401180D+01,
     *  -0.365805308591928435D+01,
     *   0.168273506377212878D+06,
     *  -0.525591948786162684D+05,
     *   0.144421037511982751D+06,
     *  -0.165331041044362937D+06,
     *   0.738001964887552749D+05,
     *  -0.157905495102099339D+05,
     *   0.144879681491836209D+04,
     *   0.589935257060922175D-01,
     *  -0.661101416330000635D+01,
     *  -0.289573768192804302D+02,
     *   0.876525683507813606D+06,
     *   0.387887406511697778D+07/

C r=r6=1.05835 A
      DATA XX6/
     *   0.255923875968635617D+01,
     *  -0.533068068602461409D+01,
     *   0.107771077313168193D+06,
     *  -0.219044839917762205D+06,
     *   0.265874971309535264D+06,
     *  -0.173347095176382485D+06,
     *   0.630010905801819827D+05,
     *  -0.131205994379011136D+05,
     *   0.148231350114034808D+04,
     *  -0.722857133485914716D+02,
     *   0.792093501325373955D-01,
     *  -0.289573768192804302D+02,
     *   0.876525683507813606D+06,
     *   0.387887406511697778D+07/



        IF(N.EQ.1) THEN
        DO I=1,14
        XX(I)=XX1(I)
        ENDDO
        ENDIF
        IF(N.EQ.2) THEN
        DO I=1,14
        XX(I)=XX2(I)
        ENDDO
        ENDIF
        IF(N.EQ.3) THEN
        DO I=1,14
        XX(I)=XX3(I)
        ENDDO
        ENDIF
        IF(N.EQ.4) THEN
        DO I=1,14
        XX(I)=XX4(I)
        ENDDO
        ENDIF
        IF(N.EQ.5) THEN
        DO I=1,14
        XX(I)=XX5(I)
        ENDDO
        ENDIF
        IF(N.EQ.6) THEN
        DO I=1,14
        XX(I)=XX6(I)
        ENDDO
        ENDIF



      TERM1 = DEXP(-XX(1)*R-XX(2))
      TERM2 = XX(3) + XX(4)*R+XX(5)*R**2+
     *        XX(6)*R**3+XX(7)*R**4+XX(8)*R**5
     *        +XX(9)*R**6+XX(10)*R**7+XX(11)*R**8
      TERM3=0.5D0*
     *   (1.D0+ DTANH(1.D0+XX(12)*R))
      TERM4 = XX(13)/(R**5)+XX(14)/(R**6)

       VPII = (TERM1*TERM2-TERM3*TERM4)
       RETURN
       END
C
      FUNCTION POT2DSG(R,rh2)
C FUNCTION FOR (R,r) DEPENDENCE OF
C F-H2 SIGMA STATE
C R,rh2 IN Angstrems
C OUTPUT IN MIKROHARTREE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension V0(3)
       CALL VlambdaSG(R,V0)
       re=0.74085D0
       z=(rh2-re)/re
       term1=V0(1)*z+VSGI(R,3)
       term2=V0(2)*z**2
       term3=V0(3)*z**3
       POT2DSG=term1+term2+term3
       RETURN
      END
C
C
C
C
      subroutine VlambdaSG(R,V0)
      implicit double precision(a-h, o-z)
      dimension T(3,3),ipvt(3)
      dimension V0(3)
      pi=dacos(-1.d0)
      re=0.74085D0
      r1=0.52918D0
      r2=0.63501D0
      r3=0.84668D0
      t1=(r1-re)/re
      t2=(r2-re)/re
      t3=(r3-re)/re
      T(1,1)=t1
      T(1,2)=t1**2
      T(1,3)=t1**3
      T(2,1)=t2
      T(2,2)=t2**2
      T(2,3)=t2**3
      T(3,1)=t3
      T(3,2)=t3**2
      T(3,3)=t3**3
      V0(1)=VSGI(R,1)-VSGI(R,3)
      V0(2)=VSGI(R,2)-VSGI(R,3)
      V0(3)=VSGI(R,4)-VSGI(R,3)
      call dgesv(3,1,T,3,ipvt,V0,3,info)
C
      return
      end
C

      FUNCTION VSGI(R,N)
CSYSTEM: F-H2 SIGMA
CLEVEL:RHF/UCCSD(T)/aug-cc-pvqz+332
C units R=A E=microEh
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(14),XX1(14),XX2(14),XX3(14),
     * XX4(14)
C r=r1=0.52918 A
      DATA XX1/
     *   0.347621098314397114D+01,
     *  -0.798710130301622012D+01,
     *  -0.467617238842874722D+05,
     *   0.491079091587328803D+05,
     *  -0.855869016233460206D+04,
     *  -0.101551258834595355D+05,
     *   0.477777240120719489D+04,
     *  -0.389575385098098886D-01,
     *  -0.487570912230031070D+03,
     *   0.125758103534125368D+03,
     *  -0.110594558700220453D+02,
     *   0.768592052646715018D+01,
     *   0.236725079764586641D+06,
     *  -0.173346530057967524D+07/

C r=r2=0.63501 A
      DATA XX2/
     *   0.333151909247473288D+01,
     *  -0.902874771605496029D+01,
     *  -0.418388904657136954D+05,
     *   0.538730352252237644D+05,
     *  -0.238232960317764409D+05,
     *  -0.612759984206418380D+03,
     *   0.404250979379386172D+04,
     *  -0.143924324941787609D+04,
     *   0.185531786841853489D+03,
     *  -0.235510589418366099D-05,
     *  -0.167298839352105655D+01,
     *   0.435400419555694447D+06,
     *   0.594406277535313508D+06,
     *  -0.432590569201788213D+07/

C r=r3=0.74085 A
      DATA XX3/
     *  0.246118595428434483D+01,
     * -0.106225896220179212D+02,
     * -0.144803130912030902D+05,
     *  0.297159808611018743D+05,
     * -0.280533690955581515D+05,
     *  0.154955056128223841D+05,
     * -0.545145928921563427D+04,
     *  0.124588772032276438D+04,
     * -0.180675232499912738D+03,
     *  0.152161667237378868D+02,
     * -0.575821782221406453D+00,
     *  0.435400419555694447D+06,
     * -0.486618092994658102D+06,
     * -0.718598600303590298D+07/

C r=r4=0.84668 A
      DATA XX4/
     *  0.340670754329443470D+01,
     * -0.144041266721032706D+02,
     *  0.747437940712431828D+04,
     * -0.190768193189354497D+05,
     *  0.223218687851440191D+05,
     * -0.151890424974961279D+05,
     *  0.654937928075852597D+04,
     * -0.182777406838075399D+04,
     *  0.323214681866157889D+03,
     * -0.332235709696224362D+02,
     *  0.154443880986455384D+01,
     *  0.435400419555694447D+06,
     * -0.542820067369836010D+07,
     *  0.393913160980229229D+08/

        IF(N.EQ.1) THEN
        DO I=1,14
        XX(I)=XX1(I)
        ENDDO
        ENDIF
        IF(N.EQ.2) THEN
        DO I=1,14
        XX(I)=XX2(I)
        ENDDO
        ENDIF
        IF(N.EQ.3) THEN
        DO I=1,14
        XX(I)=XX3(I)
        ENDDO
        ENDIF
        IF(N.EQ.4) THEN
        DO I=1,14
        XX(I)=XX4(I)
        ENDDO
        ENDIF

      TERM1 = DEXP(-XX(1)*R-XX(2))
      TERM2 = XX(3) + XX(4)*R+XX(5)*R**2+
     *        XX(6)*R**3+XX(7)*R**4+XX(8)*R**5
     *        +XX(9)*R**6+XX(10)*R**7+XX(11)*R**8
      TERM3=0.5D0*
     *   (1.D0+ DTANH(1.D0+XX(12)*R))
      TERM4 = XX(13)/(R**5)+XX(14)/(R**6)

       VSGI = (TERM1*TERM2-TERM3*TERM4)
       RETURN
       END



      
C===================================================
c  The routine supplies the value of the diabatic 
c  electronic states for the F + H2 system as a function
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
c====================================================

       subroutine diabats(R,theta,V1,V2,V3,V12)
C===================================================
c  il programma calcola gli stati elettronici diabatici
c  per il sistema F + H2 in funzione di R di
c  Jacobi e dell'angolo theta 
C===================================================
      IMPLICIT REAL*8 (A-H,O-Z)                                
      COMMON/BLK1/EPSI0,RM0,BETA,C000,X3,X4,B1,B2,B3,B4    
      COMMON/BLK2/X1,X2,A1,A2,A3,A4,AEXP,ALFA1       
      COMMON/BLK3/A202,ALFA202,C202
      COMMON/BLK4/C022,A022,ALFA022,GAMMA022           
      COMMON/BLK5/A220,ALFA220      
      COMMON/BLK6/A222,ALFA222      
      COMMON/BLK7/C224      
      parameter(pigreco = 3.1415926535897932)

c parametri termini radiali

c V000 (MSV)

	RM0 = 3.3400d0                                           
	EPSI0= 4.120d0                                
	BETA = 6.300D0                                  
	C000 = 7516.0d0                               
 	AEXP= 0.0000D0                               
 	ALFA1=0.00000D0                                    
c	AEXP= 0.3800D0                               
c	ALFA1=17.00000D0                                    
c	X1 = 0.6900D0                                  
c	X2 = 0.8600D0                                  
 	X1 = 0.0000D0                                  
 	X2 = 0.0000D0                                  
 	X3 = 1.1100D0                                            
 	X4 = 1.5000D0   

c V202
c	A202=2.0d5
c	ALFA202=4.000d0

 	A202=1.3d5
 	ALFA202=3.500d0

	C202=C000*0.0928d0

c V022
	C022=902.d0 
	A022=8.1149d5
	ALFA022=3.5d0
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
	C224=0.49d0*1.544d0*1129.0d0


      CALL SPLINE_PES                    
c     write(*,11) b1,b2,b3,b4
11    format(4f11.4)

	c1=5.d0*dsqrt(5.0d0)
	c2=2.0d0*dsqrt(2.0d0)/(5.0d0*dsqrt(7.0d0))
	c3=6.0d0*dsqrt(2.0d0)/(5.0d0*dsqrt(35.d0))
	c4=9.0d0/(10.0d0*dsqrt(70.0d0))

C=======chiamate alle routines che calcolano i coeff. Vabc===
 	call v000(r,p000)
	call v202(r,p202)
 	call v022(r,p022)
 	call v220(r,p220)
 	call v222(r,p222)
 	call v224(r,p224)

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


