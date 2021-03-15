*Ar-NO(A^2Sigma 3ssigma) PES RCCSD(T)/AVTZ+33221 scaled with 1.23 factor
* Reference:
* J. Klos, M. H. Alexander, R. Hernandez-Lamoneda and T. G. Wright
* J. Chem. Phys. 129, 244303 (2008)

      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(10)
      include "common/parpot"
      potnam='Klos et al Ar-NO(A) CCSDT PES scaled'
      print *, potnam
1      print *, ' r (bohr)'
      rshift=0.5
      xfact=0.8
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vsig',/,11(1pe16.8))
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
      potnam='Klos et al Ar-NO(A^2Sigma) CCSDT PES with 1.23 fact'
      lammin(1)=1
      lammax(1)=10
      mproj(1)=0
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
*    vvl(1-10) expansion coefficients in dl0 (l=1:11) of vsum

* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  8-oct-1993
* revised for He-NO(X) : 1-20-95 by Moonbong Yang
* revised for CCSD(T) PES: 2002.10.13 by Jacek Klos
*revised for Ar-OH(A^2Sigma+) 2005 J.Klos
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension xlam1(7),xlam2(7),r0(7),
     :          vsum(11),xsum(11),
     :          d0(121),aa(121),thta(11),cthta(11)
      dimension kpvt(11),qraux(11),work(121),rsd(11)

      common /covvl/ vvl(10)
      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
* for distances beyond rmax diffeeence potential is damped
      data rmax /13d0/
* coefficicients for d0 rotation matrices
* stored (by column) for each of 7 angles and for l=0:6
* angles are 0 20 40 60 80 90 100 120 140 160 180
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
      do i=1,11
      cthta(i)=dcos(thta(i)*dacos(-1.D0)/180.D0)
      enddo
      do 100 i=1,11
        vsum(i)=VPES(r,thta(i))
100   continue
* solve simultaneous equations for solutions
      tol=1.e-10
      call dcopy(121,d0,1,aa,1)
      call dqrank(aa,11,11,11,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,11,11,11,kr,vsum,xsum,rsd,kpvt,qraux)
* convert to hartree
      conv=1.d0/219474.6
      call dscal(11,conv,xsum,1)
      vv0=xsum(1)
      call dcopy(10,xsum(2),1,vvl,1)
      end

      FUNCTION VPES(R,theta)
C*********************************
C System: Ar-NO(A2Sigma 3ssigma )
C Method:RCCSD(T)
C Basis:aug-cc-pvtz+3s3p2d2f1g1h
C CP and BSSE corrected
C Variables:
C R in Bohr
C Theta n Degrees
C Output in cm-1
C Theta=0 for Ar---NO geom.
C********************************
C Authors:
C Dr Jacek Klos
C Department of Chemistry and Biochemistry
C University of Maryland
C College Park, MD 20742
C email:jklos@mail.umd.edu
C 
C Dr Ramon Hernandez
C Universidad Autonoma de Estado de Morelos
C Dr Tim Wright
C University of Nottingham
C**********************************
C Needs link with LAPACK
C*********************************
      implicit double precision(a-h, o-z)
      dimension V0(7)
      dimension T(7)
      pi=dacos(-1.d0)
      call Vlpes(R,V0)
       do j=1,7
       T(j)=PLGNDR((j-1),0,DCOS(theta*pi/180.D0))
       enddo
       s=0.D0
       do i=1,7
       s=s+V0(i)*T(i)
       enddo
       SCALEFACT=1.23D0
       VPES=s*SCALEFACT
       return
       end


      subroutine Vlpes(R,V0)
      implicit double precision(a-h, o-z)
      dimension T(7,7),ipvt(7),work(100)
      dimension V0(7),theta(7)
      LWORK=100
      pi=dacos(-1.d0)
      theta(1)=0.D0
      theta(2)=30.D0
      theta(3)=60.D0
      theta(4)=90.D0
      theta(5)=120.D0
      theta(6)=150.D0
      theta(7)=180.D0
      do i=1,7
       do j=1,7
       T(i,j)=PLGNDR((j-1),0,DCOS(theta(i)*pi/180.D0)) 
       enddo
      enddo
      do k=1,7
      V0(k)=VPESN(R,k)
      enddo
      call dgesv(7,1,T,7,ipvt,V0,7,info)
c      call dgels('N',11,9,1,T,11,V0,11,WORK,LWORK,info)
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



      FUNCTION VPESN(R,N)
CSYSTEM: Ar-NO(Sigma) 2ssigma Rydberg state Theta=0 for Ar-NO
CLEVEL:RHF/RCCSD(T)/aug-cc-pvtz+332211
C units R=au E=microEh(parameters)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(15),XX1(15),XX2(15),XX3(15),
     * XX4(15),XX5(15),XX6(15),XX7(15),
     * XX8(15),XX9(15),XX10(15),XX11(15)
C theta=0
      DATA XX1/
     *    .189675783208241233D+01,
     *   -.309195004777979143D+01,
     *   -.116531365637548618D+11,
     *    .957399371781109810D+10,
     *   -.270115688871362400D+10,
     *    .131225730419954672D+09,
     *    .986307719193495959D+08,
     *   -.261243938556781746D+08,
     *    .300572119932373473D+07,
     *   -.174007436735790980D+06,
     *    .418957467898137111D+04,
     *    .904264556193752966D+01,
     *    .364105970166922271D+09,
     *   -.141607057842372894D+10,
     *    .412195477942114182D+11/

C theta=30
      DATA XX2/
     *    .173107296604431515D+01,
     *    .737379343308677495D-04,
     *   -.487527598170212364D+10,
     *    .647452421583293629D+10,
     *   -.230528142992037106D+10,
     *    .185762863381997764D+09,
     *    .855809970015178770D+08,
     *   -.272566888508405536D+08,
     *    .351143216033564880D+07,
     *   -.223737637863645534D+06,
     *    .601940287135423932D+04,
     *    .904264556193752966D+01,
     *    .399563417748368144D+09,
     *   -.313457283677237797D+10,
     *    .697753079412517548D+11/

C theta=60
      DATA XX3/       
     *    .181110708541646059D+01,
     *   -.183515540527415862D-02,
     *   -.578180269769829845D+10,
     *    .111908179067941456D+11,
     *   -.506537471001978874D+10,
     *    .749487441108477473D+09,
     *    .912560717300215811D+08,
     *   -.488702694296826720D+08,
     *    .737954670387391467D+07,
     *   -.524456385242260993D+06,
     *    .155184538731292978D+05,
     *    .904264556193752966D+01,
     *    .446758503045856476D+09,
     *   -.490250745743442917D+10,
     *    .970570252295909119D+11/

C theta=90
      DATA XX4/       
     *    .176352817561121267D+01,
     *    .168184400617172859D-03,
     *    .188614417127255583D+10,
     *    .324998677034662390D+10,
     *   -.181634740407821918D+10,
     *    .208684516734100252D+09,
     *    .731053175089619607D+08,
     *   -.270482666059069224D+08,
     *    .380649017730596568D+07,
     *   -.262700418548160640D+06,
     *    .780635903118705937D+04,
     *    .904264556193752966D+01,
     *    .495355533432190239D+09,
     *   -.716921229136979103D+10,
     *    .133839401140674118D+12/

C theta=120
      DATA XX5/       
     *    .180601660651678153D+01,
     *    .316942778923712221D-03,
     *    .107778260058283787D+11,
     *   -.491758807447486687D+10,
     *    .142360015381824923D+10,
     *   -.552982434821508288D+09,
     *    .203771391833207130D+09,
     *   -.454305627726987228D+08,
     *    .579085608980125096D+07,
     *   -.397265125641557621D+06,
     *    .120254676811789905D+05,
     *    .904264556193752966D+01,
     *    .468696940920499682D+09,
     *   -.558171343202733707D+10,
     *    .112859762677947037D+12/

C theta=150
      DATA XX6/
     *    .177281436023091477D+01,
     *   -.228384099665141532D+00,
     *    .146563745478028154D+10,
     *    .283565668014795446D+10,
     *   -.166489885698309255D+10,
     *    .207670123641992271D+09,
     *    .688149277762462497D+08,
     *   -.269004237529385984D+08,
     *    .387283941660520900D+07,
     *   -.269204185034881637D+06,
     *    .785643570293914217D+04,
     *    .904264556193752966D+01,
     *    .446979793798343599D+09,
     *   -.475405897638587475D+10,
     *    .104187333074331604D+12/

C theta=180
      DATA XX7/
     *    .178920621032253790D+01,
     *    .253316893964572853D+00,
     *    .174283618682380257D+11,
     *   -.881477420972724342D+10,
     *    .214695101520432377D+10,
     *   -.572593284989805222D+09,
     *    .206461813950128168D+09,
     *   -.494426655049585029D+08,
     *    .657311348627357371D+07,
     *   -.455796041772730707D+06,
     *    .135340370295620214D+05,
     *    .904264556193752966D+01,
     *    .405092790907102644D+09,
     *   -.289843999447082996D+10,
     *    .816031187896175079D+11/


      IF(N.EQ.1) THEN
      DO I=1,15
      XX(I)=XX1(I)
      ENDDO
      ENDIF
      IF(N.EQ.2) THEN
      DO I=1,15
      XX(I)=XX2(I)
      ENDDO
      ENDIF
      IF(N.EQ.3) THEN
      DO I=1,15
      XX(I)=XX3(I)
      ENDDO
      ENDIF
      IF(N.EQ.4) THEN
      DO I=1,15
      XX(I)=XX4(I)
      ENDDO
      ENDIF
      IF(N.EQ.5) THEN
      DO I=1,15
      XX(I)=XX5(I)
      ENDDO
      ENDIF
      IF(N.EQ.6) THEN
      DO I=1,15
      XX(I)=XX6(I)
      ENDDO
      ENDIF
      IF(N.EQ.7) THEN
      DO I=1,15
      XX(I)=XX7(I)
      ENDDO
      ENDIF

      TERM1 = DEXP(-XX(1)*R-XX(2))
      TERM2 = XX(3) + XX(4)*R+XX(5)*R**2+
     *        XX(6)*R**3+XX(7)*R**4+XX(8)*R**5
     *        +XX(9)*R**6+XX(10)*R**7+XX(11)*R**8
      TERM3=0.5D0*
     *   (1.D0+ DTANH(1.D0+XX(12)*R))
      TERM4 = XX(13)/(R**6)+XX(14)/(R**7)+XX(15)/(R**8)
       TOCM=0.219474643545745D0
       VPESN = (TERM1*TERM2-TERM3*TERM4)*TOCM
       RETURN
       END

