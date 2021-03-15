* Ar-OH(A^2Sigma+) PES, r=re of OH(A),  RCCSD(T) Klos et al
* Reference
* J. Klos, M. H. Alexander, M. Brouard, C. J. Eyles and F. J. Aoiz
* J. Chem. Phys.  129, 054301 (2008)

      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(10)
      include "common/parpot"
      potnam='Klos et al Ar-OH(A^2Sigma+) CCSDT PES'
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
      potnam='Klos et al Ar-OH(A^2Sigma+) CCSDT PES'
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
      dimension kpvt(11),qraux(11),work(121),rsd(7),re(14)

      common /covvl/ vvl(10)
      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
* for distances beyond rmax difference potential is damped
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
C System: Ar-OH(A2Sigma+) r=re
C Method:RCCSD(T)
C Basis:aug-cc-pvqz+3s3p2d2f1g
C CP and BSSE corrected
C Variables:
C R in Bohr
C Theta n Degrees
C Output in cm-1
C Theta=0 for Ar---OH geom.
C********************************
C Authors:
C Dr Jacek Klos
C Department of Chemistry and Biochemistry
C University of Maryland
C College Park, MD 20742
C email:jklos@mail.umd.edu
C**********************************
C Needs link with LAPACK
C*********************************
      implicit double precision(a-h, o-z)
      dimension V0(9)
      dimension T(9)
      pi=dacos(-1.d0)
      call Vlsum(R,V0)
       do j=1,9
       T(j)=PLGNDR((j-1),0,DCOS(theta*pi/180.D0))
       enddo
       s=0.D0
       do i=1,9
       s=s+V0(i)*T(i)
       enddo
       VPES=s
       return
       end

      subroutine Vlsum(R,V0)
      implicit double precision(a-h, o-z)
      dimension T(9,9),ipvt(9),work(100)
      dimension V0(9),theta(9)
      LWORK=100
      pi=dacos(-1.d0)
      theta(1)=0.D0
      theta(2)=22.5D0
      theta(3)=45.D0
      theta(4)=67.5D0
      theta(5)=90.D0
      theta(6)=112.5D0
      theta(7)=135.D0
      theta(8)=157.5D0
      theta(9)=180.D0
      do i=1,9
       do j=1,9
       T(i,j)=PLGNDR((j-1),0,DCOS(theta(i)*pi/180.D0)) 
       enddo
      enddo
      do k=1,9
      V0(k)=VAPRIMI(R,k)
      enddo
      call dgesv(9,1,T,9,ipvt,V0,9,info)
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



      FUNCTION VAPRIMI(R,N)
CSYSTEM: Ar-OH(A^2Sigma+) 
CLEVEL:RHF/RCCSD(T)/aug-cc-pvqz+33221
C units R=au E=microEh
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(15),XX1(15),XX2(15),XX3(15),
     * XX4(15),XX5(15),XX6(15),XX7(15),XX8(15),
     * XX9(15)
C theta=0
      DATA XX1/
     *    .167041449146508958D+01,
     *   -.333696272770431901D+00,
     *    .830154856794724655D+10,
     *   -.100099884307374420D+11,
     *    .547337025560947418D+10,
     *   -.172259354337886095D+10,
     *    .336835895836324751D+09,
     *   -.416720194531489462D+08,
     *    .317891113664396340D+07,
     *   -.136708435109207756D+06,
     *    .254284630467577381D+04,
     *    .253055840830723319D+05,
     *    .317860854070319608D+08,
     *    .789667970289935350D+09,
     *    .243470736517366314D+10/
C theta=22.5
      DATA XX2/
     *    .598725710484723717D+01,
     *   -.282900604520101737D+02,
     *   -.168128221973281789D+10,
     *    .266564111293966579D+10,
     *   -.181660094054920387D+10,
     *    .691349716012006283D+09,
     *   -.159343008076246321D+09,
     *    .224379012980882861D+08,
     *   -.182999635578120523D+07,
     *    .734514931277571450D+05,
     *   -.830268264724267169D+03,
     *    .253055840830723326D+02,
     *    .350450764312567115D+08,
     *    .652086434384151459D+09,
     *    .105452144064331412D+10/
C theta=45
      DATA XX3/       
     *    .298970691593739835D+01,
     *   -.823520279736410821D+01,
     *   -.120855544625483418D+10,
     *    .210108072995053816D+10,
     *   -.159076820182764220D+10,
     *    .685286523093381763D+09,
     *   -.183737498437959671D+09,
     *    .314194940110719912D+08,
     *   -.334979935395574244D+07,
     *    .203977430750166095D+06,
     *   -.543812823904013294D+04,
     *    .253055840830723326D+02,
     *    .383130205844744146D+08,
     *    .424759986171416938D+09,
     *    .120513928592014343D+09/
C theta=67.5
      DATA XX4/       
     *    .317658357575944228D+01,
     *   -.107924406715431456D+02,
     *   -.160609000341002417D+10,
     *    .251376781010435152D+10,
     *   -.171491662636207628D+10,
     *    .666164413375604272D+09,
     *   -.161189445670164973D+09,
     *    .248831930138184652D+08,
     *   -.239359690994976461D+07,
     *    .131174869717816880D+06,
     *   -.313044605777145989D+04,
     *    .253055840830723326D+02,
     *    .370560181867787763D+08,
     *    .315461556118382454D+09,
     *    .865698576477292478D+08/
C theta=90
      DATA XX5/       
     *    .314613274697740275D+01,
     *   -.110584726982116521D+02,
     *   -.181391789407562447D+10,
     *    .270581758085410452D+10,
     *   -.176040917011285233D+10,
     *    .652578762241834164D+09,
     *   -.150789155851522595D+09,
     *    .222453349933418483D+08,
     *   -.204656896726195654D+07,
     *    .107361410512573129D+06,
     *   -.245560357800092652D+04,
     *    .253055840830723326D+02,
     *    .356043421438017786D+08,
     *    .294205008625618756D+09,
     *    .417965062810200378D+08/
C theta=112.5
      DATA XX6/       
     *    .299718958358611420D+01,
     *   -.944320165898971453D+01,
     *   -.188617825622271967D+10,
     *    .280771271202880383D+10,
     *   -.182417228479559016D+10,
     *    .675962236657286882D+09,
     *   -.156336709836749136D+09,
     *    .231240387484182119D+08,
     *   -.213760602157374192D+07,
     *    .113017113710314938D+06,
     *   -.261616809990241200D+04,
     *    .253055840830723326D+02,
     *    .379672417447080091D+08,
     *    .254540377254802853D+09,
     *    .139569671445584409D+08/
C theta=135
      DATA XX7/       
     *    .287732667967476186D+01,
     *   -.717655640868726952D+01,
     *   -.162158377765151954D+10,
     *    .249739815786805153D+10,
     *   -.168182102830656981D+10,
     *    .648958088866146684D+09,
     *   -.157239490938466758D+09,
     *    .245328686793795228D+08,
     *   -.240914410556727182D+07,
     *    .136375340894662804D+06,
     *   -.340949951121997947D+04,
     *    .253055840830723326D+02,
     *    .426163024501064494D+08,
     *    .174407852476842761D+09,
     *    .112778635971554592D+08/
C theta=157.5
      DATA XX8/       
     *    .276481498766124378D+01,
     *   -.435929502602190855D+01,
     *   -.206828381102243137D+10,
     *    .277329293996062613D+10,
     *   -.161850402204040384D+10,
     *    .569170106573991418D+09,
     *   -.138619178207336068D+09,
     *    .242046429027803130D+08,
     *   -.286503826625642274D+07,
     *    .202296473081753036D+06,
     *   -.636771236082359883D+04,
     *    .253055840830723326D+02,
     *    .466466092335424572D+08,
     *    .933458930334756523D+08,
     *   -.156843347270052527D+05/

C theta=180
      DATA XX9/
     *    .320639359884166808D+01,
     *   -.858021726750815539D+01,
     *   -.922144464833773732D+09,
     *    .186612788984572005D+10,
     *   -.162657175790944386D+10,
     *    .797040672310728073D+09,
     *   -.240061763153670579D+09,
     *    .455385079136549830D+08,
     *   -.532220943708122987D+07,
     *    .351349163891480188D+06,
     *   -.100649181380998543D+05,
     *    .253055840830723326D+02,
     *    .481128280025808588D+08,
     *    .603004105580350831D+08,
     *    .305919855630838796D+08/

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
      IF(N.EQ.8) THEN
      DO I=1,15
      XX(I)=XX8(I)
      ENDDO
      ENDIF
      IF(N.EQ.9) THEN
      DO I=1,15
      XX(I)=XX9(I)
      ENDDO
      ENDIF
      
      TERM1 = DEXP(-XX(1)*R-XX(2))
      TERM2 = XX(3) + XX(4)*R+XX(5)*R**2+
     *        XX(6)*R**3+XX(7)*R**4+XX(8)*R**5
     *        +XX(9)*R**6+XX(10)*R**7+XX(11)*R**8
      TERM3=0.5D0*
     *   (1.D0+ DTANH(1.D0+XX(12)*R))
      TERM4 = XX(13)/(R**6)+XX(14)/(R**7)
     *        +XX(15)/(R**8)
       TOCM=0.219474643545745D0
       VAPRIMI = (TERM1*TERM2-TERM3*TERM4)*TOCM

      IF(N.EQ.1.AND.R.LT.4.D0) THEN
      A=-2.41540493009711D0
      B=16.4277091522484D0
      C=-15.8635420775928D0
      VAPRIMI=DEXP(A*R**2+B*R+C)*TOCM
      ENDIF
      IF(N.EQ.2.AND.R.LT.4.D0) THEN
      A=-2.41540493009711D0
      B=16.4277091522484D0
      C=-15.8635420775928D0
      VAPRIMI=DEXP(A*R**2+B*R+C)*TOCM
      ENDIF
      IF(N.EQ.4.AND.R.LT.4.D0) THEN
      A=-0.156755402048759D0
      B=-0.611757203763026D0
      C=15.8560130076858D0
      VAPRIMI=DEXP(A*R**2+B*R+C)*TOCM
      ENDIF
      IF(N.EQ.5.AND.R.LT.4.D0) THEN
      A=-0.138883473747242D0
      B=-0.722669906228287D0
      C=16.0621883213063D0
      VAPRIMI=DEXP(A*R**2+B*R+C)*TOCM
      ENDIF
      IF(N.EQ.6.AND.R.LT.4.D0) THEN
      A=-0.176698802131906D0
      B=-0.45547722176048D0
      C=15.3897449249026D0
      VAPRIMI=DEXP(A*R**2+B*R+C)*TOCM
      ENDIF
      IF(N.EQ.7.AND.R.LT.4.D0) THEN
      A=0.0181274953590527D0
      B=-2.46880883585686D0
      C=19.7591976984636D0
      VAPRIMI=DEXP(A*R**2+B*R+C)*TOCM
      ENDIF




      RETURN
       END
