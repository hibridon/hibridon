*Kr-NO(A^2Sigma 3ssigma) PES RCCSD(T)+DKROLL
*Reference:
* J. Klos, M. H. Alexander, R. Hernandez-Lamoneda and T. G. Wright,
* J. Chem. Phys. 129, 244303 (2008)

      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(10)
      include "common/parpot"
      potnam='Klos et al Kr-NO(A^2Sigma+) DKROLL RCCSDT PES'
      print *, potnam
1      print *, ' r (bohr)'
      rshift=0.5
      xfact=0.8
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vsum',/,7(1pe16.8))
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
      potnam='Klos et al Kr-NO(A^2Sigma) DKROLL RCCSDT PES'
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
*revised for Kr-NO(A^2Sigma+) 2007 J.Klos
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension xlam1(7),xlam2(7),r0(7),
     :          vsum(11),xsum(11),
     :          d0(121),aa(121),thta(11),cthta(11)
      dimension kpvt(11),qraux(11),work(121),rsd(11),re(14)

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
C System: Kr-NO(A2Sigma 3ssigma )
C Method:RCCSD(T)
C Basis:aug-cc-pvtz+3s3p2d2f1g1h+dkroll
C Including Douglas-Kroll relativistic effect
C CP and BSSE corrected
C Variables:
C R in Bohr
C Theta n Degrees
C Output in cm-1
C Theta=0 for Kr---NO geom.
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
       VPES=s
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
CSYSTEM: Kr-NO(Sigma) 2ssigma Rydberg state Theta=0 for Kr-NO
CLEVEL:RHF/RCCSD(T)/aug-cc-pvtz+332211
C units R=au E=microEh(parameters)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(15),XX1(15),XX2(15),XX3(15),
     * XX4(15),XX5(15),XX6(15),XX7(15),
     * XX8(15),XX9(15),XX10(15),XX11(15)
C theta=0
      DATA XX1/
     *    .164712471638325830D+01,
     *   -.419332624316846481D+01,
     *    .256277414772642422D+10,
     *   -.147797650576702261D+10,
     *    .286542638737985730D+09,
     *   -.768309056340548396D+07,
     *   -.496847779646645486D+07,
     *    .767753835295920027D+06,
     *   -.457732167083797540D+05,
     *    .948319418651164028D+03,
     *    .636484379096013786D+01,
     *    .904264556193752966D+01,
     *    .617425496217356920D+09,
     *   -.606861347362439442D+10,
     *    .118904194855057205D+12/


C theta=30
      DATA XX2/
     *    .170942297067917326D+01,
     *   -.698382510008604823D-01,
     *    .662942133397109866D+09,
     *   -.279407765476401255D+06,
     *    .197669580927786529D+09,
     *   -.443529611542983800D+08,
     *   -.357557811594507610D+07,
     *    .172835478339616419D+07,
     *   -.707010082358985092D+05,
     *   -.140168629572804530D+05,
     *    .121038835318113570D+04,
     *    .904264556193752966D+01,
     *    .622069551141422153D+09,
     *   -.557276006606412220D+10,
     *    .115852804825448868D+12/

C theta=60
      DATA XX3/       
     *    .166579192815633426D+01,
     *   -.896310813086594660D+00,
     *    .414543762596067047D+10,
     *   -.173168098482884479D+10,
     *    .231485905532741845D+09,
     *    .142951126658383720D+08,
     *   -.533162300839414261D+07,
     *    .239029691393645329D+02,
     *    .109218211210223584D+06,
     *   -.121379124266172967D+05,
     *    .505041964306414570D+03,
     *    .904264556193752966D+01,
     *    .734425227336765766D+09,
     *   -.102558051766473446D+11,
     *    .184472942245054260D+12/

C theta=90
      DATA XX4/       
     *    .175420853312275726D+01,
     *   -.193505575176768851D+01,
     *    .271787121939807749D+10,
     *   -.968658906661303759D+09,
     *   -.350051443032424822D+08,
     *    .787753903276718110D+08,
     *   -.142845985378605910D+08,
     *    .331772931149782205D+06,
     *    .184523607288288156D+06,
     *   -.224574764382473695D+05,
     *    .907372295070675250D+03,
     *    .904264556193752966D+01,
     *    .771230547979266763D+09,
     *   -.119438211273232479D+11,
     *    .217928612391660919D+12/

C theta=120
      DATA XX5/       
     *    .172617409588513748D+01,
     *   -.854893473648289659D-02,
     *    .576116736428887272D+10,
     *   -.149515751916462445D+10,
     *    .188414483096892625D+09,
     *   -.104054297868861333D+08,
     *   -.215827082367856381D+07,
     *    .106098864504383714D+07,
     *   -.932981031308511738D+05,
     *   -.358370504323249042D+04,
     *    .884746142531126679D+03,
     *    .904264556193752966D+01,
     *    .719046253499636769D+09,
     *   -.892541888913840485D+10,
     *    .177896043404205261D+12/

C theta=150
      DATA XX6/
     *    .179172522260526201D+01,
     *   -.235050969511213292D+01,
     *    .108887193060283422D+10,
     *   -.808452941081363440D+09,
     *    .247651287831033379D+09,
     *   -.167560142542067375D+08,
     *   -.537727077757506259D+07,
     *    .960594050883064861D+06,
     *   -.639204020984089460D+00,
     *   -.940173027963206550D+04,
     *    .577993712733358279D+03,
     *    .904264556193752966D+01,
     *    .653118598607211828D+09,
     *   -.608682769918468571D+10,
     *    .144610083943216766D+12/


C theta=180
      DATA XX7/
     *    .189662280475808864D+01,
     *   -.467200035098823907D+01,
     *    .224045498887993670D+10,
     *   -.147138715952513194D+10,
     *    .290708726098975301D+09,
     *    .118239571499489695D+08,
     *   -.115236586789156813D+08,
     *    .131719622772722435D+07,
     *    .179661563044265472D+04,
     *   -.895808166645952770D+04,
     *    .447551041549007664D+03,
     *    .904264556193752966D+01,
     *    .616586058247982860D+09,
     *   -.468030299956171417D+10,
     *    .129912212745167404D+12/


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

