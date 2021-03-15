* System:  Cl(-)(X 1Sigma+)+H2,  ab initio CCSD(T) Surfaces
* References:
*  A. A. Buchachenko, T. A. Grinev, J. Klos, E. J. Bieske,
*  M. M. Szczesniak and G. Chalasinski,
*  J. Chem. Phys. 119, 12931, (2003)
*
*  M. H. Alexander, J. Klos and D. Manolopoulos
*  J. Chem. Phys. 128, 084312 (2008)
*
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
      potnam='Klos, Buchachenko Cl(-)-H2 PES for r=re'
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
*  Cl(-)-H2(X) potential of Klos-Buchachenko
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
* adapted for Cl(-)-H2 by Klos 2007
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension xlam1(4),xlam2(4),r0(4),c1(4),c2(4),c3(4),
     :          clr(4),vsum(4),xsum(4),
     :          vap(4), d0(16),aa(25)
      dimension kpvt(7),qraux(7),work(25),rsd(7),Theta(4)

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
       THETA(1)=0.D0
       THETA(2)=30.d0
       THETA(3)=60.d0
       THETA(4)=90.d0
       rh2=1.400D0
      do 100 i=1,4
        vap(i)=POTCLMH2(r,rh2,THETA(i))
100   continue
* solve simultaneous equations for solutions
* first for vsigma
      tol=1.e-10
      call dcopy(16,d0,1,aa,1)
      call dqrank(aa,4,4,4,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,4,4,4,kr,vap,xsum,rsd,kpvt,qraux)
* convert to hartree
      conv=1.d0/219474.6d0
      call dscal(4,conv,xsum,1)
      vv0=xsum(1)
      call dcopy(3,xsum(2),1,vvl,1)
      end


c      PROGRAM POTCHECK
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      dimension thet(91),R(86),rh2(131)
c      DO I=1,91
c      thet(I)=DFLOAT(I-1)*2.D0
c      ENDDO
c      DO I=1,86
c      R(I)=0.0D0+DFLOAT(I-1)*0.1D0
c      ENDDO
c      DO I=1,131
c      rh2(I)=0.2D0+DFLOAT(I-1)*0.01D0
c      ENDDO
c       rd=1.1D0
c      READ(5,*)THETA
c      DO I=1,86
cc       DO J=1,131
c       WRITE(6,*) R(I),POTCLMH2(R(I),1.4003D0,THETA)
cc       ENDDO
c      ENDDO 
c      END

C
      FUNCTION POTCLMH2(R,rh2,THETA)
C FUNCTION FOR (R,r,Theta) DEPENDENCE OF
C Cl(-)-H2  POTENTIAL 
C R,rh2 IN  BOHRS
C OUTPUT IN WAVENUMBERS
C AB INITIO POINTS : CCSD(T) 
C CALCULATIONS BY Dr ALEXEI BUCHACHENKO
C FIT AND EXPANSION IN z^k*P_l(cos(theta)) 
C z=(rh2-r_e)/r_e, r_e=0.7408 A /BY Dr J. KLOS      
C REQUIRES LINKING WITH LAPACK LIBRARY      
C**************************************
C* Dr Jacek Klos                      *
C*  Institute of Theoretical Chemistry*
C*  University of Nijmegen            *
C*  Toernooiveld 1, 6525 ED, Nijmegen *
C*  The Netherlands                   *
C*email: jakl@theochem.kun.nl         *      
C**************************************      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension V0(4)
       pi=dacos(-1.d0)
       CT=DCOS(THETA*pi/180.D0)
       au2Ang=0.529177249D0
       CALL Vl(R*au2Ang,rh2*au2Ang,V0)
       term1=V0(1)*PLGNDR(0,0,CT)
       term2=V0(2)*PLGNDR(2,0,CT)
       term3=V0(3)*PLGNDR(4,0,CT)
       term4=V0(4)*PLGNDR(6,0,CT)
       POTCLMH2=term1+term2+term3+term4
       RETURN
      END
C
C
C
C
      subroutine Vl(R,rsmall,V0)
C EXPANSION IN EVEN ORDINARY LEGENDRE
C POLYNOMIALS UP TO 6th ORDER      
      implicit double precision(a-h, o-z)
      dimension T(4,4),ipvt(4)
      dimension V0(4)
      pi=dacos(-1.d0)
      t1=0.D0
      t2=30.D0
      t3=60.D0
      t4=90.D0
      ct1=DCOS(t1*pi/180.D0)
      ct2=DCOS(t2*pi/180.D0)
      ct3=DCOS(t3*pi/180.D0)
      ct4=DCOS(t4*pi/180.D0)
      T(1,1)=PLGNDR(0,0,ct1)
      T(1,2)=PLGNDR(2,0,ct1)
      T(1,3)=PLGNDR(4,0,ct1)
      T(1,4)=PLGNDR(6,0,ct1)
      T(2,1)=PLGNDR(0,0,ct2)
      T(2,2)=PLGNDR(2,0,ct2)
      T(2,3)=PLGNDR(4,0,ct2)
      T(2,4)=PLGNDR(6,0,ct2)
      T(3,1)=PLGNDR(0,0,ct3)
      T(3,2)=PLGNDR(2,0,ct3)
      T(3,3)=PLGNDR(4,0,ct3)
      T(3,4)=PLGNDR(6,0,ct3)
      T(4,1)=PLGNDR(0,0,ct4)
      T(4,2)=PLGNDR(2,0,ct4)
      T(4,3)=PLGNDR(4,0,ct4)
      T(4,4)=PLGNDR(6,0,ct4)
      V0(1)=POTTH0(R,rsmall)
      V0(2)=POTTH30(R,rsmall)
      V0(3)=POTTH60(R,rsmall)
      V0(4)=POTTH90(R,rsmall)
      call dgesv(4,1,T,4,ipvt,V0,4,info)
C
      return
      end
C
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


C
      FUNCTION POTTH0(R,rh2)
C FUNCTION FOR (R,r) DEPENDENCE OF
C Cl(-)-H2 PES FOR THETA=0 (COLLINEAR)
C R,rh2 IN Angstrems
C OUTPUT IN CM-1
C AB INITIO POINTS : CCSD(T) 
C CALCULATIONS BY Dr ALEXEI BUCHACHENKO
C FIT AND EXPANSION IN z^k BY Dr JACEK KLOS
C REQUIRES LINKING WITH LAPACK LIBRARY      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension V0(4)
       CALL Vlambdath0(R,V0)
       re=0.7408D0
       z=(rh2-re)/re
       term1=V0(1)*z+VTH0(R,3)
       term2=V0(2)*z**2
       term3=V0(3)*z**3
       term4=V0(4)*z**4
       POTTH0=term1+term2+term3+term4
       RETURN
      END
C
C
C
C
      subroutine Vlambdath0(R,V0)
      implicit double precision(a-h, o-z)
      dimension T(4,4),ipvt(4)
      dimension V0(4)
      pi=dacos(-1.d0)
      re=0.7408D0
      r1=0.45D0
      r3=0.60D0
      r4=0.85D0
      r5=1.1D0
      t1=(r1-re)/re
      t3=(r3-re)/re
      t4=(r4-re)/re
      t5=(r5-re)/re
      T(1,1)=t1
      T(1,2)=t1**2
      T(1,3)=t1**3
      T(1,4)=t1**4
      T(2,1)=t3
      T(2,2)=t3**2
      T(2,3)=t3**3
      T(2,4)=t3**4
      T(3,1)=t4
      T(3,2)=t4**2
      T(3,3)=t4**3
      T(3,4)=t4**4
      T(4,1)=t5
      T(4,2)=t5**2
      T(4,3)=t5**3
      T(4,4)=t5**4
      V0(1)=VTH0(R,1)-VTH0(R,3)
      V0(2)=VTH0(R,2)-VTH0(R,3)
      V0(3)=VTH0(R,4)-VTH0(R,3)
      V0(4)=VTH0(R,5)-VTH0(R,3)
      call dgesv(4,1,T,4,ipvt,V0,4,info)
C
      return
      end
C
      FUNCTION VTH0(R,N)
CSYSTEM: Cl(-)-H2  THETA=0 (COLLINEAR)
CLEVEL:RHF/CCSD(T)/ ALEXEI BUCHACHENKO
C units R=A E=microEh
C FIT / JACEK KLOS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(17),XX1(17),XX2(17),XX3(17),
     * XX4(17),XX5(17)
C r=0.45A / RMS=0.01
      DATA XX1/
     *   0.171223729729281948D+01,
     *  -0.309894737301498635D+01,
     *   0.719940705059182801D+05,
     *  -0.860465381830779224D+05,
     *   0.288769244975569200D+05,
     *   0.521684420380041036D+00,
     *  -0.303156952168253110D+04,
     *   0.100128848569549871D+04,
     *  -0.159210422251484545D+03,
     *   0.130431107957872765D+02,
     *  -0.460830904849329204D+00,
     *   0.180741695362990236D-01,
     *   0.253315047974868648D+06,
     *  -0.185598618092404748D+07,
     *   0.605743295502562891D+06,
     *   0.311025370200709067D+07,
     *  -0.231712865302507952D+07/
C r=0.60A / RMS=0.01
      DATA XX2/
     *    0.161426758126935388D+01,
     *   -0.261537621430277545D+01,
     *    0.533551889450118979D+05,
     *   -0.886877992086724407D+05,
     *    0.302294683093666754D+05,
     *    0.430746109129674926D-01,
     *    -0.316008688106317095D+04,
     *    0.101642009973748293D+04,
     *    -0.157776887129108758D+03,
     *    0.125988928562964890D+02,
     *    -0.437862129072214445D+00,
     *    0.207282587315916374D-01,
     *     0.423046439313460025D+06,
     *   -0.333575580599352997D+07,
     *    0.194633532500908384D+07,
     *    0.477292648426634632D+07,
     *   -0.437066775569675956D+07/
C r=re=0.7408A / RMS=0.01
      DATA XX3/
     *    0.151107023389725925D+01,
     *   -0.195715972084233192D+01,
     *    0.932690477335292235D+04,
     *   -0.887917086316965433D+05,
     *     0.319622475829332325D+05,
     *    -0.777208534037589368D-03,
     *    -0.331208502931245221D+04,
     *    0.101658261533088125D+04,
     *    -0.150106625257583488D+03,
     *    0.113363838299233226D+02,
     *   -0.374782181912938051D+00,
     *    0.356578128711631040D-01,
     *   0.583033010507522384D+06,
     *  -0.471118086778729502D+07,
     *   0.346694462059843680D+07,
     *   0.604893642991827801D+07,
     *  -0.641208591852132231D+07/
C r=0.85A / RMS=0.01
      DATA XX4/
     *    0.132653945555429376D+01,
     *  -0.139653879242908663D+01,
     *  -0.179665548795161565D+06,
     *  -0.121833722895891460D+05,
     *   0.282628895079033100D+05,
     *  -0.794296071252643742D+04,
     *  -0.139109425385837787D+00,
     *   0.325089287629446801D+03,
     *  -0.646663974874380187D+02,
     *   0.534625836837256685D+01,
     *  -0.182242970559121809D+00,
     *   0.179868192506720637D+01,
     *   0.856987102935241535D+06,
     *  -0.913372313112292252D+07,
     *   0.156405215201796703D+08,
     *  -0.679474435838123132D+07,
     *  -0.180908665530709527D+07/
C r=1.1A / RMS=0.03
      DATA XX5/
     *    0.137792336618727185D+01,
     *   -0.185150817325020167D+01,
     *  0.854696797173463363D+00,
     * -0.874915556705465278D+05,
     *  0.261037052246335479D+05,
     *  0.387147611321633030D+04,
     * -0.428512445800887963D+04,
     *  0.110473599213357397D+04,
     * -0.144713653697217808D+03,
     *  0.980276216526209332D+01,
     * -0.288514675120921660D+00,
     *  0.421425441067243417D-01,
     *  0.115168149031761545D+07,
     * -0.114901657376439199D+08,
     *  0.244856121901887394D+08,
     * -0.177392548471592031D+08,
     * 0.138181825925986725D+07/
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
C
C
      TERM1 = DEXP(-XX(1)*R-XX(2))
      TERM2 = XX(3) + XX(4)*R+XX(5)*R**2+
     *        XX(6)*R**3+XX(7)*R**4+XX(8)*R**5
     *        +XX(9)*R**6+XX(10)*R**7+XX(11)*R**8
      TERM3=0.5D0*
     *   (1.D0+ DTANH(1.D0+XX(12)*R))
      TERM4 = XX(13)/(R**4)+XX(14)/(R**5)
     *        +XX(15)/(R**6)+XX(16)/(R**7)
     *        +XX(17)/(R**8)
C
       VTH0 = (TERM1*TERM2-TERM3*TERM4)
       RETURN
       END


C
      FUNCTION POTTH30(R,rh2)
C FUNCTION FOR (R,r) DEPENDENCE OF
C Cl(-)-H2 PES FOR THETA=30
C R,rh2 IN Angstrems
C OUTPUT IN CM-1
C AB INITIO POINTS : CCSD(T) 
C CALCULATIONS BY Dr ALEXEI BUCHACHENKO
C FIT AND EXPANSION in z^k BY Dr JACEK KLOS
C REQUIRES LINKING WITH LAPACK LIBRARY      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension V0(4)
       CALL Vlambdath30(R,V0)
       re=0.7408D0
       z=(rh2-re)/re
       term1=V0(1)*z+VTH30(R,3)
       term2=V0(2)*z**2
       term3=V0(3)*z**3
       term4=V0(4)*z**4
       POTTH30=term1+term2+term3+term4
       RETURN
      END
C
C
C
C
      subroutine Vlambdath30(R,V0)
      implicit double precision(a-h, o-z)
      dimension T(4,4),ipvt(4)
      dimension V0(4)
      pi=dacos(-1.d0)
      re=0.7408D0
      r1=0.45D0
      r3=0.60D0
      r4=0.85D0
      r5=1.1D0
      t1=(r1-re)/re
      t3=(r3-re)/re
      t4=(r4-re)/re
      t5=(r5-re)/re
      T(1,1)=t1
      T(1,2)=t1**2
      T(1,3)=t1**3
      T(1,4)=t1**4
      T(2,1)=t3
      T(2,2)=t3**2
      T(2,3)=t3**3
      T(2,4)=t3**4
      T(3,1)=t4
      T(3,2)=t4**2
      T(3,3)=t4**3
      T(3,4)=t4**4
      T(4,1)=t5
      T(4,2)=t5**2
      T(4,3)=t5**3
      T(4,4)=t5**4
      V0(1)=VTH30(R,1)-VTH30(R,3)
      V0(2)=VTH30(R,2)-VTH30(R,3)
      V0(3)=VTH30(R,4)-VTH30(R,3)
      V0(4)=VTH30(R,5)-VTH30(R,3)
      call dgesv(4,1,T,4,ipvt,V0,4,info)
C
      return
      end
C
      FUNCTION VTH30(R,N)
CSYSTEM: Cl(-)-H2  THETA=30
CLEVEL:RHF/CCSD(T)/ ALEXEI BUCHACHENKO
C units R=A E=microEh
C FIT / JACEK KLOS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(17),XX1(17),XX2(17),XX3(17),
     * XX4(17),XX5(17)
C r=0.45A / RMS=0.01
      DATA XX1/
     *  0.139032476927790216D+01,
     * -0.134627575579868441D+01,
     * -0.123386625039271137D+06,
     * -0.154787023018665059D+05,
     *  0.313394321999330423D+05,
     * -0.116139278215055165D+05,
     *  0.158922008485397964D+04,
     *  0.866120588174998379D-04,
     * -0.282571668801545250D+02,
     *  0.323304222649032624D+01,
     * -0.130475869453770710D+00,
     * -0.289071131452723999D-01,
     *  0.311724253479954496D+06,
     * -0.324775472054268094D+07,
     *  0.145066498652969068D+07,
     *  0.474516432812912576D+07,
     * -0.369892559975029808D+07/

C r=0.60A / RMS=0.01
      DATA XX2/
     *    0.131797415485845337D+01,
     *   -0.145661888803034367D+01,
     *   -0.127294200792077900D+06,
     *   -0.229407480365412403D+04,
     *    0.285386311017458356D+05,
     *   -0.130714693941969072D+05,
     *    0.270190158315010649D+04,
     *   -0.299900379617490898D+03,
     *    0.148198614564063718D+02,
     *    0.904053113819259802D-04,
     *   -0.244826191343795393D-01,
     *   -0.161660253032696495D-01,
     *    0.406288644407385204D+06,
     *   -0.396571931976757618D+07,
     *    0.182430928147088806D+07,
     *    0.654669079155520909D+07,
     *   -0.547452502410516702D+07/
C r=re=0.7408A / RMS=0.01
      DATA XX3/
     *    0.149227880471129848D+01,
     *   -0.201228523630604483D+01,
     *    0.307380169095503675D+01,
     *   -0.810561661510545819D+05,
     *    0.326724586433967852D+05,
     *   -0.612927095810405604D+03,
     *   -0.357240411711617480D+04,
     *    0.121009652313099264D+04,
     *   -0.192568521824217527D+03,
     *    0.154952959906180698D+02,
     *   -0.527883999206096188D+00,
     *   -0.655539927859111063D-03,
     *    0.500978340160817141D+06,
     *   -0.458252985648713913D+07,
     *    0.344024782073218655D+07,
     *    0.578442800270957407D+07,
     *   -0.598443268935817946D+07/
C r=0.85A / RMS=0.01
      DATA XX4/
     *    0.132531143082507796D+01,
     *   -0.129514685491954906D+01,
     *   -0.144118384264626337D+06,
     *   -0.167839706868529793D+05,
     *    0.313338442595869456D+05,
     *   -0.108469788017407809D+05,
     *    0.138497916232003899D+04,
     *    0.207836946673194342D-01,
     *   -0.220215370693219796D+02,
     *    0.239269633706107587D+01,
     *   -0.929225734252175500D-01,
     *    0.140428328217546104D-01,
     *    0.566354667620536522D+06,
     *   -0.519432982407136634D+07,
     *    0.358780918144316226D+07,
     *    0.770003556793260761D+07,
     *   -0.783157536088911723D+07/
C r=1.1A / RMS=0.03
      DATA XX5/
     *    0.147976102671807674D+01,
     *   -0.171423900407105001D+01,
     *    0.100238172041578188D+02,
     *   -0.909935294321298279D+05,
     *    0.307953966409183158D+05,
     *    0.307843699230599725D+03,
     *   -0.311905377714521910D+04,
     *    0.910134889076069840D+03,
     *   -0.130267299735922762D+03,
     *    0.962395290789412883D+01,
     *   -0.317532331359072162D+00,
     *    0.396394251606365655D-01,
     *    0.671437723310046946D+06,
     *   -0.525155832046271767D+07,
     *    0.457745156524579879D+07,
     *    0.715297652135270089D+07,
     *   -0.889885703211604618D+07/
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
C
C
      TERM1 = DEXP(-XX(1)*R-XX(2))
      TERM2 = XX(3) + XX(4)*R+XX(5)*R**2+
     *        XX(6)*R**3+XX(7)*R**4+XX(8)*R**5
     *        +XX(9)*R**6+XX(10)*R**7+XX(11)*R**8
      TERM3=0.5D0*
     *   (1.D0+ DTANH(1.D0+XX(12)*R))
      TERM4 = XX(13)/(R**4)+XX(14)/(R**5)
     *        +XX(15)/(R**6)+XX(16)/(R**7)
     *        +XX(17)/(R**8)
C
       VTH30 = (TERM1*TERM2-TERM3*TERM4)
       RETURN
       END


C
      FUNCTION POTTH60(R,rh2)
C FUNCTION FOR (R,r) DEPENDENCE OF
C Cl(-)-H2 PES FOR THETA=60
C R,rh2 IN Angstrems
C OUTPUT IN CM-1
C AB INITIO POINTS : CCSD(T) 
C CALCULATIONS BY Dr ALEXEI BUCHACHENKO
C FIT AND EXPANSION BY JACEK KLOS
C REQUIRES LINKING WITH LAPACK LIBRARY      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension V0(4)
       CALL Vlambdath60(R,V0)
       re=0.7408D0
       z=(rh2-re)/re
       term1=V0(1)*z+VTH60(R,3)
       term2=V0(2)*z**2
       term3=V0(3)*z**3
       term4=V0(4)*z**4
       POTTH60=term1+term2+term3+term4
       RETURN
      END
C
C
C
C
      subroutine Vlambdath60(R,V0)
      implicit double precision(a-h, o-z)
      dimension T(4,4),ipvt(4)
      dimension V0(4)
      pi=dacos(-1.d0)
      re=0.7408D0
      r1=0.45D0
      r3=0.60D0
      r4=0.85D0
      r5=1.1D0
      t1=(r1-re)/re
      t3=(r3-re)/re
      t4=(r4-re)/re
      t5=(r5-re)/re
      T(1,1)=t1
      T(1,2)=t1**2
      T(1,3)=t1**3
      T(1,4)=t1**4
      T(2,1)=t3
      T(2,2)=t3**2
      T(2,3)=t3**3
      T(2,4)=t3**4
      T(3,1)=t4
      T(3,2)=t4**2
      T(3,3)=t4**3
      T(3,4)=t4**4
      T(4,1)=t5
      T(4,2)=t5**2
      T(4,3)=t5**3
      T(4,4)=t5**4
      V0(1)=VTH60(R,1)-VTH60(R,3)
      V0(2)=VTH60(R,2)-VTH60(R,3)
      V0(3)=VTH60(R,4)-VTH60(R,3)
      V0(4)=VTH60(R,5)-VTH60(R,3)
      call dgesv(4,1,T,4,ipvt,V0,4,info)
C
      return
      end
C
      FUNCTION VTH60(R,N)
CSYSTEM: Cl(-)-H2  THETA=60
CLEVEL:RHF/CCSD(T)/ ALEXEI BUCHACHENKO
C units R=A E=microEh
C FIT / JACEK KLOS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(17),XX1(17),XX2(17),XX3(17),
     * XX4(17),XX5(17)
C r=0.45A / RMS=0.003
      DATA XX1/
     *    0.140695045880379621D+01,
     *   -0.137927006399580510D+01,
     *   -0.150198792127149820D+06,
     *    0.101001716334281500D+05,
     *    0.240462242551462878D+05,
     *   -0.125863414710397446D+05,
     *    0.273146973182596776D+04,
     *   -0.316638919194688242D+03,
     *    0.161037353220220396D+02,
     *   -0.315699122577735823D-02,
     *   -0.299762055119252233D-01,
     *   -0.103491239402443716D+00,
     *    0.250257834237943433D+06,
     *   -0.323201174300837982D+07,
     *    0.148081829264557967D+07,
     *    0.484406038100257795D+07,
     *   -0.377210811041178694D+07/
C r=0.60A / RMS=0.004
      DATA XX2/
     *    0.704203450770159112D+00,
     *    0.474545211620003204D+01,
     *   -0.209853878187665077D+04,
     *   -0.868035669773566537D+06,
     *    0.514917353458871250D+06,
     *   -0.138591782756275148D+06,
     *    0.210974740089464176D+05,
     *   -0.193953046087472512D+04,
     *    0.106617805036503071D+03,
     *   -0.322961213107566891D+01,
     *    0.412366755042961400D-01,
     *   -0.715075495234876790D+00,
     *    0.133017474983237730D+07,
     *   -0.102678564006306045D+08,
     *    0.196114868319640718D+08,
     *   -0.151459141726626735D+08,
     *    0.436975600044166576D+07/
C r=re=0.7408A / RMS=0.006
      DATA XX3/
     *    0.696962058572321963D+00,
     *    0.419914199487625339D+01,
     *   -0.191896978383664396D+04,
     *   -0.666223697928342503D+06,
     *    0.411880843047917704D+06,
     *   -0.115267041869000896D+06,
     *    0.182171944901430907D+05,
     *   -0.173200769687580737D+04,
     *    0.979956674371649257D+02,
     *   -0.303913065021967599D+01,
     *    0.396667542394386063D-01,
     *   -0.697555802992127494D+00,
     *    0.144434761614709883D+07,
     *   -0.107763702991164792D+08,
     *    0.203297450254042670D+08,
     *   -0.148565742863274831D+08,
     *    0.365067773455836857D+07/
C r=0.85A / RMS=0.004
      DATA XX4/
     *    0.658106669567382285D+00,
     *    0.399986826645696381D+01,
     *    0.201647798640282638D+07,
     *   -0.198711127476505004D+07,
     *    0.786362025550031452D+06,
     *   -0.172514869919821504D+06,
     *    0.230824264877182723D+05,
     *   -0.193185863213234961D+04,
     *    0.986336704770926787D+02,
     *   -0.280614131046438686D+01,
     *    0.340056679242773247D-01,
     *   -0.581233197314076522D+00,
     *    0.196798333109608386D+07,
     *   -0.942397431188258156D+07,
     *    0.109061379627501965D+08,
     *   -0.211613599986839972D+04,
     *   -0.400794816880642949D+07/
C r=1.1A / RMS=0.008
      DATA XX5/
     *    0.629568672496801263D+00,
     *    0.432110705501143766D+01,
     *   -0.168488396227247734D+07,
     *    0.714207928964297171D+06,
     *   -0.110383136913591123D+06,
     *   -0.637040740709343911D+01,
     *    0.228844506436141046D+04,
     *   -0.337630160512065572D+03,
     *    0.232251383592715506D+02,
     *   -0.805807400108535132D+00,
     *    0.113655322290210446D-01,
     *   -0.637705233205278721D+00,
     *    0.163121816770550655D+07,
     *   -0.112210169098399282D+08,
     *    0.210769615488709584D+08,
     *   -0.143301704185998514D+08,
     *    0.238392591009412147D+07/
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
C
C
      TERM1 = DEXP(-XX(1)*R-XX(2))
      TERM2 = XX(3) + XX(4)*R+XX(5)*R**2+
     *        XX(6)*R**3+XX(7)*R**4+XX(8)*R**5
     *        +XX(9)*R**6+XX(10)*R**7+XX(11)*R**8
      TERM3=0.5D0*
     *   (1.D0+ DTANH(1.D0+XX(12)*R))
      TERM4 = XX(13)/(R**4)+XX(14)/(R**5)
     *        +XX(15)/(R**6)+XX(16)/(R**7)
     *        +XX(17)/(R**8)
C
       VTH60 = (TERM1*TERM2-TERM3*TERM4)
       RETURN
       END


C
      FUNCTION POTTH90(R,rh2)
C FUNCTION FOR (R,r) DEPENDENCE OF
C Cl(-)-H2 PES FOR THETA=90 (T-shaped)
C R,rh2 IN Angstrems
C OUTPUT IN CM-1
C AB INITIO POINTS : CCSD(T) 
C CALCULATIONS BY Dr ALEXEI BUCHACHENKO
C FIT AND EXPANSION BY Dr JACEK KLOS
C REQUIRES LINKING WITH LAPACK LIBRARY      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension V0(4)
       CALL Vlambdath90(R,V0)
       re=0.7408D0
       z=(rh2-re)/re
       term1=V0(1)*z+VTH90(R,3)
       term2=V0(2)*z**2
       term3=V0(3)*z**3
       term4=V0(4)*z**4
       POTTH90=term1+term2+term3+term4
       RETURN
      END
C
C
C
C
      subroutine Vlambdath90(R,V0)
      implicit double precision(a-h, o-z)
      dimension T(4,4),ipvt(4)
      dimension V0(4)
      pi=dacos(-1.d0)
      re=0.7408D0
      r1=0.45D0
      r3=0.60D0
      r4=0.85D0
      r5=1.1D0
      t1=(r1-re)/re
      t3=(r3-re)/re
      t4=(r4-re)/re
      t5=(r5-re)/re
      T(1,1)=t1
      T(1,2)=t1**2
      T(1,3)=t1**3
      T(1,4)=t1**4
      T(2,1)=t3
      T(2,2)=t3**2
      T(2,3)=t3**3
      T(2,4)=t3**4
      T(3,1)=t4
      T(3,2)=t4**2
      T(3,3)=t4**3
      T(3,4)=t4**4
      T(4,1)=t5
      T(4,2)=t5**2
      T(4,3)=t5**3
      T(4,4)=t5**4
      V0(1)=VTH90(R,1)-VTH90(R,3)
      V0(2)=VTH90(R,2)-VTH90(R,3)
      V0(3)=VTH90(R,4)-VTH90(R,3)
      V0(4)=VTH90(R,5)-VTH90(R,3)
      call dgesv(4,1,T,4,ipvt,V0,4,info)
C
      return
      end
C
      FUNCTION VTH90(R,N)
CSYSTEM: Cl(-)-H2  THETA=90
CLEVEL:RHF/CCSD(T)/ ALEXEI BUCHACHENKO
C units R=A E=microEh
C FIT / JACEK KLOS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(17),XX1(17),XX2(17),XX3(17),
     * XX4(17),XX5(17)
C r=0.45A / RMS=0.003
      DATA XX1/
     *    0.613835959460779401D+00,
     *    0.309454640133259051D+01,
     *    0.135321933427209729D+02,
     *   -0.780367590852468420D+05,
     *    0.473799849239383184D+05,
     *   -0.127570737607083865D+05,
     *    0.193649592834947293D+04,
     *   -0.177093588310402367D+03,
     *    0.971509750436553610D+01,
     *   -0.295377219709907490D+00,
     *    0.386633255732044814D-02,
     *   -0.663267700255923076D+00,
     *    0.629610277442689985D+06,
     *   -0.374425968353941059D+07,
     *    0.161285273268963676D+07,
     *    0.529471940378062986D+07,
     *   -0.411741315643498441D+07/
C r=0.60A / RMS=0.003
      DATA XX2/
     *    0.587980982644684569D+00,
     *    0.315649303149274196D+01,
     *   -0.313952735611494980D+02,
     *   -0.784057317878382310D+05,
     *    0.475866833870262053D+05,
     *   -0.127083503436684514D+05,
     *    0.191183657006472549D+04,
     *   -0.173077324435009956D+03,
     *    0.939957452608550348D+01,
     *   -0.282953922830049276D+00,
     *    0.367758762211109957D-02,
     *   -0.632510991542724499D+00,
     *    0.844003780690020183D+06,
     *   -0.460999167337300815D+07,
     *    0.275957993055702047D+07,
     *    0.520336385614155047D+07,
     *   -0.455988417737672105D+07/
C r=re=0.7408A / RMS=0.008
      DATA XX3/
     *    0.512408245963489106D+00,
     *    0.416928473774896169D+01,
     *   -0.181079933390793798D+06,
     *    0.210591803309326142D+01,
     *    0.355208178830138713D+05,
     *   -0.121195527986485122D+05,
     *    0.195165609065968943D+04,
     *   -0.177336465598537302D+03,
     *    0.934194896884068271D+01,
     *   -0.266360391875546276D+00,
     *    0.322098615291943570D-02,
     *   -0.602506588414989386D+00,
     *    0.985358869511889410D+06,
     *   -0.515159241967203282D+07,
     *    0.349772546958416887D+07,
     *    0.520185178583765496D+07,
     *   -0.492940015614692122D+07/
C r=0.85A / RMS=0.006
      DATA XX4/
     *    0.569910972787980530D+00,
     *    0.308213307590109764D+01,
     *    0.132925880001034602D+01,
     *   -0.764955602578773105D+05,
     *    0.471241314729965161D+05,
     *   -0.126270954938210871D+05,
     *    0.190311723728890615D+04,
     *   -0.172319699578743666D+03,
     *    0.935928680123202383D+01,
     *   -0.281713737004089360D+00,
     *    0.366842146106331236D-02,
     *   -0.615117541262720890D+00,
     *    0.132223785853271931D+07,
     *   -0.791803779652955011D+07,
     *    0.111949722719332408D+08,
     *   -0.371366761559249228D+07,
     *   -0.115805240223872173D+07/
C r=1.1A / RMS=0.012
      DATA XX5/
     *     0.509061396510676722D+00,
     *     0.381876517132442439D+01,
     *    -0.285127834938180051D+06,
     *     0.141911773459697346D+06,
     *    -0.280332682934475961D+05,
     *     0.249710287106528949D+04,
     *     0.919026173264051449D-01,
     *    -0.201892937096076750D+02,
     *     0.185649210809466725D+01,
     *    -0.730610718028838341D-01,
     *     0.115396694229618051D-02,
     *    -0.559595761015455606D+00,
     *     0.136274837722120108D+07,
     *    -0.774679103899470996D+07,
     *     0.106861011875193939D+08,
     *    -0.270441513625251502D+07,
     *    -0.194465682715598308D+07/
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
C
C
      TERM1 = DEXP(-XX(1)*R-XX(2))
      TERM2 = XX(3) + XX(4)*R+XX(5)*R**2+
     *        XX(6)*R**3+XX(7)*R**4+XX(8)*R**5
     *        +XX(9)*R**6+XX(10)*R**7+XX(11)*R**8
      TERM3=0.5D0*
     *   (1.D0+ DTANH(1.D0+XX(12)*R))
      TERM4 = XX(13)/(R**4)+XX(14)/(R**5)
     *        +XX(15)/(R**6)+XX(16)/(R**7)
     *        +XX(17)/(R**8)
C
       VTH90 = (TERM1*TERM2-TERM3*TERM4)
       RETURN
       END
