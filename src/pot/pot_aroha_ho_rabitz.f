*Ar-OH(A^2Sigma+) PES of Ho and Rabitz
*References: 
* Tak‐San Ho, Herschel Rabitz, Seung E. Choi, 
* and Marsha I. Lester J. Chem. Phys. 102, 2282 (1995)
* also  same authors   J. Chem. Phys. 104, 1187 (1996)
      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(6)
      include "common/parpot"
      potnam='Ar-OH(A^2Sigma+) PES of Ho and Rabitz'
      print *, potnam
1      print *, ' r (bohr)'
      rshift=0.5
      xfact=0.8
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vsum',/,7(1pe16.8))
      goto 1
99    rmin=4
      dr=0.25
      do i=1,70
         r=rmin+(i-1)*dr
         call pot(vv0,r)
         write (9,101) r, vv0, (vvl(ik),ik=1,6)
101      format(f9.5,7(1pg12.4))
      enddo
      end

      include "common/syusr"
      include "common/bausr"
      include "common/ground"
      subroutine loapot(iunit,filnam)
* ------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      potnam='Ar-OH(A^2Sigma+) PES of Ho and Rabitz'
      lammin(1)=1
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
*  Ar-OH(A) potential of Ho Rabitz Choi Lester
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*  vv0        contains isotropic term (P0 term)
*  variable in common block /covvl/
*    vvl:     vector of length 6 to store r-dependence of each term
*             in potential expansion
*    vvl(1-6) expansion coefficients in Pl0 (l=1:6)
*  variables in common block /coexpan/
*    rang:    spline knots (angstroms)
*    vln (n=0:6):  expansion coefficients (hartree) from tak-sun ho

* author:  tak-sun ho
* revised by mha:  27-apr-2008
* ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (nr=50)
      dimension rang(nr),vl0(nr),vl1(nr),vl2(nr),vl3(nr),vl4(nr),
     : vl5(nr),vl6(nr),vvl(6)
* spline values for legendre term (l=0:6)
      dimension pl0(nr),pl1(nr),pl2(nr),pl3(nr),pl4(nr),pl5(nr),pl6(nr),
     :  dpl0(nr),dpl1(nr),dpl2(nr),dpl3(nr),dpl4(nr),dpl5(nr),dpl6(nr),
     :          xx(nr),yy(nr),dd(nr)
      common /coexpan/ rang,vl0,vl1,vl2,vl3,vl4,vl5,vl6
      common /covvl/ vvl
      data interp /1/
* determine spline coefficients for first pass
      kr=nr
      if (interp .eq. 1) then
         call dcopy(kr,rang,1,xx,1)
         call dcopy(kr,vl0,1,yy,1)
         dfinal=(yy(nr)-yy(nr-1))/(xx(nr)-xx(nr-1))
         dinitial=(yy(2)-yy(1))/(xx(2)-xx(1))
         call splinex(xx,yy,kr,dinitial,dfinal,dd)
         call dcopy(kr,yy,1,pl0,1)
         call dcopy(kr,dd,1,dpl0,1)

         call dcopy(kr,vl1,1,yy,1)
         dfinal=(yy(nr)-yy(nr-1))/(xx(nr)-xx(nr-1))
         dinitial=(yy(2)-yy(1))/(xx(2)-xx(1))
         call splinex(xx,yy,kr,dinitial,dfinal,dd)
         call dcopy(kr,dd,1,dpl1,1)
         dfinal=(yy(nr)-yy(nr-1))/(xx(nr)-xx(nr-1))
         dinitial=(yy(2)-yy(1))/(xx(2)-xx(1))
         call splinex(xx,yy,kr,dinitial,dfinal,dd)
         call dcopy(kr,yy,1,pl1,1)
         call dcopy(kr,dd,1,dpl1,1)

         call dcopy(kr,vl2,1,yy,1)
         dfinal=(yy(nr)-yy(nr-1))/(xx(nr)-xx(nr-1))
         dinitial=(yy(2)-yy(1))/(xx(2)-xx(1))
         call splinex(xx,yy,kr,dinitial,dfinal,dd)
         call dcopy(kr,yy,1,pl2,1)
         call dcopy(kr,dd,1,dpl2,1)

         call dcopy(kr,vl3,1,yy,1)
         dfinal=(yy(nr)-yy(nr-1))/(xx(nr)-xx(nr-1))
         dinitial=(yy(2)-yy(1))/(xx(2)-xx(1))
         call splinex(xx,yy,kr,dinitial,dfinal,dd)
         call dcopy(kr,yy,1,pl3,1)
         call dcopy(kr,dd,1,dpl3,1)

         call dcopy(kr,vl4,1,yy,1)
         dfinal=(yy(nr)-yy(nr-1))/(xx(nr)-xx(nr-1))
         dinitial=(yy(2)-yy(1))/(xx(2)-xx(1))
         call splinex(xx,yy,kr,dinitial,dfinal,dd)
         call dcopy(kr,yy,1,pl4,1)
         call dcopy(kr,dd,1,dpl4,1)

         call dcopy(kr,vl5,1,yy,1)
         dfinal=(yy(nr)-yy(nr-1))/(xx(nr)-xx(nr-1))
         dinitial=(yy(2)-yy(1))/(xx(2)-xx(1))
         call splinex(xx,yy,kr,dinitial,dfinal,dd)
         call dcopy(kr,yy,1,pl5,1)
         call dcopy(kr,dd,1,dpl5,1)

         call dcopy(kr,vl6,1,yy,1)
         dfinal=(yy(nr)-yy(nr-1))/(xx(nr)-xx(nr-1))
         dinitial=(yy(2)-yy(1))/(xx(2)-xx(1))
         call splinex(xx,yy,kr,dinitial,dfinal,dd)
         call dcopy(kr,yy,1,pl6,1)
         call dcopy(kr,dd,1,dpl6,1)
         interp=0
      endif
* convert distance to angstroms
      rr=r*0.52917715d0
* determine expansion coefficients
      if (rr .lt. rang(1) ) then
         vr0=vl0(1)
         vr1=vl1(1)
         vr2=vl2(1)
         vr3=vl3(1)
         vr4=vl4(1)
         vr5=vl5(1)
         vr6=vl6(1)
      else if (rr .gt. rang(kr)) then
* r-6 extrapolation for lambda=0 and lambda=2
         rm6=1d0/rr**6
         rm7=rm6/rr
         vr0=-1.2046*rm6
         vr2=-0.1772*rm6
* r-7 extrapolation for lambda=1 and lambda=3
         vr1=-1.8935d0*rm7
         vr3=-0.9745d0*rm7
* set others equal to zero
         vr4=0d0
         vr5=0d0
         vr6=0d0
      else
* here for spline interpolation
         call splint0(rang,pl0,dpl0,kr,rr,vr0)
         call splint0(rang,pl1,dpl1,kr,rr,vr1)
         call splint0(rang,pl2,dpl2,kr,rr,vr2)
         call splint0(rang,pl3,dpl3,kr,rr,vr3)
         call splint0(rang,pl4,dpl4,kr,rr,vr4)
         call splint0(rang,pl5,dpl5,kr,rr,vr5)
         call splint0(rang,pl6,dpl6,kr,rr,vr6)
      endif
      vvl(1)=vr1
      vvl(2)=vr2
      vvl(3)=vr3
      vvl(4)=vr4
      vvl(5)=vr5
      vvl(6)=vr6
      vv0=vr0
      return
      end
      block data ho_rabitz
      implicit real*8(a-h,o-z)
      parameter (nr=50)
      dimension rang(nr),vl0(nr),vl1(nr),vl2(nr),vl3(nr),vl4(nr),
     :                  vl5(nr),vl6(nr)
      common /coexpan/ rang,vl0,vl1,vl2,vl3,vl4,vl5,vl6
*  jacobi vector (Angstrom) and legendre expansion coefficients vllam (lam=0-6)
*  in hartree
      data rang/
     : 1.9000000d0, 1.9360701d0, 1.9728833d0, 2.0116627d0, 2.0528303d0,
     : 2.0965942d0, 2.1430813d0, 2.1923823d0, 2.2445704d0, 2.2997100d0,
     : 2.3578619d0, 2.4190857d0, 2.4834419d0, 2.5509930d0, 2.6218045d0,
     : 2.6959455d0, 2.7734896d0, 2.8545157d0, 2.9391080d0, 3.0273577d0,
     : 3.1193629d0, 3.2152298d0, 3.3150737d0, 3.4190200d0, 3.5272055d0,
     : 3.6397795d0, 3.7569061d0, 3.8787657d0, 4.0055575d0, 4.1375022d0,
     : 4.2748459d0, 4.4178635d0, 4.5668641d0, 4.7221978d0, 4.8842628d0,
     : 5.0535163d0, 5.2304871d0, 5.4157932d0, 5.6101653d0, 5.8144773d0,
     : 6.0297921d0, 6.2574271d0, 6.4990483d0, 6.7568269d0, 7.0336896d0,
     : 7.3337690d0, 7.6632828d0, 8.0324955d0, 8.4611051d0, 8.9999999d0/
      data vl0/
     : 1.2607681d-1, 1.0181059d-1, 8.3219741d-2, 6.8419479d-2,
     : 5.6377898d-2, 4.6477862d-2, 3.8225429d-2, 3.1260322d-2,
     : 2.5331108d-2, 2.0257419d-2, 1.5918347d-2, 1.2230841d-2,
     : 9.1481765d-3, 6.6182707d-3, 4.5563309d-3, 2.9331860d-3,
     : 1.7264646d-3, 8.3777236d-4, 1.9955285d-4, -2.4070280d-4,
     : -5.0774833d-4, -6.2937775d-4, -6.6622024d-4, -6.5687384d-4,
     : -6.2561154d-4, -5.7495144d-4, -5.0879966d-4, -4.3843666d-4,
     : -3.6993741d-4, -3.0738843d-4, -2.5257993d-4, -2.0600623d-4,
     : -1.6706238d-4, -1.3487276d-4, -1.0847183d-4, -8.6948570d-5,
     : -6.9483347d-5, -5.5365330d-5, -4.3988979d-5, -3.4845727d-5,
     : -2.7512814d-5, -2.1641886d-5, -1.6947943d-5, -1.3199147d-5,
     : -1.0207681d-5, -7.8216457d-6, -5.9179573d-6, -4.3958332d-6,
     : -3.1695421d-6, -2.1530930d-6/
      data vl1/
     : 6.2337068d-2, 4.7591108d-2, 3.7269165d-2, 2.9624747d-2,
     : 2.3632515d-2, 1.8771436d-2, 1.4693083d-2, 1.1218873d-2,
     : 8.2658172d-3, 5.7899159d-3, 3.7611734d-3, 2.1552040d-3,
     : 9.5100145d-4, 9.5185923d-5, -5.1017525d-4, -9.0269916d-4,
     : -1.0986212d-3, -1.1721675d-3, -1.1622529d-3, -1.1028322d-3,
     : -9.9920733d-4, -8.5534643d-4, -7.1147111d-4, -5.8145498d-4,
     : -4.7215893d-4, -3.8192287d-4, -3.0632336d-4, -2.4433462d-4,
     : -1.9312747d-4, -1.5095284d-4, -1.1683117d-4, -8.9830877d-5,
     : -6.8769523d-5, -5.2494115d-5, -3.9998484d-5, -3.0438966d-5,
     : -2.3140927d-5, -1.7573786d-5, -1.3327579d-5, -1.0088078d-5,
     : -7.6159162d-6, -5.7291759d-6, -4.2896354d-6, -3.1921149d-6,
     : -2.3563843d-6, -1.7210496d-6, -1.2389474d-6, -8.7360481d-7,
     : -5.9627309d-7, -3.8223270d-7/
      data vl2/ 
     : 1.7393770d-2, 2.5275176d-3, -5.8311630d-3, -1.0456782d-2,
     : -1.3063412d-2, -1.4465414d-2, -1.5165738d-2, -1.5371997d-2,
     : -1.5161998d-2, -1.4571782d-2, -1.3641829d-2, -1.2441135d-2,
     : -1.1047491d-2, -9.5612948d-3, -8.1115558d-3, -6.7552641d-3,
     : -5.5018790d-3, -4.3905383d-3, -3.4325067d-3, -2.6470278d-3,
     : -2.0100088d-3, -1.4864677d-3, -1.0802168d-3, -7.6391980d-4,
     : -5.2068177d-4, -3.4419469d-4, -2.2475817d-4, -1.4999897d-4,
     : -1.0391867d-4, -7.4811107d-5, -5.5344378d-5, -4.1749021d-5,
     : -3.1844832d-5, -2.4459411d-5, -1.8874535d-5, -1.4621512d-5,
     : -1.1364536d-5, -8.8564121d-6, -6.9127792d-6, -5.3969856d-6,
     : -4.2082202d-6, -3.2721075d-6, -2.5332462d-6, -1.9496783d-6,
     : -1.4889789d-6, -1.1255877d-6, -8.3899891d-7, -6.1244014d-7,
     : -4.3161547d-7, -2.8210941d-7/
      data vl3/ 
     : 1.5214532d-1, 1.0273179d-1, 7.1378058d-2, 5.0998610d-2,
     : 3.7191821d-2, 2.7584935d-2, 2.0513913d-2, 1.5051002d-2,
     : 1.0714285d-2, 7.2476051d-3, 4.5189966d-3, 2.4446832d-3,
     : 9.5612096d-4, -4.3454084d-5, -6.9561583d-4, -1.0874717d-3,
     : -1.2632002d-3, -1.2826289d-3, -1.1844718d-3, -1.0473498d-3,
     : -9.0512737d-4, -7.5817675d-4, -6.3803882d-4, -5.2772550d-4,
     : -4.1822946d-4, -3.2114675d-4, -2.4151544d-4, -1.8030674d-4,
     : -1.3469949d-4, -1.0096990d-4, -7.5909797d-5, -5.7266582d-5,
     : -4.3317181d-5, -3.2812571d-5, -2.4863470d-5, -1.8822196d-5,
     : -1.4220022d-5, -1.0710713d-5, -8.0362607d-6, -6.0015005d-6,
     : -4.4574350d-6, -3.2894980d-6, -2.4093127d-6, -1.7486608d-6,
     : -1.2549562d-6, -8.8773740d-7, -6.1591770d-7, -4.1558329d-7,
     : -2.6811763d-7, -1.5796514d-7/
      data vl4/ 
     : 8.0888133d-2, 5.2183065d-2, 3.4496506d-2, 2.3523903d-2,
     : 1.6550682d-2, 1.2084566d-2, 9.0641862d-3, 6.8841494d-3,
     : 5.2216955d-3, 3.9079652d-3, 2.8659579d-3, 2.0536899d-3,
     : 1.4520241d-3, 1.0279244d-3, 7.0780717d-4, 4.5712096d-4,
     : 2.7919470d-4, 1.7639972d-4, 1.5770130d-4, 1.5251182d-4,
     : 1.2705191d-4, 9.3800594d-5, 3.4338388d-5, -1.8997302d-5,
     : -4.8591174d-5, -6.2785397d-5, -6.5275781d-5, -6.0684204d-5,
     : -5.2590308d-5, -4.3735461d-5, -3.5262248d-5, -2.7944257d-5,
     : -2.1873387d-5, -1.6972846d-5, -1.3066662d-5, -9.9802388d-6,
     : -7.5569638d-6, -5.6675986d-6, -4.2060490d-6, -3.0854715d-6,
     : -2.2345075d-6, -1.5947380d-6, -1.1187397d-6, -7.6850833d-7,
     : -5.1402148d-7, -3.3189775d-7, -2.0417764d-7, -1.1727174d-7,
     : -6.1139992d-8, -2.9015422d-8/
      data vl5/
     : 1.4377118d-3, 1.4096455d-3, 1.3810078d-3, 1.3508445d-3,
     : 1.3187829d-3, 1.2846761d-3, 1.2482469d-3, 1.2095751d-3,
     : 1.1710141d-3, 1.1333066d-3, 1.0941225d-3, 1.0581734d-3,
     : 1.0427239d-3, 1.0452764d-3, 1.0177235d-3, 9.3954779d-4,
     : 8.2192983d-4, 7.0391484d-4, 6.5263349d-4, 6.1970048d-4,
     : 5.6545985d-4, 5.0532043d-4, 4.1993802d-4, 3.3565236d-4,
     : 2.6501995d-4, 1.9953586d-4, 1.4425151d-4, 1.0185901d-4,
     : 7.0552165d-5, 4.8010560d-5, 3.2274013d-5, 2.1455701d-5,
     : 1.4133463d-5, 9.2449181d-6, 6.0118134d-6, 3.8922208d-6,
     : 2.5107788d-6, 1.6146450d-6, 1.0350378d-6, 6.6085852d-7,
     : 4.1951850d-7, 2.6393820d-7, 1.6371558d-7, 9.9282094d-8,
     : 5.8061498d-8, 3.1976927d-8, 1.5844869d-8, 6.3525379d-9,
     : 1.4287587d-9, 0.0000000d0/
      data vl6/
     : 1.0678495d-4, 1.2736550d-4, 1.4835239d-4, 1.7044923d-4,
     : 1.9397894d-4, 2.1904228d-4, 2.4570108d-4, 2.7353557d-4,
     : 2.9977824d-4, 3.2279960d-4, 3.4281692d-4, 3.6097707d-4,
     : 3.9212596d-4, 4.3684292d-4, 4.5182586d-4, 4.0643422d-4,
     : 2.9179097d-4, 1.5343952d-4, 9.0268946d-5, 6.7406651d-5,
     : 3.9005620d-5, 2.1113254d-5, -3.3505597d-6, -1.4826673d-5,
     : -1.3824339d-5, -1.5483628d-5, -1.6216761d-5, -1.4190982d-5,
     : -1.0828174d-5, -7.6528019d-6, -4.9108677d-6, -2.8977832d-6,
     : -1.5106459d-6, -6.3148218d-7, -1.1456645d-7, 1.5729089d-7,
     : 2.7560013d-7, 3.0373137d-7, 2.8384921d-7, 2.4217164d-7,
     : 1.9410692d-7, 1.4797538d-7, 1.0772013d-7, 7.4713867d-8,
     : 4.8931506d-8, 2.9673327d-8, 1.5994693d-8, 6.9531359d-9,
     : 1.7622446d-9, 0.0000000d0/
      end
c******** Cubic splines
c
      SUBROUTINE SPLINEX(X,Y,N,YP1,YPN,Y2)
      implicit real*8(a-h,o-z)
      PARAMETER (NMAX=2000)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
      IF (YP1.GT..99d30) THEN
        Y2(1)=0.d0
        U(1)=0.d0
      ELSE
        Y2(1)=-0.5d0
        U(1)=(3.d0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.d0
        Y2(I)=(SIG-1.d0)/P
        U(I)=(6.d0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99d30) THEN
        QN=0.d0
        UN=0.d0
      ELSE
        QN=0.5d0
        UN=(3.d0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.d0)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END
c
      SUBROUTINE SPLINT0(XA,YA,Y2A,N,X,Y)
      implicit real*8(a-h,o-z)
      DIMENSION XA(N),YA(N),Y2A(N)
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.d0) PAUSE 'Bad XA input.'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.d0
      RETURN
      END
