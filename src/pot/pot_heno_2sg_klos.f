* System:  NO(A^2Sigma)+He, original ab initio RCCSD(T) PES's
* Reference:
* J. Klos, M. H. Alexander, R. Hernandez-Lamoneda and T. G. Wright, 
* J. Chem. Phys. 129, 244303 (2008)

      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(6)
      include "common/parpot"
      potnam='Klos et al He-NO(A 3ssigma) RCCSDT PES'
      print *, potnam
1      print *, ' r (bohr)'
      rshift=0.5
      xfact=0.8
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vsum',/,7(1pe16.8),/,
     :    '  vdif',/,5e16.8)
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
      potnam='Klos et al He-NO(A 3ssigma) RCCSDT PES'
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
*    vvl(1-6) expansion coefficients in dl0 (l=1:6) of vsum
*    vvl(7-11) expansion coefficients in dl2 (l=2:6) of vdif

* uses linear least squares routines from cmlib

* author:  millard alexander
* 0 degree for He-NO and 180 for He-ON
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension xlam1(14),xlam2(14),r0(14),c1(14),c2(14),c3(14),
     :          clr(14),vsum(7),xsum(7),vdif(7),xdif(7),
     :          ddif(7),vap(7),va2p(7),
     :          d0(49),aa(64),thta(7),cthta(7)
      dimension kpvt(8),qraux(7),work(55),rsd(7),re(7)

      common /covvl/ vvl(6)
      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
* for distances beyond rmax difference potential is damped
      data rmax /20d0/
* coefficicients for d0 rotation matrices
* stored (by column) for each of 7 angles and for l=0:6
* angles are  0 30 60 90 120 150 180
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

      thta(1)=0.D0
      thta(2)=30.D0
      thta(3)=60.D0
      thta(4)=90.D0
      thta(5)=120.D0
      thta(6)=150.D0
      thta(7)=180.D0
      do i=1,7
      cthta(i)=dcos(thta(i)*dacos(-1.D0)/180.D0)
      enddo
      do 100 i=1,7
        vsum(i)=VPESN(r,i)
100   continue
* solve simultaneous equations for solutions
* first for vsigma
      tol=1.e-10
      call dcopy(49,d0,1,aa,1)
      call dqrank(aa,7,7,7,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,7,7,7,kr,vsum,xsum,rsd,kpvt,qraux)
* convert to hartree
      conv=1.d0/219474.6
      call dscal(7,conv,xsum,1)
      vv0=xsum(1)
      call dcopy(6,xsum(2),1,vvl,1)
      end



      FUNCTION VPESN(R,N)
CSYSTEM: NO-HE 2Sigma 3s Rydb state Theta=0 for He---NO
C units R=au E=cm-1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(14),XX1(14),XX2(14),XX3(14),
     * XX4(14),XX5(14),XX6(14),XX7(14),
     * XX8(14),XX9(14),XX10(14),XX11(14)
C theta=0
      DATA XX1/
     *    .131423916900709514D+01,
     *    .176258927527663523D+01,
     *   -.813582302539818883D+08,
     *    .565857686758687615D+09,
     *   -.339658930992513120D+09,
     *    .911991590227617919D+08,
     *   -.138542527935384586D+08,
     *    .128408653557724995D+07,
     *   -.721908149771633034D+05,
     *    .232765234762263981D+04,
     *   -.333502783119149555D+02,
     *    .904264556193752966D+01,
     *    .430833583554335907D+08,
     *    .145446230894912958D+09/

C theta=30
      DATA XX2/
     *    .121559070173313355D+01,
     *    .243314311824340868D+01,
     *    .492146069840028882D+09,
     *   -.639903181824158132D+07,
     *   -.105391707985159844D+09,
     *    .379821367188983858D+08,
     *   -.636734978125522565D+07,
     *    .613264650978507358D+06,
     *   -.343477739806441678D+05,
     *    .106296997246434307D+04,
     *   -.143988510770430196D+02,
     *    .904264556193752966D+01,
     *    .444342984539582953D+08,
     *    .156461260822119564D+09/

C theta=60
      DATA XX3/       
     *    .126718513074087769D+01,
     *    .299041547167617949D+01,
     *    .101797407616952360D+10,
     *   -.375906920829937518D+09,
     *    .108111143785015214D+07,
     *    .249200366993475743D+08,
     *   -.639977216553394403D+07,
     *    .837033757021745783D+06,
     *   -.607898811685022883D+05,
     *    .245183244478681900D+04,
     *   -.428791701636777987D+02,
     *    .904264556193752966D+01,
     *    .427071144516928568D+08,
     *    .269938644874472260D+09/

C theta=90
      DATA XX4/       
     *    .124120178937217984D+01,
     *    .201737429086422404D+01,
     *    .590388146604815125D+09,
     *   -.359068349030894816D+09,
     *    .940792706746427417D+08,
     *   -.129421259234622847D+08,
     *    .902206084164039814D+06,
     *    .674312333011266446D+03,
     *   -.421464194175796274D+04,
     *    .276248405895754900D+03,
     *   -.617148936880593713D+01,
     *    .904264556193752966D+01,
     *    .392148021343831122D+08,
     *    .400194990152434886D+09/

C theta=120
      DATA XX5/       
     *    .125785690057248645D+01,
     *    .143439623934946070D+01,
     *    .236275413429721802D+09,
     *   -.838210605832316130D+08,
     *   -.290120441367796622D+07,
     *    .673713638180521876D+07,
     *   -.158949763355272659D+07,
     *    .194348040938935650D+06,
     *   -.132393877090626647D+05,
     *    .497674056588049439D+03,
     *   -.813664059061183309D+01,
     *    .904264556193752966D+01,
     *    .397027373591580763D+08,
     *    .398071187419760346D+09/

C theta=150
      DATA XX6/
     *    .124024779759014381D+01,
     *    .234652308329455117D+01,
     *    .515664904199745655D+09,
     *   -.156801506967336182D+08,
     *   -.116763717651413396D+09,
     *    .455496506216331124D+08,
     *   -.819598005539610237D+07,
     *    .849924761059873155D+06,
     *   -.514055280826786257D+05,
     *    .172103872880607105D+04,
     *   -.249778669371362874D+02,
     *    .904264556193752966D+01,
     *    .399179635667340383D+08,
     *    .384341745920602202D+09/

C theta=180
      DATA XX7/
     *    .121079289011488100D+01,
     *    .320591036076181712D+01,
     *    .299917195642765188D+10,
     *   -.150789668296675634D+10,
     *    .259233729875880092D+09,
     *   -.457474546157879196D+07,
     *   -.457900224571545981D+07,
     *    .749919031174719799D+06,
     *   -.545322199640776744D+05,
     *    .200998468844074091D+04,
     *   -.309028438740895552D+02,
     *    .904264556193752966D+01,
     *    .389613380060090870D+08,
     *    .397455165480571270D+09/


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
      TERM4 = XX(13)/(R**6)+XX(14)/(R**7)
       TOCM=0.219474643545745D0
       VPESN = (TERM1*TERM2-TERM3*TERM4)*TOCM
       RETURN
       END


