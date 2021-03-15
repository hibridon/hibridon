* System:  He-SH(X^2Pi) original ab initio RCCSD(T) PES's
* Reference: Cybulski et al et al J. Chem. Phys. 113, 9549 (2000).


      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(19)
      include "common/parpot"
      potnam='Cybulski et al He-SH(X) RCCSD(T) PES'
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
      potnam='Cybulski et al He-SH(X) RCCSD(T) PES'
      lammin(1)=1
      lammax(1)=10
      lammin(2)=2
      lammax(2)=10
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
*    vvl:     vector of length 6 to store r-dependence of each term
*             in potential expansion
*    vvl(1-10) expansion coefficients in dl0 (l=1:10) of vsum
*    vvl(11-19) expansion coefficients in dl2 (l=2:10) of vdif

* uses linear least squares routines from cmlib

* revised for RCCSD(T) He-SH PES: 2008.02.26 by Jacek Klos
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
      data rmax /25d0/
* coefficicients for d0 rotation matrices
* stored (by column) for each of 7 angles and for l=0:6
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
* stored (by column) for each of 5 angles and for l=2:6
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

      do i=1,11
      cthta(i)=dcos(thta(i)*dacos(-1.D0)/180.D0)
      enddo
* determine A' and A" potentials at angles
      rang=r*0.529177259D0
      do 100 i=1,11
        vsum(i)=PHEHSS(rang,cthta(i))
* don't compute vdif for colinear geometries
        if (i.ne.1 .and. i.ne.11) then
C MINUS SIGN IS TO CORRECT FOR OPPOSITE DEFINITION
C OF VDIF USED BY CYBULSKI ET AL IN THEIR CODE
          vdif(i-1)=-PHEHSD(rang,cthta(i))
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
*      conv=1.d0/219474.6
      conv=1.0D-3
      call dscal(11,conv,xsum,1)
      vv0=xsum(1)
      call dcopy(lsum-1,xsum(2),1,vvl,1)
      call dcopy(81,d2,1,aa,1)
      call dqrank(aa,9,9,ldif,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,9,9,ldif,kr,vdif,xdif,rsd,kpvt,qraux)
      call dscal(9,conv,xdif,1)
      call dcopy(ldif,xdif,1,vvl(11),1)
      end

C----He-SH POTENTIAL DECKS------


*Deck Phehss
      FUNCTION PHEHSS  ( R, C )
      IMPLICIT REAL*8  ( A-H, O-Z )
      REAL*8             A( 6 )
      PARAMETER        ( ME = 6, ML = 6, MT = 2 * ME + 4 * ML + 4 )
      REAL*8             X( MT )
*
*     Subroutine for evaluating the potential energy surface for
*     the He-HS (X state) obtained by fitting the RCCSD(T) results of
*     S.M. Cybulski and described in J. Chem. Phys. 113, 9549 (2000).
*     This function evaluates  V = 0.5 ( V(A') + V(A") )
*
*     On input:
*     =========
*     R - intermolecular distance (in angstroms)
*     C - cos( theta ), where theta is the angle between the bond axis
*         of HS and a vector from He to the center of mass of HS.
*
*     On return:
*     ==========
*     PHEHSS - interaction energy in mE_h (milihartrees)
*
      DATA               X  /
     *  0.10133614D+02, 0.14340552D+01, 0.14377761D+00, 0.32561888D+00,
     *  0.63067982D-01,-0.10810939D-01,-0.31281365D+01,-0.55839881D+00,
     *  0.10585865D+00,-0.40375325D-01,-0.17155868D-01, 0.12530605D-01,
     *  0.70039030D+00,-0.22134900D+01, 0.24431613D+00,-0.12967997D+01,
     * -0.28443694D+00, 0.12363489D+00, 0.19810241D+00, 0.16214110D+01,
     * -0.14352046D+00, 0.14679245D+01, 0.24134805D+00,-0.17565422D+00,
     * -0.55995385D-01,-0.25231721D+00, 0.22362101D-01,-0.44262752D+00,
     * -0.36384450D-02, 0.87707493D-01,-0.36374967D-02, 0.20090256D-02,
     * -0.28156997D-02, 0.39363142D-01,-0.81674589D-02,-0.12338368D-01,
     * -0.42036767D+03,-0.27433242D+02,-0.46336501D+03,-0.32584188D+03/
*
      PARAMETER        ( ONE = 1.0D0, THREE = 3.0D0, F15 = 15.0D0,
     *                   ZERO = 0.0D0, FIVE = 5.0D0, THIRTY = 30.0D0,
     *                   SEVNTY = 70.0D0, HALF = 0.5D0, T35 = 35.0D0,
     *                   S63 = 63.0D0, THRESH = 1.0D-8, O125 = 0.125D0 )
*
 2000 FORMAT( ' No convergence in F', I1, ' !' )
*
*     Start by evaluating the required Legendre polynomials
*
      T = C
      T2 = T * T
      T3 = T2 * T
      T4 = T3 * T
      T5 = T4 * T
      T6 = T5 * T
      A( 1 ) = ONE
      A( 2 ) = T
      A( 3 ) = HALF * ( THREE * T2 - ONE )
      A( 4 ) = HALF * ( FIVE * T3 - THREE * T )
      A( 5 ) = O125 * ( T35 * T4 - THIRTY * T2 + THREE )
      A( 6 ) = O125 * ( S63 * T5 - SEVNTY * T3 + F15 * T )
*
*     VSH( R, THETA ) part:
*
      D = ZERO
      B = ZERO
      DO 100 I = 1, ME
         D = D + X( I ) * A( I )
         B = B + X( I + ME ) * A( I )
 100  CONTINUE
      NE = 2 * ME
      R2 = R * R
      R3 = R2 * R
      G = ZERO
      DO 200 I = 1, ML
         IE = I + NE
         G = G + ( X( IE ) + X( IE + ML ) * R + X( IE + 2 * ML )
     *         * R2 + X( IE + 3 * ML ) * R3 ) * A( I )
 200  CONTINUE
      VSH = G * EXP( D + B * R )
*
*     VAS( R, THETA ) part:
*
      BR = B * R
      S = ONE
      T = ONE
      DBR = DABS( BR )
      DO 300 K = 1, 6
         T = T * DBR / DFLOAT( K )
         S = S + T
 300  CONTINUE
      F6 = ONE - DEXP( BR ) * S
      T2 = T * DBR / DFLOAT( 7 )
      S = S + T2
      F7 = ONE - DEXP( BR ) * S
*
*     If F6 or F7 are smaller than 1.0D-8, they should be recalculated
*     and the segment below will do the trick.
*
      IF( DABS( F6 ) .LT. THRESH ) THEN
         F6 = ZERO
         DO 400 I = 7, 1000
            T = T * DBR / DFLOAT( I )
            F6 = F6 + T
            IF( ( T / F6 ) .LT. THRESH ) GO TO 500
 400     CONTINUE
         WRITE( *, 2000 ) I - 1
         STOP
 500     CONTINUE
         F6 = F6 * DEXP( BR )
      END IF
*
      IF( DABS( F7 ) .LT. THRESH ) THEN
         F7 = ZERO
         DO 600 I = 8, 1000
            T2 = T2 * DBR / DFLOAT( I )
            F7 = F7 + T2
            IF( ( T2 / F7 ) .LT. THRESH ) GO TO 700
 600     CONTINUE
         WRITE( *, 2000 ) I - 1
         STOP
 700     CONTINUE
         F7 = F7 * DEXP( BR )
      END IF
*
      R6 = R3 * R3
      R7 = R6 * R
      NE = NE + 4 * ML
      VAS = F6 * ( X( 1 + NE ) + X( 2 + NE ) * A( 3 ) ) / R6
     *    + F7 * ( X( 3 + NE ) * A( 2 ) + X( 4 + NE ) * A( 4 ) ) / R7
      PHEHSS = VSH + VAS
      RETURN
      END


*Deck Phehsd
      FUNCTION PHEHSD  ( R, C )
      IMPLICIT REAL*8  ( A-H, O-Z )
      PARAMETER        ( ME = 6, ML = 6, MT = 2 * ME + 4 * ML + 8 )
      REAL*8             X( MT )
      REAL*8             A( 6 ), E( 6 )
*
*     Subroutine for evaluating the potential energy surface for
*     the He-HS (X state) obtained by fitting the RCCSD(T) results of
*     S.M. Cybulski and described in J. Chem. Phys. 113, 9549 (2000).
*     This function evaluates  V = 0.5 ( V(A') - V(A") )
*     Note that by definition V = 0 at 0 and 180 degrees.
*
*     On input:
*     =========
*     R - intermolecular distance (in angstroms)
*     C - cos( theta ), where theta is the angle between the bond axis
*         of HS and a vector from He to the center of mass of HS.
*
*     On return:
*     ==========
*     PHEHSD - interaction energy in mE_h (milihartrees)
*
*     CCSD(T) surface
      DATA               X  /
     *  -0.1137747129D+01, 0.8423153318D+00, 0.6019101450D+00,
     *  -0.8216049541D-01, 0.4074365680D-02, 0.6956475555D-02,
     *  -0.3844950019D+01, 0.5608945013D-01, 0.1008925382D-01,
     *   0.8768866182D-01, 0.1233946367D-01, 0.4738972898D-05,
     *  -0.4544786904D+04, 0.4336221617D+04,-0.5481080077D+04,
     *   0.2011465231D+02, 0.1309315465D+03, 0.2924552033D+02,
     *  -0.1671913334D+05,-0.1575914716D+04, 0.6425676382D+04,
     *  -0.3542149676D+03,-0.1418363273D+03,-0.1238991829D+02,
     *   0.8275404700D+04, 0.4010914349D+03,-0.2317400850D+04,
     *   0.1092215801D+03, 0.6279052921D+02, 0.4810261064D+01,
     *  -0.4299352255D+04, 0.3836118274D+03, 0.3216912310D+03,
     *  -0.1183439192D+02,-0.1068653075D+02,-0.1373588367D+01,
     *  -0.2906699905D+00, 0.4157552035D-01,-0.3161136359D+00,
     *  -0.1433557616D+01, 0.4159436308D+03,-0.3335013868D+01,
     *   0.4166439948D+01, 0.3784761062D+02/
*
      PARAMETER        ( ONE = 1.0D0, THREE = 3.0D0, EIGHT = 8.0D0,
     *                   ZERO = 0.0D0, FOUR = 4.0D0, D7P5 = 7.5D0,
     *                   SEVEN = 7.0D0, D15 = 15.0D0, D52P5 = 52.5D0,
     *                   D13125 = 13.125D0, D19 = 19.0D0, D33 = 33.0D0,
     *                   D51 = 51.0D0, D7875 = 7.875D0, D125 = 125.0D0,
     *                   D253 = 253.0D0, D143 = 143.0D0, HALF = 0.5D0,
     *                   THRESH = 1.0D-8, O125 = 0.125D0, FIVE = 5.0D0,
     *                   THIRTY = 30.0D0, SEVNTY = 70.0D0, D35 = 35.0D0,
     *                   D63 = 63.0D0 )
*
 2000 FORMAT( ' No convergence in F', I1, ' !' )
*
*     Quick return for 0 or 180.
*
*       IF( DABS( C ) .EQ. ONE ) THEN
*         PHEHSD = ZERO
*         RETURN
*      END IF
*
*     Start by evaluating the required Legendre polynomials
*
      T = C
      T2 = T * T
      T3 = T2 * T
      T4 = T3 * T
      T5 = T4 * T
      T6 = T5 * T
      T7 = T6 * T
      A( 1 ) = THREE * ( ONE - T2 )
      A( 2 ) = D15 * ( T - T3 )
      A( 3 ) = D7P5 * ( EIGHT * T2 - SEVEN * T4 - ONE )
      A( 4 ) = D52P5 * ( FOUR * T3 - THREE * T5 - T )
      A( 5 ) = D13125 * ( ONE - D19 * T2 + D51 * T4 - D33 * T6 )
      A( 6 ) = D7875 * ( D15 * T - D125 * T3 + D253 * T5 - D143 * T7 )
      E( 1 ) = ONE
      E( 2 ) = T
      E( 3 ) = HALF * ( THREE * T2 - ONE )
      E( 4 ) = HALF * ( FIVE * T3 - THREE * T )
      E( 5 ) = O125 * ( D35 * T4 - THIRTY * T2 + THREE )
      E( 6 ) = O125 * ( D63 * T5 - SEVNTY * T3 + D15 * T )
*
*     VSH( R, THETA ) part:
*
      D = ZERO
      B = ZERO
      DO 100 I = 1, ME
         D = D + X( I ) * E( I )
         B = B + X( I + ME ) * E( I )
 100  CONTINUE
      NE = 2 * ME
      R2 = R * R
      R3 = R2 * R
      G = ZERO
      DO 200 I = 1, ML
         IE = I + NE
         G = G + ( X( IE ) + X( IE + ML ) * R + X( IE + 2 * ML )
     *         * R2 + X( IE + 3 * ML ) * R3 ) * A( I )
 200  CONTINUE
      VSH = G * EXP( D + B * R )
*
*     VAS( R, THETA ) part:
*
      BR = B * R
      S = ONE
      T = ONE
      DBR = DABS( BR )
      DO 300 K = 1, 6
         T = T * DBR / DFLOAT( K )
         S = S + T
 300  CONTINUE
      F6 = ONE - DEXP( BR ) * S
      T2 = T * DBR / DFLOAT( 7 )
      S = S + T2
      F7 = ONE - DEXP( BR ) * S
      T3 = T2 * DBR / DFLOAT( 8 )
      S = S + T3
      F8 = ONE - DEXP( BR ) * S
      T4 = T3 * DBR / DFLOAT( 9 )
      S = S + T4
      F9 = ONE - DEXP( BR ) * S
*
*     If F6, ..., F9 are smaller than 1.0D-8, they should be
*     recalculated and the segment below will do the trick.
*
      IF( DABS( F6 ) .LT. THRESH ) THEN
         F6 = ZERO
         DO 400 I = 7, 1000
            T = T * DBR / DFLOAT( I )
            F6 = F6 + T
            IF( ( T / F6 ) .LT. THRESH ) GO TO 500
 400     CONTINUE
         WRITE( *, 2000 ) I - 1
         STOP
 500     CONTINUE
         F6 = F6 * DEXP( BR )
      END IF
*
      IF( DABS( F7 ) .LT. THRESH ) THEN
         F7 = ZERO
         DO 600 I = 8, 1000
            T2 = T2 * DBR / DFLOAT( I )
            F7 = F7 + T2
            IF( ( T2 / F7 ) .LT. THRESH ) GO TO 700
 600     CONTINUE
         WRITE( *, 2000 ) I - 1
         STOP
 700     CONTINUE
         F7 = F7 * DEXP( BR )
      END IF
*
      IF( DABS( F8 ) .LT. THRESH ) THEN
         F8 = ZERO
         DO 800 I = 9, 1000
            T3 = T3 * DBR / DFLOAT( I )
            F8 = F8 + T3
            IF( ( T3 / F8 ) .LT. THRESH ) GO TO 900
 800     CONTINUE
         WRITE( *, 2000 ) I - 1
         STOP
 900     CONTINUE
         F8 = F8 * DEXP( BR )
      END IF
*
      IF( DABS( F9 ) .LT. THRESH ) THEN
         F9 = ZERO
         DO 910 I = 10, 1000
            T4 = T4 * DBR / DFLOAT( I )
            F9 = F9 + T4
            IF( ( T4 / F9 ) .LT. THRESH ) GO TO 920
 910     CONTINUE
         WRITE( *, 2000 ) I - 1
         STOP
 920     CONTINUE
         F9 = F9 * DEXP( BR )
      END IF
*
      R6 = R3 * R3
      R7 = R6 * R
      R8 = R7 * R
      R9 = R8 * R
      NE = NE + 4 * ML
      VAS = F6 * ( X( 1 + NE ) * A( 1 ) + X( 2 + NE ) * A( 3 ) ) / R6
     *    + F7 * ( X( 3 + NE ) * A( 2 ) + X( 4 + NE ) * A( 4 ) ) / R7
     *    + F8 * ( X( 5 + NE ) * A( 1 ) + X( 6 + NE ) * A( 3 ) ) / R8
     *    + F9 * ( X( 7 + NE ) * A( 2 ) + X( 8 + NE ) * A( 4 ) ) / R9
      PHEHSD = VSH + VAS
      RETURN
      END


