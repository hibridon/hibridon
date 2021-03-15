* System:  NO(X 2Pi)+He, original ab initio RCCSD(T) PES's
* Reference: J. Klos, G. Chałasiński, M. T. Berry, R. Bukowski, 
* and S. M. Cybulski, J. Chem. Phys. 112, 2195 (2000).


      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(11)
      include "common/parpot"
      potnam='Klos et al He-NO(X) CCSDT PES'
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
      potnam='Klos et al He-NO(X) CCSDT PES'
** ccsdt potential using avtz+332 basis  for 7 angles **
      lammin(1)=1
      lammax(1)=6
      lammin(2)=2
      lammax(2)=6
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
*    vvl(1-6) expansion coefficients in dl0 (l=1:6) of vsum
*    vvl(7-11) expansion coefficients in dl2 (l=2:6) of vdif

* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  8-oct-1993
* revised for He-NO(X) : 1-20-95 by Moonbong Yang
* revised for CCSD(T) PES: 2002.10.13 by Jacek Klos
* 0 degree for He-NO and 180 for He-ON
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension xlam1(14),xlam2(14),r0(14),c1(14),c2(14),c3(14),
     :          clr(14),vsum(7),xsum(7),vdif(7),xdif(7),
     :          ddif(7),vap(7),va2p(7),
     :          d0(49),d2(25),aa(64),thta(7),cthta(7)
      dimension kpvt(8),qraux(7),work(55),rsd(7),re(14)

      common /covvl/ vvl(11)
      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
* for distances beyond rmax difference potential is damped
      data rmax /13d0/
* coefficicients for d0 rotation matrices
* stored (by column) for each of 7 angles and for l=0:6
* angles are 0 30 60 90 120 150 180
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
* coefficicients for d2 rotation matrices
* stored (by column) for each of 5 angles and for l=2:6
* angles are 30 60 90 120 150
      data d2/
     : 1.5309311d-1,  4.5927932d-1,  6.1237244d-1,  4.5927934d-1,
     : 1.5309312d-1,
     : 2.9646353d-1,  5.1348990d-1,  1.8279042d-8, -5.1348989d-1,
     : -2.9646355d-1,
     : 4.1999000d-1,  2.2234766d-1, -3.9528471d-1,  2.2234762d-1,
     : 4.1999002d-1,
     : 4.9023048d-1, -1.6982081d-1, -2.4180900d-8,  1.6982085d-1,
     : -4.9023049d-1,
     : 4.8532921d-1, -3.4523418d-1,  3.2021721d-1, -3.4523418d-1,
     : 4.8532920d-1/


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
* determine A' and A" potentials at angles
      do 100 i=1,7
        vap(i)=POTFUNPRIM(r,cthta(i))
        va2p(i)=POTFUNBIS(r,cthta(i))
        vsum(i)=half*(va2p(i)+vap(i))
* don't compute vdif for colinear geometries
        if (i.ne.1 .and. i.ne.7) then
          vdif(i-1)=half*(va2p(i)-vap(i))
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
      call dcopy(49,d0,1,aa,1)
      call dqrank(aa,7,7,7,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,7,7,7,kr,vsum,xsum,rsd,kpvt,qraux)
* convert to hartree
      conv=1.d0/219474.6
      call dscal(7,conv,xsum,1)
      vv0=xsum(1)
      call dcopy(6,xsum(2),1,vvl,1)
      call dcopy(25,d2,1,aa,1)
      call dqrank(aa,5,5,5,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,5,5,5,kr,vdif,xdif,rsd,kpvt,qraux)
      call dscal(5,conv,xdif,1)
      call dcopy(5,xdif,1,vvl(7),1)
      end
*Deck POTFUNBIS
      FUNCTION POTFUNBIS(R,C)
      IMPLICIT REAL*8  (A-H,O-Z)
      REAL*8 A(6)
      PARAMETER (ME = 6, ML = 6,
     *  MT = 2 * ME + 4 * ML + 4 )
      PARAMETER (ONE = 1.0D0, THREE = 3.0D0,
     * F15 = 15.0D0,
     *  ZERO = 0.0D0, FIVE = 5.0D0, THIRTY = 30.0D0,
     *  SEVNTY = 70.0D0, HALF = 0.5D0, T35 = 35.0D0,
     *   S63 = 63.0D0, THRESH = 1.0D-8) 
      REAL*8    X(MT)
*
*     Potential energy function for 
*      a rare gas atom - asymmetric
*     linear molecule complexes.
*
*     On input:
*     =========
*R - intermolecular distance (in angstroms)
*C - cos( theta ), where theta is the angle between
*     the bond axis
*    of a diatomic and a vector from Rg to the
*    center of mass
*    of a diatomic
*
*  On return:
*  ==========
*  POTFUN - interaction energy in mE_h (milihartrees)
       data X/0.93729727E+01,-0.34394679,
     *0.12075059E+01,
     *   0.87512578E-01,
     *  0.24040713E+00,0.14739815E+00,-0.20236282E+01,
     *  -0.44819367E-01,
     * -0.18935935E-01, 0.43253788E-04,-0.44905825E-01,
     *  -0.42649430E-01,
     * -0.25182060E+01, 0.21820896E+01,-0.35647710E+01,
     *  0.56402415E+00,
     *  0.60600803E+01, 0.24719691E+00,0.38623219E+01,
     *   0.26295718E+00,
     *  0.23270757E+01, -0.10360823E+01,-0.39459956E+01,
     *  -0.67839402E+00,
     * -0.53872460E+00,-0.11447268E+00,-0.37394241E+00,
     *   0.24835212E+00,
     *  0.78682216E+00,0.24061793E+00,0.90072824E-02,
     *   0.39021151E-02,
     * 0.18131326E-01,-0.17611395E-01,-0.49279888E-01,
     *   -0.19927465E-01,
     * -0.12539203E+05,-0.47580271E+04,0.15366833E+05,
     *  0.51647119E+04/

*A2 state
*INSERT THE PARAMETERS HERE
*
*
*     Start by evaluating the required Legendre polynomials
*
      T = -C
      T2 = T * T
      T3 = T2 * T
      T4 = T3 * T
      T5 = T4 * T
      T6 = T5 * T
      A( 1 ) = ONE
      A( 2 ) = T
      A( 3 ) = HALF * ( THREE * T2 - ONE )
      A( 4 ) = HALF * ( FIVE * T3 - THREE * T )
      A( 5 ) = 0.125D0 * ( T35 * T4 - THIRTY * T2 + THREE )
      A( 6 ) = 0.125D0 * ( S63 * T5 - SEVNTY * T3 + F15 * T )
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
      T = T * DBR / DFLOAT( 7 )
      S = S + T
      F7 = ONE - DEXP( BR ) * S
*
*     If F6 or f7 are smaller than 1.0D-8,
*       they should be recalculated
*     and the segment below will do the trick.
*
      IF( DABS( F6 ) .LT. THRESH ) THEN
         F6 = ZERO
         DO 400 I = 7, 1000
            T = T * BR / DFLOAT( I )
            F6 = F6 + T
            IF( ( T / F6 ) .LT. THRESH ) GO TO 500
 400     CONTINUE
         WRITE( 6, * ) 'No convergence for F6.'
         STOP
 500     CONTINUE
         F6 = F6 * DEXP( BR )
      END IF
*
      IF( DABS( F7 ) .LT. THRESH ) THEN
         F7 = ZERO
         DO 600 I = 8, 1000
            T = T * BR / DFLOAT( I )
            F7 = F7 + T
            IF( ( T / F7 ) .LT. THRESH ) GO TO 700
 600     CONTINUE
         WRITE( 6, * ) 'No convergence for F7.'
 700     CONTINUE
         F7 = F7 * DEXP( BR )
      END IF
*
      R6 = R3 * R3
      R7 = R6 * R
      NE = NE + 4 * ML
      VAS = F6 * ( X( 1 + NE ) + X( 2 + NE ) * A( 3 ) ) / R6
     *    + F7 * ( X( 3 + NE )
     *  * A( 2 ) + X( 4 + NE ) * A( 4 ) ) / R7
      POTFUNBIS = (VSH + VAS)*1000.D0/4.556335D0
      RETURN
      END
*Deck POTFUNPRIM
      FUNCTION POTFUNPRIM(R,C)
      IMPLICIT REAL*8  (A-H,O-Z)
      REAL*8 A(6)
      PARAMETER (ME = 6, ML = 6,
     *  MT = 2 * ME + 4 * ML + 4 )
      PARAMETER (ONE = 1.0D0, THREE = 3.0D0,
     * F15 = 15.0D0,
     *  ZERO = 0.0D0, FIVE = 5.0D0, THIRTY = 30.0D0,
     *  SEVNTY = 70.0D0, HALF = 0.5D0, T35 = 35.0D0,
     *   S63 = 63.0D0, THRESH = 1.0D-8) 
      REAL*8    X(MT)
*
*     Potential energy function for 
*      a rare gas atom - asymmetric ( case A' HeNO RCCSD(T))
*     linear molecule complexes.
*
*     On input:
*     =========
*R - intermolecular distance (in angstroms)
*C - cos( theta ), where theta is the angle between
*     the bond axis
*    of a diatomic and a vector from Rg to the
*    center of mass
*    of a diatomic
*
*  On return:
*  ==========
*  POTFUN - interaction energy in mE_h (milihartrees)
      data X/0.89614790D+01,-0.66184380D+00,
     *  0.78007803D+00,-0.60974193D-02,
     * -0.35808998D-01, -0.11214555D-01,
     * -0.20242288D+01, -0.34243587D-01,
     * -0.33588331D-01, -0.37086659D-02,
     *  0.77114573D-02,  0.24071004D-02,
     * -0.24565438D+00,  0.32383547D+01,
     * -0.18898698D+01,  0.13834898D+01,
     *  0.23578296D+01,  0.14503735D+01,
     *  0.43039404D+01,  0.22754183D+01,
     *  0.47018602D+01,  0.24381747D+00,
     * -0.25353828D+00, -0.49050007D+00,
     * -0.44107679D+00, -0.38335519D+00,
     * -0.55140365D+00,  0.42936877D-01,
     * -0.19769691D+00,  0.19501492D-01,
     * -0.74599920D-02, -0.40805756D-02,
     *  0.23613677D-07, -0.56473632D-02,
     *  0.11274139D-01, -0.33655160D-02,
     * -0.13246868D+05, -0.41907209D+04,
     *  0.20463763D+05,  0.58529744D+03/    
*A1 state HeNO
*INSERT THE PARAMETERS HERE
*
*
*     Start by evaluating the required Legendre polynomials
*
      T = -C
      T2 = T * T
      T3 = T2 * T
      T4 = T3 * T
      T5 = T4 * T
      T6 = T5 * T
      A( 1 ) = ONE
      A( 2 ) = T
      A( 3 ) = HALF * ( THREE * T2 - ONE )
      A( 4 ) = HALF * ( FIVE * T3 - THREE * T )
      A( 5 ) = 0.125D0 * ( T35 * T4 - THIRTY * T2 + THREE )
      A( 6 ) = 0.125D0 * ( S63 * T5 - SEVNTY * T3 + F15 * T )
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
      T = T * DBR / DFLOAT( 7 )
      S = S + T
      F7 = ONE - DEXP( BR ) * S
*
*     If F6 or f7 are smaller than 1.0D-8,
*       they should be recalculated
*     and the segment below will do the trick.
*
      IF( DABS( F6 ) .LT. THRESH ) THEN
         F6 = ZERO
         DO 400 I = 7, 1000
            T = T * BR / DFLOAT( I )
            F6 = F6 + T
            IF( ( T / F6 ) .LT. THRESH ) GO TO 500
 400     CONTINUE
         WRITE( 6, * ) 'No convergence for F6.'
         STOP
 500     CONTINUE
         F6 = F6 * DEXP( BR )
      END IF
*
      IF( DABS( F7 ) .LT. THRESH ) THEN
         F7 = ZERO
         DO 600 I = 8, 1000
            T = T * BR / DFLOAT( I )
            F7 = F7 + T
            IF( ( T / F7 ) .LT. THRESH ) GO TO 700
 600     CONTINUE
         WRITE( 6, * ) 'No convergence for F7.'
 700     CONTINUE
         F7 = F7 * DEXP( BR )
      END IF
*
      R6 = R3 * R3
      R7 = R6 * R
      NE = NE + 4 * ML
      VAS = F6 * ( X( 1 + NE ) + X( 2 + NE ) * A( 3 ) ) / R6
     *    + F7 * ( X( 3 + NE )
     *  * A( 2 ) + X( 4 + NE ) * A( 4 ) ) / R7
      POTFUNPRIM = (VSH + VAS)*1000.D0/4.556335D0
      RETURN
      END
