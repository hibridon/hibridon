* System:  NH(A 3Pi)+He, original ab initio CEPA PES's
* Reference: R. Jonas and V. Staemmler, Z. Phys. D 14, 143 (1989)

      include "common/syusr"
      include "common/bausr"
      include "common/ground"
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(9)
      common /coconv/ econv
      include "common/parpot"
      potnam='JONAS-STAEMMLER CEPA He-NH(A)'
      econv=219474.6d0
      print *, potnam
1      print *, ' r (bohr)'
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vsum',/,6(1pe16.8),/,
     :    '  vdif',/,4e16.8)
      goto 1
99    end
      subroutine loapot(iunit,filnam)
* ------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      potnam='JONAS-STAEMMLER CEPA He-NH(A)'
      lammin(1)=1
      lammax(1)=5
      lammin(2)=2
      lammax(2)=5
      mproj(1)=0
      mproj(2)=2
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* ----------------------------------------------------------------------
      subroutine pot (vv0, rz)
*  subroutine to calculate the rz-dependent coefficients in the
*  NH(3PI)-He potential from R. Jonas and V. Staemmler,
*  Z. Phys. D 14, 143 (1989)
*  in units of hartree for energy and bohr for distance
* ----------------------------------------------------------------------
*  current revision:  14-mar-1991 by mha
*  on entry:
*    rz:      interparticle distance
*  vv0 contains the isotropic term (hartree)
*  vvl(i = 1,5) contains the potential expansion coefficients for
*  the anistropic components in the sum potential for lam = 1,5
*  vvl(i = 6,9) contains the expansion coefficients for the
*  difference potential for lam = 2,5
*  variable in common block /conlam/ used here
*    nlam:    the number of angular coupling terms actually used
*  variable in common block /covvl/
*    vvl:     array to store rz-dependence of each
*             angular term in the potential
*  variable in common block /coconv/ used here
*   econv:     conversion factor from cm-1 to hartrees
* ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(9)
      common /coconv/ econv
      dimension c1(10), c2(10), c3(10), vlam(10), xlam1(10), xlam2(5)

*  these parameters are derived from m.h. alexander's fit, where
*  elements 1 to 5 correspond to v10, v40, v50, v32 and v42 in cm-1
*  (fit to [c1*exp(-lam1*R)+(c2+c3*R)exp(-lam2*R)]*D(R)) and
*  elements 6 to 10 correspond to v00, v20, v30, v22 and v52 in cm-1
*  (fit to exp(-lam*R)*(c1+c2*R+c3*R~2)*D(R))
*  (with in both cases D(R) = -0.5*(tanh(1.5*(R-7))-1) and R in bohr)
      data c1 / -5.7448844d+06, 1.8974292d+05, -1.3659195d+05,
     :          5.1199742d+04, -2.5862815d+04,
     :          2.8005552d+09, 4.8153597d+05, 4.7415214d+07,
     :          -6.2088590d+05, 1.8265885d+04/
      data c2 / 8.0120657d+02, -7.0624787d+04, 1.0697699d+03,
     :          -5.5847327d+01, 7.1874093d+01,
     :          -1.5768763d+09, -1.7910035d+05, -2.5006055d+07,
     :          6.9560799d+04, 6.3078713d+01/
      data c3 / 0., 8.1031079d+03, -158.05584, 0., 0., 2.4022923d+08,
     :          1.5734957d+04, 2.9877603d+06, -2.7406807d+02, 0./
      data xlam1 / 2.0687256, 1.1178223, 1.58125, 1.2674194, 1.3,
     :             2.964186, 1.1045815, 2.2289111, 1.3258313, 1.4666067/
      data xlam2 / 0.75860962, 0.75517578, 0.38125, 0.33920288, 0.5/
*  NOTE: since Staemmler defines vdiff=0.5*(VA'-VA"), all the signs of
*  the Vlam,2 expansion coefficients should be CHANGED from what are
*  given above, since the correct expansion is vdiff=0.5*(VA"-VA')
*  (see m.h. alexander, chem. phys. 92, 337 (1985))
*  note that theta=180 corresponds to N-H-He

*  calculate expansion coefficients from analytical expressions
*  introduce damping of all terms at large r
      rz2 = rz * rz
      drconv = - 0.5 * (tanh(1.5*(rz-7.d0)) - 1.) / econv
      do 20 i = 1, 5
        vlam(i) = drconv * (c1(i) * exp(-xlam1(i)*rz) +
     :            (c2(i) + c3(i) * rz) * exp(-xlam2(i)*rz))
        vlam(i+5) = drconv * exp(-xlam1(i+5)*rz) *
     :              (c1(i+5) + c2(i+5) * rz + c3(i+5) * rz2)
20    continue
*  here for isotropic term
      vv0 = vlam(6)
*  here for anistropic terms
      vvl(1) = vlam(1)
      vvl(2) = vlam(7)
      vvl(3) = vlam(8)
      vvl(4) = vlam(2)
      vvl(5) = vlam(3)
      vvl(6) = vlam(9)
      vvl(7) = vlam(4)
      vvl(8) = vlam(5)
      vvl(9) = vlam(10)
* now multiply by minus one to compensate for different definition of vdiff
      do 30 i=6,9
        vvl(i)=-vvl(i)
30    continue
      return
      end
