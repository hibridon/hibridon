*  H-He potential
*  Molpro calculation - avqz basis + mid-bon functions
*  RCCSD(T) calculation with BSSE correction
*
*  p.dagdigian (23-sep-2014) 
*
      include "common/syusr"
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(1)
      include "common/parpot"
      econv=219474.6d0
      potnam='H-He RCCSD(T) potential'
      print *, potnam
1      print *, ' r (bohr)'
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0 * econv
100   format(' v =',1pe16.8)
      goto 1
99    end
      include "common/bausr"
      include "common/ground"
* --------------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parpot"
      include "common/parbas"
      potnam='H-He RCCSD(T) potential'
      lammin(1)=2
      lammax(1)=2
      mproj(1)=0
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end

      subroutine pot (vv0, r)
*  -----------------------------------------------------------------------

*  subroutine to calculate the r-dependent coefficients in the
*  collision of a homonuclear diatomic with a structureless target
*  in units of hartree for energy and bohr for distance

*  on return:
*  vv0 contains the isotropic term (n=0) in the potential
*  the coefficients for each angular term in the coupling potential
*  [ vvl(i) for i = 1, nlam ] are returned in common block /covvl/
*  vvl(1) contains the anisotropic (n=2) term in the potential

*  variable in common block /conlam/
*    nlammx:    the maximum number of anisotropic terms
*    nlam:      the total number of angular coupling terms
*  variable in common block /covvl/
*    vvl        array to store r-dependence of each angular term in the
*               potential
*
*  H-He potential.  Molpro calculation - avqz basis + mid-bon functions
*  RCCSD(T) calculation with BSSE correction
*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      common /conlam/ nlam,nlammx
      include "common/parbas"
      common /covvl/ vvl(1)
      dimension rrr(29), vvv(29), cs(29)
      data rrr /
     +  2.00000,  2.25000,  2.50000,  2.75000,  3.00000, 
     +  3.25000,  3.50000,  3.75000,  4.00000,  4.25000, 
     +  4.50000,  4.75000,  5.00000,  5.25000,  5.50000, 
     +  6.00000,  6.50000,  7.00000,  7.50000,  8.00000, 
     +  8.50000,  9.00000,  9.50000, 10.00000, 11.00000, 
     + 12.00000, 13.00000, 15.00000, 20.00000 /
      data vvv /
     +  1.3700850e+04, 9.3789900e+03, 6.3572100e+03, 4.2602900e+03, 
     + 2.8201500e+03, 1.8425800e+03, 1.1872000e+03, 7.5341000e+02, 
     + 4.7004000e+02, 2.8744000e+02, 1.7145000e+02, 9.8910000e+01, 
     + 5.4340000e+01, 2.7540000e+01, 1.1850000e+01, -1.6900000e+00, 
     + -4.8100000e+00, -4.6100000e+00, -3.6200000e+00, -2.6600000e+00, 
     + -1.9000000e+00, -1.3500000e+00, -9.6000000e-01, -6.9000000e-01, 
     + -3.7000000e-01, -2.0000000e-01, -1.1000000e-01, -3.0000000e-02, 
     + 2.0000000e-02 /   
      econv=219474.6d0
*
*  spline fits
      if (ifirst .eq. 0) then
*  spline fit of the coefficients for the PE curves
*  evaluate derivative at first point
        der1=(vvv(2)-vvv(1))/(rrr(2)-rrr(1))
        call dspline(rrr,vvv,28,der1,0.d0,cs)
        ifirst = 1
      endif
*
*  determine splined coefficients for PE curves at r=R
      call dsplint(rrr,vvv,cs,28,r,vvx)
      vv0 = vvx
*  fix long-range behavior of PE curve as c6/r^6 from next to last point
      if (r.gt.15.d0) then
        c6 = vvv(28) * rrr(28)**6
        vv0 = c6 / r**6
      endif
      vv0 = vv0 / econv
      vvl(1) = 0.d0
      return
      end
