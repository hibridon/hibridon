*  H2(a3Sig(u)-) AVQZ MRCI potential
*  computed by P. Dagdigian, Jul-2013 (molpro 2010) 
*
      include "common/syusr"
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(1)
      include "common/parpot"
      econv=219474.6d0
      potnam='H2(a) MRCI AVQZ POTENTIAL'
      print *, potnam
1      print *, ' r (bohr)'
      read (5, *, end=99) r
      if (r.le.0.d0) goto 99
      call pot(vv0,r)
      write (6, 100) vv0 * econv
100   format(' v =',1pe16.8)
      goto 1
99    rr=2.6d0
      dr=0.25d0
      open (unit=12,file='h2a_pec.txt')
      write(12,209)
209   format(' %R/bohr    V(H2a)')
      do i=1,201
        call pot(vv0,rr)
        write (12,110) rr,econv*vv0
110     format(f7.2,4(1pe16.8))
        rr = rr + dr
      enddo
      close(12)
      end
      include "common/bausr"
      include "common/ground"
* --------------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parpot"
      include "common/parbas"
      potnam='H2(a) MRCI AVQZ POTENTIAL'
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

*  -----------------------------------------------------------------------
*  H2(a3Sig(u)-) AVQZ MRCI potential 
      implicit double precision (a-h,o-z)
      common /conlam/ nlam,nlammx
      include "common/parbas"
      common /covvl/ vvl(1)
      dimension rr(29),vmrci(29)
      dimension csplin(29),vec(29)
* hyperbolic tangent scaling factor
      data alph /1.2d0/
      data rmax /13d0/
*
      data ifirst /0/
*  grid of computed interaction energies
      data rr /0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8,
     :  1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 3, 3.5, 4, 4.5, 5, 6, 
     :  7, 8, 9, 10, 15, 20/
      data vmrci /98224.98, 82953.7, 71116.99, 61663.24, 53904.05, 
     :  47384.09, 41799.72, 36946.8, 32685.88, 28919.27, 25575.91, 
     :  22601.82, 19954, 17596.76, 15499.42, 13635.06, 11979.63, 
     :  6148.79, 3049.4, 1456.87, 667.37, 290.5, 42.57, 0.06, -4.01, 
     :-2.74, -1.59, -0.11, 0/
*
* spline fit of the vv0 coefficient
      if (ifirst .eq. 0) then
*    evaluate derivative at first point
        der1=(vmrci(2)-vmrci(1))/(rr(2)-rr(1))
        call dspline(rr,vmrci,29,der1,0d0,csplin)
        ifirst = 1
      end if
      call dcopy(29,vmrci,1,vec,1)
      call dsplint(rr,vec,csplin,29,r,vvx)
*  switch to 1/r^6 for long-range
*  c^ from potential at r = 15
      c6 = -1.253d06
      switch_lr=0.5*(tanh(1.0*(r - 15.d0)) + 1.d0)
* merge with asymptotic form
      vvx = vvx*(1.d0 - switch_lr) + switch_lr*c6/(r**6)
*  convert to hartree
      vv0 = vvx/219474.6d0
      vvl(1) = 0.d0
      return
      end
