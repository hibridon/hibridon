*  H2(X1Sig(g)+) AVQZ MRCI potential
*  computed by P. Dagdigian, Jul-2013 (molpro 2010) 
*
      include "common/syusr"
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(1)
      include "common/parpot"
      econv=219474.6d0
      potnam='H2(X) MRCI AVQZ POTENTIAL'
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
      open (unit=12,file='h2x_pec.txt')
      write(12,209)
209   format(' %R/bohr    V(H2X)')
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
      potnam='H2(X) MRCI AVQZ POTENTIAL'
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
*  H2(X1Sig+) AVQZ MRCI potential 
      implicit double precision (a-h,o-z)
      common /conlam/ nlam,nlammx
      include "common/parbas"
      common /covvl/ vvl(1)
      dimension rr(33),vmrci(33)
      dimension csplin(33),vec(33)
* hyperbolic tangent scaling factor
      data alph /1.2d0/
      data rmax /13d0/
*
      data ifirst /0/
*  grid of computed interaction energies
      data rr /0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 
     :  1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 3, 3.5, 4, 
     :  4.5, 5, 6, 7, 8, 9, 10, 15, 20/
      data vmrci /104285.36, 50875.56, 17363.14, -4196.65, -18181.1, 
     :  -27176.07, -32790.86, -36068.08, -37705.37, -38181.96, 
     :  -37834.47, -36903.95, -35565.75, -33949.05, -32150.01, 
     :  -30240.83, -28276.07, -26297.19, -24335.87, -22416.26, 
     :  -20556.77, -12534.14, -6958.17, -3576, -1741.82, -823.24, 
     :  -180.14, -41.97, -11.63, -4.08, -1.82, -0.11, 0/
*
* spline fit of the vv0 coefficient
      if (ifirst .eq. 0) then
*    evaluate derivative at first point
        der1=(vmrci(2)-vmrci(1))/(rr(2)-rr(1))
        call dspline(rr,vmrci,33,der1,0d0,csplin)
        ifirst = 1
      end if
      call dcopy(32,vmrci,1,vec,1)
      call dsplint(rr,vec,csplin,32,r,vvx)
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
