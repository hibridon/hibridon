*  Ar-Ar potential.  Taken from notes on Quantum Molecular
*  Collision Theory, by M.H.Alexander 
*
      include "common/syusr"
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(1)
      include "common/parpot"
      econv=219474.6d0
      potnam='AR-AR POTENTIAL (M.H.A. NOTES)'
      print *, potnam
1      print *, ' r (bohr)'
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0 * econv
100   format(' v +',1pe16.8)
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
      potnam='AR-AR POTENTIAL (M.H.A. NOTES)'
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
      implicit double precision (a-h,o-z)
      common /conlam/ nlam,nlammx

      include "common/parbas"
      common /covvl/ vvl(1)


*  Ar-Ar potential.  Taken from notes on Quantum Molecular
*  Collision Theory, by M.H.Alexander 

*  de and re in hartree and bohr
      data de, re, beta / 4.456e-4, 7.1053d0, .89305d0 /

      vv0 = de*(exp(-2.d0*beta*(r - re)) 
     :  - 2.d0*exp(-beta*(r - re)))
      vvl(1) = 0.d0
      return
      end
