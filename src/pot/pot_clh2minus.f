*system: Cl-H2 minus using Alexander avqz  potential
*
* Reference:  M.H.Alexander, J. Chem. Phys. 118, 9637 (2003)
*
      include "common/syusr"
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(2)
      include "common/parpot"
      potnam='ALEXANDER CLH2 MINUS' 
      print *, potnam 
1     print *, ' r (bohr)'
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vsum',/,7(1pe16.8))
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
      potnam='ALEXANDER CLH2 MINUS' 
      lammin(1)=2
      lammax(1)=4
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
      common /covvl/ vvl(2)
      dimension xlam1(3),xlam2(3),c1(3),c2(3),c3(3),v(3)


*  alexander clh2minus avqz-ci potential
      data xlam1 / 1.6406d+0,  2.8710d-1,  6.9625d-1/
      data xlam2 / 1.5952d-1,  1.1195d+0,  2.4741d+0/
      data c1 / 4.7867d+6, -1.9602d+3, -8.7879d+2/
      data c2 /-3.3148d+3,  2.2473d+5,  1.5736d+6/
      data c3 / 7.0924d+2, -7.0676d+4, -5.0293d+1/
* dissociation energy in cm-1
      de=835.65d0
      data cmtohar / 219474.6d0/
      
      do i=1,3
         v(i)=c1(i)*exp(-xlam1(i)*r)+(c2(i)+c3(i)*r)*exp(-xlam2(i)*r)
      enddo
      vv0=v(1)/cmtohar
* subtract dissociation energy
      vv0=vv0-de/cmtohar
      vvl(1)=v(2)/cmtohar
      vvl(2)=v(3)/cmtohar

      return
      end
