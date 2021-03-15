*system: Cl-H2 minus using Alexander avqz ccsd(t)  potential
*
* Reference:  M.H.Alexander, J. Chem. Phys. 118, 9637 (2003)
*
      include "common/syusr"
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(2)
      include "common/parpot"
      potnam='ALEXANDER CLH2 MINUS CCSD(T)' 
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
      potnam='ALEXANDER CLH2 MINUS CCSD(T)' 
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
      dimension xlam1(3),xlam2(3),r0(3),c1(3),c2(3),c3(3),v(3),cl(3)


*  alexander clh2minus avqz-ci potential (v=0)
*     data xlam1 / 1.5779d+0, 2.8537d-1, 4.9314d-1/
*     data xlam2 / 5.4004d-1, 1.1025d+0, 1.6474d+0/
*     data r0 /9.3310d+0, 0d0, 7.7767d+0/
*     data c1 / 3.8705d+6, -2.0407d+3, -1.5539d+2/
*     data c2 /-1.4550d+4, 2.0786d+5, 4.1522d+5/
*     data c3 / 1.5073d+1, -6.9288d+4, -9.2548d+4/
*     data cl / 8.2506d+6,0d0, -1.3871d+5/
*  alexander clh2minus avqz-ci potential (v=1)
      data xlam1 / 1.5025d+0, 2.8811d-1, 5.4117d-1/
      data xlam2 / 5.8193d-1, 1.1017d+0, 1.7252d+0/
      data r0 /9.1616d+0, 0d0, 7.7211d+0/
      data c1 / 2.8709d+6, -2.3640d+3, -3.5411d+2/
      data c2 /-2.2115d+4, 2.0327d+5, 8.0065d+5/
      data c3 / 7.3802d+0, -7.5382d+4, -1.7656d+5/
      data cl / 1.4273d+7, 0d0, -1.9187d+5/
      data one,half,alph /1d0,0.5d0,1.2d0/
      data cmtohar / 219474.6d0/
      rm6=r**(-6) 
      do i=1,3
         if (i.eq.1 .or. i.eq.3) then
            v(i)=c1(i)*dexp(-xlam1(i)*r)+
     :        (c2(i)+c3(i)*r)*dexp(-xlam2(i)*r)-
     :        half*(tanh(alph*(r-r0(i)))+one)*cl(i)*rm6
         else
            v(i)=c1(i)*dexp(-xlam1(i)*r)+
     :        (c2(i)+c3(i)*r)*dexp(-xlam2(i)*r)
         endif
      enddo
      vv0=v(1)/cmtohar
      vvl(1)=v(2)/cmtohar
      vvl(2)=v(3)/cmtohar

      return
      end
