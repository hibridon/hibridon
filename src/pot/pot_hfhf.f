* System:  HF-HF
* Reference:  M. H. Alexander and A. E. DePristo, J. Chem. Phys. 65, 5009 (1976) P. F. Vohralik, R. O. Watts, and M. H. Alexander, J. Chem. Phys. 91, 7563 (1989); 93, 3983 (1990).


      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(3)
      common /coconv/ econv
      include "common/parpot"
      potnam='ALEXANDER-DEPRISTO HF-HF 3-TERM'

      econv=219474.6d0
1      print *, ' r (bohr)'
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(4(1pe16.8))
      goto 1
99    r=4
      write (2,*) 'HF-HF 3-term PES'
      do i=1,40
         call pot(vv0,r)
         write (2,101) r,vv0,(vvl(j), j=1,3)
101      format(f8.4,4(pe16.8))
         r=r+0.2d0
      enddo
      end

      include "common/syusr"
      include "common/bausr"
      include "common/ground"
* --------------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parpot"
      include "common/parbas"
      potnam='ALEXANDER-DEPRISTO HF-HF 3-TERM'
      lammin(1)=2
      lammax(1)=2
      mproj(1)=0
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* --------------------------------------------------------------------
      subroutine pot (vv0, r)
* --------------------------------------------------------------------
*  subroutine to calculate the r-dependent terms or their derivatives
*  for the alexander-depristo hf-hf rigid rotor potential
*  nlam is the number of anisotropic terms included
*    if nlam = 1, then just dipole-dipole
*    if nlam = 2, then dipole-dipole and dipole-quadrupole
*    if nlam = 3, then dipole-dipole, dipole-quadrupole, and short-range
*  on return:
*   vv0 contains the isotropic term
*   the coefficients for each angular term in the coupling potential
*   [ vvl(i) for i = 1, nlam ] are returned in common block /covvl/
*   vvl(1) contains the dipole-dipole term
*   vvl(2) contains the dipole-quadrupole term (if desired)
*   vvl(3) contains the short-range term (if desired)
*  variable in common block /conlam/
*    nlam:      the total number of angular coupling terms
*  variable in common block /covvl/
*    vvl        array to store r-dependence of each angular term in the
*               potential
* --------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      common /conlam/ nlam
      common /covvl/ vvl(3)
*  factor is the conversion from kcal/mole to hartree
      factor = 1.5935636d-3
      rm1 = 1.d0 / r
      rm2 = rm1 * rm1
*  44.54662398 is (4*pi)**1.5
      fact = factor / 44.54662398d0
*  now calculate r-dependence of potential
*  first isotropic term
      e1 = exp( - r)
      e2 = e1 * e1
      ep1 = exp( - .1 * r)
      ex20 = 3.05d6 * e2
      ex11 = 2.25d4 * e1 * ep1
      vv0 = fact * (ex20 - ex11)
*   now dipole - dipole term
      ex = 490.0 * e1
      ex04 = 45.0 * exp( - .475 * r)
      p3 = 5233.0 * rm2 * rm1
      vvl(1)  =  - fact * (ex + ex04 + p3)
      if (nlam .ge. 2) then
*   here if dipole - quadrupole term included
        p1 = 4.2d5 * rm1
        ex2 = e2 * exp( - .04 * r)
        p4 = 1.462d4 * rm2 * rm2
        vvl(2) = fact * ((1.0d5 - p1) * ex2 + p4)
      end if
      if (nlam .eq. 3) then
*   here if short - range anisotropy included
        ex19 = 1.83d5 * e2 / ep1
        ex07 = 610.0 * exp( - .7885 * r)
        vvl(3) = fact * (ex07 - ex19)
      end if
      return
      end
