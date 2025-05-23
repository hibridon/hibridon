!system: N2-Ar using canonical Pattengill, LaBudde, Bernstein potential
! reference:  m.d. pattengill, r. a. labudde, r.b. bernstein,
!  and c.f. curtiss, j. chem. phys. 55, 5517 (1971)]
#include "common/syusr.F90"
subroutine driver
use mod_covvl, only: vvl
use mod_parpot, only: potnam=>pot_name, label=>pot_label
implicit double precision (a-h,o-z)
potnam='PATTENGILL-LABUDDE-BERNSTEIN AR-N2'
print *, potnam
1  print *, ' r (bohr)'
read (5, *, end=99) r
call pot(vv0,r)
write (6, 100) vv0,vvl
100 format(' vsum',/,7(1pe16.8))
goto 1
99 end
#include "common/bausr.F90"
#include "common/ground.F90"
! --------------------------------------------------------------------------
subroutine loapot(iunit,filnam)
use mod_parbas, only: ntv, ivcol, ivrow, lammin, lammax, mproj
use mod_parpot, only: potnam=>pot_name, label=>pot_label
! --------------------------------------------------------------------------
character*(*) filnam
UNUSED_DUMMY(iunit)
UNUSED_DUMMY(filnam)
potnam='PATTENGILL-LABUDDE-BERNSTEIN AR-N2'
lammin(1)=2
lammax(1)=2
mproj(1)=0
ntv(1)=1
ivcol(1,1)=0
ivrow(1,1)=0
return
end

subroutine pot (vv0, r)
!  -----------------------------------------------------------------------

!  subroutine to calculate the r-dependent coefficients in the
!  collision of a homonuclear diatomic with a structureless target
!  in units of hartree for energy and bohr for distance

!  on return:
!  vv0 contains the isotropic term (n=0) in the potential
!  the coefficients for each angular term in the coupling potential
!  [ vvl(i) for i = 1, nlam ] are returned in module mod_covvl
!  vvl(1) contains the anisotropic (n=2) term in the potential

!  -----------------------------------------------------------------------
use mod_covvl, only: vvl
use mod_conlam, only: nlam
use mod_parbas, only: lammin
implicit double precision (a-h,o-z)
real(8), intent(out) :: vv0
real(8), intent(in) :: r  ! intermolecular distance

!  pattengill ar-n2 potential [m.d. pattengill, r. a. labudde, r.b. bernstein,
!  and c.f. curtiss, j. chem. phys. 55, 5517 (1971)]
!  updated to full double precision 5/5/97 by mha

data eps, a6, a12, r0 / 83.05d0, 0.13d0, 0.5d0, 3.929d0/


!     convert r0 to bohr

rr0 = r0 / 0.52917715d0

!     convert eps to hartree

eps_au = eps / 219474.6d0

rat = rr0 / r
r2 = rat * rat
r3 = r2 * rat
r6 = r3 * r3
r12 = r6 * r6
vv0   = eps_au * (r12 - 2.d0* r6)
vvl(1) = eps_au * (a12 * r12 - 2.d0 * a6 * r6)
if (nlam .ge. 2) vvl(2)=vvl(1)*0.2d0
if (nlam .ge. 3) vvl(3)=vvl(1)*0.1d0
if (nlam .ge. 4) vvl(4)=vvl(1)*0.d0
if (nlam .ge. 5) vvl(5)=vvl(1)*0.025d0
if (nlam .ge. 6) vvl(6)=vvl(1)*0.02d0
return
end
