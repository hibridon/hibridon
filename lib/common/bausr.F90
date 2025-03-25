#include "hiutil.inc.F90"
!comdeck bausr
! --------------------------------------------------------------------
subroutine bausr (bqs, jhold, ehold, ishold, nlevel, nlevop, &
                  sc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu, &
                  csflag, clist, bastst, ihomo, nu, numin, jlpar, &
                  n, nmax, ntop, v2)
use mod_ancou, only: ancou_type, ancouma_type
use mod_hitypes, only: bqs_type
implicit none
type(bqs_type), intent(out) :: bqs
integer, intent(out), dimension(:) :: jhold
real(8), intent(out), dimension(:) :: ehold
integer, intent(out), dimension(:) :: ishold
integer, intent(out) :: nlevel
integer, intent(out) :: nlevop
real(8), intent(out), dimension(:) :: sc1
real(8), intent(out), dimension(:) :: sc2
real(8), intent(out), dimension(:) :: sc3
real(8), intent(out), dimension(:) :: sc4
real(8), intent(in) :: rcut
integer, intent(in) :: jtot
logical, intent(in) :: flaghf
logical, intent(in) :: flagsu
logical, intent(in) :: csflag
logical, intent(in) :: clist
logical, intent(in) :: bastst
logical, intent(in) :: ihomo
integer, intent(in) :: nu
integer, intent(in) :: numin
integer, intent(in) :: jlpar
integer, intent(out) :: n
integer, intent(in) :: nmax
integer, intent(out) :: ntop
type(ancou_type), intent(out), allocatable, target :: v2

UNUSED(bqs)
call bqs%init(0)  ! to silence warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.
UNUSED(jhold)
jhold(1) = 0  ! to silence warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.
UNUSED(ehold)
ehold(1) = 0.0  ! to silence warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.
UNUSED(ishold)
ishold(1) = 0  ! to silence warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.
UNUSED(nlevel)
nlevel = 0  ! to silence warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.
UNUSED(nlevop)
nlevop = 0  ! to silence warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.
UNUSED(sc1)
sc1(1) = 0.0  ! to silence warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.
UNUSED(sc2)
sc2(1) = 0.0  ! to silence warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.
UNUSED(sc3)
sc3(1) = 0.0  ! to silence warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.
UNUSED(sc4)
sc4(1) = 0.0  ! to silence warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.
UNUSED(rcut)
UNUSED(jtot)
UNUSED(flaghf)
UNUSED(flagsu)
UNUSED(csflag)
UNUSED(clist)
UNUSED(bastst)
UNUSED(ihomo)
UNUSED(nu)
UNUSED(numin)
UNUSED(jlpar)
UNUSED(n)
n = 0  ! to silence warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.
UNUSED(nmax)
UNUSED(ntop)
ntop = 0  ! to silence warning #6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.
UNUSED(v2)

stop "this bausr subroutine is a default one. It is not supposed to be called. It's here just to allow linking in case the user doesn't provide one."

return
end
! --------------------------------------------------------------------
