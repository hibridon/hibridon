!************************************************************************
!                                                                       *
!                           module basis                                *
!                                                                       *
!************************************************************************
!  Contains definition of the abstract class basis                      *
!************************************************************************
module mod_basis

   type, abstract :: ab_basis
      integer            :: id          ! Id number of the basis (1-30 & and 99 currently)
      character(len=200) :: description ! Short description of the collision type
   contains
      ! All the following procedures all called in hisystem.F90, for all the basis
      procedure(compute_angular_coupling), deferred, nopass :: compute_angular_coupling
      procedure(read_sys_dep_params), deferred, nopass :: read_sys_dep_params
      procedure(save), deferred, nopass :: save
      procedure(read_pot), deferred, nopass :: read_pot

      procedure, public :: is_twomol => basis_is_twomol
      procedure, public :: uses_j12 => basis_uses_j12
      ! procedure(init_spin_interface), deferred, nopass :: basis_init_spin
      procedure, public :: init_spin => basis_init_spin
   end type ab_basis

   interface
      ! Interface for ba subroutine
      ! subroutine to determine angular coupling potential
      subroutine compute_angular_coupling(j, l, is, jhold, ehold, ishold, nlevel, nlevop, sc1, sc2, sc3, sc4, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, ihomo, nu, numin, jlpar, twomol, n, nmax, ntop)
         integer :: j(:), l(:), is(:), nlevel, nlevop, jtot, nu, numin, jlpar, n, nmax, ntop
         integer :: jhold(5), ishold(5)
         logical :: ihomo, flaghf, csflag, clist, flagsu, bastst, twomol
         real(8) :: ehold(5), sc1(5), sc2(5), sc3(5), sc4(5)
         real(8) :: rcut
      end subroutine compute_angular_coupling
      ! Interface for read_sys_dep_params subroutine
      ! subroutine to read in system dependent parameters
      subroutine read_sys_dep_params(irpot, readpt, iread)
         integer :: iread, irpot
         logical :: readpt
      end subroutine read_sys_dep_params
      ! Interface for sav entry of sy subroutine
      subroutine save(readpt)
         logical :: readpt
      end subroutine save
      ! Interface for ptr entry of sy subroutine
      subroutine read_pot(fname,readpt)
         character*(*) fname
         logical :: readpt
      end subroutine read_pot

      ! subroutine init_spin_interface(this, flaghf, spnj2, spnj12, spntot, spin)
      !    import
      !    class(ab_basis), intent(in) :: this
      !    logical, intent(in) :: flaghf
      !    real(8), intent(out) :: spnj2
      !    real(8), intent(out) :: spnj12
      !    real(8), intent(out) :: spntot
      !    real(8), intent(inout) :: spin
      ! end subroutine init_spin_interface

   end interface
contains
   subroutine init_vib_qnumbers
#include "common/parbas.F90"
      ! set default for vibrational quantum numbers to zero for each term
      do it=1,maxtrm
         ivrow(1,it)=0
         ivcol(1,it)=0
         ntv(it)=1
      end do
   end subroutine init_vib_qnumbers

   logical function basis_is_twomol(this)
      implicit none
      class(ab_basis), intent(in) :: this
      integer :: ibasty
      ibasty = this%id
      basis_is_twomol = is_twomol(ibasty)
   end function basis_is_twomol

   logical function basis_uses_j12(this)
      implicit none
      class(ab_basis), intent(in) :: this
      integer :: ibasty
      ibasty = this%id
      basis_uses_j12 = is_j12(ibasty)
   end function basis_uses_j12

   subroutine basis_init_spin(this, flaghf, spnj2, spnj12, spntot, spin)
      implicit none
      class(ab_basis), intent(in) :: this
      logical, intent(in) :: flaghf
      real(8), intent(out) :: spnj2
      real(8), intent(out) :: spnj12
      real(8), intent(out) :: spntot
      real(8), intent(inout) :: spin
      integer :: ibasty
      ibasty = this%id
      spnj2 = 0.d0
      spnj12 = 0.d0
      spntot = 0.d0
      if (flaghf) then
        spnj12 = 0.5d0
        spntot = 0.5d0
        if (ibasty.eq.12 .or. ibasty.eq.15) then
          spin = 0.d0
          spnj2 = 0.5d0
        endif
      endif
      if (ibasty.eq.23) then
        spnj2 = 0.5d0
        spnj12 = 0.5d0
        spntot = 0.5d0
      endif
   end subroutine basis_init_spin


!  ------------------------------------------------------------------
   logical function is_twomol(ibasty)
   !     ------------------------------------------------------------------
   !
   !     checks if a basis is for molecule-molecule collision (j=10j1+j2)
   !
   !     written by q. ma
   !     current revision:  24-jul-2019 (p.dagdigian)
   !     ------------------------------------------------------------------
   implicit none
   integer, intent(in) :: ibasty
   if ((ibasty .eq. 9) .or. (ibasty .eq. 20) .or. (ibasty .eq. 21) &
   .or. (ibasty .eq. 25) .or. (ibasty .eq. 26) &
   .or. (ibasty .eq. 28) .or. (ibasty .eq. 30) &
   .or. (ibasty .eq. 100)) &
   then
   is_twomol = .true.
   else
   is_twomol = .false.
   end if
   return
   end function is_twomol
   !     ------------------------------------------------------------------
   logical function is_j12(ibasty)
   !     ------------------------------------------------------------------
   !
   !     checks if j12 is used in a basis
   !
   !     written by q. ma
   !     current revision:  17-oct-2018 (p.dagdigian)
   !     ------------------------------------------------------------------
   implicit none
   integer, intent(in) :: ibasty
   if (is_twomol(ibasty) .or. (ibasty .eq. 12) &
   .or. (ibasty .eq. 13) .or. (ibasty .eq. 15) &
   .or. (ibasty .eq. 23)) then
   is_j12 = .true.
   else
   is_j12 = .false.
   end if
   return
   end function is_j12
   
end module mod_basis


    