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
end module mod_basis
    