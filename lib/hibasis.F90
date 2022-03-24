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
      procedure, nopass :: ba
      procedure, nopass :: sy
      procedure, nopass :: sav
      procedure, nopass :: ptr
   end type ab_basis

   interface
      ! Interface for ba subroutine
      subroutine ba(j, l, is, jhold, ehold, ishold, nlevel, nlevop, sc1, sc2, sc3, sc4, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, ihomo, nu, numin, jlpar, n, nmax, ntop)
         integer :: j(:), l(:), is(:), nlevel, nlevop, jtot, nu, numin, jlpar, n, nmax, ntop
         logical :: ihomo, flaghf, csflag, clist, flagsu, bastst
         real(8) :: jhold(5), ehold(5), sc1(5), sc2(5), sc3(5), sc4(5), ishold(5)
         real(8) :: rcut
      end subroutine ba
      ! Interface for sy subroutine
      subroutine sy(irpot, readpt, iread)
         integer :: iread, irpot
         logical :: readpt
      end subroutine sy
      ! Interface for sav entry of sy subroutine
      subroutine sav(readpt)
         logical :: readpt
      end subroutine sav
      ! Interface for ptr entry of sy subroutine
      subroutine ptr(fname,readp)
         character*(*) fname
         logical :: readpt
      end subroutine ptr
   end interface
end module mod_basis
    