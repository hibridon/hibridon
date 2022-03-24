module mod_ba_dummy
    use mod_basis, only: ab_basis


    type, extends(ab_basis) :: hiba_dummy

    contains
    procedure, nopass :: ba
    procedure, nopass :: sy
    procedure, nopass :: sav
    procedure, nopass :: ptr

    end type hiba_dummy


contains
subroutine ba(j, l, is, jhold, ehold, ishold, nlevel, nlevop, sc1, sc2, sc3, sc4, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, ihomo, nu, numin, jlpar, n, nmax, ntop)
         integer :: j(:), l(:), is(:), nlevel, nlevop, jtot, nu, numin, jlpar, n, nmax, ntop
         logical :: ihomo, flaghf, csflag, clist, flagsu, bastst
         real(8) :: jhold(5), ehold(5), sc1(5), sc2(5), sc3(5), sc4(5), ishold(5)
         real(8) :: rcut
end subroutine ba

subroutine sy(irpot, readpt, iread)
         integer :: iread, irpot
         logical :: readpt
end subroutine sy

subroutine sav(readpt)
         logical :: readpt
end subroutine sav

subroutine ptr(fname,readp)
         character*(*) fname
         logical :: readpt
end subroutine ptr


end module mod_ba_dummy