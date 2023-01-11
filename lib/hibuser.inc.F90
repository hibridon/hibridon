interface
    subroutine bausr(bqs, jhold, ehold, ishold, nlevel, nlevop, &
        sc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu, &
        csflag, clist, bastst, ihomo, nu, numin, jlpar, &
        n, nmax, ntop, v2)
        !use mod_ancou, only: ancou_type
        import :: ancou_type
        import :: bqs_type

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
    end subroutine

end interface
