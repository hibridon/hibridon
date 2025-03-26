module mod_hipot
  ! provides interfaces to user-provided subroutines related to potential
  interface

    subroutine pot(vv0, r)
      real(8), intent(out) :: vv0
      real(8), intent(in) :: r
    end subroutine

    subroutine loapot(iunit, filnam)
      integer, intent(in) :: iunit  ! if a data file is used, this subroutine is expected to use this unit to open it in read mode (not used here)
      character*(*), intent(in) :: filnam  ! if a data file is used, the file name of the data file
    end subroutine

  end interface
end module mod_hipot
