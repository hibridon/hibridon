module mod_hipot
  ! provides interfaces to user-provided subroutines related to potential
  interface

    subroutine pot(vv0, r)
      real(8), intent(out) :: vv0
      real(8), intent(in) :: r  ! intermolecular distance
    end subroutine

    subroutine loapot(iunit, filnam)
      integer, intent(in) :: iunit  ! if a data file is used, this subroutine is expected to use this unit to open it in read mode (not used here)
      character*(*), intent(in) :: filnam  ! if a data file is used, the file name of the data file
    end subroutine

    subroutine ground(wf, r, nch, nphoto, mxphot)
      real(8), intent(out) :: wf(nch*nphoto) ! array of dimension nch*nphoto, containing, on return,
      ! ground state wavefunction in each of nch components
      ! nphoto is number of difference ground state wavefunctions
      real(8), intent(out) :: r  ! value of separation coordinate
      integer, intent(in) :: nch  ! total number of channels (row dimension of q)
      integer, intent(in) :: nphoto ! number of different wavefunctions calculated
      ! column index of q vector
      integer, intent(in) :: mxphot  ! maximum size of q vector (mxphot .ge. nch*nphoto)
    end subroutine

    subroutine wfintern(wf, yymin, nch, nphoto, nny, ifull)
      real(8), intent(out) :: wf(nch*nphoto) ! array of dimension nch*nphoto, containing, on return,
      ! ground state wavefunction in each of nch components
      ! nphoto is number of difference ground state wavefunctions
      real(8), intent(in) :: yymin
      integer, intent(in) :: nch  ! total number of channels (row dimension of q)
      integer, intent(in) :: nphoto ! number of different wavefunctions calculated
      ! column index of q vector
      integer, intent(in) :: nny
      logical, intent(in) :: ifull
    end subroutine

  end interface
end module mod_hipot
