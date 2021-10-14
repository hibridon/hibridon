
subroutine fassert(filename,linenum)
	character(*), intent(in) :: filename
	integer,      intent(in) :: linenum
	write(6,'(a,a,a,i4.4)') "assertion failed in file ", filename,":",linenum
	stop 1
end subroutine fassert

! subroutine FortranAssert(prepost,expression,filename,linenum)
! 	character(*), intent(in) :: prepost,expression, filename
! 	integer,      intent(in) :: linenum
! 	write(6,'(a,a,a,a,a,a,i4.4)') prepost," (", expression, ") failed in file ", filename," line ",linenum
! 	stop
! end subroutine FortranAssert
