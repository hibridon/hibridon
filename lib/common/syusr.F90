!comdeck syusr
subroutine syusr (irpot, readpt, iread)
!  dummy syusr subroutine
logical readpt
character*(*)fname
if(.not.readpt.or.iread.eq.0) then
  call loapot(1,' ')
  return
endif
entry ptrusr (fname)
entry savusr (readpt)
entry chkusr
return
end
