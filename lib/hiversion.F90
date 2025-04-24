module mod_version
save
contains
subroutine execute_command(command, command_exit_code, error_message)
  character(*), intent(in) :: command
  integer, intent(out) :: command_exit_code
  character(*), intent(out) :: error_message
  integer :: cmd_stat

  call execute_command_line(command, exitstat=command_exit_code, cmdstat=cmd_stat, cmdmsg=error_message)
  if ( cmd_stat /= 0 ) then
      write(6, *) "the command ", command, " failed"
      write(6, *) "error message :", error_message
      write(6, *) "exit code :", command_exit_code
      command_exit_code = 127
  end if

end subroutine

subroutine get_tmp_dir(tmp_dir)
  character(len=256), intent(out) :: tmp_dir
  integer :: get_env_status

  call get_environment_variable("TMPDIR", tmp_dir, status=get_env_status)
  if(get_env_status /= 0) then
    tmp_dir = '.'
    ! stop 'failed to get the value of TMPDIR environment variable'
  end if

end subroutine

subroutine get_command_stdout(command, stdout, default_stdout)
  character(*), intent(in) :: command
  character(*), intent(out) :: stdout
  character(*), intent(in) :: default_stdout
  character(len=:), allocatable :: redirected_command
  integer :: exit_stat
  integer :: iostat
  integer, parameter :: tmp_file_unit = 666
  character(len=256) :: cmd_msg
  character(len=256) :: tmp_dir
  character(len=256) :: tmp_file_path
  character(*), parameter :: redirect_operator = ' > '

  call get_tmp_dir(tmp_dir)
  tmp_file_path = trim(tmp_dir) // '/hw.info'

  allocate(character(len=len(command) + len(redirect_operator) + len(tmp_file_path)) :: redirected_command)
  redirected_command = trim(command) // redirect_operator // tmp_file_path
  call execute_command(redirected_command, &
     command_exit_code=exit_stat, error_message=cmd_msg)
  if ( exit_stat == 0 ) then
      open(unit=tmp_file_unit, file=tmp_file_path, status="old")
      read(tmp_file_unit,'(a)', iostat=iostat) stdout
      close(tmp_file_unit)
  else
      stdout = default_stdout
  end if
  deallocate(redirected_command)
end subroutine

subroutine version(out_unit)
  implicit none
  integer, intent(in) :: out_unit
  integer :: i
  character(len=100) :: txt
  character(len=256) :: tmp_dir

  call get_tmp_dir(tmp_dir)

  write(out_unit,10) _BUILD_VERS_
10   format (/' ---------------------------------------------------', &
  '-----------------------',/ &
   ,9x, &
   '  HIBRIDON SCATTERING CODE V ',f4.2,// &
   '      AUTHORS: M. ALEXANDER, D. MANOLOPOULOS,', &
   ' H.-J. WERNER, B. FOLLMEG,'/ &
   '               P. DAGDIGIAN',/ &
   ' CONTRIBUTORS: D. LEMOINE, P. VOHRALIK,', &
   ' G. COREY, R. JOHNSON, T. ORLIKOWSKI,',/ &
   '          A. BERNING, A. DEGLI-ESPOSTI,', &
   ' C. RIST, B. POUILLY, J. KLOS, Q. MA,',/, &
   '          G. VAN DER SANDEN, M. YANG, F. DE WEERD', &
   ', S. GREGURICK, F. LIQUE',/ &
   ' ---------------------------------------------------', &
   '-----------------------')

 write (out_unit,15)
15  format(11x,'BUILD CONFIGURATION:')
  write(out_unit,'(a18,a)') 'CMAKE BUILD DATE: ', _BUILD_DATE_
  write(out_unit,'(a18,a)') 'CMAKE BUILD TYPE: ', _BUILD_TYPE_
  write(out_unit,'(a19,a)') 'CODE GIT REVISION: ', _GIT_REVISION_
  write(out_unit,'(a5,a,a1,a,a1,a)')  'SYS: ', _BUILD_SYS1_ , ' ', &
                                        _BUILD_SYS2_, ' ', &
                                        _BUILD_PROC_
  write(out_unit,'(a10,a)') 'COMPILER: ', _BUILD_COMP_
  write(out_unit,'(a17,a)') 'COMPILE OPTIONS: ', _BUILD_FLAGS_
  write (out_unit,20)
20  format (' ---------------------------------------------------', &
  '-----------------------')

write (out_unit,30)
30 format(11x,' CURRENT HARDWARE CONFIGURATION:')
  write(out_unit,*)
  if (_BUILD_SYS1_ == "Darwin") then
      call execute_command_line('system_profiler SPHardwareDataType | grep -E "Model|Processor|Cores|Cache|Hyper|Memory" > ' // trim(tmp_dir) // '/hw.info')
      open(unit=666, file=(trim(tmp_dir) // "/hw.info"), status="old")
      do 
          read(666,'(a)', iostat=i) txt
          if(i.ne.0) exit
          write(out_unit,'(a)') trim(txt)
      enddo 
      close(666)
elseif (_BUILD_SYS1_ == "Linux") then

  call get_command_stdout('lscpu | grep -m 1 "Model name"', stdout=txt, default_stdout='Model name: unknown')
  write(out_unit,'(a)') trim(txt)    

  call get_command_stdout('lscpu | grep -m 1 "CPU MHz"', stdout=txt, default_stdout='CPU MHz: unknown')
  write(out_unit,'(a)') trim(txt)    

  call get_command_stdout('lscpu | grep -m 1 "Socket(s)"', stdout=txt, default_stdout='Socket(s): unknown')
  write(out_unit,'(a)') trim(txt)    

  call get_command_stdout('lscpu | grep -m 1 "Core(s)"', stdout=txt, default_stdout='Core(s): unknown')
  write(out_unit,'(a)') trim(txt)    

  call get_command_stdout('lscpu | grep -m 1 "L2"', stdout=txt, default_stdout='L2: unknown')
  write(out_unit,'(a)') trim(txt)    

  call get_command_stdout('lscpu | grep -m 1 "L3"', stdout=txt, default_stdout='L3: unknown')
  write(out_unit,'(a)') trim(txt)    

  call get_command_stdout('cat /proc/cpuinfo | grep -m 1 -o "ht"', stdout=txt, default_stdout='')
  if(trim(txt)=='ht') then
      write(out_unit,'(a35)') "Hyper-Threading Technology: Enabled"
  end if

  call get_command_stdout('cat /proc/meminfo | grep "MemTotal"', stdout=txt, default_stdout='')
  write(out_unit,'(a)') trim(txt) 
endif
write(out_unit,31)
31 format( &
 ' ---------------------------------------------------', &
  '-----------------------')

end subroutine version

subroutine acknow(iunit,ipos)
  implicit none
  integer :: iunit
  logical :: ipos
  if (ipos) write (iunit, 10)
10 format(/' --------------------------------------------------------', &
  '----------------------------------------------------------', &
  '-------------', &
  /,' All publications resulting from use of the integrators', &
  ' included in the Hibridon code must include', &
   ' the following reference:', &
 //,' D. E. Manolopoulos, J. Chem. Phys. 85, 6425 (1986);', &
  ' M. H. Alexander and D. E. Manolopoulos, J. Chem. Phys.', &
  ' 80, 2044 (1987).', &
  //,' All publications resulting from use of the', &
  ' Hibridon package must include', &
  ' the following reference:', &
 //,' HIBRIDON is a package of programs for the time-independent', &
  ' quantum treatment', &
  ' of inelastic collisions and photodissociation', &
  /,' written by', &
  ' M. H. Alexander, D. E.  Manolopoulos, H.-J. Werner, and', &
 ' B. Follmeg,', &
  ' with contributions by', &
  ' P. F. Vohralik, D. Lemoine,', &
  /,' G. Corey, B. Johnson, T. Orlikowski, A. Berning,', &
  ' A. Degli-Esposti, C. Rist, P. Dagdigian, B. Pouilly,', &
  ' G. van der Sanden, M. Yang, F. de Weerd, S. Gregurick,', &
  ' J. Klos, and F. Lique', &
 /,' --------------------------------------------------------', &
  '-------------', &
  '----------------------------------------------------------')
if (.not.ipos) write (iunit, 20)
20 format(/' --------------------------------------------------------', &
  '---------------------', &
  /,' All publications resulting from use of the integrators', &
  ' included in the',/,' Hibridon code must include', &
   ' the following reference:', &
 //,' D. E. Manolopoulos, J. Chem. Phys. 85, 6425 (1986);', &
  ' M. H. Alexander and D. E.',/, &
    ' Manolopoulos, J. Chem. Phys.', &
  ' 80, 2044 (1987).', &
  //,' All publications involving the determination of', &
     'photodissociation cross sections',/,' must also include', &
     ' the following reference:', &
  //,' M. H. Alexander, Comput. Phys. Commun, 75, 87 (1993).', &
  //,' All publications investigating flux redistribution must', &
  ' include the following references:', &
  //, ' M. H. Alexander, J. Chem. Phys. 95, 8931 (1991); D. E.', &
      ' Manolopoulos and',/,' M. H. Alexander, J. Chem. Phys.', &
      '  97, 2527 (1992).', &
  //,' All publications resulting from use of the', &
  ' Hibridon package must include', &
  /,' the following reference:', &
 //,' HIBRIDON is a package of programs for the time-independent', &
  ' quantum treatment', &
  /,' of inelastic collisions and photodissociation', &
  ' written by', &
  ' M. H. Alexander,',/,' D. E.  Manolopoulos, H.-J. Werner, and', &
 ' B. Follmeg,', &
  ' with contributions by', &
  /,' P. F. Vohralik, D. Lemoine,', &
  ' G. Corey, B. Johnson, T. Orlikowski, A. Berning,', &
  /,' A. Degli-Esposti, C. Rist, P. Dagdigian, B. Pouilly,', &
  ' G. van der Sanden, M. Yang, F. de Weerd, S. Gregurick,', &
  ' J. Klos, and F. Lique', &
 /,' --------------------------------------------------------', &
  '---------------------')
end subroutine acknow
end module mod_version
