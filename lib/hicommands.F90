#include "assert.h"
#include "command.inc.F90"

module mod_hicommands
use mod_command, only: command_type, command_mgr_type
implicit none
  class(command_mgr_type), pointer :: command_mgr  ! singleton

  ! run
  type, extends(command_type) :: run_command_type
  contains
    procedure :: execute => run_execute
  end type run_command_type

  ! showpot
  type, extends(command_type) :: showpot_command_type
  contains
    procedure :: execute => showpot_execute
  end type showpot_command_type

  ! show
  type, extends(command_type) :: show_command_type
  contains
    procedure :: execute => show_execute
  end type show_command_type

  ! prsbr
  type, extends(command_type) :: prsbr_command_type
  contains
    procedure :: execute => prsbr_execute
  end type prsbr_command_type

contains

  subroutine update_nu_params()
    !  numin and numax should be 0 if cc calculation, if not, then set them
    !  equal to zero
    use mod_par, only: lpar, ipar
    use ipar_enum
    use lpar_enum

    common /coselb/ ibasty
    integer ibasty

    if (.not. lpar(LPAR_CSFLAG)) then
      lpar(LPAR_NUCROS)=.false.
      if (ipar(IPAR_NUMAX) .ne. 0) then
        write (6, 105)
        105     format ('  CC calculation, numax set to zero')
        ipar(IPAR_NUMAX) = 0
      end if
      ! NB this is disabled currently for 2P atom + homonuclear
      if (ipar(IPAR_NUMIN) .ne. 0.and.ibasty.ne.12) then
        write (6, 106)
        106     format ('  CC calculation, numin set to zero')
        ipar(IPAR_NUMIN) = 0
      end if
    end if
  end subroutine update_nu_params

  subroutine print_main_params(out_unit)
    use mod_hibasis, only: basknd
    use mod_cosys, only: scod
    use mod_cosysl, only: islcod, lspar
    use mod_cosysi, only: nscode, isicod, ispar
    use mod_cosysr, only: isrcod, rspar
    use mod_hinput_state, only: lindx
    use mod_si_params, only: lcode
    use mod_par, only: lpar
    use lpar_enum
    implicit none
    integer, intent(in) :: out_unit

#include "common/parbas.F90"

    common /coselb/ ibasty
    integer ibasty

    common /cotwo/ numj,nj1j2(5)
    integer :: numj
    integer :: nj1j2

    ! fcod = Flags CODes : stores the name of system independent parameters of type logical
    common /cofcod/ fcod
    character*8 :: fcod(lcode)

    integer :: j
    integer :: length

    if (ibasty .lt. 99) then
      length = index(basknd(ibasty),' ') - 1
      if (length .eq. -1) length=9
#if defined(HIB_UNIX) || defined(HIB_CRAY) || defined(HIB_MAC)
      write(out_unit,710) basknd(ibasty)(1:length)//' system parameters:', &
                 (scod(j),ispar(j),j = 1,isicod)
#endif
    else
      write(out_unit,710) 'user defined system parameters:', &
                   (scod(j),ispar(j),j = 1,isicod)
    endif
    if(isrcod.gt.0) &
      write(out_unit,720) (scod(isicod+j),rspar(j),j = 1,isrcod)
    if(islcod.gt.0) &
      write(out_unit,735) (scod(isicod+isrcod+j),lspar(j),j = 1,islcod)
    if (.not. lpar(LPAR_TWOMOL) ) then
      write(out_unit,736) 'LAMMIN: ',(lammin(j),j=1,ispar(1))
      write(out_unit,736) 'LAMMAX: ',(lammax(j),j=1,ispar(1))
      write(out_unit,736) 'MPROJ:  ',(mproj(j),j=1,ispar(1))
    else if (lpar(LPAR_TWOMOL)) then
      write (out_unit, 738)'J1/J2: ',(nj1j2(j)/10,mod(nj1j2(j),10), &
                        j=1,numj)
    end if
    write(out_unit,730) 'Flags:',(fcod(j),lpar(lindx(j)),j = 1,lcode)
    710 format(5x,'*** ',(a)/(4(1x,a7,'=',i4,7x)))
    720 format(4(1x,a7,'=',1pg11.4))
    730 format(5x,'*** ',(a)/(6(1x,a6,'=',l2,3x)))
    735 format(3(1x,a7,'=',l2,9x))
    736   format(1x,(a),10i4,/,9x,10i4)
    738   format (1x,(a),1x,20(2i1,'  ') )

  end subroutine print_main_params

  subroutine run_execute(this, statements, bofargs, next_statement, post_action)
    use mod_par, only: ipar, lpar, rpar
    use mod_si_params, only: iicode, ircode, icode, lcode, set_param_names
    use ipar_enum
    use lpar_enum
    use mod_coiout, only: niout, indout
    use mod_cosout, only: nnout, jout
    use mod_coener, only: energ
    use mod_command, only: k_post_action_interpret_next_statement, k_post_action_read_new_line
    use mod_hibrid2, only: enord
    use mod_hinput_state, only: lindx, irpot, irinp
    use mod_command, only: k_post_action_write_cr_and_exit, k_post_action_exit_hinput, k_post_action_exit_hibridon
    use rpar_enum
    !
    ! label:execute_run
    !
    ! start execution,run
    use mod_command, only: k_post_action_read_new_line
    implicit none
    class(run_command_type) :: this
    character(len=K_MAX_USER_LINE_LENGTH), intent(in) :: statements
    integer, intent(in) :: bofargs
    integer, intent(out) :: next_statement
    integer, intent(out) :: post_action

#include "common/parpot.F90"

    common /cofile/ input, output, jobnam, savfil
    character*40 :: input
    character*40 :: output
    character*40 :: jobnam
    character*40 :: savfil

    common /coselb/ ibasty
    integer ibasty

    ! pcod = Parameters CODes : stores the name of system independent parameters of type integer and real
    common /copcod/ pcod
    character*8 :: pcod(icode)

    common /cosavi/ iipar, ixpar(iicode)
    integer :: iipar
    integer :: ixpar

    common /cosavr/ irpar, junks, rxpar(ircode)
    integer :: irpar
    integer :: junks
    real(8) :: rxpar

    integer :: i, j
    integer :: nerg

    call update_nu_params()
    call set_param_names(lpar(LPAR_BOUNDC),pcod,icode)
    if (lpar(LPAR_CSFLAG).and.ipar(IPAR_NUD).ne.1) lpar(LPAR_NUCROS)=.true.
    nerg=ipar(IPAR_NERG)
    ! check to see if flags are ok if wavefunction desired or
    ! photodissociation calculation
    call genchk
    call enord(energ,nerg)
    do i = 1,ircode
      rxpar(i) = rpar(i)
    end do
    do i = 1,iicode
      ixpar(i) = ipar(i)
    end do
    if(irinp.eq.0) then
      write(6,505)
      505 format (/,' ** SAVE DEFAULT VARIABLES OR SPECIFY INPUT', &
              ' FILE WITH INP = filename')
      if(lpar(LPAR_BATCH)) then
        post_action = k_post_action_exit_hibridon
        return
      else
        post_action = k_post_action_read_new_line
        return
      end if
    end if
    if(rpar(RPAR_SCAT_TOLAI).eq.0) then  ! graffy: todo : shouldn't it be RPAR_XMU instead of RPAR_SCAT_TOLAI here ?
      write(6,507)
      507   format(/,' ** SPECIFY COLLISION REDUCED MASS WITH XMU = mass')
      post_action = k_post_action_read_new_line
      return
    end if
    if(irpot.ne.0.or..not.lpar(LPAR_READPT)) then
    ! open output file
    ! first make sure it is lower case
      call lower(output)
      call upper(output(1:1))
      call openf(9,output,'sf',0)
    !     write input data to file outpt
      write (9, 508) label
      508   format (1x,a)
      write (9, 508) potnam
      write (9, 240)
      240   format(1h ,30('='))
      write(9,710) 'Parameters:',(pcod(j),ipar(j),j = 1,iicode)
      write(9,720) (pcod(iicode+j),rpar(j),j = 1,ircode)
      call print_main_params(out_unit=9)
      710 format(5x,'*** ',(a)/(4(1x,a7,'=',i4,7x)))
      720 format(4(1x,a7,'=',1pg11.4))

      736   format(1x,(a),10i4,/,9x,10i4)
      737   format(1x,(a),i2,(a),20(i3,1x))
      739   format (1x,(a),i2,(a),20(2i1,'  ') )
      call enord(energ,ipar(IPAR_NERG))
      if(ipar(IPAR_NERG).gt.0) write(9,740) (energ(j),j = 1,ipar(IPAR_NERG))
      740 format(1x,'** Energies:',(t15,5f15.6))
      if(nnout.ne.0) then
        if (.not.lpar(LPAR_TWOMOL) ) then
          write (9,737) 'NOUT: ',nnout, &
               '; JOUT:',(jout(j), j=1,iabs(nnout) )
        else
          write (9,739) 'NOUT: ',nnout, &
             '; J1/J2-OUT: ', &
            (jout(j)/10, mod (jout(j),10), j = 1,iabs(nnout))
        end if
      end if
      if(niout.ne.0) write (9,736) 'INDOUT: ',(indout(j), j=1,niout)
      write (9, 240)
      post_action = k_post_action_exit_hinput
      return
    else
      write(6,510)
    510   format(' Potential not yet defined!')
      post_action = k_post_action_read_new_line
      return
    end if

    post_action = k_post_action_write_cr_and_exit
  end subroutine run_execute


  subroutine showpot_execute(this, statements, bofargs, next_statement, post_action)
    use mod_command, only: k_post_action_read_new_line
    class(showpot_command_type) :: this
    character(len=K_MAX_USER_LINE_LENGTH), intent(in) :: statements
    integer, intent(in) :: bofargs
    integer, intent(out) :: next_statement
    integer, intent(out) :: post_action
    write(6,*) "************************************************************"
    write(6,*) "Entering the DRIVER subroutine of the potential"
    write(6,*) "Press Ctrl+D to go back to Hibridon's console"
    write(6,*) "************************************************************"
    call driver
    post_action = k_post_action_read_new_line
  end subroutine

  subroutine show_execute(this, statements, bofargs, next_statement, post_action)
    use mod_par, only: ipar, lpar, rpar
    use mod_hinput_state, only: batch
    use mod_si_params, only: iicode, ircode, icode, lcode, set_param_names
    use ipar_enum
    use lpar_enum
    use mod_coiout, only: niout, indout
    use mod_conlam, only: nlammx
    use mod_codim, only: nmax => mmax
    use mod_cosout, only: nnout, jout
    use mod_coener, only: energ
    use mod_command, only: k_post_action_interpret_next_statement, k_post_action_read_new_line
    use mod_hibrid2, only: enord


    ! show all parameters and flags
    ! show
    class(show_command_type) :: this
    character(len=K_MAX_USER_LINE_LENGTH), intent(in) :: statements
    integer, intent(in) :: bofargs
    integer, intent(out) :: next_statement
    integer, intent(out) :: post_action

#include "common/parbas.F90"
#include "common/parpot.F90"


    common /cofile/ input, output, jobnam, savfil
    character*40 :: input
    character*40 :: output
    character*40 :: jobnam
    character*40 :: savfil

    common /coselb/ ibasty
    integer ibasty

    ! pcod = Parameters CODes : stores the name of system independent parameters of type integer and real
    common /copcod/ pcod
    character*8 :: pcod(icode)

    ! fcod = Flags CODes : stores the name of system independent parameters of type logical
    common /cofcod/ fcod
    character*8 :: fcod(lcode)

    common /cotwo/ numj,nj1j2(5)
    integer :: numj
    integer :: nj1j2

    integer :: j
    integer :: l, lc, length, leninp, lenjob, lenout
    logical :: jtrunc
    character*40 :: code
    character(len=K_MAX_USER_LINE_LENGTH) :: answer

    l = bofargs
    next_statement = l
    call set_param_names(lpar(LPAR_BOUNDC),pcod,icode)
    call get_token(statements,l,code,lc)

    if (.not.lpar(LPAR_BOUNDC)) then
      write(6,710) &
       'Parameters (scattering):',(pcod(j),ipar(j),j = 1,iicode)
    else
      write(6,710) &
       'Parameters (bound-state):',(pcod(j),ipar(j),j = 1,iicode)
    endif
    write(6,720)   (pcod(iicode+j),rpar(j),j = 1,ircode-1)
    write(6,1720)  pcod(iicode+ircode), rpar(ircode)
    if(nnout.ne.0) then
      if (.not.lpar(LPAR_TWOMOL) ) then
        write (6,737) 'NOUT: ',nnout, &
               '; JOUT:',(jout(j), j=1,iabs(nnout) )
      else
        write (6,739) 'NOUT: ',nnout, &
             '; J1/J2-OUT: ', &
            (jout(j)/10, mod (jout(j),10), j = 1,iabs(nnout))
      end if
    end if
    if(niout.ne.0) write (6,701) 'INDOUT:',(indout(j), j=1,niout)
    10 format((a))
    701   format(1x,(a),10i5,/,8x,10i5,/,5x,10i5)
    call print_main_params(out_unit=6)
    write(6,731)  nmax, nlammx
    if (ipar(IPAR_LSCREEN) .le. 24 .and. .not. batch) then
      write (6, 703)
    703   format (6x,'enter <return> to continue,', &
                 ' <q> for prompt')
      read (5, 10) answer
      if (answer(1:1) .eq. 'q' .or. answer(1:1) .eq. 'q') then  ! fixme: both 'or' conditions are identical
        post_action = k_post_action_read_new_line
        return
      end if
    end if
    call enord(energ,ipar(IPAR_NERG))
    if(ipar(IPAR_NERG).gt.0) write(6,740) (energ(j),j = 1,ipar(IPAR_NERG))
    710 format(5x,'*** ',(a)/(4(1x,a7,'=',i4,7x)))
    720 format(4(1x,a7,'=',1pg11.4))
    1720 format(1x,a7,'=',f10.5)
    731 format(1x,'** Maximum Channels: ', i4, '; ', &
      'Maximum Anisotropic Terms: ',i5)
    737   format(1x,(a),i2,(a),20(i3,1x))
    739   format (1x,(a),i2,(a),20(2i1,'  ') )
    740 format(1x,'** Energies:',(t15,5f15.6))
    lenout = index(output,' ')-1
    lenjob = index(jobnam,' ')-1
    jtrunc = .false.
    if (lenjob.gt.8) then
       jobnam=jobnam(1:8)
       jtrunc = .true.
    endif
    leninp=index(input,' ')-1
    write (6, 751) label
    751 format (1x,'** Label:      ',(a))
    write (6, 752) potnam
    752 format (1x,'** Pot name:      ',(a))
    if (.not. jtrunc) then
       write (6, 753) input(1:leninp), &
         output(1:lenout), jobnam(1:lenjob)
    753 format(1x,'** Input File:  ',(a),/ &
           1x,'** Output file: ',(a),/,1x,'** Jobname:     ',(a))
    else
       write (6, 754) input(1:leninp), &
         output(1:lenout), jobnam(1:lenjob)
    754 format(1x,'** Input File:  ',(a),/ &
           1x,'** Output file: ',(a),/,1x,'** Jobname:     ',(a), &
           ' (** TRUNCATED TO 8 CHARACTERS **)')
    endif
    next_statement = l
    post_action = k_post_action_interpret_next_statement
  end subroutine show_execute


  subroutine prsbr_execute(this, statements, bofargs, next_statement, post_action)
    ! pressure broadening cross sections - added by p. dagdigian
    use mod_hiutil, only: assignment_parse
    use mod_command, only: k_post_action_read_new_line
    class(prsbr_command_type) :: this
    character(len=K_MAX_USER_LINE_LENGTH), intent(in) :: statements
    integer, intent(in) :: bofargs
    integer, intent(out) :: next_statement
    integer, intent(out) :: post_action
    integer :: l
    character(len=40) :: fnam1, fnam2
    character*8 empty_var_list(0)
    character*40 :: code
    integer :: i, j, lc
    integer, parameter :: k_num_args = 12
    real(8) :: a(k_num_args)

    common /cofile/ input, output, jobnam, savfil
    character*40 :: input
    character*40 :: output
    character*40 :: jobnam
    character*40 :: savfil


    l = bofargs

    call get_token(statements,l,fnam1,lc)
    if(fnam1 .eq. ' ') fnam1 = jobnam
    call lower(fnam1)
    call upper(fnam1(1:1))
    ! get iener for 1st smt file
    a(1) = 0.d0
    if(l .eq. 0) goto 3205
    call get_token(statements,l,code,lc)
    call assignment_parse(code(1:lc),empty_var_list,j,a(1))
    3205 call get_token(statements,l,fnam2,lc)
    if(fnam2 .eq. ' ') fnam2 = jobnam
    call lower(fnam2)
    call upper(fnam2(1:1))
    ! get iener for 2nd smt file
    a(2) = 0.d0
    if(l .eq. 0) goto 3210
    call get_token(statements,l,code,lc)
    call assignment_parse(code(1:lc),empty_var_list,j,a(2))
    ! get k, j1, in1, j2, in2, diag, j1p, in1p, j2p, in2p
    3210 do 3220 i = 3, 12
      a(i) = 0.d0
      if(l .eq. 0) goto 3220
      call get_token(statements,l,code,lc)
      call assignment_parse(code(1:lc),empty_var_list,j,a(i))
    3220 continue
    call prsbr(fnam1,fnam2,a)
    post_action = k_post_action_read_new_line
  end subroutine prsbr_execute



  subroutine init()
    class(command_type), allocatable :: com
    if (.not. associated(command_mgr)) then
      allocate(command_mgr)
      command_mgr%num_commands = 0
    end if

    com = run_command_type()
    call command_mgr%register_command('RUN', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    com = showpot_command_type()
    call command_mgr%register_command('SHOWPOT', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    com = show_command_type()
    call command_mgr%register_command('SHOW', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    com = prsbr_command_type()
    call command_mgr%register_command('PRSBR', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    ASSERT(associated(command_mgr))
  end subroutine init

end module mod_hicommands
