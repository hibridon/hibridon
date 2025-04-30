#include "assert.h"
#include "command.inc.F90"
#if defined(HIB_UNIX_IBM) || defined(HIB_UNIX_AIX)
@proc ss noopt
#endif
#if defined(HIB_UNIX_HP)
!$hp$optimize off
#endif


module mod_candidates

  ! dimension of codex, ihold, lhold, should be equal to largest number
  ! of identical strings of 1:nnn characters in names of all variables
  ! (probably 'p' is the most recurring string:  12 times in
  !  pcod, fcod, and bcod)
  integer, parameter :: k_max_candidates = 15

  type candidates_type
  ! this class represents a list of candidate codexes that match the shortcut codex the user has inputed.
  
  ! For example, if the user inputs the statement 'PR=T', then the candidates would be any codex that starts with 'PR'
  ! - system independent flags:
  !   - 'PRAIRY'
  !   - 'PRLOGD'
  !   - 'PRPART'
  !   - 'PRSMAT'
  !   - 'PRT2'
  !   - 'PRXSEC'
  ! - commands:
  !   - 'PRINTC'
  !   - 'PRINTS'
  !   - 'PRSBR'
  private
    character(len=8) codex(k_max_candidates)
    integer :: ihold(k_max_candidates)
    integer :: lhold(k_max_candidates)
    integer :: num_candidates = 0
  contains
    procedure, public :: add_candidate => candidates_add_candidate
    procedure, public :: get_num_candidates => candidates_get_num_candidates
    procedure, public :: get_codex => candidates_get_codex
    procedure, public :: get_codex_index => candidates_get_codex_index
    procedure, public :: get_bost => candidates_get_bost
    procedure, public :: empty => candidates_empty
  end type candidates_type

contains

!
! candidates_type implementation
!

subroutine candidates_add_candidate(this, codex, codex_index, bost)
  class(candidates_type), intent(inout) :: this
  character(len=8), intent(in) :: codex
  integer, intent(in) :: codex_index
  integer, intent(in) :: bost  ! index of line that points to the beginnning of the statement (eg 'JTOT=42')

  this%num_candidates = this%num_candidates + 1
  if (this%num_candidates > k_max_candidates) then
    stop 'error : k_max_candidates is too small'
  end if
  this%lhold(this%num_candidates) = bost
  this%ihold(this%num_candidates) = codex_index
  this%codex(this%num_candidates) = codex
end subroutine candidates_add_candidate

function candidates_get_num_candidates(this) result(num_candidates)
  class(candidates_type), intent(in) :: this
  integer :: num_candidates
  num_candidates = this%num_candidates
end function candidates_get_num_candidates

function candidates_get_codex(this, candidate_index) result(codex)
  use mod_assert, only: fassert
  class(candidates_type), intent(in) :: this
  integer, intent(in) :: candidate_index
  character(len=8) :: codex
  ASSERT(candidate_index <= this%num_candidates)
  codex = this%codex(candidate_index)
end function candidates_get_codex

function candidates_get_codex_index(this, candidate_index) result(codex_index)
  use mod_assert, only: fassert
  class(candidates_type), intent(in) :: this
  integer, intent(in) :: candidate_index
  integer :: codex_index
  ASSERT(candidate_index <= this%num_candidates)
  codex_index = this%ihold(candidate_index)
end function candidates_get_codex_index

function candidates_get_bost(this, candidate_index) result(bost)
  use mod_assert, only: fassert
  class(candidates_type), intent(in) :: this
  integer, intent(in) :: candidate_index
  integer :: bost
  ASSERT(candidate_index <= this%num_candidates)
  bost = this%lhold(candidate_index)
end function candidates_get_bost

subroutine candidates_empty(this)
  class(candidates_type), intent(inout) :: this
  this%num_candidates = 0
end subroutine candidates_empty

end module mod_candidates

module mod_hinput
  use mod_assert, only: fassert
  implicit none

  enum, bind( C )
  enumerator :: &
    k_keyword_execute_command     =  1, &   !   40 label:execute_command(i)
    k_keyword_set_si_ir_param     =  2, &   !  100 label:set_si_ir_param(line, l)
    k_keyword_set_si_l_param      =  3, &   !  200 label:set_si_l_param(line, l)
    k_keyword_set_sd_param        =  4, &   ! 1400 label:set_sd_param(line, l)
    k_keyword_set_ibasty          =  5, &   !   50 label:set_ibasty(line,l)
    k_keyword_execute_command_mgr_command     =  6   !   45 label:execute_command_mgr_command(i)
  end enum

  integer, parameter :: ncode = 26  !  ncode is the number of bcod's
  character(len=8), parameter :: bcod(ncode) = [ &  ! bcod stores hibridon's commands
    'DEBROGLI', &
    'DIFFER  ', &
    'DIFCRS  ', &
    'ENERGY  ', &
    'EXIT    ', &
    'HELP    ', &
    'JOUT    ', &
    'MINPOT  ', &
    'MRCRS   ', &
    'NNOUT   ', &
    'OPTIMIZE', &
    'POT     ', &
    'PRINTC  ', &
    'PRINTS  ', &
    'PSI     ', &
    'QUIT    ', &
    'SAVE    ', &
    'TENXSC  ', &
    'TESTPOT ', &
    'TURN    ', &
    'INDOUT  ', &
    'PARTC   ', &
    'FLUX    ', &
    'J1J2    ', &
    'EADIAB  ', &
    'SYSCONF ']

  character(len=8), parameter :: bascod(1) = ['BASISTYP']

contains
subroutine hinput(first_time)
!  subroutine to redefine system independent input parameters for
!  hibridon code
!  author:  h.-j. werner
!  revision: 13-nov-1996 by mby to include bound states
!  revised further:  24-feb-1998 by mha
!  maxbas extended to 28:  30-jul-2018
!  jobname limited to 8 characters:   04-oct-2001, format corrected 2-dec-2007
!  added command sysconf 24-dec-2007
!  added command to compute cross sections for mixed singlet-triplet states
!    e.g. CH3 X3B1 - a1A1 collisions involving perturbed levels - p. dagdigian
!  altered input list for hypxsc:  5-jan-2011 by p. dagdigian
!  added bastp1 basis subroutine (for symmetric top without inversion doubling:
!     p. dagdigian
!  added calculation of transport cross sections - p. dagdigian
!  added 2sigma-2pi (no perts.) basis routine - p. dagdigian
!  last line "quit/exit" now optional in a .com file - q. ma
!  allow calculation of differential cross sections for ibasty=9 - p.dagdigian
!  eliminate special call for ibasty=4 for printc command - p. dagdigian
!
!  revision: 15-mar-2012 by p. dagdigian
!  revision: 20-apr-2012 by q. ma (new eadiab command)
!  revision: 8-oct-2012 by q. ma (allow calculation for user-defined ibasty 100,
!     whose channels have j12's)
!  revision: 27-apr-2013 by p. dagdigian (add prsbr [pressure broadening
!     cross sections] command)
!  revision: 10-jul-2013 by q. ma (new 2pi--1sigma basis)
!  revision: 18-jul-2013 by q. ma (new symtop--1sigma basis)
!  revision: 18-sep-2014 by p. dagdigian (new 3P atom + 2S atom basis)
!  revision: 21-jul-2015 by p. dagdigian (new spherical top-atom basis)
!  revision: 10-may-2107 by p. dagdigian (hypxsc for mol-mol collisions implemented)
!  revision:  1-jun-2017 by p. dagdigian (new 1sig + 1sig basis)
!  revision:  8-jun-2017 by p. dagdigian (new 2sig + 1sig basis)
!  revision: 19-sep-2017 by p. dagdigian (new C2V asym top basis)
!  revision: 30-jul-2018 by p. dagdigian (new 3sig + 1sig basis)
!  revision: 16-jan-2019 by p. dagdigian (new chiral astop + atom basis)
!  revision: 20-jun-2019 by p. dagdigian (new C2V asym - diat molecule basis)
!
!  current revision: 16-jan-2019 by q. dagdigian
! ---------------------------------------------------------------------
use mod_com, only: com_file, com
use mod_cosout, only: nnout, jout
use mod_coiout, only: niout, indout
use mod_codim, only: nmax => mmax
use mod_coamat, only: scmat => toto ! scmat(1)
use mod_coener, only: energ, max_en
use mod_cosys, only: scod
use mod_cosysi, only: nscode, isicod, ispar
use mod_cosysl, only: islcod, lspar
use mod_cosysr, only: isrcod, rspar
use mod_version, only : version
use mod_hibrid5, only : readpc
use mod_difcrs, only: difcrs
use mod_hibasis, only: is_twomol
use mod_hibrid2, only: enord, prsg
use mod_hibrid3, only: potmin
use mod_hiutil, only: assignment_parse
use fcod_enum
use lpar_enum
use ipar_enum
use rpar_enum
use mod_par, only: lpar, ipar, rpar, fcod, pcod
use mod_candidates, only: candidates_type
use mod_hicommands, only: command_init => init, command_mgr, command_type, update_nu_params
use mod_hinput_state, only: batch
use mod_si_params, only: iicode, ircode, icode, lcode, set_param_names
use mod_hinput_state, only: lindx, irpot, irinp
use mod_command, only: k_post_action_interpret_next_statement, k_post_action_read_new_line, k_post_action_exit_hibridon, k_post_action_exit_hinput, k_post_action_write_cr_and_exit
use mod_selb, only: ibasty
use mod_file, only: input, jobnam
use mod_sav, only: iipar, ixpar, irpar, rxpar
use mod_tensor, only: tenopa, mrcrs
use mod_hitestptn, only: testptn
use mod_two, only: numj, nj1j2
use mod_opti, only: optifl
use mod_hiutil, only: get_token, lower, upper, lenstr, vaxhlp, sys_conf
use mod_hibrid1, only: difs, turn
use mod_hibrid4, only: psi, eadiab1, sprint
use mod_hypxsc, only: hypxsc
use mod_hiiolib1, only: openf, gendat, savdat, genchk
use mod_hisystem, only: baschk, sysdat, syssav, ptread
use mod_histmix, only: stmix
implicit none
logical, intent(inout) :: first_time
character(len=K_MAX_USER_LINE_LENGTH) line
character(len=40) :: fnam1
character(len=40) :: fnam2
character*40 :: code
character*8 empty_var_list(0)
integer nerg
logical logp, opti
real(8) :: a(15)  ! real arguments
integer :: ia(10)  ! integer arguments of commands

! these were members of cokeyl common block
integer :: nncode  
integer :: llcode
integer :: ijcode


integer :: ipr, istep, inam, i, ienerg, iflux, ii, im, imx, inew, iprint, iskip, itx, ityp, izero
integer :: j, jm, jmx, jtot2x, l, l1, l2, lc, lcc, ld, len
integer :: nde
real(8) :: optacm, r, thrs, val, waveve, xmu
real(8) :: a1, acc, acclas, optval, optacc, accmx, delt_e, e, e1
type(candidates_type) :: candidates
! class(command_type), pointer :: command
integer :: post_action
integer :: next_statement  ! index of the next statement to be interpreted
save ipr, opti, a, a1, acc, acclas, optval, optacc, istep, inam, &
     code, lc, jtot2x

if (first_time) then
  call command_init()
end if
ASSERT(associated(command_mgr))

nerg = 0
lindx(FCOD_AIRYFL) = LPAR_AIRYFL
lindx(FCOD_BASTST) = LPAR_BASTST
lindx(FCOD_BATCH) = LPAR_BATCH
lindx(FCOD_CHLIST) = LPAR_CHLIST
lindx(FCOD_CSFLAG) = LPAR_CSFLAG
lindx(FCOD_FLAGHF) = LPAR_FLAGHF
lindx(FCOD_FLAGSU) = LPAR_FLAGSU
lindx(FCOD_IHOMO) = LPAR_IHOMO
lindx(FCOD_IPOS) = LPAR_IPOS
lindx(FCOD_LOGDFL) = LPAR_LOGDFL
lindx(FCOD_NOPRIN) = LPAR_NOPRIN
lindx(FCOD_NUCROS) = LPAR_NUCROS
lindx(FCOD_PHOTOF) = LPAR_PHOTOF
lindx(FCOD_PRAIRY) = LPAR_PRAIRY
lindx(FCOD_PRLOGD) = LPAR_PRLOGD
lindx(FCOD_PRPART) = LPAR_PRPART
lindx(FCOD_PRSMAT) = LPAR_PRSMAT
lindx(FCOD_PRT2) = LPAR_PRT2
lindx(FCOD_PRXSEC) = LPAR_PRXSEC
lindx(FCOD_READPT) = LPAR_READPT
lindx(FCOD_RSFLAG) = LPAR_RSFLAG
lindx(FCOD_T2TEST) = LPAR_T2TEST
lindx(FCOD_TWOMOL) = LPAR_TWOMOL
lindx(FCOD_WAVEFL) = LPAR_WAVEFL
lindx(FCOD_WRPART) = LPAR_WRPART
lindx(FCOD_WRSMAT) = LPAR_WRSMAT
lindx(FCOD_WRXSEC) = LPAR_WRXSEC
lindx(FCOD_BOUNDC) = LPAR_BOUNDC

! addresses for commands
! debrogli: 1800
! differ: 1500
! difcrs: 2000
! energy: 300
! exit: 600
! help: 75
! jout: 400
! minpot: 1700
! mrcrs: 2500
! nnout: 2400
! optimimize: 2100
! pot: 1000
! printc: 2600
! prints: 1900
! psi: 2800
! quit: 600
! save: 1300
! tenxsc: 2300
! testpot: 1200
! turn: 1600
! indout: 430
! partc: 2650
! flux: 2800
! eadiab: 2850
! j1j2:  460
! sysconf:  2900
! nb after changing the following list, check that all the variables "incode"
! that follow after address 900 are changed accordingly
!
iipar=iicode
irpar=ircode
nncode=ncode
llcode=lcode
ijcode=icode
! Open command file given by user
if(com) open(unit=1312, status='old', file=trim(com_file))

!   define system dependent parameter codes
if(first_time) then
   islcod=0
   isrcod=0
   isicod=0
   izero=0
   call sysdat(irpot, lpar(LPAR_READPT), izero)
   first_time = .false.
   call version(6)
!  in this next statement the $ sign implies no line feed
!  replace this with an equivalent formatting character if your system
!  doesn't accept this extension
   goto 1  ! label:read_new_line
else
  do 3 i = 1, ircode
3   rpar(i) = rxpar(i)
  do 4 i = 1, iicode
4   ipar(i) = ixpar(i)
  if(opti) goto 2160
end if
! label:read_new_line
1 if (.not. lpar(LPAR_BATCH) .and. .not. batch) write (6, 2)
optifl = .false.
!  in this next statement the $ sign implies no line feed
!  replace this with an equivalent formatting character if your system
!  doesn't accept this extension
#if defined(HIB_UNIX) || defined(HIB_MAC)
2 format(' Hibridon> ',$)
#endif
#if defined(HIB_CRAY)
2  format(' Hibridon> ')
#endif
call set_param_names(lpar(LPAR_BOUNDC),pcod,icode)

! read the next command
if(com) then  
  read(1312, 10, end=599) line  ! label:write_cr_and_exit in case of error
else
  read(5, 10, end=599) line  ! label:write_cr_and_exit in case of error
endif
10 format((a))
if(line .eq. ' ') goto 1  ! label:read_new_line
if (line(1:1) .eq. '?') then
    code='help '//line(2:)
    call vaxhlp(code)
!         call helppr(line)
    goto 1  ! label:read_new_line
else if (line (1:4).eq.'help' .or. line (1:4).eq.'HELP') then
    call vaxhlp(line)
!         line = '?intro'
!         call helppr(line)
    goto 1  ! label:read_new_line
else if (line(1:3) .eq.'BAT' .or. line(1:3) .eq. 'bat' .or. &
         line(1:4) .eq.' BAT' .or. line(1:4) .eq.' bat') then
    lpar(LPAR_BATCH)=.true.
    batch = .true.
    goto 1  ! label:read_new_line
end if
call upper(line)
l1 = 1
!
! label:interpret_next_statement(line, l1)
!
! interpret the statement starting at index l1 of line
! examples of statements:
!  - 'VMAX(1)=4'
!  - 'RUN'
! warning! line can contain multiple statements separated with comma or semicolon (eg 'VMAX(1)=4;VMAX(2)=4')
15 continue
if(l1 .eq. 0) goto 1  ! label:read_new_line
ASSERT(l1 >= 0)  ! graffy: I suspect that l1 can never be negative
! read the next token, which is expected to be either a command or a parameter
l = -iabs(l1) ! consider '=' as a token delimiter
ASSERT(l <= 0)
call get_token(line, l, code, lc)
ASSERT(l >= 0)
if(lc .eq. 0) goto 1  ! label:read_new_line
call candidates%empty()
! search in commands
do i = 1,ncode
  len = lenstr(bcod(i))
  if(bcod(i)(1:lc) .eq. code(1:lc)) then
    if (lc .eq. len) goto 40  ! label:execute_command(i)
    call candidates%add_candidate(codex=bcod(i), codex_index=i, bost=l)
    iskip = k_keyword_execute_command
  end if
end do
! search in command_mgr
do i = 1, command_mgr%num_commands
  len = lenstr(command_mgr%commands(i)%codex)
  if(command_mgr%commands(i)%codex(1:lc) .eq. code(1:lc)) then
    if (lc .eq. len) goto 45  ! label:execute_command_mgr_command(i)
    call candidates%add_candidate(codex=bcod(i), codex_index=i, bost=l)
    iskip = k_keyword_execute_command_mgr_command
  end if
end do

ASSERT(l1 >= 0) ! graffy: I suspect that l1 can never be negative
l = iabs(l1)  ! reset l at the start of a statement because the different keyword handling entry points will expect it to be pointing to the beginning of the statement
ASSERT(l >= 0)  ! make sure l is positive (even if one day we remove the suspected unneeded iabs above)
! search in system independent parameters of type integer and real
do i = 1,icode
  len = index(pcod(i),' ') - 1
  if(pcod(i)(1:lc) .eq. code(1:lc)) then
    if (lc .eq. len) goto 100  ! label:set_si_ir_param(line, l)
    call candidates%add_candidate(codex=pcod(i), codex_index=i, bost=l)

    iskip = k_keyword_set_si_ir_param
  end if
end do
! search in system independent parameters of type logical
do i = 1,lcode
len=lenstr(fcod(i))
  if(fcod(i)(1:lc) .eq. code(1:lc)) then
    if (lc .eq. len) goto 200  ! label:set_si_l_param(line, l)
    call candidates%add_candidate(codex=fcod(i), codex_index=i, bost=l)
    iskip = k_keyword_set_si_l_param
  end if
end do
! search in legacy base parameters (system dependent parameters)
do i = 1,nscode
  len = index(scod(i),' ') - 1
  if(scod(i)(1:lc) .eq. code(1:lc)) then
    if (lc .eq. len) goto 1400  ! label:set_sd_param(line, l)

    call candidates%add_candidate(codex=scod(i), codex_index=i, bost=l)
    iskip = k_keyword_set_sd_param
  end if
end do
! search in bascod parameters (only contains BASISTYP at the moment)
len = 8
if(bascod(1)(1:lc) .eq. code(1:lc)) then
  if (lc .eq. len) goto 50  ! label:set_ibasty(line,l)
  call candidates%add_candidate(codex=bascod(1), codex_index=1, bost=l)
  iskip = k_keyword_set_ibasty
end if
if (candidates%get_num_candidates() == 0) then
  ! the input string matched none of the parameters
  write(6, 27) code(1:lc),(bcod(j),j = 1,ncode)
27   format( &
    /' *** invalid keyword "',(a),'"; valid request keys are:'// &
   (1x,6(a8,5x)))
  write(6,31) bascod, (pcod(j),j = 1,icode)
31   format (/(1x,6(a8,5x)))
  write(6,31) (fcod(j),j = 1,lcode)
  write(6,31) (scod(j),j = 1,nscode)
  goto 1  ! label:read_new_line
else if (candidates%get_num_candidates() > 1) then
  ! more than one parameter matched the input string : we can't decide which one to choose
  write (6, 28) code(1:lc), (candidates%get_codex(i), i = 1, candidates%get_num_candidates())
28   format (' *** ambiguity between input string ',(a), &
           ' and request keys:',/,7(a10) )
  goto 1  ! label:read_new_line
else
  ! exactly one parameter matched the input string, even if the input string was shorter than the parameter
  ! the input string is then considered a valid (ie non ambiguous) shortcut for the parameter
  ! set the parameter accordingly
  i = candidates%get_codex_index(1)
  l = candidates%get_bost(1)
!  goto (40, 100, 200, 1400, 50, 1410, 1420, 1430), iskip
  goto (40, 100, 200, 1400, 50, 45), iskip
end if
!
! label:execute_command(i)
!
40 goto ( &
      1800,1500,2000,300,600, &
      75, &
      400,1700, &
      2500, &
      2400,2100,1000,2600, &
      1900,2800,600, &
      1300,2300, &
      1200,1600,430,2650,2800, &
      460,2850,2900),i
!
! label:execute_command_mgr_command(i)
!
45 continue
  next_statement = HIB_UNINITIALIZED_INTEGER_VALUE
  call command_mgr%commands(i)%item%execute(statements=line, bofargs=l, next_statement=next_statement, post_action=post_action)
#ifndef DISABLE_HIB_ASSERT
#if (HIB_UNINITIALIZED_INTEGER_VALUE != 0)
      ! make sure that next_statement has been initialized by the command's execute method
      ! note: this test can be removed when the warnings -Wunused-dummy-argument and -Wunused-variable trigger errors (when done, if a command forgets to set the output dummy argumen next_statement, then the compiler will detect a -Wunused-dummy-argument and trigger an error at compile time)
      if (next_statement == HIB_UNINITIALIZED_INTEGER_VALUE) then
        write(6,*) 'this command forgot to set next_statement :', command_mgr%commands(i)%codex
        ASSERT(next_statement /= HIB_UNINITIALIZED_INTEGER_VALUE)  
      end if
#endif
#endif
  
! label:on_execute_command_completion(post_action)
!
47 continue
  if (post_action == k_post_action_read_new_line) then
    goto 1  ! label:read_new_line
  else if (post_action == k_post_action_interpret_next_statement) then
    l1 = next_statement
    goto 15  ! label:interpret_next_statement(line, l1)
  else if (post_action == k_post_action_exit_hibridon) then
    call exit()
  else if (post_action == k_post_action_write_cr_and_exit) then
    goto 599  ! label:write_cr_and_exit
  else if (post_action == k_post_action_exit_hinput) then
    return
  else
    ASSERT( .false. )  ! unexpected value for post_action
  end if 
!
! label:set_ibasty(line,l)
! basis type and kind of calculation
! 
50 if(l.eq.0) goto 1  ! label:read_new_line
l1 = l
call get_token(line,l,code,lc)
call assignment_parse(code(1:lc),bascod,j,val)
if(j .eq. 0) goto 15  ! label:interpret_next_statement(line, l1)
if(j .lt. 0) goto 1  ! label:read_new_line
ibasty=int(val)
call baschk(ibasty)
! set twomolecule true
if (is_twomol(ibasty)) then
  lpar(LPAR_TWOMOL)=.true.
else
  lpar(LPAR_TWOMOL)=.false.
endif
call sysdat(irpot, lpar(LPAR_READPT), izero)
l1 = l
goto 15  ! label:interpret_next_statement(line, l1)
!  request help
!75    line = '?intro'
!      call helppr(line)
75 call vaxhlp(line)
goto 1  ! label:read_new_line
!
! label:set_si_ir_param(line, l)
!
! set system independent parameters (integer and real)
!  specify parameters in the form cod1=val1, cod2=val2, etc.
100 if(l .eq. 0) goto 1  ! label:read_new_line
ASSERT(l > 0)
call set_param_names(lpar(LPAR_BOUNDC),pcod,icode)
l1 = l
call get_token(line,l,code,lc)
call assignment_parse(code(1:lc),pcod,j,val)
if(j .eq. 0) goto 15  ! label:interpret_next_statement(line, l1)
if(j .lt. 0) goto 1  ! label:read_new_line
if (j .eq. 5) then
  if (ipar(j) .lt. val) write (6, 101)
101   format &
  (1x,'** NERG INCREASED; VERIFY ARRAY OF COLLISION ENERGIES')
end if
if(j.le.iicode) then
  ipar(j) = val
  if(ipar(IPAR_NUD).ne.1) lpar(LPAR_NUCROS)=.true.
else
  rpar(j-iicode) = val
end if
call enord(energ,ipar(IPAR_NERG))
call update_nu_params()
goto 100  ! label:set_si_ir_param(line, l)
!
! label:set_si_l_param(line, l)
!
! set system independent parameters (flags)
! specify flags in the form cod1=val1,cod2=val2, etc
! where val(i) must be either "t(rue)" or "f(alse)"
200 if(l .eq. 0) goto 1  ! label:read_new_line
l1 = l
call get_token(line,l,code,lc)
call assignment_parse(code(1:lc),fcod,j,val)
if(j .eq. 0) goto 15  ! label:interpret_next_statement(line, l1)
if(j .lt.0) goto 1  ! label:read_new_line
logp = .false.
if(val .eq. 1) logp = .true.
if (j .eq. 3) then
  write (6, 201)
201   format (' ** BATCH FLAG CAN NOT BE SET INTERACTIVELY!')
  lpar(LPAR_BATCH) = batch
  goto 1  ! label:read_new_line
end if
lpar(lindx(j)) = logp
goto 200  ! label:set_si_l_param(line, l)
! energies
! specify energies in the form
! energ=e1,e2,e3...
! terminate the string with a semicolon if other parameters will follow
! on the same card, e.g. energ=e1,e2,e3;jtot1=0,jtot2=2....
300 i = 0
310 if(l .eq. 0) goto 320
if(line(l-1:l-1) .eq. ';') goto 320
i = i+1
call get_token(line,l,code,lc)
call assignment_parse(code(1:lc),empty_var_list,j,energ(i))
goto 310
320 if (energ(1) .gt. 0) then
   if (i .ne. nerg) then
     if (i .le. max_en) then
       write (6, 321) i
321        format (' ** NERG HAS BEEN RESET TO',i3)
       nerg = i
     else
       write (6, 322) max_en
322        format (' ** NERG RESET TO ',i0,' (MAXIMUM VALUE)')
       nerg = max_en
     end if
   end if
elseif (energ(1) .lt. 0) then
   nerg=energ(4)+0.001d0
   if (nerg > max_en) then
      nerg = max_en
      write (6,322)
   endif
   e1=energ(2)
   delt_e=(energ(3)-e1)/(nerg-1)
   nde=delt_e*100d0
   delt_e=nde/100d0
   write (6,323) e1,nerg, delt_e,e1+(nerg-1)*delt_e
323    format(' ** GRID OF ENERGIES:  E-FIRST = ', f9.2, &
     '; NERG = ',i2, '; DELTA_E = ',f7.2, &
        '; E-LAST = ',f9.2)
   do ii=1,nerg
      energ(ii)=e1+(ii-1)*delt_e
   enddo
endif

call enord(energ,nerg)
ipar(IPAR_NERG) = nerg
l1 = l
goto 15  ! label:interpret_next_statement(line, l1)
! jout values
! specify jout values in the form
! jout,nnout,jout(1),...,jout(iabs(nnout))
! terminate the string with a semicolon if other parameters will follow
! on the same card, e.g. jout,-3,0,2,4;energ=e1,e2,e3;jtot1=0,jtot2=2....
400 i = 0
call get_token(line,l,code,lc)
call assignment_parse(code(1:lc),empty_var_list,j,val)
l1 = l
nnout=val
410 if(l .eq. 0) goto 420
if(line(l-1:l-1) .eq. ';') goto 420
call get_token(line,l,code,lc)
call assignment_parse(code(1:lc),empty_var_list,j,val)
i = i+1
jout(i) = val
goto 410
420 if(nnout.ge.0) nnout = i
if(nnout.lt.0) nnout = -i
l1 = l
goto 15  ! label:interpret_next_statement(line, l1)
! indout values
! specify indout values in the form
! indout,niout,indout(1),...,indout(niout)
! terminate the string with a semicolon if other parameters will follow
! on the same card, e.g. indout,2,1,-1;energ=e1,e2,e3;jtot1=0,jtot2=2....
430 i = 0
call get_token(line,l,code,lc)
call assignment_parse(code(1:lc),empty_var_list,j,val)
l1 = l
niout=val
if(niout.eq.0) goto 15  ! label:interpret_next_statement(line, l1)
440 if(l .eq. 0) goto 450
if(line(l-1:l-1) .eq. ';') goto 450
call get_token(line,l,code,lc)
call assignment_parse(code(1:lc),empty_var_list,j,val)
i = i+1
indout(i) = val
goto 440
450 if(niout.ge.0) niout = i
l1 = l
goto 15  ! label:interpret_next_statement(line, l1)
! j1j2 values
! specify j1j2 values in the form
! j1j2,numj,j1j2(1),...,j1j2(numj)
! terminate the string with a semicolon if other parameters will follow
! on the same card, e.g. j1j2,2,00,10;energ=e1,e2,e3;jtot1=0,jtot2=2....
460 if (.not.lpar(LPAR_TWOMOL)) then
  write (6, 465)
465   format(' ** NUMJ CAN ONLY BE DEFINED IF TWOMOL = .TRUE.')
  goto 15  ! label:interpret_next_statement(line, l1)
endif
i = 0
call get_token(line,l,code,lc)
call assignment_parse(code(1:lc),empty_var_list,j,val)
l1 = l
numj=val
if(niout.eq.0) goto 15  ! label:interpret_next_statement(line, l1)
470 if(l .eq. 0) goto 480
if(line(l-1:l-1) .eq. ';') goto 480
call get_token(line,l,code,lc)
call assignment_parse(code(1:lc),empty_var_list,j,val)
i = i+1
nj1j2(i) = val
goto 470
480 if(numj.ge.0) then
  numj = i
else
  write (6, 485)
485   format(' ** YOU MUST SPECIFY A VALUE OF NUMJ')
endif
l1 = l
goto 15  ! label:interpret_next_statement(line, l1)

!
! label:execute_run
!
! start execution,run
!  numin and numax should be 0 if cc calculation, if not, then set them
!  equal to zero
500 continue
call command_mgr%execute_command('RUN', post_action)
ASSERT(post_action /= k_post_action_interpret_next_statement)  ! make sure that this command doesn't output the position of the next statement in the current line, as it will be ignored
goto 47 ! label:on_execute_command_completion(post_action)

510   format(' Potential not yet defined!')

! no more commands
! label:write_cr_and_exit
599 write (6, *)
! exit
600 continue
call exit
! read parameters for potential
!     pot=potfile
1000 call get_token(line,l,code,lc)
call ptread(code(1:lc),lpar(LPAR_READPT))
l1 = l
goto 15  ! label:interpret_next_statement(line, l1)
! test potential
! testpot
! you will be prompted for r and theta. to exit, specify r=0
1200 call testptn(lpar(LPAR_IHOMO))
goto 1  ! label:read_new_line
! save input parameters
!     save=filename
!     if filename is not specified, the inputfile is overwritten
1300 inew=0
call set_param_names(lpar(LPAR_BOUNDC),pcod,icode)
if(l.ne.0) then
  call get_token(line,l,code,lc)
  if(lc .eq. 0) then
    code = input
  else
    inew=1
! convert input files to lower case
    call lower(code)
    call upper(code(1:1))
  endif
  
else
  code = input
end if
call savdat(inew,code)
call syssav()
close(8)
l1 = l
irinp = 1
goto 15  ! label:interpret_next_statement(line, l1)
!
! label:set_sd_param(line, l)

!
! redefine system dependent parameters
! specify in the same way as other parameters, e.g.
!     jmin=0,jmax=4,brot=2.2...
1400 continue
ASSERT(l >= 0)
if(l.eq.0) goto 1  ! label:read_new_line
l1=l
! read the next token, which is expected to be an assignment (eg 'VMAX(1)=4' if line(l:)=='VMAX(1)=4;VMAX(2)=4') 
call get_token(line,l,code,lc)
! at this point, l points to the beginning of the next statement
call assignment_parse(code(1:lc),scod,j,val)
if(j .eq. 0) then
  ! parameter not found
  goto 15  ! label:interpret_next_statement(line, l1)
end if
if(j .lt. 0) then
  ! malformed assignment (eg 'JTOT;')
  goto 1  ! label:read_new_line
end if
if(j.eq.1 .and. .not.lpar(LPAR_TWOMOL)) then
  write(6,'(1x,a,"CAN NOT BE MODIFIED")') scod(j)
  goto 1400  ! label:set_sd_param(line, l)
end if
if(j.le.isicod) then
  ispar(j) = val
else if(j.le.isrcod+isicod) then
  rspar(j-isicod) = val
else
  if(val .eq. 1) then
    lspar(j-isrcod-isicod) = .true.
  else
    lspar(j-isrcod-isicod) = .false.
  end if
end if
goto 1400  ! label:set_sd_param(line, l)
! !
! ! label:set_sd_i_param(l)
! ! 
! ! set system dependent integer parameter
! 1410 if(l.eq.0) goto 1  ! label:read_new_line
! l1=l
! ! l is expected to point at the beginning of a 'JTOT=42' assignment inside the line string
! call get_token(line,l,code,lc)
! ! code(1:lc) is expected to contain the param name eg 'JTOT'
! ! l is now expected to point on the character after the equal sign eg in 'JTOT=42'
! call assignment_parse(code(1:lc),scod,j,val)
! if(j .eq. 0) goto 15  ! label:interpret_next_statement(line, l1)
! if(j .lt. 0) goto 1  ! label:read_new_line
! if(j.eq.1 .and. .not.lpar(LPAR_TWOMOL)) then
!   write(6,'(1x,a,"CAN NOT BE MODIFIED")') scod(j)
!   goto 1410  ! label:set_sd_i_param()
! end if
! if(j.le.isicod) then
!   ispar(j) = val
! else if(j.le.isrcod+isicod) then
!   rspar(j-isicod) = val
! else
!   if(val .eq. 1) then
!     lspar(j-isrcod-isicod) = .true.
!   else
!     lspar(j-isrcod-isicod) = .false.
!   end if
! end if
! goto 1410  ! label:set_sd_i_param()
!.....differences of s-matrices
!  dif,jobfile1,jobfile2,iprint,ienerg,thrs
!  this compares s-matrices in the files jobfile1.smt and jobfile2.smt
!  if iprint.eq.0 largest average and absolute differences are printed
!  if iprint.eq.1 these values are given for all individual s-matrices
!  if iprint.ge.2 s-matrices are printed
!  thrs: threshold for neglect of small s-matrix elements in comparison
1500 call get_token(line,l,code,lc)
if(code.ne.' ') fnam1 = code
call lower(fnam1)
call upper(fnam1(1:1))
call get_token(line,l,code,lc)
if(code.ne.' ') fnam2 = code
call lower(fnam2)
call upper(fnam2(1:1))
if(fnam1 .eq. ' '.or.fnam2 .eq. ' ') goto 1  ! label:read_new_line
iprint = 0
ienerg = 1
thrs = 1.e-5
if(l.ne.0) then
  call get_token(line,l,code,lc)
  call assignment_parse(code(1:lc),empty_var_list,j,val)
  iprint = val
end if
if(l.ne.0) then
  call get_token(line,l,code,lc)
  call assignment_parse(code(1:lc),empty_var_list,j,val)
  ienerg = val
  ienerg = max0(1,ienerg)
end if
if(l.ne.0) then
  call get_token(line,l,code,lc)
  call assignment_parse(code(1:lc),empty_var_list,j,val)
  thrs = val
end if
call difs(fnam1,fnam2,ienerg,iprint,acc,accmx,thrs,imx,jmx,ityp)
code = '?'
lc = 1
if (thrs .lt. 0.) then
  if(ityp .eq. 1) then
    code = 'S real'
    lc = 6
  end if
  if(ityp .eq. 2) then
    code = 'S imaginary'
    lc = 11
  end if
else
  code = 'S modulus'
  lcc = 9
end if
if(iprint .eq. 0) write(6,1510) acc,accmx,imx,jmx,code(1:lc), &
  abs(thrs)
1510 format(' Average relative difference:',f10.2,'%'/ &
       ' Largest relative difference:',f10.2,'%'/ &
 ' (i = ',i2,'  j = ',i2,') element of ',(a)/ &
 ' Inspection threshold is ',1pg8.1)
goto 1  ! label:read_new_line
!.....determine turning point from isotropic potential
!     turn
1600 if(irpot .eq. 0) then
  write(6,510)  ! potentiel not yet defined
  goto 1  ! label:read_new_line
end if
e = 0
do 1605 i = 1,ipar(IPAR_NERG)
1605 e = max(e,energ(i))
if(e .eq. 0) then
  write(6,1610)
1610   format(' Total energy has not been given a value !')
  goto 1  ! label:read_new_line
end if
r = turn(e)
write(6,1620) r
1620 format(' Turning point for isotropic potential at R = ', &
         f5.2, ' bohr')
l1 = l
goto 15  ! label:interpret_next_statement(line, l1)
!  determine minimum of isotropic potential
!  minpot
1700 if(irpot .eq. 0) then
  write(6,510)
  goto 1  ! label:read_new_line
end if
r = potmin()
write(6,1710) r
1710 format(' Minimum of isotropic potential at r = ',f5.2, ' bohr')
l1 = l
goto 15  ! label:interpret_next_statement(line, l1)
!  calculate de broglie wavelength in bohr (defined as 2pi/k)
!  debrogli
1800 e = 0
do 1810 i = 1,ipar(IPAR_NERG)
1810 e = max(e,energ(i))
if(e .eq. 0) then
  write(6,1610)
  goto 1  ! label:read_new_line
end if
xmu = rpar(RPAR_XMU)
if(xmu .eq. 0) then
  write(6,1820)
1820   format(' Collision reduced mass has not been given a value !')
  goto 1  ! label:read_new_line
end if
r = 48.75/sqrt(xmu*e)
write(6,1830) r,0.2*r
1830 format('   de Broglie wavelength  = ',f6.3,' bohr / 5  = ',f6.3)
waveve = 6.283/r
write (6, 1831) waveve, 1.8897*waveve
1831 format ('   wavevector =', g12.5,' Bohr^-1 = ', &
   g12.5,' Angstroms^-1')
l1 = l
goto 15  ! label:interpret_next_statement(line, l1)
!
! label:execute_command_prints(line, l)
!
!.....print s-matrices:
!  prints,jobfile,j1,j2,jd,jlp,ienerg
!  first matrix printed is for jtot=j1
!  last matrix printed is for jtot=j2
!  increment of jtot is jd
!  default values are:
!  j1 = jtot1, j2=jtot1, jd=jtotd, jlp=jlpar,ienerg=1
!  jlp: parity; zero means both values
1900 continue
ASSERT(l >= 0)  ! verify that '=' is not treated as a delimiter in get_token
! read the job file name
call get_token(line,l,fnam1,lc)
if(fnam1 .eq. ' ') fnam1 = jobnam
call lower(fnam1)
call upper(fnam1(1:1))
! read the remaining arguments
do i = 1,5
  ia(i) = 0
  if(l .eq. 0) exit
  call get_token(line,l,code,lc)
  call assignment_parse(code(1:lc),empty_var_list,j,a(i))
  ia(i)=a(i)
end do
if(ia(1).eq.0) ia(1)=ipar(IPAR_JTOT1)
if(ia(2).eq.0) ia(2)=ipar(IPAR_JTOT2)
if(ia(3).eq.0) ia(3)=ipar(IPAR_JTOTD)
if(ia(4).eq.0) ia(4)=ipar(IPAR_JLPAR)
call lower(fnam1)
call upper(fnam1(1:1))
call sprint(fnam1,ia)
goto 1  ! label:read_new_line
!.....differential cross sections:
!  diffc,jobfile,j1,in1,j2,in2,ang1,ang2,dang,ienerg,jtotend
!
!  differential cross sections computed for atom-molecule
!  collisions or symmetric top-linear molecule collisions
2000 if (.not. lpar(LPAR_TWOMOL) .or. is_twomol(ibasty)) then
  call get_token(line,l,fnam1,lc)
  if(fnam1 .eq. ' ') fnam1 = jobnam
  call lower(fnam1)
  call upper(fnam1(1:1))
  do 2010 i = 1,15
     a(i) = 0.d0
     if(l .eq. 0) goto 2010
     call get_token(line,l,code,lc)
     call assignment_parse(code(1:lc),empty_var_list,j,a(i))
2010   continue
  write (6,*) 'hinput : lpar = ' 
  call difcrs(fnam1,a,lpar(LPAR_FLAGHF))
else
  write (6, 2011)
2011   format(' Sorry, differential cross sections not yet', &
         /,'  implemented for most molecule-molecule ', &
         'collisions')
end if
goto 1  ! label:read_new_line
!.....optimize
!  opt,code,start,end,fak,add,accav,accmx,thrs
!  code is name of variable to be optimized
!  start: start value
!  end:   end value
!  the value is updated in each step according to
!  val(n)=val(n-1)*fak+add
!  if(fak.eq.0) fak=1.0
!  optimization is stopped if the relative difference in the
!  s-matrices of subsequent runs are smaller than accav
!  (average error) and accmx (maximum error)
!  optimization is performed for jtot=jtot1 only
!  thrs:  threshold for check for s-matrix elements. all s-matrix
!  elements which are smaller than thrs are not compared
!  comparison is of s modulus if thrs .ge. 0, otherwise of s matrix
2100 if (ipar(IPAR_NERG) .gt. 1) then
! no optimization if more than one energy requested
  write (6, 2101)
2101   format(' ** NERG SET EQUAL TO 1 FOR OPTIMIZATION')
  ipar(IPAR_NERG)=1
endif
call get_token(line,l,code,lc)
do 2110 ipr = iicode+1,icode
2110 if(code(1:lc) .eq. pcod(ipr)(1:lc)) goto 2120
write(6,2115) code(1:lc)
2115 format(' Invalid parameter: "',(a),'"')
goto 1  ! label:read_new_line
2120 do 2130 i = 1,7
a(i) = 0.d0
if(l .eq. 0) goto 2130
call get_token(line,l,code,ld)
call assignment_parse(code(1:ld),empty_var_list,j,a(i))
2130 continue
if(a(1) .eq. 0.or.a(2) .eq. 0) then
  write(6,2140)
2140   format(' Initial or final value of parameter to optimize', &
         ' has not been defined')
  goto 1  ! label:read_new_line
end if
if (.not.lpar(LPAR_WRSMAT)) then
  write (6,2141)
2141   format (' Flag WRSMAT set to .TRUE. for optimization')
  lpar(LPAR_WRSMAT) = .true.
end if
if (nnout .lt. 0) then
  write (6, 2142)
2142   format (' NNOUT set positive for optimization')
  nnout = iabs (nnout)
else if (nnout .eq. 0) then
  write (6, 2143)
2143   format (' NNOUT=0; optimization not possible; reset NNOUT')
  goto 1  ! label:read_new_line
end if
if (.not.lpar(LPAR_NOPRIN)) then
   write (6,2144)
2144    format(' Flag NOPRIN set to .TRUE. for optimization')
   lpar(LPAR_NOPRIN)=.true.
endif
!
if(a(3) .eq. 0.and.a(4) .eq. 0.and.a(2).gt.a(1)) a(3) = 2.0d0
if(a(3) .eq. 0.and.a(4) .eq. 0.and.a(2).lt.a(1)) a(3) = 0.5d0
if(a(3) .eq. 0) a(3) = 1.0d0
if(a(5) .eq. 0) a(5) = 1.0d0
if(a(6) .eq. 0) a(6) = 5.0d0
thrs = 1.d-5
if(a(7).ne.0) thrs=a(7)
a1 = a(1)
rpar(ipr-iicode) = a1
if(a(1).lt.a(2).and.(a(1)*a(3)+a(4).lt.a(1)).or. &
   a(1).gt.a(2).and.(a(1)*a(3)+a(4).gt.a(1)).or. &
   a(1) .eq. a(2)) then
   write(6,2150)
2150    format(' Invalid step parameters for OPTIMIZE')
   goto 1  ! label:read_new_line
end if
jtot2x = ipar(IPAR_JTOT2)
write(6,2151) pcod(ipr)(1:lc),ipar(IPAR_JTOT1),(a(i),i = 1,6), &
              abs(thrs)
2151 format(' Optimization of ',(a),' for Jtot = ',i3,/, &
  ' Start:',f7.3,'  End:',f7.3,'  Factor:',f5.2, &
  ' Increment:',f7.3,/, &
  ' Average error limit:',f4.1,'%',/, &
  ' Maximum error limit:',f4.1,'%',/, &
  ' Threshold for S-matrix elements:',e10.2)
ipar(IPAR_JTOT2) = ipar(IPAR_JTOT1)
fnam1 = 'Joba'
fnam2 = 'Jobb'
! avoid too long filenames on cray (cos)
inam = 1
istep = 1
jobnam = fnam1
acclas = 1.d10
optval = a1
opti = .true.
optifl = .true.
write(6,255) pcod(ipr)(1:lc),a(1)
255 format(1x,(a),' = ',f7.3)
goto 500
2160 if(istep.ge.2) then
  call difs(fnam1,fnam2,1,0,acc,accmx,thrs,im,jm,ityp)
  code = '?'
  lcc = 1
  if(ityp .eq. 1) code = 'S real'
  if(ityp .eq. 1) lcc = 6
  if(ityp .eq. 2) code = 'S imaginary'
  if(ityp .eq. 2) lcc = 11
  if (thrs .ge. 0.) then
    code = 'S modulus'
    lcc  = 9
  end if
  write(6,2165) code(1:lcc), &
         acc,code(1:lcc),accmx,im,jm,code(1:lcc)
2165   format(' average difference between old and new ',(a), &
    ' = ',f10.2,'%',/, &
         ' Largest difference between old and new ',(a), &
    ' = ', f10.2,'%',/, &
     ' in (i = ',i2,' j = ',i2,') element of ',(a))
  if(acc.lt.acclas) then
    optval = a(1)
    optacc = acc
    optacm = accmx
    imx = im
    jmx = jm
    itx = ityp
  end if
  acclas = acc
end if
a(1) = a(1)*a(3)+a(4)
if((acclas.lt.a(5).and.accmx.lt.a(6)) &
   .or.(a(2).gt.a1.and.a(1).gt.a(2)).or. &
   (a(2).lt.a1.and.a(1).lt.a(2))) then
   code = '?'
   lcc = 1
   if(itx .eq. 1) code = 'S real'
   if(itx .eq. 1) lcc = 6
   if(itx .eq. 2) code = 'S imaginary'
   if(itx .eq. 2) lcc = 11
    if (thrs .ge. 0.) then
      code = 'S modulus'
      lcc  = 9
    end if
   rpar(ipr-iicode) = optval
   write(6,2170) pcod(ipr)(1:lc),optval,code(1:lcc), &
      optacc,code(1:lcc),optacm,imx,jmx,code(1:lcc)
2170    format(' optimized value for ',a,' = ',g11.4,/, &
   ' average difference in old and new ',(a),' is', &
     f10.2,'%',/, &
   ' Largest difference in old and new ',(a),' is', &
     f10.2,'%',/, &
     ' in (i = ',i2,' j = ',i2,') element of ',(a))
   opti = .false.
   ipar(IPAR_JTOT2) = jtot2x
   goto 1  ! label:read_new_line
end if
rpar(ipr-iicode) = a(1)
if(inam .eq. 1) then
  jobnam = fnam2
  inam = 2
else
  jobnam = fnam1
  inam = 1
end if
istep = istep+1
write(6,255) pcod(ipr)(1:lc),a(1)
goto 500
!.....tensor cross sections
!  tenxsc,jobfile,maxn,iframe,in1,in2,ienerg,jtotend,minj,maxj
2300 call get_token(line,l,fnam1,lc)
if(fnam1 .eq. ' ') fnam1 = jobnam
call lower(fnam1)
call upper(fnam1(1:1))
do 2310 i = 1,9
a(i) = 0.d0
if(l .eq. 0) goto 2310
call get_token(line,l,code,lc)
call assignment_parse(code(1:lc),empty_var_list,j,a(i))
2310 continue
#if defined(HIB_MAC)
call exit
write(6,2320)
2320 format(' Sorry, tensor cross sections not yet implemented')
#endif
#if defined(HIB_UNIX)
call tenopa(fnam1,a)
#endif
goto 1  ! label:read_new_line
!....nnout must be preceded by jout
2400 write (6, 2410)
2410 format(' To change NNOUT, enter the command line',/, &
  '    jout,nnout,jout(1),...,jout(iabs(nnout))' )
goto 1  ! label:read_new_line
!.....m-resolved cross sections
!  mrcrs,jobfile,ienerg
2500 call get_token(line,l,fnam1,lc)
if(fnam1 .eq. ' ') fnam1 = jobnam
call lower(fnam1)
call upper(fnam1(1:1))
do 2510 i = 1,1
a(i) = 0
if(l .eq. 0) goto 2510
call get_token(line,l,code,lc)
call assignment_parse(code(1:lc),empty_var_list,j,a(i))
2510 continue
#if defined(HIB_MAC)
call exit
write(6,2520)
2520 format(' M-resolved cross sections not yet implemented')
#endif
#if defined(HIB_UNIX) || defined(HIB_CRAY)
call mrcrs(fnam1,a)
#endif
goto 1  ! label:read_new_line
! printc : print selected integral cross sections from ics file
2600 call get_token(line,l,fnam1,lc)
if(fnam1 .eq. ' ') fnam1 = jobnam
call lower(fnam1)
call upper(fnam1(1:1))
do 2610 i = 1,8
a(i) = 0
if(l .eq. 0) goto 2610
call get_token(line,l,code,lc)
call assignment_parse(code(1:lc),empty_var_list,j,a(i))
2610 continue
if(ibasty.eq. 5) write (6, 2611)
2611 format &
 (' *** printed triplet pi cross sections should be checked')
if(ibasty.eq. 6) write (6, 2612)
2612 format &
 (' *** printed symmetric top cross sections should be checked')
!
!  eliminate special call for ibasty = 4
call prsg(fnam1,a)
!      if(ibasty.ne.4) call prsg(fnam1,a)
!      if(ibasty.eq.4) call prsgpi(fnam1,a)
goto 1  ! label:read_new_line
! print selected partial cross sections from pcs file
2650 call get_token(line,l,fnam1,lc)
if(fnam1 .eq. ' ') fnam1 = jobnam
call lower(fnam1)
call upper(fnam1(1:1))
do 2660 i = 1,8
a(i) = 0
if(l .eq. 0) goto 2660
call get_token(line,l,code,lc)
call assignment_parse(code(1:lc),empty_var_list,j,a(i))
2660 continue
call readpc(fnam1, a, scmat, nmax)
goto 1  ! label:read_new_line
!  psi(wavefunction calculation),jobfile,mchannel
!  flux calculation,jobfile,mchannel,iflux,thresh,iprint
2800 call get_token(line,l,fnam1,lc)
call lower(fnam1)
call upper(fnam1(1:1))
if(fnam1 .eq. ' ') fnam1 = jobnam
do 2810 i = 1,10
a(i) = 0
if(l .eq. 0) goto 2810
call get_token(line,l,code,lc)
call assignment_parse(code(1:lc),empty_var_list,j,a(i))
2810 continue
iflux=a(1)
if (a(2) .eq. 0.d0) iflux=2
call psi(fnam1,a)
goto 1  ! label:read_new_line
! adiabatic energy calculation, jobfile
2850 call get_token(line,l,fnam1,lc)
call lower(fnam1)
call upper(fnam1(1:1))
if(fnam1 .eq. ' ') fnam1 = jobnam
call get_token(line,l,code,lc)
if (code .eq. ' ') then
   l1 = 1
   l2 = 10
else
   read (code, *, err=2860, end=2860) l2
   call get_token(line,l,code,lc)
   if (code .eq. ' ') then
      l1 = 1
   else
      l1 = l2
      read (code, *, err=2860, end=2860) l2
   end if
end if
call eadiab1(fnam1,l1,l2)
goto 1  ! label:read_new_line
2860 write (6, *) 'Parameters to EADIAB cannot be recognized'
goto 1  ! label:read_new_line
!  print out system parameters
2900 call sys_conf
goto 1  ! label:read_new_line

end

end module mod_hinput
