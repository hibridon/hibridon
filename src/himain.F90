#include "assert.h"
program logair
use mod_com, only: com_file, com
use mod_clseg, only: allocate_clseg
use mod_cobuf, only: allocate_cobuf
use mod_cofil, only: allocate_cofil
use mod_conlam, only: allocate_conlam

use mod_cdbf, only: allocate_cdbf

use mod_cosout, only: allocate_cosout
use mod_coiout, only: allocate_coiout
use mod_cocent, only: allocate_cocent
use mod_coeint, only: allocate_coeint
use mod_coj12, only: allocate_coj12
use mod_coj12p, only: allocate_coj12p
use mod_covvl, only: allocate_covvl
use mod_cofact, only: allocate_cofact

use mod_coener, only: allocate_coener

use mod_coatpi, only: allocate_coatpi
use mod_coatpr, only: allocate_coatpr
use mod_coatp1, only: allocate_coatp1
use mod_coatp2, only: allocate_coatp2
use mod_coatp3, only: allocate_coatp3
use mod_coz, only: allocate_coz, z
use mod_cow, only: allocate_cow, w
use mod_cozmat, only: allocate_cozmat, zmat
use mod_coamat, only: allocate_coamat, amat
use mod_cobmat, only: allocate_cobmat, bmat
use mod_cotq1, only: allocate_cotq1, tq1
use mod_cotq2, only: allocate_cotq2, tq2
use mod_cotq3, only: allocate_cotq3, tq3
use mod_cojq, only: allocate_cojq, jq
use mod_colq, only: allocate_colq, lq
use mod_coinq, only: allocate_coinq, inq
use mod_cojhld, only: allocate_cojhld, jhold
use mod_coehld, only: allocate_coehld, ehold
use mod_coinhl, only: allocate_coinhl, inhold
use mod_coisc1, only: allocate_coisc1, isc1
use mod_coisc2, only: allocate_coisc2, isc2
use mod_coisc3, only: allocate_coisc3, isc3
use mod_coisc4, only: allocate_coisc4, isc4
use mod_coisc5, only: allocate_coisc5
use mod_coisc6, only: allocate_coisc6
use mod_coisc7, only: allocate_coisc7
use mod_coisc8, only: allocate_coisc8
use mod_coisc9, only: allocate_coisc9
use mod_coisc10, only: allocate_coisc10
use mod_coisc11, only: allocate_coisc11
use mod_coisc12, only: allocate_coisc12
use mod_colsc1, only: allocate_colsc1, lsc1
use mod_cosc1, only: allocate_cosc1, sc1
use mod_cosc2, only: allocate_cosc2, sc2
use mod_cosc3, only: allocate_cosc3, sc3
use mod_cosc4, only: allocate_cosc4, sc4
use mod_cosc5, only: allocate_cosc5, sc5
use mod_cosc6, only: allocate_cosc6, sc6
use mod_cosc7, only: allocate_cosc7, sc7
use mod_cosc8, only: allocate_cosc8, sc8
use mod_cosc9, only: allocate_cosc9, sc9
use mod_cosc10, only: allocate_cosc10
use mod_cosc11, only: allocate_cosc11
use mod_coeig, only: allocate_coeig
use mod_coeig2, only: allocate_coeig2
#if defined(HIB_UNIX)
use mod_cokaux, only: allocate_cokaux
#endif
use mod_cotble, only: allocate_cotble
use mod_coqvec, only: allocate_coqvec
use mod_coqvec2, only: allocate_coqvec2
use mod_codim, only: allocate_codim, mairy, mmax
use mod_comxbs, only: allocate_comxbs
use mod_comxm, only: allocate_comxm

use mod_cosyr, only: allocate_cosyr
use mod_cosys, only: allocate_cosys
use mod_cosysl, only: allocate_cosysl
use mod_cosysi, only: allocate_cosysi
use mod_cosysr, only: allocate_cosysr
use mod_par, only: allocate_par
use mod_flow, only: flow
!**********************************************************************
!***   this code is not released for general public use            ****
!***   all use must be by specific prior arrangement with:         ****
!***     millard alexander, department of chemistry,               ****
!***     university of maryland, college park, md, 20742           ****
!***     tel: 1.301.405.1823; email: mha@umd.edu                   ****
!***   no part of this program may be copied or used for other     ****
!***   purposes without the author's permission.                   ****
!**********************************************************************
!  ***  driver for log-derivative integration ***
!  author:  millard alexander
!  revised:  24-jan-2012 by p. dagdigian (added common blocks /coj12/, /co12p/)
!  revised:  replaced statements before call to propag to accommodate
!    q.ma's revised version of bound (28-jun-2013, p.dagdigian)
!  revised: 16-jun-2019 by p. dagdigian (increased kmxbas to 30)
!  added common block from ba3p2s basis routine
!  current revision date:  19-jun-2019 (p. dagdigian)
!
implicit none
#if defined(HIB_UNIX_IBM) || defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
character *40 test
#endif
!  ----------------------------------------------------------
!  in the following parameter statements:
!     kmxbas is maximum number of included basis routines, this should
!     be updatated as basis routines are added
!     kmax is maximum number of channels set at run time
!     klammx is maximum number of anisotropic terms in potential
!     kfact is maximum value for which log(n!) is computed
!     ken is number of total energies allowed
!     kmxpho is maximum number of different initial states allowed in
!     photodissociation calculation
!     kout is number of different values of rotational quantum number
!     for which s-matrix will be stored on disk

!  revised by p. dagdigian (13-dec-2019)
integer :: kmax
integer :: kairy
integer :: ktri
integer, parameter :: kbig = 10

! size of scratch array for matrix inversion
integer :: kaux

! size of scratch array for matrix factorization and generalized
!   eigenvalue determination
integer :: kaux3 

integer :: kq2
integer :: kqmax

! krotmx set the maximum size of the asym top e.fn expansion
integer :: krotmx

! modified klammx (pjd - 17-jan-2019)
! increased kfact (pjd - 15-dec-2020)
integer, parameter :: klammx = 2000
integer, parameter :: kfact = 3500
integer, parameter :: kout = 50
integer, parameter :: ken = 25

integer, parameter :: kmxpho = 3
integer, parameter :: knphot = 1
integer, parameter :: kmxbas = 30 ! number of base types
integer, parameter :: kmaxpar = 150 ! max number of base specific parameters
integer :: men
!
!  ----------------------------------------------------------
common /comom/  xmom(3), imom(13)
real(8) :: xmom
integer :: imom

integer :: arg_index
integer :: stat
character(len=32) :: arg


kmax = 0

arg_index = 1
do while( arg_index <= command_argument_count() )
      call get_command_argument(arg_index, arg)

      select case (arg)
          case ('-k', '--kmax')
              arg_index = arg_index + 1
              if ( arg_index > command_argument_count() ) then
                  print '(2a, /)', 'missing value for <max_num_channels> ', arg
                  call print_help()
                  stop 1
              end if
              call get_command_argument(arg_index, arg)
              read(arg,*,iostat=stat) kmax 

          case ('-h', '--help')
              call print_help()
              stop 1

          case ('-c', '--com')
            arg_index = arg_index + 1
            if ( arg_index > command_argument_count() ) then
                print '(2a, /)', 'missing value for <command_file> ', arg
                call print_help()
                stop 1
            end if
            call get_command_argument(arg_index, com_file)
            com = .true.

          case default
              print '(2a, /)', 'unrecognised command-line option: ', arg
              call print_help()
              stop 1
      end select
      arg_index = arg_index + 1
end do
if (kmax == 0) then
      print '(a, /)', 'error : --kmax argument is mandatory.'
      call print_help()
      stop 1
end if

kairy = kmax
ktri = kmax * (kmax + 1) / 2
kq2 = 2 * kmax
kqmax = kmxpho * kmax
! krotmx set the maximum size of the asym top e.fn expansion
! increased factor multiplying kmax below and narray from 12 to 100 (pjd:  19-sep-2017)
krotmx = 500 * kmax

! set size of scratch array for matrix inversion
#if defined(HIB_UNIX_IBM)
kaux = 100 * kmax
#elif defined(HIB_UNIX) && !defined(HIB_UNIX_IBM)
! set size of scratch array for matrix inversion with lapack routines
! warning, this assumes a blocksize of 64
kaux = 128 * kmax
#else
#error "hibridon doesn't know what value to choose for kaux on this architecture"
#endif

! set size of scratch array for matrix factorization and generalized
!   eigenvalue determination
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
kaux3 = max( 66 * kmax, 2 * kmax * kmax + 1)
#else
#error "hibridon doesn't know what value to choose for kaux on this architecture"
#endif

!   square matrices and vectors
call allocate_coener(ken)
call allocate_clseg()
call allocate_cobuf()
call allocate_cofil()
call allocate_conlam(klammx)
call allocate_cdbf()
call allocate_coatpi(krotmx, 100)
call allocate_coatpr(krotmx)
call allocate_coatp1(krotmx)
call allocate_coatp2(krotmx)
call allocate_coatp3(krotmx)
call allocate_coz(kmax)
call allocate_cow(kmax)
call allocate_cozmat(kmax)
call allocate_coamat(kmax)
call allocate_cobmat(kairy)
call allocate_cotq1(T_MATRIX_SIZE)
call allocate_cotq2(T_MATRIX_SIZE)
call allocate_cotq3(T_MATRIX_SIZE)
call allocate_cojq(kmax)
call allocate_colq(kmax)
call allocate_coinq(kmax)
call allocate_cojhld(kmax)
call allocate_coehld(kmax)
call allocate_coinhl(kmax)
call allocate_coisc1(kmax)
call allocate_coisc2(kmax)
call allocate_coisc3(kmax)
call allocate_coisc4(kmax)
call allocate_coisc5(kmax)
call allocate_coisc6(kmax)
call allocate_coisc7(kmax)
call allocate_coisc8(kmax)
call allocate_coisc9(kmax)
call allocate_coisc10(kmax)
call allocate_coisc11(kmax)
call allocate_coisc12(kmax)
call allocate_colsc1(kmax)
call allocate_cosc1(kmax)
call allocate_cosc2(kmax)
call allocate_cosc3(kmax)
call allocate_cosc4(kmax)
call allocate_cosc5(kmax)
call allocate_cosc6(kmax)
call allocate_cosc7(kmax)
call allocate_cosc8(kmax)
call allocate_cosc9(kmax)
call allocate_cosc10(kmax)
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
call allocate_cosc11(kaux3)
#endif
#if defined(HIB_UNIX) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
call allocate_cosc11(kaux)
#endif
!   total matrix and vector storage required is:
!     7 kmax**2 + 25 kmax + v2 storage + kfact -- without airy integration
!     8 kmax**2 + 25 kmax + v2 storage + kfact -- with airy integration
!   if linked with -b option, then storage requirements drop to
!     5 kmax**2 + 25 kmax + v2 storage + kfact -- with airy integration
!
!  parameter below sets maximum size of asymmetric top basis fn expansion
call allocate_cosout(kout)
call allocate_coiout(kout)
call allocate_cocent(kmax)
call allocate_coeint(kmax)
call allocate_coj12(kmax)
call allocate_coj12p(kmax)
call allocate_covvl(klammx)
call allocate_cofact(kfact)
call allocate_coeig()
call allocate_coeig2()
#if defined(HIB_UNIX_IBM)
call allocate_cokaux(anaux=max(kaux,1800))
#endif
call allocate_cotble(kfact)
call allocate_coqvec(kqmax, kmxpho, knphot)
call allocate_coqvec2(kq2)
call allocate_codim(kairy, kmax, kbig)
call allocate_comxbs(kmxbas)
call allocate_comxm()

call allocate_cosysi(kmaxpar)
call allocate_cosysl(kmaxpar)
call allocate_cosysr(kmaxpar)
call allocate_cosyr(kmaxpar)
call allocate_cosys(2*kmaxpar+3)
call allocate_par()
!
men = ken
!  calculate array containing logs of factorials for use in vector
!  coupling coefficient routines:  factlg(i)=log(fact(i-1))
call factlg(kfact)
! start scattering calculation
call flow (z, w, zmat, amat, bmat, jq, lq, inq, jhold, ehold, &
           inhold, isc1, isc2, isc3, isc4, lsc1, &
           sc2, sc1, sc3, sc4, sc5, &
           sc6, sc7, sc8, sc9, tq1, tq2, tq3, men, mmax, mairy)
call finit ! Closes all I/O units
contains
subroutine print_help()
print '(a, /)', 'command-line options:'
print '(a)',    '  -k, --kmax <max_num_channels>  defines the max number of channels'
print '(a)',    '  -c, --com <command_file>  specifies a command file instead of input redirection'
print '(a, /)', '  -h, --help        print usage information and exit'
end subroutine print_help
end