#include "assert.h"
module mod_com
   implicit none
   character(len=300) :: com_file 
   logical :: com = .false.
end module mod_com

module mod_cosout
   !  variables in common block /cosout/
   !    nnout:     number of different rotational levels for which s-matrix
   !               elements are to be saved in file 13,
   !    jout(i):   values of rotational angular momentum for these lvels
   ! *  variables in common block /cosout/
   ! *    nnout:     number of different rotational levels for which s-matrix
   ! *               elements are to be saved in file nfile
   ! *    jout(i):   values of rotational angular momentum for these levels
   ! *               if nnout is positive, than an s-matrix element will be saved
   ! *               only if both the initial and final quantum numbers correspond
   ! *               to one of the values of jout(i)
   ! *               if nnout is negative, then every column of the s-matrix for
   ! *               which the initial quantum numbers correspond to one of the
   ! *               values of jout(i) will be printed
   ! *  variables in common block /cosout/
   ! *    nnout:     number of different rotational levels for which s-matrix
   ! *               elements are to be saved in files smat1, smat2, ...
   ! *    jout(i):   values of rotational angular momentum for these lvels
   ! *  variables in common block /cosout/
   ! *    nnout:     number of different rotational levels for which s-matrix
   ! *               elements are to be saved in file nfile
   ! *    jout(i):   values of rotational angular momentum for these levels
   ! *               if nnout is positive, than an s-matrix element will be saved
   ! *               only if both the initial and final quantum numbers correspond
   ! *               to one of the values of jout(i)
   ! *               if nnout is negative, then every column of the s-matrix for
   ! *               which the initial quantum numbers correspond to one of the
   ! *               values of jout(i) will be printed
   ! from hiiolib_c.F90
   ! *  variables in common block /cosout/
   ! *    nnout:     number of different rotational levels for which s-matrix
   ! *               elements are to be saved in files smat1, smat2, ...
   ! *    jout(i):   values of rotational angular momentum for these lvels
   ! *  variables in common block /cosout/
   ! *    nnout:     number of different rotational levels for which s-matrix
   ! *               elements are to be saved in file nfile
   ! *    jout(i):   values of rotational angular momentum for these levels
   ! *               if nnout is positive, than an s-matrix element will be saved
   ! *               only if both the initial and final quantum numbers correspond
   ! *               to one of the values of jout(i)
   ! *               if nnout is negative, then every column of the s-matrix for
   ! *               which the initial quantum numbers correspond to one of the
   ! *               values of jout(i) will be printed
   implicit none
   integer, dimension(:), allocatable :: jout
   integer, allocatable               :: nnout
   contains
   subroutine allocate_cosout(aout)
      integer, intent(in) :: aout
      allocate(jout(aout)) ; allocate(nnout)
   end subroutine allocate_cosout
end module mod_cosout

module mod_coiout
   implicit none
   integer, dimension(:), allocatable :: indout
   integer, allocatable               :: niout
   contains
   subroutine allocate_coiout(aout)
      integer, intent(in) :: aout
      allocate(indout(aout)) ; allocate(niout)
   end subroutine allocate_coiout
end module mod_coiout

module mod_conlam
   ! *  variables
   ! *    nlam:      the number of angular coupling terms actually used
   ! *    nlammx:    the maximum number of angular coupling terms allowed
   ! *    lamnum:    number of non-zero v2 matrix elements for each lambda
   ! *               lamnum is an array of dimension nlammx
   implicit none
   integer, dimension(:), allocatable :: lamnum
   integer, allocatable               :: nlam, nlammx
   contains
   subroutine allocate_conlam(n)
      integer, intent(in) :: n
      allocate(lamnum(n)) ; allocate(nlam) ; allocate(nlammx)
      nlammx = n
      nlam = 0
   end subroutine allocate_conlam
end module mod_conlam

module mod_cocent
   ! variables in this module
   !   cent:      array containing centrifugal barrier of each channel
   implicit none
   real(8), dimension(:), allocatable :: cent
   contains
   subroutine allocate_cocent(amax)
      integer, intent(in) :: amax
      allocate(cent(amax)) ;
   end subroutine allocate_cocent
end module mod_cocent

module mod_coeint
   ! variables in this module
   !   eint:     array containing channel energies (in hartree)
   !             the zero of energy is assumed to be the 0(0,0) level
   !             + energy of j2min level of linear molecule
   implicit none
   real(8), dimension(:), allocatable :: eint
   contains
   subroutine allocate_coeint(amax)
      integer, intent(in) :: amax
      allocate(eint(amax)) ;
   end subroutine allocate_coeint
end module mod_coeint

module mod_coj12
   ! variables in this module
   !    j12:   array containing vector sum of ja + j (similar
   !           situation with molecule-molecule scattering,
   !           see hibastpln basis routine)
   !    j12:      array containing vector sum of ja + j (similar
   !              here j12 is lower-case j in alexander, dagdigian, klos notes
   !              situation with molecule-molecule scattering,
   !              see hibastpln basis routine)
   !    j12:       vector resultant of j1 + j2
   !    j12:       array containing the j12 quantum number of each channel  
   !    j12:       array containing vector sum of j1 + j2 for molecule-molecule
   !               systems and open-shell atom - molecule systems
   
   implicit none
   integer, dimension(:), allocatable :: j12
   contains
   subroutine allocate_coj12(amax)
      integer, intent(in) :: amax
      allocate(j12(amax)) ;
   end subroutine allocate_coj12
end module mod_coj12

module mod_coj12p
   !  variables in this module
   !    j12p:  arrays containing vector sum of j1 + j2 for molecule-molecule 
   !               systems and open-shell atom - molecule systems
   implicit none
   integer, dimension(:), allocatable :: j12pk
   contains
   subroutine allocate_coj12p(amax)
      integer, intent(in) :: amax
      allocate(j12pk(amax)) ;
   end subroutine allocate_coj12p
end module mod_coj12p

module mod_covvl
   !  variable in this module
   !    vvl:       array to store r-dependence of each angular term in the
   !               potential   
   implicit none
   real(8), dimension(:), allocatable :: vvl
   contains
   subroutine allocate_covvl(alammx)
      integer, intent(in) :: alammx
      allocate(vvl(alammx)) ;
   end subroutine allocate_covvl
end module mod_covvl

module mod_cofact
   !  variable in this module
   !    si:        table of logarithms of factorials (1-based index)
   !       si(n) = ln((n-1)!) 
   !       si(0) = ln(0!) = 0.0
   !       si(1) = ln(1!) = ln(1) + ln(0)
   !       si(2) = ln(2!) = ln(2) + ln(1) + ln(0)

   implicit none
   real(8), dimension(:), allocatable :: si
   contains
   subroutine allocate_cofact(afact)
      integer, intent(in) :: afact
      allocate(si(0:afact-1)) ;
   end subroutine allocate_cofact
end module mod_cofact

module mod_coener
   ! *  variable in this module
   ! *    energ:     array containing total energies (in cm-1) at which scattering
   ! *               calculations are to be performed
   implicit none
   real(8), dimension(:), allocatable :: energ
   contains
   subroutine allocate_coener(n)
      integer, intent(in) :: n
      allocate(energ(n)) ;
   end subroutine allocate_coener
end module mod_coener

module mod_cdbf
   ! buffer of 32 bits words (integers)
   implicit none
   integer, allocatable :: ldbuf
   integer, allocatable :: libuf                ! length of ibuf
   integer, allocatable :: ibfil                ! current fileid 
   integer, allocatable :: ibrec                ! current record
   integer, allocatable :: ibof                 ! current position in the buffer (idbuf[1:ibof-1] is expected to be already filled)
   integer, allocatable :: ibstat               ! 0: buffer already saved to disc, 1: buffer contains unsaved data
   integer, dimension(:), allocatable :: idbuf  ! buffer of llbuf integers
   integer, allocatable :: llbuf
   integer :: ibadr                             ! where to write idbuf on disc (word position relative to the current record)
   contains
   subroutine allocate_cdbf()
      allocate(llbuf)
#if defined(HIB_CRAY) || defined(HIB_FPS)
      llbuf=512
#elif defined(HIB_VAX) || defined(HIB_UNIX) || defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IBM)
      llbuf=1024
#elif defined(HIB_UNIVAC)
      llbuf=1792
#else
#error "unknow architecture for setting llbuf"
#endif
      allocate(ldbuf)
      allocate(libuf)
      allocate(ibfil)
      allocate(ibrec)
      allocate(ibof)
      allocate(ibstat)
      allocate(idbuf(llbuf))
      ibadr = 0
   end subroutine allocate_cdbf
end module mod_cdbf


module mod_clseg
   !  variables in this module
   !    lseg:      number of integer words per disc sector (length of a segment)
   !    intrel:    number of integer words per real words
   !    lchar:     number of characters per integer word
   implicit none
   integer, allocatable :: lseg
   integer, allocatable :: intrel ! size of real(8) compared to the size of an integer
   integer, allocatable :: lchar
   contains
   subroutine allocate_clseg()
      allocate(lseg) ; allocate(intrel) ; allocate(lchar) ;
#if defined(HIB_UNIX_DARWIN64)
      lseg = 1024
      intrel = 1
      lchar = 8
#endif
#if defined(HIB_UNIX) && !defined(HIB_UNIX_DARWIN)
      lseg = 1024
      intrel = 2
      lchar = 4
#endif
#if defined(HIB_CRAY)
      lseg = 512
      intrel = 1
      lchar = 8
#endif
   end subroutine allocate_clseg
end module mod_clseg

module mod_cobuf
   ! *  variables 
   ! *    lbuf:      length of i/o buffer (should be multiple of lseg)
   ! *    ibuf:      buffer used in several i/o routines
   implicit none
   integer, allocatable               :: lbuf
   integer, dimension(:), allocatable :: ibuf
   contains
   subroutine allocate_cobuf()
      allocate(lbuf)
#if defined(HIB_CRAY)
   lbuf = 512
#else
   lbuf = 1024
#endif
      allocate(ibuf(lbuf))
   end subroutine allocate_cobuf
end module mod_cobuf

module mod_cofil
   ! *  variables
   ! *     nfl:      logical unit number
   ! *     iofbuf:   buffer offset
   ! *     maxrec:   number of records on file associated with nfl
   ! *     iofrec:   offset in file
   implicit none
   integer, dimension(:), allocatable :: maxrec, iofrec
   integer, allocatable               :: nfl, iofbuf, nwrec
   contains
   subroutine allocate_cofil()
      allocate(maxrec(60)) ; allocate(iofrec(60)) ;
      allocate(nfl)        ; allocate(iofbuf)      ; allocate(nwrec)
   end subroutine allocate_cofil
end module mod_cofil


module mod_coatpi
!  variables in this module
!     narray:   maximum size of asymmetric top basis fn expansion
!               set to 500 (see krotmx above) in himain 
!     isiz:    length of eigenfunction expansion for each rot. level
   implicit none
   integer, dimension(:), allocatable :: isiz
   integer, allocatable               :: narray
   contains
   subroutine allocate_coatpi(n, anarray)
      integer, intent(in) :: n
      integer, intent(in) :: anarray
      allocate(isiz(n)) ; allocate(narray)
      narray = anarray
   end subroutine allocate_coatpi
end module mod_coatpi

module mod_coatp3
   ! variable in this module
   !   isizh:      temporary storage for length (related to isiz)
   implicit none
   integer, dimension(:), allocatable :: isizh
   contains
   subroutine allocate_coatp3(n)
      integer, intent(in) :: n
      allocate(isizh(n))
   end subroutine allocate_coatp3
end module mod_coatp3

module mod_coatpr
   ! variable in this module
   !      c:        expansion coefficients for asymmetric top rotor wave fns.   
   ! *   c:          expansion coefficienfs for the rotational/fine-structure
   ! *               wave functions in the symmetrized basis:
   ! *               (1) sigma=1, eps=+1; (2) sigma=1, eps=-1, (3) sigma = 0, eps=+1
   ! *  variable in common block /coatpr/
   ! *   c:          expansion coefficients for spherical top rotor wave fns.
   ! *               in a signed-k symmmetric top basis |j,k,m>.
   ! *               if c(+k) = c(-k), then eps = +1, or if c(+k) = -c(-k), eps = -1.
   ! *               the expansion coefficients for a given eps stored in the array c are:
   ! *               c(k=-j), c(k=-j+2), ... c(k=j-10, c(k=j).
   ! *  variable in common block /coatpr/
   ! *   c:          expansion coefficients for asymmetric top rotor wave fns.
   ! *               in a symmetrized symmetric top basis.  these basis functions
   ! *               are given by green, jcp 64, 3463 (1976) as:
   ! *                  |j,k,m,s> = [2*(1+delta(k,0)]^-1 *
   ! *                      (|j,k,m> + (-1)*s * |j,-k,m>)
   ! *               (note that green uses eps = (-1)*s for the symmetry index.)
   ! *               in this basis, the asymmetric top hamiltonian block diagonalizes
   ! *               into 4 groups:  (1) k even, s = +1, (2) k even, s= -1,
   ! *               (3) k odd, s = =1, and (4) k odd, s = _1.
   ! *               the expansion coefficients stored in the array c are:
   ! *               c(k=0), c(k=2), ... c(k=j) for even k, and
   ! *               c(k=1), c(k=3), ... c(k=j) for odd k.
   ! *               the levels are denoted in the channel list by j and the
   ! *               symmetry index is.  the prolate-limit projection quantum number kp
   ! *               equals abs(is), and s is the sign of is.  the oblate limit projection
   ! *               quantum number ko can be obtained as follows:
   ! *                  if is >= 0, then ko = j - kp
   ! *                  if is < 0, then ko = j + 1 - kp
   ! *               if bastst=t, the values of j, is, kp, and ko are printed out, as
   ! *               well as the values of the expansion coefficients

   
   implicit none
   real(8), dimension(:), allocatable :: c
   contains
   subroutine allocate_coatpr(n)
      integer, intent(in) :: n
      allocate(c(n))
   end subroutine allocate_coatpr
end module mod_coatpr

module mod_coatp1
   ! *  variable in this module
   ! *   ctemp:      temporary storage for rot. e.fn coeffs.   
   implicit none
   real(8), dimension(:), allocatable :: ctemp
   contains
   subroutine allocate_coatp1(n)
      integer, intent(in) :: n
      allocate(ctemp(n))
   end subroutine allocate_coatp1
end module mod_coatp1

module mod_coatp2
   ! *  variable in this module
   ! *   chold:      ditto   
   implicit none
   real(8), dimension(:), allocatable :: chold
   contains
   subroutine allocate_coatp2(n)
      integer, intent(in) :: n
      allocate(chold(n))
   end subroutine allocate_coatp2
end module mod_coatp2

module mod_coz
   implicit none
   real(8), dimension(:,:), allocatable, target :: z
   real(8), dimension(:), pointer :: z_as_vec  ! z matrix viewed as a vector
   contains
   subroutine allocate_coz(n)
      use, intrinsic :: ISO_C_BINDING
      integer, intent(in) :: n
      allocate(z(n,n))

      call C_F_POINTER (C_LOC(z), z_as_vec, [n*n])
   end subroutine allocate_coz
end module mod_coz

module mod_cow
   implicit none
   real(8), dimension(:,:), allocatable, target :: w
   real(8), dimension(:), pointer :: w_as_vec  ! w matrix viewed as a vector
   contains
   subroutine allocate_cow(n)
      use, intrinsic :: ISO_C_BINDING
      integer, intent(in) :: n
      allocate(w(n,n))

      call C_F_POINTER (C_LOC(w), w_as_vec, [n*n])
   end subroutine allocate_cow
end module mod_cow

module mod_cozmat
   implicit none
   real(8), dimension(:,:), allocatable, target :: zmat
   real(8), dimension(:), pointer :: zmat_as_vec  ! zmat matrix viewed as a vector
   contains
   subroutine allocate_cozmat(n)
      use, intrinsic :: ISO_C_BINDING

      integer, intent(in) :: n
      allocate(zmat(n,n))
      call C_F_POINTER (C_LOC(zmat), zmat_as_vec, [n*n])
   end subroutine allocate_cozmat
end module mod_cozmat

module mod_coamat
   implicit none
   ! note: these members are not used at the same time
   ! (they used to be a different view of the same area in memory
   ! when this module was a common block)
   ! todo : reduce memory usage by allocating the data only while being used
   real(8), dimension(:,:), allocatable, target :: amat
   real(8), dimension(:), allocatable :: toto
   integer, dimension(:), allocatable :: ietmp
   real(8), dimension(:), allocatable :: zbuf
   real(8), dimension(:), allocatable :: simag2
   real(8), dimension(:), allocatable :: psir
   integer, dimension(:), allocatable :: labadr
   contains
   subroutine allocate_coamat(n)
      integer, intent(in) :: n
      allocate(amat(n,n))
      allocate(toto(n*n))
      allocate(ietmp(n))
      allocate(zbuf(n))
      allocate(simag2(n*n))
      allocate(psir(n*n))
      allocate(labadr(MAX_NJTOT))
   end subroutine allocate_coamat
end module mod_coamat

module mod_cobmat
   implicit none
   real(8), dimension(:,:), allocatable :: bmat
   real(8), dimension(:), allocatable   :: psii
   real(8), dimension(:), allocatable   :: sigma
   contains
   subroutine allocate_cobmat(n)
      integer, intent(in) :: n
      allocate(bmat(n,n))
      allocate(psii(n*n))
      allocate(sigma(n*n))
   end subroutine allocate_cobmat
end module mod_cobmat

module mod_cotq1
   ! stores t matrix
   implicit none
   ! note: these members are not used at the same time
   ! (they used to be a different view of the same area in memory
   ! when this module was a common block)
   real(8), dimension(:,:), allocatable :: tq1
   real(8), dimension(:,:,:), allocatable :: vec
   real(8), dimension(:), allocatable :: dpsir
   real(8), dimension(:), allocatable :: scmat
   contains
   subroutine allocate_cotq1(n)
      integer, intent(in) :: n
      allocate(tq1(n,n))
      allocate(vec(3,3,1))
      allocate(dpsir(n*n))  ! note : the size has been found by trial and error (with all tests passing)
      allocate(scmat(n*n))  ! note : the size has been found by trial and error (with all tests passing)
   end subroutine allocate_cotq1
end module mod_cotq1

module mod_cotq2
   ! note: these members are not used at the same time
   ! (they used to be a different view of the same area in memory
   ! when this module was a common block)
   ! stores t matrix
   implicit none
   real(8), dimension(:,:), allocatable :: tq2
   real(8), dimension(:), allocatable :: dpsii
   real(8), dimension(:), allocatable :: scmat
   contains
   subroutine allocate_cotq2(n)
      integer, intent(in) :: n
      allocate(tq2(n,n))
      allocate(dpsii(n*n))  ! note : the size has been found by trial and error (with all tests passing)
      allocate(scmat(n*n))  ! note : the size has been found by trial and error (with all tests passing)
   end subroutine allocate_cotq2
end module mod_cotq2

module mod_cotq3
   ! note: these members are not used at the same time
   ! (they used to be a different view of the same area in memory
   ! when this module was a common block)
   ! stores t matrix
   implicit none
   real(8), dimension(:,:), allocatable :: tq3
   real(8), dimension(:), allocatable :: scmat
   contains
   subroutine allocate_cotq3(n)
      integer, intent(in) :: n
      allocate(tq3(n,n))
      allocate(scmat(n*n))  ! note : the size has been found by trial and error (with all tests passing)
   end subroutine allocate_cotq3
end module mod_cotq3

module mod_cojq
   implicit none
   integer, dimension(:), allocatable :: jq
   contains
   subroutine allocate_cojq(n)
      integer, intent(in) :: n
      allocate(jq(n))
   end subroutine allocate_cojq
end module mod_cojq

module mod_colq
   implicit none
   integer, dimension(:), allocatable :: lq
   contains
   subroutine allocate_colq(n)
      integer, intent(in) :: n
      allocate(lq(n))
   end subroutine allocate_colq
end module mod_colq

module mod_coinq
   implicit none
   integer, dimension(:), allocatable :: inq
   contains
   subroutine allocate_coinq(n)
      integer, intent(in) :: n
      allocate(inq(n))
   end subroutine allocate_coinq
end module mod_coinq

module mod_cojhld
   implicit none
   integer, dimension(:), allocatable :: jhold
   contains
   subroutine allocate_cojhld(n)
      integer, intent(in) :: n
      allocate(jhold(n))
   end subroutine allocate_cojhld
end module mod_cojhld

module mod_coehld
   implicit none
   real(8), dimension(:), allocatable :: ehold
   integer, dimension(:), allocatable :: eholdint
   contains
   subroutine allocate_coehld(n)
      integer, intent(in) :: n
      allocate(ehold(n))
      allocate(eholdint(n))
   end subroutine allocate_coehld
end module mod_coehld

module mod_coinhl
   implicit none
   integer, dimension(:), allocatable :: inhold
   contains
   subroutine allocate_coinhl(n)
      integer, intent(in) :: n
      allocate(inhold(n))
   end subroutine allocate_coinhl
end module mod_coinhl

module mod_coisc1
   implicit none
   integer, dimension(:), allocatable :: isc1
   contains
   subroutine allocate_coisc1(n)
      integer, intent(in) :: n
      allocate(isc1(n))
   end subroutine allocate_coisc1
end module mod_coisc1

module mod_coisc2
   implicit none
   integer, dimension(:), allocatable :: isc2
   integer, allocatable :: nj
   contains
   subroutine allocate_coisc2(n)
      integer, intent(in) :: n
      allocate(isc2(n))
      allocate(nj)
   end subroutine allocate_coisc2
end module mod_coisc2

module mod_coisc3
   implicit none
   integer, dimension(:), allocatable :: isc3
   contains
   subroutine allocate_coisc3(n)
      integer, intent(in) :: n
      allocate(isc3(n))
   end subroutine allocate_coisc3
end module mod_coisc3

module mod_coisc4
   implicit none
   integer, dimension(:), allocatable :: isc4
   contains
   subroutine allocate_coisc4(n)
      integer, intent(in) :: n
      allocate(isc4(n))
   end subroutine allocate_coisc4
end module mod_coisc4

module mod_coisc5
   implicit none
   integer, dimension(:), allocatable :: isc5
   contains
   subroutine allocate_coisc5(n)
      integer, intent(in) :: n
      allocate(isc5(n))
   end subroutine allocate_coisc5
end module mod_coisc5

module mod_coisc6
   implicit none
   integer, dimension(:), allocatable :: isc6
   contains
   subroutine allocate_coisc6(n)
      integer, intent(in) :: n
      allocate(isc6(n))
   end subroutine allocate_coisc6
end module mod_coisc6

module mod_coisc7
   implicit none
   integer, dimension(:), allocatable :: isc7
   contains
   subroutine allocate_coisc7(n)
      integer, intent(in) :: n
      allocate(isc7(n))
   end subroutine allocate_coisc7
end module mod_coisc7

module mod_coisc8
   !  variable in this module
   !    isc8:     for CC calculations, on return contains atom+diatom angular mom
   !              for CD calculations, on return contains diatom angular momentum

   implicit none
   integer, dimension(:), allocatable :: isc8
   contains
   subroutine allocate_coisc8(n)
      integer, intent(in) :: n
      allocate(isc8(n))
   end subroutine allocate_coisc8
end module mod_coisc8

module mod_coisc9
   implicit none
   integer, dimension(:), allocatable :: isc9
   contains
   subroutine allocate_coisc9(n)
      integer, intent(in) :: n
      allocate(isc9(n))
   end subroutine allocate_coisc9
end module mod_coisc9

module mod_coisc10
   implicit none
   integer, dimension(:), allocatable :: isc10
   contains
   subroutine allocate_coisc10(n)
      integer, intent(in) :: n
      allocate(isc10(n))
   end subroutine allocate_coisc10
end module mod_coisc10

module mod_coisc11
   implicit none
   integer, dimension(:), allocatable :: isc11
   contains
   subroutine allocate_coisc11(n)
      integer, intent(in) :: n
      allocate(isc11(n))
   end subroutine allocate_coisc11
end module mod_coisc11

module mod_coisc12
   implicit none
   integer, dimension(:), allocatable :: isc12
   contains
   subroutine allocate_coisc12(n)
      integer, intent(in) :: n
      allocate(isc12(n))
   end subroutine allocate_coisc12
end module mod_coisc12

module mod_colsc1
   implicit none
   logical, dimension(:), allocatable :: lsc1
   contains
   subroutine allocate_colsc1(n)
      integer, intent(in) :: n
      allocate(lsc1(n))
   end subroutine allocate_colsc1
end module mod_colsc1

module mod_cosc1
   implicit none
   real(8), dimension(:), allocatable :: sc1
   contains
   subroutine allocate_cosc1(n)
      integer, intent(in) :: n
      allocate(sc1(n))
   end subroutine allocate_cosc1
end module mod_cosc1

! general purpose array of scalars
! note : sc2 and sc2int are never used at the same time, as they initially were 2 views of the same data
module mod_cosc2
   implicit none
   real(8), dimension(:), allocatable :: sc2
   integer, dimension(:), allocatable :: sc2int
   contains
   subroutine allocate_cosc2(n)
      integer, intent(in) :: n
      allocate(sc2(n))
      allocate(sc2int(n))
   end subroutine allocate_cosc2
end module mod_cosc2

! general purpose array of scalars
! note : sc3 and sc3int are never used at the same time, as they initially were 2 views of the same data
module mod_cosc3
   implicit none
   real(8), dimension(:), allocatable :: sc3
   integer, dimension(:), allocatable :: sc3int
   contains
   subroutine allocate_cosc3(n)
      integer, intent(in) :: n
      allocate(sc3(n))
      allocate(sc3int(n))
   end subroutine allocate_cosc3
end module mod_cosc3

! general purpose array of scalars
! note : sc4 and sc4int are never used at the same time, as they initially were 2 views of the same data
module mod_cosc4
   implicit none
   real(8), dimension(:), allocatable :: sc4
   integer, dimension(:), allocatable :: sc4int
   contains
   subroutine allocate_cosc4(n)
      integer, intent(in) :: n
      allocate(sc4(n))
      allocate(sc4int(n))
   end subroutine allocate_cosc4
end module mod_cosc4

! general purpose array of scalars
! note : sc5 and sc5int are never used at the same time, as they initially were 2 views of the same data
module mod_cosc5
   implicit none
   real(8), dimension(:), allocatable :: sc5
   integer, dimension(:), allocatable :: sc5int
   contains
   subroutine allocate_cosc5(n)
      integer, intent(in) :: n
      allocate(sc5(n))
      allocate(sc5int(n))
   end subroutine allocate_cosc5
end module mod_cosc5

module mod_cosc6
   implicit none
   real(8), dimension(:), allocatable :: sc6
   contains
   subroutine allocate_cosc6(n)
      integer, intent(in) :: n
      allocate(sc6(n))
   end subroutine allocate_cosc6
end module mod_cosc6

module mod_cosc7
   implicit none
   real(8), dimension(:), allocatable :: sc7
   contains
   subroutine allocate_cosc7(n)
      integer, intent(in) :: n
      allocate(sc7(n))
   end subroutine allocate_cosc7
end module mod_cosc7

module mod_cosc8
   implicit none
   real(8), dimension(:), allocatable :: sc8
   contains
   subroutine allocate_cosc8(n)
      integer, intent(in) :: n
      allocate(sc8(n))
   end subroutine allocate_cosc8
end module mod_cosc8

module mod_cosc9
   implicit none
   real(8), dimension(:), allocatable :: sc9
   contains
   subroutine allocate_cosc9(n)
      integer, intent(in) :: n
      allocate(sc9(n))
   end subroutine allocate_cosc9
end module mod_cosc9

module mod_cosc10
   implicit none
   real(8), dimension(:), allocatable :: sc10
   contains
   subroutine allocate_cosc10(amax)
      integer, intent(in) :: amax
      allocate(sc10(amax))
   end subroutine allocate_cosc10
end module mod_cosc10

module mod_cosc11
   implicit none
   real(8), dimension(:), allocatable :: sc11
   contains
   subroutine allocate_cosc11(aaux)
      integer, intent(in) :: aaux
      allocate(sc11(aaux))
   end subroutine allocate_cosc11
end module mod_cosc11

module mod_coeig
   ! module for ba3d1p basis routine
   implicit none
   real(8), dimension(:,:), allocatable :: c0, c1, c2
   contains
   subroutine allocate_coeig()
      allocate(c0(4,4)) ; allocate(c1(3,3)) ; allocate(c2(2,2))
   end subroutine allocate_coeig
end module mod_coeig

module mod_coeig2
   ! * module for ba3p2s basis routine
   ! matrices for transformation between atomic and molecular BF functions
   implicit none
   real(8), dimension(:,:), allocatable :: t12, t32
   contains
   subroutine allocate_coeig2()
      allocate(t12(5,5)) ; allocate(t32(3,3))
   end subroutine allocate_coeig2
end module mod_coeig2

module mod_cokaux
   implicit none
   integer, allocatable :: naux
   contains
   subroutine allocate_cokaux(anaux)
      integer, intent(in) :: anaux
      allocate(naux)
      naux = anaux
   end subroutine allocate_cokaux
end module mod_cokaux

module mod_cotble
   ! *     npnt:     max. number of pointer
   ! *     jttble:   array containing pointer to records in s-matrix   
   implicit none
   integer, dimension(:), allocatable :: jttble
   integer, allocatable               :: npnt
   contains
   subroutine allocate_cotble(afact)
      integer, intent(in) :: afact
      allocate(jttble(afact)) ; allocate(npnt)
      npnt = afact
   end subroutine allocate_cotble
end module mod_cotble


module mod_coqvec
   ! *  vectors dimensioned nphoto*nch
   !       q

   ! *  variables in common block /coqvec/
   ! *     mxphot        maximum column dimension of q matrix
   ! *     nphoto        actual column dimension of q matrix
   ! *     q             accumulated gamma2 inhomogeneous propagator
   ! *                   only calculated if photof = .true.
   ! *                   this is stored as a column vector for each separate
   ! *                   ground-state wavefunction
   
   ! *     variables in common block /coqvec/
   ! *     mxphot        maximum dimension of q matrix
   ! *     nphoto        actual collumn dimension of q matrix
   ! *     q             accumulated overlap matrix with ground state
   ! *                   only calculated if photof = .true.
   ! *                   this is stored with each wavefunction as a column vector

   ! *     variables in common block /coqvec/
   ! *     mxphot        maximum column dimension of q matrix
   ! *     nphoto        actual collumn dimension of q matrix
   ! *     q             accumulated overlap matrix with ground state
   ! *                   only calculated if photof = .true.
   ! *                   this is stored with each wavefunction as a column vector
   
   ! *     variables in common block /coqvec/
   ! *     mxphot        maximum column dimension of q matrix
   ! *     nphoto        actual column dimension of q matrix (gamma2)
   ! *     q             accumulated overlap matrix with ground state
   ! *                   only calculated if photof = .true.
   ! *                   this is stored with each wavefunction as a column vector

   ! *     variables in common block /coqvec/
   ! *     mxphot        maximum column dimension of q matrix
   ! *     nphoto        actual column dimension of q matrix
   ! *     q             accumulated overlap matrix with ground state
   ! *                   only calculated if photof = .true.
   ! *                   this is stored with each wavefunction as a column vector

   ! *     variables in common block /coqvec/
   ! *     mxphot        maximum dimension of q matrix
   ! *     nphoto        actual collumn dimension of q matrix
   ! *     q             accumulated overlap matrix with ground state
   ! *                   only calculated if photof = .true.
   ! *                   this is stored with each wavefunction as a column vector


   ! mxphot:   maximum number of initial states (maximum columns in q)
   ! nphoto:   actual number of initial states used
   ! q:        accumulated overlap matrix with ground state
   !           only calculated if photof = .true.
   implicit none
   real(8), dimension(:), allocatable :: q
   integer, allocatable               :: mxphot
   integer, allocatable               :: nphoto
   contains
   subroutine allocate_coqvec(aqmax, amxpho, anphot)
      integer, intent(in) :: aqmax
      integer, intent(in) :: amxpho
      integer, intent(in) :: anphot
      allocate(q(aqmax)) ; allocate(mxphot) ; allocate(nphoto)
      mxphot = amxpho
      nphoto = anphot
   end subroutine allocate_coqvec
end module mod_coqvec

module mod_coqvec2
   implicit none
   real(8), dimension(:), allocatable :: q2
   contains
   subroutine allocate_coqvec2(aq2)
      integer, intent(in) :: aq2
      allocate(q2(aq2))
   end subroutine allocate_coqvec2
end module mod_coqvec2

module mod_codim
   implicit none
   save
   integer, allocatable :: mairy, mmax, mbig
   contains
   subroutine allocate_codim(aairy, amax, abig)
      integer, intent(in) :: aairy
      integer, intent(in) :: amax
      integer, intent(in) :: abig
      allocate(mairy) ; allocate(mmax) ; allocate(mbig)
      mairy = aairy
      mmax = amax
      mbig = abig
   end subroutine allocate_codim
end module mod_codim

! cotwo
! stores data related to systems with 2 molecules
!    numj:     number of j1-j2 values
!    nj1j2:    specific j1-j2 values (up to a maximum of 50)
!              N.B. this dimension is set here
module mod_two
   integer :: numj
   integer :: nj1j2(50)
end module mod_two

! coopti
!    optifl:      flag, signals if the calculation is an optimization
module mod_opti
   logical :: optifl
end module mod_opti

module mod_comxbs
   implicit none
   save
   integer, allocatable :: maxbas  ! maximum number of allowed basis types 
   ! THIS SHOULD BE INCREASED AS BASIS ROUTINES ARE ADDED
   contains
   subroutine allocate_comxbs(n)
      integer, intent(in) :: n
      allocate(maxbas) ; maxbas = n
   end subroutine allocate_comxbs
end module mod_comxbs

module mod_comxm
   implicit none
   save
   integer, allocatable :: ncache, mxmblk
   contains
   subroutine allocate_comxm()
      allocate(ncache)  ! cache size in words
      allocate(mxmblk)
      !  determine cache and block sizes for matrix multiply
      ncache=4096
#if defined(HIB_UNIX_HP)
      ncache=16000
#endif
      mxmblk=64
   end subroutine allocate_comxm
end module mod_comxm

module mod_cosys
   implicit none
   save
   integer :: lencod
   character(len=8), dimension(:), allocatable :: scod
   contains
   subroutine allocate_cosys(alencod)
      integer, intent(in) :: alencod
      lencod = alencod
      allocate(scod(lencod))
   end subroutine allocate_cosys
end module mod_cosys

module mod_cosyr
   implicit none
   save
   integer :: maxpar
   character(len=8), dimension(:), allocatable :: rcod
   contains
   subroutine allocate_cosyr(amaxpar)
      integer, intent(in) :: amaxpar
      maxpar = amaxpar
      allocate(rcod(maxpar))
   end subroutine allocate_cosyr
end module mod_cosyr

module mod_cosysi
   implicit none
   save
   integer :: maxpar
   integer :: nscode
   integer :: isicod
   integer, dimension(:), allocatable, target :: ispar
   contains
   subroutine allocate_cosysi(amaxpar)
      integer, intent(in) :: amaxpar
      maxpar = amaxpar
      allocate(ispar(amaxpar))
   end subroutine allocate_cosysi
   subroutine convert_ispar_to_mat(nlines,ncols,istart,new)
      integer, intent(in) :: nlines, ncols, istart
      integer, dimension(nlines,ncols) :: new
      integer i, j
      do i=1,nlines
         do j=1,ncols
           new(i,j) = ispar(istart-1+i+(j-1)*nlines)
         enddo
      enddo
   end subroutine convert_ispar_to_mat
end module mod_cosysi


module mod_cosysl
   implicit none
   save
   integer :: maxpar
   integer :: islcod
   integer, dimension(:), allocatable, target :: lspar
   contains
   subroutine allocate_cosysl(amaxpar)
      integer, intent(in) :: amaxpar
      maxpar = amaxpar
      islcod = 0
      allocate(lspar(amaxpar))
   end subroutine allocate_cosysl
end module mod_cosysl

module mod_cosysr
   implicit none
   save
   integer :: maxpar
   integer :: isrcod
   integer :: junkr
   real(8), dimension(:), allocatable, target :: rspar
   contains
   subroutine allocate_cosysr(amaxpar)
      integer, intent(in) :: amaxpar
      maxpar = amaxpar
      allocate(rspar(maxpar))
   end subroutine allocate_cosysr
   subroutine convert_rspar_to_mat(nlines,ncols,new)
      integer, intent(in) :: nlines, ncols
      real(8), dimension(nlines,ncols) :: new
      integer i, j
      do i=1,nlines
         do j=1,ncols
           new(i,j) = rspar(i+(j-1)*nlines)
         enddo
      enddo
   end subroutine convert_rspar_to_mat
end module mod_cosysr


module mod_par
   use mod_hiparcst, only: LPAR_COUNT, IPAR_COUNT, RPAR_COUNT
   implicit none
   save

   !  fcod stores logical flags (length = lcode)
   ! fcod = Flags CODes : stores the name of system independent parameters of type logical
   character(len=8), parameter :: fcod(LPAR_COUNT) = [ &
      'AIRYFL', &
      'BASTST', &
      'BATCH ', &
      'CHLIST', &
      'CSFLAG', &
      'FLAGHF', &
      'FLAGSU', &
      'IHOMO ', &
      'IPOS  ', &
      'LOGDFL', &
      'NOPRIN', &
      'NUCROS', &
      'PHOTOF', &
      'PRAIRY', &
      'PRLOGD', &
      'PRPART', &
      'PRSMAT', &
      'PRT2  ', &
      'PRXSEC', &
      'READPT', &
      'RSFLAG', &
      'T2TEST', &
      'TWOMOL', &
      'WAVEFL', &
      'WRPART', &
      'WRSMAT', &
      'WRXSEC', &
      'BOUNDC']

   ! pcod = Parameters CODes : stores the name of system independent parameters of type integer and real
   character(len=8) :: pcod(IPAR_COUNT+RPAR_COUNT) = [ &
      ! parameters of type integer
      'JTOT1   ', &  
      'JTOT2   ', &
      'JTOTD   ', &
      'JLPAR   ', &
      'NERG    ', &
      'NUMAX   ', &
      'NUMIN   ', &
      'NUD     ', &
      'LSCREEN ', &
      'IPRINT  ', &
      ! parameters of type real
      'FSTFAC  ', &  
      'RINCR   ', &
      'RCUT    ', &
      'RENDAI  ', &
      'RENDLD  ', &
      'RSTART  ', &
      'SPAC    ', &
      'TOLAI   ', &
      'XMU     ']


   logical, dimension(:), allocatable, target :: lpar
   !  variables in common block /colpar/
   !    airyfl:      if .true., then airy propagation will take place
   !    prairy:      if .true., then step-by-step information is printed out in
   !                 airy propagation
   !    bastst:      if .true., then execution terminates after the first call
   !                 to basis
   !    batch:       if .true., then the job is run as a batch job
   !                 if .false., then the job is assumed to be interactive
   !    chlist:      if .true., then the channel quantum numbers and energies are
   !                 printed out at at each total-j
   !                 if .false., then  this is done only at first total-j
   !    csflag:      if .true., then coupled-states calculation is desired
   !                 if .false., then close-coupled calculation is desired
   !    flaghf:      if .true., then the system has even multiplicity
   !                 (half-integer total angular momentum)
   !    flagsu:      if .true., then the problem is assumed to a molecule
   !                 scattering off a surface, in which case the diagonal
   !                 elements of the transition probabilities are equal to the
   !                 modulus squared of the s-matrix (not t-matrix elements)
   !    ihomo:       if .true., then the molecule is assumed to be homonuclear
   !    ipos:        if .true., then printout is suited for a 132-position printe
   !                 if .false., then printout is suited for a  80 -position
   !                 printer
   !    logdfl:      if .true., then logd propagation will take place
   !    prlogd:       if .true., then the lower triangle of log-derivative matrix
   !                 is printed at the end of the logd and the end of the airy
   !                 integrations
   !    noprin:      if .true., then most printing is suppressed
   !    prpart:       if .true, then the full matrix of partial cross sections
   !                 (summed over m-states) is printed
   !    readpt:      if .true., then potential parameters are expected
   !    rsflag:      if .true., then calculation is to be restarted
   !                 a check will be made to see if all required files
   !                 are present:  these may include
   !                    trstrt, tmp10, tmp11, xsecn (or tmpxn), smatn,
   !                    psecn, tmp35, ...
   !    prsmat:       if .true., then the upper triangle of real and imaginary
   !                 parts of s-matrix are printed
   !    t2test:      if .true., then the first two columns of the square modulus
   !                 of the t-matrix are printed
   !    prt2:      if .true., then the upper triangle of square modulus of
   !                 t-matrix is printed
   !    twomol:      .true. if molecule-molecule collision
   !    wrsmat:       if .true., and nnout is > 0, then those s-matrix elements
   !                 for which both the initial and final rotational quantum
   !                 numbers are in the array jout (input line 12) are written to
   !                 files smat1, smat2, ...
   !                 if nnout < 0, then each column of the s-matrix whose initial
   !                 index is in the array jout is written to files smat1, smat2,
   !    wrpart:      if .true., then input data and the matrix of partial cross
   !                 sections (summed over m-states) is written to file pxsec
   !    wrxsec:      if .true., then some input data and the full matrix of
   !                 integral cross sections ((summed over m-states and summed
   !                 from jtot1 to jtot2) is written to file xsec1, xsec2, ....
   !    prxsec:      if .true., then the full matrix of integral cross sections
   !                 ((summed over m-states and summed from jtot1 to jtot2) is
   !                 printed
   !


   !    airyfl:      if .true., then airy propagation will take place
   logical, pointer :: airyfl

   !    prairy:      if .true., then step-by-step information is printed out in
   !                 airy propagation
   logical, pointer :: prairy

   !    bastst:      if .true., then execution terminates after the first call
   !                 to basis
   logical, pointer :: bastst

   logical, pointer :: batch

   !     chlist:  if .true., then the channel quantum numbers and energies are
   !              printed out at at each total-j
   !              if .false., then  this is done only at first total-j
   logical, pointer :: chlist

   !    csflag:   if .true., then coupled-states calculation is desired
   !              if .false., then close-coupled calculation is desired
   logical, pointer :: csflag

   !    flaghf:   if .true., then the system has even multiplicity (half-integer
   !              total angular momentum)
   logical, pointer :: flaghf

   !    flagsu:   if .true., then the problem is assumed to a molecule scattering
   !              of a surface, in which case the diagonal elements of the
   !              transition probabilities are equal to the modulus squared of
   !              the s-matrix (not t-matrix elements)
   logical, pointer :: flagsu

   !    ihomo:    if .true., then the molecule is assumed to be homonuclear
   logical, pointer :: ihomo

   !     ipos:    if .true., then printout is suited for a 132-position printer
   !              if .false., then printout is suited for a  80 -position printer
   logical, pointer :: ipos

   !    logdfl:      if .true., then logd propagation will take place
   logical, pointer :: logdfl

   !     prlogd:  if .true., then the lower triangle of log-derivative matrix
   !              is printed at the end of the logd and the end of the airy
   !              integrations
   logical, pointer :: prlogd

   !     noprin:  if .true., then most printing is suppressed
   logical, pointer :: noprin

   !     prpart:  if .true, then the full matrix of partial cross sections (summe
   !              over m-states) is printed
   logical, pointer :: prpart

   !    readpt:      if .true., then potential parameters are expected
   logical, pointer :: readpt

   !    rsflag:      if .true., then calculation is to be restarted
   !                 a check will be made to see if all required files
   !                 are present:  these may include
   !                    trstrt, tmp10, tmp11, xsecn (or tmpxn), smatn,
   !                    psecn, tmp35, ...  
   logical, pointer :: rsflag

   !     prsmat:  if .true., then the upper triangle of real and imaginary parts
   !              of s-matrix are printed
   logical, pointer :: prsmat

   !     t2test:  if .true., then the first two columns of the square modulus of
   !              the t-matrix are printed
   logical, pointer :: t2test

   !     prt2:    if .true., then the upper triangle of square modulus of t-matri
   !              is printed
   logical, pointer :: prt2

   !    twomol:   if .true., then molecule-molecule collision is assumed
   logical, pointer :: twomol

   !     wrsmat:  if .true., and nnout is > 0, then those s-matrix elements
   !              for which both the initial and final rotational quantum
   !              numbers are in the array jout (input line 12) are written to
   !              files smat1, smat2, ...
   !              if nnout < 0, then each column of the s-matrix whose initial
   !              index is in the array jout is written to files smat1, smat2, ..
   logical, pointer :: wrsmat

   !     wrpart:  if .true., then input data and the matrix of partial cross
   !              sections (summed over m-states) is written to file pxsec
   logical, pointer :: wrpart

   !     wrxsec:  if .true., then some input data and the full matrix of integral
   !              cross sections ((summed over m-states and summed from
   !              jtot1 to jtot2) is written to file xsec1, xsec2, ....
   logical, pointer :: wrxsec

   !     prxsec:  if .true., then the full matrix of integral cross sections
   !              ((summed over m-states and summed from jtot1 to jtot2) is print
   logical, pointer :: prxsec

   !     nucros:  parameter to control how CS integral cross sections are
   !              computed
   logical, pointer :: nucros

   !     photof:  if .true. then photodissociation calculation
   logical, pointer :: photof

   !     wavefl:  if .true. then information is written to calculate,
   !              subsequently, wavefunctions, fluxes, and adiabatic energies
   logical, pointer :: wavefl

   !     boundc:  if .true. then susan gregurick's bound state calculation
   !              is implemented
   logical, pointer :: boundc

   ! ipar (integer parameters)

   integer, dimension(:), allocatable, target :: ipar


   integer, pointer :: jtot1
   integer, pointer :: jtot2
   integer, pointer :: jtotd
   integer, pointer :: jlpar
   integer, pointer :: nerg
   integer, pointer :: numax
   integer, pointer :: numin
   integer, pointer :: nud
   
   ! lscreen is the number of lines available on your terminal screen
   integer, pointer :: lscreen

   ! iprint controls degree of print output in some routines
   !     iprint=-1 (no print); iprint=0 (min print); iprint=1 some print, etc
   integer, pointer :: iprint

   ! rpar (real number parameters)

   real(8), dimension(:), allocatable, target :: rpar

   ! parameters for scattering mode
   real(8), pointer :: scat_fstfac
   real(8), pointer :: scat_rincr
   real(8), pointer :: scat_rcut
   real(8), pointer :: scat_rendai
   real(8), pointer :: scat_rendld
   real(8), pointer :: scat_rstart
   real(8), pointer :: scat_spac
   real(8), pointer :: scat_tolai

   ! parameters for bound state mode
   real(8), pointer :: bound_r1
   real(8), pointer :: bound_r2
   real(8), pointer :: bound_c
   real(8), pointer :: bound_spac
   real(8), pointer :: bound_delr
   real(8), pointer :: bound_hsimp
   real(8), pointer :: bound_eigmin
   real(8), pointer :: bound_tolai

   ! parameters common to scattering mode and bound state mode
   real(8), pointer :: xmu

   contains
   subroutine allocate_par()
      use mod_hiparcst, only: LPAR_COUNT, IPAR_COUNT, RPAR_COUNT
      use lpar_enum
      use ipar_enum
      use rpar_enum
      implicit none
      integer :: num_lpar = LPAR_COUNT
      integer :: num_ipar = IPAR_COUNT
      integer :: num_rpar = RPAR_COUNT

      allocate(lpar(num_lpar))

      airyfl => lpar(LPAR_AIRYFL)
      prairy => lpar(LPAR_PRAIRY)
      bastst => lpar(LPAR_BASTST)
      batch => lpar(LPAR_BATCH)
      chlist => lpar(LPAR_CHLIST)
      csflag => lpar(LPAR_CSFLAG)
      flaghf => lpar(LPAR_FLAGHF)
      flagsu => lpar(LPAR_FLAGSU)
      ihomo => lpar(LPAR_IHOMO)
      ipos => lpar(LPAR_IPOS)
      logdfl => lpar(LPAR_LOGDFL)
      prlogd => lpar(LPAR_PRLOGD)
      noprin => lpar(LPAR_NOPRIN)
      prpart => lpar(LPAR_PRPART)
      readpt => lpar(LPAR_READPT)
      rsflag => lpar(LPAR_RSFLAG)
      prsmat => lpar(LPAR_PRSMAT)
      t2test => lpar(LPAR_T2TEST)
      prt2 => lpar(LPAR_PRT2)
      twomol => lpar(LPAR_TWOMOL)
      wrsmat => lpar(LPAR_WRSMAT)
      wrpart => lpar(LPAR_WRPART)
      wrxsec => lpar(LPAR_WRXSEC)
      prxsec => lpar(LPAR_PRXSEC)
      nucros => lpar(LPAR_NUCROS)
      photof => lpar(LPAR_PHOTOF)
      wavefl => lpar(LPAR_WAVEFL)
      boundc => lpar(LPAR_BOUNDC)

      allocate(ipar(num_ipar))

      jtot1 => ipar(IPAR_JTOT1)
      jtot2 => ipar(IPAR_JTOT2)
      jtotd => ipar(IPAR_JTOTD)
      jlpar => ipar(IPAR_JLPAR)
      nerg => ipar(IPAR_NERG)
      numax => ipar(IPAR_NUMAX)
      numin => ipar(IPAR_NUMIN)
      nud => ipar(IPAR_NUD)
      lscreen => ipar(IPAR_LSCREEN)
      iprint => ipar(IPAR_IPRINT)

      allocate(rpar(num_rpar))

      scat_fstfac => rpar(RPAR_SCAT_FSTFAC)
      scat_rincr => rpar(RPAR_SCAT_RINCR)
      scat_rcut => rpar(RPAR_SCAT_RCUT)
      scat_rendai => rpar(RPAR_SCAT_RENDAI)
      scat_rendld => rpar(RPAR_SCAT_RENDLD)
      scat_rstart => rpar(RPAR_SCAT_RSTART)
      scat_spac => rpar(RPAR_SCAT_SPAC)
      scat_tolai => rpar(RPAR_SCAT_TOLAI)

      bound_r1 => rpar(RPAR_BOUND_R1)
      bound_r2 => rpar(RPAR_BOUND_R2)
      bound_c => rpar(RPAR_BOUND_C)
      bound_spac => rpar(RPAR_BOUND_SPAC)
      bound_delr => rpar(RPAR_BOUND_DELR)
      bound_hsimp => rpar(RPAR_BOUND_HSIMP)
      bound_eigmin => rpar(RPAR_BOUND_EIGMIN)
      bound_tolai => rpar(RPAR_BOUND_TOLAI)


      xmu => rpar(RPAR_XMU)

   end subroutine allocate_par

   subroutine set_param_names(boundc, param_names, param_names_size)
      !  subroutine to change param_names's for bound state or scattering
      use mod_hiparcst, only: IPAR_COUNT
      use rpar_enum
      implicit none
      logical, intent(in) :: boundc
      integer, intent(in) :: param_names_size  ! size of param_names array
      character*8, intent(out) :: param_names(param_names_size)  ! array containing the name of each parameter (old name: pcod)
      if (boundc) then
        param_names(IPAR_COUNT + RPAR_BOUND_R1)      = 'R1'
        param_names(IPAR_COUNT + RPAR_BOUND_R2)      = 'R2'
        param_names(IPAR_COUNT + RPAR_BOUND_C)       = 'C' 
        param_names(IPAR_COUNT + RPAR_BOUND_SPAC)    = 'SPAC' 
        param_names(IPAR_COUNT + RPAR_BOUND_DELR)    = 'DELR' 
        param_names(IPAR_COUNT + RPAR_BOUND_HSIMP)   = 'HSIMP' 
        param_names(IPAR_COUNT + RPAR_BOUND_EIGMIN)  = 'EIGMIN' 
        param_names(IPAR_COUNT + RPAR_BOUND_TOLAI)   = 'TOLAI' 
      else
        param_names(IPAR_COUNT + RPAR_SCAT_FSTFAC)  = 'FSTFAC'
        param_names(IPAR_COUNT + RPAR_SCAT_RINCR)   = 'RINCR'
        param_names(IPAR_COUNT + RPAR_SCAT_RCUT)    = 'RCUT'
        param_names(IPAR_COUNT + RPAR_SCAT_RENDAI)  = 'RENDAI'
        param_names(IPAR_COUNT + RPAR_SCAT_RENDLD)  = 'RENDLD'
        param_names(IPAR_COUNT + RPAR_SCAT_RSTART)  = 'RSTART'
        param_names(IPAR_COUNT + RPAR_SCAT_SPAC)    = 'SPAC'
        param_names(IPAR_COUNT + RPAR_SCAT_TOLAI)   = 'TOLAI'
      endif
      return
   end


end module mod_par


! used to be common block coselb
module mod_selb
   integer :: ibasty  ! base type
end module mod_selb


! used to be common block coered
!    ered:      collision energy in atomic units (hartrees)
!    rmu:       collision reduced mass in atomic units (mass of electron = 1)
module mod_ered
real(8) :: ered
real(8) :: rmu
end module mod_ered

! used to be common block coskip
!   nskip  for a homonuclear molecule lamda is running in steps of nskip=2
!          for a heteronuclear molecule nskip=1
!
!   iskip   same as nskip, used for consistency check
module mod_skip
integer :: nskip
integer :: iskip
end module mod_skip

! used to be common block covib
module mod_vib
use mod_parbas, only: maxvib
integer :: nvibs = -1
integer :: ivibs(maxvib) = -1
integer :: nvibp = -1
integer :: ivibp(maxvib) = -1
end module mod_vib

! used to be common block cophot
!     photof        true if photodissociation calculation
!                   false if scattering calculation
!     wavefn        true if G(a,b) transformation matrices are saved
!                   to be used later in computing the wavefunction
module mod_phot
logical :: photof
logical :: wavefn
logical :: boundf
logical :: writs
end module mod_phot

! used to be common block cospbf
module mod_spbf
   integer :: lnbufs
   integer :: lnbufl
   integer :: nbuf
   integer :: maxlsp
   integer :: maxllb
   integer :: ihibuf
   integer :: igjtp
end module mod_spbf

! used to be common block comom
module mod_mom
  real(8) :: spin
  real(8) :: xj1
  real(8) :: xj2
  integer :: j1
  integer :: in1
  integer :: j2
  integer :: in2
  integer :: maxjt
  integer :: maxjot
  integer :: nwaves
  integer :: jfsts
  integer :: jlparf
  integer :: jlpars
  integer :: njmax
  integer :: j1min
  integer :: j2max
end module mod_mom

! used to be common block cofile
module mod_file
  character(40) :: input
  character(40) :: output
  character(40) :: jobnam
  character(40) :: savfil
end module mod_file

! used to be common block cosurf
!    flagsu:    if .true., then molecule-surface collisons
!                 this variable is set equal to flagsu, it is held in a
!                 separate common block for compatability with subroutines
!                 smatop, soutpt, and xwrite
!                 if .true., then the problem is assumed to a molecule
!                 scattering off a surface, in which case the diagonal
!                 elements of the transition probabilities are equal to the
!                 modulus squared of the s-matrix (not t-matrix elements)

module mod_surf
  logical :: flagsu
end module mod_surf


! used to be common block cojtot
module mod_jtot
  integer :: jjtot
  integer :: jjlpar
end module mod_jtot

! used to be common block coja
module mod_ja
  integer :: jja(9)
end module mod_ja

! used to be common block coel
module mod_el
  integer :: ll(9)
end module mod_el

! used to be common block cosavi and cosavr
module mod_sav
  use mod_hiparcst, only: IPAR_COUNT, RPAR_COUNT

  integer :: iipar
  integer :: ixpar(IPAR_COUNT)

  integer :: irpar
  real(8) :: rxpar(RPAR_COUNT)
end module mod_sav

! used to be common block copmat
module mod_pmat
  real(8) :: rtmn
  real(8) :: rtmx
  integer :: iflag
end module mod_pmat

! used to be common block cputim
module mod_cputim
  real(8) :: cpuld
  real(8) :: cpuai
  real(8) :: cpupot
  real(8) :: cpusmt
  real(8) :: cpupht
end module mod_cputim

 ! all the commons blocks from hiiolib_f.F90:
    !!   common/cdbf/ ldbuf,libuf,ibfil,ibrec,ibof,ibstat,idbuf(llbuf)

 ! all the commons blocks from himain.t:
    !!   common /comom/  xmom(3), imom(13)
    !!   common /cosout/ nnout, jout(kout)
    !!   common /coiout/ niout, indout(kout)
    !!   common /cocent/ cent(kmax)
    !!   common /coeint/ eint(kmax)
    !!   common /coj12/ j12(kmax)
    !!   common /coj12p/ j12pk(kmax)
    !!   common /covvl/  vvl(klammx)
    !!   common /cofact/ si(kfact)
    !!  common /coener/ energ(ken)
    !!   common /clseg/  lseg,intrel,lchar
    !!   common /cobuf/  lbuf,ibuf(1024)
    !!   common /cofil/  nfl,iofbuf,maxrec(60),iofrec(60),nwrec
    !!   common /conlam/ nlam, nlammx, lamnum(klammx)
    !!   common /coatpi/ narray, isiz(krotmx)
    !!   common /coatp3/ isizh(krotmx)
    !!   common /coatpr/ c(krotmx)
    !!   common /coatp1/ ctemp(krotmx)
    !!   common /coatp2/ chold(krotmx)
    !!   common /coz/ z(kmax,kmax)
    !!   common /cow/ w(kmax,kmax)
    !!   common /cozmat/ zmat(kmax,kmax)
    !!   common /coamat/ amat(kmax,kmax)
    !!   common /cobmat/ bmat(kairy,kairy)
    !!   common /cotq1/ tq1(kmax,kmax)
    !!   common /cotq2/ tq2(kmax,kmax)
    !!   common /cotq3/ tq3(kmax,kmax)
    !!   common /cojq/ jq(kmax)
    !!   common /colq/ lq(kmax)
    !!   common /coinq/ inq(kmax)
    !!   common /cojhld/ jhold(kmax)
    !!   common /coehld/ ehold(kmax)
    !!   common /coinhl/ inhold(kmax)
    !!   common /coisc1/ isc1(kmax)
    !!   common /coisc2/ isc2(kmax)
    !!   common /coisc3/ isc3(kmax)
    !!   common /coisc4/ isc4(kmax)
    !!   common /coisc5/ isc5(kmax)
    !!   common /coisc6/ isc6(kmax)
    !!   common /coisc7/ isc7(kmax)
    !!   common /coisc8/ isc8(kmax)
    !!   common /coisc9/ isc9(kmax)
    !!   common /coisc10/ isc10(kmax)
    !!   common /coisc11/ isc11(kmax)
    !!   common /coisc12/ isc12(kmax)
    !!   common /colsc1/ lsc1(kmax)
    !!   common /cosc1/ sc1(kmax)
    !!   common /cosc2/ sc2(kmax)
    !!   common /cosc3/ sc3(kmax)
    !!   common /cosc4/ sc4(kmax)
    !!   common /cosc5/ sc5(kmax)
    !!   common /cosc6/ sc6(kmax)
    !!   common /cosc7/ sc7(kmax)
    !!   common /cosc8/ sc8(kmax)
    !!   common /cosc9/ sc9(kmax)
    !!   common /cosc10/ sc10(kmax)
    !!   common /coeig2/  t12(5,5), t32(3,3)
    !!   common /coeig/  c0(4,4), c1(3,3), c2(2,2)
    !!   common /cosc11/ sc11(kaux3)
    !!   common /cosc11/ sc11(kaux)
    !!   common /cokaux/ naux
    !!   common /cotble/ npnt, jttble(kfact)
    !!   common /coqvec/ mxphot, nphoto, q(kqmax)
    !!   common /coqvec2/ q2(kq2)
    !!   common /codim/ mairy,mmax,mbig
    !!   common /comxbs/ maxbas
    !!   common /comxm/ ncache, mxmblk