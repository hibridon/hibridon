!     This file contains various modules that replace common blocks in
!     hibridon.
!
!     Physical constants
module constants
implicit none
real(8), parameter :: tolerance = 1.0d-10
real(8), parameter :: zero = 0.0d0
real(8), parameter :: half = 0.5d0
real(8), parameter :: one  = 1.0d0
real(8), parameter :: two  = 2.0d0
real(8), parameter :: four = 4.0d0
real(8), parameter :: rad = 57.29577951308232d0  ! number of angular degrees in a radian

real(8), parameter :: pi=dacos(-1d0)
real(8), parameter :: s4pi = sqrt( 4.d0 * pi )  ! normalization factor for isotropic potential
!     Below are constants used in Hibridon 4.4
!
!     All physical constants below are adapted from 2010 values from
!     http://physics.nist.gov/cuu/Constants/index.html
!
!     econv is the conversion from hartree to cm-1
!     hartree to wavenumber
!     The value of econv is the rydberg constant in m-1 devided by 50
!$$$      real(8), parameter :: econv=10973731.568539/50d0
!     xmconv converts mass in atom units (electron mass) to atom mass
!     unit (1/12 of the mass of 12C)
!$$$      real(8), parameter :: xmconv=5.4857990946d-4
!     ang2c converts bohr^2 to angstrom^2
!     The value of ang2 is (10^10 Bohr Radius)^2
!$$$      real(8), parameter :: ang2c=0.52917721092d0**2
real(8), parameter :: econv = 219474.6d0
real(8), parameter :: xmconv = 5.485930d-4
real(8), parameter :: ang2c = 0.280002846d0
real(8), parameter :: amu_to_emu = 1822.88848477d0  ! amu to electrom mass (1.0 / xmconv)

end module constants

! file units used in hibridon
module funit
implicit none

   enum, bind( C )
   enumerator ::  &
      FUNIT_CS             =  1, &   ! cross sections 
      FUNIT_ICS            =  1, &   ! <job-name>.ics
      FUNIT_EADIAB         =  2, &   ! <job-name>.eadiab
      FUNIT_TENS_OUTPUT    =  2, &   ! <job-name>.tcs or <job-name>.dch or <job-name>.dcga
      FUNIT_MCS            =  2, &   ! <job-name>.mcs
      FUNIT_PSI            =  2, &   ! <job-name>.psi
      FUNIT_SAV            =  3, &   ! <job-name>.sav
      FUNIT_XSC            =  3, &   ! <job-name>.xsc
      FUNIT_FLX            =  3, &   ! <job-name>.flx
      FUNIT_TCB            =  4, &   ! <job-name>.tcb (result of tenxsc)
      FUNIT_STDOUT         =  6, &   ! standard output
      FUNIT_INP            =  8, &   ! input file (*.inp)
      FUNIT_OUT            =  9, &   ! output file (*.out)
      FUNIT_TRANS_MAT      = 10, &   ! temporary storage for transformation matrices (t<unit><pid>.tmp)
      FUNIT_QUAD_MAT       = 11, &   ! temporary storage for quadrature matrices (t<unit><pid>.tmp)
      FUNIT_CHANNEL_PARAMS = 12, &   ! temporary storage for channel parameters (t<unit><pid>.tmp)
      FUNIT_WFU            = 22, &   ! <job-name>.wfu
      FUNIT_PCS_START      = 25, &   ! <job-name><energy_index>.pcs files (partial cross sections; one for each energy)
      FUNIT_APCS_START     = 35, &   ! accumulation of cs partial cross sections at each projection index; one for each energy. (t<unit><pid>.tmp)
      FUNIT_SMT_START      = 45, &   ! <job-name><energy_index>.smt storage of selected s-matrix elements; one for each energy
      FUNIT_ICS_START      = 70      ! <job-name><energy_index>.ics files (integral cross sections; one for each energy)
   end enum


end module funit

module mod_fileid
implicit none

   enum, bind( C )
   enumerator ::  &
      FILEID_SAV           =  1, &   ! <job-name>.sav
      FILEID_TMP           =  1      ! <job-name>.tmp
   end enum


end module mod_fileid

module mod_hitypes
   implicit none
   private

   type, public :: bqs_type  ! base quantum states
      ! stores an array of base quantum states
      ! rotational angular momenta, orbital angular momenta, and
      ! additional quantum index for each channel
      ! if the calculation involves the collisions of two diatomic
      ! molecules, the bqs%jq = j1 + 10000 j2, where j1 and j2 are the
      ! rotational quantum numbers of each molecule      

      ! the vector jq will hold the rotational quantum numbers in the
      ! packed bas
      !
      ! the vector lq will hold the orbital angular momenta of each
      ! channel in the packed basis
      !
      ! the vector inq will hold the symmetry indices in the packed basis
      ! first sum over the unpacked states
      
      integer, allocatable :: jq(:)  ! combined rotational quantum number. This integer stores an encoded representation of the set of rotational quantum numbers: for most bases, it only contains the value of j (rotational quantum number), but in case of 2 molecule collisions it encodes both j1 (rotational quantum number of modlecule 1) and j2 (rotational quantum number of molecule 2). The base itself is responsible for encoding and decoding this number.
      integer, allocatable :: lq(:)  ! orbital angular momentum quantum number
      integer, allocatable :: inq(:)  ! symmetry index
      integer, allocatable :: j12(:)   ! rotational quantum number (only for 2 molecules systems)
      integer :: length  ! number of elements used in all arrays
      integer :: max_length
   contains
      procedure :: init => bqs_type_init
      procedure :: deinit => bqs_type_deinit
   end type bqs_type

   type, public :: rbesself_type  ! ricatti bessel functions
      real(8), allocatable :: fj(:)
      real(8), allocatable :: fpj(:)
      real(8), allocatable :: fn(:)
      real(8), allocatable :: fpn(:)
      integer, allocatable :: length  ! number of elements used in all arrays
   contains
      procedure :: init => rbesself_type_init
      procedure :: deinit => rbesself_type_deinit
   end type rbesself_type


contains
   subroutine bqs_type_init(this, nopen)
      class(bqs_type) :: this
      integer, intent(in) :: nopen

      call bqs_type_deinit(this)
      
      allocate(this%jq(nopen))
      allocate(this%lq(nopen))
      allocate(this%inq(nopen))
      allocate(this%j12(nopen))
      this%length = 0
      this%max_length = nopen
   end subroutine 

   subroutine bqs_type_deinit(this)
      class(bqs_type) :: this
      if (allocated(this%jq)) deallocate(this%jq)
      if (allocated(this%lq)) deallocate(this%lq)
      if (allocated(this%inq)) deallocate(this%inq)
      if (allocated(this%j12)) deallocate(this%j12)
      this%length = 0
      this%max_length = 0
   end subroutine 

   subroutine rbesself_type_init(this, nopen)
      class(rbesself_type) :: this
      integer, intent(in) :: nopen

      call rbesself_type_deinit(this)
      
      allocate(this%fj(nopen))
      allocate(this%fpj(nopen))
      allocate(this%fn(nopen))
      allocate(this%fpn(nopen))
      this%length = 0
   end subroutine 

   subroutine rbesself_type_deinit(this)
      class(rbesself_type) :: this
      if (allocated(this%fj)) deallocate(this%fj)
      if (allocated(this%fpj)) deallocate(this%fpj)
      if (allocated(this%fn)) deallocate(this%fn)
      if (allocated(this%fpn)) deallocate(this%fpn)
      this%length = 0
   end subroutine 

end module mod_hitypes