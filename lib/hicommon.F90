!     This file contains various modules that replace common blocks in
!     hibridon.
!
!     Physical constants
module constants
implicit none
real(8), parameter :: pi=dacos(-1d0)
real(8), parameter :: s4pi = sqrt( 4.d0 * pi )  ! normalization factor for isotropic potential
!     Below are constants used in Hibridon 4.4
real(8), parameter :: econv=219474.6d0, xmconv=5.485930d-4, &
     ang2c=0.280002846d0
!
!     All physical constants below are adapted from 2010 values from
!     http://physics.nist.gov/cuu/Constants/index.html
!
!     econv is the conversion from hartree to cm-1
!     The value of econv is the rydberg constant in m-1 devided by 50
!$$$      real(8), parameter :: econv=10973731.568539/50d0
!     xmconv converts mass in atom units (electron mass) to atom mass
!     unit (1/12 of the mass of 12C)
!$$$      real(8), parameter :: xmconv=5.4857990946d-4
!     ang2c converts bohr^2 to angstrom^2
!     The value of ang2 is (10^10 Bohr Radius)^2
!$$$      real(8), parameter :: ang2c=0.52917721092d0**2
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
      FUNIT_SAV            =  3, &   ! <job-name>.sav
      FUNIT_XSC            =  3, &   ! <job-name>.xsc
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

   type, public :: packed_base_type
      ! the vector iorder will point to the position in the unpacked basis
      ! of each state in the packed basis
      !
      ! the vector jpack will hold the rotational quantum numbers in the
      ! packed bas
      !
      ! the vector lpack will hold the orbital angular momenta of each
      ! channel in the packed basis
      !
      ! the vector epack will hold the channel energies in the packed basis
      !
      ! the vector inpack will hold the symmetry indices in the packed basis
      ! first sum over the unpacked states
      
      integer, allocatable :: iorder(:)
      integer, allocatable :: jpack(:)
      integer, allocatable :: lpack(:)
      real(8), allocatable :: epack(:)
      integer, allocatable :: inpack(:)
      integer, allocatable :: j12pk(:)
      integer, allocatable :: length  ! number of elements used in all arrays
   contains
      procedure :: init => packed_base_type_init
      procedure :: deinit => packed_base_type_deinit
   end type packed_base_type

contains
   subroutine packed_base_type_init(this, nopen)
      class(packed_base_type) :: this
      integer, intent(in) :: nopen

      call packed_base_type_deinit(this)
      
      allocate(this%iorder(nopen))
      allocate(this%jpack(nopen))
      allocate(this%lpack(nopen))
      allocate(this%epack(nopen))
      allocate(this%inpack(nopen))
      allocate(this%j12pk(nopen))
      this%length = 0
   end subroutine 

   subroutine packed_base_type_deinit(this)
      class(packed_base_type) :: this
      if (allocated(this%iorder)) deallocate(this%iorder)
      if (allocated(this%jpack)) deallocate(this%jpack)
      if (allocated(this%lpack)) deallocate(this%lpack)
      if (allocated(this%epack)) deallocate(this%epack)
      if (allocated(this%inpack)) deallocate(this%inpack)
      if (allocated(this%j12pk)) deallocate(this%j12pk)
      this%length = 0
   end subroutine 

end module mod_hitypes