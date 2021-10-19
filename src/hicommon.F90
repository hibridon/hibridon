!     This file contains various modules that replace common blocks in
!     hibridon.
!
!     Physical constants
module constants
implicit none
real(8), parameter :: pi=dacos(-1d0)
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
