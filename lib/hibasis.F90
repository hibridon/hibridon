!************************************************************************
!                                                                       *
!                            basis  library                             *
!                                                                       *
!************************************************************************
!                          routines included:                           *
!                                                                       *
!  1. basis      dispatcher routine to select basis                     *
!                see in header for basis subroutine for list of basis   *
!                routines currently available                           *
!  4. is_twomol check if a basis is for mol-mol collision (j=10j1+j2)   *
!  5. is_j12    check if j12 is used in a basis                         *
!                                                                       *
!     current revision:  24-jul-2019 (p.dagdigian)                       *
!                                                                       *
!************************************************************************
module mod_hibasis
#include "assert.h"
contains
! --------------------------------------------------------------------
subroutine basis (j, l, is, jhold, ehold, ishold, nlevel, nlevop, &
                  sc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu, &
                  csflag, clist, bastst, ihomo, nu, numin, jlpar, &
                  twomol, n, nmax, ntop, v2)
! --------------------------------------------------------------------
!
!  dispatcher to select the correct basis routine for given problem
!  the variable in common block /coselb/ ibasty determines the
!  basis used in the current calculation
!
!  basis routines currently available are:
!  ibasty         basis routine          kind of problem
!    1              ba1sg             singlet sigma scattering
!    2              ba2sg             doublet sigma scattering
!    3              ba2pi             doublet pi scattering
!    4              basgpi            sigma/pi scattering
!    5              bapi              general pi scattering
!    6              bastp             symmetric top scattering
!    7              ba13p             1/3 P atom scattering
!    8              ba2mol            1sigma + 1sigma
!    9              bastpln           symmetric top + linear molecule
!    10             ba22p             2/2 P atom scattering
!    11             ba1del            singlet delta scattering
!    12             bah2p             homonuclear + 2P atom
!    13             bah3p             homonuclear + 3P atom
!    14             ba2del            doublet delta scattering
!    15             badiat2p          heteronuclear + 2P atom  **to check
!    16             baastp            asymmetric top scattering
!    17             bach2x            CH2(X B1) (0,v2,0) bender levels
!    18             bastp1            symmetric top - no inversion doubling
!    19             basgpi1           2sigma | 2pi + atom (no pertubations)
!    20             ba2pi1sg          2pi molecule + 1sigma molecules
!    21             bastp1sg          symmetric top + 1sigma molecules
!    22             ba1d3p            1D/3P atom + closed-shell atom
!    23             ba3p2s            3P atom + 2S atom
!    24             basphtp           spherical top + atom scattering
!    25             ba1sg1sg           1sigma + 1sigma (different molecules)
!    26             ba2sg1sg          2sigma + 1sigma molecules
!    27             baastp1           C2v asymmetric top scattering, w body-frame
!                                     quant axis along C2 axis (compatible w MOLSCAT)
!    28             ba3sg1sg          3sigma + 1sigma molecules
!    29             baastp2           chiral asymmetric top + atom scattering
!    30             baastp3           C2v asymmetric top + linear molecule
!    99 or higher   bausr             user defined basis
!  author: b. follmeg
!  current revision of list:  20-jun-2019 (p.dagdigian)
! --------------------------------------------------------------------
!  variables in call list:
!    j:        on return contains rotational quantum numbers for each
!              channel
!    l:        on return contains orbital angular momentum for each
!              channel
!    is:       on return contains symmetry index for each channel
!    jhold:    on return contains rotational quantum numbers for each
!              rotational level
!    ehold:    on return contains energy in hartrees of each rotational
!              level
!    ishold:   on return contains symmetry index of each rotational level
!    nlevel:   on return contains number of energetically distinct
!              rotational levels used in channel basis
!    nlevop:   on return contains number of energetically distinct
!              rotational levels used in channel basis which are open
!              asymptotically
!    sc1,sc2:  scratch vectors of length at least nmax
!    sc3,sc4:  scratch vectors of length at least nmax
!              these scratch vectors are not used here
!    rcut:     cut-off point for keeping higher energy channels
!              if any open channel is still closed at r=rcut, then
!              all closed channels as well any open channels which are
!              still closed at r=rcut are dropped from basis
!              note that this is ignored in molecule-surface collisions!!!
!    jtot:     total angular momentum
!    flaghf:   if .true., then system with half-integer spin
!              if .false., then system with integer spin (this is the case
!              here)
!    flagsu:   if .true., then molecule-surface collisons
!    csflag:   if .true., then coupled-states calculation
!              if .false., the close-coupled calculation
!    clist:    if .true., then quantum numbers and energies listed for
!              each channel
!    bastst:   if .true., then execution terminates after the first call
!              to basis
!    ihomo:    if .true. , then homonuclear molecule
!              if .false., then heteronuclear molecule
!              if the molecule is homonuclear (ihomo = .true.), the
!              rotational levels included go from jmin to jmax in steps
!              of 2 and only even lambda terms in the anisotropy are
!              included
!    nu:       coupled-states projection index
!    numin:    minimum coupled states projection index
!              for cc calculations nu and numin are both set = 0 by calling
!              program
!    jlpar:    total parity of included channels in cc calculation
!              parity=(-1)**jlpar (by definition parity=(-1)**(j+l+jtot) )
!    n:        on return equals number of channels
!    nmax:     maximum dimension of arrays
!    ntop:     maximum row dimension of all matrices passed to subroutines
!              propag and soutpt.  ntop is set in basis only if nu = numin
!              otherwise it is unchanged from the value supplied by the
!              calling program
!    note!!!   if flaghf = .true., then the true values of the rotational
!    quantum numbers, the total angular momentum, and the coupled-states
!    projection index are equal to the values stored in j, jtot, and nu
!    plus 1/2
use mod_ancou, only: ancou_type
use mod_hiba01_1sg, only: ba1sg
use mod_hiba02_2sg, only: ba2sg
use mod_hiba03_2pi, only: ba2pi
use mod_hiba04_sgpi, only: basgpi
use mod_hiba05_pi, only: bapi
use mod_hiba06_stp, only: bastp
use mod_hiba07_13p, only: ba13p
use mod_hiba08_2mol, only: ba2mol
use mod_hiba09_stpln, only: bastpln
use mod_hiba10_22p, only: ba22p
use mod_hiba11_1del, only: ba1del
use mod_hiba12_h2p, only: bah2p
use mod_hiba13_h3p, only: bah3p
use mod_hiba14_2del, only: ba2del
use mod_hiba15_diat2p, only: badiat2p
use mod_hiba16_astp, only: baastp
use mod_hiba17_ch2x, only: bach2x
use mod_hiba18_stp1, only: bastp1
use mod_hiba19_sgpi1, only: basgpi1
use mod_hiba20_2pi1sg, only: ba2pi1sg
use mod_hiba21_stp1sg, only: bastp1sg
use mod_hiba22_1d3p, only: ba1d3p
use mod_hiba23_3p2s, only: ba3p2s
use mod_hiba24_sphtp, only: basphtp
use mod_hiba25_1sg1sg, only: ba1sg1sg
use mod_hiba26_2sg1sg, only: ba2sg1sg
use mod_hiba27_astp1, only: baastp1
use mod_hiba28_3sg1sg, only: ba3sg1sg
use mod_hiba29_astp2, only: baastp2
use mod_hiba30_astp3, only: baastp3
!use mod_bausr, only: bausr
!use mod_hibuser, only: bausr
use, intrinsic :: ISO_C_BINDING   ! for C_LOC and C_F_POINTER

implicit double precision (a-h,o-z)
type(ancou_type), intent(out), allocatable :: v2
integer :: j(:)
integer :: l(:)
integer :: is(:)
integer :: jhold(:)
real(8) :: ehold(:)
integer :: ishold(:)
real(8), target :: sc1(:), sc2(:), sc3(:), sc4(:)
integer nlevel, nlevop, jtot, nu, &
        jlpar, n, nmax
!      real ehold, sc1, sc2, sc3, sc4, rcut
logical flaghf, flagsu, csflag, clist, bastst, ihomo, twomol
integer ibasty
integer, pointer :: sc1_as_int(:), sc2_as_int(:), sc3_as_int(:), sc4_as_int(:)
#include "hibuser.inc.F90"
common /coselb/ ibasty
!
ASSERT(.not. allocated(v2))

call C_F_POINTER (C_LOC(sc1), sc1_as_int, [nmax])
call C_F_POINTER (C_LOC(sc2), sc2_as_int, [nmax])
call C_F_POINTER (C_LOC(sc3), sc3_as_int, [nmax])
call C_F_POINTER (C_LOC(sc4), sc4_as_int, [nmax])

!  select basis routine according to value of ibasty
if (ibasty .ge. 99) then
!  user supplied routine
  call bausr(j, l, is, jhold, ehold, ishold, nlevel, nlevop, &
                  sc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu, &
                  csflag, clist, bastst, ihomo, nu, numin, jlpar, &
                  n, nmax, ntop, v2)
  return
endif
goto (100,200,300,400,500,600,700,800,900,1000,1100,1200,1300, &
      1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400, &
      2500,2600,2700,2800,2900,3000) &
      ibasty
!  singlet sigma basis
100   call ba1sg (j, l, is, jhold, ehold, ishold, nlevel, nlevop, &
                  sc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu, &
                  csflag, clist, bastst, ihomo, nu, numin, jlpar, &
                  n, nmax, ntop, v2)
return
!  doublet sigma basis
200   call ba2sg (j, l, is, jhold, ehold, ishold, nlevel, &
                  nlevop, sc1_as_int, sc2, sc3, sc4, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, &
                  ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
return
!  doublet pi basis
300   call ba2pi (j, l, is, jhold, ehold, ishold, nlevel, &
                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, &
                  ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
return
!  sigma/pi basis
400   call basgpi (j, l, is, jhold, ehold, ishold, nlevel, &
                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, &
                  ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
return
!   general pi basis
500    call bapi (j, l, is, jhold, ehold, ishold, nlevel, &
                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, &
                  ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
return
!  symmetric top basis, with inversion doubling
600   call bastp (j, l, is, jhold, ehold, ishold, nlevel, &
                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, &
                  ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
return
!  1/3 P atom basis
700   call ba13p (j, l, is, jhold, ehold, ishold, nlevel, &
                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, &
                  ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
return
! 1sigma + 1sigma basis
800   call ba2mol (j, l, is, jhold, ehold, ishold, nlevel, &
                  nlevop, sc1_as_int, sc2_as_int, sc3_as_int, sc4, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, &
                  ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
return
! symmetric top + 1 sigma basis
900   call bastpln (j, l, is, jhold, ehold, ishold, nlevel, &
                  nlevop, sc1_as_int, sc2_as_int, sc3_as_int, sc4_as_int, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, &
                  ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
return
! 2/2 P atom basis
1000   call ba22p (j, l, is, jhold, ehold, ishold, nlevel, &
                  nlevop, sc1_as_int, sc2, sc3, sc4, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, &
                  ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
return
!  singlet delta basis
1100  call ba1del (j, l, is, jhold, ehold, ishold, nlevel, &
                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, &
                  ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
return
!  homonuclear + 2P atom basis
1200  call bah2p (j, l, is, jhold, ehold, ishold, nlevel, &
                  nlevop, sc2, sc3, sc4, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, &
                  ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
return
!  homonuclear + 3P atom basis
1300  call bah3p (j, l, is, jhold, ehold, ishold, nlevel, &
                  nlevop, sc1_as_int, sc2, sc3, sc4, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, &
                  ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
return
!  doublet delta basis
1400  call ba2del (j, l, is, jhold, ehold, ishold, nlevel, &
                  nlevop, sc1_as_int, sc2, sc3, sc4, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, &
                  ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
return
!  heteronuclear + 2P atom basis
!1500  call bah2p (j, l, is, jhold, ehold, ishold, nlevel, &
!    :                  nlevop, sc2, sc3, sc4, rcut, jtot, &
!    :                  flaghf, flagsu, csflag, clist, bastst, &
!    :                  ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
1500  call badiat2p (j, l, is, jhold, ehold, ishold, nlevel, &
                  nlevop, sc1_as_int, sc2, sc3, sc4, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, &
                  ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
!  asymmetric top basis
1600   call baastp (j, l, is, jhold, ehold, ishold, nlevel, &
                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, &
                  ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
return
!  CH2(X 3B1) (0,v2,0) bender levels
1700   call bach2x (j, l, is, jhold, ehold, ishold, nlevel, &
                  nlevop, sc1, sc2, sc3, sc4, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, &
                  ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
return
!  symmetric top basis, with no inversion doubling
1800   call bastp1 (j, l, is, jhold, ehold, ishold, nlevel, &
                  nlevop, sc1_as_int, sc2_as_int, sc3_as_int, sc4_as_int, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, &
                  ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
return
!  2sig-2pi + atom scattering (one 2sigma state and one or more 2pi
!   vibrational levels, no sigma-pi spectroscopic perturbations)
1900   call basgpi1 (j, l, is, jhold, ehold, ishold, nlevel, &
                  nlevop, sc1_as_int, sc2, sc3, sc4, rcut, jtot, &
                  flaghf, flagsu, csflag, clist, bastst, &
                  ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
return
!  2pi + 1sigma molecules
 2000 call ba2pi1sg(j, l, is, jhold, ehold, ishold, nlevel, &
     nlevop, rcut, jtot, flaghf, flagsu, csflag, clist, &
     bastst, ihomo, nu, numin, jlpar, twomol, n, nmax, ntop, v2)
return
!  symmetric top + 1sigma molecules
 2100 call bastp1sg(j, l, is, jhold, ehold, ishold, nlevel, &
     nlevop, rcut, jtot, flaghf, flagsu, csflag, clist, &
     bastst, ihomo, nu, numin, jlpar, twomol, n, nmax, ntop, v2)
return
!  1D/3P atom + closed-shell atom
 2200 call ba1d3p(j, l, is, jhold, ehold, ishold, nlevel, &
     nlevop, rcut, jtot, flaghf, flagsu, csflag, clist, &
     bastst, ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
return
!   3P atom + 2S atom
 2300 call ba3p2s (j, l, is, jhold, ehold, ishold, nlevel, &
     nlevop, rcut, jtot, flaghf, flagsu, csflag, clist, &
     bastst, ihomo, nu, numin, jlpar, n, nmax, ntop, v2)
return
!   spherical top + atom
 2400 call basphtp (j, l, is, jhold, ehold, ishold, nlevel, &
     nlevop, sc1, sc2, sc3, sc4, rcut, jtot, &
     flaghf, flagsu, csflag, clist, bastst, ihomo, &
     nu, numin, jlpar, n, nmax, ntop, v2)
return
!   1sigma + 1sigma basis (different molecules)
 2500 call ba1sg1sg (j, l, is, jhold, ehold, ishold, nlevel, &
     nlevop, sc1, sc2, sc3, sc4, rcut, jtot, &
     flaghf, flagsu, csflag, clist, bastst, ihomo, &
     nu, numin, jlpar, n, nmax, ntop, v2)
return
!   2sigma + 1sigma basis
 2600 call ba2sg1sg (j, l, is, jhold, ehold, ishold, nlevel, &
     nlevop, sc1, sc2, sc3, sc4, rcut, jtot, &
     flaghf, flagsu, csflag, clist, bastst, ihomo, &
     nu, numin, jlpar, n, nmax, ntop, v2)
return
!   C2v asymmetric top + atom scattering
 2700 call baastp1 (j, l, is, jhold, ehold, ishold, nlevel, &
     nlevop, sc1, sc2, sc3, sc4, rcut, jtot, &
     flaghf, flagsu, csflag, clist, bastst, ihomo, &
     nu, numin, jlpar, n, nmax, ntop, v2)
return
!   3sigma + 1sigma basis
 2800 call ba3sg1sg (j, l, is, jhold, ehold, ishold, nlevel, &
     nlevop, sc1, sc2, sc3, sc4, rcut, jtot, &
     flaghf, flagsu, csflag, clist, bastst, ihomo, &
     nu, numin, jlpar, n, nmax, ntop, v2)
return
!   chiral asymmetric top + atom scattering
 2900 call baastp2 (j, l, is, jhold, ehold, ishold, nlevel, &
     nlevop, sc1, sc2, sc3, sc4, rcut, jtot, &
     flaghf, flagsu, csflag, clist, bastst, ihomo, &
     nu, numin, jlpar, n, nmax, ntop, v2)
return
!   C2v asymmetric top + linear molecule scattering
 3000 call baastp3 (j, l, is, jhold, ehold, ishold, nlevel, &
     nlevop, sc1, sc2, sc3, sc4, rcut, jtot, &
     flaghf, flagsu, csflag, clist, bastst, ihomo, &
     nu, numin, jlpar, n, nmax, ntop, v2)
return
end
!     ------------------------------------------------------------------
logical function is_twomol(ibasty)
!     ------------------------------------------------------------------
!
!     checks if a basis is for molecule-molecule collision (j=10j1+j2)
!
!     written by q. ma
!     current revision:  24-jul-2019 (p.dagdigian)
!     ------------------------------------------------------------------
implicit none
integer, intent(in) :: ibasty
if ((ibasty .eq. 9) .or. (ibasty .eq. 20) .or. (ibasty .eq. 21) &
      .or. (ibasty .eq. 25) .or. (ibasty .eq. 26) &
      .or. (ibasty .eq. 28) .or. (ibasty .eq. 30) &
      .or. (ibasty .eq. 100)) &
     then
   is_twomol = .true.
else
   is_twomol = .false.
end if
return
end function is_twomol
!     ------------------------------------------------------------------
logical function is_j12(ibasty)
!     ------------------------------------------------------------------
!
!     checks if j12 is used in a basis
!
!     written by q. ma
!     current revision:  17-oct-2018 (p.dagdigian)
!     ------------------------------------------------------------------
implicit none
integer, intent(in) :: ibasty
if (is_twomol(ibasty) .or. (ibasty .eq. 12) &
     .or. (ibasty .eq. 13) .or. (ibasty .eq. 15) &
     .or. (ibasty .eq. 23)) then
   is_j12 = .true.
else
   is_j12 = .false.
end if
return
end function is_j12
end module mod_hibasis
