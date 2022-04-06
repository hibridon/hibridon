#include "assert.h"
!   cov2 module
!       nv2max: maximum core memory allocated for the v2 matrix
!       v2: lower triangle of nonzero angular coupling matrix elements stored in packed column form that is (1,1), (2,1), (3,1) ... (n,1), (2,2), (3,2) ... (n,2), etc.
use mod_cov2, only: nv2max, v2
!
!   coiv2 block
!       iv2: row+column index of v2 matrix for each non-zero element
use mod_coiv2, only: iv2
!   coiout block
!       niout: number of level indeces included in the output of hibridon
!       indout: level indeces included in the output of hibridon
use mod_coiout, only: niout, indout
!
!   cocent block
!       cent: array containing centrifugal barrier of each channel
use mod_cocent, only: cent

!   coeint block
!       eint: array containing channel energies (in hartree). The zero of energy is assumed to be the v2=0, j=0, k=0 level
use mod_coeint, only: eint
!
!   covvl block
!       vvl: r-dependence of each term in potential expansion
!   represent <v_2'|v_{\lambda\mu}(Q_2, R)|v_2> in the following order:
!       <0|v_00|0>, <0|v_20|0>, ..., <0|v_99|0>, <0|v_10|1>, ...
!   where v2 >= v2'
use mod_covvl, only: vvl
!     size of vvl : NVVL

use mod_conlam, only: nlam, nlammx, lamnum
use mod_cosys, only: scod
use mod_cosysr, only: isrcod, junkr, rspar
use mod_cosysi, only: nscode, isicod, ispar
implicit double precision (a-h, o-z)
!   Define the sizes of grids
!       V2MAX: maximum value of v2 included in the pot file
!       V2TMAX: number of (v2, v2') combination, C(V2MAX+1, 2)
!       NVLM: number of v_lm terms for each (v2, v2') combination
!       NVVL: total number of v_lm terms, for all (v2, v2') blocks
!       NTHETA, NPHI: number of theta/phi's in the ab initio calculation
!       NANGLE: number of (theta, phi) tuples
!       NDIST: number of distances included in the ab initio calculation
integer V2MAX, V2TMAX, NVLM, NTHETA, NPHI, NANGLE, NVVL, NDIST
parameter (V2MAX=3, V2TMAX=(V2MAX+1)*(V2MAX+2)/2)
parameter (NVLM=12, NVVL=NVLM*V2TMAX)
parameter (NTHETA=19, NPHI=7, NANGLE=NTHETA*NPHI)
parameter (NDIST=19)
!
!   Conversion factor
!       ECONV: hartree to wavenumber
!       XMCONV: amu to electrom mass
double precision ECONV, XMCONV
parameter (ECONV=219474.6d0, XMCONV=1822.88848477d0)
!
!   Max number of channels
!      integer KKMAX
!      parameter (KMAX=10000)
!
!   Lengths of cod array, 
!       ICOD, IRCOD: lenghts of cod array
integer ICOD, IRCOD
parameter (ICOD=5, IRCOD=4)
!
!   ch3he block: data used only by this pot/basis combination
!       brot, crot: rotational constants of CH3 for each vibrational level
!       evib: vibrational level energies
!       nlamsi: number of v_lm terms for each (v2, v2') combination
!       lamsym, musym: list of lambda/mu's used for the coupling potential symmetric to theta = 90 deg
!       lamasy, muasy: list of lambda/mu's used for the coupling potential anti-symmetric to theta = 90 deg
!
!     Source of rotational constants:
!     Yamada, C., et. al., JCP, 75, 5256
!     Amano, T., et. al., JCP, 77, 5284
!
common /ch3he/ brot, crot, evib, &
               lamsym, musym, lamasy, muasy
double precision brot(V2MAX+1), crot(V2MAX+1), evib(V2MAX+1)
integer lamsym(NVLM), musym(NVLM), lamasy(NVLM), muasy(NVLM)
data brot /9.57789d0, 9.25814d0, 8.93320d0, 8.60974d0/
data crot /4.74275d0, 4.811643d0, 4.871213d0, 4.92655d0/
data evib /0d0, 6.064531d2, 1.28809d3, 2.0191657d3/
!   lambda/mu's used in this pot routine
!
!   When <v2'|V|v2> is symmetric about theta=90 deg (even v2 + v2')
!       lambda =  0  2  4  6  8  3  5  7  9  6  8  9
!       mu     =  0  0  0  0  0  3  3  3  3  6  6  9
!   When <v2'|V|v2> is anti-symmetric about theta=90 deg (odd v2 + v2')
!       lambda =  1  3  5  7  9  4  6  8 10  7  9 10
!       mu     =  0  0  0  0  0  3  3  3  3  6  6  9
data lamsym /0, 2, 4, 6, 8, 3, 5, 7, 9, 6, 8, 9/
data musym /0, 0, 0, 0, 0, 3, 3, 3, 3, 6, 6, 9/
data lamasy /1, 3, 5, 7, 9, 4, 6, 8, 10, 7, 9, 10/
data muasy /0, 0, 0, 0, 0, 3, 3, 3, 3, 6, 6, 9/
!
!
!   cosysi block
!       nscod: total number of variable names which are passed to HINPUT, nscod must equal isrcod + isicod + 3
!       isicod: total number of integer system dependent variables
!       nterm: number of different associated legendre terms in expansion of potential
!       numpot: the number of the potential used, this variable is passed to the pot subroutine
!       ipotsy: cylindrical symmetry of potential. Should be set to 3 for CH3.
!       iop: ortho/para label for molecular states. Only para states are included if iop=1 and only ortho states if iop=-1.
!

!
!   cosysr block
!       isrcod: total number of real system dependent variables
!       junkr: junk variable (required by hibridon)
!       vmax: maximum value of v2 (starts from zero) included in the calculation
!       emax0, emax1, emax2, emax3: maximum total energy of a level to be included in the channel basis, for four vibrational levels

!
!   coered block
!       ered: collision energy in atomic units (hartrees)
!       rmu: collision reduced mass in atomic units
common /coered/ ered, rmu
double precision ered, rmu
!
!
!   coipar block
!       junkip: not refered to
!       iprint: level of printing in calculation
common /coipar/ junkip, iprint
integer junkip(9), iprint
real(8), pointer :: emax0, emax1, emax2, emax3
integer, pointer :: nterm, ipotsy, iop, jmax, vmax