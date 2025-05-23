!  system:  Ar-CH4 PES
!  written by P. Dagdigian, Jul 2015
!
!  PES calculated by Heijmen et al. (JCP 107, 902 (1997).  Fortran code
!  to generate the tetrahedrl expansion coefficients provided by
!  A. van der Avoird (July 2015)
!
!  expansion coefficients determined by quadrature, as described in
!  Heijmen et al. paper (NOTE:  these coefficients slightly different than
!  those plotted in Fig. 4.)
!
!  the PES is fitted with 5 angular terms
!
#include "unused.h"
#include "common/syusr.F90"
#include "common/ground.F90"
#include "common/bausr.F90"
! ------------------------------------------------------------------------
subroutine driver
use mod_covvl, only: vvl
use mod_parpot, only: potnam=>pot_name
use mod_hipot, only: pot
implicit double precision (a-h,o-z)
econv=219474.6d0
potnam='Ar-CH4 Nijmegen 1997'
print *, potnam
1 print *, 'R (bohr):'
read (5, *, end=99) r
call pot(vv0,r)
if (r.le.0.d0) goto 99
!     vlm coefficient is returned in atomic units (hartree)
!     convert from atomic units for printout
write (6, 100) vv0*econv, (econv*vvl(i), i=1,4)
100 format(6(1pe16.8))
goto 1
99 rr=3.5d0
dr=0.5d0
open (unit=12,file='arch4_vlms.dat')
write(12,109)
109 format(' %R/bohr V0   V3    V4     V6    V7')
do i=1,40
  call pot(vv0,rr)
  write (12,110) rr,vv0*econv, (econv*vvl(j),j=1,4)
110   format(f7.2,5(1pe16.8))
  rr = rr + dr
enddo
close(12)
end
! ------------------------------------------------------------------------
subroutine loapot(iunit,filnam)
use mod_conlam, only: nlam, nlammx, lamnum
use mod_cosysi, only: ispar
use mod_parbas, only: ntv, ivcol, ivrow, lammin, lammax, mproj
use mod_parpot, only: potnam=>pot_name
use mod_selb, only: ibasty
implicit double precision (a-h,o-z)
integer, intent(in) :: iunit  ! if a data file is used, this subroutine is expected to use this unit to open it in read mode (not used here)
character*(*), intent(in) :: filnam  ! if a data file is used, the file name of the data file (not used here)    

integer, pointer :: nterm
nterm=> ispar(1)
UNUSED_DUMMY(iunit)
UNUSED_DUMMY(filnam)
potnam='Ar-CH4 Nijmegen 1997'
ibasty = 24
!
nterm = 4
lammin(1) = 3
lammin(2) = 4
lammin(3) = 6
lammin(4) = 7
do i=1,nterm
  mproj(i) = 0
  lammax(i) = lammin(i)
  lamnum(i) = 1
enddo
!  total number of anisotropic terms
nlammx = nterm
nlam = nterm
!
ntv(1)=1
ivcol(1,1)=0
ivrow(1,1)=0
return
end
! ----------------------------------------------------------------------
subroutine pot (vv0, r)
!  subroutine to calculate the r-dependent coefficients
!  in atomic units (distance and energy)
! ----------------------------------------------------------------------
!  on entry:
!    r:          interparticle distance
!  on return:
!    vv0         contains isotropic term (Y00)
!  variable in module mod_covvl
!    vvl:        vector of length 4 to store r-dependence of each term
!                in potential expansion
!  uses linear least squares routines from lapack
!
! author:  paul dagdigian
! latest revision date:  27-jul-2015
! ----------------------------------------------------------------------
! this pes is a fit to 19 values of R and 190 orientations, namely
! R = [3:0.5:10 11 12 13 15 20]
!
use mod_covvl, only: vvl
use mod_hiblas, only: dscal, dcopy
use mod_hipotutil, only: dspline, dsplint
implicit double precision (a-h,o-z)
real(8), intent(out) :: vv0
real(8), intent(in) :: r  ! intermolecular distance

real(8) :: v(5)
real(8), save :: csplin(69,5)
real(8) :: rr(69), vl(345),vec(69)
integer, save :: ifirst=0
! here are the 69 values of R
v=0d0 ; csplin=0  ; vec=0d0
data rr/ &
  3.000,  3.250,  3.500,  3.750,  4.000, &
  4.250,  4.500,  4.750,  5.000,  5.250, &
  5.500,  5.750,  6.000,  6.250,  6.500, &
  6.750,  7.000,  7.250,  7.500,  7.750, &
  8.000,  8.250,  8.500,  8.750,  9.000, &
  9.250,  9.500,  9.750, 10.000, 10.250, &
 10.500, 10.750, 11.000, 11.250, 11.500, &
 11.750, 12.000, 12.250, 12.500, 12.750, &
 13.000, 13.250, 13.500, 13.750, 14.000, &
 14.250, 14.500, 14.750, 15.000, 15.250, &
 15.500, 15.750, 16.000, 16.250, 16.500, &
 16.750, 17.000, 17.250, 17.500, 17.750, &
 18.000, 18.250, 18.500, 18.750, 19.000, &
 19.250, 19.500, 19.750, 20.000 / 
! here are column ordered angular expansion coefficients (5 total) for the
! potential at each of 69 values of R (345 values)
data vl/ &
  1.6287369e+04,  4.4046281e+04,  4.6359940e+04,  3.9291904e+04, &
  3.0089847e+04,  2.1677293e+04,  1.4960028e+04,  9.9775120e+03, &
  6.4547340e+03,  4.0502760e+03,  2.4555580e+03,  1.4247960e+03, &
  7.7534900e+02,  3.7745800e+02,  1.4186100e+02,  8.6910000e+00, &
 -6.1351000e+01, -9.3562000e+01, -1.0391400e+02, -1.0228500e+02, &
 -9.4641000e+01, -8.4475000e+01, -7.3754000e+01, -6.3512000e+01, &
 -5.4234000e+01, -4.6094000e+01, -3.9091000e+01, -3.3143000e+01, &
 -2.8128000e+01, -2.3920000e+01, -2.0395000e+01, -1.7443000e+01, &
 -1.4969000e+01, -1.2890000e+01, -1.1140000e+01, -9.6620000e+00, &
 -8.4090000e+00, -7.3440000e+00, -6.4350000e+00, -5.6560000e+00, &
 -4.9860000e+00, -4.4080000e+00, -3.9080000e+00, -3.4740000e+00, &
 -3.0960000e+00, -2.7650000e+00, -2.4760000e+00, -2.2210000e+00, &
 -1.9970000e+00, -1.7990000e+00, -1.6240000e+00, -1.4680000e+00, &
 -1.3300000e+00, -1.2070000e+00, -1.0970000e+00, -9.9900000e-01, &
 -9.1000000e-01, -8.3100000e-01, -7.6000000e-01, -6.9600000e-01, &
 -6.3800000e-01, -5.8500000e-01, -5.3800000e-01, -4.9500000e-01, &
 -4.5600000e-01, -4.2100000e-01, -3.8900000e-01, -3.5900000e-01, &
 -3.3200000e-01, &
 -2.3045716e+04,  1.6147254e+04,  2.7963537e+04,  2.7878401e+04, &
  2.3437793e+04,  1.8095866e+04,  1.3263332e+04,  9.3791520e+03, &
  6.4562090e+03,  4.3485290e+03,  2.8744320e+03,  1.8674630e+03, &
  1.1927960e+03,  7.4835800e+02,  4.6013900e+02,  2.7609200e+02, &
  1.6044700e+02,  8.9077000e+01,  4.5963000e+01,  2.0618000e+01, &
  6.2640000e+00, -1.4190000e+00, -5.1490000e+00, -6.6090000e+00, &
 -6.8270000e+00, -6.4120000e+00, -5.7200000e+00, -4.9470000e+00, &
 -4.1960000e+00, -3.5180000e+00, -2.9280000e+00, -2.4280000e+00, &
 -2.0120000e+00, -1.6680000e+00, -1.3870000e+00, -1.1560000e+00, &
 -9.6800000e-01, -8.1400000e-01, -6.8700000e-01, -5.8300000e-01, &
 -4.9700000e-01, -4.2600000e-01, -3.6600000e-01, -3.1600000e-01, &
 -2.7400000e-01, -2.3900000e-01, -2.0900000e-01, -1.8300000e-01, &
 -1.6100000e-01, -1.4100000e-01, -1.2500000e-01, -1.1100000e-01, &
 -9.8000000e-02, -8.7000000e-02, -7.8000000e-02, -7.0000000e-02, &
 -6.2000000e-02, -5.6000000e-02, -5.0000000e-02, -4.5000000e-02, &
 -4.0000000e-02, -3.7000000e-02, -3.3000000e-02, -3.0000000e-02, &
 -2.7000000e-02, -2.4000000e-02, -2.2000000e-02, -2.0000000e-02, &
 -1.8000000e-02, &
  4.6837293e+04,  1.2599032e+04, -1.4628790e+03, -6.1199430e+03, &
 -6.7105150e+03, -5.7813780e+03, -4.4774900e+03, -3.2620860e+03, &
 -2.2824190e+03, -1.5505090e+03, -1.0290000e+03, -6.6947700e+02, &
 -4.2774100e+02, -2.6846700e+02, -1.6535500e+02, -9.9679000e+01, &
 -5.8517000e+01, -3.3155000e+01, -1.7830000e+01, -8.7850000e+00, &
 -3.6100000e+00, -7.7800000e-01,  6.6600000e-01,  1.3110000e+00, &
  1.5130000e+00,  1.4830000e+00,  1.3440000e+00,  1.1650000e+00, &
  9.8300000e-01,  8.1400000e-01,  6.6700000e-01,  5.4300000e-01, &
  4.4100000e-01,  3.5700000e-01,  2.9000000e-01,  2.3600000e-01, &
  1.9300000e-01,  1.5800000e-01,  1.3000000e-01,  1.0800000e-01, &
  9.0000000e-02,  7.5000000e-02,  6.3000000e-02,  5.4000000e-02, &
  4.6000000e-02,  3.9000000e-02,  3.4000000e-02,  2.9000000e-02, &
  2.5000000e-02,  2.2000000e-02,  1.9000000e-02,  1.7000000e-02, &
  1.5000000e-02,  1.3000000e-02,  1.1000000e-02,  1.0000000e-02, &
  9.0000000e-03,  8.0000000e-03,  7.0000000e-03,  6.0000000e-03, &
  6.0000000e-03,  5.0000000e-03,  4.0000000e-03,  4.0000000e-03, &
  4.0000000e-03,  3.0000000e-03,  3.0000000e-03,  3.0000000e-03, &
  2.0000000e-03, &
 -2.6455581e+04, -1.0512445e+04, -3.3217010e+03, -3.3239100e+02, &
  7.2266000e+02,  9.4458000e+02,  8.4921200e+02,  6.6519400e+02, &
  4.8459000e+02,  3.3748600e+02,  2.2792700e+02,  1.5054400e+02, &
  9.7788000e+01,  6.2727000e+01,  3.9873000e+01,  2.5203000e+01, &
  1.5902000e+01,  1.0063000e+01,  6.4250000e+00,  4.1680000e+00, &
  2.7690000e+00,  1.8980000e+00,  1.3500000e+00,  9.9900000e-01, &
  7.6600000e-01,  6.0600000e-01,  4.9100000e-01,  4.0500000e-01, &
  3.3700000e-01,  2.8200000e-01,  2.3600000e-01,  1.9700000e-01, &
  1.6500000e-01,  1.3700000e-01,  1.1400000e-01,  9.4000000e-02, &
  7.8000000e-02,  6.4000000e-02,  5.2000000e-02,  4.2000000e-02, &
  3.4000000e-02,  2.8000000e-02,  2.2000000e-02,  1.8000000e-02, &
  1.4000000e-02,  1.1000000e-02,  9.0000000e-03,  7.0000000e-03, &
  6.0000000e-03,  5.0000000e-03,  4.0000000e-03,  3.0000000e-03, &
  2.0000000e-03,  2.0000000e-03,  1.0000000e-03,  1.0000000e-03, &
  1.0000000e-03,  1.0000000e-03,  0.0000000e+00,  0.0000000e+00, &
  0.0000000e+00,  0.0000000e+00,  0.0000000e+00,  0.0000000e+00, &
  0.0000000e+00,  0.0000000e+00,  0.0000000e+00,  0.0000000e+00, &
  0.0000000e+00, &
  1.9051572e+04,  1.0157816e+04,  5.3310390e+03,  2.7402960e+03, &
  1.3690300e+03,  6.5626900e+02,  2.9470900e+02,  1.1749800e+02, &
  3.5041000e+01, -1.0500000e-01, -1.2605000e+01, -1.4974000e+01, &
 -1.3339000e+01, -1.0549000e+01, -7.8110000e+00, -5.5460000e+00, &
 -3.8300000e+00, -2.5950000e+00, -1.7390000e+00, -1.1610000e+00, &
 -7.8000000e-01, -5.3100000e-01, -3.7000000e-01, -2.6700000e-01, &
 -2.0000000e-01, -1.5600000e-01, -1.2600000e-01, -1.0500000e-01, &
 -8.9000000e-02, -7.6000000e-02, -6.6000000e-02, -5.7000000e-02, &
 -4.9000000e-02, -4.2000000e-02, -3.6000000e-02, -3.1000000e-02, &
 -2.6000000e-02, -2.2000000e-02, -1.9000000e-02, -1.6000000e-02, &
 -1.3000000e-02, -1.1000000e-02, -9.0000000e-03, -7.0000000e-03, &
 -6.0000000e-03, -5.0000000e-03, -4.0000000e-03, -3.0000000e-03, &
 -3.0000000e-03, -2.0000000e-03, -2.0000000e-03, -1.0000000e-03, &
 -1.0000000e-03, -1.0000000e-03, -1.0000000e-03, -1.0000000e-03, &
 -0.0000000e+00, -0.0000000e+00, -0.0000000e+00, -0.0000000e+00, &
 -0.0000000e+00, -0.0000000e+00, -0.0000000e+00, -0.0000000e+00, &
 -0.0000000e+00, -0.0000000e+00, -0.0000000e+00, -0.0000000e+00, &
 -0.0000000e+00 /
!
! spline fit
if (ifirst .eq. 0) then
! spline fit of the vl coefficients
   ind=1
   do ilam=1,5
     call dcopy(69,vl(ind),1,vec,1)
!    evaluate derivative at first point
     der1=(vec(2)-vec(1))/(rr(2)-rr(1))
     call dspline(rr,vec,69,der1,0d0,csplin(1,ilam))
     ind = ind + 69
   enddo
   ifirst = 1
 end if
! r^-6 fit to isotropic part of potential
 c6 = vl(63)*rr(63)**6
! switching function for long-range
 switch_lr=0.5*(tanh(0.5*(r - 18.d0)) + 1.d0)
! determine splined coefficients at R=r
 ind=1
 do ilam=1,5
   call dcopy(69,vl(ind),1,vec,1)
   call dsplint(rr,vec,csplin(1,ilam),69,r,vvx)
! kill anisotropic terms at large R
   vvx = (1.d0 - switch_lr)*vvx
   if (ilam.eq.1) then
! merge with asymptotic form
      vvx = vvx + switch_lr*c6/(r**6)
   endif
   v(ilam)=vvx
   call dcopy(69,vl(ind),1,vec,1)
   ind = ind + 69
 enddo
 call dcopy(4,v(2),1,vvl(1),1)
! convert to hartree
 econv=1.d0/219474.6d0
 call dscal(4,econv,vvl,1)
 vv0 = v(1) * econv
!
 return
 end
