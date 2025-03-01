#include "unused.h"
! System:  CN(X 2Sigma)+He, ab initio RCCSD(T) PES's
! BASIS AVQZ Lique, Spielfiedel
! 
! Reference: François Lique, Annie Spielfiedel,
! Nicole Feautrier,Ioan F. Schneider, Jacek Klos,
! and Millard H. Alexander,JCP 132, 024303, 2010
!
!  Note:  this pot routine requires a data file to be in hibxx/bin/progs/potdata:
!         hecn_fitmlv.dat
!
subroutine driver
use mod_covvl, only: vvl
use mod_cosysr, only: rspar
use mod_parpot, only: potnam=>pot_name, label=>pot_label
implicit double precision (a-h,o-z)
real(8), pointer :: rshift, xfact
rshift=>rspar(1); xfact=>rspar(2)
econv=219474.6d0
potnam='He-CN(2Sigma) CCSDT PES'
print *, potnam
1  print *, ' r (bohr)'
rshift=0.5
xfact=0.8
read (5, *, end=99) r
!
!  modification to print out table of vlm's (12-aug-2014, p.dagdigian)
if (r.le.0.d0) goto 50
call pot(vv0,r)
write (6, 100) vv0,vvl
100 format(' vsum',/,7(1pe16.8))
goto 1
!
!  save table of vlm's (20-may-2013, p.dagdigian)
50 write(6,53)
53 format('Enter Rmin,Rmax, dR:')
read (5,*) rmin, rmax, dr
npts = nint((rmax - rmin)/dr) + 1
open (unit=22,file='hecnx_vlms.txt')
write(22,52)
52 format(' %R/bohr  V00  V10  V20  V30  V40  V50  V60  V70', &
  '  V80  V90  V10,0  V11,0  V12,0')
do 55 ir=1,npts
  rr = rmin + (ir - 1)*dr
  call pot(vv0,rr)
  write(22,110) rr,vv0*econv,(econv*vvl(j),j=1,12)
110   format(f7.2,12(1pe16.8))
55 continue
close(22)
!
99 end
#include "common/syusr.F90"
#include "common/bausr.F90"
#include "common/ground.F90"
! ----------------------------------------------------------------
subroutine loapot(iunit,filnam)
use mod_parbas, only: ntv, ivcol, ivrow, lammin, lammax, mproj
use mod_parpot, only: potnam=>pot_name, label=>pot_label
character*(*) filnam

potnam='He-CN(2Sigma) CCSDT PES'
lammin(1)=1
lammax(1)=12
mproj(1)=0
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
!    r:      interparticle distance
!  on return:
!  vv0        contains isotropic term (d00 term in vsum)
!  variable in module mod_covvl
!    vvl:     vector of length 6 to store r-dependence of each term
!             in potential expansion
!    vvl(1-12) expansion coefficients in dl0 (l=1:12) of vsum
! uses linear least squares routines from cmlib
! 0 degree for He-NP and 180 for He-PN
! ----------------------------------------------------------------------
use mod_covvl, only: vvl
implicit double precision (a-h,o-z)
dimension &
          vsum(13), &
          d0(169),thta(13), &
          dd0(13,13)
dimension work(1053)
!     dimension of work is determined by setting LWORK to -1 
!     and looking at the value of work(1)
data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
data lwork /1053/
!     lwork is size of work array for dgels
! coefficicients for d0 rotation matrices
! stored (by column) for each of 13 angles and for l=0:12
! angles are 0:15:180 
data d0/ &
     1.d0, &
     1.d0, &
     1.d0, &
     1.d0, &
     1.d0, &
     1.d0, &
     1.d0, &
     1.d0, &
     1.d0, &
     1.d0, &
     1.d0, &
     1.d0, &
     1.d0, & !l=0
     1.000000000000000d0, &
     0.965925826289069d0, &
     0.866025403784439d0, &
     0.707106781186547d0, &
                   0.5d0, &
     0.258819045102521d0, &
                     0d0, &
    -0.258819045102521d0, &
                  -0.5d0, &
    -0.707106781186548d0, &
    -0.866025403784439d0, &
    -0.965925826289069d0, &
    -1.000000000000000d0, & !l=1
                     1d0, &
     0.899519052838329d0, &
                 0.625d0, &
                  0.25d0, &
                -0.125d0, &
    -0.399519052838329d0, &
                  -0.5d0, &
    -0.399519052838329d0, &
                -0.125d0, &
                  0.25d0, &
                 0.625d0, &
     0.899519052838329d0, &
                     1d0, & !l=2
                     1d0, &
     0.804163923099993d0, &
     0.324759526419165d0, &
    -0.176776695296637d0, &
               -0.4375d0, &
    -0.344884596328147d0, &
                     0d0, &
     0.344884596328147d0, &
                0.4375d0, &
     0.176776695296637d0, &
    -0.324759526419165d0, &
    -0.804163923099993d0, &
                    -1d0, & !l=3
                     1d0, &
     0.684695438682637d0, &
    0.0234375000000004d0, &
              -0.40625d0, &
            -0.2890625d0, &
     0.143429561317363d0, &
                 0.375d0, &
     0.143429561317363d0, &
            -0.2890625d0, &
              -0.40625d0, &
    0.0234375000000004d0, &
     0.684695438682637d0, &
                     1d0, & !l=4
                     1d0, &
     0.547125874778594d0, &
    -0.223272174413175d0, &
    -0.375650477505353d0, &
            0.08984375d0, &
     0.342727820841858d0, &
                     0d0, &
    -0.342727820841858d0, &
           -0.08984375d0, &
     0.375650477505353d0, &
     0.223272174413175d0, &
    -0.547125874778594d0, &
                    -1d0, & !l=5
                     1d0, &
     0.398305991010481d0, &
         -0.3740234375d0, &
            -0.1484375d0, &
          0.3232421875d0, &
    0.0431002589895194d0, &
               -0.3125d0, &
    0.0431002589895194d0, &
          0.3232421875d0, &
            -0.1484375d0, &
         -0.3740234375d0, &
     0.398305991010481d0, &
                     1d0, & !l=6
                     1d0, &
     0.245541045229049d0, &
    -0.410178047690872d0, &
     0.127058249744458d0, &
         0.22314453125d0, &
     -0.27304996323882d0, &
                     0d0, &
      0.27304996323882d0, &
        -0.22314453125d0, &
    -0.127058249744457d0, &
     0.410178047690872d0, &
    -0.245541045229049d0, &
                    -1d0, & !l=7
                     1d0, &
    0.0961843272422365d0, &
    -0.338775634765625d0, &
         0.29833984375d0, &
    -0.073638916015625d0, &
    -0.170219971773485d0, &
             0.2734375d0, &
    -0.170219971773485d0, &
    -0.073638916015625d0, &
         0.29833984375d0, &
    -0.338775634765625d0, &
    0.0961843272422365d0, &
                     1d0, & !l=8
                     1d0, &
   -0.0427678470871814d0, &
    -0.189575202067438d0, &
     0.285535794942029d0, &
    -0.267898559570312d0, &
     0.159493867392234d0, &
                     0d0, &
    -0.159493867392234d0, &
     0.267898559570312d0, &
    -0.285535794942029d0, &
     0.189575202067438d0, &
    0.0427678470871814d0, &
                    -1d0, & !l=9
                     1d0, &
    -0.165055973786964d0, &
  -0.00703811645507864d0, &
       0.1151123046875d0, &
    -0.188228607177734d0, &
     0.231630070466652d0, &
           -0.24609375d0, &
     0.231630070466652d0, &
    -0.188228607177734d0, &
       0.1151123046875d0, &
  -0.00703811645507864d0, &
    -0.165055973786964d0, &
                     1d0, & !l=10
                     1d0, &
    -0.265489992206791d0, &
     0.160704825466514d0, &
    -0.104184312120626d0, &
    0.0638713836669922d0, &
     -0.03054390246936d0, &
                     0d0, &
      0.03054390246936d0, &
   -0.0638713836669922d0, &
     0.104184312120625d0, &
    -0.160704825466514d0, &
     0.265489992206791d0, &
                    -1d0, & !l=11
                     1d0, &
    -0.340215667541778d0, &
     0.273202657699585d0, &
    -0.246719360351563d0, &
     0.233752965927124d0, &
    -0.227479473296845d0, &
          0.2255859375d0, &
    -0.227479473296845d0, &
     0.233752965927124d0, &
    -0.246719360351563d0, &
     0.273202657699585d0, &
    -0.340215667541778d0, &
                     1d0/ !l=12
data ifirst /0/
if (ifirst .eq. 0) then
   call init_potential()
   ifirst=1
endif
icount=0
do i=1,13
do j=1,13
 icount=icount+1
 dd0(i,j)=d0(icount)
enddo
enddo
do i=1,13
thta(i)=0.D0+(i-1)*15.D0
enddo
! detemine potential at angles
do 100 i=1,13
 call potential(r,thta(i),V)
  vsum(i)=V
100 continue
! solve simultaneous equations for solutions
! first for vsigma
tol=1.e-10
!      call dcopy(169,d0,1,aa,1)
!      call dqrank(aa,13,13,13,tol,kr,kpvt,qraux,work)
!      call dqrlss(aa,13,13,13,kr,vsum,xsum,rsd,kpvt,qraux)
call dgels('T',13,13,1,dd0,13,vsum,13,work,lwork,info)
! convert to hartree
conv=1.d0/1.d6
call dscal(13,conv,vsum,1)
vv0=vsum(1)
call dcopy(12,vsum(2),1,vvl,1)
end

module pot_hech_dgels
implicit none
save
! max_b: maximum number of blocks
! max_p
! max_r
! max_rr
! max_par : maximum number of parameters which can be fitted
!max_theta
integer*4, parameter :: max_b = 90
integer*4, parameter :: max_p = 10
integer*4, parameter :: max_r = 37
integer*4, parameter :: max_rr = 70
integer*4, parameter :: max_par = 50
integer*4, parameter :: max_theta = 13


type, public ::  fitmlv
character*80 :: label
integer*4 :: nterms
integer*4, dimension(max_p) :: npa
real*8 a(max_p,max_par,max_b)
integer*4, dimension(max_p) :: ntotal
integer*4, dimension(max_p) :: nblk
integer*4, dimension(max_p) :: nangle
integer*4, dimension(max_p) :: nr
integer*4, dimension(max_p) :: isym
integer*4  npblk(max_p,max_b)
real*8, dimension(max_p,max_b) :: angle
real*8 r(max_p,max_b)
real*8 rr(max_p,max_rr,max_b)
real*8 v(max_p,max_rr,max_b)
real*8 pinv(max_p,max_theta,max_theta) 
integer*4 maxpws(max_p)
integer*4 minmps(max_p)
integer*4 maxmps(max_p)
integer*4 mpsstp(max_p)
real*8 updmax(max_p)
real*8 wfac(max_p)
real*8 ang(max_p,max_theta)
real*8 rmol(max_p,max_r)
integer*4 mld(max_p)
real*8 re(max_p)
logical fitflg(max_p,max_b)
real*8 rinv(max_p,max_r,max_r)
integer*4 maxit(max_p)
real*8 eps(max_p)
real*8 xr(max_p,max_b)
integer*4 :: ntblk
integer*4 :: ntot
integer*4 :: ntpa

contains

procedure :: read_from_file

end type fitmlv

type(fitmlv) :: mlv

contains
subroutine read_from_file(this, fitmlv_filepath)
implicit none
class(fitmlv) :: this
character(len=255), intent(in) :: fitmlv_filepath


real*8 bidon1(max_par,max_b),bidon2(max_b)
real*8 bidon3(max_rr,max_b),bidon4(max_theta,max_theta)
real*8 bidon5(max_r)
logical bidon6(max_b)
real*8 bidon7(max_r,max_r),bidon8(max_theta)
integer*4 bidon9(max_b)
integer*4 input

integer*4 :: ierr
real*8 theta,det,rd
real*8 dlmm0
integer :: i, j, l, n
integer :: s1, s2
integer :: lmax
integer*4 :: ipot

this%ntblk = 0
this%ntpa = 0
this%ntot = 0
ierr=0
input = 6

open(unit=1,file= &
    fitmlv_filepath, &
     status='OLD',form='FORMATTED')
read(1,"(a80)") this%label
write(input,*) 'label : ', this%label
read(1,"(i10)") this%nterms
write(input,*) 'nterms : ',this%nterms
do ipot=1,this%nterms
  read(1,"(i10)") this%npa(ipot)
  read(1,"(e30.20)") bidon1
  do s1=1,max_par
    do s2=1,max_b
      this%a(ipot,s1,s2)=bidon1(s1,s2)
    enddo
  enddo
  read(1,"(i10)") this%ntotal(ipot)
  read(1,"(i10)") this%nblk(ipot)
  read(1,"(i10)") this%nangle(ipot)
  read(1,"(i10)") this%nr(ipot)
  read(1,"(i10)") this%isym(ipot)
  read(1,"(i10)") bidon9
  do s1=1,max_b
    this%npblk(ipot,s1)=bidon9(s1)
  enddo
  read(1,"(e30.20)") bidon2
  do s1=1,max_b
    this%angle(ipot,s1)=bidon2(s1)
  enddo
  read(1,"(e30.20)") bidon2
  do s1=1,max_b
    this%r(ipot,s1)=bidon2(s1)
  enddo
  read(1,"(e30.20)") bidon3
  do s1=1,max_rr
    do s2=1,max_b
      this%rr(ipot,s1,s2)=bidon3(s1,s2)
    enddo
  enddo
  read(1,"(e30.20)") bidon3
  do s1=1,max_rr
    do s2=1,max_b
      this%v(ipot,s1,s2)=bidon3(s1,s2)
    enddo
  enddo
  read(1,"(i10)") this%maxpws(ipot)
  read(1,"(i10)") this%minmps(ipot)
  read(1,"(i10)") this%maxmps(ipot)
  read(1,"(i10)") this%mpsstp(ipot)
  read(1,"(e30.20)") this%updmax(ipot)
  read(1,"(e30.20)") this%wfac(ipot)
  read(1,"(e30.20)") bidon8
  do s1=1,max_theta
    this%ang(ipot,s1)=bidon8(s1)
  enddo
  read(1,"(e30.20)") bidon5
  do s1=1,max_r
    this%rmol(ipot,s1)=bidon5(s1)
  enddo
  read(1,"(e30.20)") this%re(ipot)
  read(1,"(i10)") this%mld(ipot)
  read(1,"(l1)") bidon6
  do s1=1,max_b
    this%fitflg(ipot,s1)=bidon6(s1)
  enddo
  read(1,"(e30.20)") bidon4
  do s1=1,max_theta
    do s2=1,max_theta
      this%pinv(ipot,s1,s2)=bidon4(s1,s2)
    enddo
  enddo
  read(1,"(e30.20)") bidon2
  do s1=1,max_b
    this%xr(ipot,s1)=bidon2(s1)
  enddo
enddo
close(1)

  do ipot=1,this%nterms
    this%ntblk = this%ntblk + this%nblk(ipot)
    write(input,*) 'ntblk', this%ntblk
    this%ntot = this%ntot + this%ntotal(ipot)
    this%ntpa = this%ntpa + this%npa(ipot)
  enddo
!       initialize matrices
  do ipot=1,this%nterms
! DEBUT INIMAT
!...subroutine to initialize transformation matrices
!                   -1
!...job 1: compute p   maxtrix
    lmax=(this%nangle(ipot)-1)* this%isym(ipot)+ this%mld(ipot)
!
!..built up p matrix for calculated angles
    do i=1,this%nangle(ipot)
      j=0
      theta=this%ang(ipot,i)
      do l=this%mld(ipot),lmax,this%isym(ipot)
        j=j+1
        this%pinv(ipot,i,j)=dlmm0(l, this%mld(ipot), theta)
      enddo
    enddo
!...now invert matrix
    do s1=1,max_theta
      do s2=1,max_theta
        bidon4(s1,s2)=this%pinv(ipot, s1, s2)
      enddo
    enddo
    call matinv(bidon4, det, max_theta, this%nangle(ipot))
    do s1=1,max_theta
      do s2=1,max_theta
        this%pinv(ipot,s1,s2) = bidon4(s1,s2)
      enddo
    enddo
!                   -1
!...job 2: compute r
!..built up r matrix
    do i=1,this%nr(ipot)
      j = 1
      rd = this%rmol(ipot,i) - this%re(ipot)
      this%rinv(ipot,1,i) = 1d0
      do n=1,this%nr(ipot)
        j = j + 1
        this%rinv(ipot,j,i) = rd
        rd = rd * (this%rmol(ipot,i) - this%re(ipot))
      enddo
    enddo
!...matrix inversion
    do s1=1,max_r
      do s2=1,max_r
        bidon7(s1,s2) = this%rinv(ipot,s1,s2)
      enddo
    enddo
    call matinv(bidon7, det, max_r, this%nr(ipot))
    do s1=1,max_r
      do s2=1,max_r
        this%rinv(ipot,s1,s2) = bidon7(s1,s2)
      enddo
    enddo
! FIN INIMAT
  enddo
  ipot=1
  write(input,'(a,a,/)') ' SURFACE LOADED FROM FILE ',trim(fitmlv_filepath)
  write(input,*) &
' LABEL:',this%label(:72), &
'  NUMBER OF POTENTIALS (NTERMS): ',this%nterms, &
'  TOTAL NUMBER OF BLOCKS:        ',this%ntblk, &
'  TOTAL NUMBER OF POINTS:       ',this%ntot, &
'  TOTAL NUMBER OF COEFFICIENTS:  ',this%ntpa
  write(input,*) 'nangle',this%nangle(ipot)

end subroutine
end module 



subroutine init_potential()
!...subroutine to load fitted surface on unformatted files (filnam.srf).
!
!   on entry:   filnam  filename
! NOM DU FICHIER ICI ET NUMERO DU CANAL D'AFFICHAGE DES RESULTATS

use pot_hech_dgels, only: mlv
implicit none
character(len=255), parameter :: file_path='potdata/hecn_fitmlv.dat'
call mlv%read_from_file(file_path)
end subroutine 

subroutine potential(R_,angl,V_)
use pot_hech_dgels, only: mlv, max_b, max_p, max_r, max_rr, max_par, max_theta
implicit none
real*8 R_,angl,V_

! DEBUT DECLARATION FONCTION dlmm0
! FIN DECLARATION FONCTION dlmm0

integer :: i, j, l, n, ipot
integer :: lmax
! DEBUT TYPE POT
! FIN TYPE POT
! DEBUT DECLARATION DES VARIABLES
integer*4 iblk
logical exstfl,openfl
real*8 vcalc,vvlsum,vvl(max_theta),vvlp(max_theta)
real*8 pvec(max_theta)
real*8 rvecp(max_theta),rvecpp(max_theta),b(max_theta,max_r),dr
real*8 fct,tanhyb,fex,pfac,rpowi,pexfac,fmfac,jst,powmj
real*8 dlmm0
real(8) :: rd

! FIN DECLARATION DES VARIABLES
! DEBUT INITIALISATION
!

!      write(input,*) 'nangle',nangle(ipot)
! potentiel
! vcalc veut un angle en degree
! cosang is the angle
! conversion
rd=0.0d0
ipot=1
!...function to calculate potential for any given geometry
lmax=(mlv%nangle(ipot)-1)*mlv%isym(ipot)+mlv%mld(ipot)
vcalc=0.
!                            mld,0                  mld,0
!...calculate pvec: pvec(1)=d     (angle), pvec(2)=d       (angle), ...
!                            l                      l+isym
!
!   where isym = 2 if homonuclear molecule and isym =1 if heteronuclear
!   molecule
i=0
do l=mlv%mld(ipot),lmax,mlv%isym(ipot)
  i=i+1
  pvec(i)=dlmm0(l,mlv%mld(ipot),angl)
enddo
!   now obtain vector of vvl(nangle) coefficients and calculate energy.
!   note that vvl(1) = vv0 is the isotropic term if mld (mproj) is 0.
!
!                      -1       -1
!   v(r,angle,r)=pvec*p  *b(r)*r  *rvec = pvec * vvl * rvec
!
!   first built up b(r) matrix
iblk=0
do l=1,mlv%nr(ipot)
  do n=1,mlv%nangle(ipot)
    iblk=iblk+1
!...this function is defined as follows:
!                                          2                    maxpws
!     b(r)=exp(-a *(r-r1))*{ a  + a *r + a *r + ... + a         *r       }
!                1          2    3      4            maxpws+2
!
!                                       minmps                minmps+mpsstp
!       -(1/2)*(1+tanh((r-r0)/rref)){a         /r       + a          /r
!                            maxpws+3             maxpws+4
!
!                        maxmps
!       + . . . + a    /r       }
!                  npa
! r0 is set to 0
!         a(npa+2,iblk,ipot)=0.0
    tanhyb=(tanh((r_-mlv%a(ipot,mlv%npa(ipot)+2,iblk))/ &
mlv%a(ipot,mlv%npa(ipot)+1,iblk))+1.d0)*0.5d0
    fex=exp(-mlv%a(ipot,1,iblk)*(r_-mlv%a(ipot,mlv%npa(ipot)+3,iblk)))
    if(fex.gt.1.d10) fex=1.d10
    pfac=0.
    rpowi=1.
    do i=2,mlv%maxpws(ipot)
      pfac=pfac+mlv%a(ipot,i,iblk)*rpowi
      rpowi=rpowi*r_
    enddo
    pexfac=fex*pfac
    if(pexfac.gt.1.d20) pexfac=1.d20
    fmfac=0.
    jst=mlv%maxpws(ipot)+1
    if (.not. (jst.gt.mlv%npa(ipot))) then
      powmj=r_**mlv%minmps(ipot)
      do j=jst,mlv%npa(ipot)
        fmfac=fmfac+mlv%a(ipot,j,iblk)/powmj
        powmj=powmj*r_**mlv%mpsstp(ipot)
      enddo
    endif
    fct=pexfac-tanhyb*fmfac
    b(n,l)=fct
  enddo
enddo
!...fast loop if nr = 1 .and. novib
if (mlv%nr(ipot).gt.1) then
!...here if normal treatment of r - dependent potential
!                                             2
!...calculate rvec: rvec = { 1., (r-re), (r-re) , ... }
!                       -1     -1
!    and multiply with r  :   r  * rvec = rvec'
  do i=1,mlv%nr(ipot)
    rvecp(i)=0.
    dr=1.
    do j=1,mlv%nr(ipot)
      rvecp(i)=rvecp(i)+mlv%rinv(ipot,i,j)*dr
      dr=dr*rd
    enddo
  enddo
!  calculate  b(r) * rvec' = rvec''
  do i=1,mlv%nangle(ipot)
    rvecpp(i)=0.
    do j=1,mlv%nr(ipot)
      rvecpp(i)=rvecpp(i)+b(i,j)*rvecp(j)
    enddo
  enddo
!              -1
!  calculate  p  * rvec'' = vvl  ( note that vvl(1) = vv0 if mld (mproj) is 0 )
  do i=1,mlv%nangle(ipot)
    vvl(i)=0.
    do j=1,mlv%nangle(ipot)
      vvl(i)=vvl(i)+mlv%pinv(ipot,i,j)*rvecpp(j)
    enddo
  enddo
else
!                                                  -1
!...here if number of bond distances is 1:  vvl = p   * b(r) where
!   b(r) is a now column vector rather than a matrix
  do i=1,mlv%nangle(ipot)
    vvl(i)=0.
    do j=1,mlv%nangle(ipot)
      vvl(i)=vvl(i)+mlv%pinv(ipot,i,j)*b(j,1)
    enddo
  enddo
endif
!  sum over all terms: v(r,angle,r) = pvec(1)*vvl(1) + pvec(2)*vvl(2) +...
vvlsum=0.
do i=1,mlv%nangle(ipot)
  vvlsum=vvlsum+vvl(i)
  vvlp(i)=pvec(i)*vvl(i)
  vcalc=vcalc+vvlp(i)
enddo
v_=vcalc
end
!***************************************************************
!* set of subroutines/functions to be used with vfitdriver.f90 *
!***************************************************************
subroutine matinv(a,det,n,nc)
implicit none
integer*4 n,nc

integer*4 i,j,k,l
real*8 det
real*8 amax,saave
real*8 a(n,n)
integer*4 ik(100),jk(100)


10 det=1.
11 do 100 k=1,nc
amax=0.
21 do 30 i=k,nc
do 30 j=k,nc
23    if (dabs(amax)-dabs(a(i,j))) 24,24,30
24    amax=a(i,j)
ik(k)=i
jk(k)=j
30 continue
31 if (amax) 41,32,41
32 det=0.
goto 140
41 i=ik(k)
if (i-k) 21,51,43
43 do 50 j=1,nc
saave=a(k,j)
a(k,j)=a(i,j)
50    a(i,j)=-saave
51 j=jk(k)
if (j-k) 21,61,53
53 do 60 i=1,nc
  saave=a(i,k)
  a(i,k)=a(i,j)
60   a(i,j)=-saave
61 do 70 i=1,nc
  if (i-k) 63,70,63
63   a(i,k)=-a(i,k)/amax
70 continue
71 do 80 i=1,nc
  do 80 j=1,nc
    if (i-k) 74,80,74
74     if (j-k) 75,80,75
75     a(i,j)=a(i,j)+a(i,k)*a(k,j)
80 continue
81 do 90 j=1,nc
  if (j-k) 83,90,83
83   a(k,j)=a(k,j)/amax
90 continue
a(k,k)=1./amax
100 det=det*amax
101 do 130 l=1,nc
  k=nc-l+1
  j=ik(k)
  if (j-k) 111,111,105
105 do 110 i=1,nc
  saave=a(i,k)
  a(i,k)=-a(i,j)
110   a(i,j)=saave
111   i=jk(k)
if (i-k) 130,130,113
113    do 120 j=1,nc
 saave=a(k,j)
 a(k,j)=-a(i,j)
120    a(i,j)=saave
130 continue
140 return
end 

function dlmm0(l,m,theta)
implicit none
integer*4  l,m
real*8 theta
real*8 dlmm0
integer*4 i,l_
real*8 x,y,thedeg,pm1,pm2,pp,rat,ai,al,al2
real*8 pi
pi=acos(-1d0)
!
!...function to calculate dlmm0(cos(theta)) as defined in "brink and satchler"
!
thedeg=(theta*pi)/180.

!  if m>l pm1=0 !

if (m.gt.l) then
  pm1=0.
  dlmm0=0
  return
endif
x=cos(thedeg)
if (m.lt.0) then 
  write (6,*) '	NEGATIVE M IN LEGENDRE ROUTINE:	 ABORT'
  stop
endif
if (m.eq.0) then 
!  here for regular legendre polynomials
  pm1=1.
  pm2=0.
  do l_=1,l
  pp=((2.*l_-1.)*x*pm1-(l_-1.)*pm2)/float(l_)
  pm2=pm1
  pm1=pp
  enddo
else
!  here for alexander-legendre polynomials
  rat=1.
  do i=2,2*m,2
    ai=i
    rat=rat*((ai-1.)/ai)
  enddo
  y=sin(thedeg)
  pm1=sqrt(rat)*(y**m)
  pm2=0.
  do l_=m+1,l
    al=(l_+m)*(l_-m)
    al=1./al
    al2=((l_+m-1)*(l_-m-1))*al
    al=sqrt(al)
    al2=sqrt(al2)
    pp=(2.*l_-1.)*x*pm1*al-pm2*al2
    pm2=pm1
    pm1=pp
  enddo
! correct phase
  pm1=((-1.)**m)*pm1
endif
dlmm0=pm1
return
end

! ------------------------------------------------------------------------

