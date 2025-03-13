!*******************************************************************
!                                                                   *
!                 potential routines for using vfit                 *
!                                                                   *
!********************************************************************
!                         routines included:                        *
!                                                                   *
!  1 loapot (vcalc, testpt, pot) loads and calculates pot.     *
!                                                                   *
!********************************************************************
!
!  Note:  this pot routine requires a data file to be in hibxx/bin/progs/potdata:
!         filnam (cetermined in last line of input file)
!
#include "unused.h"
#include "common/syusr.F90"
#include "common/bausr.F90"
#include "common/ground.F90"
! --------------------------------------------------------------------------
subroutine loapot(iunit,filnam)
! --------------------------------------------------------------------------
!
!   subroutine to load potential provided from "vfit.prg"
!   author: b. follmeg
!   current revision date: 6-may-1997 by mha
!
!   on entry: iunit -> unit assigned to file containing the potential
!             filnam -> filename
!
!   variables in common blocks:
!   the settings used in the scattering program
!
!
!   nc         -> array containing the number of coefficients for the
!                 r - dependent potential for each term
!   a          -> array containing the coefficients for each term
!   maxpw      -> array containing the maximum exponent in powerseries
!                 expansion for each term
!   minms      -> array containig minimum exponent in multipolseries
!                 expansion for each term
!
!   maxms      -> array containig maximum exponent in multipolseries
!                 expansion for each term
!
!   msstp      -> array containing step in multipol expansion
!
!   mdimp      -> array containing the dimensions of the transformation
!                 matrices used to fit the angular dependent potential
!
!   pinv       -> the transformation matrices for each term
!
!   mdimr      -> array containing the length of the vectors vibmat,
!                 where the elements of vibmat are vibrational matrix
!                 elements
!
!   ntv        -> array containing the number of different vectors
!                 for each term
!
!   vibmat     -> array containing vectors of vibrational matrix elements
!
!
!   re         -> equilibrium bond distance of the molecule used in the
!                 calculation of the vibrational matrix elements
!
!  the common block /cofit/ is used to pass arguments
!
! --------------------------------------------------------------------------
use mod_cosysi, only: ispar
use mod_par, only: readpt
use mod_parbas, only: maxtrm, maxvb2, ntv, ivcol, ivrow, lammin, lammax, mproj
use mod_parpot, only: potnam=>pot_name
use mod_skip, only: nskip, iskip
use mod_hiiolib1, only: openf
implicit double precision(a-h,o-z)
character*(*) filnam
character*80 potlab
character*68 filnm1
common /coptx/ nblkx,maxpwx,minmpx,maxmpx,mpsstx,junk, &
  rex(20),rin(20),tanhy(250),fex(250),a1(250),ah(250), &
  art(250),are(250)
#include "common/parvfit.F90"
common /copot/ nc(maxtrm),a(20,50,maxtrm),maxpw(maxtrm), &
               minms(maxtrm),maxms(maxtrm),msstp(maxtrm), &
               mdimp(maxtrm),pinv(maxang,maxang,maxtrm), &
               mdimr(maxtrm),vibmat(maxnr,maxvb2,maxtrm), &
               avec(maxang)
common /cofit/ npa,maxpws,minmps,maxmps,mpsstp,idimp,idimr
common /core/  re
integer, pointer :: nterm
nterm=>ispar(1)
potnam='WERNER-FOLLMEG VFIT'
! default data
lammin(1)=2
lammax(1)=2
mproj(1)=0
ntv(1)=1
ivcol(1,1)=0
ivrow(1,1)=0
!      ibasty=4
! return at this point if readpt is false
if (.not.readpt) return
filnm1 = 'potdata/'//filnam
call openf(iunit,filnm1,'sf',0)
!aber
rewind iunit
!aber
read(iunit,'(a)',err=888,iostat=ierr) potlab
read(iunit,*,err=888,iostat=ierr) nterm, iskip, re
nblkx=0
maxpwx=0
maxmpx=0
minmpx=100
mpsstx=100
if(nterm.gt.maxtrm) stop 'nterm'
do 100 it = 1, nterm
   read(iunit,*,err=888,iostat=ierr) maxpw(it),minms(it), &
     maxms(it),msstp(it),nc(it),nblk,mdimp(it),lammin(it), &
     lammax(it),mproj(it),mdimr(it),ntv(it)
   maxpwx=max(maxpwx,maxpw(it))
   if(maxms(it).ne.0) maxmpx=max(maxmpx,maxms(it))
   if(minms(it).ne.0) minmpx=min(minmpx,minms(it))
   if(msstp(it).ne.0) mpsstx=min(mpsstx,msstp(it))
   if(nblk.gt.50) stop 'nblk'
   if(nc(it)+3.gt.20) stop 'nc'
   do 30 iblk = 1, nblk
   nblkx=nblkx+1
   read(iunit,*,err=888,iostat=ierr) (a(i,iblk,it),i=1,nc(it)+3)
   a1(nblkx)=a(1,iblk,it)
   art(nblkx)=a(nc(it)+2,iblk,it)
   are(nblkx)=a(nc(it)+3,iblk,it)
30    ah(nblkx)=1.0d0/a(nc(it)+1,iblk,it)
   nangle = mdimp(it)
   if(nangle.gt.maxang) stop 'nangle'
   if(ntv(it).gt.maxvb2) stop 'maxvb2'
   if(mdimr(it).gt.maxnr) stop 'maxnr'
   do 40 i = 1, nangle
40    read(iunit,*,err=888,iostat=ierr) (pinv(i,j,it),j=1,nangle)
   do 50 j = 1, ntv(it)
   read(iunit,*,err=888,iostat=ierr) ivrow(j,it),ivcol(j,it)
50    read(iunit,*,err=888,iostat=ierr)(vibmat(i,j,it),i=1,mdimr(it))
100 continue
close(iunit)
write(6,105) filnam,potlab(:72)
105 format(/,' POTENTIAL LOADED FROM FILE ',(a),//,' LABEL: ',(a))
potlab=' '
if(iskip.eq.1) then
  potlab=' Molecule is heteronuclear'
else if(iskip.eq.2) then
  potlab=' Molecule is homonuclear'
else
  potlab=' *** INVALID SYMMETRY, SOMETHING IS WRONG'
end if
nskip=iskip
write(6,110) potlab(:73),nterm
110 format((a),/,' NUMBER OF TERMS (NTRM) READ IN: ',i2,/)
write(6,115)
115 format('  NTRM    LMMIN    LMMAX    MLD    NVBLOCKS')
do 150 i=1,nterm
write(6,120) i,lammin(i),lammax(i),mproj(i),ntv(i)
120 format(t4,i2,t12,i2,t21,i2,t29,i2,t36,i2)
150 continue
return
888 write(6,300) ierr
300 format(' *** READ ERROR IN LOAPOT, IOSTAT: ',i4,', ABORT ***')
close(iunit)
return
end

subroutine driver
use mod_covvl, only: vvl
use mod_par, only: readpt
implicit double precision (a-h,o-z)
character *48 potnam
character *40 filnam
readpt=.true.
potnam='WERNER-FOLLMEG VFIT'
print *, 'potential subroutine:  ',potnam
print *, 'filename ?  '
read (5, *)filnam
call loapot(1,filnam)
2 print *, ' r (bohr) '
read (5, *, end=93) r
if (r .lt. 0.d0) go to 93
call pot(vv0,r)
write (6, 100) vvl
100 format(' vvl(1:7) vsigma',7(1pe16.8),/, &
       ' vvl(8:14)   vpi',7(1pe16.8),/, &
       ' vvl(15:19)   v1',5(1pe16.8),/, &
       ' vvl(20:24)   v2',5(1pe16.8))
goto 2

93 do j=1,30
   r=(j-1)*.25+4
   call pot(vv0,r)
   write (6,101) r,vv0,vvl
101    format(f5.2,30(1pg10.2))
enddo
end
! -------------------------------------------------------------------
subroutine pot(vv0,r)
! -------------------------------------------------------------------
!  subroutine to calculate the r-dependent coefficients
!
!  author: b. follmeg
!  current revision date: 18-aug-1998 by mha
! -------------------------------------------------------------------
use mod_covvl, only: vvl
use mod_cosysi, only: ispar
use mod_parbas, only: maxtrm, maxvb2, ntv, lammin, lammax
use mod_selb, only: ibasty
use mod_skip, only: nskip
implicit double precision(a-h,o-z)
#include "common/parvfit.F90"
common /copot/ nc(maxtrm),a(20,50,maxtrm),maxpw(maxtrm), &
               minms(maxtrm),maxms(maxtrm),msstp(maxtrm), &
               mdimp(maxtrm),pinv(maxang,maxang,maxtrm), &
               mdimr(maxtrm),vibmat(maxnr,maxvb2,maxtrm), &
               avec(maxang)
common /cofit/ npa,maxpws,minmps,maxmps,mpsstp,idimp,idimr
common /coptx/ nblkx,maxpwx,minmpx,maxmpx,mpsstx,junk, &
  rex(20),rin(20),tanhy(250),fex(250),a1(250),ah(250), &
  art(250),are(250)
common /core/  re
integer, pointer :: nterm
nterm=>ispar(1)
rr = r
vv0 = 0.d0
jj = 0
do 10 i=1,maxpwx-1
10 rex(i)=rr**(i-1)
ri=1.0d0/rr
do 20 i=minmpx,maxmpx
20 rin(i)=ri**i
do 30 i=1,nblkx
  fex(i)=exp(-a1(i)*(rr-are(i)))
!        if (art(i) .eq. 0.d0) then
! here for follmeg form: long range part just multiplied by
! tanh(R)
!          tanhy(i)=tanh(rr)
!        else
! here for Berning form:  long range part multiplied as in
! degli-esposti and werner
    tanhy(i)=(1.d0+tanh(ah(i)*(rr-art(i))))*0.5d0
!        endif
30 continue

!
! loop over all terms starts here
!
iblkxx=0
ishift=0
do 100 iterm = 1, nterm
   npa = nc(iterm)
   maxpws = maxpw(iterm)
   minmps = minms(iterm)
   maxmps = maxms(iterm)
   mpsstp = msstp(iterm)
   idimp = mdimp(iterm)
   idimr = mdimr(iterm)
   lmin = lammin(iterm)
   lmax = lammax(iterm)
   ivmax = ntv(iterm)
!
! loop over all vibrational terms
!
   do 80 iv = 1, ivmax
   iblkx=iblkxx
     call vcalc(rr,a(1,1,iterm),pinv(1,1,iterm), &
                vibmat(1,iv,iterm),avec,iblkx)
     ii = 0
!
! store isotropic term in vv0 in certain cases
!
     if (ibasty .ne. 4 .and. ibasty .ne. 7 &
        .and. ibasty .ne. 10) then
       if(iterm .eq.1 .and. iv .eq. 1) &
       then
         vv0 = avec(1)
         ii = 1
       else
         if (ishift.eq.0) then
! shift jj down by one
            jj=jj-1
            ishift=1
         endif
       end if
     endif
!
! store anisotropic terms in vvl array
!
     do 50 l = lmin, lmax, nskip
       ii = ii + 1
       jj = jj + 1
       vvl(jj) = avec(ii)
!             print *, l,jj,vvl(jj)
50     continue
80   continue
iblkxx=iblkx
100 continue
return
end
! --------------------------------------------------------------------
subroutine vcalc(rr,a,pinv,rvecp,avec,iblkx)
! --------------------------------------------------------------------
!
! author: b.follmeg
! current revision date: 26-may-1991 by mha
!
! --------------------------------------------------------------------
use mod_himatrix, only: mxva
implicit double precision(a-h,o-z)
parameter (maxang=10)
dimension pinv(maxang,maxang),a(20,50),rvecp(10),avec(10)
common /cofit/ npa,maxpws,minmps,maxmps,mpsstp,idimp,idimr
!
! the arrays in common block /cowork/ are used as scratch arrays
!
common /cowork/ pvec(10),rvecpp(10),b(10,10)
common /coptx/ nblkx,maxpwx,minmpx,maxmpx,mpsstx,junk, &
  rex(20),rin(20),tanhy(250),fex(250),a1(250),ah(250), &
  art(250),are(250)
!
!    obtain vector of a(r) coefficients:
!
!               -1          -1          -1
!              p  * b(r) * r  * rvec = p  * b(r) * rvec' = a(r)
!
!                       -1
!   the vector rvec' = r  * rvec must be provided by vfit !
!
!   first built up b(r) matrix
!
UNUSED_DUMMY(rr)
iblk = 0
ne=maxpws-1
m1=maxpws+1
nm=0
if(mpsstp.ne.0.and.maxmps.ne.0.and.minmps.ne.0) &
      nm=(maxmps-minmps)/mpsstp+1
do 200 l = 1, idimr
do 100 n = 1, idimp
iblk =iblk +1
iblkx=iblkx+1
b(n,l)=fex(iblkx)*ddot(ne,a(2,iblk),1,rex,1)
if(nm.ne.0) b(n,l)=b(n,l) &
    -tanhy(iblkx)*ddot(nm,a(m1,iblk),1,rin(minmps),mpsstp)
100 continue
200 continue
!
! fast loop if idimr = 1 and rvecp(1) = 1
!
if((idimr.gt.1).or.(rvecp(1).ne.1.)) then
!
!  calculate  b(r) * rvec' = rvec''
!
  call mxva(b,1,10,rvecp,1,rvecpp,1,idimp,idimr)
!              -1
!  calculate  p  * rvec'' = a(r)
!
  call mxva(pinv,1,maxang,rvecpp,1,avec,1,idimp,idimp)
!                                                  -1
!  here if number of bond distances is 1: avec = p   * b(r) where
!  b(r) is a now column vector rather than a matrix
!
else
  call mxva(pinv,1,maxang,b,1,avec,1,idimp,idimp)
!
end if
!
!  avec now contains
!  avec(1) = vv0, avec(2) = vvl(1), ... ( if iterm=1 and lammin=0) or
!  avec(1) = vvl(1), avec(2) = vvl(2), ... ( if iterm > 1 or lammin > 0)
!
return
end
