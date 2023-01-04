module mod_hibrid4
contains
#include "assert.h"
!***********************************************************************
!                                                                       *
!                         hibridon 4  library                           *
!                                                                       *
!************************************************************************
!                          routines included:                           *
!                                                                       *
!   2. sprint         prints s-matrices on the screen                   *
!   6. turn           function, determines classical turning point      *
!   9. waverd         writes and reads header file for wavefunction     *
!  10. psi            to determine wavefunction
!  11. flux           to determine fluxes
!  12. transmt        print out transformation matrix at rout
!************************************************************************
! ----------------------------------------------------------------------
subroutine sprint (fname, ia)
!   reads and prints s-matrices
!   author: h.-j. werner
!   now print out j12 for molecule-molecule systems
!   and for 2P atom + molecule system
!   current revision date: 19-jun-2015 by p.dagdigian
! ----------------------------------------------------------------------
use mod_cosout, only : nnout, jout
use mod_coj12p, only: j12pk
use constants
use mod_codim, only: nmax => mmax
use mod_cojq, only: jq ! jq(1)
use mod_colq, only: lq ! lq(1)
use mod_coinq, only: inq ! inq(1)
use mod_coinhl, only: jlev => inhold ! jlev(1)
use mod_coisc1, only: inlev => isc1 ! inlev(1)
use mod_cosc1, only: elev => sc1 ! elev(1)
use mod_coz, only: sreal => z_as_vec ! sreal(1)
use mod_cozmat, only: simag => zmat_as_vec ! simag(1)
use mod_hibasis, only: is_j12
use mod_parpot, only: potnam=>pot_name, label=>pot_label
use mod_selb, only: ibasty
use mod_hiutil, only: gennam
use mod_hismat, only: smatread, rdhead
use mod_hitypes, only: packed_base_type
implicit none
character*(*), intent(in) :: fname
integer, intent(in) :: ia(4)

real(8) :: ered, rmu
integer :: i, iadr, ienerg, ierr, ij
integer :: j, ja, je, jfinal, jfirst, jj1, jj2, jlp, jlpar, jtot, jtota, jtotb, jtotd
integer :: lenx, ncol, nlevel, nlevop, nopen, nu,nud, numax, numin

type(packed_base_type) :: packed_base
character*20 cdate
character*40 xname
logical  existf, csflag, flaghf, flagsu, twomol, nucros

!
!.....jtota: first jtot to be printed
!.....jtotb: last jtot to be printed
!.....jlp:   jlpar to be printed. if zero, both parities are searched
!.....ienerg: number of energy
jtota = ia(1)
jtotb = ia(2)
jlp = ia(3)
ienerg = ia(4)
if(ienerg.eq.0) ienerg = 1
if(jtotb.eq.0) jtotb = jtota
call gennam(xname, fname, ienerg, 'smt', lenx)
inquire (file = xname,  exist = existf)
if (.not. existf) then
  write (6,  20) xname(1:lenx)
20   format(/' FILE ', (a), ' NOT FOUND')
  return
end if
call openf(1, xname, 'tu', 0)
call rdhead(1,cdate, &
   ered,rmu,csflag,flaghf,flagsu,twomol, &
   nucros,jfirst,jfinal,jtotd,numin,numax,nud,nlevel,nlevop,nnout, &
   jlev,inlev,elev,jout)
  if(csflag) jlp=0
  write (6, 70) xname, cdate, label, potnam, ered*econv
70   format(/' S-MATRICES READ FROM FILE ', (a)/ &
          ' WRITTEN: ', (a), /, ' LABEL:   ', (a),/, &
          ' POT NAME:  ', (a) &
         /,' ENERGY: ', f10.3)
  if (.not. twomol) then
    write (6, 80)
80     format ('   N   J   INDEX   EINT(CM-1)')
    do  95  i  =  1,  nlevop
      if (ibasty.ne.12) then
        if (.not. flaghf) then
          write (6,  85) i,  jlev(i),  inlev(i),  elev(i) * econv
85           format (i4,  i5,  i6,  f11.3)
        else
          write (6,  90) i,  (jlev(i)+0.5d0),  inlev(i), &
                           elev(i) * econv
90           format (i4,  f5.1,  i6,  f11.3)
        end if
      else
        write (6,91) i, jlev(i),inlev(i)+0.5,elev(i)*econv
91         format (i4,i5,f6.1,f11.3)
      endif
95     continue

  else
    write (6, 100)
100     format('   N   J1   J2  INDEX  EINT(CM-1)')
      do 135  i  =  1,  nlevop
        jj1 = jlev(i) / 10
        jj2 = mod(jlev(i), 10)
        write (6,  130) i,  jj1,  jj2,  inlev(i), &
                        elev(i) * econv
130         format (i4,  2i5,  i6,  f11.3)
135       continue
  end if
jtotb = min0(jtotb, jfinal)
!
iadr=0
30 nopen = 0
call smatread (iadr, sreal, simag, jtot, jlpar, nu, &
                  jq, lq, inq, packed_base, &
                   1, nmax, nopen, ierr)
if(csflag) jlpar=0
if(ierr.eq.-1) goto 400
if(ierr.lt.-1) then
  write(6,35) xname
35   format(' ERROR READING FILE ',(a))
  goto 400
end if
!.....assume that jlpar=1 is stored first
220 if (jlpar .eq.1.and.jlp.eq.-1) goto 30
if(jtot.lt.jtota) goto 30
if(jtot.gt.jtotb) then
   if(jlp.eq.jlpar.or.jlp.eq.-1) goto 400
   jlp=-1
   goto 220
end if
if (.not.flaghf) then
  write (6, 230) jtot, jlpar, nu, nnout
230   format(/' JTOT=', i4, '  JLPAR=', i2, '  NU=', i3, &
          '  NNOUT=', i3)
else
  write (6, 240) jtot+0.5d0, jlpar, nu+0.5d0, nnout
240   format(/' JTOT=', f6.1, '  JLPAR=', i2, '  NU=', f5.1, &
          '  NNOUT=', i3)
end if
if (nnout.le.0) then
  write (6, 250)
250   format(/' COLUMN INDICES:')
  write (6, 290) 'N    ', (j, j=1, nopen)
  if (.not. twomol) then
    if (flaghf) then
      write (6, 260) 'J    ', (jq(j)+0.5d0, j=1, nopen)
260       format (1x, (a), (t10, 20(f6.1)) )
    else
      write (6, 290) 'J    ', (jq(j), j=1, nopen)
    end if
  else
    write (6, 270) 'J1/J2', (jq(j), j=1, nopen)
270     format(1x, (a), (t10, 10i6) )
  end if
  write (6, 290) 'L    ', (lq(j), j=1, nopen)
  write (6, 290) 'INDEX', (inq(j), j=1, nopen)
end if
  write (6, 280)
280   format(/' ROW INDICES:')
  write (6, 290) 'N    ', (j, j=1, packed_base%length)
  if (.not. is_j12(ibasty)) then
    if (flaghf) then
      write (6, 260) 'J    ', (packed_base%jpack(j)+0.5d0, j=1, packed_base%length)
    else
      write (6, 290) 'J    ', (packed_base%jpack(j), j=1, packed_base%length)
    end if
    write (6, 290) 'IS   ', (packed_base%inpack(j), j=1, packed_base%length)
  else
    if (ibasty.eq.12 .or. ibasty.eq.15) then
      ASSERT(.false.)
      write (6, 290) 'J    ', (packed_base%jpack(j), j=1, packed_base%length)
      write (6, 260) 'JA   ', (packed_base%inpack(j)+0.5d0, j=1, packed_base%length)
      write (6, 260) 'J12  ', (j12pk(j)+0.5d0, j=1, packed_base%length)
    else
      write (6, 270) 'J1/J2', (packed_base%jpack(j), j=1, packed_base%length)
      write (6, 270) 'J12', (j12pk(j), j=1, packed_base%length)
      write (6, 270) 'IS', (packed_base%inpack(j), j=1, packed_base%length)
    endif
  end if
  write (6, 290) 'L    ', (packed_base%lpack(j), j=1, packed_base%length)
290   format(1x, (a), (t10, 10i6))
ncol = nopen
if(nnout.gt.0) ncol = packed_base%length
write (6, 300) 'REAL PART OF THE S-MATRIX'
300 format(/1x, (a))
do 330 ja = 1, packed_base%length, 10
  je=min0(ja+9,packed_base%length)
  write (6, 310) (j, j = ja, je)
310   format(10i12)
  ij=1-nmax
  do 320 i=1,ncol
    if(nnout.gt.0) je=min0(ja+9,i)
    if (je.ge.ja) write (6, 315) i, (sreal(ij+j*nmax), j=ja,je)
315     format(1x, i3, 10(1pe12.4))
320   ij=ij+1
330 continue
write (6, 300) 'IMAGINARY PART OF THE S-MATRIX'
do 350 ja = 1, packed_base%length, 10
  je=min0(ja+9,packed_base%length)
  write (6, 310) (j, j = ja, je)
  ij=1-nmax
  do 340 i = 1, ncol
    if(nnout.gt.0) je=min0(ja+9,i)
    if (je.ge.ja) write (6, 315) i, (simag(ij+j*nmax), j=ja,je)
340   ij=ij+1
350 continue
goto 30
400 call closf(1)
close(1)
return
end
!
!     ------------------------------------------------------------------
function iwavsk(irecr)
!     ------------------------------------------------------------------
!     Function to return offset of wfu file for recrod #irec (stream IO)
!
!     Written by: Qianli Ma
!     Latest revision: 20-apr-2012
!
!     This function needs nchwfu, ipos2 and ipos3 from the mod_wave
!     module.  These variables are set by waverd.
!
!     The stream IO counterpart for `dbrr(,,,irec)` is `read
!     (,pos=wavesk(irec))...
!     ------------------------------------------------------------------
use mod_wave, only: nchwfu, ipos2, ipos3, nrlogd, inflev, get_wfu_rec1_length, get_wfu_logd_rec_length, get_wfu_airy_rec_length
implicit none
integer, intent(in) :: irecr
integer(8) :: iwavsk

!     The following variables are for size-determination of (machine
!     dependent) built-in types
integer int_t
double precision dble_t
character char_t

integer(8) :: lr1  ! length of record 1 in bytes
integer(8) :: lrlogd  ! length of a logd record
integer(8) :: lrairy  ! length of an airy record
integer, parameter :: char_size = int(sizeof(char_t), kind(int_t))
integer, parameter :: int_size = int(sizeof(int_t), kind(int_t))
integer, parameter :: dbl_size = int(sizeof(dble_t), kind(int_t))
!
if (irecr .le. 0) then
   iwavsk = -1
   goto 100
end if
if (irecr .eq. 1) then
   iwavsk = 1
   goto 100
end if
if (irecr .eq. 2) then
   iwavsk = ipos2
   goto 100
end if
if (irecr .eq. 3) then
   iwavsk = ipos3
   goto 100
end if
!     Length for record 1 (at the beginning of the file)
lr1 = get_wfu_rec1_length(nchwfu)
!     Length for each block written by the Airy and LOGDpropagator
lrairy = get_wfu_airy_rec_length(nchwfu, inflev)
lrlogd = get_wfu_logd_rec_length(nchwfu, inflev)
!
if ((irecr - 3) .le. nrlogd) then
!     within the logd range of the file
   iwavsk = lr1 + (irecr - 4) * lrlogd + 1
   goto 100
else
!     airy range of the file
   iwavsk = lr1 + nrlogd * lrlogd &
        + (irecr - 4 - nrlogd) * lrairy + 1
   goto 100
end if
!     should never reach here if function called properly
write (0, *) '*** OOPS! ERROR SEEKING WFU FILE. ABORT.'
call exit()
100 continue
ASSERT(iwavsk > 0)
end
!
! -------------------------------------------------------------------------
subroutine wavewr(jtot,jlpar,nu,nch,rstart,rendld)
! -------------------------------------------------------
!  subroutine to write initial header information on wavefunction file
!  (file jobname.WFU, logical unit 22), unit is opened in subroutine openfi
!     written by:  millard alexander
!     latest revision:  11-dec-2011
!     common blocks amat and bmat are used to store real and
!     imaginary parts of asymptotic wavefunction (only used in
!     read of wavefunction from saved file)
!
!     Major revision: 16-mar-2012 by Q. Ma
!     Use stream I/O for smaller file size and better compatibility
!
!     current revision: 20-apr-2012 by q. ma
!
! -------------------------------------------------------
#define AMAT_AS_VEC_METHOD_DISTINCT 1
#define AMAT_AS_VEC_METHOD_POINTER 2
#define AMAT_AS_VEC_METHOD_NOVEC 3
#define AMAT_AS_VEC_METHOD AMAT_AS_VEC_METHOD_DISTINCT
use mod_coeint, only: eint
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_POINTER)
use, intrinsic :: ISO_C_BINDING
use mod_coamat, only: amat ! amat(25)
#endif
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_NOVEC)
use mod_coamat, only: amat ! amat(25)
#endif
use mod_cojq, only: jq ! jq(1)
use mod_colq, only: lq ! lq(1)
use mod_coinq, only: inq ! inq(1)
use mod_par, only: csflag, flaghf, wrsmat, photof
use funit
use mod_wave, only: irec, ifil, nchwfu, ipos2, ipos3, nrlogd, iendwv, get_wfu_rec1_length, wfu_format_version
use mod_parpot, only: potnam=>pot_name, label=>pot_label
use mod_ered, only: ered, rmu
use mod_hiutil, only: dater
implicit none
integer, intent(in) :: jtot
integer, intent(in) :: jlpar
integer, intent(in) :: nu
integer, intent(in) :: nch
real(8), intent(in) :: rstart
real(8), intent(in) :: rendld

character*20 :: cdate

integer :: i
integer(8) :: end_of_rec1_pos
!
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_POINTER)
real, pointer :: amat_as_vec(:)
#endif

ifil = FUNIT_WFU

nchwfu = nch
ipos2 = -1
ipos3 = -1
nrlogd = 0
!     Mark the position of the EOF of the WFU file in order to by pass
!     (likely) a bug in the intel compiler that INQUIRE does not return
!     the proper offset
iendwv = 1
!     Write magic number
write (ifil, err=950) char(128), 'WFU'

if (wrsmat) then
   write (ifil, err=950) char(0), wfu_format_version, char(0), char(0)
else
   write (ifil, err=950) char(1), wfu_format_version, char(0), char(0)
end if
!
write (ifil, err=950) ipos2, ipos3, nrlogd
call dater(cdate)
write (ifil, err=950) cdate, label, potnam
!     Four zero-bytes for alignment / C struct compatibility
write (ifil, err=950) char(0), char(0), char(0), char(0)
!
write (ifil, err=950) jtot, jlpar, nu, nch, csflag, flaghf, photof
write (ifil, err=950) ered, rmu, rstart, rendld
!
write (ifil, err=950) (jq(i), i=1, nch), (lq(i), i=1, nch), &
     (inq(i), i=1, nch)
write (ifil, err=950) (eint(i), i=1, nch)
!
write (ifil, err=950) 'ENDWFUR', char(1)
inquire(ifil, pos=end_of_rec1_pos)
ASSERT(end_of_rec1_pos == (get_wfu_rec1_length(nchwfu) + 1))
!
iendwv = iendwv + get_wfu_rec1_length(nchwfu)
irec=3  ! irec=2 and irec=3 are reserved and their position in the file are stored in ipos2 and ipos3
return

950 write (0, *) '*** ERROR WRITING WFU FILE. ABORT.'
call exit
return

end subroutine wavewr
!
!     ------------------------------------------------------------------
!     reads header file for wavefunction (wfu file)
subroutine waverd(jtot,jlpar,nu,nch,npts,nopen,nphoto,jflux, &
     rstart,rendld,rinf)
use mod_wave, only: irec, ifil, nchwfu, ipos2, ipos3, nrlogd, inflev, get_wfu_rec1_length, wfu_format_version
use mod_coeint, only: eint
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_DISTINCT)
use mod_coamat, only: amat => psir ! amat(25) psir(nopen, nopen)
#endif
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_POINTER)
use, intrinsic :: ISO_C_BINDING
use mod_coamat, only: amat ! amat(25)
#endif
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_NOVEC)
use mod_coamat, only: amat ! amat(25)
#endif
use mod_cobmat, only: bmat => psii ! bmat(25), here bmat is used as a vector 
use mod_cotq1, only: dpsir ! dpsir(25)
use mod_cotq2, only: dpsii ! dpsii(25)
use mod_cojq, only: jq ! jq(1)
use mod_colq, only: lq ! lq(1)
use mod_coinq, only: inq ! inq(1)
use mod_coisc1, only: isc1 ! isc1(25)
use mod_cosc1, only: sc1 ! sc1(10)
use mod_cosc2, only: sc2 ! sc2(10)
use mod_cosc3, only: sc3 ! sc3(10)
use mod_cosc4, only: sc4 ! sc4(10)
use mod_cosc5, only: sc5 ! sc5(10)
use mod_cow, only: w => w_as_vec ! w(25)
use mod_cozmat, only: zmat => zmat_as_vec ! zmat(25)
use mod_par, only: csflag, flaghf, photof
use funit
use mod_parpot, only: potnam=>pot_name, label=>pot_label
use mod_ered, only: ered, rmu
use mod_hivector, only: dset
implicit none
integer, intent(out) :: jtot
integer, intent(out) :: jlpar
integer, intent(out) :: nu
integer, intent(out) :: nch
integer, intent(out) :: npts
integer, intent(out) :: nopen
integer, intent(out) :: nphoto
integer, intent(out) :: jflux
real(8), intent(out) :: rstart
real(8), intent(out) :: rendld
real(8), intent(out) :: rinf


character*48 :: oldlab, oldpot
character*20 :: olddat

character :: csize8(8), csize4(4)
integer :: i
integer :: nopsq
integer :: nrecs
real(8), parameter :: zero=0.d0
! integer, parameter :: izero=0

!
ifil = FUNIT_WFU ! the wfu file is expected to be open using this unit
!     Read the magic number (from the start of the file)
read (ifil, pos=1, end=900, err=950) csize8
inflev = ichar(csize8(5))
if (csize8(6) /= wfu_format_version) then
  write (0,'(a,i3,a,i3,a)') '*** UNHANDLED VERSION OF WFU FORMAT : ', ichar(csize8(6)), ' (THIS VERSION OF HIBRIDON ONLY HANDLES WFU FORMAT VERSION ',  ichar(wfu_format_version),'). ABORT.'
  call exit()
end if
!
read (ifil, end=900, err=950) ipos2, ipos3, nrlogd
!
read (ifil, end=900, err=950) olddat, oldlab, oldpot
label = oldlab
potnam = oldpot
!     Read four zero bytes
read (ifil, end=900, err=950) csize4
!
read (ifil, end=900, err=950) jtot, jlpar, nu, nch, csflag, &
     flaghf, photof
!     nchwfu is used in locating the position for records
nchwfu = nch
read (ifil, end=900, err=950) ered, rmu, rstart, rendld
write (6, 245) olddat
if (jflux .ne. 0) write (3, 245) olddat
if (jflux .eq. 0) write (2, 245) olddat
245 format('    FROM CALCULATION ON: ',(a))
if (jflux .ne. 0) write (3, 250) oldlab
if (jflux .eq. 0) write (2, 250) oldlab
write (6, 250) oldlab
250 format('    INITIAL JOB LABEL: ', (a))
if (jflux .ne. 0) write (3, 251) oldpot
if (jflux .eq. 0) write (2, 251) oldpot
write (6, 251) oldpot
251 format('    INITIAL POT NAME: ', (a))
!
!     Read in channel labels
read (ifil, end=900, err=950) (jq(i), i=1, nch), &
     (lq(i), i=1, nch), (inq(i), i=1, nch), &
     (eint(i), i=1, nch)
!
! start reading in information from record 2 here
read (ifil, end=900, err=950, pos=iwavsk(2)) nrecs, nopen, &
     nphoto, rinf
npts = nrecs - 3
! read in wavevectors, bessel functions j, j', n, n'
! first initialize to zero for all channels
call dset(nch,zero,sc1,1)
call dset(nch,zero,sc2,1)
call dset(nch,zero,sc3,1)
call dset(nch,zero,sc4,1)
call dset(nch,zero,sc5,1)
read (ifil, end=900, err=950) (sc1(i), i=1, nopen), &
     (sc2(i), i=1, nopen), (sc3(i), i=1, nopen), &
     (sc4(i), i=1, nopen), (sc5(i), i=1, nopen)
nopsq = nopen ** 2
! read in sreal and simag, store in w and zmat
read (ifil, end=900, err=950) (w(i), i=1, nopsq), &
     (zmat(i), i=1, nopsq)
if (photof) then
! read in number of initial photodissociation states
!        call dbri(mphoto,1,ifil,REC_LAST_USED)
!        nphoto=mphoto
! read in real part of photodissociation amplitude
! overlay sreal which is not needed for photodissociation problem
   read (ifil, end=900, err=950) (w(i), i=1, nphoto * nopen)
! read in imaginary part of photodissociation amplitude
! overlay simag which is not needed for photodissociation problem
   read (ifil, end=900, err=950) (zmat(i), i=1, nphoto * nopen)
endif
! read in channel packing list and real and imaginary parts
! of scattering wavefunction and derivative
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_DISTINCT)
read (ifil, end=900, err=950, pos=iwavsk(3)) &
     (isc1(i), i=1, nopen), (amat(i), i=1, nopsq), &
     (bmat(i), i=1, nopsq), (dpsir(i), i=1, nopsq), &
     (dpsii(i), i=1, nopsq)
#endif
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_POINTER)
! amat_as_vec is a view of the matrix amat(nopen, nopen) as a vector(nopen*nopen)
call C_F_POINTER (C_LOC(amat), amat_as_vec, [nopsq])
read (ifil, end=900, err=950, pos=iwavsk(3)) &
     (isc1(i), i=1, nopen), (amat_as_vec(i), i=1, nopsq), &
     (bmat(i), i=1, nopsq), (dpsir(i), i=1, nopsq), &
     (dpsii(i), i=1, nopsq)
#endif
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_NOVEC)
read (ifil, end=900, err=950, pos=iwavsk(3)) &
     (isc1(i), i=1, nopen), amat, &
     (bmat(i), i=1, nopsq), (dpsir(i), i=1, nopsq), &
     (dpsii(i), i=1, nopsq)
#endif
irec = 3
return
!
900 continue
950 write (0, *) '*** ERROR READING WFU FILE. ABORT.'
call exit
return
!
end
! ------------------------------------------------------------------
subroutine psi(filnam,a)
!
! driver subroutine to calculate scattering wavefunction and fluxes
! from information stored in direct access file
!
! author: millard alexander
! current revision date (algorithm): 15-apr-1997 by mha
! revised on 30-mar-2012 by q. ma for stream I/O of wfu files
! current revision: 20-apr-2012 by q. ma
!
! special version for 13p collisions
!
! ------------------------------------------------------------------
use mod_cosout, only: nnout, jout
use mod_coiout, only: niout, indout
use constants
use mod_coqvec, only: nphoto
use mod_coeint, only: eint
use mod_coamat, only: psir ! psir(100) psir(nopen,nopen)
use mod_cobmat, only: psii ! psii(100) 
use mod_cotq1, only: dpsir ! dpsir(100)
use mod_cotq2, only: dpsii ! dpsii(100)
use mod_cojq, only: jq ! jq(60)
use mod_colq, only: lq ! lq(10)
use mod_coinq, only: inq ! inq(60)
use mod_coisc1, only: ipack => isc1 ! ipack(10)
use mod_coisc2, only: nlist => isc2 ! nlist(50)
use mod_coisc3, only: nalist => isc3 ! nalist(60)
use mod_coisc5, only: nblist  => isc5   ! nblist(60)
use mod_cosc2, only: fj  => sc2   ! fj(10)
use mod_cosc3, only: fjp => sc3   ! fjp(10)
use mod_cosc4, only: fn  => sc4   ! fn(10)
use mod_cosc5, only: fnp => sc5   ! fnp(10)
use mod_cosc6, only: sc  => sc6   ! sc(100)
use mod_cosc7, only: sc1  => sc7   ! sc1(100)
use mod_cosysi, only: ispar
use mod_coz, only: scmat => z_as_vec ! scmat(100)
use mod_cow, only: sr => w_as_vec ! sr(100)
use mod_cozmat, only: si => zmat_as_vec ! si(100)
use mod_version, only : version
use mod_hibrid3, only: expand
use mod_hiba07_13p, only: tcasea
use mod_par, only: batch, csflag, photof
use mod_wave, only: irec, inflev
use funit
use mod_selb, only: ibasty
use mod_ered, only: ered, rmu
use mod_hiutil, only: gennam, mtime, gettim, dater
use mod_hiutil, only: daxpy_wrapper
use mod_hivector, only: dset, matmov, dsum
implicit double precision (a-h,o-z)
character*(*) filnam
character*40  psifil, wavfil, flxfil
character*20  cdate
character*10  elaps, cpu
character*5   s13p(12)
logical exstfl, adiab, &
                kill,propf, sumf, &
                coordf
! common for y1, y2, y4
dimension a(7)  ! arguments
data s13p /'3SG0f','3SG1f','3PI0f','3PI1f','3PI2f','1PI1f', &
           '3SG1e','3PI0e','3PI1e','3PI2e','1SG0e','1PI1e'/
!

integer, pointer :: ipol
integer, parameter :: psifil_unit = 2
ipol=>ispar(3)


one=1.d0
onemin=-1.d0
zero=0.d0
! initialize timer
call mtime(cpu0,ela0)
! input
iflux=a(1)
if (iflux .gt. 4 .or. iflux .lt. -3) then
  write (6, 2) iflux
2   format (' *** IFLUX =',i3,' MUST BE -3 ... 4  ***')
  return
endif
iprint=a(2)
inchj = a(5)
inchl = a(6)
inchi = a(7)
thresh=a(3)
factr=a(4)
coordf=.false.
sumf=.false.
adiab=.false.
jflux=iabs(iflux)
if (iflux .eq. -1 .or. iflux .eq. 2) adiab = .true.
if (iflux .eq. -2) then
  sumf=.true.
  jflux=1
  adiab = .false.
endif
if (iflux .eq. 3) then
  jflux=1
  adiab=.false.
  coordf=.true.
  ny=a(2)
  ymin=a(3)
  dy=a(4)
  iprint=a(5)
  thresh=a(6)
  factr=a(7)
endif
if (jflux .eq. 4) then
  rout=a(2)
  adiab=.true.
  sumf=.false.
  coordf=.true.
endif
if (photof) then
   if (thresh .eq. 0.d0) thresh=-1.d9
   if (iprint .eq. 0) iprint= 1
else
   if (factr .eq. 0.d0) factr=1.d0
endif
if (iprint .ne. 0) then
  kill=.false.
else
  kill=.true.
endif
!
! generate filename and check if it is present
ien=0
wavfil=filnam//'.wfu'
call gennam(wavfil,filnam,ien,'wfu',lenfs)
inquire(file = wavfil, exist = exstfl)
if (.not. exstfl) then
    write(6,10) wavfil(1:lenfs)
10     format(' ** WAVEFUNCTION INFORMATION FILE ',(a), &
           ' NOT FOUND **')
    return
end if
! open file which holds transformation data and asymptotic wavefunction
call openf(22, wavfil, 'TU', 0)
call dater(cdate)
! open file to save generated wavefunction
if (jflux .eq. 0) then
  call gennam(psifil,filnam,ien,'psi',lenft)
  call openf(psifil_unit, psifil,'sf',0)
! write a header
  call version(2)
  write(2,11)
  write(6,11)
11   format(/' ** WAVEFUNCTION DETERMINATION ***',/)
  write (psifil_unit, 12) wavfil
12   format('    INFORMATION FROM FILE: ',(a))
  write (psifil_unit,13) cdate
13 format('    THIS CALCULATION ON: ',(a))
endif
if (jflux .ne. 0) then
! open file to save generated flux
  call gennam(flxfil,filnam,ien,'flx',lenft)
  call openf(3,flxfil,'sf',0)
! write a header
  call version(3)
  if (jflux .eq. 1) then
    if (photof) then
      write(3,14)
      write(6,14)
14       format(/, &
   ' ** DETERMINATION OF OUTGOING FLUX ***')
    else
      write(3,15)
      write(6,15)
15       format(/, &
   ' ** DETERMINATION OF INCOMING AND OUTGOING FLUX ***')
    endif
  endif
  if (jflux .eq. 2) then
    write(3,16)
    write(6,16)
16     format(/' ** ADIABATIC ENERGIES ***',/)
  endif
  if (jflux .eq. 4) then
    write(3,17)
    write(6,17)
17     format(/' ** TRANSFORMATION MATRIX **',/)
  endif
  write (3,18) cdate
18   format('    THIS CALCULATION ON: ',(a))
endif
! read header information, s matrix, and asymptotic wavefunction and
! derivative
call waverd(jtot,jlpar,nu,nch,npts,nopen,nphoto, &
            jflux,rstart,rendld,rinf)
if (inflev .ne. 0) then
   write (6, *) '** CALCULATION WITH WRSMAT=.T. REQUIRED.'
   goto 700
end if
if (photof) then
  write (6, 19)
  if (jflux .eq. 0) write(psifil_unit, 19)
  if (jflux .ne. 0) write (3, 19)
19   format('    PHOTODISSOCIATION BOUNDARY CONDITIONS')
else
  write (6, 20)
  if (jflux .eq. 0) write(psifil_unit, 20)
  if (jflux .ne. 0) write (3, 20)
20   format('    SCATTERING BOUNDARY CONDITIONS')
  photof=.false.
endif
if (adiab) then
  if (jflux .eq. 0)  write(psifil_unit,21)
  if (jflux .ne. 0)  write (3,21)
  write (6,21)
21   format ('    ADIABATIC BASIS')
endif
if (.not.adiab .and. .not. sumf) then
  if (.not. coordf) then
    if (ibasty .ne. 7) then
      if (jflux .eq. 0) write(psifil_unit,22)
      if (jflux .ne. 0)  write (3,22)
      write (6,22)
22       format ('    DIABATIC (ASYMPTOTIC) BASIS')
    else
      if (jflux .eq. 0) write(psifil_unit,23)
      if (jflux .ne. 0)  write (3,23)
!  print flux even inside of closed region in molecular basis
      kill = .false.
      write (6,23)
23       format ('    MOLECULAR (CASE A) BASIS')
    endif
  else
    if (ny .gt. 0) then
      if (jflux .eq. 0) write(psifil_unit,24)
      if (jflux .ne. 0) write (3,24)
      write (6,24)
24       format &
       ('    COORDINATE SPACE FLUX; POSITIVE INDEX CHOSEN')
    else
      if (jflux .eq. 0) write(psifil_unit,25)
      if (jflux .ne. 0) write (3,25)
      write (6,25)
25       format &
       ('    COORDINATE SPACE FLUX; NEGATIVE INDEX CHOSEN')
    endif
  endif
endif
if (sumf) then
      if (jflux .eq. 0) write(psifil_unit,26)
      if (jflux .ne. 0) write (3,26)
      write (6,26)
26       format &
   ('    DIABATIC (ASYMPTOTIC) BASIS;', &
    ' INELASTIC FLUX SUMMED OVER ROTATIONAL LEVELS')
endif
if (jflux .ne. 0) write (3, 12) wavfil
write (6, 12) wavfil
if (csflag) then
  if (jflux .ne. 0) &
    write(3,27) ered*econv, rmu*xmconv, jtot, nu
  if (jflux .eq.0) &
    write(2,27) ered*econv, rmu*xmconv, jtot, nu
  write(6,27) ered*econv, rmu*xmconv, jtot, nu
27   format('    ENERGY = ',f10.3,' cm(-1);  MASS = ',f9.4, &
     ' amu',/,'    CS CALCULATION:  JTOT = ',i3,'; NU =',i3)
else
  if (jflux .ne. 0) &
     write(3,29) ered*econv, rmu*xmconv, jtot, jlpar
  if (jflux .eq. 0) &
      write(2,29) ered*econv, rmu*xmconv, jtot, jlpar
  write(6,29) ered*econv, rmu*xmconv, jtot, jlpar
29   format('    ENERGY = ',f10.3,' cm(-1);  MASS = ',f9.4, &
     /,'    CC CALCULATION:  JTOT = ',i3,'; JLPAR =',i3)
endif
if (iabs(jflux).eq.1) then
  if (rendld .ge. rinf) then
    write (6, 30) rendld, rinf
30     format (' *** FLUX DESIRED; BUT RENDLD=',f7.3, &
            ' .GE. RINF=',f7.3)
    return
  endif
  write(3, 31) rendld, rinf
  write (6, 31) rendld, rinf
31   format ('    FLUXES DETERMINED FROM R = ', &
        f7.3,' TO R = ',f7.3)
  if (coordf) then
    write (6,34) ymin, dy, ymin+(iabs(ny)-1)*dy
    write (3,34) ymin, dy, ymin+(iabs(ny)-1)*dy
34     format ('                           R-INT = ',f5.2,':', &
          f5.2,':',f5.2)
  endif
  write(3,32) thresh
  write (6,32) thresh
32   format ('    CLOSED CHANNEL THRESHOLD = ',1pg10.3)
  write(3,33) factr
  write (6,33) factr
33   format ('    FACTOR FOR CLOSED CHANNEL DAMP = ',f5.2)
else if(jflux.eq.2) then
  write(3, 35) rendld, rinf
  write (6, 35) rendld, rinf
35   format ('    ADIABATIC ENERGIES DETERMINED FROM R = ', &
        f7.3,' TO R = ',f7.3)
else if(jflux.eq.4) then
  write (6,36) rout
  write (3,36) rout
36   format ( &
 '    DIABATIC->ADIABATIC TRANSFORMATION REQUESTED AT R = ', &
     f7.3)
else if(jflux.eq.0) then
  write(2, 37) rstart, rinf, npts
  write (6, 37) rstart, rinf, npts
37   format ( &
  /,'    WF DEFINED FROM R = ',f7.3,' TO R = ',f7.3,' AT ', &
   i4,' POINTS',/)
endif
inch=0
! check if initial channel is in list of channels
! not for photodissociation
if (.not. photof) then
  if (ibasty .ne. 7) then
    do 40 nn=1, nch
      j1 = jq(nn)
      l1 = lq(nn)
      i1 = inq(nn)



 write(6,443) nn,j1,l1,il,inchj,inchl,inchi
443  format(' ch#',i3,3i5,'  req:',3i5)



      if(j1.eq.inchj.and.l1.eq.inchl.and.i1.eq.inchi) then
        inch=nn
        goto 41
      endif
40     continue
  else
    inch=inchj
  endif
41   if ((jflux .ne. 2 .and. jflux .ne. 4).and. inch .eq. 0) then
    if (jflux.ne.0) write(3, 43) inchj, inchl, inchi
    write (6, 43) inchj, inchl, inchi
43     format( /,' ** CHANNEL (J, l, IN = ',2i3,i4,') NOT IN LIST')
    return
  endif
endif
do 45 i = 1, nch
45   nalist(i)=i
! reorder channels in increasing energy since this is eispack ordering
if (nch .gt. 1) then
  call dcopy(nch,eint,1,sc1,1)
  do 50 i = 1,nch-1
    do 48 j = i+1,nch
      if (sc1(j).lt.sc1(i)) then
!  switch
        ehold=sc1(i)
        sc1(i)=sc1(j)
        sc1(j)=ehold
        inhold=nalist(i)
        nalist(i)=nalist(j)
        nalist(j)=inhold
      endif
48     continue
50   continue
endif
do 51 i = 1, nch
  if (inch .eq. nalist(i)) then
    incha=i
    goto 52
  endif
51 continue
52 continue
! nalist now contains ordering of channels in energy
!   nalist(i) is the number in original list of the channel of ith
!   lowest energy
if (.not.coordf)  then
  if (jflux.ne.0) write(3, 55)
  write (6, 55)
55   format('    CHANNEL PARAMETERS:', &
       '   N   J   L   IN    ENERGY(CM-1)    SQRT(K)')
  if (jflux .ne. 2) then
    inchc=inch
    if (adiab) inchc=incha
    if (.not.photof) then
      if (ibasty .ne. 7 .or. &
          (ibasty .eq. 7 .and. (ipol .eq. 0 .or. adiab &
           .and. jlpar .eq. 1))) then
        if (jflux.eq.0) &
          write(2, 57) inchc, jq(inch), lq(inch), inq(inch), &
            econv*eint(inch)
        if (jflux.ne.0) &
          write(3, 57) inchc, jq(inch), lq(inch), inq(inch), &
            econv*eint(inch)
        write (6, 57) inchc, jq(inch), lq(inch), inq(inch), &
            econv*eint(inch)
57         format(/,15x,'INITIAL:',3i4,i5,f13.3)
      else if (ibasty .eq.7 .and. jlpar .eq. -1 &
              .and. ipol.ne.0) then
        if (jflux.eq.0) &
          write(2, 58) inchc, s13p(inch+6), &
            econv*eint(inch)
        if (jflux.ne.0) &
          write(3, 58) inchc, s13p(inch+6), &
            econv*eint(inch)
        write (6, 58) inchc, s13p(inch+6), &
            econv*eint(inch)
58         format(/,15x,'INITIAL:  ',i2,3x,a5,f18.3)
      endif
    else
! initial channel always 1 for photodissociation, since only one
! column in wavefunction
      inch=1
      incha=1
      inchc=1
    endif
  endif
endif
if (jflux .eq. 4) then
  write (3,55)
  do 60  i=1, nch
    write (3, 59) i, jq(i), lq(i), inq(i), econv*eint(i)
59     format(10x,4i4,f13.3)
60   continue
endif
if (photof) then
  inch=1
  incha=1
  inchc=1
endif
nchsq=nch*nch
nopsq=nopen*nopen
! make a list of pointers
! nlist is pointer to desired probed channels in full channel list
nj = 0
if (ibasty .ne. 7) then
  if (nnout .lt. 0) then
    write (6, 65) nnout
65     format ('  ** WARNING: NNOUT = ',i3, &
          ' .LE. O IN SUBROUTINE PSI')
  else if (nnout .eq. 0) then
    if(jflux.ne. 4) then
      write (6, 68)
68       format ('  ** NO PROBE STATES REQUESTED?  ABORT ***')
      go to 700
    endif
  endif
endif
do 120 i=1, iabs(nnout)
   jo = jout(i)
   if (jflux .eq. 2) then
     if (ibasty .ne. 4) then
       llo = iabs(indout(i))/100
       io = sign(iabs(indout(i))-100*llo,indout(i))
     else
       llo = iabs(indout(i))/1000
       io = sign(iabs(indout(i))-1000*llo,indout(i))
     endif
   endif
   do 100 nn=1, nch
     j1 = jq(nn)
     if(j1.ne.jo) goto 100
       if (jflux .ne. 2) then
         do 75 in = 1, iabs(niout)
           io = indout(in)
           if(inq(nn).ne.io) goto 75
           nj = nj + 1
           nlist(nj) =nn
75          continue
       else
          innq=inq(nn)
          if(innq.ne.io .or. lq(nn) .ne. llo) goto 100
! check to see if state has already been found
           ifound=0
           do 80 if=1,nj
             if (nn .eq. nlist(if)) ifound=1
80            continue
           if (ifound .eq. 0) then
             nj = nj + 1
             nlist(nj) = nn
           endif
       endif
100    continue
120 continue
! check if there had been any match
! if coordinate space flux or 13p scattering, include all states
if (coordf .or. ibasty .eq. 7 .or. sumf) then
  nj=nch
  do 121 i=1, nch
    nlist(i)=i
121   continue
endif
if(nj.eq.0) then
  if(jflux.ne.4) then
    if (jflux .ne. 0) write(3,130)
    if (jflux .eq. 0) write(2,130)
    write(6,130)
    if(.not. batch) write(6,130)
130       format(' *** NO PROBE STATES FOUND, ABORT ***')
    goto 700
  endif
end if
! reorder list in terms of energy
if (adiab) then
  if (nj .gt. 1) then
    do 135 i=1, nj-1
      do 134 j=i+1,nj
        if (eint(nlist(j)) .lt. eint(nlist(i))) then
          nhold=nlist(i)
          nlist(i)=nlist(j)
          nlist(j)=nhold
        endif
134       continue
135     continue
  endif
endif
! establish equivalence with adiabatic level list
do 138 i = 1, nj
  do 136 j =1,nch
    if(nalist(j) .eq. nlist(i)) then
      nblist(i)=j
      goto 138
    endif
136   continue
138 continue
if (.not. coordf) then
  if (.not.sumf) then
    do 142 i=1, nj
      nn=nlist(i)
      nnn=nn
      if (adiab) nnn=nblist(i)
      if (eint(nn) .le. ered) then
        sq=sqrt(2.d0*rmu*(ered-eint(nn)))
      else
        sq=-sqrt(2.d0*rmu*(-ered+eint(nn)))
      endif
      if (ibasty .ne. 7 .or. adiab) then
        if (jflux.eq.0) &
        write(2, 140) nnn, jq(nn), lq(nn), inq(nn), &
                 eint(nn)*econv, sq
        if (jflux.ne.0) write(3,140) nnn,jq(nn),lq(nn),inq(nn), &
                 eint(nn)*econv, sq
        write (6, 140) nnn, jq(nn), lq(nn), inq(nn), &
                 eint(nn)*econv, sq
140         format(16x,'PROBED:',3i4,i5,f13.3,f13.4)
      else
        if (jlpar .eq. -1) nn=nn+6
        if (jflux.eq.0) &
        write(2, 141) nnn, s13p(nn), &
                 eint(nnn)*econv, sq
        if (jflux.ne.0) write(3,141) nnn,s13p(nn), &
                 eint(nnn)*econv, sq
        write (6, 141) nnn, s13p(nn), &
                 eint(nnn)*econv, sq
141         format(16x,'PROBED:  ',i2,3x,a5,f18.3,f13.4)
      endif
142     continue
  else
    do 145 ni = 1, niout
      write (3, 143) indout(ni)
      write (6, 143) indout(ni)
143       format(16x,'PROBED:  INDEX =',i4)
145     continue
  endif
endif
if (jflux.eq.1) then
  if(.not.coordf) then
     write(3, 146)
     write (6, 146)
146      format(17x,'TOTAL:  LAST COLUMN')
  else
     write(3, 147)
     write (6, 147)
147      format('    TOTAL FLUX IN LAST COLUMN')
  endif
endif
! reverse order of adiabatic states, since eispack routines return
! highest energy first
do 148 i=1,nj
  nalist(i)=nch-nblist(i)+1
148 continue
npoint=(inch-1)*nch + 1
! npoint points to top of column inch of full wavefunction matrix
if (ibasty .eq.7 .and. jlpar.eq.-1 .and. ipol.eq.1) then
!         if (adiab) then
!           write (6, 149)
!           if (jflux.ne.0) write(3, 149)
!149        format('    STATE 5 IS PARALLEL POLARIZATION, STATE 6'm
!     :            ' IS PERPENDICULAR')
!         endif
   call waverot(jtot,nch)
 endif
! if 13p  or 2s-2p scattering, determine the matrix for case (a) -> case (e)
if (ibasty .eq.7) call tcasea(jtot,jlpar)
if (jflux .eq. 0) then
! here if wavefunction calculation
  if (.not.photof) then
!        if (photof) then
! here for scattering
! now expand psir, psii, dpsir, dpsii into nch x nch matrix, putting zeros
! as closed channel components
    if (nch .gt. nopen) then
      call expand(nopen,nopen,nch,nch,ipack, &
                psir,psii,scmat)
      call expand(nopen,nopen,nch,nch,ipack, &
                dpsir,dpsii,scmat)
    endif
! store desired wavefunction column (real and imaginary part) in 1st and 2nd
! columns of psir
    if (inch .ne. 1) call dcopy(nch,psir(npoint),1,psir,1)
    call dcopy(nch,psii(npoint),1,psir(nch+1),1)
    write(2, 150)
150     format(/' R (BOHR) AND REAL PART OF WAVEFUNCTION', &
          ' (R < 0 INDICATES AIRY PROPAGATION)',/)
    call psicalc(npts, nch, nchsq, nj, psifil_unit)
! copy imaginary part of asymptotic wavefunction into first column of psir
! so we can use same loop as above
    call dcopy(nch,psir(nch+1),1,psir,1)
    write(2, 185)
185     format(/' R (BOHR) AND IMAGINARY PART OF WAVEFUNCTION', &
          '(R < 0 INDICATES AIRY PROPAGATION)',/)
    call psicalc(npts,nch,nchsq,nj, psifil_unit)
  else
! here for photdissociation, in which case outgoing wavefunction is a
! given column of chi
! npoint points to which of the nphot ground state wavefunctions are to
! be selected
    if (nch .gt. nopen) then
      call expand(nopen,nopen,nch,nch,ipack, &
                  psir,psii,scmat)
      call expand(nopen,nopen,nch,nch,ipack, &
                  dpsir,dpsii,scmat)
    endif
! store desired wavefunction column (real and imaginary part) in 1st and
! 2nd columns of psir
    call dcopy(nch,psir(npoint),1,psir,1)
    call dcopy(nch,psii(npoint),1,psir(nch+1),1)
! store derivatives (real and imaginary)  in 3rd and 4th columns of psir
    call dcopy(nch,dpsir(npoint),1,psir(2*nch+1),1)
    call dcopy(nch,dpsii(npoint),1,psir(3*nch+1),1)
    ipoint=2*nch+1
    irec=npts+4
    write(2, 200)
200     format(/' R (BOHR) AND REAL PART OF CHI')
    iwf = 1
    propf=.true.
    call flux(npts,nch,nchsq,ipoint,nj,adiab,thresh,factr,kill, &
            photof,propf,sumf,iwf,coordf,ny,ymin,dy,psifil_unit)
    write(2, 210)
210     format(/' R (BOHR) AND IMAGINARY PART OF CHI')
! reread asymptotic information
    call waverd(jtot,jlpar,nu,nch,npts,nopen,nphoto, &
            jflux,rstart,rendld,rinf)
    if (nch .gt. nopen) then
      call expand(nopen,nopen,nch,nch,ipack, &
                  psir,psii,scmat)
      call expand(nopen,nopen,nch,nch,ipack, &
                  dpsir,dpsii,scmat)
    endif
! store desired wavefunction column (real and imaginary part) in 1st and
! 2nd columns of psir
    call dcopy(nch,psir(npoint),1,psir,1)
    call dcopy(nch,psii(npoint),1,psir(nch+1),1)
! store derivatives (real and imaginary)  in 3rd and 4th columns of psir
    call dcopy(nch,dpsir(npoint),1,psir(2*nch+1),1)
    call dcopy(nch,dpsii(npoint),1,psir(3*nch+1),1)
    iwf = -1
    irec=npts+4
    call flux(npts,nch,nchsq,ipoint,nj,adiab,thresh,factr,kill, &
            photof,propf,sumf,iwf,coordf,ny,ymin,dy,psifil_unit)
  endif
else if (jflux .eq. 2) then
  write(3, 300)
300   format(/' R (BOHR) AND ADIABATIC ENERGIES (CM-1)',/)
  irec=npts+4
  call eadiab(npts,nch,nj)
else if (jflux .eq. 4) then
  call transmt(npts,nch,nchsq,rout)
else if (jflux .eq. 1) then
! here for flux calculation
! first for total flux (only for scattering)
! now expand psir, psii, dpsir, dpsii into nch x nch matrix, putting zeroz
! as closed channel components
  if (.not. photof) then
    if (nch .gt. nopen) then
      call expand(nopen,nopen,nch,nch,ipack, &
           psir,psii,scmat)
      call expand(nopen,nopen,nch,nch,ipack, &
           dpsir,dpsii,scmat)
    endif
! store desired wavefunction column (real and imaginary part) in 1st and 2nd
! columns of psir
    call dcopy(nch,psir(npoint),1,psir,1)
    call dcopy(nch,psii(npoint),1,psir(nch+1),1)
! store derivatives (real and imaginary)  in 3rd and 4th columns of psir
    call dcopy(nch,dpsir(npoint),1,psir(2*nch+1),1)
    call dcopy(nch,dpsii(npoint),1,psir(3*nch+1),1)
! now compute desired total fluxes
    write(3, 305)
305     format(/' R (BOHR) AND TOTAL FLUXES (UNITS OF HBAR/MU; ', &
            'NO DAMPING)'/)
! first initial flux (not if adiabatic or not if 13p scattering or
! not summed fluxes)
    scsum=0.d0
    if (.not. adiab .and. ibasty .ne. 7 .and. .not.sumf) then
      do 310 i=1, nj
        nni=nlist(i)
        sc(i)=psir(nni)*psir(3*nch+nni)-psir(nch+nni) &
              *psir(2*nch+nni)
      scsum=scsum+sc(i)
310       continue
      write(3, 320) rinf, (sc(i), i=1,nj), scsum
320       format(f10.4,30(1pe11.3))
    endif

    if ((.not. adiab .and. .not. coordf) .and. ibasty .ne.7) then
    scsum=0.d0
    if (.not.sumf) then
      do 321 i=1, nj
        nni=nlist(i)
        sc(i)=psir(nni)*psir(3*nch+nni)-psir(nch+nni) &
            *psir(2*nch+nni)
        scsum=scsum+sc(i)
321       continue
      nout=nj
    else
! here for sum of fluxes over desired indices
      do 322 i=1,niout
        sc(i)=0.d0
322       continue
      scsum=0.d0
      do 330 i=1,nch
        nni=nalist(i)
        scc=psir(nni)*psir(3*nch+nni)-psir(nch+nni) &
            *psir(2*nch+nni)
        do 327 ni=1, niout
          if (inq(nni) .eq. indout(ni)) sc(ni)=sc(ni)+scc
327         continue
330       continue
      scsum=dsum(niout,sc,1)
      nout=niout
    endif
    write(3, 320) rinf, (sc(i), i=1,nout), scsum
    endif
    irec=npts+4
    ipoint=2*nch+1
    iwf = 0
    propf=.true.
! plot out all fluxes for total flux which is numerically well behaved
    tthresh=-1.e9
    call flux(npts,nch,nchsq,ipoint,nj,adiab,thresh,factr,.false., &
            photof,propf,sumf,iwf,coordf,ny,ymin,dy,psifil_unit)
  endif
  if (.not. photof) then
! now for incoming flux (only for scattering)
! first construct real and imaginary parts of incoming waves
! here for scattering, in which case incoming flux is a diagonal matrix
    call dset(nchsq,zero,psir,1)
    call dset(nchsq,zero,psii,1)
    call dset(nchsq,zero,dpsir,1)
    call dset(nchsq,zero,dpsii,1)
    ipoint=1
    do 335 i=1, nch
      psir(ipoint)=-fn(i)
      psii(ipoint)=-fj(i)
      dpsir(ipoint)=-fnp(i)
      dpsii(ipoint)=-fjp(i)
      ipoint=ipoint+nch+1
335     continue
    write(3, 340)
340     format(/' R (BOHR) AND INCOMING FLUXES (UNITS OF HBAR/MU)',/)
    if (ibasty .eq.7 .and. jlpar.eq.-1 .and. ipol.eq.1) &
        call waverot(jtot,nch)
! store desired wavefunction column (real and imaginary part) in 1st and 2nd
! columns of psir
    call dcopy(nch,psir(npoint),1,psir,1)
    call dcopy(nch,psii(npoint),1,psir(nch+1),1)
! store derivatives (real and imaginary)  in 3rd and 4th columns of psir
    call dcopy(nch,dpsir(npoint),1,psir(2*nch+1),1)
    call dcopy(nch,dpsii(npoint),1,psir(3*nch+1),1)
! now compute desired incoming fluxes
! first initial flux (not if adiabatic or not if 13p scattering or
! not summed fluxes)
    if (.not. adiab .and. ibasty .ne. 7 ) then
      if (.not.sumf) then
        scsum=0.d0
        do 360 i=1, nj
          nni=nlist(i)
          sc(i)=psir(nni)*psir(3*nch+nni)-psir(nch+nni) &
              *psir(2*nch+nni)
        scsum=scsum+sc(i)
360         continue
        write(3, 320) rinf, (sc(i), i=1,nj), scsum
      else
! here for sum of fluxes over desired indices
        do 362 i=1,niout
          sc(i)=0.d0
362         continue
        scsum=0.d0
        do 365 i=1,nch
          nni=nalist(i)
          scc=psir(nni)*psir(3*nch+nni)-psir(nch+nni) &
            *psir(2*nch+nni)
          do 364 ni=1, niout
            if (inq(nni) .eq. indout(ni)) sc(ni)=sc(ni)+scc
364          continue
365         continue
        scsum=dsum(niout,sc,1)
        nout=niout
      endif
      write(3, 320) rinf, (sc(i), i=1,nout), scsum
    endif
    irec=npts+4
    ipoint=2*nch+1
    iwf = 0
    propf=.false.
    call flux(npts,nch,nchsq,ipoint,nj,adiab,thresh,factr,kill, &
            photof,propf,sumf,iwf,coordf,ny,ymin,dy,psifil_unit)
  endif
! now for outgoing flux
  if (.not.photof) then
! outgoing wave: move Sreal into psir and Simag into psii
    call matmov (sr, psir, nopen, nopen, nopen, nopen)
    call matmov (si, psii, nopen, nopen, nopen, nopen)
! premultiply Sr by diagonal matrix yl(kr) and Si by jl(kr)
    do 430 irow = 1, nopen
      fac1=fn(irow)
      fac2=fj(irow)
      call dscal(nopen, fac1, psir(irow), nopen)
      call dscal(nopen, fac2, psii(irow), nopen)
430     continue
! add together, resave in psir, this is real part of outgoing wave
    call daxpy_wrapper(nopsq, one, psii, 1, psir, 1)
! repeat for derivative of part of outgoing wave
    call matmov (sr, dpsir, nopen, nopen, nopen, nopen)
    call matmov (si, dpsii, nopen, nopen, nopen, nopen)
! premultiply Sr by diagonal matrix -ylp(kr) and Si by jlp(kr)
    do 440 irow = 1, nopen
      fac1=fnp(irow)
      fac2=fjp(irow)
      call dscal(nopen, fac1, dpsir(irow), nopen)
      call dscal(nopen, fac2, dpsii(irow), nopen)
440     continue
! add together, resave in dpsir, this is derivative of
! real part of outgoing wave
    call daxpy_wrapper(nopsq, one, dpsii, 1, dpsir, 1)
! repeat for imaginary part of outgoing wave
    call matmov (si, psii, nopen, nopen, nopen, nopen)
    call matmov (sr, scmat, nopen, nopen, nopen, nopen)
! premultiply Sr by diagonal matrix -jl(kr) and Si by yl(kr)
    do 450 irow = 1, nopen
      fac1=fn(irow)
      fac2=-fj(irow)
      call dscal(nopen, fac1, psii(irow), nopen)
      call dscal(nopen, fac2, scmat(irow), nopen)
450     continue
! add together, resave in psii, this is imaginary part of outgoing wave
    call daxpy_wrapper(nopsq, one,scmat, 1, psii, 1)
! repeat for derivative of imaginary part of outgoing wave
    call matmov (si, dpsii, nopen, nopen, nopen, nopen)
    call matmov (sr, scmat, nopen, nopen, nopen, nopen)
! premultiply Sr by diagonal matrix jl(kr) and Si by yl(kr)
    do 460 irow = 1, nopen
      fac1=fnp(irow)
      fac2=-fjp(irow)
      call dscal(nopen, fac1, dpsii(irow), nopen)
      call dscal(nopen, fac2, scmat(irow), nopen)
460     continue
! add together, resave in dpsii, this is imaginary part of derivative
! of outgoing wave
    call daxpy_wrapper(nopsq, one,scmat, 1, dpsii, 1)
! now expand psir, psii, dpsir, dpsii into nch x nch matrix, putting zeroz
! as closed channel components
    if (nch .gt. nopen) then
      call expand(nopen,nopen,nch,nch,ipack, &
                  psir,psii,scmat)
      call expand(nopen,nopen,nch,nch,ipack, &
                  dpsir,dpsii,scmat)
    endif
    if (ibasty .eq.7 .and. jlpar.eq.-1 .and. ipol.eq.1 &
      .and. .not.photof) &
       call waverot(jtot,nch)
  else
! here for photdissociation, in which case outgoing wavefunction is a
! given column of chi
! npoint points to which of the nphot ground state wavefunctions are to
! be selected
    if (nch .gt. nopen) then
      call expand(nopen,nopen,nch,nch,ipack, &
                  psir,psii,scmat)
      call expand(nopen,nopen,nch,nch,ipack, &
                  dpsir,dpsii,scmat)
    endif
  endif
! store desired wavefunction column (real and imaginary part) in 1st and 2nd
! columns of psir
  call dcopy(nch,psir(npoint),1,psir,1)
  call dcopy(nch,psii(npoint),1,psir(nch+1),1)
! store derivatives (real and imaginary)  in 3rd and 4th columns of psir
  call dcopy(nch,dpsir(npoint),1,psir(2*nch+1),1)
  call dcopy(nch,dpsii(npoint),1,psir(3*nch+1),1)


! now compute desired outgoing fluxes
  if (.not. photof) write(3, 550)
550   format(/' R (BOHR) AND OUTGOING FLUXES (UNITS OF HBAR/MU)',/)
  if (photof) write(3, 551)
551   format(/' R (BOHR) AND OUTGOING FLUXES (ATOMIC UNITS)',/)
! first initial flux
  if (coordf) then
    write (3, 554) (ymin+dy*(iy-1), iy=1,iabs(ny))
554     format('    R-INT:',f9.3,30(f11.3))
  endif
  if ((.not. adiab .and. .not. coordf) .and. ibasty .ne.7) then
    scsum=0.d0
    if (.not.sumf) then
      do 560 i=1, nj
       nni=nlist(i)
        sc(i)=psir(nni)*psir(3*nch+nni)-psir(nch+nni) &
            *psir(2*nch+nni)
        scsum=scsum+sc(i)
560       continue
      nout=nj
    else
! here for sum of fluxes over desired indices
        do 575 i=1,niout
          sc(i)=0.d0
575         continue
        scsum=0.d0
        do 580 i=1,nch
          nni=nalist(i)
          scc=psir(nni)*psir(3*nch+nni)-psir(nch+nni) &
            *psir(2*nch+nni)
          do 579 ni=1, niout
            if (inq(nni) .eq. indout(ni)) sc(ni)=sc(ni)+scc
579          continue
580         continue
        scsum=dsum(niout,sc,1)
        nout=niout
    endif
    if (photof) then
      call dscal(nj,1.d0/rmu,sc,1)
      scsum=scsum/rmu
    endif
    write(3, 320) rinf, (sc(i), i=1,nout), scsum
  endif
  ipoint=2*nch+1
  irec=npts+4
  iwf = 0
  if (photof) propf=.true.
  if (.not. photof) propf=.false.
  call flux(npts,nch,nchsq,ipoint,nj,adiab,thresh,factr,kill, &
            photof,propf,sumf,iwf,coordf,ny,ymin,dy,psifil_unit)
endif
700 if (photof .or. jflux .eq. 0) close (psifil_unit)
if (jflux .ne. 0) close (3)
close (FUNIT_WFU)
call mtime(cpu1,ela1)
cpu1 = cpu1 - cpu0
ela1 = ela1 - ela0
call gettim(cpu1,cpu)
call gettim(ela1,elaps)
if(.not. batch) then
  if (iflux .eq. 0) write(6,720) elaps, cpu
720   format(/,' ** WAVEFUNCTION CALCULATION FINISHED:', &
       /,'    ELAPSED TIME:',(a),'  CPU TIME: ',(a))
  if (iflux .ne. 0) write(6,730) elaps, cpu
730   format(/,' ** FLUX CALCULATION FINISHED:', &
       /,'    ELAPSED TIME:',(a),'  CPU TIME: ',(a))
endif
return
end
! ------------------------------------------------------------------
subroutine flux(npts,nch,nchsq,ipoint,nj,adiab,thresh,factr,kill, &
                photof, propf, sumf,iwf,coordf,nny,ymin,dy,psifil_unit)
!
! subroutine to calculate fluxes
!
! author: millard alexander
! current revision date (algorithm): 19-apr-1996 by mha
! revised on 30-mar-2012 by q. ma for stream I/O of wfu files
! current revision: 20-apr-2012 by q. ma
!
! ------------------------------------------------------------------
use mod_coiout, only: niout, indout
use mod_cocent, only: sc2 => cent
use mod_coamat, only: psir ! psir(100) (4,nch)
use mod_cobmat, only: psii ! psii(100) Here psii is used as a vector
use mod_cotq2, only: scmat2 => dpsii ! scmat2(100)
use mod_cotq3, only: scmat3 => scmat ! scmat3(100)
use mod_coinq, only: inq ! inq(60)
use mod_cosc1, only: pk => sc1 ! pk(6)
use mod_coisc2, only: nlist => isc2 ! nlist(60)
use mod_coisc3, only: nalist => isc3 ! nalist(60)
use mod_cosc6, only: sc => sc6 ! sc(60)
use mod_cosc7, only: sc1 => sc7 ! sc1(6)
use mod_cosc8, only: sc8
use mod_cosc9, only: sc9
use mod_coz, only: scmat => z_as_vec ! scmat(100)
use mod_cozmat, only: tcoord => zmat_as_vec ! tcoord(100)
use mod_wave, only: irec, ifil, nrlogd
use mod_coqvec, only: nphoto
use mod_selb, only: ibasty
use mod_ered, only: ered, rmu
use mod_hiba07_13p, only: ttrans
use mod_hiutil, only: daxpy_wrapper
use mod_himatrix, only: mxma
use mod_hivector, only: dset, vadd, vmul, dsum
! steve, you may need more space, but i doubt it since tcoord is dimensioned n
implicit double precision (a-h,o-z)
logical adiab, kill, photof, propf, sumf, coordf

dimension scc(100)
data zero, one, onemin /0.d0, 1.d0, -1.d0/
data ione, mone /1,-1/
integer :: psifil_unit
! if propf = true then true back-subsititution for flux
! if propf = false then inward propagation
! noffset is start of 5th column of psir
  noffset=4*nch+1
!        open (unit=23, file='propagators.txt',status='unknown')
! determine coordinate matrix if coordf true
  if (coordf) then
    ind=1
    ny=iabs(nny)
    do 50 i= 1, nch
      sc(i)=zero
      if(nny .gt. 0 .and. inq(i).gt. 0) sc(i)=one
      if(nny .lt. 0 .and. inq(i).lt. 0) sc(i)=one
50     continue
! sc is now mask for those states for which index is desired
    do 100 iy = 1, ny
      y=ymin+(iy-1)*dy
      call wfintern(scmat, y, nch, nphoto, 1, .false.)
! steve, you'll need to modify wfintern so that scmat returns both function an
! scmat is a vector of length nch containing the nch internal states
! evaluated at internal coordinate y
      call vmul(sc,1,scmat,1,tcoord(ind:),1,nch)
! this masks the internal states depending on whether the index is
! positive or negative
      ind=ind+nch
100     continue
  endif
! tcoord now contains as rows internal states and as columns
! values of internal coordinate

! here beings loop over sectors (R), starting from outermost sector
  do 420 kstep=1, npts
    irec=irec-1
    read (ifil, end=900, err=950, pos=iwavsk(irec)) r
! branch out if we've gone beyond airy propagation region
    if (r .gt. 0) goto 450
    read (ifil, end=900, err=950) drnow
! read in local wavevectors
    read (ifil, end=900, err=950) (sc8(i), i=1, nch)
! read in first G(n,n+1) matrix,  then Tn transformation
! matrix into local interval
    read (ifil, end=900, err=950) (scmat3(i), i=1, nchsq), &
         (scmat(i), i=1, nchsq)
! read in propagators (y1=pk, y2=sc1, y4=sc2, gam1=sc9,
!     muab= 5th column of psi)
    read (ifil, end=900, err=950) (pk(i), i=1, nch), &
         (sc1(i), i=1, nch), (sc2(i), i=1, nch), &
         (sc9(i), i=1, nch), (psir(noffset - 1 + i), i=1, nch)
!          write (23, 299) -r, drnow, (scmat(ii), ii=1,4),
!     :      (pk(ii),ii=1,2),(sc1(ii),ii=1,2),
!     :      (sc2(ii),ii=1,2),(sc9(ii),ii=1,2),(psir(noffset+ii),ii=0,1)

! 299     format (2f16.12,14(1pe22.12e3))
! transform wave function into local basis
    if (propf) then
! scmat2(nch, 2) = scmat(nch, nch) * psir(nch, 2)
      call mxma(scmat,1,nch,psir,1,nch,scmat2,1,nch,nch,nch,2)
      call dcopy(2*nch,scmat2,1,psii(ipoint),1)
! here for back substitution for wavefunction
! 3rd and 4th columns of psii contain real and imaginary parts of function
! propagate functions
! evaluate G-tilde* psi-tilde(b)
       call mxma(scmat3,1,nch,psii(ipoint),1,nch,psii,1, &
               nch,nch,nch,2)
! 1st two column of psii now contain G-tilde(A-B)*psi-tilde(b)
! if photodissociation, subtract off mu-tilde(a,b) from real part
      if (photof) call vadd(mone,psii,1 &
                          ,psir(noffset),1,nch)
! at this point 1st two columns of psii contain real and imaginary
!   psi-tilde(a)
! 3rd and 4th columns still contain psi-tilde(b)
! propagate derivatives
      call vmul(pk,1,psii,1,scmat2,1,nch)
      call vmul(pk,1,psii(1+nch:),1,scmat2(1+nch:),1,nch)
! first and second columns of scmat2 now contain y1 Fa
      call vmul(sc1,1,psii(ipoint:),1,scmat2(ipoint:),1,nch)
      call vmul(sc1,1,psii(ipoint+nch:),1, &
                scmat2(ipoint+nch:),1,nch)
! 3rd and 4th columns of scmat2 now contain y2 Fb
      call daxpy_wrapper(2*nch,onemin,scmat2,1,scmat2(ipoint),1)
! 3rd and 4th columns of scmat2 now contains (-y1 Fa + y2 Fb)
! this is psip-tilde(a)  move to 2nd and 3rd columns of psir
      call dcopy(2*nch,scmat2(ipoint),1,psir(ipoint),1)
! if photodissociation subtract off gamma1 from real part of derivative
      if (photof) call vadd(mone,psir(ipoint),1,sc9,1,nch)
! for compatibility with previous version, copy derivative into 3rd and
! 4th columns of psii
      call dcopy(2*nch,psir(ipoint),1,psii(ipoint),1)
    else
      do 300 ii=1, nch
        sc1(ii)=onemin/sc1(ii)
300       continue
! scmat2(nch, 4) = scmat(nch, nch) * psir(nch, 4)
      call mxma(scmat,1,nch,psir,1,nch,scmat2,1,nch,nch,nch,4)
      call dcopy(2*nch,scmat2,1,psii(ipoint),1)
      call dcopy(2*nch,scmat2(ipoint),1,psii,1)
! 1st two columns of psii now contain real and imaginary part of derivative
! 3rd and 4th columns contain real and imaginary parts of function
! propagate functions
      call vmul(sc2,1,psii(ipoint:),1,scmat2,1,nch)
      call vmul(sc2,1,psii(ipoint+nch:),1,scmat2(1+nch:),1,nch)
! first and second columns of scmat2 now contain y4 Fb
      call daxpy_wrapper(2*nch,onemin,scmat2,1,scmat2(ipoint),1)
! 3rd and 4th columns of scmat2 now contains (Fb' - y4 Fb)
! if photodissociation, subtract off gamma2 from real part
      if (photof) call vadd(mone,scmat2(2*nch+1),1 &
                          ,psir(noffset),1,nch)
! if photodissociation, 3rd column of scmat2 now contains
!     Re(Fb' - y4 Fb - gamma 2)
! multiply by - y2^-1 to get Fa
      call vmul(sc1,1,scmat2(ipoint:),1,psii,1,nch)
      call vmul(sc1,1,scmat2(ipoint+nch:),1,psii(1+nch:),1,nch)
! 1st and 2nd columns of psii now contain Fa
      do 320 ii=1, nch
        sc1(ii)=onemin/sc1(ii)
320       continue
      call vmul(sc1,1,psii(ipoint:),1,scmat2,1,nch)
      call vmul(sc1,1,psii(ipoint+nch:),1,scmat2(1+nch:),1,nch)
! 1st and 2nd columns of scmat2 now contain y2 Fb
! if photodissociation subtract off gamma1 from real part
      if (photof) call vadd(mone,scmat2,1,sc9,1,nch)
      call dscal(nch,onemin,pk,1)
      call vmul(pk,1,psii,1,psii(ipoint:),1,nch)
      call vmul(pk,1,psii(1+nch:),1,psii(ipoint+nch:),1,nch)
! 3rd and 4th columns of psii now contain -y1 Fa
! add on y2 Fb and store in 3rd and 4th columns of psii
      call daxpy_wrapper(2*nch,one,scmat2,1,psii(ipoint),1)
! if wavevector matrix positive (lambda negative) channel is closed
! kill the corresponding components of psi and psi'
! thresh is the threshold for killing closed channel components
    do 330 ii=1, nch
      if (sc8(ii) .lt. thresh) then
        fact=exp(-sqrt(abs(sc8(ii)))*drnow*factr)
        psii(ii)=fact*psii(ii)
        psii(ii+nch)=fact*psii(ii+nch)
        psii(ii+2*nch)=fact*psii(ii+2*nch)
        psii(ii+3*nch)=fact*psii(ii+3*nch)
      endif
330     continue
    endif
! 1st two columns of psii now contain Fa
! 3rd and 4th columns of psii now contain Fa'
    if (adiab) then
! here for flux calculation in locally adiabatic basis
      scsum=0.d0
      do 360 i=1, nj
        nni=nalist(i)
        sc(i)=psii(nni)*psii(3*nch+nni)- &
              psii(nch+nni)*psii(2*nch+nni)
!  store real or imaginary parts of wf if desired
        if (iwf .eq. 1) pk(i)=psii(nni)
        if (iwf .eq. -1) pk(i)=psii(nch+nni)
! if wavevector matrix positive (lambda negative) channel is closed
! kill the corresponding components of psi and psi'
!              if (sc8(nni) .lt. thresh.and.kill) then
        if (sc8(nni) .lt. 0.d0.and.kill) then
          sc(i)=zero
          if (iwf .ne. 0) pk(i)=0.d0
        endif
        scsum=scsum+sc(i)
360       continue
      nout=nj
    endif
! transform wavefunction and derivative into asymptotic basis
    call mxma(scmat,nch,1,psii,1,nch,psir,1,nch,nch,nch,4)
! psir now contains this information
!          write (23, 299) -r, drnow, (psir(ii), ii=1,8)
    if (.not.adiab) then
! here for flux calculation in asymptotic basis
! calculate flux
      if (coordf) then
! here for coordinate space calculation of fluxes
        scsum=zero
! fluxes will be stored in vector sc, initialize to zero
        call dset(nch,zero,sc,1)
        do i=1,ny
          ind=(i-1)*nch+1
          scc1=ddot(nch,psir,1,tcoord(ind),1)
! scc1 contains sum(psi-real*phi_internal)
          scc2=ddot(nch,psir(nch+1),1,tcoord(ind),1)
! scc2 contains sum(psi-imag*phi_internal)
          scc3=ddot(nch,psir(2*nch+1),1,tcoord(ind),1)
! scc3 contains sum(dpsi-real*phi_internal)
          scc4=ddot(nch,psir(3*nch+1),1,tcoord(ind),1)
! scc4 contains sum(dpsi-imag*phi_internal)
          sc(i)=scc1*scc4-scc2*scc3
!steve, you'll need to append to this to calculate r-component of current dens
! try using sc9 as scratch storage for these
        enddo
        call dscal(ny,dy,sc,1)
! multiply by step width to recover J.R at ri times dri
        do i=1, ny
          scsum=scsum+sc(i)
        enddo
        nout=ny
      endif
      if (.not. coordf) then
! here for channel fluxes in diabatic basis
! transform wavefunction and derivative into molecular basis (only for
! 13p or 2s-2p
        if (sumf) call dset(niout,zero,scc,1)
        if (ibasty .eq. 7) then
! psii(nch, 4) = ttrans(nch, nch) * psir(nch, 4)
          call mxma(ttrans,nch,1,psir,1,nch, &
                    psii,1,nch,nch,nch,4)
! N. B.  as written, mxma multiplies psir by the transpose of ttrans
! psii now contains this transformed wavefunction and derivative
! copy these back to psir
          call dcopy(4*nch,psii,1,psir,1)
        endif
        do 370 i=1, nj
          mmi=nalist(i)
          nni=nlist(i)
          sc(i)=psir(nni)*psir(3*nch+nni)- &
            psir(nch+nni)*psir(2*nch+nni)
!  store real or imaginary parts of wf if desired
          if (iwf .eq. 1) pk(i)=psii(nni)
          if (iwf .eq. -1) pk(i)=psii(nch+nni)
! if wavevector matrix positive (lambda negative) channel is closed
! kill the corresponding components of psi and psi'
          if (sc8(mmi) .lt. 0.d0.and.kill) then
!              if (sc8(mmi) .lt. thresh.and.kill) then
            sc(i)=zero
            if (iwf .ne. 0) pk(i)=0.d0
          endif
370         continue
        if (sumf) then
          do 390 i=1, nj
            mmi=nalist(i)
            do 375 ni=1, niout
              if (inq(mmi) .eq. indout(ni)) then
                scc(ni)=scc(ni)+sc(mmi)
              endif
375             continue
390           continue
        endif
        if (.not.sumf) then
          nout=nj
          scsum=dsum(nj,sc,1)
        else
          nout=niout
          scsum=dsum(niout,scc,1)
        endif
! transform wavefunction and derivative back into asymptotic basis (only for
! 13p or 2s-2p
        if (ibasty .eq. 7) then
          call mxma(ttrans,1,nch,psir,1,nch, &
                    psii,1,nch,nch,nch,4)
! psii now contains this transformed wavefunction and derivative
! copy these back to psir
          call dcopy(4*nch,psii,1,psir,1)
        endif
      endif
    endif
    if (iwf .ne. 0) &
      write(psifil_unit, 400) -r, (pk(i), i=1,nj)
    if (iwf .eq. 0) then
      if (photof) then
! for photodissociation, so, steve, you'll need to scale sc9 also
        call dscal(nout,1.d0/rmu,sc,1)
        scsum=scsum/rmu
        if (scsum .lt. 1.d-13) scsum=zero
      endif
      if (sumf) then
        write (3, 400) -r, (scc(i), i=1,nout), scsum
      else
        write (3, 400) -r, (sc(i), i=1,nout), scsum
! steve you'll also have to write out sc9
      endif
400       format(f10.4,30(1pe11.3))
    endif
420   continue

450   continue
!        close (23)
  return
!
900   continue
950   write (0, *) '*** ERROR READING WFU FILE (FLUX). ABORT'
  call exit()
  end
! ------------------------------------------------------------------
subroutine eadiab(npts,nch,nj)
!
! subroutine to readin and print out adiabatic energies
!
! author: millard alexander
! current revision date (algorithm): 12-may-1997 by mha
! revised on 30-mar-2012 by q. ma for stream I/O of wfu files
!
! ------------------------------------------------------------------
use mod_coisc3, only: nalist => isc3 ! nalist(10)
use mod_cosc6, only: sc => sc6 ! sc(6)
use mod_cosc8, only: sc8
use mod_wave, only: irec, ifil
use mod_ered, only: ered, rmu
implicit none
integer, intent(in) :: npts
integer, intent(in) :: nch
integer, intent(in) :: nj
integer :: i, nni
integer :: kstep
real(8) :: r
real(8) :: drnow
! common for y1, y2, y4
real(8), parameter :: two = 2.d0
real(8), parameter :: conv = 219474.6d0
  do 420 kstep=1, npts
    irec=irec-1
    read (ifil, end=900, err=950, pos=iwavsk(irec)) r
! branch out if we've gone beyond airy propagation region
    if (r .gt. 0) return
    read (ifil, end=900, err=950) drnow
! read in local wavevectors
    read (ifil, end=900, err=950) (sc8(i), i=1, nch)
    do 370 i=1, nj
      nni=nalist(i)
      sc(i)=-conv*(sc8(nni)/(two*rmu)-ered)
370     continue
!          write (3, 170) -r+0.5*drnow, (nalist(i), i=1,nj),
    write (3, 170) -r+0.5*drnow, &
        (sc(i), i=1,nj)
!170       format(f10.4,2i4,20(1pe12.4))
170     format(f10.4,20(1pe12.4))
420   continue
  return
!
900   continue
950   write (0, *) '*** ERROR READING WFU FILE (EADIAB). ABORT'
  call exit()
  end
! ------------------------------------------------------------------
subroutine waverot(jtot,nch)
! special for singlet-triplet mixing;
! rearrange wavefunction to correspond to state assignment:
! psi-parallel = (Jtot+1)^1/2 |el=J+1> - Jtot^1/2 |el=J-1>
! psi-perpendicular = Jtot^1/2 |el=J+1> + (Jtot+1)^1/2 |el=J-1>
! with assignment that original column five is singlet with el=J-1 and
! column 6 is singlet with el=J+1
! use orthogonal plane rotation
use mod_coamat, only: psir ! psir(100) (4,nch)
use mod_cobmat, only: psii ! psii(100) Here psii is used as a vector
use mod_cotq1, only: dpsir ! dpsir(100) 
use mod_cotq2, only: dpsii ! dpsii(100)
implicit none
integer, intent(in) :: jtot
integer, intent(in) :: nch
integer :: jpoint
real(8) :: cs, sn ! cosine and sine of rotation angle
real(8) :: jtot_as_real
real(8) :: norm
jtot_as_real=jtot
cs=sqrt(jtot_as_real+1.d0)
sn=sqrt(jtot_as_real)
norm=sqrt(2.d0*jtot_as_real+1.d0)
cs=cs/norm
sn=sn/norm
jpoint=4*nch+1
call drot(nch,psir(jpoint),1,psir(jpoint+nch),1,cs,sn)
call drot(nch,psii(jpoint),1,psii(jpoint+nch),1,cs,sn)
call drot(nch,dpsir(jpoint),1,dpsir(jpoint+nch),1,cs,sn)
call drot(nch,dpsii(jpoint),1,dpsii(jpoint+nch),1,cs,sn)
return
end
! ------------------------------------------------------------------
subroutine psicalc(npts, nch, nchsq, nj, psifil_unit)
!
! subroutine to propagate wavefunctions inward
!
! author: millard alexander
! current revision date (algorithm): 4-oct-1991 by mha
! revised on 30-mar-2012 by q. ma for stream I/O of wfu files
!
! ------------------------------------------------------------------
use mod_coamat, only: psir ! psir(100)  psir(nch, 1)
use mod_coisc2, only: nlist => isc2 ! nlist(6)
use mod_cosc6, only: sc => sc6 ! sc(6)
use mod_cosc7, only: sc1 => sc7 ! sc1(6)
use mod_cow, only: sr => w_as_vec ! sr(100)
use mod_cozmat, only: si => zmat_as_vec ! si(100)
use mod_wave, only: irec, ifil
use mod_himatrix, only: mxma
implicit none
integer, intent(in) :: npts
integer, intent(in) :: nch
integer, intent(in) :: nchsq
integer, intent(in) :: nj
integer, intent(in) :: psifil_unit
integer :: i, kstep
real(8) :: drnow, r
  irec=npts+4
  do 180 kstep=1, npts
    irec=irec-1
    read (ifil, end=900, err=950, pos=iwavsk(irec)) r, drnow
! read in adiabtic energies (which we don't need)
    read (ifil, end=900, err=950) (sc1(i), i=1, nch)
! read in transformation matrices
    read (ifil, end=900, err=950) (sr(i), i=1, nchsq)
! si(nch, 1) = sr(nch, nch) * psir(nch, 1)
    call mxma(sr,1,nch,psir,1,nch,si,1,nch,nch,nch,1)
    do 160 i=1, nj
      sc(i)=si(nlist(i))
160     continue
    call dcopy(nch,si,1,psir,1)
    write(psifil_unit, 170) r, (sc(i), i=1,nj)
170     format(f10.4,12(1pe10.2))
180   continue
return
!
900 continue
950 write (0, *) '*** ERROR READING WFU FILE (PSICALC). ABORT.'
call exit()
end
! ------------------------------------------------------------------
subroutine transmt(npts,nch,nchsq,rout)
!
! subroutine to print out transformation matrix at rout
!
! author: millard alexander
! current revision date (algorithm): 18-may-2008 by mha
! revised on 30-mar-2012 by q. ma for stream I/O of wfu files
! current revision: 20-apr-2012 by q. ma
!
! ------------------------------------------------------------------
use mod_cosc7, only: sc1 => sc7 ! sc1(6)
use mod_cow, only: sr => w_as_vec ! sr(100)
use mod_cozmat, only: si => zmat_as_vec ! si(100)
use mod_wave, only: irec, ifil
use mod_himatrix, only: transp
implicit double precision (a-h,o-z)
logical renormf
dimension scrvec(64)
irec=npts+4
delold=1.d+18
do 200 kstep=1, npts
    irec=irec-1
    read (ifil, end=900, err=950, pos=iwavsk(irec)) r, drnow
    del=abs(abs(r)-rout)
! read in adiabtic energies (which we don't need)
    read (ifil, end=900, err=950) (sc1(i), i=1, nch)
! read in first G(n,n+1) matrix,  then Tn transformation
! matrix into local interval
    read (ifil, end=900, err=950) (sr(i), i=1, nchsq), &
         (si(i), i=1, nchsq)
    if (rout .gt. 0.d0) then
      if (del .gt.delold) then
! transpose matrix is stored column by column
! after transposition, columns correspond to eigenvectors
        call transp (si, nch, nch)
        write (3,160) -r+0.5*drnow
160    format('    TRANSFORMATION MATRIX AT R = ',f10.6,' IS:',/)
        write (6,161) abs(r)
161    format('    TRANSFORMATION MATRIX DETERMINED AT R = ',f10.6)
! print out transpose matrix (reverse order since eispack
! determines highest eigenvalue first
        ind=nchsq-nch+1
! if more than 64 channels, then eliminate all elements with i>64, and renormalize vector
        nnch=nch
        if (nch.gt.64) then
           renormf=.true.
           write (3,162) nch
           write (6,162) nch
162    format('    CHANNEL EXPANSION TRUNCATED FROM NCH =',i4, &
             ' TO NCH = 64;')
           write (3,163)
           write (6,163)
163    format('       EIGENVECTOR RENORMALIZED')
           nnch=64
        else
           renormf=.false.
        endif
        call openf(4,'tmatrix.dat','sf',0)
        nstate=min(12,nch)
        do 190 ii=1,nstate
!             do 190 ii=1,nch
          if (renormf) then
             call dcopy(64,si(ind),1,scrvec,1)
             cnorm=ddot(64,scrvec,1,scrvec,1)
             cnorm=1d0/sqrt(cnorm)
             call dscal(64,cnorm,scrvec,1)
             call dcopy(64,scrvec,1,si(ind),1)
          endif
          write (3,175) -r+0.5*drnow, ii, &
            (si(ij), ij=ind,ind+nnch-1)
          write (4,175) -r+0.5*drnow, ii, &
            (si(ij), ij=ind,ind+nnch-1)
175         format(f7.4,i3,64f10.6)
          ind=ind-nch
190         continue
        close(4)
        return
      endif
    else
        call transp (si, nch, nch)
! print out transpose matrix (reverse order since eispack
! determines highest eigenvalue first
        ind=nchsq-nch+1

        call openf(4,'tmatrix.dat','sf',0)
        do 195 ii=1,nstate
          if (renormf) then
             call dcopy(64,si(ind),1,scrvec,1)
             cnorm=ddot(64,scrvec,1,scrvec,1)
             cnorm=1d0/sqrt(cnorm)
             call dscal(64,cnorm,scrvec,1)
             call dcopy(64,scrvec,1,si(ind),1)
          endif
!              do 195 ii=1,nch
          write (4,175) -r+0.5*drnow, ii, &
            (si(ij), ij=ind,ind+nnch-1)
          ind=ind-nch
195         continue
        close(4)
    endif
    call dcopy(nchsq,sr,1,si,1)
    delold=del
200 continue
! here if rout not reached at end of airy propagation region, print out
! last transformation matrix
call transp (si, nch, nch)
write (3,160)-r+0.5*drnow
write (6,160)-r+0.5*drnow
! print out transpose matrix
ind=nchsq-nch+1
do 250 ii=1,12
!      do 250 ii=1,nch
  write (3,175) -r+0.5*drnow,  ii, (si(ij), ij=ind,ind+nch-1)
  ind=ind-nch
250 continue
return
!
900 continue
950 write (0, *) '*** ERROR READING WFU FILE (TRANSMT). ABORT'
call exit()
end
!
!     ------------------------------------------------------------------
subroutine eadiab1(filnam, nchmin, nchmax)
!
!     Subroutine to readin and print out adiabatic energies from a wfu
!     file.  wfu files with only adiabatic energy information
!     (wrsmat=.f. when performaing wavefunction calculations) can be
!     supported.
!
!     Adiabatic energies as a function of r will be saved to
!     filnam.eadiab, with columns the adiabatic bender curves.  r is
!     stored in the first column.
!
!     author: qianli ma
!     current revision: 20-apr-2012
!
!     ------------------------------------------------------------------
use constants
use mod_cosc8, only: sc8
use mod_wave, only: ifil, nrlogd
use funit
use mod_ered, only: ered, rmu
use mod_hiutil, only: gennam
implicit none
character*(*), intent(in) :: filnam
integer, intent(in) :: nchmin
integer, intent(inout) :: nchmax

integer :: i, j
integer :: jtot, jlpar, nu, nch, npts, nopen, nphoto
integer(8) :: noffst
integer :: lenfs, lenft, jflux
integer :: nchpr
real(8) :: drnow, rstart, rendld, rinf, r
logical :: exstfl
character*40 :: wavfil, eadfil
!
integer(8) :: seek_pos
!
double precision :: dble_t
integer, parameter :: eadfil_unit = FUNIT_EADIAB
integer, parameter :: ien = 0
!
if (nchmax .ne. 0 .and. nchmax .lt. nchmin) goto 990
wavfil = filnam // '.wfu'
call gennam(wavfil, filnam, ien, 'wfu', lenfs)
inquire (file=wavfil, exist=exstfl)
if (.not. exstfl) then
   write (6, 10) wavfil(1:lenfs)
10    format (' ** WAVEFUNCTION INFORMATION FILE ', (a), &
        ' NOT FOUND **')
   return
end if
call openf(FUNIT_WFU, wavfil, 'TU', 0)

eadfil = filnam // '.eadiab'
call gennam(eadfil, filnam, ien, 'eadiab', lenft)
call openf(eadfil_unit, eadfil, 'sf', 0)
write (6, 15) eadfil(1:lenft)
15 format (' ** WRITING ADIABATIC ENERGIES TO ', (a))
!
call waverd(jtot, jlpar, nu, nch, npts, nopen, nphoto, &
     jflux, rstart, rendld, rinf)
if (nchmin .gt. nch) goto 990
if (nchmax .eq. 0 .or. nchmax .gt. nch) nchmax = nch
nchpr = nchmax - nchmin + 1
noffst = (nch - nchmax + 2) * sizeof(dble_t)
!
write(eadfil_unit, 17) nchmin, nchmax
17 format (' ** ADIABATIC ENERGIES FROM NO.', i5, ' TO NO.', &
     i5, ' REQUESTED')
do i = 4 + nrlogd, npts + 3
   seek_pos = iwavsk(i)
   read (ifil, end=900, err=950, pos=seek_pos) r, drnow
   read (ifil, end=900, err=950, pos=seek_pos+noffst) &
        (sc8(j), j=1, nchpr)
   write(eadfil_unit, 20) -r + 0.5 * drnow
20    format (f10.5, 1x, $)
   do j = nchpr, 1, -1
      write(eadfil_unit, 30) -econv * (sc8(j) / (2d0 * rmu) - ered)
   end do
30    format (f13.5, 1x, $)
   write(eadfil_unit, 20)
end do
write (eadfil_unit, 20)
goto 990
!
900 continue
950 write (0, *) '*** ERROR READING WFU FILE'
990 close(ifil)
close(eadfil_unit)
return
end
end module mod_hibrid4