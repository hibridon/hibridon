#include "assert.h"
module mod_hismat
contains

! ---------------------------------------------------------------------------
! hismat library (hibridon s-matrix library)
!
! subroutines included:
!
!   sread (sinqr/readrc/rdhead/saddr) stream s-matrix i/o
!   swrite (wrhead) stream s-matrix i/o
!
! ---------------------------------------------------------------------------
!
!     ------------------------------------------------------------
subroutine readrc(iadr, smt_file_unit, lrec, jtot, jlpar, nu, nopen, &
     length, nnout)
!     author: h.j. werner
!     revision: 14-jun-1990 by mha
!     rewritten: 07-jan-2012 by q. ma
!
!     Read the header of the s-matrix for a single partial wave
!
!     INPUT ARGUMENTS:
!     iadr: 0 if read next record, otherwize read from byte iadr (one
!     based) of the file.
!
!     RETURNED ARGUMENTS:
!     lrec: on return contains the size in byte current s-matrix,
!     or -1 on end of file
!     Other parameters have their traditional meanings.
!     ------------------------------------------------------------
implicit none
integer :: iadr
integer, intent(in) :: smt_file_unit
integer, intent(out) :: lrec
integer, intent(out) :: jtot
integer, intent(out) :: jlpar
integer, intent(out) :: nu
integer, intent(out) :: nopen
integer, intent(out) :: length
integer, intent(out) :: nnout
if (iadr .eq. 0) inquire (smt_file_unit, pos=iadr)
read (smt_file_unit, pos=iadr, end=900, err=950) &
     lrec, jtot, jlpar, nu, nopen, length ,nnout
return
!
900 lrec = -1
return
950 write (0, *) '*** ERROR READING S-MATRIX FILE (readrc). ABORT.'
call exit()
end
!     ------------------------------------------------------------
!
!     ------------------------------------------------------------
subroutine rdhead(smt_file_unit, cdate, ered, rmu, csflag, flaghf, &
     flagsu, twomol, nucros, jfirst, jfinal, jtotd, numin, numax, nud, &
     nlevel, nlevop, nnout, jlev, inlev, elev, jout)
!
!     subroutine to read header from file smt_file_unit
!     authors: h.j. werner and b. follmeg
!     revision: 27-oct-95
!     major revision: 07-jan-2012 by q. ma
!     ------------------------------------------------------------
use mod_parpot, only: potnam=>pot_name, label=>pot_label
implicit none
integer, intent(in) :: smt_file_unit ! logical unit used to read smt file 
character*20, intent(out) :: cdate
double precision, intent(out) :: ered
double precision, intent(out) :: rmu
logical, intent(out) :: csflag, flaghf, flagsu, twomol, nucros
integer, intent(out) :: jfirst, jfinal, jtotd
integer, intent(out) :: numin, numax, nud
integer, intent(out) :: nlevel, nlevop, nnout
integer, dimension(1), intent(out) :: jlev  ! dimension(1:nlevel)
integer, dimension(1), intent(out) :: inlev ! dimension(1:nlevel)
double precision, dimension(1), intent(out) :: elev ! dimension(1:nlevel)
integer, dimension(1), intent(out) :: jout  ! dimension(1:iabs(nnout))
character*8 csize8
character*4 csize4
integer :: lenhd
integer :: ibasty
integer :: i

!     Read the magic number (from the start of the file)
read (smt_file_unit, pos=1, end=900, err=950) csize8
!     first word on file contains the length of the header
read (smt_file_unit, end=900, err=950) lenhd
!
read (smt_file_unit, end=900, err=950) cdate, label, potnam
!     Four bytes for alignment
read (smt_file_unit, end=900, err=950) csize4
!
read (smt_file_unit, end=900, err=950) ered, rmu, csflag, flaghf, flagsu, &
     twomol, nucros, jfirst, jfinal, jtotd, numin, numax, &
     nud, nlevel, nlevop, nnout, ibasty
read (smt_file_unit, end=900, err=950) (jlev(i), i=1, nlevel), &
     (inlev(i), i=1, nlevel), (elev(i), i=1, nlevel)
read (smt_file_unit, end=900, err=950) (jout(i), i=1, iabs(nnout))
!
read (smt_file_unit, end=900, err=950) csize8
if (csize8 .ne. 'ENDOFHDR') goto 950
return
!
900 continue
950 write (0, *) '*** ERROR READING S-MATRIX FILE. ABORT.'
call exit
end
!     ------------------------------------------------------------------
subroutine sinqr(smt_file_unit, njtot, nchmx)
!
!     Subroutine to scan an s-matrix file for the number of partial
!     waves and the maximum number of channels.
!
!     If any error occured, njtot is set to -1.
!
!     The file pointer will be pointed to the end of the file on return.
!     It is recommended to call this subroutine prior to calling rdhead.
!
!     Author: Qianli Ma
!     ------------------------------------------------------------------
implicit none
integer smt_file_unit, njtot, nchmx, iaddr, lrec
integer jtot, jlpar, nu, nopen, length, nnout
!
!     Read the length of the header
read (smt_file_unit, pos=9, end=900, err=950) iaddr
!
njtot = 0
nchmx = 0
iaddr = iaddr + 1
20 call readrc(iaddr, smt_file_unit, lrec, jtot, jlpar, nu, nopen, length, &
     nnout)
if (lrec .lt. 0) return
njtot = njtot + 1
if (nchmx .lt. length) nchmx = length
iaddr = iaddr + lrec
goto 20
900 continue
950 njtot = -1
return
end
subroutine sread (iadr, sreal, simag, jtot, jlpar, nu, &
                  row_bqs, packed_bqs, &
                  smt_file_unit, nmax, nopen, ierr)
!     authors: h.j. werner and b. follmeg
!     revision: 21-feb-2006 by mha
!     major revision: 07-feb-2012 by q. ma
!     current revision:  8-oct-2012 by q. ma
!
!     read real and imaginary parts of s-matrix together with
!     other information as written in soutpt, swrite
!     if iadr = 0 read sequential (next record)
!     if iadr > 0 read absolute
!     if nopen = -1, the lower triangle is filled
!     ------------------------------------------------------------
use mod_hibasis, only: is_j12
use mod_selb, only: ibasty
use mod_hitypes, only: bqs_type
implicit none
integer, intent(in) :: iadr
real(8), intent(out) :: sreal(nmax,1)
real(8), intent(out) :: simag(nmax,1)
integer, intent(out) :: jtot
integer, intent(out) :: jlpar
integer, intent(out) :: nu
type(bqs_type), intent(out) :: row_bqs
type(bqs_type), intent(out) :: packed_bqs
integer, intent(in) :: smt_file_unit
integer, intent(in) :: nmax
integer, intent(inout) :: nopen
integer, intent(out) :: ierr
logical :: triang
character*8 :: csize8
integer :: iaddr
integer :: lrec, i, nnout, irow, icol, ioff

call row_bqs%init(nmax)
call packed_bqs%init(nmax)

!
ierr=0
triang =.false.
if (nopen < 0) then
   triang = .true.
   nopen = iabs(nopen)
end if
!     on return:
!
!     the packed_bqs will hold the column indices of
!     the packed basis (dimension length)
!
!     the vectors row_bqs%jq, row_bqs%lq, row_bqs%inq will hold the row indices of the packed
!     basis (dimension nopen)
!
!     read next s-matrix header
iaddr = iadr
call readrc(iaddr,smt_file_unit,lrec,jtot,jlpar,nu,nopen,packed_bqs%length,nnout)
if (lrec < 0) goto 900

read (smt_file_unit, end=900, err=950) &
     (packed_bqs%jq(i), i=1, packed_bqs%length), &
     (packed_bqs%lq(i), i=1, packed_bqs%length), &
     (packed_bqs%inq(i), i=1, packed_bqs%length)
if (is_j12(ibasty)) then
     read (smt_file_unit, end=900, err=950) (packed_bqs%j12(i), i=1, packed_bqs%length)
end if
!
if (nnout .gt. 0) then
   do 50 i = 1, packed_bqs%length
      row_bqs%jq(i) = packed_bqs%jq(i)
      row_bqs%lq(i) = packed_bqs%lq(i)
      row_bqs%inq(i)= packed_bqs%inq(i)
      if (is_j12(ibasty)) then
         row_bqs%j12(i) = packed_bqs%j12(i)
      end if
50    continue
   nopen = packed_bqs%length
   if (triang) then
      ioff = 1
      do 70 irow = 1, packed_bqs%length
         read (smt_file_unit, end=900, err=950) &
              (sreal(ioff + i, 1), i=0, irow - 1), &
              (simag(ioff + i, 1), i=0, irow - 1)
         ioff = ioff + irow
70       continue
      goto 800
   end if
!     read s-matrix
   do 80  icol = 1, packed_bqs%length
      read (smt_file_unit, end=900, err=950) &
           (sreal(i, icol), i=1, icol)
      read (smt_file_unit, end=900, err=950) &
           (simag(i, icol), i=1, icol)
80    continue
!     fill lower triangle
   do icol=1,packed_bqs%length
      do irow=1,icol
         sreal(icol,irow)=sreal(irow,icol)
         simag(icol,irow)=simag(irow,icol)
      end do
   end do
!
else if (nnout .le. 0) then
!     here if you have written out columns of the s-matrix
   read (smt_file_unit, end=900, err=950) (row_bqs%jq(i), i=1, nopen), &
        (row_bqs%lq(i), i=1, nopen), (row_bqs%inq(i), i=1, nopen)
   if (is_j12(ibasty)) then
      read (smt_file_unit, end=900, err=950) (row_bqs%j12(i), i=1, nopen)
   end if
!     now read columns of the s-matrix
   do 140 icol = 1, packed_bqs%length
      read (smt_file_unit, end=900, err=950) &
           (sreal(i, icol), i=1, nopen), &
           (simag(i, icol), i=1, nopen)
140    continue
end if
row_bqs%length = nopen
!
!     Read eight bytes "ENDOFSMT"
800 read (smt_file_unit, end=900, err=950) csize8
if (csize8 .ne. 'ENDOFSMT') goto 950
return
!
!     End-of-file
900 continue
!     Read error
950 ierr = -1
return
end
!     ------------------------------------------------------------
!
!     ------------------------------------------------------------
subroutine swrite (sreal, simag, jtot, jlpar, nu, &
                   row_bqs, iorder, packed_bqs, &
                   epack, nfile, nmax, nopen)
!  subroutine to write selected elements of s-matrix to file nfile
!  author:  millard alexander
!  modified by  h.j. werner and b. follmeg
!  revision: 21-feb-2006 by mha
!  major revision: 07-jan-2012 by q. ma
!  added /coj12/ common block (p.dagdigian)
!  current revision:  8-oct-2012 by q. ma
!  ------------------------------------------------------------------
!  variables in call list:
!    sreal:     on entry: contains real part of open-channel s-matrix
!               on return: contains real part of packed s-matrix
!    simag:     on entry: contains imaginary part of open-channel s-matrix
!               on return: contains imaginary part of packed s-matrix
!    jtot:      total angular momentum
!    csflag:    if .true. coupled-states calculation
!               if .false. close-coupled calculation
!    flaghf:    if .true., then system with half-integer spin
!                if .false., then system with integer spin
!    nu:        cs projection index (not used in cc calculation)
!    row_bqs%jq:        channel rotational angular momenta
!    row_bqs%lq:        channel orbital angular momenta
!    row_bqs%inq:       additional quantum index of each channel
!    note!!!   if flaghf = .true., then the true values
!    of the rotational quantum numbers, the total angular momentum,
!    and the coupled-states projection index are equal to the values
!    stored in row_bqs%jq, jtot, and nu plus 1/2
!    nfile:     logical unit for output of s-matrices
!    nmax:      maximum row dimension of matrices
!    nopen:     number of channels

!  ------------------------------------------------------------------
use mod_cosout, only: nnout, jout
use mod_coeint, only: eint
use mod_hibasis, only: is_j12
use mod_selb, only: ibasty
use mod_hitypes, only: bqs_type
implicit none
real(8), intent(inout) :: sreal(nmax,nmax)
real(8), intent(inout) :: simag(nmax,nmax)
integer, intent(in) :: jtot
integer, intent(in) :: jlpar
integer, intent(in) :: nu
type(bqs_type), intent(in) :: row_bqs
integer, intent(out) :: iorder(nopen)
type(bqs_type), intent(out) :: packed_bqs
real(8), intent(out) :: epack(nopen)
integer, intent(in) :: nfile
integer, intent(in) :: nmax
integer, intent(in) :: nopen

integer :: ic, icol, i, ii, ir, irow, length, mmout
integer :: nchnid, nrecw

integer int_t
double precision dble_t
!

! by default, each channel is uses three parameters: row_bqs%j, row_bqs%l, row_bqs%inq
nchnid = 3
if (is_j12(ibasty)) then
   ! some basis have an additional channel parameter row_bqs%j12
   nchnid = nchnid + 1
end if
!
! the vector iorder will point to the position in the unpacked basis
! of each state in the packed basis
!
! the vector epack will hold the channel energies in the packed basis
! first sum over the unpacked states
!
mmout = iabs(nnout)
call packed_bqs%init(nopen)
length = 0
do icol = 1, nopen
   ! now sum over the packed states, find labels
   do ii = 1, mmout
      if (row_bqs%jq(icol) == jout(ii) ) then
         ! here if match
         length = length + 1
         packed_bqs%jq(length) = row_bqs%jq(icol)
         packed_bqs%lq(length) = row_bqs%lq(icol)
         epack(length) = eint(icol)
         packed_bqs%inq(length) = row_bqs%inq(icol)
         if (is_j12(ibasty)) packed_bqs%j12(length) = row_bqs%j12(icol)
         iorder(length) = icol
         exit
      end if
   end do
end do
packed_bqs%length = length
! calculate number of words that will be written
nrecw = sizeof(int_t) * 7 + sizeof(int_t) * nchnid * length
if(nnout > 0) then
   nrecw = nrecw + sizeof(dble_t) * length * (length + 1)
else
   nrecw = nrecw + sizeof(int_t) * nchnid * nopen + &
        sizeof(dble_t) * 2 * length * nopen
end if
nrecw = nrecw + 8
! write out general information on next record
write (nfile, err=950) nrecw, jtot, jlpar, nu, nopen, &
     packed_bqs%length ,nnout
write (nfile, err=950) (packed_bqs%jq(i), i=1, length)
write (nfile, err=950) (packed_bqs%lq(i), i=1, length)
write (nfile, err=950) (packed_bqs%inq(i), i=1, length)
if (is_j12(ibasty)) then
   write (nfile, err=950) (packed_bqs%j12(i), i=1, length)
end if

! here we pack the s-matrix and print out just those elements for
! which the initial and final rotational quantum numbers correspond
! to an element in the array jout
if (nnout > 0) then
   ! the dimension of the packed s-matrix is length x length now pack
   ! the real part of the s-matrix
   do icol = 1, length
      ic = iorder(icol)
      do irow = 1, length
         ir = iorder(irow)
         sreal(irow,icol) = sreal(ir,ic)
         simag(irow,icol) = simag(ir,ic)
      end do
   end do
   ! write s-matrix into buffer
   do icol = 1, length
      write (nfile, err=950) (sreal(i, icol), i=1, icol)
      write (nfile, err=950) (simag(i, icol), i=1, icol)
   end do
   ! here if you want to print out columns of the s-matrix
else
   ASSERT (nnout <= 0)
   write (nfile, err=950) (row_bqs%jq(i), i=1, nopen)
   write (nfile, err=950) (row_bqs%lq(i), i=1, nopen)
   write (nfile, err=950) (row_bqs%inq(i), i=1, nopen)
   if (is_j12(ibasty)) write (nfile, err=950) (row_bqs%j12(i), i=1, nopen)
   ! now write out columns of the s-matrix into buffer length is the
   ! number of columns of the s-matrix to save
   do ii = 1, length
      icol = iorder(ii)
      write (nfile, err=950) (sreal(i, icol), i=1, nopen), &
           (simag(i, icol), i=1, nopen)
   end do
end if
write (nfile, err=950) 'ENDOFSMT'
return
!
! On error
950 write (0, *) '*** ERROR WRITING S-MATRIX FILE. ABORT.'
call exit()
end
! ------------------------------------------------------------
!
! ------------------------------------------------------------
subroutine wrhead(nfile,cdate, &
     ered,rmu,csflag,flaghf, &
     flagsu,twomol,nucros,jfirst,jfinal,jtotd,numin,numax,nud, &
     nlevel,nlevop,nnout,jlev,inlev,elev,jout)
!
!     initialize buffering for file nfile and write general information
!     (header)
!
!     author: h.j. werner
!     revision: 27-oct-1995 by mha
!     major revision: 07-jan-2012 by q.ma (stream I/O, write ibasty)
!     ------------------------------------------------------------
use mod_parpot, only: potnam=>pot_name, label=>pot_label
use mod_selb, only: ibasty
implicit double precision (a-h,o-z)
logical csflag, flaghf, flagsu, twomol, nucros
character*20 cdate
dimension jlev(1),inlev(1),elev(1),jout(1)
integer int_t
!     Calculate the length of the header (in bytes)
lenhd = sizeof(cdate) + sizeof(label) + sizeof(potnam) + &
     sizeof(dble_t) * 2 + sizeof(int_t) * 16 + &
     (sizeof(int_t) * 2 + sizeof(dble_t)) * nlevel + &
     sizeof(int_t) * iabs(nnout) + 20
!
!     Write eight-byte file magic number
write (nfile, err=950) char(128), 'SMT', char(0), char(2), &
     char(0), char(0)
!
write (nfile, err=950) lenhd
write (nfile, err=950) cdate, label, potnam
!     Four zero-bytes for alignment / C struct compatibility
write (nfile, err=950) char(0), char(0), char(0), char(0)
!
write (nfile, err=950) ered, rmu
write (nfile, err=950) csflag, flaghf, flagsu, twomol, nucros, &
     jfirst, jfinal, jtotd, numin, numax, nud, nlevel, nlevop, &
     nnout, ibasty
!
write (nfile, err=950) (jlev(i), i=1, nlevel)
write (nfile, err=950) (inlev(i), i=1, nlevel)
write (nfile, err=950) (elev(i), i=1, nlevel)
!
write (nfile, err=950) (jout(i), i=1, iabs(nnout))
!
write (nfile, err=950) 'ENDOFHDR'
!
return
!
!     On error:
950 write (0, *) '*** ERROR WRITING S-MATRIX FILE. ABORT.'
call exit()
return
end
end module mod_hismat
