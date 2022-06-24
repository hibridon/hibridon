! ------------------------------------------------------------------
!  in this header, arrays are set up as allocatable.
!  in addition, procedures to allocate and deallocate these
!  arrays are included here.  see below for full details about
!  hitensor
!
!  revision:  6-jun-2013 (q. ma)
! ------------------------------------------------------------------
module tensor
implicit none
integer :: jmx, kmx, lmx, kkmx, lbufs, lbuflb
parameter (jmx=90, kmx=2*jmx+1, lmx=kmx, kkmx=3*kmx)
parameter (lbufs=13504500, lbuflb=9000)
real(8), dimension(:), allocatable :: srealp, simagp
integer, dimension(:), allocatable :: ipackp, jpackp, lpackp
integer :: ialloc
!
contains
!
subroutine tensor_allocate()
allocate(srealp(lbufs), stat=ialloc)
allocate(simagp(lbufs), stat=ialloc)
allocate(ipackp(lbuflb), stat=ialloc)
allocate(jpackp(lbuflb), stat=ialloc)
allocate(lpackp(lbuflb), stat=ialloc)
if (ialloc .ne. 0) then
   print *, '  *** MEMORY ALLOCATION FAILS, EXITING ***'
   call tensor_free()
   call exit
end if
end subroutine tensor_allocate
!
subroutine tensor_free()
if (allocated(srealp)) deallocate(srealp)
if (allocated(simagp)) deallocate(simagp)
if (allocated(ipackp)) deallocate(ipackp)
if (allocated(jpackp)) deallocate(jpackp)
if (allocated(lpackp)) deallocate(lpackp)
end subroutine tensor_free
end module tensor
! ------------------------------------------------------------------
subroutine tenopa(filnam,a)
!
! subroutine to calculate tensor opacities from s-matrix elements

! see B. Follmeg, P. Rosmus, and H.-J. Werner,
! J. Chem. Phys. 93, 4687 (1990)

! this subroutine returns the sigma(lambda,ki,kf,ji,jf) of Eq. 30
! of this paper, or, for lambda=0, sigma(ji,jf,K) which is Eq. 31
!
! this subroutine requires CC s-matrix for both parities
!
! author: b. follmeg
! revision date:  13-oct-2011 by pjd
! changed open command to open smt files for streaming reads
! replaced call closf(1) with close(1)  (2 places)
! revision date:  8-aug-2012 by q. ma
!
! revised by pj dagdigian - to get SIGK working, etc.
!  with half-integral j and even multiplicity open-shell molecules (12-nov-2008)
! SIGKP not working.  Do not set N (also called lambda) GT 0
! added calculation of tensor cross sections with quantization
!  axis along initial relative velocity vector (lambda = 0 only) (5-mar-2010)
! added calculation of m-resolved differential cross sections for an
!  elastic transition in the helicity frame
!  (quant. axis for initial levels along initial velocity,
!  quant. axis for final level along final vel. vector.  this is followed
!  by evaluation of the tensor diff. cross sections (started 8-mar-2010, pjd)
! added calculation of m-resolved differential cross sections for an
!  elastic transition in the geometric apse frame.  this is followed
!  by evaluation of tensor diff. cross sections (started 7-oct-2011, pjd).
!  see statements after "input:" for details of types of cross sections
!  that can be calculated
! took care of special problem with is lables for ibasty = 4 and 19
!  in sigk (NOTE:  this is still a problem in other types of calculations)
!
! revision:  5-oct-2012 (p.j.dagdigian)
! revision:  6-jun-2013 make common/partens as a dynamically allocated
!     module (q. ma)
!---------------------------------------------------------------*
!                                                               *
!  IMPORTANT NOTE:  ONE SHOULD READ THE INPUT FILE FOR THE      *
!  SYSTEM BEFORE RUNNING DIFCRS, SO THAT THE LOGICAL PARAMETERS *
!  FOR THE BASIS TYPE ARE IN MEMORY                             *
!                                                               *
!---------------------------------------------------------------*
use mod_cosout
use mod_codim, only: mmax
use mod_cotble, only: npnt, jttble
use mod_coamat, only: labadr ! labadr(1)
use mod_cobmat, only: sigma ! sigma(kkmx*jmx*jmx)
use mod_cojq, only: jq ! jq(1)
use mod_colq, only: lq ! lq(1)
use mod_coinq, only: inq ! inq(1)
use mod_coisc2, only: inlev => isc2 ! inlev(1)
use mod_coisc3, only: jlev => isc3 ! jlev(1)
use mod_coisc4, only: jpack => isc4 ! jpack(1)
use mod_coisc5, only: lpack => isc5 ! lpack(1)
use mod_coisc6, only: ipack => isc6 ! ipack(1)
use mod_coisc7, only: matel => isc7 ! matel(1)
use mod_coisc8, only: jlist => isc8 ! jlist(1)
! added two common blocks - levels for which xs's to be computed (pjd)
use mod_coisc9, only: jslist => isc9 ! jslist(1)
use mod_coisc10, only: inlist => isc10 ! inlist(1)
use mod_cosc1, only: elev => sc1 ! elev(1)
use mod_cosc2, only: prefac => sc2 ! prefac(1)
use mod_coz, only: sreal => z_as_vec ! sreal(1)
use mod_cow, only: simag => w_as_vec ! simag(1)
use mod_cozmat, only: jtotpa => zmat_as_vec ! jtotpa(1)
use mod_hibrid5, only: sread
use tensor
use constants, only: econv, xmconv, ang2c
use mod_par, only: batch
implicit double precision (a-h,o-z)
character*(*) filnam
character*40  tcsfil, smtfil, tcbfil, dchfil
character*20  cdate
character*10  elaps, cpu
logical csflag, flaghf, flagsu, twomol, exstfl, &
        fast, nucros
!
#include "common/parpot.F90"
common /cospbf/ lnbufs, lnbufl, nbuf, maxlsp, maxllb, ihibuf, &
                igjtp
common /comom/  spin, xj1,xj2, j1, in1, j2, in2, maxjt, maxjot, &
            nwaves, jfsts, jlparf, jlpars, njmax, j1min, j2max
!
common /colnlb/ lenlab(MAX_NJTOT)
common /ckli/  kplist(0:kmx)
common /cf9a/ f9pha(kkmx)
save nout
dimension a(9)
data  tol,   zero,   nstep &
    /1.d-7,  0.d0,     2/
!
call tensor_allocate()
! initialize timer and arays
call mtime(cpu0,ela0)
srealp(lbufs)=zero
simagp(lbufs)=zero
ipackp(lbuflb)=0
jpackp(lbuflb)=0
lpackp(lbuflb)=0
! these 2 parameters are set in partens
lnbufs = lbufs
lnbufl = lbuflb
!
! input:
!  iframe = 0 for quantization axis set to lab Z-axis - isotropic
!    velocity distribution
!  iframe = 1 for quant. axis along initial relative velocity vector
!    (collision frame)
!  iframe = 2 for helicity frame (initial level quantized along initial
!    velocity vector, final level along final velocity vector).  In this
!    case, m-resolved differential cross sections computed for (j1min,in1)
!    -> (j1min,in1) elastic transition, and then tensor cross sections
!    determined and printed out.  The m-resolved differential cross sections
!    are also printed out if iframe is set equal to -2.
!    NOTE:  collision frame differential cross sections (both m-resolved and
!    tensor) can be calculated by setting iframe = +or- 2 and setting ihel
!    equal to zero in amplih subroutine.
!  iframe = 3 for geometric apse frame.  Calculations similar to those for
!   iframe = 2 are carried out.  Likewise, m-resolved differential cross
!   sections are printed out if iframe is set equal to -3.
!
iframe = a(2)
n = a(3)
! NB this variable is called lambda in Follmeg's paper and thesis
! it is the moment of the velocity distribution
if(n.lt.0) then
   minn = iabs(n)
   maxn = minn
else
   minn = 0
   maxn = n
end if
n = minn
maxk=a(4)
in1 = a(5)
!      if(in1.eq.0) in1=1
! previous statement deleted, mha 10/27/95
in2 = a(6)
!      if(in2.eq.0) in2=1
! previous statement deleted, mha 10/27/95
iener = a(1)
maxjot = a(7)
j1min = a(8)
j2max = a(9)
if (iener.le. 0) iener = 1
!
! generate filename and check if it is present
!
call gennam(smtfil,filnam,iener,'smt',lenfs)
inquire(file = smtfil, exist = exstfl)
if (.not. exstfl) then
    write(6,10) smtfil(1:lenfs)
10     format(' ** FILE ',(a),' NOT FOUND **')
    goto 4000
end if
!
! open smatrix-file
!
call openf(1,smtfil,'tu',0)
!
if (abs(iframe) .le. 1) then
! for iframe = 0 and 1, open file for tensor opacities
  call gennam(tcsfil,filnam,iener,'tcs',lenft)
  call openf(2,tcsfil,'sf',0)
  call gennam(tcbfil,filnam,iener,'tcb',lenfb)
! open and rewind .tcb file
  open(form='formatted',unit=4,file=tcbfil)
  rewind (4)
! for iframe = 2 or 3, instead open file for differential cross sections
else
  if (abs(iframe).eq.2) then
    call gennam(dchfil,filnam,iener,'dch',lenft)
    call openf(2,dchfil,'sf',0)
  end if
  if (abs(iframe).eq.3) then
    call gennam(dchfil,filnam,iener,'dcga',lenft)
    call openf(2,dchfil,'sf',0)
  end if
end if
!
! read header of smatrix-file
!
call rdhead(1,cdate,ered,rmu,csflag,flaghf,flagsu,twomol, &
   nucros,jfirst,jfinal,jtotd,numin,numax,nud,nlevel,nlevop,nnout, &
   jlev,inlev,elev,jout)
!
! we need the s-matrices as lower triangles, so nnout  m u s t  be > 0
!
if (nnout.lt.0) then
   write(2,11)
   if(.not.batch) write(6,11)
11    format(' ** NNOUT < 0, ABORT **')
   goto 4000
end if
nout = nnout
!
!  molecule-molecule cross sections not yet implemented
!
if (twomol) then
   write(2,12)
   if (.not. batch) write(6,12)
12    format(' *** TENSOR OPACITIES FOR MOLECULE -', &
          ' MOLECULE COLLISIONS NOT YET IMPLEMENTED ***')
   goto 300
end if
if (flagsu) then
   write(2,14)
   if(.not. batch) write(6,14)
14    format(' *** TENSOR OPACITIES FOR SURFACE -', &
          ' COLLISIONS NOT YET IMPLEMENTED ***')
   goto 300
end if
if (csflag) then
   write(2,16)
   if(.not. batch) write(6,16)
16    format(' *** CS TENSOR OPACITIES', &
          ' NOT YET IMPLEMENTED ***')
   goto 300
end if
!
fast = .true.
spin = 0.
if(flaghf) then
   spin = 0.5d0
! deleted next line - all pairs of jtot and jtotp must be included
! in sums in sigk and sigkkp (pjd)
!         fast = .false.
end if
!
write (2, 20) smtfil, cdate, label, potnam
if(.not. batch) write (6, 20) smtfil, cdate, label,potnam
20 format(/' CLOSE COUPLED TENSOR OPACITIES',/, &
        ' S-MATRICES READ FROM FILE ',(a),/, &
 '      WRITTEN:   ',(a),/, &
 '      LABEL:     ',(a)/, &
 '      POT NAME:  ',(a))
!
! obtain new date
!
call dater(cdate)
write(2, 22) cdate
if(.not. batch) write(6, 22) cdate
22 format(' DATE:    ',(a))
! reset maxjt to jfinal if necessary
if (maxjot.eq.0) maxjot = jfinal
maxjt = min0(jfinal,maxjot)
if (maxjt .lt. matjot) then
  write (2, 23) maxjot, jfinal
  if (.not. batch) write (6, 23) maxjot, jfinal
23   format (' MAX(JTOT) RESET TO JFINAL = ',i4, ' IN TENOPA')
endif
maxjot = maxjt
! nwaves is the number of partial waves
nwaves = jfinal - jfirst + 1
! calculate prefactors and check if j' s are in jout list
! save pointer in array jlist
! also save in values in array inlist (pjd)
nj = 0
jmax = -1
if (j1min.eq.0) j1min = jout(1)
if (j2max.eq.0) j2max = jout(iabs(nnout))
do 40 i=1, iabs(nout)
   jo = jout(i)
   if (jo .gt. j2max) goto 45
   if (jo .lt. j1min) goto 40
   do 30 j=1, nlevel
      j1 = jlev(j)
      in = inlev(j)
! test for both in1 and in2 (pjd)
      if (j1.ne.jo) goto 30
      if (in.ne.in1 .and. in.ne.in2) goto 30
      if (abs(iframe) .ne. 2) then
! test that j1 lies between j1min and j2max for iframe = 0 and 1
        if (j1.lt.j1min .or. j1.gt.j2max) goto 30
! for iframe = 2, make sure (j1min,in1) level is in list
      else
        if (j1 .ne. j1min) goto 30
        if (in .ne. in1) goto 30
      end if
      if (j1 .gt. jmax) jmax = j1
! include only open channels
      if (ered .gt. elev(j)) then
        nj = nj + 1
        jlist(nj) = j
! next two statements added to save state j values and index (pjd)
        jslist(nj) = j1
        inlist(nj) = in
! changed next 2 statements - should index with j, not j1 (pjd)
!              prefac(j1+1) = ang2c * 3.1415926535897932d0 /
!     :                  (2.d0* rmu * (ered - elev(j)))
!              matel(j1+1) = nj
        prefac(nj) = ang2c * 3.1415926535897932d0 / &
                  (2.d0* rmu * (ered - elev(j)))
        matel(nj) = nj
! end of changed statements (pjd)
      endif
30   continue
40 continue
! check if there had been any match
45 if(nj.eq.0) then
   write(2,50)
   if(.not. batch) write(6,50)
50    format(' *** NO TRANSITIONS FOUND, ABORT ***')
   goto 300
end if
!      maxk = 2 * jmax
!      maxk=0
njmax = nj
if(nj.gt.jmx) then
  write(6,51) nj,jmx
51   format(/' TENOPA: NJ =',i4, ' .GT. JMX = ',i4)
  goto 4000
end if
! build up pointer table for direct access i/o
call saddr(1,jfirst,jfinal,numin,numax,csflag,jttble,npnt, &
           jfsts,jlparf,jlpars,ierr)
if (ierr .ne. 0) goto 300
! read s-matrix for jfinal to find out the maximum buffer length
jlp = 1
if (jlparf .eq. 1) jlp = 0
jj = jlp * nwaves + jfinal + 1
iaddr = jttble(jj)
nopen = -1
call sread ( iaddr, sreal, simag, jtot, jlpar, nu, &
            jq, lq, inq, ipack, jpack, lpack, &
             1, mmax, nopen, length, ierr)
maxlsp = (length*(length+1))/2
maxllb = length
nbuf1 = lnbufs/maxlsp
nbuf2 = lnbufl/maxllb
nbuf = min(nbuf1,nbuf2)
! store information on job.tcb file for later transformation into
! m-resolved cross sections
!
if (abs(iframe) .le. 1) then
! write *.tcb as formatted file
! output label and cdate in separate write statements
  write(4, *, err=999) label
  write(4, *, err=999) cdate
  write(4, *, err=999) ered, rmu, flaghf, nlevel, &
                  nlevop, njmax, minn, maxn, nstep, maxk
  write(4, *, err=999) (jlev(i),i=1, nlevel)
  write(4, *, err=999) (inlev(i),i=1, nlevel)
  write(4, *, err=999) (elev(i),i=1, nlevel)
  write(4, *, err=999) (jlist(i),i=1, njmax)
end if
! write header
write(2,60) ered*econv,rmu*xmconv,jfirst,maxjt,maxk,minn,maxn,nbuf
if(.not.batch) write(6,60) &
  ered*econv,rmu*xmconv,jfirst,maxjt,maxk,minn,maxn,nbuf
60 format(/,' ENERGY: ',f11.3,' cm(-1)    MASS: ',f11.3,/, &
     ' SUMMING PARTIAL WAVES FROM JTOT=',i3,' TO JTOT=',i3,/ &
     ' MAX(KI,KF) = ',i3,';',i3,' .LE. LAMDA .LE.',i3, &
     '; MAXIMUM NUMBER OF BUFFERS =',i5)
!
! write level list
write(2,65)
if(.not. batch) write(6,65)
if (abs(iframe) .le. 1)  write(6,64)
64 format(/,' COLUMNS ARE INITIAL STATES; ROWS ARE FINAL STATES')
! switch row and column headings to match integral cross section output
! CORRESPONDING CHANGE IN CREATING XS ARRAYS MADE ONLY IN sigk AND sigkkp SUBRS
65 format (/,' LEVEL LIST FOR TENSOR OPACITIES (OPEN CHANNELS)', &
        /,'   N     J   INDEX  EINT(cm-1)')
do 68 i = 1, nj
   ii=jlist(i)
   if (.not. flaghf) then
     write(2,66) i,jlev(ii),inlev(ii),elev(ii)*econv
     if(.not.batch) write(6,66) i,jlev(ii),inlev(ii), &
                                elev(ii)*econv
66      format (i4, 1x, i5, i6, f11.3)
   else
     write(2,67) i,jlev(ii)+spin,inlev(ii),elev(ii)*econv
     if(.not.batch) write(6,67) i,jlev(ii)+spin,inlev(ii), &
                                elev(ii)*econv
67      format (i4, 1x, f5.1, i6, f11.3)
   endif

68 continue

if (abs(iframe) .le. 1) write(2,64)
! loop over n
n = minn
! fast algorithm if n = 0
70 if (n.eq.0) then
   ifr = abs(iframe) + 1
   goto (171, 172, 173, 174), ifr
     write (6,271) iframe
271      format(/' IFRAME =',i3,' NOT IMPLEMENTED.  ABORT  ***')
     goto 300
! iframe = 0
171        continue
       call sigk(maxk,nnout,jfirst,jfinal,jtotd,nj,mmax,jpack, &
             lpack,ipack,jttble,prefac,sigma, &
             sreal,simag,matel,lenlab,labadr, &
             jtotpa,fast,ierr)
       if(ierr.ne.0) goto 4000
       goto 300
! iframe = 1
172        call sigkc(maxk,nnout,jfirst,jfinal,jtotd,nj,mmax,jpack, &
             lpack,ipack,jttble,prefac, &
             sreal,simag,matel,lenlab,labadr, &
             jtotpa,fast,ierr)
       if(ierr.ne.0) goto 4000
       goto 300
! iframe = 2
173        jlevel = jlev(jlist(1))
       inlevel = inlev(jlist(1))
       call dsigh(maxk,nnout,jfirst,jfinal,jtotd,jpack,lpack, &
             ipack,jlevel,inlevel,elev(jlist(1)),flaghf, &
             iframe,ierr)
       if(ierr.ne.0) goto 4000
       goto 300
! iframe = 3
174        jlevel = jlev(jlist(1))
       inlevel = inlev(jlist(1))
       call dsigga(maxk,nnout,jfirst,jfinal,jtotd,jpack,lpack, &
             ipack,jlevel,inlevel,elev(jlist(1)),flaghf, &
             iframe,ierr)
       if(ierr.ne.0) goto 4000
else
! here for n > 0
   if (abs(iframe) .eq. 1) then
     write (6,72)
72      format ('/, ** N.GT.0 CALCULATION FOR COLL. FRAME NOT', &
       ' IMPLEMENTED.  ABORT  **')
     goto 300
   end if
   nk = 0
   xn = n
   do 75 k = 0, maxk
   xk = k
   do 75 kp = 0, maxk
   xkp = kp
   x = xf3jm0(xkp,xk,xn)
   if (abs(x) .lt. tol) goto 75
   nk = nk + 1
75    continue
   if(nk.gt.kkmx) then
     write(6,76) kkmx,nk
76      format(/' TENOPA: kkmx too small:',2i5)
     goto 4000
   end if
   call sigkkp(n,maxk,nk,nnout,jfirst,jfinal,jtotd,nj,mmax, &
        jpack,lpack,ipack,jttble,prefac,sigma, &
        sreal,simag,matel,lenlab,labadr, &
        jtotpa,kplist,f9pha,fast,ierr)
   if(ierr.ne.0) goto 4000
end if
! next n
n = n + nstep
if (n .le. maxn) goto 70
300 call mtime(cpu1,ela1)
ela1 = ela1 - ela0
cpu1 = cpu1 - cpu0
call gettim(ela1,elaps)
call gettim(cpu1,cpu)
write(2,400) elaps, cpu
if(.not. batch) write(6,400) elaps, cpu
400 format(/,' ** TENXSC FINAL TIMING, ELAPSED: ',a,'  CPU: ',a,' **')
close (1)
close (2)
close (3)
close (4)
call tensor_free()
return
999 write(2,1000)
if(.not.batch) write(6,1000)
1000 format(' *** I/O ERROR IN TENOPA; ABORT')
close (1)
close (2)
close (3)
rewind (4)
close (4)
4000 call tensor_free()
return
end
! ------------------------------------------------------------------
subroutine addsp(jtmin,jtmax,jlp, &
                 labadr,lenlab,jtotpa,jttble)
!
! current revision date: 12-nov-2008 by pj dagdigian
!
use ISO_FORTRAN_ENV, only : ERROR_UNIT
use tensor
use mod_cojq, only: jqp => jq ! jqp(1)
use mod_colq, only: lqp => lq ! lqp(1)
use mod_coinq, only: inqp => inq ! inqp(1)
use mod_coisc4, only: jpack => isc4 ! jpack(1)
use mod_coisc5, only: lpack => isc5 ! lpack(1)
use mod_coisc6, only: ipack => isc6 ! ipack(1)
use mod_hibrid5, only: sread
use mod_par, only: batch
implicit double precision (a-h,o-z)
logical lprnt
common /comom/  spin, xj1,xj2, j1, in1, j2, in2, maxjt, maxjot, &
             nwaves, jfsts, jlparf, jlpars, njmax
common /cospbf/ lnbufs, lnbufl, nbuf, maxlsp, maxllb, ihibuf, &
                igjtp
common /coipar/ ipar(9), iprnt
! add these three common blocks (mha 9/30/08)
dimension labadr(1),lenlab(1),jtotpa(1),jttble(1)
!
! FLAG FOR DIAGNOSTIC PRINTING
lprnt = .false.
if (iprnt .eq. 2) lprnt = .true.
!
ierr = 0
ibuf = ihibuf
! not sure of the function of igjtp.  ignore (pjd)
!      jtpmin = max((igjtp+1),jtmin)
jtpmin = jtmin
jtpmax = jtmax
if(jtmax-jtmin+1.gt.nbuf) then
  write(6,10) jtmax-jtmin+1,nbuf
10   format(/' NOT ENOUGH BUFFERS IN ADDSP:',2i5)
  stop
end if
! set offsets differently than in follmeg's version (pjd)
! by adding lengths of arrays for preceding jtp's
ioffs = 1
ioff = 1
!* DIAGNOSTIC PRINT
   if (lprnt) write (6,1001) j1p,jfsts,nwaves,jlparf, &
      jlpars
1001    format ('entered addsp:  j1p=',i3,' jfsts=',i4, &
     ' nwaves=',i3, ' jlparf=',i3,' jlpars=',i3)
!
do 100 jtp=jtpmin,jtpmax
! calculate address of s-matrix in file
   jjp = jlp * nwaves + jtp + 1
! alternate parity for delta(jtotp-jtot) even/odd (pjd)
   if (((jtp-jtpmin)/2)*2 .ne. (jtp-jtpmin)) then
! delta=odd:  switch parity of jtotp
! first check whether jlpar=1 or -1
      if (jtp+1 .eq. jjp) then
! jlpar=1, skip forward to index for jtotp, jlpar=-1
         jjp = jjp + nwaves
      else
! jlpar=-1:  skip backward to index for jtotp, jlpar=1
         jjp = jjp - nwaves
      end if
   end if
! end of added/modified code (pjd)
   iaddrp = jttble(jjp)
   if(iaddrp .lt. 0) goto 100
! calculate address of s-matrix in buffer
   ibuf = ibuf + 1
   if (ibuf .gt. nbuf) ibuf = 1
! set offsets differently (pjd)
!         ioffs = (ibuf-1)*maxlsp + 1
! replace 35535 below with global parameter lbufs (pjd)
   if (ioffs .gt. lbufs) write (6,11) jtp, ioffs
11    format (' *** JTP = ',i3, '; IOFFS = ', i8)
! set offsets differently (pjd)
!         ioff = (ibuf-1)*maxllb + 1
! replace 1421 below with global parameter lbuflb (pjd)
   if (ioff .gt. lbuflb) write (6,12) jtp, ioff
12    format (' *** JTP = ',i3, '; IOFF = ', i5)
!
! read s-matrix for jtot' into buffer
   nopenp = -1
   call sread ( iaddrp, srealp(ioffs), simagp(ioffs), jtotp, &
                jlparp, nu, jqp, lqp, inqp, ipackp(ioff), &
                jpackp(ioff), lpackp(ioff), &
                1, maxlsp, nopenp, lengtp, ierr)
   if(ierr.eq.-1) goto 999
   if(ierr.lt.-1) then
      write(2,20)
      if(.not.batch) write(6,20)
20       format(' ** READ ERROR IN ADDSP, ABORT **')
      return
   end if
!* DIAGNOSTIC PRINT
   if (lprnt) write (6,2001) jtp,j1p,jjp,jtotp,jlparp, &
       iaddrp
2001    format ('read jtp=',i4,' jlp=',i3,' jjp=',i4, &
      ' jtotp=',i4,' jlpar=',i3,' iaddrp=',i4)
!
! save pointer for s-matrix and for label arrays
   if (jtotp+1 > MAX_NJTOT) then
      write(ERROR_UNIT,*) 'error: max number of jtot values exceeded (please increase MAX_NJTOT)'
      stop 1
   end if

   jtotpa(jtotp+1) = ioffs - 1
   labadr(jtotp+1) = ioff - 1
   lenlab(jtotp+1) = lengtp
! set offsets for next jtp (pjd)
   ioffs = ioffs + lengtp*(lengtp+1)/2
   ioff = ioff + lengtp
100 continue
ihibuf = ibuf
! ignore igjtp (pjd)
!      igjtp = jtotp
!
! here if end of file detected
!
999 return
end
! ------------------------------------------------------------------
subroutine mrcrs(filnam,a)
!
!
! author: b. follmeg
! current revision date:  9-oct-2008 by pjd
!
! subroutine to compute integral fully m-resolved cross sections
! sigma(j1,j2,m1,m2) from tensorial cross sections sigma(n,k,k',j1,j2),
!
!                                    j1+j2-m1-m2                1/2
! sigma(lam,j1,j2,m1,m2) =  sum  (-1)            [(2k+1)(2k'+1)]
!                           k,k'
!
!                   ! j1  j1  k ! ! j2  j2  k'!
!                   !           ! !           ! sigma(lam,k,k',j1,j2)  ,
!                   ! m1 -m1  0 ! ! m2 -m2  0 !
!
!
! partially degeneracy averaged cross sections sigma(lam,j1,j2,m1),
!
! sigma(lam,j1,j2,m1) = sum sigma(lam,j1,j2,m1,m2) ,
!                      m2
!
!
! and fully degeneracy averaged cross sections sigma(lam,j1,j2),
!
!                          -1
! sigma(lam,j1,j2) = (2j1+1)  sum  sigma(lam,j1,j2,m1,m2) .
!                            m1,m2
!
! ------------------------------------------------------------------
use tensor
use mod_codim, only: mmax
use mod_coamat, only: sigmam => toto ! sigmam(1)
use mod_cobmat, only: sigmak => bmat ! sigmak(1)
use mod_coinhl, only: jlev => inhold ! jlev(1)
use mod_coisc2, only: inlev => isc2 ! inlev(1)
use mod_coisc3, only: jlist => isc3 ! jlist(1)
use mod_cosc1, only: elev => sc1 ! elev(1)
use mod_coz, only: xm1lab => z_as_vec ! xm1lab(1)
use mod_cow, only: xm2lab => w_as_vec ! xm2lab(1)
use mod_cozmat, only: sigma => zmat_as_vec ! sigma(1)
use constants, only: econv, xmconv, ang2c
use mod_par, only: batch
implicit double precision(a-h,o-z)
character*(*) filnam
character*40  tcbfil, mcsfil
character*20  cdate
character*10  elaps, cpu
logical flaghf, exstfl
#include "common/parpot.F90"
common /coipar/ ipar(9),iprnt
common /comom/  spin, xj1,xj2, j1, in1, j2, in2, maxjt, maxjot, &
             nwaves, jfsts, jlparf, jlpars, njmax
data zero /0.d0/
!
dimension a(1)
!
! initialize timer
call mtime(cpu0,ela0)
zero =  zero
! input
iener = a(1)
if (iener .le. 0) iener = 1
!
! generate filename and check if it is present
call gennam(tcbfil,filnam,iener,'tcb',lenft)
inquire(file = tcbfil, exist = exstfl)
if (.not. exstfl) then
    write(6,10) tcbfil(1:lenft)
10     format(' ** FILE ',(a),' NOT FOUND **')
    return
end if
!
! converted *.tcb to formatted file (pjd)
!      call openf(1,tcbfil,'su',0)
!      call openf(1,tcbfil,'sf',0)
open(form='formatted',unit=4,file=tcbfil)
print *, 'tcbfil:  ', tcbfil
call gennam(mcsfil,filnam,iener,'mcs',lenm)
print *, 'mcsfil:  ', mcsfil
call openf(2,mcsfil,'sf',0)
! converted *.tcb to formatted file and use unit=4 (pjd)
!      read(4, err=999) label, cdate, ered, rmu, flaghf, nlevel,
!     :                 nlevop, njmax, minn, maxn, nstep, maxk
! read label and cdate in separate statements (pjd)
read(4, 12, err=999) label
12 format(a40)
read(4, 13, err=999) cdate
13 format(a21)
!      read(1, *, err=999) label, cdate, ered, rmu, flaghf, nlevel,
!     :                 nlevop, njmax, minn, maxn, nstep, maxk
read(4, *, err=999) ered, rmu, flaghf, nlevel, &
                 nlevop, njmax, minn, maxn, nstep, maxk
! converted *.tcb to formatted file (pjd)
!      read(1 ,err=999) (jlev(i),i=1, nlevel)
!      read(1 ,err=999) (inlev(i),i=1, nlevel)
!      read(1 ,err=999) (elev(i),i=1, nlevel)
!      read(4, err=999) (jlist(i),i=1, njmax)
read(4, * ,err=999) (jlev(i),i=1, nlevel)
read(4, * ,err=999) (inlev(i),i=1, nlevel)
read(4, * ,err=999) (elev(i),i=1, nlevel)
read(4, *, err=999) (jlist(i),i=1, njmax)
!
spin = 0.
if(flaghf) spin = 0.5
!
write (2, 20) tcbfil, cdate, label, maxk, maxk
if(.not. batch) write (6, 20) tcbfil, cdate, label, maxk, maxk
20 format(/' CLOSE COUPLED M-RESOLVED CROSS SECTIONS',/, &
        ' K K''-MATRICES READ FROM FILE ',(a),/, &
        ' WRITTEN: ',(a),/,' LABEL:   ',(a),/, &
        ' MAX(K, K'') IN PREVIOUS CALCULATION = ',i3,/,/, &
        ' WARNING: M-DEPENDENCE CORRECT ONLY FOR J+J'' .LE.',i3,/, &
        /,' ROWS ARE INITIAL STATES, COLUMNS ARE FINAL STATES',/, &
        ' LAST COLUMN IS SUM OF THE ROW')
! NB the variable n, (minn .le. n .le. maxn in even steps)
! is called lambda in Follmeg's paper and thesis
! it is the moment of the velocity distribution
! loop over transitions
   do 400 i = 1, njmax
      ii = jlist(i)
      j1 = jlev(ii)
      xj1 = j1 + spin
      m1comp = nint(2.d0*xj1 + 1.d0)
      in1 = inlev(ii)
      do 300 j = 1, njmax
         jj = jlist(j)
         j2 = jlev(jj)
         xj2 = j2 + spin
         m2comp = nint(2.d0*xj2 + 1.d0)
         in2 = inlev(jj)
         call sigms(numk,i,j,m1comp,m2comp,minn,maxn, &
                nstep,xm1lab,xm2lab,sigmak,sigmam,mmax,ierr)
         if (ierr .ne. 0) return
300       continue
400     continue
! obtain timing information
call mtime(cpu1,ela1)
cpu1 = cpu1 - cpu0
ela1 = ela1 - ela0
call gettim(ela1,elaps)
call gettim(cpu1,cpu)
write(2,600) elaps, cpu
if(.not. batch) write(6,600) elaps, cpu
600 format(/,' ** MRCRS FINAL TIMING ELAPSED: ', &
         (a),' CPU: ',(a),/)
close(4)
close(2)
return
999 write(2,1000)
write(6,1000)
1000 format(' ** READ ERROR IN MRCRS, ABORT')
close(4)
return
end
! -------------------------------------------------------------------
subroutine sigms(numk,ii,jj,m1comp,m2comp,minn,maxn,nstep, &
                 xm1lab,xm2lab,sigmak,sigmam,mmax,ierr)
! ------------------
! current revision date: 9-oct-1997 by pjd
! ------------------
use mod_coisc4, only: isc1 => isc4 ! isc1(1)
use mod_coisc5, only: isc2 => isc5 ! isc2(1)
use mod_coisc6, only: isc3 => isc6 ! isc3(1)
use mod_coisc7, only: isc4 => isc7 ! isc4(1)
use mod_cosc2, only: sc1 => sc2 ! sc1(1)
use mod_par, only: batch, flaghf
implicit double precision(a-h,o-z)
character*20  cdate
!
#include "common/parpot.F90"
common /comom/  spin, xj1,xj2, j1, in1, j2, in2, maxjt, maxjot, &
             nwaves, jfsts, jlparf, jlpars, njmax
dimension xm1lab(1),xm2lab(1),sigmak(mmax,1),sigmam(mmax,1)
data tol /1.0d-10/
data zero /0.0d0/
ierr = 0
rewind (4)
! dummy read
read(4, 11,err=999) label
11 format(a40)
read(4, 12, err=999) cdate
12 format(a21)
!     read(4, err=999) label, cdate, ered, rmu, flaghf, nlevel,
!    :                 id2, id3, id4, id5, id6, id7
read(4,*,err=999) ered, rmu, flaghf, nlevel, &
                 id2, id3, id4, id5, id6, id7
read(4,*,err=999) (isc1(i),i=1, nlevel)
read(4,*,err=999) (isc2(i),i=1, nlevel)
read(4,*,err=999) (sc1(i),i=1, nlevel)
read(4,*,err=999) (isc3(i),i=1, id3)
! loop over n
do 600 n = minn, maxn, nstep
! clear array
do 30 i = 1, m1comp
do 20 j = 1, m2comp+1
   sigmam(i,j) = zero
20 continue
30 continue
sigmat = 0.
! loop over k and k'
read(4,*, err=999, end=400) numk
do 300 kk = 1, numk
   read(4, *, err = 999, end = 400) k,kp
   xk = k
   xkp = kp
   fack = sqrt((2.d0*xk+1.d0)*(2.d0*xkp+1.d0))
! read matrix for k,k' values
   do 50 i = 1, njmax
! switched row and column index in SIGK and SIGKKP (pjd)
!50       read(4,*, err=999, end = 400) (sigmak(i,j), j = 1, njmax)
50    read(4,*, err=999, end = 400) (sigmak(j,i), j = 1, njmax)
   trans = sigmak(ii,jj)
   if (abs(trans) .lt. tol) goto 300
! loop over components
   xm1 = -xj1
   do 200 m1 = 1, m1comp
      xm1lab(m1) = xm1
      xm2 = -xj2
      do 100 m2 = 1, m2comp
          xm2lab(m2) = xm2
! evaluate phase
           phase = 1.d0
           ipower = nint(xj1+xj2-xm1-xm2)
           if((ipower/2)*2 .ne. ipower) phase = -1.d0
           f3ja = xf3j(xj1,xj1,xk,xm1,-xm1,zero)
           f3jb = xf3j(xj2,xj2,xkp,xm2,-xm2,zero)
           s = phase * fack * f3ja * f3jb * trans
           sigmam(m1,m2) = sigmam(m1,m2) + s
           sigmam(m1,m2comp+1) = sigmam(m1,m2comp+1) + s
           sigmat = sigmat + s
100       xm2 = xm2 + 1.d0
200    xm1 = xm1 + 1.d0
300 continue
! print out results
400 continue
if (flaghf) then
   write(2, 410) xj1, in1, xj2, in2, n, sigmat/(2.d0*xj1+1.d0)
   if (.not. batch) &
    write(6, 410) xj1, in1, xj2, in2, n, sigmat/(2.d0*xj1+1.d0)
410     format &
      (/' TRANSITION J1 = ',f4.1,i3,' -> J2 = ',f4.1,i3, &
      ', LAM = ',i2,', TOTAL CROSS SECTION = ',1pe10.3,/)
else
   write(2, 415) nint(xj1), nint(xj2),n, sigmat/(2.d0*xj1+1.d0)
   if (.not. batch) &
   write(6, 415) nint(xj1), nint(xj2),n, sigmat/(2.d0*xj1+1.d0)
415    format &
   (/' TRANSITION J1 = ',i2,' -> J2 = ',i2,', LAM = ',i2, &
       ' ,TOTAL CROSS SECTION = ',1pe10.3,/)
endif
lmax = 0
420 lmin = lmax + 1
lmax = lmax + 9
lmax = min0(lmax,m2comp)
if (flaghf) then
  write(2,430) (xm2lab(l),l=lmin,lmax)
  if(.not. batch) write(6,430) (xm2lab(l),l=lmin,lmax)
430   format(10x,9(f5.1,6x))
else
  write(2,435) (nint(xm2lab(l)),l=lmin,lmax)
  if(.not. batch) write(6,435) (nint(xm2lab(l)),l=lmin,lmax)
435   format(8x,9(i5,6x))
endif
do 500 m1 = 1, m1comp
    xm1 = xm1lab(m1)
    if (lmax .eq. m2comp) then
       if (flaghf) then
         write(2,440) xm1,(sigmam(m1,l),l=lmin,lmax), &
                         sigmam(m1,m2comp+1)
         if (.not. batch) write(6,440) &
                        xm1,(sigmam(m1,l),l=lmin,lmax), &
                         sigmam(m1,m2comp+1)
440          format(1x,f4.1,12(1pe11.3))
       else
         write(2,445) nint(xm1),(sigmam(m1,l),l=lmin,lmax), &
                         sigmam(m1,m2comp+1)
         if (.not. batch) write(6,445) &
                        nint(xm1),(sigmam(m1,l),l=lmin,lmax), &
                         sigmam(m1,m2comp+1)
445          format(1x,i4,12(1pe11.3))
       endif
    else
       if (flaghf) then
         write(2,440) xm1,(sigmam(m1,l),l=lmin,lmax)
         if (.not. batch) write(6,440) &
                    xm1,(sigmam(m1,l),l=lmin,lmax)
       else
         write(2,445) nint(xm1),(sigmam(m1,l),l=lmin,lmax)
         if (.not. batch) write(6,445) &
                    nint(xm1),(sigmam(m1,l),l=lmin,lmax)
       endif
    end if
500  continue
 write(2,'(a)') ' '
 if(.not. batch) write(6,'(a)') ' '
 if((m2comp-lmax)) 600,600,420
600  continue
 return
999 write(2,1000)
write(6,1000)
1000 format(' ** READ ERROR IN SIGMS, ABORT')
ierr = 1
return
end
!------- -----------------------------------------------------------------
subroutine sigk(maxk,nnout,jfirst,jfinal,jtotd,nj,mmax,jpack, &
                lpack,ipack,jttble,prefac,sigma, &
                sreal,simag,matel,lenlab,labadr, &
                jtotpa,fast,ierr)
!
! subroutine to calculate sigma(k,j1,j2) cross sections:
! ( see also " m.h. alexander and s.l. davis, jcp 78(11),6748(1983)"
!   equation 23 )
!
!         k          pi                                l1+l2-j1-j2+2j
!   sigma       =  -----      sum   (2j+1)(2j'+1)  (-1)
!         j1,j2    k(j1)^2  j,j',l1,l2
!
!             ! j1 j1 k ! ! j2 j2 k !    j             * j'
!             {         } {         }  t              t
!             ! j  j' l1! ! J  J' l2!    j1,l1,j2,l2    j1,l1,j2,l2
!
! author: b. follmeg
! current revision date: 5-oct-2012 by pj dagdigian
!
!------------------------------------------------------------------------
use tensor
use mod_cojq, only: jq ! jq(1)
use mod_colq, only: lq ! lq(1)
use mod_coinq, only: inq ! inq(1)
! modules for levels for which xs's to be computed
use mod_coisc9, only: jslist => isc9 ! jslist(1)
use mod_coisc10, only: inlist => isc10 ! inlist(1)
use mod_hibrid2, only: mxoutd
use mod_hibrid5, only: sread
use mod_par, only: batch, ipos
implicit double precision (a-h,o-z)
complex*8 t, tp
logical diag, diagj, diagin, &
        twopar, fast

!* flags for diagnostic printing
logical lprnt,lprntf

character*10 elaps, cpu
common /comom/  spin, xj1,xj2, j1, in1, j2, in2, maxjt, maxjot, &
            nwaves, jfsts, jlparf, jlpars, njmax, j1min, j2max
common /cospbf/ lnbufs, lnbufl, nbuf, maxlsp, maxllb, ihibuf, &
                igjtp
common /coipar/ ipar(9),iprnt
common /coselb/ ibasty
! 3rd subscript is for state index (subscript = 5 + IN)
! states with up to 9 state indices allowed
common/cadr/ iadr(0:2*jmx,lmx,9)
common/c6jt/ f6a(kmx,0:2*jmx,lmx),f6p(kmx)
!
dimension jpack(1),ipack(1),lpack(1)
dimension sreal(1), simag(1)
dimension prefac(1), matel(1), labadr(1), jtotpa(1)
dimension lenlab(1), jttble(1)
dimension sigma(nj,nj,maxk+1)
! array to store partial-j cross sections (K=0-2 only)
dimension psig0(1001),psig1(1001),psig2(1001)
!
! FLAG FOR DIAGNOSTIC PRINTING
lprnt = .false.
if (iprnt .eq. 2) lprnt = .true.

! FLAG FOR PARTIAL-Jtot CROSS SECTIONS
lprntf = .false.
if (iprnt .eq. 3) lprntf = .true.
!
! initialize timer
call mtime(cpu0,ela0)
ierr = 0
! clear sigma array
do 6 k = 0, maxk
do 5 j = 1, njmax
do 4 i = 1, njmax
  sigma(i,j,k+1) = 0.d0
4 continue
5 continue
6 continue
if(maxk.gt.kmx.or.j2max.gt.jmx) then
  write(6,7) maxk,kmx,j2max,jmx
7   format(/' SIGK: jmx or kmx too small:',4i5)
  ierr=-1
  return
end if
diagin = .false.
if (in1 .eq. in2) diagin = .true.
twopar = .false.
if (jlpars.ne.0) twopar = .true.
jt = jfirst
jstart = jfirst
jlp = 1
if (jlparf .eq. 1) jlp = 0
jjoff = jlp * nwaves + 1
ihibuf = 0
igjtp = -1
jtpstp = 1
!
! initialize arrays to store partial-jtot cross sections (K=0-2 only)
if (lprntf) then
  do 1006 jj = 1, 1001
  psig0(jj) = 0.d0
  psig1(jj) = 0.d0
  psig2(jj) = 0.d0
1006   continue
end if
!
! sum over jtot and over parities (if twopar = .true.)
10 continue
!
! find address of jtot in smatrix file
jj = jjoff + jt
iaddr = jttble(jj)
if(iaddr .lt. 0) goto 700
! read s-matrix for jtot
nopen = -1
call sread ( iaddr, sreal, simag, jtot, jlpar, nu, &
            jq, lq, inq, ipack, jpack, lpack, &
             1, mmax, nopen, length, ierr)
if(ierr.eq.-1) goto 999
if(ierr.lt.-1) then
  write(2,20)
  if(.not.batch) write(6,20)
20   format(' *** READ ERROR, ABORT')
  return
end if
xjtot = jtot + spin
facj = (2.d0 * xjtot + 1.d0)
!
! boundaries for sum over jtot'
jtpmin = jtot
jtpmax = min((jtot + maxk), jfinal)
! set parity to be considered for jtot' (pjd)
jlparp = jlpar
!
! fill buffer with required s' matrices
! parity for each jtot' needs to be kept track of in addsp
call addsp(jtpmin,jtpmax,jlp, &
           labadr,lenlab,jtotpa,jttble)
if (srealp(lbufs) .ne. 0.d0) print *, 'srealp error in sigk'
if (simagp(lbufs) .ne. 0.d0) print *, 'simagp error in sigk'
if (ipackp(lbuflb) .ne. 0) print *, 'ipackp error in sigk'
if (lpackp(lbuflb) .ne. 0) print *, 'lpackp error in sigk'
if (jpackp(lbuflb) .ne. 0) print *, 'jpackp error in sigk'
! prepare for sum over jtot'
jtotp = jtpmin
!
! start of loop for sum over jtot'
60 continue
! set parity to be considered for jtotp:
!   jlparp = jlpar for delta(jtot - jtotp) = even
!   jlparp = -jlpar for delta () = odd
jlparp = jlpar
if (((jtotp-jtot)/2) *2 .ne. (jtotp-jtot)) jlparp = -jlpar
! switch=1 to include (jtot,jtotp) and (jtotp,jtot) terms for jtot .ne. jtotp
switch = 1.d0
if (jtotp .eq. jtot .or. jtotp .gt. maxjt) switch = 0.
! find address of jtot' in s-matrix buffer
jj = jtotp + 1
! alternate parity for delta(jtotp-jtot) even/odd
if (((jtotp-jtot)/2)*2 .ne. (jtotp-jtot)) then
! delta=odd:  switch parity of jtotp
! first check whether jlpar=1 or -1
  if (jjoff .eq. 1) then
! jlpar=1, skip forward to index for jtotp, jlpar=-1
     jjp = jjp + nwaves
  else
! jlpar=-1:  skip backward to index for jtotp, jlpar=1
     jjp = jjp - nwaves
  end if
end if

ioffs = jtotpa(jj)
ioff  = labadr(jj)
lengtp = lenlab(jj)
xjtotp = jtotp + spin
facjjp = (2.d0 * xjtotp + 1.d0) * facj
jminjp = iabs(jtot-jtotp)
jplujp = nint(xjtot + xjtotp)
kmin=jminjp
jpmax = 0
do irowp = 1,lengtp
  jpmax = max(jpackp(ioff+irowp),jpmax)
end do
! now loop over all transitions
lmax=jtotp+jpmax+1
! clear iadr array
do 69 j=0,jmx
do 68 l=1,lmx
do 67 index=1,9
  iadr(j,l,index) = 0
67 continue
68 continue
69 continue
!
do 70 irowp = 1, lengtp
  irp = ioff + irowp
  j1p = jpackp(irp)
  l1p = lpackp(irp)
  indp = ipackp(irp)
  if (ipackp(irp).ne.in1 .and. ipackp(irp).ne.in2) goto 70
!
! special treatment for ibasty = 4 and 19
  if (iabsty.eq.4 .or. ibasty.eq.19) then
    if (mod(abs(indp),100).eq.10) then
      write (6,1075)
1075       format(' *** v=1 IN PI STATE NOT IMPLEMENTED.', &
        ' ABORT ***')
      return
    end if
    if (indp.eq. 100) isub = 1
    if (indp.eq.-100) isub = 2
    if (indp.eq. 200) isub = 3
    if (indp.eq.-200) isub = 4
    if (indp.eq. 300) isub = 5
    if (indp.eq.-300) isub = 6
    iadr(j1p,lmax-l1p,isub) = irowp
! all other basis types
  else
! check that kp < 4
    if (abs(indp).gt.4) then
      write (6,1074) abs(indp)
1074       format(' *** KP =',i3,' > 4. ABORT ***')
      return
    end if
    iadr(j1p,lmax-l1p,5+indp) = irowp
  end if
70 continue
!
jmx1=0
!
! sum over row index for jtot
do 400 irow = 1, length
  j1 = jpack(irow)
  xj1 = j1 + spin
  if (j1 .gt. j2max ) goto 400
  if (j1 .gt. jpmax ) goto 400
  if (j1 .lt. j1min ) goto 400
  l1 = lpack(irow)
  xl1 = l1
  if (xl1.lt.abs(xjtotp-xj1) .or. xl1.gt.xjtotp+xj1) goto 400
  j1t2 = nint(2.d0*xj1)
  do 74 ij = 1, njmax
  if (j1.ne.jslist(ij)) goto 74
  if (ipack(irow).ne.inlist(ij)) goto 74
  ir = ij
  denrow = prefac(ij)
  goto 75
74   continue
! level not in list
  ir = 100
75   continue
  jmx1=max(j1,jmx1)
  xjmx1=jmx1+spin
  kmax=min(jplujp,int(2*xjmx1))
! kmax must be <= maxk of tensor xs's to be calculated (pjd)
  kmax=min(kmax,maxk)
  kk=0
  l1m=lmax-l1
  do 80 k=kmin,kmax
    xk=k
    kk=kk+1
    if (kk.gt.kmx) write (6,*) 'kk error in sigk'
    if (j1.gt.jmx) write (6,*) 'j1 error in sigk'
    if (l1m.gt.lmx) write (6,*) 'l1m error in sigk'
    f6a(kk,j1,l1m) = xf6j(xj1,xj1,xk,xjtot,xjtotp,xl1)
80   continue
!
! sum over column index for jtot
  do 200 icol = 1, irow
    if (ipack(irow).eq.in1 .and. ipack(icol).eq.in2) &
       go to 82
    if (ipack(irow).eq.in2 .and. ipack(icol).eq.in1) &
       go to 82
    go to 200
82     continue
    j2 = jpack(icol)
    xj2 = j2 + spin
    if (j2 .gt. j2max) goto 200
    if (j2 .gt. jpmax) goto 200
    if (j2 .lt. j1min) goto 200
    l2 = lpack(icol)
    xl2 = l2
    if (xl2.lt.abs(xjtotp-xj2) .or. xl2.gt.xjtotp+xj2) goto 200
    diagj = j1.eq.j2
    j2t2 = nint(2.d0*xj2)
    min2j = min(j1t2,j2t2)
    if (min2j .lt. jminjp) goto 200
    do 84 ij = 1, njmax
      if (j2.ne.jslist(ij)) goto 84
      if (ipack(icol).ne.inlist(ij)) goto 84
      ic = ij
      dencol = prefac(ij)
      goto 85
84     continue
! level not in list
    ic = 100
85     continue
    ii = (irow*(irow-1))/2 + icol
    l2m=lmax-l2
! special treament for iabsty = 4 and 19
    if (ibasty.eq.4 .or. ibasty.eq.19) then
      if (ipack(irow).eq. 100) isubr = 1
      if (ipack(irow).eq.-100) isubr = 2
      if (ipack(irow).eq. 200) isubr = 3
      if (ipack(irow).eq.-200) isubr = 4
      if (ipack(irow).eq. 300) isubr = 5
      if (ipack(irow).eq.-300) isubr = 6
      if (ipack(icol).eq. 100) isubc = 1
      if (ipack(icol).eq.-100) isubc = 2
      if (ipack(icol).eq. 200) isubc = 3
      if (ipack(icol).eq.-200) isubc = 4
      if (ipack(icol).eq. 300) isubc = 5
      if (ipack(icol).eq.-300) isubc = 6
      irowp=iadr(j1,l1m,isubr)
      icolp=iadr(j2,l2m,isubc)
! all other basis types
    else
      irowp=iadr(j1,l1m,5+ipack(irow))
      icolp=iadr(j2,l2m,5+ipack(icol))
    end if
! if l1p <> l1 or l2p <> l2, do not include this (jtot,jtotp) term
    if (irowp.eq.0 .or. icolp.eq.0) goto 200
    if(irowp.ge.icolp) then
      iip = ioffs + (irowp*(irowp-1))/2 + icolp
    else
      iip = ioffs + (icolp*(icolp-1))/2 + irowp
    end if
    diag = diagj .and. diagin
    if (diag) dencol = 0.
    wd = 1.d0
    if (diag .and. l1.ne.l2) wd = 2.d0
! phase factor
    phasea = 1.d0
    ipower= nint(xl1+xl2-xj1-xj2+2.d0*xjtot)
    if ((ipower/2)*2 .ne. ipower) phasea = -1.d0
    phaseb = 1.d0
    ipower= ipower+2*nint(xjtotp-xjtot)
    if ((ipower/2)*2 .ne. ipower) phaseb = -1.d0
! convert s-matrix to t-matrix: t(j,j1,in1,l1 ; j2,in2,l2) =
!     delta(j1,in1,l1 ; j2,in2,l2) - s(j,j1,in1,l1 ; j2,in2,l2)
    t  = -cmplx(sreal(ii),simag(ii))
    tp = -cmplx(srealp(iip),simagp(iip))
    if (iip .gt. lbufs) write (6,11)iip
11     format (' *** IIP = ', i8)
    if (diag .and. l1.eq.l2) then
      t = t + (1.d0,0)
      tp = tp + (1.d0,0)
    end if
! note below that real(t*conjg(tp)) = real(tp*conjg(t))
! and imag(t*conjg(tp)) = -imag(tp*conjg(t))
    factor = wd * facjjp * real(t*conjg(tp))
! terms in parentheses on RHS of next line is to include (jtot,jtotp) and (jtotp,jtot) terms
    factor = factor * (phasea + phaseb*switch)
! loop over k
    kmax = min(min2j,jplujp)
! kmax must be <= maxk of tensor xs's to be calculated
    kmax=min(kmax,maxk)
    kk=0
    do 90 k = kmin, kmax
      kk=kk+1
      f6p(kk) = f6a(kk,j1,l1m)*f6a(kk,j2,l2m)*factor
90     continue
    kk=0
    do 100 k = kmin, kmax
      kk=kk+1
      sigma(ir,ic,k+1) = sigma(ir,ic,k+1) + f6p(kk)*dencol
      sigma(ic,ir,k+1) = sigma(ic,ir,k+1) + f6p(kk)*denrow
!
! accummulate partial cross sections
      if (lprntf) then
! print values for only one transition
! the number of the states are obtained from the list of states
! in tenxsc output
        istatei = 4
        istatef = 1
        if ((ir.eq.istatei .and. ic.eq.istatef) .or. &
             (ir.eq.istatef .and. ic.eq.istatei)) then
          if (k .eq. 0) psig0(jtot+1) = psig0(jtot+1) &
             + f6p(kk)*denrow
          if (k .eq. 1) psig1(jtot+1) = psig1(jtot+1) &
             + f6p(kk)*denrow
          if (k .eq. 2) psig2(jtot+1) = psig2(jtot+1) &
             + f6p(kk)*denrow
        end if
      end if
100     continue
!
200   continue
400 continue
!
! next jtot'
500 jtotp = jtotp + jtpstp
if(jtotp.le.jtpmax) goto 60
!
! next jtot
700 jt = jt + 1
if (jt.le.maxjt) goto 10
!
! next parity
if (twopar) then
  twopar = .false.
  jt = jfsts
  jstart = jfsts
  jlp = 1
  jjoff = jlp * nwaves + 1
  ihibuf = 0
  igjtp = -1
! loop back for next parity
  goto 10
end if
! here if calculation has been done
! write *.tcb as formatted file
999 write(4, *) maxk + 1
do 1100 k = 0, maxk
  write(4, *) k,k
  write(2,800) k
  if(.not.batch) write(6,800) k
800   format(/' TENSOR RANK K =',i3)
  do 1000 i = 1, njmax
! write *.tcb as formattted file (pjd)
    write(4, *) (sigma(i,j,k+1),j=1,njmax)
1000   continue
  call mxoutd (2, sigma(1,1,k+1), njmax, njmax, 0, ipos)
! use diag as scratch variable (ipos = .false. for screen output)
  diag = .false.
  if(.not. batch) &
    call mxoutd (6, sigma(1,1,k+1), njmax, njmax, 0, diag)
1100 continue

! PRINT OUT PARTIAL CROSS SECTIONS
if (lprntf) then
  write (6,1011) istatei,istatef
1011   format (/' PARTIAL-JTOT CROSS SECTIONS for transition: ', &
     '  level = ',i2,' -> ',' level = ',i2/ &
     '   jtot',7x,'P-XS(0)',8x,'P-XS(1)',8x,'P-XS(2)')
  do jj = 0,maxjt
    write (6,1012) jj,psig0(jj+1),psig1(jj+1), &
       psig2(jj+1)
1012     format(i6,2x,3e15.6)
  end do
end if

call mtime(cpu1,ela1)
cpu1 = cpu1 - cpu0
ela1 = ela1 - ela0
call gettim(ela1,elaps)
call gettim(cpu1,cpu)
write(2,1200) elaps, cpu
if(.not. batch) write(6,1200) elaps, cpu
1200 format(/' ** N = 0 COMPLETED, TIMING ELAPSED: ',a, &
        ' CPU: ',a,/)
return
end
!------------------------------------------------------------------------
subroutine sigkkp(n,maxk,nk,nnout,jfirst,jfinal,jtotd,nj,mmax, &
         jpack,lpack,ipack,jttble,prefac,sigma, &
         sreal,simag,matel,lenlab,labadr, &
         jtotpa,kplist,f9pha,fast,ierr)
!
! subroutine to calculate sigma(lambda; j1, j2, ki, kf) cross section
! defined by follmeg et al., jcp 93(7), 4687 (1990).
! NOTE:  lambda = n in subr.
!
! author: b. follmeg
! previous revision date: 18-apr-1997 by mha
! current revision date: 9-oct-2008 by pj dagdigian
!
!** subr only partially revised to calculate cross sections correctly
!** for half-integer j and even multiplicity open-shell molecules
!** DO NOT USE TO CALCULATE sigma(k,kp) TENSOR CROSS SECTIONS
!** REVISIONS NOT COMPLETED (pjd)
!
!------------------------------------------------------------------------
use tensor
use mod_cojq, only: jq ! jq(1)
use mod_colq, only: lq ! lq(1)
use mod_coinq, only: inq ! inq(1)
! added two common blocks - levels for which xs's to be computed (pjd)
use mod_coisc9, only: jslist => isc9 ! jslist(1)
use mod_coisc10, only: inlist => isc10 ! inlist(1)
use mod_hibrid5, only: sread
use mod_par, only: batch, ipos
implicit double precision (a-h,o-z)
complex*8 t, tp, ai, cphase
logical diag, diagj, diagin, &
        twopar, fast

!* flag for diagnostic printing
logical lprnt2

character*10 elaps, cpu
common /comom/  spin, xj1,xj2, j1, in1, j2, in2, maxjt, maxjot, &
            nwaves, jfsts, jlparf, jlpars, njmax, j1min, j2max
common /cospbf/ lnbufs, lnbufl, nbuf, ihibuf,  maxlsp, maxllb, &
               igjtp
common /coipar/ ipar(9),iprnt
common/cccpu/ tchi,tlog,tsetup,tdelt,lenk,max1,max2,max3,maxkk
! add 3rd subscript for state index (subscript = 5 - IN) (pjd)
! states with up to 9 state indices allowed
common/cadr/ iadr(0:2*jmx,lmx,9)
!      common/cadr/ iadr(0:jmx,lmx)
common/c6jt/ f6a(kmx,0:2*jmx,lmx),f9a(3*kmx,lmx)

dimension jpack(1),ipack(1),lpack(1)
dimension sreal(1), simag(1)
dimension prefac(1), matel(1), labadr(1), jtotpa(1)
dimension lenlab(1), jttble(1), kplist(0:maxk)
dimension sigma(nk,nj,nj), f9pha(nk)
data tol /1.d-10/
data zero /0.d0/
#if defined(HIB_UNIX)
data ai / (0.d0,1.d0)/
#endif
#if defined(HIB_CRAY)
data ai / (0.d0,1.d0)/
#endif

!* flag for diagnostic printing
lprnt2 = .false.
if (iprnt.eq.4) lprnt2 = .true.

! initialize timer
call mtime(cpu0,ela0)
t6j=0
t9j=0
ierr = 0
icount = 0
xn = n
! check limits of sigma array size
if(maxk.gt.kmx-1) then
  write(6,4) kmx,maxk
4   format(/' SIGKKP: KMX too small:',i4,'  need:',i4)
  ierr=-1
  return
end if
if(j2max.gt.jmx) then
  write(6,5) jmx,j2max
5   format(/' SIGKKP: JMX too small:',i4,'  need:',i4)
  ierr=-1
  return
end if
nkk=0
iof=0
do 6 kp=0,maxk
kplist(kp)=nkk
do 6 k=abs(kp-n),min(maxk,kp+n),n
nkk=nkk+1
if (nkk .gt. nk) print *, 'nkk error',kp, k, nk
!
! clear sigma array
do 6 i=1,nj
do 6 j=1,nj
6 sigma(nkk,j,i)=zero
!
!      if(nkk.ne.nk) then
!        print *, maxk, n,k,kp,nkk
!        write(6,7) nk,nkk
!7       format(/' nk ne nkk:',2i4)
!        ierr=-1
!        return
!      end if
if(nkk.gt.kmx*3) then
  write(6,8) kmx**2,nkk
8   format(/' SIGKKP: KMX**2 too small:',i4,'  need:',i4)
  ierr=-1
  return
end if
diagin = .false.
if (in1 .eq. in2) diagin = .true.
twopar = .false.
if(jlpars.ne.0) twopar = .true.
jt = jfirst
jstart = jfirst
jlp = 1
if (jlparf .eq. 1) jlp = 0
jjoff = jlp * nwaves + 1
ihibuf = 0
igjtp = -1
! for integer jtot we need only even/even or odd/odd combinations of
! jtot/jtotp
jtpstp = 1
if (fast) jtpstp = 2
! sum over jtot and over parities (if twopar = .true.)
10 continue
! find address of jtot in smatrix file
jj = jjoff + jt
iaddr = jttble(jj)
if(iaddr .lt. 0) goto 700
! read s-matrix for jtot
nopen = -1
call sread ( iaddr, sreal, simag, jtot, jlpar, nu, &
            jq, lq, inq, ipack, jpack, lpack, &
             1, mmax, nopen, length, ierr)
if(ierr.eq.-1) goto 999
if(ierr.lt.-1) then
  write(2,20)
  if(.not.batch) write(6,20)
20   format(' *** READ ERROR IN SIGKKP, ABORT')
  return
end if
xjtot = jtot + spin
facj = (2.d0* xjtot + 1.d0)
!
! boundaries for sum over jtot'
jtpmin = jtot
jtpmax = min((jtot + maxk), jfinal)
! fill buffer with required s' matrices
call addsp(jtpmin,jtpmax,jlp, &
           labadr,lenlab,jtotpa,jttble)
!
! prepare for sum over jtot'
jtotp = jtpmin
!
! start of loop for sum over jtot'
60 continue
switch = 1.d0
if (jtot .eq. jtotp .or. jtotp .gt. maxjt) switch = 0.
! find address of jtot' in s-matrix buffer
jj = jtotp + 1
ioffs = jtotpa(jj)
ioff  = labadr(jj)
lengtp = lenlab(jj)
xjtotp = jtotp + spin
facjjp = (2.d0* xjtotp + 1.d0) * facj
jminjp = iabs(jtot-jtotp)
jplujp = nint(xjtot + xjtotp)
! precompute 6j symbols
!
! next 2 statements assume last basis level for jtot' has highest j
! replace with scan over all jtot' levels (pjd)
!      j2mxp=jpackp(ioff+lengtp)
!      j2mx =jpack(length)
j2mxp = 0
do irowp = 1,lengtp
  j2mxp = max(jpackp(ioff+irowp),j2mxp)
end do
io = labadr(jtot + 1)
j2mx = 0
do irow = 1,length
  j2mx = max(jpack(io+irow),j2mx)
end do
! end of addition to set j2mxp and j2mx (pjd)
! modify next 2 statements for half-integral j (pjd)
!      lmax =jtot+j2mx+1
!      lmaxp=jtotp+j2mxp+1
lmax = nint(xjtot+j2mx+spin+1.d0)
lmaxp = nint(xjtotp+j2mxp+spin+1.d0)
kpmin=jminjp
kmin=abs(kpmin-n)
t1=second()
do 70 icol=1,length
! modify next if statement in include test for in2 (pjd)
!      if(ipack(icol).ne.in2) goto 70
if(ipack(icol).ne.in2 .and. ipack(icol).ne.in1) goto 70
j2=jpack(icol)
l2=lpack(icol)
xj2=j2+spin
xl2=l2
kpmax=min(jplujp,nint(2*xj2))
kkp=0
do 65 kp=kpmin,kpmax
kkp=kkp+1
xkp=kp
65 f6a(kkp,j2,lmax-l2) = xf6j(xkp,xjtotp,xjtot,xl2,xj2,xj2)
70 continue
t6j=t6j+second()-t1
kpmx=kpmax
! clear iadr array (pjd)
do 79 j=0,jmx
do 78 l=1,lmx
do 77 index=1,9
iadr(j,l,index) = 0
77 continue
78 continue
79 continue
do 80 icolp=1,lengtp
irp = ioff + icolp
ind = ipackp(irp)
! moved next 2 statements before if statement (pjd)
j1=jpackp(irp)
l1=lpackp(irp)
! modify if statement to include test for in2 (pjd)
!      if(ipackp(irp).ne.in1) goto 80
if(ipackp(irp).ne.in1 .and. ipackp(irp).ne.in2) goto 80
! added 3rd subscript in next statement (pjd)
iadr(j1,lmaxp-l1,5+ind)=icolp
80 continue
!
! now loop over all transitions
do 400 irow = 1, length
   if(ipack(irow).ne.in1) goto 400
   j1 = jpack(irow)
! moved next statement ahead of if statements (pjd)
   xj1 = j1 + spin
   if (j1 .gt. j2max) goto 500
   if (j1 .gt. j2mxp) goto 500
   if (j1 .lt. j1min) goto 400
   l1 = lpack(irow)
! moved next statement ahead of if statements (pjd)
   xl1 = l1
   l1pmin=max(abs(l1-n),abs(j1-jtotp))
! use actual angular momenta, not integer values, in next statement (pjd)
!         if((-1)**(l1pmin+j1-jtotp)).ne.jlpar) l1pmin=l1pmin+1
   if((-1)**nint((l1pmin+xj1-xjtotp)).ne.jlpar) l1pmin=l1pmin+1
   l1pmax=min(l1+n,j1+jtotp)
   if(l1pmax.lt.l1pmin) goto 400
   j1t2 = nint(2.d0*xj1)
! change index for matel and prefac (pjd)
!         ir = matel(j1+1)
!         denrow = prefac(j1+1)
! replace previous 2 statements (pjd)
   do 84 ij=1,njmax
   if (j1.ne.jslist(ij)) goto 84
   if (ipack(irow).ne.inlist(ij)) goto 84
   ir = ij
   denrow = prefac(ij)
   goto 85
84    continue
! level not in list (pjd)
   ir = 100
85    continue
! end of replacement (pjd)
   xxj=xjtot+xj1+xjtotp+xj1
   kmax=nint(2*xj1)
   factor=denrow*facjjp
   ll=0
   t1=second()
do 100 l1p=l1pmin,l1pmax,2
   ll=ll+1
   xl1p = l1p
   cphase = 1.d0
   ipower = l1 + l1p
   if (ipower .ne. 0) cphase = ai ** ipower
   f3a = real(cphase)*xf3jm0(xl1,xl1p,xn) &
      *sqrt((2.d0*xl1+1.d0)*(2.d0*xl1p+1.d0))*factor
! precalculate 9j symbols and phase factors
   do 90 kp=kpmin,kpmx
   kmin1=abs(kp-n)
   kmax1=min(kp+n,kmax)
   ik=kplist(kp)
   xkp = kp
   do 90 k=kmin1,kmax1,n
   xk = k
   ik=ik+1
   f9a(ik,ll) =  f3a*(-1)**(k+kp) &
            * sqrt((2.d0*xk+1.d0)*(2.d0*xkp+1.d0)) &
            * xf3jm0(xkp,xk,xn) &
!     :            * xf9j(xkp,xk,xn,xjtot,xj1,xl1,xjtotp,xj1,xl1p)
             * xf9j(xn,xkp,xk,xl1,xjtot,xj1,xl1p,xjtotp,xj1)
   f9pha(ik) = switch
   ipower = nint(xk+xkp+xxj)
   if(mod(ipower,2).ne.0) f9pha(ik) = -switch
90 continue
100 continue
   t9j=t9j+second()-t1
!
! compute all transitions separately, unlike in SIGK
do 300 icol = 1, length
   if(ipack(icol).ne.in2) goto 300
   j2 = jpack(icol)
   xj2 = j2 + spin
   if (j2 .gt. j2max) goto 400
   if (j2 .gt. j2mxp) goto 400
   if (j2 .lt. j1min) goto 300
   l2 = lpack(icol)
   xl2 = l2
! next if statement tested on actual angular momenta, not integer values (pjd)
   if(xl2.gt.xjtotp+xj2 .or. xl2.lt.abs(xjtotp-xj2)) goto 300
   kpmax=min(jplujp,nint(2*xj2))
   isrow = max(irow,icol)
   iscol = min(irow,icol)
   ii = (isrow * (isrow-1)) / 2 + iscol
   t  = -cmplx(sreal(ii),simag(ii))
   diagj = j1 .eq. j2
   diag = diagj .and. diagin
   if (diag .and. l1.eq.l2) t = t + (1.d0,0.d0)
   l2m=lmax-l2
! added 3rd subscript in next statement (pjd)
   icolp=iadr(j2,lmaxp-l2,5+ipack(icol))
   if(icolp.eq.0) goto 300
   potenz = ((2.d0* xj1) - xj2 + xl2)
   ipower = nint(xjtotp+potenz)
   ifaka=1
   if(mod(ipower,2).ne.0) ifaka=-1
   ipower = nint(xjtot+potenz)
   ifakb=1
   if(mod(ipower,2).ne.0) ifakb=-1
! change index for matel (pjd)
!         ic = matel(j2+1)
! replace previous statement with following code (pjd)
   do 104 ij = 1,njmax
   if (j2.ne.jslist(ij)) goto 104
   if (ipack(icol).ne.inlist(ij)) goto 104
   ic = ij
   goto 105
104    continue
! level not in list (pjd)
   ic = 100
105    continue
! end of added code (pjd)
   ll=0
do 200 l1p=l1pmin,l1pmax,2
   ll=ll+1
! added 3rd subscript in next statement (pjd)
   irowp=iadr(j1,lmaxp-l1p,5+ipack(irow))
! do not include terms if both levels not ln list (pjd)
!         if(irowp.eq.0) goto 200
   if(irowp.eq.0 .or. icolp.eq.0) goto 200
   isrowp = max(irowp,icolp)
   iscolp = min(irowp,icolp)
   iip = ioffs + (isrowp * (isrowp-1)) / 2 + iscolp
!
! convert s-matrix to t-matrix: t(j,j1,in1,l1 ; j2,in2,l2) =
!     delta(j1,in1,l1 ; j2,in2,l2) - s(j,j1,in1,l1 ; j2,in2,l2)
!
! replace 35525 below with global parameter lbufs (pjd)
   if (iip .gt. lbufs) write (6,11)iip
11    format (' *** IIP = ', i8)
   tp = -cmplx(srealp(iip),simagp(iip))
   if (diag .and. l1p.eq.l2) tp = tp + (1.d0,0d0)
   fac = real(t*conjg(tp))
   factra = fac*ifaka
   factrb = fac*ifakb
! loop over k
   kkp=0
!      print*,'jtot,jtotp,j1,l1,j2,l2,l1p,icolp,irowp',
!    1         jtot,jtotp,j1,l1,j2,l2,l1p,icolp,irowp
 do 160 kp=kpmin,kpmax
   kmin1=abs(kp-n)
   kmax1=min(kp+n,kmax)
   ikk=kplist(kp)
   kkp=kkp+1
   fac=f6a(kkp,j2,l2m)
   faca=factra*fac
   facb=factrb*fac
 do 150 k=kmin1,kmax1,n
   ikk=ikk+1
! switch row and column headings to match integral cross section output (pjd)
!         sigma(ikk,ir,ic)=sigma(ikk,ir,ic)
!     1       +f9a(ikk,ll)*(faca+f9pha(ikk)*facb)
   sigma(ikk,ic,ir)=sigma(ikk,ic,ir) &
       +f9a(ikk,ll)*(faca+f9pha(ikk)*facb)
150 continue
160 continue
200 continue
300 continue
400 continue
! next jtot'
500 jtotp = jtotp + jtpstp
if(jtotp.le.jtpmax) goto 60
! next jtot
700 jt = jt + 1
if (jt.le.maxjt) goto 10
! next parity
if (twopar) then
    twopar = .false.
    jt = jfsts
    jstart = jfsts
    jlp = 1
    jjoff = jlp * nwaves + 1
    ihibuf = 0
    igjtp = -1
    goto 10
end if
! here if calculation has been done
! write *.tcb as formattted file (pjd)
!999   write(4) nk
999 write(4, *) nk
ikk=0
do 1000 kp= 0, maxk
do 1000 k = abs(kp-n),min(maxk,kp+n),n
ikk=ikk+1
! write *.tcb as formattted file (pjd)
!      write(4) k,kp
write(4, *) k,kp
write(2,800) n,k,kp
if(.not.batch) write(6,800) n,k,kp
800 format(/' LAMBDA =',i2,'; TENSOR RANK KI =',i3,' KF =',i3/)
do 1000 i = 1, njmax
write(2,1010) (sigma(ikk,i,j),j=1,njmax)
if(.not. batch) write(6,1010) (sigma(ikk,i,j),j=1,njmax)
! write *.tcb as formattted file (pjd)
!1000  write(4) (sigma(ikk,i,j),j=1,njmax)
1000 write(4, *) (sigma(ikk,i,j),j=1,njmax)
1010 format(1x,10(1pd12.4))
! use diag as scratch variable (ipos = .false. for screen output)
diag = .false.
call mtime(cpu1,ela1)
cpu1 = cpu1 - cpu0
ela1 = ela1 - ela0
call gettim(ela1,elaps)
call gettim(cpu1,cpu)
!	      write(2,1200) n,elaps, cpu
!      if(.not. batch) write(6,1200) n,elaps, cpu
!1200  format(/' ** N =',i2,' COMPLETED, TIMING ELAPSED: ',a,
!     :        ' CPU: ',a)
!	      write(2,1200) n,elaps, cpu,t6j,t9j,
!     :      lenk,max1,max2,max3,maxkk,tsetup,tdelt,tchi,tlog
!      if(.not. batch) write(6,1200) n,elaps, cpu,t6j,t9j,
!     :      lenk,max1,max2,max3,maxkk,tsetup,tdelt,tchi,tlog
!1200  format(/' ** N =',i2,' COMPLETED, TIMING ELAPSED: ',a,
!     :        ' CPU: ',a/
!     :        ' 6J Symbols:',f10.2,'  9J Symbols:',f10.2/
!     : ' LK=',i4,'  L1=',i4,'  L2=',i4,'  L3=',i4,'  LKK',i4/
!     : '  SETUP:',f10.2,'  TDELT:',f10.2,'  TCHI:',f10.2,
!     : '  TLOG:',f10.2)
return
end
!-------------------------------------------------------------------------
subroutine sigkc(maxk,nnout,jfirst,jfinal,jtotd,nj,mmax,jpack, &
                lpack,ipack,jttble,prefac, &
                sreal,simag,matel,lenlab,labadr, &
                jtotpa,fast,ierr)
!
! subroutine to calculate tensor cross sections with the
! quantization axis along the initial relative velocity vector
! (derivation by pj dagdigian)
!
! set k = kp and q = qp = 0, to make connection with sigk values
!
! author: pj dagdigian
! current revision date: 5-mar-2010 by pj dagdigian
!------------------------------------------------------------------------
use tensor
use mod_cojq, only: jq ! jq(1)
use mod_colq, only: lq ! lq(1)
use mod_coinq, only: inq ! inq(1)
! common blocks for levels for which xs's to be computed
use mod_coisc9, only: jslist => isc9 ! jslist(1)
use mod_coisc10, only: inlist => isc10 ! inlist(1)
use mod_hibrid2, only: mxoutd
use mod_hibrid5, only: sread
use mod_par, only: batch, ipos
implicit double precision (a-h,o-z)
complex*8 t, tp
logical diag, diagj, diagin, &
        twopar, fast

!* flags for diagnostic printing
logical lprnt,lprntf

character*10 elaps, cpu
common /comom/  spin, xj1,xj2, j1, in1, j2, in2, maxjt, maxjot, &
            nwaves, jfsts, jlparf, jlpars, njmax, j1min, j2max

common /cospbf/ lnbufs, lnbufl, nbuf, maxlsp, maxllb, ihibuf, &
                igjtp
common /coipar/ ipar(9),iprnt
! add 3rd subscript for state index (subscript = 5 + IN)
! states with up to 9 state indices allowed
common/cadr/ iadr(0:jmx,lmx,9)
!      common/cadr/ iadr(0:jmx,lmx)
common/c6jt/ f6a(kmx,0:jmx,lmx),f6p(kmx)
!
dimension jpack(1),ipack(1),lpack(1)
dimension sreal(1), simag(1)
dimension prefac(1), matel(1), labadr(1), jtotpa(1)
dimension lenlab(1), jttble(1)
dimension sigma(nj,nj,maxk+1)
! array to store partial-j cross sections (K=0-2 only)
dimension psig0(1001),psig1(1001),psig2(1001)
!
! FLAG FOR DIAGNOSTIC PRINTING
lprnt = .false.
if (iprnt .eq. 2) lprnt = .true.

! FLAG FOR PARTIAL-Jtot CROSS SECTIONS
lprntf = .false.
if (iprnt .eq. 3) lprntf = .true.
!
! initialize timer
call mtime(cpu0,ela0)
ierr = 0
! clear sigma array
do 6 k = 0, maxk
do 5 j = 1, njmax
do 4 i = 1, njmax
sigma(i,j,k+1) = 0.d0
4 continue
5 continue
6 continue
if(maxk.gt.kmx.or.j2max.gt.jmx) then
  write(6,7) maxk,kmx,j2max,jmx
7   format(/' SIGKC: jmx or kmx too small:',4i5)
  ierr=-1
  return
end if
diagin = .false.
if (in1 .eq. in2) diagin = .true.
twopar = .false.
if (jlpars.ne.0) twopar = .true.
jt = jfirst
jstart = jfirst
jlp = 1
if (jlparf .eq. 1) jlp = 0
jjoff = jlp * nwaves + 1
ihibuf = 0
igjtp = -1
jtpstp = 1
!
! initialize arrays to store partial-jtot cross sections (K=0-2 only)
if (lprntf) then
  do 1006 jj = 1, 301
  psig0(jj) = 0.d0
  psig1(jj) = 0.d0
  psig2(jj) = 0.d0
1006   continue
end if
!
! sum over jtot and over parities (if twopar = .true.)
10 continue
!
! find address of jtot in smatrix file
jj = jjoff + jt
iaddr = jttble(jj)
if(iaddr .lt. 0) goto 700
! read s-matrix for jtot
nopen = -1
call sread ( iaddr, sreal, simag, jtot, jlpar, nu, &
            jq, lq, inq, ipack, jpack, lpack, &
             1, mmax, nopen, length, ierr)
if(ierr.eq.-1) goto 999
if(ierr.lt.-1) then
  write(2,20)
  if(.not.batch) write(6,20)
20   format(' *** READ ERROR, ABORT')
  return
end if
xjtot = jtot + spin
facj = (2.d0 * xjtot + 1.d0)
!
! boundaries for sum over jtot'
jtpmin = jtot
jtpmax = min((jtot + maxk), jfinal)
! set parity to be considered for jtot' (pjd)
jlparp = jlpar
!
! fill buffer with required s' matrices
! parity for each jtot' needs to be kept track of in addsp
call addsp(jtpmin,jtpmax,jlp, &
           labadr,lenlab,jtotpa,jttble)
if (srealp(lbufs) .ne. 0.d0) print *, 'srealp error in sigkc'
if (simagp(lbufs) .ne. 0.d0) print *, 'simagp error in sigkc'
if (ipackp(lbuflb) .ne. 0) print *, 'ipackp error in sigkc'
if (lpackp(lbuflb) .ne. 0) print *, 'lpackp error in sigkc'
if (jpackp(lbuflb) .ne. 0) print *, 'jpackp error in sigkc'
! prepare for sum over jtot'
jtotp = jtpmin
!
! start of loop for sum over jtot'
60 continue
!
! set parity to be considered for jtotp:
!   jlparp = jlpar for delta(jtot - jtotp) = even
!   jlparp = -jlpar for delta () = odd
jlparp = jlpar
if (((jtotp-jtot)/2)*2 .ne. (jtotp-jtot)) jlparp = -jlpar
! switch=1 to include (jtot,jtotp) and (jtotp,jtot) terms for jtot .ne. jtotp
switch = 1.d0
if (jtotp .eq. jtot .or. jtotp .gt. maxjt) switch = 0.
! find address of jtot' in s-matrix buffer
jj = jtotp + 1
! alternate parity for delta(jtotp-jtot) even/odd
   if (((jtotp-jtot)/2)*2 .ne. (jtotp-jtot)) then
! delta=odd:  switch parity of jtotp
! first check whether jlpar=1 or -1
      if (jjoff .eq. 1) then
! jlpar=1, skip forward to index for jtotp, jlpar=-1
         jjp = jjp + nwaves
      else
! jlpar=-1:  skip backward to index for jtotp, jlpar=1
         jjp = jjp - nwaves
      end if
   end if
!
ioffs = jtotpa(jj)
ioff  = labadr(jj)
lengtp = lenlab(jj)
xjtotp = jtotp + spin
facjjp = (2.d0 * xjtotp + 1.d0) * facj
jminjp = iabs(jtot-jtotp)
jplujp = nint(xjtot + xjtotp)
kmin=jminjp
jpmax = 0
do irowp = 1,lengtp
  jpmax = max(jpackp(ioff+irowp),jpmax)
end do
! now loop over all transitions
lmax=jtotp+jpmax+1
! clear iadr array
do 69 j=0,jmx
do 68 l=1,lmx
do 67 index=1,9
iadr(j,l,index) = 0
67 continue
68 continue
69 continue
!
do 70 irowp = 1, lengtp
  irp = ioff + irowp
  j1p = jpackp(irp)
  l1p = lpackp(irp)
  indp = ipackp(irp)
  if (ipackp(irp).ne.in1 .and. ipackp(irp).ne.in2) goto 70
  iadr(j1p,lmax-l1p,5+indp) = irowp
70 continue
!
jmx1=0
!
! sum over row index for jtot
do 400 irow = 1, length
  j1 = jpack(irow)
  xj1 = j1 + spin
  if (j1 .gt. j2max ) goto 400
  if (j1 .gt. jpmax ) goto 400
  if (j1 .lt. j1min ) goto 400
  l1 = lpack(irow)
  xl1 = l1
! delete next statement
!        if (xl1.lt.abs(xjtotp-xj1) .or. xl1.gt.xjtotp+xj1) goto 400

  j1t2 = nint(2.d0*xj1)
  do 74 ij = 1, njmax
  if (j1.ne.jslist(ij)) goto 74
  if (ipack(irow).ne.inlist(ij)) goto 74
  ir = ij
  denrow = prefac(ij)
  goto 75
74   continue
! level not in list
  ir = 100
75   continue
  jmx1=max(j1,jmx1)
  xjmx1=jmx1+spin
  kmax=min(jplujp,int(2*xjmx1))
! kmax must be <= maxk of tensor xs's to be calculated (pjd)
  kmax=min(kmax,maxk)
  kk=0
  l1m=lmax-l1
!
! sum over column index for jtot - scan over full set of (irow,icol) values
   do 200 icol = 1,length
    if (ipack(irow).eq.in1 .and. ipack(icol).eq.in2) &
       go to 82
    if (ipack(irow).eq.in2 .and. ipack(icol).eq.in1) &
       go to 82
    go to 200
82     continue
    j2 = jpack(icol)
    xj2 = j2 + spin
    if (j2 .gt. j2max) goto 200
    if (j2 .gt. jpmax) goto 200
    if (j2 .lt. j1min) goto 200
    l2 = lpack(icol)
    xl2 = l2
    if (xl2.lt.abs(xjtotp-xj2) .or. xl2.gt.xjtotp+xj2) goto 200
    diagj = j1.eq.j2
    j2t2 = nint(2.d0*xj2)
    min2j = min(j1t2,j2t2)
    if (min2j .lt. jminjp) goto 200
    do 84 ij = 1, njmax
      if (j2.ne.jslist(ij)) goto 84
      if (ipack(icol).ne.inlist(ij)) goto 84
      ic = ij
      dencol = prefac(ij)
      goto 85
84     continue
! level not in list
    ic = 100
85     continue
! determine index of element in s matrix
    if(irow .ge. icol) then
      ii = (irow*(irow-1))/2 + icol
    else
      ii = (icol*(icol-1))/2 + irow
    end if
    l2m=lmax-l2
!
! sum over row index for jtotp
    do 150 irowp = 1, lengtp
      irpp = ioff + irowp
      j1pp = jpackp(irpp)
      l1pp = lpackp(irpp)
! basis function must be the same as for irow
      if (j1pp.ne.j1) goto 150
      if (ipackp(irpp) .ne. ipack(irow)) goto 150
      xl1pp = l1pp
! icolp is index for (j2,l2,ind.jtotp) basis function
      icolp =iadr(j2,l2m,5+ipack(icol))
      if (icolp .eq. 0) goto 150
! determine index of element in s' matrix
      if(irowp.ge.icolp) then
        iip = ioffs + (irowp*(irowp-1))/2 + icolp
      else
        iip = ioffs + (icolp*(icolp-1))/2 + irowp
      end if
      diag = diagj .and. diagin
      if (diag) dencol = 0.
!* phase factor - the same for (jtot,jtotp) and (jtot,jtotp) terms
      phasea = 1.d0
      ipower= nint(xj1 + xj2 - xjtot - xjtotp + xl2)
      if ((ipower/2)*2 .ne. ipower) phasea = -1.d0
!
! convert s-matrix to t-matrix: t(j,j1,in1,l1 ; j2,in2,l2) =
!     delta(j1,in1,l1 ; j2,in2,l2) - s(j,j1,in1,l1 ; j2,in2,l2)
      t  = -cmplx(sreal(ii),simag(ii))
      tp = -cmplx(srealp(iip),simagp(iip))
      if (iip .gt. lbufs) write (6,11)iip
11       format (' *** IIP = ', i8, ' .ge. LBUFS')
      if (diag .and. l1.eq.l2) then
        t = t + (1.d0,0)
      end if
      if (diag .and. l1pp.eq.l2) then
        tp = tp + (1.d0,0)
      end if

!** DIAGNOSTIC PRINT
!          write (6,*) t,tp

! iph below is i ** (l1 - l1pp) factor
      iph = nint((xl1pp - xl1)/2.d0)
      phasel = 1.d0
      if ((iph/2)*2 .ne. iph) phasel = -1.d0
      factl = sqrt((2.d0 * xl1 + 1.d0) &
        * (2.d0 * xl1pp + 1.d0))
      factor = facjjp * phasea * phasel &
        * factl * real(t*conjg(tp))
!
! loop over k
      kmax = min(min2j,jplujp)
! kmax must be <= maxk of tensor xs's to be calculated
      kmax = min(kmax,maxk)
      do 100 k = kmin, kmax
        xk = k
        factk = (2.d0 * xk + 1.d0)
! next statements are to include (jtot, jtotp) term
        fmsum = 0.d0
        mterms = nint(2.d0 * xj1 + 1.d0)
        xmmin = -xj1
        do 95 m = 1,mterms
          xm = xmmin + (m - 1)
          fmsum = fmsum + xf3j(xj1, xj1, xk, -xm, xm, 0.d0) &
            * xf3j(xj1, xl1, xjtot, xm, 0.d0, -xm) &
            * xf3j(xj1, xl1pp, xjtotp, xm, 0.d0, -xm) &
            * xf3j(xjtotp, xjtot, xk, xm, -xm, 0.d0)
95         continue
        factr = factk * fmsum * &
          xf6j(xj2,xj2,xk,xjtot,xjtotp,xl2) * factor

!  below for diagnostic print
       xx6j = xf6j(xj2,xj2,xk,xjtot,xjtotp,xl2)

! next statements are to include (jtotp, jtot) term for jtotp .ne. jtot
! the contribution of the (jtot,jtotp,l1,l1pp) and (jtotp,jtot,l1pp,l1)
! terms are equal
        if (switch .eq. 1.d0) then
          factr = 2.d0 * factr
        end if
        sigma(ir,ic,k+1) = sigma(ir,ic,k+1) + factr*dencol
        sigma(ic,ir,k+1) = sigma(ic,ir,k+1) + factr*denrow
!
! accummulate partial cross sections
        if (lprntf) then
! print values for only one transition
          istate = 1
          if (ir.eq.istate .and. ic.eq.istate) then
            if (k .eq. 0) psig0(jtot+1) = psig0(jtot+1) &
              + factr*denrow
            if (k .eq. 1) psig1(jtot+1) = psig1(jtot+1) &
              + factr*denrow
            if (k .eq. 2) psig2(jtot+1) = psig2(jtot+1) &
              + factr*denrow
          end if
        end if
100       continue
!
150     continue
200   continue
400 continue
!
! next jtot'
500 jtotp = jtotp + jtpstp
if(jtotp.le.jtpmax) goto 60
!
! next jtot
700 jt = jt + 1
if (jt.le.maxjt) goto 10
!
! next parity
if (twopar) then
  twopar = .false.
  jt = jfsts
  jstart = jfsts
  jlp = 1
  jjoff = jlp * nwaves + 1
  ihibuf = 0
  igjtp = -1
! loop back for next parity
  goto 10
end if
! here if calculation has been done
!
! write *.tcb as formatted file
999 write(4,798)
798 format('cross sections in the collision frame')
write(4, *) maxk + 1
write (2,799)
if(.not.batch) write(6,799)
799 format(/' CROSS SECTIONS IN THE COLLISION FRAME')
do 1100 k = 0, maxk
  write(4, *) k,k
  write(2,800) k
  if(.not.batch) write(6,800) k
800   format(/' TENSOR RANK K =',i3)
  do 1000 i = 1, njmax
! write *.tcb as formattted file (pjd)
    write(4, *) (sigma(i,j,k+1),j=1,njmax)
1000   continue
  call mxoutd (2, sigma(1,1,k+1), njmax, njmax, 0, ipos)
! use diag as scratch variable (ipos = .false. for screen output)
  diag = .false.
  if(.not. batch) &
    call mxoutd (6, sigma(1,1,k+1), njmax, njmax, 0, diag)
1100 continue

! PRINT OUT PARTIAL CROSS SECTIONS
if (lprntf) then
  write (6,1011) istate,istate
1011   format (/' PARTIAL-JTOT CROSS SECTIONS for transition: ', &
     '  level = ',i2,' -> ',' level = ',i2/ &
     '   jtot',7x,'P-XS(0)',8x,'P-XS(1)',8x,'P-XS(2)')
  do jj = 0,maxjt
    write (6,1012) jj,psig0(jj+1),psig1(jj+1), &
       psig2(jj+1)
1012     format(i6,2x,3e15.6)
  end do
end if

call mtime(cpu1,ela1)
cpu1 = cpu1 - cpu0
ela1 = ela1 - ela0
call gettim(ela1,elaps)
call gettim(cpu1,cpu)
write(2,1200) elaps, cpu
if(.not. batch) write(6,1200) elaps, cpu
1200 format(/' ** N = 0 COMPLETED, TIMING ELAPSED: ',a, &
        ' CPU: ',a,/)
return
end
!-------------------------------------------------------------------------
subroutine dsigh(maxk,nnout,jfirst,jfinal,jtotd,jpack, &
                lpack,ipack,jlevel,inlevel,elevel,flaghf, &
                iframe,ierr)
!
! subroutine to calculate m-resolved differential cross sections
! for the elastic (j1,in1) -> (j1,in1) transition in the
! helicity frame.  these are then used to compute the
! corresponding tensor cross sections
!
! author: pj dagdigian
! current revision date: 7-oct-2011 by pj dagdigian
!------------------------------------------------------------------------
use mod_codim, only: mmax
use mod_cov2, only: nv2max, ndummy, y => v2
use mod_cojq, only: jq ! jq(1)
use mod_colq, only: lq ! lq(1)
use mod_coinq, only: inq ! inq(1)
use mod_cojhld, only: jlev => jhold ! jlev(1)
use mod_coisc1, only: inlev => isc1 ! inlev(1)
! common blocks for levels for which xs's to be computed
use mod_coisc9, only: jslist => isc9 ! jslist(1)
use mod_coisc10, only: inlist => isc10 ! inlist(1)
use mod_cosc1, only: elev => sc1 ! elev(1)
use mod_coz, only: sreal => z_as_vec ! sreal(1)
use mod_cow, only: simag => w_as_vec ! simag(1)
use mod_hibrid5, only: sread
use mod_difcrs, only: sphn
use constants, only: econv, xmconv, ang2c
use mod_par, only: batch, ipos
implicit double precision (a-h,o-z)
! size of q for j <= 5 and 0.5 deg angle increment
complex*16 q(43681)
logical diag, diagj, diagin, &
        twopar, fast, flaghf
logical existf,csflg1,flghf1,flgsu1,twomol, &
        nucros, iprint
character*10 elaps, cpu
character*20 cdate1
common /comom/  spin, xj1,xj2, j1, in1, j2, in2, maxjt, maxjot, &
            nwaves, jfsts, jlparf, jlpars, njmax, j1min, j2max
common /cospbf/ lnbufs, lnbufl, nbuf, maxlsp, maxllb, ihibuf, &
                igjtp
common /coipar/ ipar(9),iprnt
common /coisc2 / jout1(1)
common /coang/ ang1, ang2, dang
!
dimension jpack(1),ipack(1),lpack(1)
! size of arrays in 2 dimension statements below are for j.le.5
dimension sigkto(21),sigk(21)
dimension sigmmp(361,11,11),mivals(11,11), &
  fmivals(11,11),mfvals(11,11),fmfvals(11,11)
data pi/3.141592653589793d0/
data rad/57.295779513082323d0/
! compute dif cross sections over 0 to 180 degrees, in 0.5 deg increments
data ang1, ang2, dang /0.d0, 180.d0, 0.5d0/
!
maxq = mmax*mmax/2
maxy = nv2max
!
! print out m-resolved differential cross sections? (immprt = 1)
immprt = 0
if (iframe .le. 0) immprt = 1
! flag for printing in amplih subr
ihfst = 0
!
if (jlevel .gt. 5) then
  write(6,12) jlevel
12   format(/' *** JLEVEL =',i3,' .GT. 5. TOO HIGH FOR DSIGH. ABORT')
  return
end if
l2max = jfinal + jlevel + 1
mlmax = jlevel + jlevel + 1
! below for half-integer spin
if (flaghf) then
  l2max  = l2max + 1
  mlmax = mlmax + 1
end if
nangle = nint((ang2 - ang1)/dang) + 1
mx = maxy/(l2max*mlmax)
jlevlp = jlevel
! below for half-integer spin
if (flaghf) jlevlp = jlevel + 1
! zero out amplitudes (complete array)
do 210 ii = 1, 43681
  q(ii) = cmplx(0.d0, 0.d0)
210 continue
!
! precalculate all required spherical harmonics
ii = 0
do 230 ml = 0, mlmax-1
  angle = ang1
  do 220 i = 1, nangle
    call sphn(ml, l2max-1, angle, y(ii+i), nangle)
    angle = angle + dang
220   continue
  ii = ii + l2max*nangle
230 continue
jtlast = -1
jplast = 0
!
!.....read header of s-matrix file
!
call rdhead(1,cdate1,ered1,rmu1,csflg1,flghf1,flgsu1, &
   twomol,nucros,jfirst,jfinal,jtotd,numin,numax,nud, &
   nlevel,nlevop,nnout,jlev,inlev,elev,jout1)
!
!.....read next s-matrix
!
250 nopen = 0
call sread (0, sreal, simag, jtot, jlpar, nu, &
            jq, lq, inq, ipack, jpack, lpack, &
            1, mmax, nopen, length, ierr)
if(ierr .eq. -1) then
   write(6,260) xnam1,jtlast,jplast
260    format(' END OF FILE DETECTED READING FILE ',(a), &
     ' LAST JTOT,JLPAR PROCESSED:',2i5)
   goto 310
end if
if(ierr .lt. -1) then
  write(6,270) xnam1,jtlast,jplast
270   format(' ERROR READING FILE ',(a), &
     ' LAST JTOT,JLPAR PROCESSED:',2i5)
  goto 310
end if
!
!.....this assumes that jlpar=1 is stored first
!
if(jtot .gt. jfinal) goto 300
!
!.....copy row labels into column labels if s-matrices are stored
!.....triangular
!
if(jlpar.eq.jplast .and. jtot.ne.jtlast+1) write(6,275) jtot,jtlast
275 format(' *** WARNING: JTOT.NE.JTLAST+1:',2i4)
jtlast=jtot
jplast=jlpar
if(nnout.gt.0) then
  do 290 i=1,length
    inq(i)=ipack(i)
    jq(i)=jpack(i)
    lq(i)=lpack(i)
290   continue
  nopen = length
end if
!
!.....calculate contributions to amplitudes for present jtot
!     for elastic (jlevel,inlevel) -> (jlevel,inlevel) transition
!
call amplih(jlevel,inlevel,jlevel,inlevel,jtot,mmax, &
  jpack,lpack,ipack,length,jq,lq,inq,nopen, &
  l2max,ihfst,nangle,flaghf,sreal,simag,y,q)
!
!.....loop back to next jtot/jlpar
!
300 if(jtot.lt.jfinal .or. jlpar.eq.1) goto 250
!
!.....ca is wavevector for initial state, ecol is collision energy
ecol = ered1 - elevel
ca=sqrt(2.d0*rmu1*ecol)
fak=ang2c/ca**2
!
!.....print m -> m' differential cross sections for this batch of angles
!
numjm = jlevlp + jlevel + 1
!
!     SUPPRESS PRINTING OF M-RESOLVED CROSS SECTIONS IN PRODUCTION RUNS (immprt = 0)
if (immprt .ne. 0) then
write(2,302)
302 format(/'%   M-DEPENDENT HELICITY FRAME ELASTIC DIFFERENTIAL', &
  ' CROSS SECTIONS (ANG^2/STR)'/)
write(6,304)
304 format(/'    M-DEPENDENT HELICITY FRAME ELASTIC DIFFERENTIAL', &
  ' CROSS SECTIONS (ANG^2/STR)'/)
do 308 mj1 = -jlevlp, jlevel
do 308 mj2 = -jlevlp, jlevel
  isub1=mj1+jlevlp+1
  isub2=mj2+jlevlp+1
  if (flaghf) then
     xmj1=dble(mj1)+spin
     xmj2=dble(mj2)+spin
     fmivals(isub1,isub2) = xmj1
     fmfvals(isub1,isub2) = xmj2
  else
     mivals(isub1,isub2) = mj1
     mfvals(isub1,isub2) = mj2
  end if
308 continue
!
if (flaghf) then
  write(2,312) ((fmivals(i1,i2), &
    fmfvals(i1,i2),i2 = 1,numjm), &
    i1 = 1,numjm)
312   format('%   ANGLE',20x,'(M -> M'')'/ &
    '%',10x,121(1x,f4.1,'->',f4.1,2x))
  write(6,313) ((fmivals(i1,i2), &
    fmfvals(i1,i2),i2 = 1,numjm), &
    i1 = 1,numjm)
313   format('    ANGLE',20x,'(M -> M'')'/ &
    10x,121(1x,f4.1,'->',f4.1,2x))
else
  write(2,314) ((mivals(i1,i2), &
    mfvals(i1,i2),i2 = 1,numjm), &
    i1 = 1,numjm)
314   format('%   ANGLE',20x,'(M -> M'')'/ &
    '%',8x,121(4x,i2,'->',i2,4x))
  write(6,315) ((mivals(i1,i2), &
    mfvals(i1,i2),i2 = 1,numjm), &
    i1 = 1,numjm)
315   format('    ANGLE',20x,'(M -> M'')'/ &
    9x,121(4x,i2,'->',i2,4x))
end if
write(2,303)
303 format (' sigmmp_hel=[')
end if
!     END OF SUPPRESSED PRINTING
!
angle = ang1
do 350 iang=1,nangle
  do 320 mj1 = -jlevlp, jlevel
  do 320 mj2 = -jlevlp, jlevel
    ii = (mj1 + jlevlp)*numjm*nangle &
         + (mj2 + jlevlp)*nangle + iang
    isub1=mj1+jlevlp+1
    isub2=mj2+jlevlp+1
    sigmmp(iang,isub1,isub2) = &
      fak*dreal(q(ii)*conjg(q(ii)))
320   continue
!     SUPPRESS THIS PRINTING
  if (immprt .ne. 0) then
  write(2,355) angle,((sigmmp(iang,i1,i2), &
    i2 = 1,numjm), i1 = 1,numjm)
  write(6,355) angle,((sigmmp(iang,i1,i2), &
    i2 = 1,numjm), i1 = 1,numjm)
355   format(f8.2,121e14.6)
  end if
  angle = angle + dang
350 continue
if (immprt .ne. 0) write(2,352)
352 format(' ];')
!     END OF SUPPRESSED PRINTING OF CROSS SECTION VALUES
!
! now compute differential tensor cross sections
!
write(2,362) (k,k=0,maxk)
362 format(/'%   HELICITY FRAME ELASTIC DIFFERENTIAL', &
  ' TENSOR CROSS SECTIONS (ANG^2/STR)'//'%   ANGLE', &
   12x,'TENSOR RANK K'/ &
  '%',6x,21(6x,i3,4x))
write(6,364) (k,k=0,maxk)
364 format(/'    HELICITY FRAME ELASTIC DIFFERENTIAL', &
  ' TENSOR CROSS SECTIONS (ANG^2/STR)'//'    ANGLE', &
   12x,'TENSOR RANK K'/ &
  ' ',6x,21(6x,i3,4x))
write(2,363)
363 format (' sigk_hel = [')
! zero integral tensor cross sections
do 360 k1=1,maxk
  sigkto(k1) = 0.d0
360 continue
angle = ang1
xj1=dble(jlevel)+spin
do 400 iang=1,nangle
  sn = sin(angle/rad)
  maxk1 = maxk+1
  do 390 k1=1,maxk1
    xk = k1 - 1
    fack = 2.d0*xk + 1.0d0
    sigk(k1)=0.d0
    do 380 mj1=-jlevlp, jlevel
      xmj1=dble(mj1)+spin
      do 375 mj2=-jlevlp, jlevel
        xmj2=dble(mj2)+spin
        ipower=nint(xj1+xj1-xmj1-xmj2)
        iph=1.d0
        if ((ipower/2)*2 .ne. ipower) iph=-1.d0
        isub1=mj1+jlevlp+1
        isub2=mj2+jlevlp+1
        sigk(k1) = sigk(k1) + iph*fack &
          *xf3j(xj1,xj1,xk,xmj1,-xmj1,0.d0) &
          *xf3j(xj1,xj1,xk,xmj2,-xmj2,0.d0) &
          *sigmmp(iang,isub1,isub2)
375       continue
380     continue
    sigkto(k1) = sigkto(k1) + sn*sigk(k1)
390   continue
  write(2,355) angle,(sigk(i),i=1,maxk1)
  write(6,355) angle,(sigk(i),i=1,maxk1)
  angle = angle + dang
400 continue
do 401 k1=1,maxk1
  sigkto(k1) = sigkto(k1)*(dang/rad)*2.d0*pi
401 continue
write(2,352)
write(2,405) (k, sigkto(k+1),k=0, maxk)
write(6,406) (k, sigkto(k+1),k=0, maxk)
405 format(/,'%  INTEGRAL TENSOR CROSS SECTIONS'/ &
  '%     RANK  XS:  ',11(i6,1pe14.5))
406 format(/,'   INTEGRAL TENSOR CROSS SECTIONS'/ &
  '      RANK  XS:',11(i6,1pe14.5))
!
310 continue
!
return
end
!===
subroutine amplih(j1,inlev1,j2,inlev2,jtot,mmax,jpack, &
  lpack,ipack,length,jq,lq,inq,nopen,maxl2,ihfst, &
  nangle,flaghf,sreal,simag,y,q)
!
! calculates scattering amplitudes for given jtot and set
! of angles
!
!     author of original ampli program:  h.-j. werner
!     revised by p.j.dagdigian for helicity frame calculation
!     for the transition (j1,inlev1) -> (j2,inlev2)
!
!     revision date: 13-oct-2011
!
!.....jpack,lpack,ipack: labels for rows
!.....jq,lq,inq:         labels for columns
implicit double precision (a-h,o-z)
complex*16 q,ai,fak2,fak3,yy,tmat
parameter (zero=0.0d0,one=1.0d0)
logical flaghf,elastc
common /coang/ ang1, ang2, dang
dimension jpack(1),lpack(1),ipack(1),jq(1),lq(1),inq(1),q(1)
dimension sreal(mmax,1),simag(mmax,1),y(1)
dimension fak1(400),fak2(400),fak3(400),ilab1(400)
sqpi=1.772453850905516d0
rad=57.295779513082323d0
!
!.....ai is sqrt(-1)
ai=cmplx(zero, one)
elastc = j1.eq.j2 .and. inlev1.eq.inlev2
!
if (flaghf) then
!.....here for half-integer spin
  fakj=sqpi*(2.d0*jtot + 2.d0)*(-1)**(j1+j2+1)
  spin=0.5d0
  j1p=j1+1
  j2p=j2+1
else
!.....here for integer spin
  fakj=sqpi*(2.d0*jtot + 1.d0)*(-1)**(j1+j2)
  spin=0.0d0
  j1p=j1
  j2p=j2
end if
xjtot=dble(jtot)+spin
xj1=dble(j1)+spin
xj2=dble(j2)+spin
ll=0
do 50 ilab=1,length
  if(jpack(ilab).ne.j1 .or. ipack(ilab).ne.inlev1) &
     goto 50
  l1=lpack(ilab)
  ll=ll+1
  fak1(ll)=fakj*sqrt(2.d0*l1 + 1.0d0)
  ilab1(ll)=ilab
50 continue
llmax=ll
!
do 500 jlab=1,nopen
  if(jq(jlab).ne.j2 .or. inq(jlab).ne.inlev2) goto 500
  l2=lq(jlab)
  xl2=l2
  do 60 ll=1,llmax
    ilab=ilab1(ll)
    l1=lpack(ilab)
!.....convert to t-matrix
    tmat=-cmplx(sreal(jlab,ilab),simag(jlab,ilab))
    if(elastc .and. l1.eq.l2) tmat = tmat + 1.0d0
    fak2(ll)=cmplx(fak1(ll)*(-1)**(l1+l2),zero) &
       *(ai**(l1-l2))*tmat
60   continue
!
  ii=0
  do 400 mj1=-j1p,j1
    xmj1=dble(mj1)+spin
!
    do 70 ll=1,llmax
      ilab=ilab1(ll)
      xl1=lpack(ilab)
      fak3(ll)=fak2(ll)*xf3j(xj1,xl1,xjtot,xmj1,zero,-xmj1)
70     continue

! set ihel = 1 for helicity-frame calculations, ihel = 0 for collision-frame
    ihel = 1

!helicity frame calculation
    if (ihel .eq. 1) then
! mj2 is helicity-frame final m quantum number
    do 300 mj2=-j2p,j2
      xmj2=dble(mj2)+spin
      angle = ang1
      do 200 iang=1,nangle
        yy=0.d0
! mj2p is collision-frame final m quantum number
        do 350 mj2p=-j2p,j2
          xmj2p=dble(mj2p)+spin
          xml2p=xmj1 - xmj2p
          ml2p=xml2p
          iyof=(iabs(ml2p)*maxl2+l2)*nangle
! redrot requires rotation angle in radians
          beta=angle/rad
          fakp = redrot(xj2,xmj2p,xmj2,beta) &
            *xf3j(xj2,xl2,xjtot,xmj2p,xml2p,-xmj1) &
            *y(iyof+iang)
          if (ml2p.gt.0) fakp = fakp*(-1)**ml2p
          yy = yy + fakp
350         continue
        ii = ii + 1
        do 100 ll=1,llmax
          q(ii) = q(ii) + fak3(ll)*yy
100         continue
        angle = angle + dang
200       continue
300     continue

    else
!collision-frame calculation

    if (ihfst.eq.0) write(6,336) ihel
336     format(/'** ihel =',i2,' - THIS IS A COLLISION-FRAME', &
      ' CALCULATION')
    ihfst = 1

    do 1300 mj2=-j2p,j2
      xmj2=dble(mj2)+spin
      ml2=mj1-mj2
      iyof=(iabs(ml2)*maxl2+l2)*nangle
      if(ml2.gt.0) fak=fak*(-1)**ml2
      xml2=ml2
      fak=xf3j(xj2,xl2,xjtot,xmj2,xml2,-xmj1)
      if(ml2.gt.0) fak=fak*(-1)**ml2
!
      angle = ang1
      do 1200 iang=1,nangle
        yy=cmplx(fak*y(iyof+iang), 0.d0)
        ii = ii + 1
        do 1100 ll=1,llmax
          q(ii)=q(ii)+fak3(ll)*yy
1100         continue
        angle = angle + dang
1200       continue
1300     continue
    end if

400   continue
500 continue
return
end
!=============
subroutine dsigga(maxk,nnout,jfirst,jfinal,jtotd,jpack, &
                lpack,ipack,jlevel,inlevel,elevel,flaghf, &
                iframe,ierr)
!
! subroutine to calculate m-resolved differential cross sections
! for the elastic (j1,in1) -> (j1,in1) transition in the
! geometric apse frame.  these are then used to compute the
! corresponding tensor cross sections
!
! author: pj dagdigian
! current revision date: 7-oct-2011 by pj dagdigian
!------------------------------------------------------------------------
use mod_codim, only: mmax
use mod_cov2, only: nv2max, ndummy, y => v2
use mod_cojq, only: jq ! jq(1)
use mod_colq, only: lq ! lq(1)
use mod_coinq, only: inq ! inq(1)
use mod_cojhld, only: jlev => jhold ! jlev(1)
use mod_coisc1, only: inlev => isc1 ! inlev(1)
! common blocks for levels for which xs's to be computed
use mod_coisc9, only: jslist => isc9 ! jslist(1)
use mod_coisc10, only: inlist => isc10 ! inlist(1)
use mod_cosc1, only: elev => sc1 ! elev(1)
use mod_coz, only: sreal => z_as_vec ! sreal(1)
use mod_cow, only: simag => w_as_vec ! simag(1)
use mod_hibrid5, only: sread
use mod_difcrs, only: sphn
use constants, only: econv, xmconv, ang2c
use mod_par, only: batch, ipos

implicit double precision (a-h,o-z)
! size of q for j <= 5 and 0.5 deg angle increment
complex*16 q(43681)
logical diag, diagj, diagin, &
        twopar, fast, flaghf
logical existf,csflg1,flghf1,flgsu1,twomol, &
        nucros, iprint
character*10 elaps, cpu
character*20 cdate1
common /comom/  spin, xj1,xj2, j1, in1, j2, in2, maxjt, maxjot, &
            nwaves, jfsts, jlparf, jlpars, njmax, j1min, j2max
common /cospbf/ lnbufs, lnbufl, nbuf, maxlsp, maxllb, ihibuf, &
                igjtp
common /coipar/ ipar(9),iprnt
common /coisc2 / jout1(1)
common /coang/ ang1, ang2, dang
!
dimension jpack(1),ipack(1),lpack(1)
! size of arrays in 2 dimension statements below are for j.le.5
dimension sigkto(21),sigk(21)
dimension sigmmp(361,11,11),mivals(11,11), &
  fmivals(11,11),mfvals(11,11),fmfvals(11,11)
data pi/3.141592653589793d0/
data rad/57.295779513082323d0/
! compute dif cross sections over 0 to 180 degrees, in 0.5 deg increments
data ang1, ang2, dang /0.d0, 180.d0, 0.5d0/
!
maxq = mmax*mmax/2
maxy = nv2max
!
! print out m-resolved differential cross sections (iframe < 0)? (immprt = 1)
immprt = 0
if (iframe .le. 0) immprt = 1
!
if (jlevel .gt. 5) then
  write(6,12) jlevel
12   format (' *** JLEVEL =',i3,' .GT. 5. TOO HIGH FOR DSIGGA.', &
    ' ABORT')
  return
end if
l2max = jfinal + jlevel + 1
mlmax = jlevel + jlevel + 1
! below for half-integer spin
if (flaghf) then
  l2max  = l2max + 1
  mlmax = mlmax + 1
end if
nangle = nint((ang2 - ang1)/dang) + 1
mx = maxy/(l2max*mlmax)
jlevlp = jlevel
! below for half-integer spin
if (flaghf) jlevlp = jlevel + 1
! zero out amplitudes (complete array)
do 210 ii = 1, 43681
  q(ii) = cmplx(0.d0, 0.d0)
210 continue
!
! precalculate all required spherical harmonics
ii = 0
do 230 ml = 0, mlmax-1
  angle = ang1
  do 220 i = 1, nangle
    call sphn(ml, l2max-1, angle, y(ii+i), nangle)
    angle = angle + dang
220   continue
  ii = ii + l2max*nangle
230 continue
jtlast=-1
jplast=0
!
!.....read header of s-matrix file
!
call rdhead(1,cdate1,ered1,rmu1,csflg1,flghf1,flgsu1, &
   twomol,nucros,jfirst,jfinal,jtotd,numin,numax,nud, &
   nlevel,nlevop,nnout,jlev,inlev,elev,jout1)
!
!.....read next s-matrix
!
250 nopen = 0
call sread (0, sreal, simag, jtot, jlpar, nu, &
                  jq, lq, inq, ipack, jpack, lpack, &
                  1, mmax, nopen, length, ierr)
if(ierr .eq. -1) then
   write(6,260) xnam1,jtlast,jplast
260    format(' END OF FILE DETECTED READING FILE ',(a), &
     ' LAST JTOT,JLPAR PROCESSED:',2i5)
   goto 310
end if
if(ierr .lt. -1) then
  write(6,270) xnam1,jtlast,jplast
270   format(' ERROR READING FILE ',(a), &
     ' LAST JTOT,JLPAR PROCESSED:',2i5)
  goto 310
end if
!
!.....this assumes that jlpar=1 is stored first
!
if(jtot .gt. jfinal) goto 300
!
!.....copy row labels into column labels if s-matrices are stored
!.....triangular
!
if(jlpar.eq.jplast .and. jtot.ne.jtlast+1) write(6,275) jtot,jtlast
275 format(' *** WARNING: JTOT.NE.JTLAST+1:',2i4)
jtlast=jtot
jplast=jlpar
if(nnout.gt.0) then
  do 290 i=1,length
    inq(i)=ipack(i)
    jq(i)=jpack(i)
    lq(i)=lpack(i)
290   continue
  nopen = length
end if
!
!.....calculate contributions to amplitudes for present jtot
!     for elastic (jlevel,inlevel) -> (jlevel,inlevel) transition
!
call amplga(jlevel,inlevel,jlevel,inlevel,jtot,mmax, &
  jpack,lpack,ipack,length,jq,lq,inq,nopen, &
  l2max,nangle,flaghf,sreal,simag,y,q)
!
!.....loop back to next jtot/jlpar
!
300 if(jtot.lt.jfinal .or. jlpar.eq.1) goto 250
!
!.....ca is wavevector for initial state, ecol is collision energy
ecol = ered1 - elevel
ca=sqrt(2.d0*rmu1*ecol)
fak=ang2c/ca**2
!
!.....print m -> m' differential cross sections for this batch of angles
!
numjm = jlevlp + jlevel + 1
!
!     SUPPRESS PRINTING OF M-RESOLVED CROSS SECTIONS IN PRODUCTION RUNS (immprt = 0)
if (immprt .ne. 0) then
write(2,302)
302 format(/'%   M-DEPENDENT GEOMETRIC APSE FRAME ELASTIC', &
  ' DIFFERENTIAL CROSS SECTIONS (ANG^2/STR)'/)
write(6,304)
304 format(/' M-DEPENDENT GEOMETRIC APSE FRAME ELASTIC', &
  ' DIFFERENTIAL CROSS SECTIONS (ANG^2/STR)'/)
do 308 mj1 = -jlevlp, jlevel
do 308 mj2 = -jlevlp, jlevel
  isub1=mj1+jlevlp+1
  isub2=mj2+jlevlp+1
  if (flaghf) then
     xmj1=dble(mj1)+spin
     xmj2=dble(mj2)+spin
     fmivals(isub1,isub2) = xmj1
     fmfvals(isub1,isub2) = xmj2
  else
     mivals(isub1,isub2) = mj1
     mfvals(isub1,isub2) = mj2
  end if
308 continue
!
if (flaghf) then
  write(2,312) ((fmivals(i1,i2), &
    fmfvals(i1,i2),i2 = 1,numjm), &
    i1 = 1,numjm)
312   format('%   ANGLE',20x,'(M -> M'')'/ &
    '%',10x,121(1x,f4.1,'->',f4.1,2x))
  write(6,313) ((fmivals(i1,i2), &
    fmfvals(i1,i2),i2 = 1,numjm), &
    i1 = 1,numjm)
313   format('    ANGLE',20x,'(M -> M'')'/ &
    10x,121(1x,f4.1,'->',f4.1,2x))
else
  write(2,314) ((mivals(i1,i2), &
    mfvals(i1,i2),i2 = 1,numjm), &
    i1 = 1,numjm)
314   format('%   ANGLE',20x,'(M -> M'')'/ &
    '%',8x,121(4x,i2,'->',i2,4x))
  write(6,315) ((mivals(i1,i2), &
    mfvals(i1,i2),i2 = 1,numjm), &
    i1 = 1,numjm)
315   format('    ANGLE',20x,'(M -> M'')'/ &
    9x,121(4x,i2,'->',i2,4x))
end if
write(2,303)
303 format (' sigmmp_hel=[')
end if
!     END OF SUPPRESSED PRINTING
!
angle = ang1
do 350 iang=1,nangle
  do 320 mj1 = -jlevlp, jlevel
  do 320 mj2 = -jlevlp, jlevel
    ii = (mj1 + jlevlp)*numjm*nangle &
         + (mj2 + jlevlp)*nangle + iang
    isub1=mj1+jlevlp+1
    isub2=mj2+jlevlp+1
    sigmmp(iang,isub1,isub2) = &
      fak*dreal(q(ii)*conjg(q(ii)))
320   continue
!     SUPPRESS THIS PRINTING
  if (immprt .ne. 0) then
  write(2,355) angle,((sigmmp(iang,i1,i2), &
    i2 = 1,numjm), i1 = 1,numjm)
  write(6,355) angle,((sigmmp(iang,i1,i2), &
    i2 = 1,numjm), i1 = 1,numjm)
355   format(f8.2,121e14.6)
  end if
  angle = angle + dang
350 continue
if (immprt .ne. 0) write(2,352)
352 format(' ];')
!     END OF SUPPRESSED PRINTING OF CROSS SECTION VALUES
!
! now compute differential tensor cross sections
!
write(2,362) (k,k=0,maxk)
362 format(/'%  GEOMETRIC APSE FRAME ELASTIC DIFFERENTIAL', &
  ' TENSOR CROSS SECTIONS (ANG^2/STR)'//'%   ANGLE', &
   12x,'TENSOR RANK K'/ &
  '%',6x,21(6x,i3,4x))
write(6,364) (k,k=0,maxk)
364 format(/'    GEOMETRIC APSE FRAME ELASTIC DIFFERENTIAL', &
  ' TENSOR CROSS SECTIONS (ANG^2/STR)'//'    ANGLE', &
   12x,'TENSOR RANK K'/ &
  ' ',6x,21(6x,i3,4x))
write(2,363)
363 format (' sigk_hel = [')
! zero integral tensor cross sections
do 360 k1=1,maxk
  sigkto(k1) = 0.d0
360 continue
angle = ang1
xj1=dble(jlevel)+spin
do 400 iang=1,nangle
  sn = sin(angle/rad)
  maxk1 = maxk+1
  do 390 k1=1,maxk1
    xk = k1 - 1
    fack = 2.d0*xk + 1.0d0
    sigk(k1)=0.d0
    do 380 mj1=-jlevlp, jlevel
      xmj1=dble(mj1)+spin
      do 375 mj2=-jlevlp, jlevel
        xmj2=dble(mj2)+spin
        ipower=nint(xj1+xj1-xmj1-xmj2)
        iph=1.d0
        if ((ipower/2)*2 .ne. ipower) iph=-1.d0
        isub1=mj1+jlevlp+1
        isub2=mj2+jlevlp+1
        sigk(k1) = sigk(k1) + iph*fack &
          *xf3j(xj1,xj1,xk,xmj1,-xmj1,0.d0) &
          *xf3j(xj1,xj1,xk,xmj2,-xmj2,0.d0) &
          *sigmmp(iang,isub1,isub2)
375       continue
380     continue
    sigkto(k1) = sigkto(k1) + sn*sigk(k1)
390   continue
  write(2,355) angle,(sigk(i),i=1,maxk1)
  write(6,355) angle,(sigk(i),i=1,maxk1)
  angle = angle + dang
400 continue
do 401 k1=1,maxk1
  sigkto(k1) = sigkto(k1)*(dang/rad)*2.d0*pi
401 continue
write(2,352)
write(2,405) (k, sigkto(k+1),k=0, maxk)
write(6,406) (k, sigkto(k+1),k=0, maxk)
405 format(/,'%  INTEGRAL TENSOR CROSS SECTIONS'/ &
  '%     RANK  XS:  ',11(i6,1pe14.5))
406 format(/,'   INTEGRAL TENSOR CROSS SECTIONS'/ &
  '      RANK  XS:',11(i6,1pe14.5))
!
310 continue
!
return
end
!======
subroutine amplga(j1,inlev1,j2,inlev2,jtot,mmax,jpack, &
  lpack,ipack,length,jq,lq,inq,nopen,maxl2,nangle, &
  flaghf,sreal,simag,y,q)
!
! calculates scattering amplitudes for given jtot and set
! of angles
!
!     author of original ampli program:  h.-j. werner
!     revised by p.j.dagdigian for geomatric apse frame calculation
!     for the transition (j1,inlev1) -> (j2,inlev2)
!
!     revision date: 13-oct-2011
!
!.....jpack,lpack,ipack: labels for rows
!.....jq,lq,inq:         labels for columns
implicit double precision (a-h,o-z)
complex*16 q,ai,fak2,fak3,yy,tmat
parameter (zero=0.0d0,one=1.0d0)
logical flaghf,elastc
common /coang/ ang1, ang2, dang
dimension jpack(1),lpack(1),ipack(1),jq(1),lq(1),inq(1),q(1)
dimension sreal(mmax,1),simag(mmax,1),y(1)
dimension fak1(400),fak2(400),fak3(400),ilab1(400)
sqpi = 1.772453850905516d0
rad = 57.295779513082323d0
piov2 = 1.570796326794897d0
!
!.....ai is sqrt(-1)
ai=cmplx(zero, one)
elastc = j1.eq.j2 .and. inlev1.eq.inlev2
!
if (flaghf) then
!.....here for half-integer spin
  fakj=sqpi*(2.d0*jtot + 2.d0)*(-1)**(j1+j2+1)
  spin=0.5d0
  j1p=j1+1
  j2p=j2+1
else
!.....here for integer spin
  fakj=sqpi*(2.d0*jtot + 1.d0)*(-1)**(j1+j2)
  spin=0.0d0
  j1p=j1
  j2p=j2
end if
xjtot=dble(jtot)+spin
xj1=dble(j1)+spin
xj2=dble(j2)+spin
ll=0
do 50 ilab=1,length
  if(jpack(ilab).ne.j1 .or. ipack(ilab).ne.inlev1) &
     goto 50
  l1=lpack(ilab)
  ll=ll+1
  fak1(ll)=fakj*sqrt(2.d0*l1 + 1.0d0)
  ilab1(ll)=ilab
50 continue
llmax=ll
!
do 500 jlab=1,nopen
  if(jq(jlab).ne.j2 .or. inq(jlab).ne.inlev2) goto 500
  l2=lq(jlab)
  xl2=l2
  do 60 ll=1,llmax
    ilab=ilab1(ll)
    l1=lpack(ilab)
!.....convert to t-matrix
    tmat = -cmplx(sreal(jlab,ilab),simag(jlab,ilab))
    if(elastc .and. l1.eq.l2) tmat = tmat + 1.0d0
    fak2(ll)=cmplx(fak1(ll)*(-1)**(l1+l2),zero) &
       *(ai**(l1-l2))*tmat
60   continue
!
  ii=0
!
! mj1 is GA-frame initial m quantum number
  do 400 mj1=-j1p,j1
    xmj1=dble(mj1)+spin
! mj2 is GA-frame final m quantum number
    do 300 mj2=-j2p,j2
      xmj2=dble(mj2)+spin
      angle = ang1
      do 200 iang=1,nangle
        ii = ii + 1
! redrot below requires rotation angle in radians
        betaga = piov2 + 0.5d0*angle/rad
! mj1p is collision-frame initial m quantum number
        do 380 mj1p=-j1p,j1
          xmj1p=dble(mj1p)+spin
          yy = 0.d0
          do 70 ll=1,llmax
            ilab=ilab1(ll)
            xl1=lpack(ilab)
            fak3(ll)=fak2(ll) &
              *xf3j(xj1,xl1,xjtot,xmj1p,zero,-xmj1p)
70           continue
! mj2p is collision-frame final m quantum number
          do 350 mj2p=-j2p,j2
            xmj2p=dble(mj2p)+spin
            xml2p=xmj1p - xmj2p
            ml2p=xml2p
            iyof=(iabs(ml2p)*maxl2+l2)*nangle
            fakp = redrot(xj1,xmj1p,xmj1,betaga) &
              *redrot(xj2,xmj2p,xmj2,betaga) &
              *xf3j(xj2,xl2,xjtot,xmj2p,xml2p,-xmj1p) &
              *y(iyof+iang)
            if (ml2p.gt.0) fakp = fakp*(-1)**ml2p
            do 100 ll=1,llmax
              yy = yy + fak3(ll)*fakp
100             continue
370             continue
350           continue
          q(ii) = q(ii) + yy
380         continue
        angle = angle + dang
200       continue
300     continue
400   continue
500 continue
return
end
!=============
function redrot (rj,rk,rm,beta)
implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This function uses eq. (4.1.23) of Edmonds
!     to calculate the reduced rotation matrix element
!     d(j,k,m;beta) = <jk|exp(+i*beta*Jy/hbar)|jm>.
!
!     The angle beta is in radians
!     -----------------------------------------------------------------
!
parameter (zero = 0.0d0)
parameter (half = 0.5d0)
parameter (one  = 1.0d0)
parameter (two  = 2.0d0)
!
!     half integer angular momenta
!
sj = half*nint(two*rj)
sk = half*nint(two*rk)
sm = half*nint(two*rm)
!
!     projection ranges
!
redrot = zero
if (sk.gt.sj .or. sk.lt.-sj)  return
if (sm.gt.sj .or. sm.lt.-sj)  return
if (mod(sj-sk,one) .ne. zero) return
if (mod(sj-sm,one) .ne. zero) return
!
!     reflection symmetries
!
if (sk+sm .ge. zero) then
  if (sk-sm .ge. zero) then
    tk = sk
    tm = sm
    isign = 0
  else
    tk = sm
    tm = sk
    isign = sk-sm
  endif
else
  if (sk-sm .ge. zero) then
    tk = -sm
    tm = -sk
    isign = 0
  else
    tk = -sk
    tm = -sm
    isign = sk-sm
  endif
endif
!
!     evaluation
!
n = sj-tk
ia = tk-tm
ib = tk+tm
a = ia
b = ib
beta2 = half*beta
cosb2 = cos(beta2)
sinb2 = sin(beta2)
cosb = (cosb2-sinb2)*(cosb2+sinb2)
d1 = pjacob(n,a,b,cosb)
d2 = cosb2**ib*sinb2**ia
d3 = d1*d2
d4 = d3*d3
ti = tm
do i = 1,ia
   ti = ti+one
   d4 = d4*(sj+ti)/(sj-ti+one)
enddo
d4 = sqrt(d4)
redrot = sign(d4,d3)
if (mod(isign,2) .ne. 0) redrot = -redrot
return
end
!=============
function pjacob (n,a,b,x)
implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     Jacobi polynomial p(n,a,b;x)
!     Abramowitz and Stegun eq. (22.7.1)
!     -----------------------------------------------------------------
!
parameter (zero = 0.0d0)
parameter (half = 0.5d0)
parameter (one  = 1.0d0)
parameter (two  = 2.0d0)
!
if (n .eq. 0) then
  fp = one
else
  f = one
  apa = a+a
  apb = a+b
  amb = a-b
  apbamb = apb*amb
  apbp1 = apb+one
  apbp2 = apb+two
  onek = zero
  twok = zero
  fp = half*(amb+apbp2*x)
  do k = 1,n-1
    onek = onek+one
    twok = twok+two
    a1 = (twok+two)*(onek+apbp1)*(twok+apb)
    a2 = (twok+apbp1)*apbamb
    a3 = (twok+apb)*(twok+apbp1)*(twok+apbp2)
    a4 = (twok+apa)*(onek+b)*(twok+apbp2)
    fm = f
    f = fp
    fp = ((a2+a3*x)*f-a4*fm)/a1
  enddo
endif
pjacob = fp
return
end
!=====eof=====
