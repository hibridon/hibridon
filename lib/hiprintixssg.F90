!************************************************************************
!                                                                       *
!                      integral cross section print                     *
!                                                                       *
!************************************************************************
!                          routines included:                           *
!                                                                       *
!   prsg/aver1/xscpr1  integral cross section print
!                                                                       *
!************************************************************************
module mod_hiprintixssg
contains

subroutine prsg (fname, a)
!  subroutine to write out selected integral cross sections
!  from file {fname1}.ics
!  author:  millard alexander
!  last bug fix:17-may-1993 by mha
!  latest revision date:  23-feb-2013 by p. dagdigian
!  ------------------------------------------------------------------
!  variables in call list:
!    zmat:    on return:  contains the nlevop x nlevop matrix of integral
!                         cross sections
!    jlev:   rotational angular momenta of each energetically open level
!    elev:   energy (in hartree) of each energetically open level
!    inlev:  additonal quantum index for each energetically open level
!    jfirst:  initial value of total angular momentum
!    jfinal:  final value of total angular momentum
!    jtotd:   step size for total angular momentum
!    ipos:    if .true., 132 column printer
!             if .false., 80 column printer
!    csflag:  if .true. coupled-states calculation
!             if .false., close-coupled calculation
!    flaghf:  if .true., then system with half-integer spin
!             if .false., then system with integer spin
!    twomol:  if .true., then molecule-molecule collision
!             if .false., then atom-molecule or molecule-surface collision
!    nlevop:  number of energetically distinct levels in channel basis which
!             are open asymptotically
!    numin:   initial value of the coupled-states projection index
!    numax:   final value of the coupled-states projection index
!    nud:     step in coupled-states projection index
!    jlpar:   parity of channels
!               eps * (-1)**(j+l-jtot) = jlpar
!             if jlpar = 0, then integral cross sections include contributions
!             of both parities
!    note!!!   if flaghf = .true.( see below), then the true values
!    of the rotational quantum numbers, the total angular momentum,
!    and the coupled-states projection index are equal to the values
!    stored in jlev, jtot, jfirst, jfinal, nu, numin, and numax plus 1/2
!    flagsu:    if .true., then molecule-surface collisons
!    xmu:       collision reduced mass in (c12) atomic mass units
!    econv:     conversion factor from cm-1 to hartrees
!  ------------------------------------------------------------------
use mod_cosout, only: nnout, jout
use mod_coiout, only: niout, indout
use mod_cojhld, only: jlev => jhold ! jlev(5)
use mod_coisc1, only: inlev => isc1 ! inlev(5)
use mod_coisc2, only: jpoint => isc2 ! jpoint(5)
use mod_cosc1, only: elev => sc1 ! elev(5)
use mod_coz, only: zmat => z_as_vec ! zmat(1)
use mod_cow, only: scmat => w_as_vec ! scmat(1)
use mod_version, only : version
use mod_par, only: csflag, flaghf, flagsu, ipos
use mod_parpot, only: potnam=>pot_name, label=>pot_label
use mod_selb, only: ibasty
use mod_hiutil, only: gennam
use mod_hivector, only: matcopy
use funit, only: FUNIT_ICS, FUNIT_XSC
implicit double precision (a-h,o-z)
character*(*) fname
character*20 cdate
character*3 stat
character*12 accs
character*40 xnam1, xnam2
logical iprint, twomol, existf, openfl, eprint
!     real econv, ener, xmu
!     real a, elev, scmat, zmat
integer i, ienerg, iout, isize, j, jbegin, jend, jfinal, &
        jfirst, jj1, jj2, jlpar, jtemp, jtotd, lenx, n, nlevel, &
        nlevop, nout, numax, numin, nud, iaver
dimension  a(3)
data econv / 219474.6d0/
integer, allocatable :: ipoint(:)
!  input parameters
iprint=.true.
eprint=.false.
if (a(1) .lt. 0.0) iprint =  .false.
if (a(1) .gt. 0.0) eprint = .true.
iaver=nint(a(2))
ienerg = a(3) + 0.1
if (ienerg .le. 0) ienerg=1
!  open integral cross section file
call gennam (xnam1, fname, ienerg, 'ics', lenx)
inquire (file = xnam1, exist = existf)
if (.not. existf) then
  write (6, 10) xnam1(1:lenx)
10   format(/' integral cross section file ',(a),' not found',/)
  return
end if
open (FUNIT_ICS, file = xnam1,status = 'old')
!  open output file for integral cross sections
call gennam (xnam2, fname, ienerg, 'xsc', lenx)
if (.not.iprint) write(6,20) xnam2
20 format(' PRINTING SELECTED CROSS SECTIONS; OUTPUT IN FILE ', &
        (a))
! inquire file specifications
inquire(file=xnam2, exist=existf, opened=openfl)
accs='sequential'

if (.not. openfl) then
  if (existf) then
    stat='old'
! make sure sequential formatted files are appended not overwritten
#if defined(HIB_UNIX_IFORT) || defined(HIB_UNIX_IFX) || defined(HIB_UNIX_PGI)
    accs='append'
#endif
#if defined(HIB_UNIX_HP) || defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IRIS) || defined(HIB_UNIX_SUN)
    accs='append'
#endif
  else
    stat='new'
  end if
  open (FUNIT_XSC, file = xnam2, status = stat, access = accs)
endif
call version(FUNIT_XSC)
read (FUNIT_ICS, 40) cdate
40 format (1x, a)
read (FUNIT_ICS, 40) label
read (FUNIT_ICS, 40) potnam
!  print job information
write (FUNIT_XSC, 50) xnam1, cdate, label, potnam
if (iprint) write (6, 50) xnam1, cdate, label, potnam
50 format(/' INTEGRAL CROSS SECTIONS READ FROM FILE ',(a)/ &
        ' WRITTEN:   ',(a)/ &
        ' LABEL:     ',(a)/, &
        ' POT NAME:  ',(a) )
read (FUNIT_ICS, 60) ener, xmu
60 format (f10.3, f15.11)
read (FUNIT_ICS, 70) csflag, flaghf, flagsu, twomol
70 format (4l3)
read (FUNIT_ICS, 75) jfirst, jfinal, jtotd, numin, numax, &
             nud, jlpar
75 format (7i5)
if (ipos) then
80   format (24i5)
  read (FUNIT_ICS, 80)  nlevel, nlevop
  read (FUNIT_ICS, 80) (jlev(i), inlev(i), i = 1, nlevel)
  read (FUNIT_ICS, 90) (elev(i), i = 1, nlevel)
90   format (8(1pe15.8))
!90      format (8f16.9)
else
  read (FUNIT_ICS, *)  nlevel, nlevop
  read (FUNIT_ICS, *) (jlev(i), inlev(i), i = 1, nlevel)
  read (FUNIT_ICS, *) (elev(i), i = 1, nlevel)
endif
!  read in matrix of cross sections, column by column
do 95  i = 1, nlevop
!  jbegin and jend point to first element and last element of column i of
!  matrix packed column by column
  jbegin = (i - 1) * nlevop + 1
  jend = i * nlevop
  if (ipos) then
    read (FUNIT_ICS, 90) (zmat(j), j = jbegin, jend)
  else
    read (FUNIT_ICS, *) (zmat(j), j = jbegin, jend)
  endif
95 continue
if (.not. flagsu) then
  if (.not.flaghf .or. ibasty.eq.12) then
    if (iprint) write (6, 100) ienerg, xmu, ener, jlpar, &
                                jfirst, jfinal, jtotd
    write (FUNIT_XSC, 100) ienerg, xmu, ener, jlpar, jfirst, &
                   jfinal, jtotd
100     format (/' ** INTEGRAL CROSS SECTIONS; IEN=', i2,' **', &
            /'    RMU=', f9.4, '  E=', f9.2,'  JLPAR=', i2, &
             '  JTOT-1=', i3, &
             '  JTOT-2=', i4,'  JTOT-D=', i3)
  else
    if (iprint) write (6, 105) ienerg, xmu, ener, jlpar, &
                      (jfirst+0.5), (jfinal+0.5), jtotd
    write (FUNIT_XSC, 105) ienerg, xmu, ener, jlpar, (jfirst+0.5), &
                   (jfinal+0.5), jtotd
105     format (/' ** INTEGRAL CROSS SECTIONS; IEN=', i2,' **', &
            /'    RMU=', f9.4, '  E=', f9.2,'  JLPAR=', i2, &
             '  JTOT-1=', f5.1, &
             '  JTOT-2=', f6.1,'  JTOT-D=', i3)
  end if
else
  if (.not. flaghf) then
    if (iprint) write (6, 110) ienerg, xmu, ener, numin, &
                               numax, nud
    write (FUNIT_XSC, 110) ienerg, xmu, ener, numin, numax, nud
110     format (/' ** SUMMED DEGENERACY AVERAGED TRANSITION', &
             ' PROBABILITIES;  IEN=', i2,' **', &
            /'    RMU=', f9.4, '  E=', f8.2,'  M-MIN=', i3, &
             '  M-MAX=', i4, '  M-STEP=', i2)
  else
    if (iprint) write (6, 115) ienerg, xmu, ener, (numin+0.5), &
                               (numax+0.5), nud
    write (FUNIT_XSC, 115) ienerg, xmu, ener, (numin+0.5), &
                   (numax+0.5), nud
115     format (/' ** SUMMED DEGENERACY AVERAGED TRANSITION', &
             ' PROBABILITIES;  IEN=', i2,' **', &
            /'    RMU=', f9.4, '  E=', f8.2,'  M-MIN=', f5.1, &
             '  M-MAX=', f6.1, '  M-STEP=', i2)
  end if
end if
if (.not. csflag) then
  if (jlpar .eq. 0) then
    write (FUNIT_XSC, 120)
    if (iprint) write (6, 120)
120     format (' ** CC CALCULATION, BOTH PARITIES **')
  else
    write (FUNIT_XSC, 125) jlpar
    if (iprint) write (6, 125) jlpar
125     format (' ** CC CALCULATION, JLPAR=', i2, ' **')
  end if
else
  if (.not. flaghf) then
    write (FUNIT_XSC, 130) numin, numax, nud
    if (iprint) write (6, 130) numin, numax, nud
130     format (' ** CS CALCULATION, NUMIN=', i2,', NUMAX=', &
               i2,' NUD=', i2, ' **')
  else
    write (FUNIT_XSC, 140) numin + 0.5, numax + 0.5, nud
    if (iprint) write (6, 140) numin+0.5, numax+ 0.5, nud
140     format (' ** CS CALCULATION, NUMIN=', f4.1, ', NUMAX=', &
            f4.1,' NUD=', i2, ' **')
  end if
end if
if (eprint) then
  if (.not. twomol) then
    if (.not. flagsu) then
      if (iprint) write (6, 145)
      write (FUNIT_XSC, 145)
145       format &
       (' LEVEL LIST FOR INTEGRAL CROSS SECTIONS', &
        /'   N   J  INDEX  EINT(CM-1)',/)
    else
      if (iprint) write (6, 150)
      write (FUNIT_XSC, 150)
150       format &
       (' LEVEL LIST FOR DEGENERACY AVERAGED', &
         ' TRANSITION PROBABILITIES', &
        /'   N   J  INDEX  EINT(CM-1)',/)
    end if
    do 170  i = 1, nlevop
      if (.not. flaghf) then
        if (iprint) &
        write (6, 160) i, jlev(i), inlev(i), elev(i) * econv
        write (FUNIT_XSC, 160) i, jlev(i), inlev(i), elev(i) * econv
160         format (i4, i5, i6, f11.3)
      else
        if (iprint) write (6, 165) i, (jlev(i)+0.5), inlev(i), &
                       elev(i) * econv
        write (FUNIT_XSC, 165) i, (jlev(i)+0.5), inlev(i), &
                       elev(i) * econv
165         format (i4, f5.1, i6, f11.3)
      end if
170     continue
  else
    if (iprint) write (6, 175)
    write (FUNIT_XSC, 175)
175       format (' LEVEL LIST FOR INTEGRAL CROSS SECTIONS', &
              /'   N   J1   J2  INDEX  EINT(CM-1)'/)
    do 190  i = 1, nlevop
      jj2 = mod( jlev(i), 10)
      jj1 = jlev(i) / 10
      if (iprint) &
      write (6, 180) i, jj1, jj2, inlev(i), elev(i) * econv
      write (FUNIT_XSC, 180) i, jj1, jj2, inlev(i), elev(i) * econv
180       format (i4, 2i5, i6, f11.3)
190     continue
  end if
endif
!  now sum and average over positive and negative values of index
!     check that number of levels with index negative equals number of levels
!     with index positive, if not abort
!     this test works by (1) checking that no index equals zero
!     and, if so, (2) adding all the values of the index.  this sum
!     should equal zero if there are as many levels with index negative
!     as index positive
if (iaver .gt. 0) then
  isum = 0
  do  200 i = 1, nlevop
    if (inlev(i) .eq. 0) then
      write (6, 195) i
      write (6, 195) i
195       format(' *** INLEV(',i3,')=0;', &
               ' AVERAGING MAY NOT WORK ***')
    else if (inlev(i) .ne. 0) then
      isum = isum + inlev(i)
    end if
200   continue
  if (isum .ne. 0) then
    write (6, 230) isum
    write (FUNIT_XSC, 230) isum
230     format (' *** SUM OF INDICES =',i3, &
               ' AVERAGING MAY NOT WORK ***')
  end if
  if  (iaver .eq. 2) then
    if (iprint) write (6, 245)
    write (FUNIT_XSC, 245)
245     format &
      (' ** CROSS SECTIONS SUMMED AND AVERAGED OVER INDEX **')
    call aver1 (zmat, scmat, nlevop)
   else if (iaver .eq. 1) then
    if (iprint) write (6, 250)
    write (FUNIT_XSC, 250)
250     format &
      (' ** CROSS SECTIONS SUMMED OVER FINAL STATE INDEX **')
   end if
end if
!  find all rows of cross sections matrix for which initial rotational
!  quantum number is equal to one of the values of jout
isize = 0
allocate(ipoint(nlevop))
nout = abs (nnout)
do  280  iout = 1, nout
  jtemp = jout(iout)
  do 270 n = 1, nlevop
    if (jlev(n) .eq. jtemp) then
      isize = isize + 1
      jpoint(isize) = n
      ipoint(isize) = n
    end if
270   continue
280 continue
!aber
insize=0
!aber
if (niout .gt. 0) then
  nout=abs(niout)
  do  282  iout = 1, nout
    indtemp=indout(iout)
    do 281 n = 1, isize
      if (inlev(jpoint(n)) .eq. indtemp) then
        insize = insize + 1
        ipoint(insize) = jpoint(n)
      end if
281     continue
282   continue
  isize=insize
endif
!  isize is the number of cross sections to be printed
if (isize .eq. 0) then
!  here if no initial states found
  write (6, 283)
  write (FUNIT_XSC, 283)
283   format (' ** NO INITIAL STATES FOUND; ABORT')
else
!  now print out desired columns of cross section matrix
  write (FUNIT_XSC, 260)
  if (iprint) write (6, 260)
260   format (/' ** COLUMN HEADINGS ARE INITIAL STATES, ROW', &
        ' HEADINGS ARE FINAL STATES **')
  call xscpr1(zmat, nlevop, isize, iaver, ipos, iprint, flaghf, FUNIT_XSC, ipoint)
endif
deallocate(ipoint)
close (1)
close (3)
return
end
subroutine aver1 (zmat, scmat, n)
use mod_cojhld, only: jlev => jhold ! jlev(1)
use mod_coisc1, only: inlev => isc1 ! inlev(1)
use mod_himatrix, only: transp
use mod_hivector, only: matcopy
implicit double precision (a-h,o-z)
!  subroutine to sum and average cross section matrix over positive
!  and negative values of index
integer i, ind, index, j, jnd, n, nn
!     real scmat, zmat
dimension zmat(1), scmat(n,n)
!  first copy cross section matrix into scmat
call matcopy (zmat, scmat, n, n, n, n)
nn = n / 2
index = 0
do  30  i = 1, nn
  ind = 2 * i
  do  20  j = 1, nn
    jnd = 2 * j
    index = index + 1
    zmat(index) = ( scmat(ind - 1, jnd - 1) + &
                   scmat(ind - 1, jnd) + &
                   scmat(ind, jnd - 1) + &
                   scmat(ind, jnd) ) * 0.5
20   continue
   jlev (i) = jlev (ind - 1)
  inlev(i) = iabs(inlev(ind - 1))
30 continue
!  because of indexing zmat need to be transposed
call transp (zmat, nn, nn)
!  reduce the size of the cross section matrix
n = nn
return
end
! -------------------------------------------------
subroutine xscpr1 (zmat, nlevop, isize, iaver, &
                   ipos, iprint, flaghf, xsc_funit, ipoint)
use mod_cojhld, only: jlev => jhold ! jlev(4)
use mod_coisc1, only: inlev => isc1 ! inlev(4)
use mod_selb, only: ibasty
use mod_himatrix, only: transp

implicit double precision (a-h,o-z)
integer, intent(in) :: xsc_funit  ! the file unit for xsc file (it's expected to be open in write mode)
integer, intent(in) :: ipoint(*)
!      current revision date: 16-dec-2007
!  subroutine to print out specified columns of cross section matrix
!  if iaver = 1, then nth and (n+1) st rows are added before printing
integer i, isize, iskip, j, jcol, jhigh, jj, jlow, jmax, &
        jrow, ncol, nlevop, iaver
integer ind
!     real elev, zmat
logical ipos, iprint, flaghf
dimension zmat(nlevop,nlevop), ind(50)
!     first transpose the cross section matrix so that initial
!     states are columns and final states are rows
call transp (zmat, nlevop, nlevop)
!  if 132 line printer, then 13 columns of matrix are printed simultaneously
!  if  80 line printer, then  7 columns of matrix are printed simultaneously
iskip = 13
if (.not. ipos) iskip = 7
irowsk = 1
if (iaver .eq. 1) irowsk = 2
!  jmax is the total number of the groups of columns to be printed
if (mod(isize, iskip) .eq. 0) then
  jmax = isize / iskip
else
  jmax = isize / iskip + 1
end if
jlow = 1
jhigh = min (iskip, isize)
!  loop over columns by groups of 6 or 10
do  150   j = 1, jmax
!  ncol is the number of columns to be printed
  ncol = jhigh - jlow + 1
  do  10   jj = 1, ncol
!  the array ind contains the index of each column which will be printed
    ind(jj) = ipoint(jj - 1 + jlow)
10   continue
!  write as a heading the column index of each column to be printed
  if (.not. flaghf .or. ibasty.eq.12) then
    if (iprint) write (6, 15) ( jlev(ind(i)), i = 1,ncol)
    write (xsc_funit, 15) ( jlev(ind(i)), i = 1, ncol )
15     format (/12x,'J=', i5, 2x, 12 (2x, i6, 2x) )
  else
    if (iprint) write (6, 30) ( jlev(ind(i))+0.5, &
                                 i = 1,ncol)
    write (xsc_funit, 30) ( jlev(ind(i))+0.5, i = 1,ncol)
30     format (/11x,'J= ', f5.1, 2x, 12 (2x, f6.1, 2x) )
  end if
  if (iprint) write (6, 40) ( inlev(ind(i)), i = 1,ncol)
  write (xsc_funit, 40) ( inlev(ind(i)), i = 1,ncol)
40   format ('   J    I | I=', i5, 2x, 12 (2x, i6, 2x))
  if (iprint) write (6, 50)
  write (xsc_funit, 50)
50   format (1h )
!  now loop through the rows of the matrix, which will be printed out
!  in groups of 10 with a blank line in between
!  loop over the rows of the matrix which will be printed
  do  95   jrow =  1 , nlevop, irowsk
!  now write out the row index followed by the desired matrix elements
    jj=jrow/10
!  inrow holds additional quantum index for this row
!  if iaver = 1, then this is positive, since both indices are summed
    inrow = inlev(jrow)
    if (iaver .ne. 1) then
      if (.not. flaghf .or. ibasty.eq.12) then
        if (iprint) &
          write (6, 60) jlev(jrow),inrow, &
            ( zmat(jrow,ind(jcol)), jcol = 1,ncol)
        write (xsc_funit, 60) jlev(jrow),inrow, &
          ( zmat(jrow,ind(jcol)), jcol = 1,ncol)
60         format (i5, i5, 2x, 13 (1pe10.3,1x) )
      else
        if (iprint) &
          write (6, 70) jlev(jrow)+0.5, inrow, &
            ( zmat(jrow,ind(jcol)), jcol = 1,ncol)
        write (xsc_funit, 70) jlev(jrow)+0.5, inrow, &
          ( zmat(jrow,ind(jcol)), jcol = 1,ncol)
70         format (f5.1, i5, 2x, 13 (1pe10.3,1x) )
      end if
    else
      if (.not. flaghf .or. ibasty.eq.12) then
        if (iprint) &
          write (6, 60) jlev(jrow),inrow, &
            (zmat(jrow,ind(jcol)) + zmat(jrow + 1,ind(jcol)), &
                jcol = 1,ncol)
        write (xsc_funit, 60) jlev(jrow),inrow, &
            (zmat(jrow,ind(jcol)) + zmat(jrow + 1,ind(jcol)), &
                jcol = 1,ncol)
      else
        if (iprint) &
          write (6, 70) jlev(jrow)+0.5, inrow, &
            (zmat(jrow,ind(jcol)) + zmat(jrow + 1,ind(jcol)), &
                jcol = 1,ncol)
        write (xsc_funit, 70) jlev(jrow)+0.5, inrow, &
            (zmat(jrow,ind(jcol)) + zmat(jrow + 1,ind(jcol)), &
                jcol = 1,ncol)
      end if
    end if
!  if this row is an integer multiple of 10, add a blank line
!  to output on unit 6 only
    jj = jrow / 10
    if (10 * jj .eq. jrow .and. jrow .ne. nlevop) then
      if (iprint) write (6, 80)
    end if
80     format(1h )
95   continue
!  all rows have been printed out for this group of columns
!  move to the next group of columns
!  jlow becomes old jhigh
!  jhigh is set equal to jhigh plus skip distance
  if (jhigh .eq. nlevop) go to 160
  jlow = jhigh + 1
  jhigh = min ( (jlow + iskip - 1), isize)
150 continue
!  entire matrix has been printed, return
160 return
end

end module mod_hiprintixssg
