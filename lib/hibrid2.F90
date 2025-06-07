!***********************************************************************
!                                                                       *
!                         hibridon 2  library                           *
!                                                                       *
!************************************************************************
!                          routines included:                           *
!                                                                       *
!   1. hiblk       block data, default settings                         *
!   2. hinput      input driver                                         *
!   4. mxoutd (mxoutr) print_integral_cross_sections      matrix print utility
!   5. prsg/aver1/xscpr1  integral cross section print
!   6. prsgpi/aver2/xscpr2 integral cross section print for sigma/pi
!                                                                       *
!************************************************************************
module mod_hibrid2
contains
! ------------------------------------------------------------------
subroutine set_default_params
!  current revision date:  22-jan-2008 by mha
! ------------------------------------------------------------------
use mod_cosout, only: nnout, jout
use mod_coiout, only: niout, indout
use mod_coener, only: energ
use mod_cosysi, only: nscode
use mod_cosysl, only: islcod
use mod_cosysr, only: isrcod
use mod_par, only: airyfl, prairy, bastst, batch, chlist, &
                csflag, flaghf, flagsu, ihomo, ipos, logdfl, &
                prlogd, noprin, prpart, readpt, rsflag, prsmat, &
                t2test, prt2, twomol, wrsmat, wrpart, wrxsec, &
                prxsec, nucros, photof, wavefl, boundc, &
                jtot1, jtot2, jtotd, jlpar, nerg, numax, numin, nud, lscreen, iprint, &
                fstfac=>scat_fstfac, rincr=>scat_rincr, rcut=>scat_rcut, rendai=>scat_rendai, rendld=>scat_rendld, rstart=>scat_rstart, spac=>scat_spac, xmu
use mod_parpot, only: potnam=>pot_name, label=>pot_label
use mod_selb, only: ibasty
use mod_file, only: input, output, jobnam
use mod_com, only: com
implicit double precision (a-h,o-z)
! nb if the nextcommon is changed, it should be also changed in common/parsys
!  this sets the maximum number of energies
!
jtot1=20
jtot2=20
jtotd=20
jlpar=1
nerg=1
numax=0
numin=0
nud=1
lscreen=48
iprint=0
energ(1)=208.509d0
do i=2,25
  energ(i)=0.d0
enddo
fstfac=15d0
rincr=8d0
rcut=30d0
rendai=25d0
rendld=8d0
rstart=5.6d0
spac=0.15d0
tolai=1.15d0
xmu=16.47d0
! econv is the conversion from hatree->cm-1
! xmconv is the conversion from au (electron mass)-> amu, i.e. the mass of the
! in amu
! ang2c is bohr^2 in angstroms^2
! default basis type is 1 (singlet sigma)
ibasty=1
airyfl=.true.
prairy=.false.
! If the command file is provided through the -c or --com command line option
! then batch mode is set to true.
batch = com
chlist=.true.
csflag=.false.
flaghf=.false.
flagsu=.false.
ihomo=.true.
ipos=.false.
logdfl=.true.
prlogd=.false.
noprin=.false.
prpart=.false.
readpt=.false.
rsflag=.false.
prsmat=.false.
t2test=.false.
prt2=.true.
wrsmat=.true.
wrpart=.false.
wrxsec=.false.
nucros=.false.
prxsec=.true.
bastst=.false.
photof=.false.
wavefl=.false.
boundc=.false.
twomol=.false.
nnout=-2
do i=1,20
  jout(i)=0
enddo
jout(2)=2
! the following values are for a singlet sigma molecule in v=0
niout=1
indout(1)=0
! the following values are for a doublet sigma type molecule
!      data niout/2/,indout(1),indout(2)/1,-1/
nscode=0
isicode=0
islcod=0
isrcod=0
!  set up default names for input and output files and jobname
input='Inpt'
output='Outpt'
jobnam='Job'
! default labels
label=' LABEL NOT YET ASSIGNED'
potnam=' POTNAM NOT YET ASSIGNED'

end
! ----------------------------------------------------------------------
subroutine enord(energ,nerg)
implicit double precision (a-h,o-z)
dimension energ(1)
do 502 i=1,nerg
emax=energ(i)
jmax=i
do 501 j=i+1,nerg
if(energ(j).gt.emax) then
  emax=energ(j)
  jmax=j
end if
501 continue
if(jmax.ne.i) then
  energ(jmax)=energ(i)
  energ(i)=emax
end if
502 continue
return
end
! ----------------------------------------------------------------------
subroutine mxoutd (iunit, xmat, nn, nmax, isym, ipos)
!  to print out nn*nn matrix xmat of maximum dimension nmax use mxoutd
!  to print out nn*m matrix xmat of maximum dimension nmax use mxoutr
!  author:  millard alexander
!  current revision date: 1-may-90 by mha
!  -------------------------------------------------------------------------
!  variables in call list:
!    iunit:     logical unit number
!    xmat:      matrix to be printed out
!    nn:         actual size of matrix (nn x nn)
!    nmax:      maximum row dimension of matrix
!    isym:      if isym = 1,  matrix is symmetrical, only lower triangle
!               is printed
!               if isym = 0 full matrix is printed
!               note that isym must be either 0 or 1 in calling sequence in ma
!               program
!    ipos:      if .true., printer is assumed to have 132 positions/line
!               if .false., printer is assumed to have 80 positions/line
!  -------------------------------------------------------------------------
implicit double precision (a-h,o-z)
logical ipos
dimension xmat(nmax,nmax), ind(10)
m=nn
n=nn
goto 1
entry mxoutr (iunit, xmat, nrow, ncol, nmax, isym, ipos)
m=nrow
n=ncol
!  abort if isym not equal to zero or one
1 if (isym .ne. 0 .and. isym .ne. 1) then
  write (iunit, 5)
  if(iunit.ne.6) write (6, 5)
5   format (/' ***  ISYM .NE. 0 OR 1 IN MXOUTD; ABORT ***')
  call exit
end if
!  if 132 line printer, then 10 columns of matrix are printed simultaneously
!  if  80 line printer, then  6 columns of matrix are printed simultaneously
iskip = 10
if (.not. ipos) iskip = 6
!  jmax is the total number of the groups of columns to be printed
jmax = n / iskip + 1
jlow = 1
jhigh = min (iskip, n)
!  loop over the groups of columns from column jlow to jhigh
do  50   j = 1, jmax
!  ncol is the number of columns to be printed
  nc = jhigh - jlow + 1
  do  10   jj = 1, nc
!  the array ind contains the index of each column which will be printed
    ind(jj) = jj - 1 + jlow
10   continue
!  write as a heading the column index of each column to be printed
  write (iunit, 15) ( ind(i), i = 1, nc )
15   format (/5x, 10 (3x, i3, 6x) )
!  now loop through the rows of the matrix, which will be printed out
!  in groups of 10 with a blank line in between
!  lowrow is the index of the first row to be printed
!  if the full matrix is desired (isym = 0), then lowrow = 1
!  if the lower triangle is desire (isym = 1), then lowrow = jlow
  lowrow = isym * jlow + (1 - isym)
!  loop over the rows of the matrix which will be printed
  do  35   jrow =  lowrow, m
!  jtop points to the maximum column which will be printed in this row
    jtop = isym * min (jhigh, jrow) + (1 - isym) * jhigh
!  now write out the row index followed by the desired matrix elements
    jj=jrow/10
    write (iunit,20) jrow,( xmat(jrow,jcol),jcol = jlow,jtop)
20     format (i4, 1x, 10 (1pe12.4) )
!  if this row is an integer multiple of 10, add a blank line
    jj = jrow / 10
    if (10 * jj .eq. jrow .and. jrow .ne. m) write (9, 25)
25     format(1h )
35   continue
!  all rows have been printed out for this group of columns
!  move to the next group of columns
!  jlow becomes old jhigh
!  jhigh is set equal to jhigh plus skip distance
  if (jhigh .eq. n) go to 60
  jlow = jhigh + 1
  jhigh = min ( (jlow + iskip - 1), n)
50 continue
!  entire matrix has been printed, return
60 return
end
! --------------------------------------------
! routine to print integral cross sections
subroutine print_integral_cross_sections(ifile, zmat, nlevop, nmax, ipos, csflag, flaghf, twomol, numax, jlev, inlev)
use constants
use mod_coisc2, only: nj, jlist => isc2 ! nj,jlist(1)
use mod_par, only: iprint
implicit none
integer, intent(in) :: ifile  ! the output file's unit (expected to be open)
real(8), intent(in) :: zmat(nmax, nlevop)
integer, intent(in) :: nlevop
integer, intent(in) :: nmax
logical, intent(in) :: ipos
logical, intent(in) :: csflag
logical, intent(in) :: flaghf
logical, intent(in) :: twomol
integer, intent(in) :: numax
integer, intent(in) :: inlev(nlevop)
integer, intent(in) :: jlev(nlevop)

integer :: i, ii, j, j1, j2
real(8) :: spin
logical csff
character*40 form
spin=0
csff=csflag
if (csflag .and. iprint .ge. 2) then
! here for full print of cs cross sections, even if not converged
  write (ifile,1) numax
1   format(/'** COUPLED STATES CROSS SECTIONS NOT CONVERGED', &
         ' FOR J .GT. NUMAX = ',i3)
  csff=.false.
endif
if(twomol) then
  if(ipos) write(ifile,5) (j,j=1,nlevop)
  if(.not.ipos) write(ifile,6) (j,j=1,nlevop)
5   format(/1x,'  INITIAL STATE',t27,'FINAL STATE'/ &
          1x,'  N   J1   J2  INDEX',(t21,i7,9i11))
6   format(/1x,'  INITIAL STATE',t27,'FINAL STATE'/ &
          1x,'  N   J1   J2  INDEX',(t21,i7,4i11))
! don't remove blank in this format!!
  form='(1x,i3,2i5,i5, (t21,10(1pd11.3)))'
else
  if(ipos) write(ifile,10) (j,j=1,nlevop)
  if(.not.ipos) write(ifile,11) (j,j=1,nlevop)
10   format(/1x,'  INITIAL STATE',t24,'FINAL STATE'/ &
          1x,'  N   J   INDEX ',(t18,i7,9i11))
11   format(/1x,'  INITIAL STATE',t24,'FINAL STATE'/ &
          1x,'  N   J   INDEX ',(t18,i7,4i11))
  if(flaghf) then
    spin=0.5d0
    form='(1x,i3,f5.1,i6,(t18,10(1pd11.3)))'
  else
! don't remove blank in this format!!
    form='(1x,i3,i4,i7,  (t18,10(1pd11.3)))'
  end if
end if
if(.not.ipos) form(21:22)=' 5'
if(csff) then
  nj=0
  do 100 i=1,nlevop
  if(jlev(i).le.numax) then
    nj=nj+1
    jlist(nj)=i
  end if
100   continue
  do 120 ii=1,nj
  i=jlist(ii)
  if(twomol) then
    j2 = mod( jlev(i), 10)
    j1 = jlev(i) / 10
    write(ifile,form) i,j1,j2,inlev(i),(zmat(i,j),j=1,nlevop)
  else if(flaghf) then
    write(ifile,form) i,jlev(i)+spin,inlev(i), &
                     (zmat(i,j),j=1,nlevop)
  else
    write(ifile,form) i,jlev(i),inlev(i),(zmat(i,j),j=1,nlevop)
  end if
120   continue
  if(ipos) write(ifile,10) (jlist(j),j=1,nj)
  if(.not.ipos) write(ifile,11) (jlist(j),j=1,nj)
  do 140 i=1,nlevop
  if(jlev(i).le.numax) goto 140
  if(twomol) then
    j2 = mod( jlev(i), 10)
    j1 = jlev(i) / 10
    write(ifile,form) i,j1,j2,inlev(i),(zmat(i,jlist(j)),j=1,nj)
  else if(flaghf) then
    write(ifile,form) i,jlev(i)+spin,inlev(i), &
                     (zmat(i,jlist(j)),j=1,nj)
  else
    write(ifile,form) i,jlev(i),inlev(i),(zmat(i,jlist(j)),j=1,nj)
  end if
140   continue
else
  do 150 i=1,nlevop
  if(twomol) then
    j2 = mod( jlev(i), 10)
    j1 = jlev(i) / 10
    write(ifile,form) i,j1,j2,inlev(i),(zmat(i,j),j=1,nlevop)
  else if(flaghf) then
    write(ifile,form) i,jlev(i)+spin,inlev(i), &
                      (zmat(i,j),j=1,nlevop)
  else
    write(ifile,form) i,jlev(i),inlev(i),(zmat(i,j),j=1,nlevop)
  end if
150   continue
end if
return
end
!  ------------------------------------------------------------------
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
!  ------------------------------------------------------------------
subroutine prsgpi (fname, a)
!  subroutine to write out selected integral cross sections
!  from file {fname1}.ics for sigma - pi transitions
!  author:  millard alexander
!  latest revision date:  5-apr-2004 by mha
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
!      implicit none
use mod_cosout, only: nnout, jout
use mod_coiout, only: niout, indout
use constants
use mod_coamat, only: zbuf ! zbuf(1)
use mod_cojhld, only: jlev => jhold ! jlev(1)
use mod_coisc1, only: inlev => isc1 ! inlev(1)
!mha
use mod_coisc2, only: jpoint => isc2 ! jpoint(5)
!mha
use mod_cosc1, only: elev => sc1 ! elev(1)
use mod_coener, only: energ
use mod_coz, only: zmat => z_as_vec ! zmat(1)
use mod_cow, only: scmat => w_as_vec ! scmat(1)
use mod_par, only: csflag, flaghf, flagsu, ihomo, ipos
use mod_parpot, only: potnam=>pot_name, label=>pot_label
use mod_hiutil, only: gennam
use mod_hivector, only: dset
use funit, only: FUNIT_ICS, FUNIT_XSC
implicit double precision (a-h,o-z)
character*(*) fname
character*20 cdate
character*40 xnam1, xnam2
character*80 line
!mha
character*3 stat
character*12 accs
!mha
logical iprint, twomol, existf, &
        openfl, eprint
dimension  a(4)
integer, allocatable :: ipoint(:)
!  input parameters
!mha
iprint=.true.
eprint=.false.
iaver = 0
!mha
!ABER
if (a(1) .lt. 0.0d0) iprint = .false.
!mha
if (a(1) .gt. 0.d0) eprint = .true.
if (a(2) .lt. 0.5) iaver = 0
if (a(2).gt. 0.d0) iaver=1
if (a(2) .gt. 1.5) iaver = 2
!ABER
ienerg = a(3) + 0.1
xthresh=a(4)
iener=energ(1)
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
! inquire file specifications
inquire(file=xnam2, exist=existf, opened=openfl)
!mha .xsc file is appended if it already exists
if (.not. openfl) then
  accs='sequential'
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
  open (FUNIT_XSC, file = xnam2, status = stat,access = accs)
endif
!mha
!mha (print out message only if no other print out)
if (.not.iprint) write(6,20) xnam2
20 format( &
  ' PRINTING SELECTED SIGMA-PI CROSS SECTIONS; OUTPUT IN FILE ', &
        (a))
!mha
!      call version(3)
read (FUNIT_ICS, 40) cdate
40 format (1x, a)
read (FUNIT_ICS, 40) label
read (FUNIT_ICS, 40) potnam
!  print job information
write (FUNIT_XSC, 50) xnam1, cdate, label, potnam
if (iprint) write (6, 50) xnam1, cdate, label, potnam
50 format( &
!mha
  /'% INTEGRAL SIGMA-PI CROSS SECTIONS READ FROM FILE ',(a)/ &
!mha
        '% WRITTEN:   ',(a)/ &
        '% LABEL:     ',(a)/, &
        '% POT NAME:  ',(a) )
read (FUNIT_ICS, 60) ener, xmu
60 format (f10.3, f15.11)
read (FUNIT_ICS, 70) csflag, flaghf, flagsu, twomol, ihomo
70 format (5l3)
read (FUNIT_ICS, 80) jfirst, jfinal, jtotd, numin, numax, &
             nud, jlpar, isa
80 format (24i5)
if (ipos) then
  read (FUNIT_ICS, 80)  nlevel, nlevop
  read (FUNIT_ICS, 80) (jlev(i), inlev(i), i = 1, nlevel)
  read (FUNIT_ICS, 90) (elev(i), i = 1, nlevel)
90   format (8f16.9)
else
  read (FUNIT_ICS, *)  nlevel, nlevop
  read (FUNIT_ICS, *) (jlev(i), inlev(i), i = 1, nlevel)
  read (FUNIT_ICS, *) (elev(i), i = 1, nlevel)
endif
!aber_begin
!  zero out zmat
call dset(nlevop*nlevop,0.d0, zmat,1)
!aber_end
!  read in matrix of cross sections, column by column
99 do 95  i = 1, nlevop
!  jbegin and jend point to first element and last element of column i of
!  matrix packed column by column
  jbegin = (i - 1) * nlevop + 1
  jend = i * nlevop
!            write (nxfile, 90) (zmat(j,i), j = 1, nlevop)
!  chnage format of reac statement (p. dagdigian, 7-mar-2012)
!        read (FUNIT_ICS, 90) (zbuf(j), j = jbegin, jend)
  read (FUNIT_ICS, *) (zbuf(j), j = jbegin, jend)
  do 96 j=jbegin,jend
    zmat(j)=zmat(j)+zbuf(j)
96   continue
95 continue
read (FUNIT_ICS,'(a)',end=999) line
if (line(1:14).eq.' ** RESTART **') then
  read (FUNIT_ICS,40) cdate
  read (FUNIT_ICS,40) label
  read (FUNIT_ICS,40) potnam
  read (FUNIT_ICS, 60) ener, xmu
  read (FUNIT_ICS, 70) csflag, flaghf, flagsu, twomol, ihomo
  read (FUNIT_ICS, 80) idum, jfinal, jtotd, numin, numax, &
               jlpar, isa
  read (FUNIT_ICS, 80)  nlevel, nlevop
  read (FUNIT_ICS, 80) (jlev(i), inlev(i), i = 1, nlevel)
  read (FUNIT_ICS, 90) (elev(i), i = 1, nlevel)
  goto 99
end if
999 if (.not. flagsu) then
  if (.not. flaghf) then
    if (iprint) write (6, 100) ienerg, xmu, ener, jlpar, &
                                jfirst, jfinal, jtotd
    write (FUNIT_XSC, 100) ienerg, xmu, ener, jlpar, jfirst, &
                   jfinal, jtotd
100     format ('% ** INTEGRAL CROSS SECTIONS; IEN=', i2,' **', &
            /'%    RMU=', f9.4, '  E=', f9.2,'  JLPAR=', i2, &
             '  JTOT-1=', i3, &
             '  JTOT-2=', i4,'  JTOT-D=', i3)
  else
    if (iprint) write (6, 105) ienerg, xmu, ener, jlpar, &
                      (jfirst+0.5), (jfinal+0.5), jtotd
    write (FUNIT_XSC, 105) ienerg, xmu, ener, jlpar, (jfirst+0.5), &
                   (jfinal+0.5), jtotd
105     format ('% ** INTEGRAL CROSS SECTIONS; IEN=', i2,' **', &
            /'%    RMU=', f9.4, '  E=', f9.2,'  JLPAR=', i2, &
             '  JTOT-1=', f5.1, &
             '  JTOT-2=', f6.1,'  JTOT-D=', i3)
  end if
else
  if (.not. flaghf) then
    if (iprint) write (6, 110) ienerg, xmu, ener, numin, &
                               numax, nud
    write (FUNIT_XSC, 110) ienerg, xmu, ener, numin, numax, nud
110     format ('% ** SUMMED DEGENERACY AVERAGED TRANSITION', &
             ' PROBABILITIES;  IEN=', i2,' **', &
            /'    RMU=', f9.4, '  E=', f8.2,'  M-MIN=', i3, &
             '  M-MAX=', i4, '  M-STEP=', i2)
  else
    if (iprint) write (6, 115) ienerg, xmu, ener, (numin+0.5), &
                               (numax+0.5), nud
    write (FUNIT_XSC, 115) ienerg, xmu, ener, (numin+0.5), &
                   (numax+0.5), nud
115     format ('% ** SUMMED DEGENERACY AVERAGED TRANSITION', &
             ' PROBABILITIES;  IEN=', i2,' **', &
            /'    RMU=', f9.4, '  E=', f8.2,'  M-MIN=', f5.1, &
             '  M-MAX=', f6.1, '  M-STEP=', i2)
  end if
end if
if (.not. csflag) then
  if (jlpar .eq. 0) then
    write (FUNIT_XSC, 120)
    if (iprint) write (6, 120)
120     format ('$ ** CC CALCULATION, BOTH PARITIES **')
  else
    write (FUNIT_XSC, 125) jlpar
    if (iprint) write (6, 125) jlpar
125     format ('% ** CC CALCULATION, JLPAR=', i2, ' **')
  end if
else
  if (.not. flaghf) then
    write (FUNIT_XSC, 130) numin, numax
    if (iprint) write (6, 130) numin, numax, nud
130     format ('% ** CS CALCULATION, NUMIN=', i2,', NUMAX=', &
               i2,', NUD=',i2,' **')
  else
    write (FUNIT_XSC, 140) numin + 0.5, numax + 0.5
    if (iprint) write (6, 140) numin+0.5, numax+ 0.5,nud
140     format ('% ** CS CALCULATION, NUMIN=', f4.1, ', NUMAX=', &
            f4.1,', NUD=',i2,' **')
  end if
end if
!mha
! print out of level list only for iprint .ge. 1
if (eprint) then
  if (.not. twomol) then
    if (.not. flagsu) then
      if (iprint) write (6, 145)
      write (FUNIT_XSC, 145)
145       format &
       (/'% LEVEL LIST FOR INTEGRAL CROSS SECTIONS', &
        /'   N   J  INDEX  EINT(CM-1)',/)
    else
      if (iprint) write (6, 150)
      write (FUNIT_XSC, 150)
150       format &
       (/'% LEVEL LIST FOR DEGENERACY AVERAGED', &
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
165         format (1x,'%',i4, f5.1, i6, f11.3)
      end if
170     continue
  else
    if (iprint) write (6, 175)
    write (FUNIT_XSC, 175)
175       format (/'% LEVEL LIST FOR INTEGRAL CROSS SECTIONS', &
              /'%   N   J1   J2  INDEX  EINT(CM-1)'/)
    do 190  i = 1, nlevop
      jj2 = mod( jlev(i), 10)
      jj1 = jlev(i) / 10
      if (iprint) &
      write (6, 180) i, jj1, jj2, inlev(i), elev(i) * econv
      write (FUNIT_XSC, 180) i, jj1, jj2, inlev(i), elev(i) * econv
180       format ('%',i4, 2i5, i6, f11.3)
190     continue
  end if
endif
!mha
!  now sum and average over positive and negative values of index
!     check that number of levels with index negative equals number of levels
!     with index positive, if not abort
!     this test works by (1) checking that no index equals zero
!     and, if so, (2) adding all the values of the index.  this sum
!     should equal zero if there are as many levels with index negative
!     as index positive
if (iaver .gt. 0) then
!ABER
!mha
  if (iaver.gt.1.d0.and.ihomo) then
!mha
    write(6,194)
194     format(' *** HOMONUCLEAR MOLECULE : AVERAGING', &
            ' MAKES NO SENSE ! ***')
    iaver=1
  end if
!ABER
  isum = 0
  do  200 i = 1, nlevop
    if (inlev(i) .eq. 0) then
      write (6, 195) i
      write (6, 195) i
195       format(/' *** INLEV(',i4,')=0;', &
               ' AVERAGING MAY NOT WORK ***')
    else if (inlev(i) .ne. 0) then
      isum = isum + inlev(i)
    end if
200   continue
  if (isum .ne. 0 .and. .not.flaghf) then
    write (6, 230) isum
    write (FUNIT_XSC, 230) isum
230     format ('% *** SUM OF INDICES =',i4, &
               ' AVERAGING MAY NOT WORK ***')
  end if
  if  (iaver .eq. 2) then
    if (iprint) write (6, 245)
    write (FUNIT_XSC, 245)
245     format &
      ('% ** CROSS SECTIONS SUMMED AND AVERAGED OVER INDEX **')
    call aver2 (zmat, scmat, nlevop)
  else if (iaver .eq. 1) then
    if (iprint) write (6, 250)
    write (FUNIT_XSC, 250)
250     format &
      ('% ** CROSS SECTIONS SUMMED FINAL STATE INDEX **')
  end if
end if
!mha
! 4 lines moved from here 2/27/92
!mha
!  find all rows of cross sections matrix for which initial rotational
!  quantum number is equal to one of the values of jout
allocate(ipoint(nlevop))
isize = 0
insize=0
nout = abs (nnout)
do  280  iout = 1, nout
  jtemp = jout(iout)
  do 270 n = 1, nlevop
    if (jlev(n) .eq. jtemp) then
      isize = isize + 1
      jpoint(isize) = n
!mha
      ipoint(isize) = n
!mha
    end if
270   continue
280 continue
!mha
if (niout .gt. 0) then
  nout=abs(niout)
  do  282  iout = 1, nout
    indtemp=indout(iout)
    do 281 n = 1, isize
      if (inlev(jpoint(n)) .eq. indtemp) then
        insize = insize + 1
        ipoint(insize) = jpoint(n)
      endif
281     continue
282   continue
  isize=insize
endif
!mha
!  isize is the number of cross sections to be printed
if (isize .eq. 0) then
!  here if no initial states found
  write (6, 283)
  write (FUNIT_XSC, 283)
283   format ('% ** NO INITIAL STATES FOUND; ABORT')
else
  write (FUNIT_XSC, 290) xthresh
  if (iprint) write (6, 290) xthresh
290   format (/'% ** COLUMN HEADINGS ARE INITIAL STATES, ROW', &
        ' HEADINGS ARE FINAL STATES **', &
      /'%      CROSS SECTION PRINT THRESHOLD=',1pd8.1)
  if (iener.lt.10) then
      write (FUNIT_XSC,295) iener
  elseif (iener.lt.100) then
      write (FUNIT_XSC,296) iener
  elseif (iener.lt.1000) then
      write (FUNIT_XSC,297) iener
  elseif (iener.lt.10000) then
      write (FUNIT_XSC,298) iener
  elseif (iener.lt.100000) then
      write (FUNIT_XSC,299) iener
  endif
295   format('x',i1,'=[')
296   format('x',i2,'=[')
297   format('x',i3,'=[')
298   format('x',i4,'=[')
299   format('x',i5,'=[')
!  now print out desired columns of cross section matrix
  call xscpr2(zmat, xthresh, nlevop, isize, iaver, iprint, isa, FUNIT_XSC, ipoint)
  write (FUNIT_XSC,300)
300   format('];')
endif
deallocate(ipoint)
close (1)
close (3)
return
end
subroutine aver2 (zmat, scmat, n)
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
                   scmat(ind, jnd) ) * 0.5d0
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
! ------------------------------------------------------
subroutine xscpr2 (zmat, xthresh, nlevop, isize, iaver, &
                   iprint, isa, xsc_funit, ipoint)
use mod_cojhld, only: jlev => jhold ! jlev(1)
use mod_coisc1, only: inlev => isc1 ! inlev(1)
use mod_par, only: flaghf, ihomo, ipos
use mod_himatrix, only: transp
implicit double precision (a-h,o-z)
integer, intent(in) :: xsc_funit  ! the file unit for xsc file (it's expected to be open in write mode)
integer, intent(in) :: ipoint(*)

!  current revision date:  10-oct-2001 by ab
!  subroutine to print out specified columns of cross section matrix
!  if iaver = 1, then nth and (n+1) st rows are added before printing
integer i, isize, iskip, j, jcol, jhigh, jj, jlow, jmax, &
        jrow, ncol, nlevop, iaver
integer ind, isa
logical iprint
dimension zmat(nlevop,nlevop), ind(50)
!  first transpose cross section matrix so that initial states are
!  columns and final states are rows
call transp (zmat, nlevop, nlevop)
!  if 132 line printer, then 13 columns of matrix are printed simultaneously
!  if  80 line printer, then  7 columns of matrix are printed simultaneously
iskip = 13
if (.not. ipos) iskip = 7
!  jamx is the total number of the groups of columns to be printed
if (mod(isize, iskip) .eq. 0) then
  jmax = isize / iskip
else
  jmax = isize / iskip + 1
end if
jlow = 1
jhigh = min (iskip, isize)
!  loop ovr columns by groups of 6 or 10
do  150   j = 1, jmax
!  ncol is the number of columns to be printed
  ncol = jhigh - jlow + 1
  do  10   jj = 1, ncol
!  the array ind contains the index of each column which will be printed
    ind(jj) = ipoint(jj - 1 + jlow)
10   continue
!  write as a heading the column index of each column to be printed
  if (.not. flaghf) then
    if (iprint) write (6, 15) ( jlev(ind(i)),i = 1,ncol)
    write (xsc_funit, 15) ( jlev(ind(i)), i = 1,ncol)
15     format (/11x,'J= ', i4, 2x, 12 (2x, i5, 2x) )
  else
    if (iprint) write (6, 30) ( jlev(ind(i))+0.5, &
                                 i = 1,ncol)
    write (xsc_funit, 30) ( jlev(ind(i))+0.5, i = 1,ncol)
30     format (/'%',11x,'J= ', f4.1, 2x, 12 (2x, f5.1, 2x) )
  end if
  if (iprint) write (6, 40) ( inlev(ind(i)), i = 1,ncol)
  write (xsc_funit, 40) ( inlev(ind(i)), i = 1,ncol)
40   format ('%   J    I | I=', i4, 2x, 12 (2x, i5, 2x))
  if (iprint) write (6, 50)
  write (xsc_funit, 50)
50   format (1h )
!  now loop through the rows of the matrix, which will be printed out
!  in groups of 10 with a blank line in between
!  lopp over the rows of the matrix which will be printed
  nlev = 0
!  nlev is the n quantum number for sigma states
  do 95 jrow = 1 , nlevop
!  inrow holds additional quantum index for this row
    inrow = inlev(jrow)
!ABER  summing will be done for HETERONUCLEAR case (PI and SIGMA)
!ABER  and HOMONUCLEAR case (only SIGMA)
    if (iaver.eq.1.and. &
        (.not.ihomo.or.aint(iabs(inrow)/100.0).eq.3)) then
!ABER  find beginning of a new vibrational branch for SIGMA ( ==> nlev=0 )
      if (aint(iabs(inrow)/100.0).eq.3) then
        if (abs(inlev(jrow))-abs(inlev(jrow-1)).ne.0) nlev = 0
!ABER  initialize SWITCH parameter depending on symmetry (ISA) and "parity"
!ABER  of row index (jrow)
        if (nlev.eq.0) then
          if (isa.eq.1.or.isa.eq.0) then
            if (mod(jrow,2).eq.0) iswtch=0
            if (mod(jrow,2).ne.0) iswtch=1
          else if (isa.eq.-1) then
            if (mod(jrow,2).eq.0) iswtch=1
            if (mod(jrow,2).ne.0) iswtch=0
          end if
          if (isa.ne.-1) go to 94
        end if
      end if
!ABER  if sum over final states desired, then sum and skip rows (depending on
!ABER  SWITCH parameter)
      if (iswtch .eq. 0) then
        if (mod (jrow,2) .eq. 0) go to 95
      else if (iswtch .eq. 1) then
        if (mod (jrow,2) .ne. 0) go to 95
      end if
!ABER  now write out the row index followed by the desired matrix elements
!ABER  here for PI states
94       if (aint(iabs(inrow)/100.0).ne.3) then
        if (.not. flaghf) then
          if (iprint) &
          write (6, 60) jlev(jrow),inrow, &
            (zmat(jrow,ind(jcol)) + zmat(jrow + 1,ind(jcol)), &
                jcol = 1,ncol)
          write (xsc_funit, 60) jlev(jrow),inrow, &
            (zmat(jrow,ind(jcol)) + zmat(jrow + 1,ind(jcol)), &
                jcol = 1,ncol)
60           format (i5, i5, 2x, 13 (1pe9.2,1x) )
        else
          if (iprint) &
            write (6, 70) jlev(jrow)+0.5, inrow, &
            (zmat(jrow,ind(jcol)) + zmat(jrow + 1,ind(jcol)), &
                jcol = 1,ncol)
          write (xsc_funit, 70) jlev(jrow)+0.5, inrow, &
            (zmat(jrow,ind(jcol)) + zmat(jrow + 1,ind(jcol)), &
                jcol = 1,ncol)
70           format (f5.1, i5, 2x, 13 (1pe9.2,1x) )
        end if
!ABER  here for SIGMA states
      else if (aint(iabs(inrow)/100.0).eq.3) then
        if (nlev.eq.0.and.isa.ne.-1) then
          if (iprint) &
            write (6, 60) jlev(jrow),inrow, &
            ( zmat(jrow,ind(jcol)), jcol = 1,ncol)
          write (xsc_funit, 60) jlev(jrow),inrow, &
           ( zmat(jrow,ind(jcol)), jcol = 1,ncol)
          nlev = 1
        else
          if (iprint) &
          write (6, 60) jlev(jrow)+1,inrow, &
            (zmat(jrow,ind(jcol)) + zmat(jrow + 1,ind(jcol)), &
                jcol = 1,ncol)
          write (xsc_funit, 60) jlev(jrow)+1,inrow, &
            (zmat(jrow,ind(jcol)) + zmat(jrow + 1,ind(jcol)), &
                jcol = 1,ncol)
          nlev = nlev + 1
        end if
      end if
!ABER  end of summing over final states encountered
    else
      zmax=0d0
      do 75 jcol=1,ncol
         zz=zmat(jrow,ind(jcol))
         if (zz.gt.zmax) zmax=zz
75       continue
      if (zmax.gt.xthresh) then
         if (.not. flaghf) then
           if (iprint) &
             write (6, 60) jlev(jrow),inrow, &
               ( zmat(jrow,ind(jcol)), jcol = 1,ncol)
           write (xsc_funit, 60) jlev(jrow),inrow, &
             ( zmat(jrow,ind(jcol)), jcol = 1,ncol)
         else
           if (iprint) &
             write (6, 70) jlev(jrow)+0.5, inrow, &
               ( zmat(jrow,ind(jcol)), jcol = 1,ncol)
           write (xsc_funit, 70) jlev(jrow)+0.5, inrow, &
             ( zmat(jrow,ind(jcol)), jcol = 1,ncol)
         end if
      endif
    end if
    jj = jrow / 10
    if (10 * jj .eq. jrow .and. jrow .ne. nlevop) then
      if (iprint) write (6, 80)
    end if
80     format(1h )
95   continue
  if (jhigh .eq. nlevop) go to 160
  jlow = jhigh + 1
  jhigh = min ( (jlow + iskip - 1), isize)
150 continue
160 return
end
end module mod_hibrid2
