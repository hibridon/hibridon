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

end module mod_hibrid2
