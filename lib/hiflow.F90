#include "assert.h"
!  -------------------------------------------------------------
module mod_flow
contains
subroutine flow (z, w, zmat, amat, bmat, jq, lq, inq, jlev, &
            elev, inlev, isc1, isc2, isc3, isc4, lsc1, &
            sc2, sc1, sc3, sc4, &
            sc5, sc6, sc7, sc8, sc9, tq1, tq2, tq3, men, &
            nmax, nairy)
!  -------------------------------------------------------------
!  program to control the log-derivative/airy integration
!  written by:  millard alexander
!  additions by: b. follmeg, h-j werner
!  current revision date:  1-oct-2001 by mha
!  -------------------------------------------------------------

!  variables in common block /copmat/
!    rtmn,rtmx: minimum and maximum turning points
!    iflag:     variable used in determination of turning points (not used her
!           iflag = 0 if all channels are in classically forbidden region
!           iflag = 1 if some channels are open
!           iflag = 2 if all asymptotically open channels are open at r
!  variable in common block /cojlpo/
!    jlpold:      parity used in xwrite subroutine to insure correct
!                 accumulation of partial waves in cases where jlpar=0
!
!  variables in module constants
!    econv:       conversion factor from cm-1 to hartree
!    xmconv:      conversion factor from amu to atomic units
!
!  variable in common block /coopti/
!    optifl:      flag, signals if the calculation is an optimization
!
use mod_cosout, only: nnout, jout
use mod_cocent, only: cent
use mod_coeint, only: eint
use mod_coener, only: energ
use mod_hibasis, only : basis
use mod_version, only : version, acknow
use mod_hibrid5, only : soutpt, nusum, xwrite, wrhead, restrt, rsave
use mod_ancou, only: ancou_type
use constants
use mod_hibrid2, only: set_default_params
use mod_hibrid3, only: propag
use mod_par, only: airyfl, prairy, bastst, chlist, &
                csflag, flaghf, flagsu, ihomo, ipos, logdfl, &
                prlogd, noprin, prpart, rsflag, prsmat, &
                t2test, prt2, twomol, wrsmat, wrpart, wrxsec, &
                prxsec, nucros, photof, wavefl, boundc, &
                jtot1, jtot2, jtotd, jlpar, nerg, numax, numin, nud, &
                lscreen, iprint, &
                fstfac=>scat_fstfac, rincr=>scat_rincr, rcut=>scat_rcut, rendai=>scat_rendai, rendld=>scat_rendld, rstart=>scat_rstart, spac=>scat_spac, tolhi=>scat_tolai, xmu
use funit
use ipar_enum
use rpar_enum
use mod_hinput, only:hinput
use mod_parpot, only: potnam=>pot_name, label=>pot_label
use mod_selb, only: ibasty
use mod_ered, only: ered, rmu
use mod_phot, only: phot_photof => photof, wavefn, boundf, writs
use mod_surf, only: surf_flagsu => flagsu
use mod_sav, only: iipar, ixpar, irpar, rxpar
implicit none
real(8), intent(out) :: z(nmax,nmax)
real(8), intent(out) :: w(nmax,nmax)
real(8), intent(out) :: zmat(nmax,nmax)
real(8), intent(out) :: amat(nmax,nmax)
real(8), intent(out) :: bmat(nairy,nairy)
integer, intent(out) :: jq(10)
integer, intent(out) :: lq(10)
integer, intent(out) :: inq(10)
integer, intent(out) :: jlev(1)
real(8), intent(out) :: elev(1)
integer, intent(out) :: inlev(1)
integer, intent(out) :: isc1(9)
integer, intent(out) :: isc2(1)
integer, intent(out) :: isc3(1)
integer, intent(out) :: isc4(1)
logical, intent(out) :: lsc1(5)
real(8), intent(out) :: sc2(1)
real(8), intent(out) :: sc1(2)
real(8), intent(out) :: sc3(1)
real(8), intent(out) :: sc4(1)
real(8), intent(out) :: sc5(1)
real(8), intent(out) :: sc6(1)
real(8), intent(out) :: sc7(1)
real(8), intent(out) :: sc8(1)
real(8), intent(out) :: sc9(1)
real(8), intent(out) :: tq1(1)
real(8), intent(out) :: tq2(1)
real(8), intent(out) :: tq3(1)
integer, intent(in) :: men
integer, intent(in) :: nmax
integer, intent(in) :: nairy

type(ancou_type), allocatable :: v2
#if defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IRIS)
real(8) secnds
common /codec/ ttim(2)
real(8) :: ttim
#endif

common /cputim/ cpuld,cpuai,cpupot,cpusmt,cpupht
real(8) :: cpuld, cpuai, cpupot, cpusmt, cpupht
common /copmat/ rtmn, rtmx, iflag
real(8) :: rtmn, rtmx
integer :: iflag

common /cojlpo/ jlpold
integer :: jlpold

common /coopti/ optifl
logical :: optifl

common /constp/ nsteps, isteps
integer :: nsteps, isteps

integer :: nlev(25)

integer :: jtotmx
character*20 :: cdate
character*10 :: time
character*10 :: timew, cpubaw, cpuptw, cpuaiw, cpuldw, cpusmw, cpuouw, &
             cpuphw, timew1, timew2, time1, time2
logical :: clist, firstj, ready
!  -------------------------------------------------------------
logical :: first

logical :: twojlp
data twojlp / .false. /

real(8), parameter :: mtime_granularity = 0.5d0
integer, parameter :: tmp_file = 1
!   to obtain timing information calls are made to system-specific
!   time subroutine:  call mtime(cpu,wallt) where cpu is clock cpu
!   time in seconds and wallt is wall clock time in seconds
!   user will have to change this subroutine for his own installation
!  subroutine to get input data

real(8) :: cpubas, cpuout, cpupt, dinsid, dlogd, ener, eshift
integer :: i, ien, ierr, ifile, ii, irec
integer :: jfirst, jfrest, jj, jlprsv, jtop, jtopo, jtot, jtoto
integer :: nch, nchmax, nchop, nchtop, nfile, nlevel, nlevop, nopen, ntop, nu, nufirs, nulast, numj, nutop
real(8) :: rstrt0, rtmn1, rtmnla, rtmx1, rtmxla
real(8) :: t1, t11, t2, t22, tb, tbm, tcpu0, tcpu1, tcpu2, tcpuf, twall0, twall1, twall2, twallf
real(8) :: xjtot

real(8) :: second


first=.true.
!  get default data
call set_default_params
1 call hinput(first)
cpupt=0
#if defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IRIS)
ttim(1)=0.d0
ttim(2)=secnds(0.0)
#endif
call mtime (tcpu0, twall0)
phot_photof = photof
wavefn = wavefl
boundf = boundc
writs = wrsmat
surf_flagsu = flagsu
!  subroutine to open required i/o files
if (.not.bastst) then
      call openfi (nerg)
end if
call version(9)
call acknow(9,ipos)
call dater(cdate)
write (9, 15) label
15 format (/' ** LABEL:     ', (a))
write (9, 16) potnam
16 format (' ** POT NAME:  ', (a))
write (9, 20) cdate
20 format ( ' ** DATE:      ', (a))
if (airyfl .and. (nairy .ne. nmax)) then
  write (9, 35)
  write (6, 35)
35   format( &
     /' *** NAIRY .NE. NMAX FOR AIRY INTEGRATION; ABORT ***')
  call exit
end if
if (.not. airyfl .and. (nairy .ne. 1)) write (9, 45)
45 format( &
  /' *** WARNING:  NAIRY .NE. 1 BUT AIRYFL .EQ. .FALSE.')
!  check for proper choice of output files and dimension limits
if (nerg .gt. men) then
  write (9, 50)
  write (6, 50)
50   format (/' *** TOO MANY TOTAL ENERGIES; ABORT ***')
  call exit
end if
if (nerg.gt.19) then
   if (wrsmat.and.(prxsec.or.wrxsec))  then
      write (6,55) nerg
55       format( &
    /' NERG = ',i2,' > 19', &
    ' FOR WRITS .AND. (WRXSEC .OR. WRPART) ALL TRUE; ABORT ***')
      call exit
   endif
endif
write (9, 60) nmax,nairy
60 format (' ** NMAX=',i4,'  NAIRY=',i4)
!  convert collision energy, rotational constant, and reduced mass
!  to atomic units
rmu = xmu / xmconv
jtotmx = 0
!  reduce maximum jtot if all channels are closed (not for bound states)
if (.not.boundc) then
  jtop = sqrt (energ(1) / econv * rcut * rcut * 2.d0 * rmu ) &
        - 0.5d0
  numj = (jtot2 - jtot1) / jtotd + 1
  jtotmx = jtot1 - jtotd
  do 75  jj = 1, numj
    jtotmx = jtotmx + jtotd
    if (jtotmx .ge. jtop) then
      jtotmx = jtotmx - jtotd
      go to 76
    end if
75   continue
endif
76 firstj=.true.
jfirst=jtot1
jtop =jtotmx
cpubas=0
cpuout=0
cpuai=0
cpuld=0
cpupot=0
cpupht=0
cpusmt=0
jtotd=ixpar(IPAR_JTOTD)
jlpar=ixpar(IPAR_JLPAR)
nerg =ixpar(IPAR_NERG)
numax=ixpar(IPAR_NUMAX)
numin=ixpar(IPAR_NUMIN)
nud=  ixpar(IPAR_NUD)
nufirs=numin
nulast=numax
nutop=numax
jlprsv=jlpar
fstfac=rxpar(RPAR_SCAT_FSTFAC)
rincr= rxpar(RPAR_SCAT_RINCR)
rcut=  rxpar(RPAR_SCAT_RCUT)
rendai=rxpar(RPAR_SCAT_RENDAI)
rendld=rxpar(RPAR_SCAT_RENDLD)
rstrt0=rxpar(RPAR_SCAT_RSTART)
spac=  rxpar(RPAR_SCAT_SPAC)
tolhi= rxpar(RPAR_SCAT_TOLAI)
isteps=0
dlogd = rendld - rstart
xmu=rxpar(RPAR_XMU)
rtmnla=rstrt0
dinsid=0
twojlp=jlpar.eq.0.and..not.csflag
if (twojlp) jlpar = 1
jlpold = jlpar
jfrest=0
74 jfirst=ixpar(IPAR_JTOT1)
if(jfrest.gt.0) jfirst=jfrest
jtotmx=jtop
if(nucros) nulast=numin
if (jtot1.eq.jfirst .or. jtot1.eq.jfirst+1) rstart=rxpar(RPAR_SCAT_RSTART)
80 jtot1 = jfirst
if (.not.boundc) jtot2 = jtotmx
nchmax = 0
!  this is beginning of loop over total angular momentum
!  (partial wave index)

100 jtot = jtot1
!  get restart parameters if run is to be restarted
ready=.false.
if (rsflag) then
!  read in restart information
  call restrt (jtot,jtopo,jtotd,jlpar,nu,nutop,nud,nerg,nlev, &
         nchmax,rtmn1,rtmx1,dinsid,wrsmat,csflag,nucros)
  rtmnla=rtmn1
  jlpold = jlpar
!  move s-matrix files to last partial wave done
  if(wrsmat) then
     do 88 ifile=1,nerg
        nfile = FUNIT_SMT_START + ifile - 1
        ered = energ(ifile)/econv
        nlevop=nlev(ifile)
        call wrhead(nfile, cdate, &
                  sc1(1),  sc1(2), lsc1(1), lsc1(2), lsc1(3), &
                 lsc1(4), lsc1(5), isc1(1), isc1(2), isc1(3), &
                 isc1(4), isc1(5), isc1(6), isc1(7), isc1(8), &
                 isc1(9), isc2, isc3, sc2, isc4)
        call fimovs(nfile,jtot,jlpar,nu,ifile,ierr)
        if(ierr.ne.0) then
           write(6,87) jtot,jlpar
           write(9,87) jtot,jlpar
87            format(/'  *** JTOT=',i5,', JLPAR=',i3,' NOT FOUND', &
                   ' IN S-MATRIX FILE, RESTART IMPOSSIBLE')
           call exit
        end if
88      continue
  end if
  firstj=.false.
  rsflag=.false.
  if(.not. csflag) then
    jtoto=jtot
    jtot=jtot+jtotd
    if(jtot.gt.jtot2) then
      if(jlpar.lt.0.or.jlprsv.gt.0) then
        write(6,89)
        write(9,89)
89         format(/'  *** RESTART, BUT NO MORE PARTIAL WAVES', &
                ' REQUESTED ')
        ready=.true.
        goto 105
      end if
      jlpar=-1
      jtot=ixpar(IPAR_JTOT1)
    else if(jlprsv.eq.0.and.jlpold.eq.-1.and.jtoto.eq.jtopo) then
      twojlp=.true.
      jlpar=1
      jlpold=1
      jfrest=jtot
    end if
    write (6, 90) jtot,jlpar
    write (9, 90) jtot,jlpar
90     format (' ** CONTINUE CC CALCULATION AT JTOT=',i5, &
            ', JLPAR=',i5)
  else
    if(nucros) then
      nu=nu+nud
      if(nu.gt.numax) then
        write(6,89)
        write(9,89)
        ready=.true.
        goto 105
      end if
      numin=nu
      write(6,92) nu
      write(9,92) nu
92       format(/' ** CONTINUE CS CALCULATION AT NU=',i5)
      goto 74
    else
      nu=ixpar(IPAR_NUMIN)
      jtot=jtot+jtotd
      if(jtot.gt.jtot2) then
        write(6,89)
        write(9,89)
        ready=.true.
        goto 105
      end if
      write(6,93) jtot
      write(9,93) jtot
93       format (' ** CONTINUE CS CALCULATION AT LBAR=',i5)
    end if
  end if
  jtot1=jtot
end if
!
!      write (9, 95)
!95    format (1h ,79('='))
xjtot = jtot
if (flaghf) xjtot = xjtot + 0.5d0
105 continue
rtmn1 = max (rendld, rendai)
rtmx1 = 0.
!  this is beginning of loop over coupled states projection index
!  in the case of cc calculation this loop is executed only once
nu = numin -  1
110 nu = nu + 1
rtmx = 0.
rtmn = max (rendld, rendai)
ien = 1
!  this is beginning of loop over collision energies
115 ener = energ(ien)
ered = ener/econv
if(.not.ready) then
if (csflag) then
  if (.not. flaghf) then
    if(prpart.and..not.noprin) then
      write (9, 120) jtot, nu, ien, ener
      write (6, 120) jtot, nu, ien, ener
120       format (/' ** LBAR=',i4,'  NU=',i2,'  IEN =',i3, &
             '  ENERGY (CM-1) =', f11.4)
    endif
  else
    if(prpart.and..not.noprin) then
      write (9, 125) jtot,nu+0.5d0, &
              ien,ener
      write (6, 125) jtot, nu+0.5d0, ien, ener
125       format (/' ** LBAR=',i4, '  NU=', f5.1, '  IEN =',i3, &
             '  ENERGY (CM-1) =', f11.4)
    endif
  end if
else
  if (.not. flaghf) then
    if(prpart.and..not.noprin) then
      write (9, 130) jtot,jlpar,ien,ener
      write (6, 130) jtot, jlpar, ien, ener
130       format (/' ** JTOT =',i4,'  JLPAR =', i2,'  IEN =',i3, &
             '  ENERGY (cm-1) =', f11.4)
    endif
  else
    if(prpart.and..not.noprin) then
      write (9, 135) jtot+0.5d0, &
      jlpar, ien, ener
      write (6, 135) jtot+0.5d0, jlpar, ien, ener
135       format (/' ** JTOT =',f7.1, '  JLPAR =', i2,'  IEN =',i3, &
             '  ENERGY (cm-1) =', f11.4)
    endif
  end if
end if
end if
!  ered is collision energy in hartree
eshift = 2.d0 * rmu * (ered - energ(1) / econv)
!  if first energy (ien = 1), then set up channel indices, array of
!  internal energies, and angular coupling matrices for each vlambda(r)
call mtime (t1,t2)
if (ien .eq. 1) then
  clist = .false.
  if (chlist .and.(jtot .eq. jfirst)) clist = .true.
  if (noprin) clist = .false.


!*        write(6,5554) nch, nmax, nchtop
!*5554    format('BEFORE BASIS:',3i6)
#ifdef DEBUG
#define ENSURE_BASIS_SCRATCHS_ARE_REAL_SCRATCHS
#endif

#ifdef ENSURE_BASIS_SCRATCHS_ARE_REAL_SCRATCHS
#define DUMMY_REAL_VALUE 42.0
  do i = 1, nmax
    sc1(i) = DUMMY_VALUE
    sc2(i) = DUMMY_VALUE
    sc3(i) = DUMMY_VALUE
    sc4(i) = DUMMY_VALUE
  end do
#endif

  call basis (jq, lq, inq, jlev, elev, inlev, nlevel, nlevop, &
              sc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu, &
              csflag, clist, bastst, ihomo, nu, numin, jlpar, &
              twomol, nch, nmax, nchtop, v2)

  ! if nch == 0, then v2 is usually not allocated at all
  if ( nch > 0 ) then
    ASSERT(allocated(v2))  ! if this fails, this means that the used base doesn't yet support v2 as growable array
  end if

#ifdef ENSURE_BASIS_SCRATCHS_ARE_REAL_SCRATCHS
  do i = 1, nmax
    sc1(i) = DUMMY_VALUE
    sc2(i) = DUMMY_VALUE
    sc3(i) = DUMMY_VALUE
    sc4(i) = DUMMY_VALUE
  end do
#endif

!*        write(6,5556) nch, nmax, nchtop
!*5556    format('AFTER BASIS: ',3i6)


  if (ready) goto 370
  if (bastst) then
     write(6,140)
     write(9,140)
140      format(' ** BASTST=.TRUE.; TEST OF SUBROUTINE BASIS')
     goto 420
  end if
!  on return from basis:
!  nch is number of channels
!  nchtop is the maximum row dimension of all matrices passed to the
!  subroutines propag and soutpt
!  inq is the additional quantum index of each channel
!  jq is array of rotational quantum numbers
!  lq is array of orbital angular momenta
  if (nch .gt. nchmax) nchmax = nch
! in rare cases one might come back from basis with no open channels (but
! with closed channels present).  to deal with this case set nchop to the
! number of open channels
  nchop = 0
  if (nch .ge. 1) then
    do 144 ii = 1 , nch
      if (ered - eint(ii) .gt. 0.) nchop = nchop + 1
144     continue
  end if
  if (nchop .eq. 0) then
    if (csflag) then
      if (nu .eq. numin) then
        write (6, 145)
        write (9, 145)
145         format( &
           /' *** NCH = 0 FOR NU=NUMIN IN CS CALCULATION;', &
            ' END CALCULATION')
!  reset jtot2 to reflect absence of all channels in cs calculations
        jtot2 = jtot - jtotd
        write (6, 150) jtot2
        write (9, 150) jtot2
150         format (/' ** RESET JTOT2=', i4)
        go to 370
      else
        if (nu - 1 .lt. nutop) then
          nutop = nu - 1
          nulast=nutop
          write (6, 155) nu, nutop
155           format &
           (' ** NCH = 0 FOR NU=',i3,' IN CS CALCULATION;', &
            ' SET NUTOP=', i3, ' AT NEXT JTOT'/, &
            '    MOVE ON TO NEXT PARTIAL WAVE')
        go to 300
        end if
      end if
    else
!  here for cc calculation, if no channels and jtot .le. 2 move on to
!  next partial wave, otherwise reset jtot2 and quit calculation
      if (jtot .le. 2) then
        write (9, 160)
        write (6, 160)
160         format (/' ** NCH = 0, MOVE ON TO NEXT PARTIAL WAVE')
!--------------------------------------------------------------------
! restart with next jtot
         jfirst = jfirst + 1
         if(jtotd.gt.1) jtotmx = jtotmx + 1
         goto 80
!--------------------------------------------------------------------
      else
        write (6, 165)
        write (9, 165)
165         format (/' *** NCH = 0 AND JTOT2.GE.2; END CALCULATION')
        jtot2 = jtot2 - jtotd
        write (6, 150) jtot2
        write (9, 150) jtot2
! Claire's modification: if nch=0 then eventually go to next parity:
        if (twojlp .and. jlpar .gt. 0) then
          jlpar = -1
          goto 74
        else
          go to 370
        end if
      end if
    end if
  end if
!  store channel parameters on unit FUNIT_CHANNEL_PARAMS if this calculation is to be
!  performed at a second energy
  if (nerg .gt. 1) then
! open file for storage of transformation matrices
    rewind (FUNIT_CHANNEL_PARAMS)
    write (FUNIT_CHANNEL_PARAMS, 170) nch
170     format (i4)
    if (nch .gt. 0) then
      do 180  i = 1, nch
        write (FUNIT_CHANNEL_PARAMS, 175) jq(i), lq(i), inq(i), cent(i), eint(i)
175         format (3i6, 2e25.15)
180       continue
    end if
  end if
else if(ien.gt.1) then
  rewind (FUNIT_CHANNEL_PARAMS)
  read (FUNIT_CHANNEL_PARAMS, 170) nch
  if (nch .gt. 0) then
    do 250  i = 1, nch
      read (FUNIT_CHANNEL_PARAMS, 175) jq(i), lq(i), inq(i), cent(i), eint(i)
250     continue
  end if
end if
irec=(ien-1)*5+2
nlevop=0
do 255 i=1,nlevel
255 if(elev(i).le.ered) nlevop=nlevop+1
nlev(ien)=nlevop
if (ien .gt. 1) then
  if (nlev(ien) .ne. nlev(ien-1)) then
    write (6, 256) ien, nlev(ien), ien-1, nlev(ien-1)
256     format (' *** NLEV(',i2,') = ',i3,' .NE. NLEV(',i2,') = ',i3, &
           ' ABORT ***')
    call exit
  endif
endif
if (wrxsec .or. prxsec .or. prpart .or. wrpart) then
  if (.not. wavefl .and. .not. photof) then
    call dres(nlevop**2+4,1,irec)
    call dres(nlevop**2+4,1,irec+1)
    call dres(nlevop**2+4,1,irec+2)
    call dres(nlevop**2+4,1,irec+3)
    call dres(nlevop**2+4,1,irec+4)
    if(nucros) then
      irec=(nerg+ien-1)*5+2
      call dres(nlevop**2+4,1,irec)
      call dres(nlevop**2+4,1,irec+1)
      call dres(nlevop**2+4,1,irec+2)
      call dres(nlevop**2+4,1,irec+3)
      call dres(nlevop**2+4,1,irec+4)
    end if
    call dsave(1)
  endif
endif
! store header and reserve space on direct access file 2 if wavefunction
! is desired (not if bound state calculation)
if (wavefl .and. .not.boundc) &
  call wavewr(jtot,jlpar,nu,nch,rstart,rendld)
call mtime (t11, t22)
tb =  t11 - t1
tbm = t22 - t2
cpubas=cpubas + tb
!  at subsequent partial waves, adjust starting point for integration for
!  next partial wave to be dinsid inside of innermost classical turning
!  point but no larger than rendld
!  10/14/08 mha:  but don't allow an increase greater than spac/2 in rstart
if (jtot .gt. jfirst) then
  rstart = max(rtmnla-dinsid,rstrt0)
  rstart = min(rstart,rstrt0+0.5*spac)
  rstrt0 = rstart
!  adjust starting point of airy integration to be no less than rstart + dlogd
!  but no larger than rendai
  rendld = min (rstart + dlogd, rendai)
  if (.not. logdfl) rendld = rstart
end if

!ger (next 2 lines)
write (9,'(/" ** J =",i5," JLPAR =",i2," STARTED")') jtot,jlpar
#if defined(HIB_UNIX_IBM) || defined(HIB_UNIX_AIX)  || defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
!      call flush_(9)
call flush(9)
#endif
#if defined(HIB_UNIX_DEC) || defined(HIB_UNIX_CONVEX) || defined(HIB_UNIX_IRIS)
call flush(9)
#endif
! replace statement below to accommodate q.ma's revised
! version of bound (28-jun-2013, p.dagdigian)
! do not reduce maximum size of matrices if bound state calculation
!      if (boundc) then
!        ntop=nmax
!      else
!        ntop=nchtop
!      endif
ntop = nchtop
!

call propag (z, w, zmat, amat, bmat, &
             jq, lq, inq, &
             ien, nerg, ered, eshift, rstart, rendld, spac, &
             tolhi, rendai, rincr, fstfac, tb, tbm, &
             ipos, prlogd, noprin, airyfl, prairy, &
             nch, nopen, nairy, ntop, v2)

! if bound state calculation, end it now
if (boundc) then
  endfile (9)
  close (9)
  goto 1
endif
!  now print out s-matrix and t-matrix squared, and calculate partial
!  cross sections and print them out, if desired
call soutpt (z, w, zmat, amat, &
             lq, jq, inq, isc1, isc2, bmat, tq1, &
             jlev, elev, inlev, jtot, jfirst, &
             jtot2, jtotd, nu, numin, nulast, nud, jlpar, ien, &
             ipos, csflag, flaghf, prsmat, prt2, t2test, &
             wrsmat, wrpart, prpart, wrxsec, prxsec, twomol, &
             nucros, firstj, nlevel, nlevop, nopen, nchtop, &
             twojlp)
cpuout = cpuout + second() - t11

!  on return from soutpt:
!     if wrxsec,prxsec, wrpart, and prpart are all .false., the upper-left
!     nopen x nopen block of z contains the modulus squared of the t-matrix
!     otherwise, the upper nlevop x nlevop block of z contains the partial
!     cross sections
!     the upper-left nopen x nopen block of w contains the real part of
!     the s-matrix
!     the upper-left nopen x nopen block of zmat contains the imaginary part
!     of the s-matrix
!     the arrays eint, cent, jq, lq, inq have been packed to eliminate
!     the closed-channel components
if (ien.eq.1 .and. nchop.gt.0) then
  rtmn1 = min(rtmn1,rtmn)
  rtmx1 = max(rtmx1,rtmx)
end if
ien = ien + 1
call mtime (tcpuf, twallf)
if (ien.eq.2) then
   tcpu1=tcpuf-tcpu0
   twall1=twallf-twall0
endif
!  go back to start calculation at another energy
if (ien .le. nerg) go to 115
!  first partial wave has been calculated for all energies
if (firstj) firstj = .false.
!  go back to start calculation at another value of coupled-states projection
!  index
300 if (nu .lt. nulast .and. nu .lt. nutop) go to 110
if( .not.nucros) nulast = nutop
!  if first partial wave, then set distance inside turning point at which
!  logd integration starts
if (jtot .eq. jfirst)  then
  dinsid = rtmn1 - rstart
  write (9, 330) dinsid
330   format (/' ** INTEGRATION WILL START', f6.3, &
           ' BOHR INSIDE INNER TURNING POINT')
end if
!  save last min and max turning points for next partial wave
rtmnla = rtmn1
rtmxla = rtmx1
if(prpart.and..not.nucros) write (9, 350)
350 format (1h ,79('='))
if(.not.nucros) then
!.....save restart information
  if (wrxsec .or. prxsec .or. prpart .or. wrpart) then
    if (.not. wavefl .and. .not. photof) then
     call rsave (jtot,jtop,jtotd,jlpar,nu,nutop,nud,nerg,nlev, &
                  nchmax,rtmn1,rtmx1,dinsid,wrsmat,csflag, &
                  nucros)
    endif
  endif
  call mtime (tcpuf, twallf)
  tcpuf = tcpuf - tcpu0
  twallf = twallf - twall0
  if (nerg.gt.1) then
      tcpu2=(tcpuf-tcpu1)/(nerg-1)
      twall2=(twallf-twall1)/(nerg-1)
  endif
  call gettim(tcpuf,time)
  call gettim(tcpu1,time1)
  call gettim(twallf,timew)
  call gettim(twall1,timew1)
  if (nerg.gt.1) then
     call gettim(tcpu2,time2)
     call gettim(twall2,timew2)
  endif
  call gettim(cpubas,cpubaw)
  call gettim(cpuld,cpuldw)
  call gettim(cpuai,cpuaiw)
  call gettim(cpupht,cpuphw)
  call gettim(cpusmt,cpusmw)
  call gettim(cpupot,cpuptw)
  call gettim(cpuout,cpuouw)
  call dater(cdate)
  if (.not. noprin) then
    write (6, 360) jtot,jlpar,cpubaw,cpuptw,cpuldw,cpuaiw,cpuphw, &
       cpusmw,cpuouw,time,timew,cdate
360     format (/' ** J =', i5,' JLPAR =', i2,' COMPLETED'/1x, &
   'CPU-TIMES:', &
   '  BASIS:',a,'  POT:',a,'  LOGD: ',a,'  AIRY: ',a/11x, &
   '  PSI0:',a,'  SMAT:',a,'  SOUT: ',a,'  CUMULATIVE:',a,/11x, &
   '  ELAPSED: ',a,'  CURRENT DATE:  ',a)
    write (6, 361) rtmn, rtmx
361     format (' TURNING POINTS (MIN/MAX) =', 2(f8.2) )
  else
    write (6, 362) jtot, jlpar, time,timew,cdate
    write (9, 362) jtot, jlpar, time,timew,cdate
362     format (' ** J =', i5,' JLPAR =', i2, &
     ' FINISHED; ', &
     ' CPU:',a,'  WALL:',a,'  DATE: ',a)
  endif
#if defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IRIS) || defined(HIB_UNIX_HP)
  call flush(6)
#endif
  cpubas=0
  cpuout=0
  cpuai=0
  cpupht=0
  cpuld=0
  cpupot=0
  cpusmt=0
end if
!  go back to start calculation at the next partial wave
jtot1 = jtot + jtotd
if (jtot1 .le. jtot2) go to 100
if (twojlp .and. jlpar .gt. 0) then
  jlpar = -1
  rstrt0=rxpar(RPAR_SCAT_RSTART)
  rstart=rstrt0
  rtmnla=rstrt0
  dinsid=0
  goto 74
end if
!.....next nu value if nu runs in outer loop
if (nucros .and. .not. wavefl .and. .not. photof) then
  call nusum (z, tq1, tq2, tq3, &
              jlev,elev, inlev, jtot, jfirst, &
              jtop, jtotd, nu, nufirs, numax, nud, jlpar, &
              nerg, ipos, csflag, flaghf, wrpart, prpart, &
              twomol, nucros, nlevel, nlev, nopen, nmax, tmp_file)
!.....save restart information
  if (wrxsec .or. prxsec .or. prpart .or. wrpart) then
    if (.not. wavefl .and. .not. photof) then
      call rsave (jtot,jtop,jtotd,jlpar,nu,nutop,nud,nerg,nlev, &
                  nchmax,rtmn1,rtmx1,dinsid,wrsmat,csflag, &
                  nucros)
    endif
  endif
  call mtime (tcpuf, twallf)
  tcpuf = tcpuf - tcpu0
  twallf = twallf - twall0
  call gettim(tcpuf,time)
  call gettim(twallf,timew)
  call gettim(cpubas,cpubaw)
  call gettim(cpuld,cpuldw)
  call gettim(cpuai,cpuaiw)
  call gettim(cpupht,cpuphw)
  call gettim(cpusmt,cpusmw)
  call gettim(cpupot,cpuptw)
  call gettim(cpuout,cpuouw)
  call dater(cdate)
  if (.not. noprin) then
    write (6, 366) nu,cpubaw,cpuptw,cpuldw,cpuaiw,cpuphw, &
   cpusmw, &
   cpuouw,time,timew,cdate
366     format (/' ** NU =', i3,' COMPLETED'/1x, &
   'CPU-TIMES:', &
 '  BASIS:',a,'  POT: ',a,'  LOGD: ',a,'  AIRY:    ',a/11x, &
 '  PHOTO:',a,'  SMAT: ',a,'  SOUT:',a,'  CUMULATIVE:',a, &
   '  ELAPSED: ',a,'  CURRENT DATE:  ',a)
    write (6, 361) rtmn, rtmx
  else
    write (6, 367) nu, time,timew,cdate
    write (9, 367) nu, time,timew,cdate
367     format (' ** NU =', i2, &
     ' COMPLETED; ', &
      ' CPU:',a,'  ELAPSED:',a,'  CURRENT DATE: ',a)
  endif
#if defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IRIS) || defined(HIB_UNIX_HP)
  call flush(6)
#endif
  cpubas=0
  cpuld=0
  cpuai=0
  cpupht=0
  cpusmt=0
  cpupot=0
  cpuout=0
  if(nulast.lt.numax) then
    numin = numin + nud
    goto 74
  end if
end if
!  calculation has now been done at all partial waves and, if
!  desired, at both values of jlpar
!  write out integral cross sections if desired
370 if (prxsec .or. wrxsec) then
  if(.not.bastst) &
  call xwrite (amat, tq3, jlev, elev, inlev, nerg, energ, &
             jfirst, jtot2, jtotd, csflag, flaghf, &
             wrxsec, prxsec, ipos, twomol, nucros, nlevel, &
             nlev, nufirs, nulast, nud, jlpar, nchtop, nmax, &
             ihomo)
 
endif
if (.not. bastst .and. &
   prxsec .or. wrxsec .or. wrpart .or. prpart) then
  do 400 ien = 1, nerg
    nfile = FUNIT_ICS_START + ien - 1
    close (nfile)
    if (wrpart) then
      nfile = FUNIT_PCS_START + ien - 1
      close (nfile)
      if (csflag .and. &
          (wrpart .or. prpart .or. wrxsec .or. prxsec)) then
        nfile = FUNIT_APCS_START + ien - 1
        close (nfile)
      end if
    end if
400   continue
end if
if (.not. bastst .and. &
   prxsec .or. wrxsec .or. wrpart .or. prpart) then
  if (.not. wavefl .and. .not. photof) then
    call dclos(1)
   endif
endif
if (wavefl .and. .not. boundc) close(FUNIT_WFU)
if (wrsmat) then
  do 410 ien = 1, nerg
  nfile = FUNIT_SMT_START + ien - 1
  call closf (nfile)
  close (nfile)
410   continue
end if
if (nerg .gt. 1) then
  if (airyfl) close (FUNIT_TRANS_MAT)
  close (FUNIT_QUAD_MAT)
end if
420 call dater (cdate)
if (.not. optifl) then
   write (6, 350)
   write (6,500) nchmax
   write (9,500) nchmax
500    format(' **** END OF CALCULATION ****',/, &
    '      MAXIMUM NUMBER OF CHANNELS USED WAS:  ',i4)
   if (nerg.eq.1) then
      twallf = max(twallf, mtime_granularity)
      write (6,505) timew, time, 100d0*(tcpuf/twallf)
      write (9,505) timew, time, 100d0*(tcpuf/twallf)
505       format('      TIMING:  ELAPSED',(a),'/ CPU',(a), &
             '/ MP RATIO',f7.1,' %')
   endif
   if (nerg.gt.1) then
      twall1 = max(twall1, mtime_granularity)
      write(6,510) timew1, time1, 100d0*(tcpu1/twall1)
      write(9,510) timew1, time1, 100d0*(tcpu1/twall1)
510       format('      TIMING (FIRST ENERGY):  ELAPSED',(a), &
            '/ CPU',(a), &
           '/ MP RATIO',f7.1)
      twall2 = max(twall2, mtime_granularity)
      write(6,515) timew2, time2, 100d0*(tcpu2/twall2)
      write(9,515) timew2, time2, 100d0*(tcpu2/twall2)
515       format('      TIMING (SECOND ENERGY): ELAPSED',(a), &
           '/ CPU',(a), &
           '/ MP RATIO',f7.1)
   endif
   write (6,520) cdate
   write (9,520) cdate
520    format('      CURRENT DATE:  ',(a))
   write (6, 350)
   write (9, 350)
#if defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IRIS) || defined(HIB_UNIX_HP)
   call flush(6)
#endif
end if
close (23)
endfile (9)
close (9)
goto 1
end
end module mod_flow
