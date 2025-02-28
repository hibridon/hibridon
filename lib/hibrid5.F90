#include "assert.h"

module mod_savfile
implicit none

  integer, parameter :: num_recs_per_energy = 5

  ! records for each energy
  enum, bind( C )
  enumerator ::  &
    EN_REC_PRESENT_INTEGRAL_XS  =   0, &  ! present integral cross section
    EN_REC_PREVIOUS_INTEGRAL_XS =   1, &  ! previous integral cross section
    EN_REC_PRESENT_PARTIAL_XS   =   2, &  ! present partial cross section
    EN_REC_PREVIOUS_PARTIAL_XS  =   3, &  ! previous partial cross section (if available)
    EN_REC_2NDLAST_PARTIAL_XS   =   4, &  ! the second last partial cross sections (if available)
    EN_REC_COUNT                =   5     ! number of records per energy
  end enum

  enum, bind( C )
  enumerator ::  &
    REC_LAST_USED         =   0, &  ! not actually a valid record index; it is used as a convention by db[wr][ri] subroutines to reuse the last used record index 
    REC_PARAMS            =   1, &  ! this record stores an array of real parameters, followed by an array of integer parameters
    REC_EN_START          =   2     ! start of energy records
  end enum


  integer, parameter :: LABEL_SIZE = 40  ! in characters
  integer, parameter :: NUM_CHARS_IN_INT = 4
  integer, parameter :: NUM_INTS_IN_LABEL = LABEL_SIZE / NUM_CHARS_IN_INT


  enum, bind( C )
  enumerator ::  &
    SAV_INTPARAMS_LABEL          =   1, &  ! the first 10 integers are used to store the 40 characters of the label
    SAV_INTPARAMS_JTOT           =  SAV_INTPARAMS_LABEL + NUM_INTS_IN_LABEL, &
    SAV_INTPARAMS_JTOTD          =  12, &
    SAV_INTPARAMS_JLPAR          =  13, &
    SAV_INTPARAMS_NU             =  14, &
    SAV_INTPARAMS_NUTOP          =  15, &
    SAV_INTPARAMS_NUD            =  16, &
    SAV_INTPARAMS_NERG           =  17, &
    SAV_INTPARAMS_NCHMAX         =  18, &
    SAV_INTPARAMS_CSFLAG         =  19, &
    SAV_INTPARAMS_WRITS          =  20, &
    SAV_INTPARAMS_NUCROS         =  21, &
    SAV_INTPARAMS_JTOP           =  22, &
    SAV_INTPARAMS_NLEV_START     =  23
  end enum

  enum, bind( C )
  enumerator ::  &
    SAV_REALPARAMS_RTMN1         =   1, &
    SAV_REALPARAMS_RTMX1         =   2, &
    SAV_REALPARAMS_DINSID        =   3, &
    SAV_REALPARAMS_ENERG_START   =  11

  end enum

contains
  function get_sav_int_params_count()
  use mod_coener, only: max_en
  integer :: get_sav_int_params_count
  get_sav_int_params_count = SAV_INTPARAMS_NLEV_START + max_en -1
  end function

  function get_sav_real_params_count()
  use mod_coener, only: max_en
  integer :: get_sav_real_params_count
  get_sav_real_params_count = SAV_REALPARAMS_ENERG_START + max_en -1
  end function

end module mod_savfile

module mod_hibrid5



contains
! ---------------------------------------------------------------------------
! hibrid5 library
!
! subroutines included:
!
!   soutpt   driver to sum up partial cross sections and to write out s-matric
!   prpartr   writes out partial cross sections
!   nusum    sums up partial cross sections over nu in case nucros=true
!   intpol   interpolates and sums up partial cross sections
!   dbout    buffers out partial cross sections and their labels
!   xwrite   writes out integral cross sections
!   readpc   to write out selected partial cross sections
!   intchk   checks consistency of restart file
!   restrt (rsave)   restart utility (useful when system crash)
!   intcrs   driver to compute integral cross sections from s-matrices
!   intcr    computes integral cross sections from s-matrices
!   tsqmat   computes squared t-matrix from s-matrix
!   partcr   computes partial cross sections from squared t-matrix
!
! ---------------------------------------------------------------------------
subroutine soutpt (tsq, sr, si, scmat, &
                   bqs, isc1, isc2, sc1, sc2, &
                   jlev, elev, inlev, jtot, jfirst, &
                   jfinal, jtotd, nu, numin, numax, nud,jlpar,ien, &
                   ipos, csflag, flaghf, prsmat, prt2, t2test, &
                   writs, wrpart, prpart, wrxsec, prxsec, twomol, &
                   nucros, firstj,nlevel, nlevop, nopen, nmax, &
                   twojlp, jlpold)
! ---------------------------------------------------------------------------
!  subroutine to:
!                1. write out the elements of the s-matrix and modulus squared
!                   t-matrix (file 9)
!                2. write out selected elements of the s-matrix to disk
!                   (units 45 - (44+ien) )
!                3. compute partial cross sections.  if desired write them
!                   to units 25 - (24+ien)
!                4. in the case of coupled-states calculations accumulate
!                   partial cross sections for all values of the projection
!                   index [units 35 - (34+ien) ]
!                5. add the partial cross sections onto the integral cross
!                   sections
!  author:  millard alexander
!  heavily modified by h.-j. werner
!  current revision date: 24-jan-2012 by p.dagdigian
! ---------------------------------------------------------------------------
!  variables in call list
!    tsq:     on input:  tsq contains the nopen x nopen square of the t-matrix
!             on return:  nlevop x nlevop block of integral cross sections
!                         (in cs calculation only for nu=numax)
!    scmat:   on return:  nlevop x nlevop block partial cross sections
!                         (in cs calculation summed up to present nu)
!    sr:      on input:   the upper-left nopen x nopen block of sr
!                         contains the real part of the s-matrix
!    si:      on input:   the upper-left nopen x nopen block of si
!                         contains the imaginary part of the s-matrix
!    bqs:     rotational angular momenta, orbital angular momenta, and
!             additional quantum index for each channel
!             if the calculation involves the collisions of two diatomic
!             molecules, the bqs%jq = j1 + 10000 j2, where j1 and j2 are the
!             rotational quantum numbers of each molecule
!    isc1,isc2: scratch vectors (min length nopen)
!    sc1, sc2:  scratch  matrices (min length nopen x nopen)
!    jlev:   rotational angular momenta of each energetically open level
!    elev:   energy (in hartree) of each energetically open level
!    inlev:  additional quantum index for each energetically open level
!             if the calculation involves the collisions of two diatomic
!             molecules, the jlev = j1 + 10000 j2, where j1 and j2 are the
!             rotational quantum numbers of each molecule
!
!    jtot:    current value of total angular momentum
!    jfirst:  initial value of total angular momentum
!    jfinal:  final value of total angular momentum
!    jtotd:   step size for total angular momentum
!    nu:      current value of the coupled-states projection index
!    numin:   initial value of the coupled-states projection index
!    numax:   final value of the coupled-states projection index
!    jlpar:   parity of channels
!               eps * (-1)**(j+l-jtot) = jlpar
!    note!!!   if flaghf = .true.( see below), then the true values
!    of the rotational quantum numbers, the total angular momentum,
!    and the coupled-states projection index are equal to the values
!    stored in bqs%jq, jtot, jfirst, jfinal, nu, numin, and numax plus 1/2
!    ien:     ordinal number of total energy at which this routine is
!             called.  i.e. if ien = 2, we are now at second value of
!             the total energy
!    ipos     if .true., 132 column printer
!             if .false., 80 column printer
!    csflag   if .true. coupled-states calculation
!             if .false., close-coupled calculation
!    flaghf:  if .true., then system with half-integer spin
!             if .false., then system with integer spin
!    prsmat    if .true., real and imaginary parts of s-matrix are printed
!                        to file 9
!    prt2   if .true., modulus squared of t-matrix is printed to file 9
!    t2test   if .true., first two columns of modulus squared of t-matrix are
!                        printed to file 9
!    writs    if .true., real and imaginary parts of selected s-matrix element
!                        are saved in files 45,46,47 ...
!    prpart    if .true., partial cross sections are printed to file 9
!    wrpart   if .true., partial cross sections are save in files 25,26,27,...
!    prxsec   if .true., integral cross sections are printed to file 9
!    wrxsec   if .true., integral cross sections are save in files 15,16,17 ..
!    twomol   if .true., then molecule-molecule cross section
!    nlevel   number of energetically distinct levels included in channel basi
!    nlevop   number of energetically distinct levels included in channel basi
!             which are open asymptotically
!    nopen    on entry:  number of open channels
!    nmax     on entry:  maximum row dimension of matrices
!    firstj   if .true. on entry, header is written to s-matrix file

!    nlevel   number of energetically distinct levels included in channel basi
!    variables in module constants
!    econv:    conversion factor from cm-1 to hartrees
!    xmconv:   converson factor from amu to atomic units
!    ang2c:     conversion factor from square bohr to square angstroms
!  ---------------------------------------------------------------------------
use mod_cosout
use constants
use mod_coqvec, only: nphoto
use mod_cocent, only: cent
use mod_coeint, only: eint
use mod_coener, only: ener => energ
use mod_hibrid2, only: mxoutd, mxoutr
use funit
use mod_savfile, only: REC_EN_START, EN_REC_COUNT
use mod_parpot, only: label=>pot_label
use mod_ered, only: ered, rmu
use mod_phot, only: photof, wavefn, boundf
use mod_surf, only: flagsu
use mod_hiutil, only: dater
use mod_himatrix, only: transp
use mod_hismat, only: wrhead, swrite
use mod_hitypes, only: bqs_type
implicit double precision (a-h,o-z)
real(8), intent(inout) :: tsq(nmax,nmax)
real(8), intent(inout) :: sr(nmax,nmax)
real(8), intent(inout) :: si(nmax,nmax)
real(8), intent(out) :: scmat(nmax,nmax)
type(bqs_type), intent(in) :: bqs
integer, intent(out) :: isc1(nopen)
integer, intent(out) :: isc2(nopen)
real(8), intent(out) :: sc1(nmax, nmax)
real(8), intent(out) :: sc2(nmax, nmax)
integer, intent(in) :: jlev(nlevop)
real(8), intent(in) :: elev(nlevop)
integer, intent(in) :: inlev(nlevop)
integer, intent(in) :: jlpold ! old value of parity, used to insure correct accumulation
                              ! of all partial waves in cases where jlpar=0 .

logical ipos, csflag, prsmat, prt2, writs, wrpart, prpart, &
        wrxsec, prxsec, flaghf, t2test, firstj, twomol, &
        nucros, faux, twojlp
type(bqs_type) :: packed_bqs

character*20 cdate
integer :: soutpt_sc_file = 1
!
data izero, ione /0, 1/
xjtot = jtot
if (flaghf .and. .not. csflag) xjtot = jtot + 0.5d0
!     print the s-matrix
if (prsmat .and. nopen .gt. 0) then
  if (photof) then
    call transp(sr, nopen, nmax)
    call transp(si,nopen,nmax)
  endif
  if (photof) then
     write (9, 100) ener(ien)
100      format(1h /, &
       ' ** REAL PART OF ASYMPTOTIC WAVEFUNCTION; ENERGY =', &
       1pe12.4,/)
     call mxoutr (9, sr, nphoto, nopen, nmax, 0, ipos)
  else
    if (.not. wavefn) then
      write (9,101) ener(ien)
101       format(1h /' ** REAL PART OF THE S MATRIX; ENERGY =', &
       1pe12.4)
    else
      write (9, 100) ener(ien)
    endif
!          call mxoutd (9, sr, nopen, nmax, 1, ipos)
    call mxoutd (9, sr, nopen, nmax, 0, ipos)
  endif
  if (photof) then
     write (9, 110) ener(ien)
110      format(1h /, &
       ' ** IMAGINARY PART OF ASYMPTOTIC WAVEFUNCTION; ENERGY =', &
       1pe12.4,/)
     call mxoutr (9, si, nphoto,nopen, nmax, 0, ipos)
  else
    if (.not. wavefn) then
      write (9,111) ener(ien)
111       format(1h /' ** IMAGINARY PART OF THE S MATRIX; ENERGY =', &
       1pe12.4)
    else
      write (9, 110) ener(ien)
    endif
!          call mxoutd (9, si, nopen, nmax, 1, ipos)
    call mxoutd (9, si, nopen, nmax, 0, ipos)
  endif
end if
if (prt2) then
  if (.not. photof) then
    write (9,120) ener(ien)
120     format(1h /' ** MODULUS SQUARED T-MATRIX; ENERGY =', &
       1pe12.4)
    call mxoutd (9, tsq, nopen, nmax, 1, ipos)
  else
    write (9, 121) ener(ien)
121     format(1h /, &
      ' ** NORMALIZED TRANSITION PROBABILITIES; ENERGY =', &
       1pe12.4)
    call mxoutr(9, tsq, nphoto, nopen, nopen, izero, ipos)
  endif
end if
if (t2test) then
  write (9, 130) ener(ien)
130   format (1h /, &
   ' ** FIRST 2 COLUMNS OF MODULUS SQUARED T-MATRIX; ENERGY =', &
       1pe12.4)
  nlow = 1
  nhigh = min (nopen, 10)
140   write (9, 145)
145   format (1h )
  do 155  i = nlow, nhigh
    write (9, 150)  i, tsq(i,1), tsq(i,2)
150     format (i4, 1x, 2 (1pe12.4) )
155   continue
  if (nhigh .lt. nopen) then
    nlow = nhigh + 1
    nhigh = min (nlow + 9, nopen)
    go to 140
  end if
end if
! if photodissociation calculation or wavefunction desired, return here
if (photof .or. wavefn) return
if (writs .and. nopen .gt. 0) then
!  here if real and imaginary parts of s-matrix for selected transitions
!  are to be written out to unit (FUNIT_SMT_START+ien-1)
    nfile = FUNIT_SMT_START + ien - 1
    if (firstj) then
       call dater (cdate)
       call wrhead(nfile, cdate, &
                ered, rmu, csflag, &
                flaghf, flagsu, twomol, nucros, jfirst, jfinal, &
                jtotd, numin, numax, nud, nlevel, nlevop, nnout, &
                jlev, inlev, elev, jout)
    end if
    call swrite (sr, si, jtot, jlpar, nu, bqs, isc1, &
                 packed_bqs, sc2, nfile, nmax, nopen)
end if

if (.not. prxsec .and. .not. wrxsec .and. .not.prpart &
    .and. .not.wrpart) return
call  partcr (tsq,  scmat, nopen, &
                   bqs, bqs, &
                   inlev, jlev, elev, jtot, nu, &
                   csflag, flaghf,twomol,flagsu, &
                   nlevop,nmax)
!  here if coupled-states calculation
if (csflag.and..not.nucros) then
!  here if partial cross sections are desired
!  if first value of projection index, then initialize matrix of
!  partial cross sections which will be stored on unit (FUNIT_APCS_START+ien-1)
!  tsq is used as scratch matrix here
  nfile = FUNIT_APCS_START + ien - 1 
  if(nu.gt.numin) then
    rewind nfile
    read(nfile) ((sc1(j,i),j=1,nlevop),i=1,nlevop)
    do 240 i=1,nlevop
    do 240 j=1,nlevop
240     scmat(j,i)=scmat(j,i)+sc1(j,i)
  end if
  if(nu.lt.numax) then
    rewind nfile
    write(nfile) ((scmat(j,i),j=1,nlevop),i=1,nlevop)
  end if
end if
if((prpart .or. wrpart) .and. .not. nucros .and. &
          (.not.csflag.or.nu.eq.numax)) then
! write out partial cross sections
   call prpartr (scmat,jlev, elev, inlev, jtot, jfirst, &
               jfinal, jtotd, nu, numin, numax, nud, jlpar, ien, &
               ipos, csflag, flaghf, wrpart, prpart, twomol, &
               twojlp, &
               nucros,nlevel, nlevop, nopen, nmax)
end if
!  sum up partial cross sections to get integral cross sections
if (.not. csflag .or. (csflag .and. (nu .eq. numax) ) ) then
   irec=(ien-1)* EN_REC_COUNT + REC_EN_START
   faux = .false.
   call intpol(irec,jtot,jfirst,jfinal,jtotd,jlpar,jlpold, &
                jlev,nmax,nlevop,tsq,scmat,faux,soutpt_sc_file)
end if
return
end
! ---------------------------------------------------------------------------
subroutine prpartr(scmat,jlev, elev, inlev, jtot, jfirst, &
                  jfinal, jtotd, nu, numin, numax, nud, jlpar,ien, &
                  ipos, csflag, flaghf, wrpart, prpart, twomol, &
                  twojlp,nucros,nlevel, nlevop, nopen, nmax)
! ---------------------------------------------------------------------------
!
!  subroutine to print partial cross sections
!
!  latest revision date: 27-oct-1995 by mha
! ---------------------------------------------------------------------------
use mod_cosout, only: nnout, jout
use mod_coiout, only: niout, indout
use mod_coisc2, only: nj, jlist => isc2 ! nj,jlist(10)
use constants
use mod_parpot, only: potnam=>pot_name, label=>pot_label
use mod_ered, only: ered, rmu
use mod_surf, only: flagsu
use mod_hiutil, only: dater
implicit double precision (a-h,o-z)
real(8), intent(in) :: scmat(nmax, nlevop)
integer, intent(in) :: jlev(nlevop)
real(8), intent(in) :: elev(nlevop)
integer, intent(in) :: inlev(nlevop)
character*20 cdate
character*40 form
logical ipos, csflag, wrpart, prpart, flaghf, twomol, nucros, &
        twojlp, headf
!  write partial opacity to unit (24+ien) if desired
!  in cs calculation this is only done if nu = numax, in which
!  case the partial opacity has been summed over all projection indices
if(csflag.and.nu.lt.numax) return
! make a list of pointers
nj = 0
do 2 i=1, iabs(nnout)
   jo = jout(i)
   in = 0
   do 2 ii=1,max(1,niout)
   if(niout.gt.0) in = indout(ii)
   do 1 j=1, nlevop
     j1 = jlev(j)
!           if(j1.ne.jo.or.(in.ne.0.and.in.ne.inlev(j))) goto 1
     if(j1.ne.jo.or.in.ne.inlev(j)) goto 1
     nj = nj + 1
     jlist(nj) = j
1  continue
2 continue
if (wrpart) then
  nfile = 24 + ien
  headf=.false.
  if (nucros .and. nu.eq.numin) headf=.true.
  if (twojlp) then
    if (jtot .eq. jfirst .and. jlpar .eq. 1) headf=.true.
  else
    if (jtot .eq. jfirst) headf=.true.
  endif
  if (headf) then
    call dater (cdate)
    write (nfile, 10) cdate
    write (nfile, 10) label
    write (nfile, 10) potnam
10     format (1x, a)
    write (nfile, 20) ered * econv, rmu * xmconv
20     format (f10.3,f15.11)
    write (nfile, 30) csflag, flaghf, flagsu, twomol, nucros
30     format (5l3)
    write (nfile, 40) jfirst, jfinal, jtotd, &
                       numin, numax, jlpar, nud
    if(nnout.lt.0) then
      write (nfile, 40) nj, nlevop
      write (nfile, 40) (jlev(i), inlev(i), i = 1, nlevop)
40       format (24i5)
      write (nfile, 50) (elev(i), i = 1, nlevop)
50       format (8g16.9)
    else
      write (nfile, 40) nj, nj
      write (nfile, 40) (jlev(jlist(i)), inlev(jlist(i)), &
                        i = 1, nj)
      write (nfile, 50) (elev(jlist(i)), i = 1, nj)
    end if
  end if
  if(nucros) then
    write (nfile, 40) nu
  else
    write (nfile, 40) jtot
  end if
! loop over initial states as specified in jout
  form='(1x,i2,i5,(t9,10d12.5))'
  if(nnout.lt.0) then
    do 60 i=1,nj
    ii=jlist(i)
    write (nfile, form) jlev(ii),inlev(ii),(scmat(ii,j), &
                      j = 1, nlevop)
60     continue
  else
    do 61 i=1,nj
    ii=jlist(i)
    write (nfile, form) jlev(ii),inlev(ii),(scmat(ii,jlist(j)), &
                      j = 1, nj)
61     continue

  end if
end if
! if cs calculation, write out partial cross sections summed over
! projection index
if (prpart) then
  if (jtot .eq. jfirst.or.(nucros.and.nu.eq.numin))  then
!  here if output of partial cross sections requested
    if( .not. twomol) then
       write (9, 100) ered * econv
100        format (/' LEVEL LIST FOR PARTIAL CROSS SECTIONS AT E=', &
        f8.2,' CM-1'/'   N    J  INDEX  EINT(CM-1)',/)
        ie=nlevop
        if(nnout.gt.0) ie=nj
        do 120  ii = 1, ie
        i=ii
        if(nnout.gt.0) i=jlist(ii)
        if (.not. flaghf) then
          write (9, 110) i, jlev(i), inlev(i), elev(i) * econv
110           format (i4, i5, i6, f11.3)
        else
          write (9, 115) i, (jlev(i)+0.5d0), inlev(i), &
                         elev(i) * econv
115           format (i4, f5.1, i6, f11.3)
        end if
120       continue
    else
      write (9, 125) ered * econv
125       format (/' LEVEL LIST FOR PARTIAL CROSS SECTIONS AT E=', &
        f8.2,' CM-1'/'   N   J1   J2  INDEX  EINT(CM-1)'/)
        ie=nlevop
        if(nnout.gt.0) ie=nj
        do 135  ii = 1, ie
        i=ii
        if(nnout.gt.0) i=jlist(ii)
        jj1 = jlev(i) / 10
        jj2 = mod(jlev(i), 10)
        write (9, 130) i, jj1, jj2, inlev(i), elev(i) * econv
130         format (i4, 2i5, i6, f11.3)
135       continue
    end if
  end if
  if (csflag) then
    xnu = nu
    if (flaghf) xnu = xnu + 0.5d0
    if (flagsu) then
      write (9, 205) xnu, ered * econv
205       format(/' DEGENERACY AVERAGED TRANSITION', &
         ' PROBABILITIES MULTIPLIED BY 2;  M=', f5.1,' E=',f8.2)
    else
      if(.not.nucros) write (9, 210) jtot, ered * econv
      if(nucros) write (9, 215) nu,ered * econv
210       format(/' CS PARTIAL CROSS SECTIONS FOR', &
         ' LBAR=',i3,', SUMMED OVER NU. E=',f8.2)
215       format(/' CS PARTIAL CROSS SECTIONS FOR', &
         ' NU=',i3,', SUMMED OVER LBAR. E=',f8.2)
    end if
  else
    write (9,220) jtot,jlpar,ered * econv
220     format (/' CC PARTIAL CROSS SECTIONS FOR', &
      ' JTOT=',i3,', JLPAR=',i2,', E=',f8.2)
  end if
  if(twomol) then
    if(nnout.lt.0) then
      if(ipos) write(9,230) (j,j=1,nlevop)
      if(.not.ipos) write(9,231) (j,j=1,nlevop)
    else
      if(ipos) write(9,230) (jlist(j),j=1,nj)
      if(.not.ipos) write(9,231) (jlist(j),j=1,nj)
    end if
230     format(/1x,'  INITIAL STATE',t27,'FINAL STATE'/ &
            1x,'  N   J1   J2  INDEX',(t21,i7,9i11))
231     format(/1x,'  INITIAL STATE',t27,'FINAL STATE'/ &
            1x,'  N   J1   J2  INDEX',(t21,i7,4i11))
! don't remove blank in this format!!
    form='(1x,i3,2i5,i6, (t21,10d11.4))'
  else
    if(nnout.lt.0) then
      if(ipos) write(9,232) (j,j=1,nlevop)
      if(.not.ipos) write(9,233) (j,j=1,nlevop)
    else
      if(ipos) write(9,232) (jlist(j),j=1,nj)
      if(.not.ipos) write(9,233) (jlist(j),j=1,nj)
    end if
232     format(/1x,'  INITIAL STATE',t24,'FINAL STATE'/ &
            1x,'  N   J   INDEX ',(t18,i7,9i11))
233     format(/1x,'  INITIAL STATE',t24,'FINAL STATE'/ &
            1x,'  N   J   INDEX ',(t18,i7,4i11))
    if(flaghf) then
      spin=0.5d0
      form='(1x,i3,f5.1,i6,(t18,10d11.4))'
    else
! don't remove blank in this format!!
      form='(1x,i3,i4,i7,  (t18,10d11.4))'
    end if
  end if
  if(.not.ipos) form(21:22)=' 5'
  if(nnout.lt.0) then
    if(twomol) then
      do 221 i=1,nj
      ii=jlist(i)
      jj1 = jlev(ii) / 10
      jj2 = mod(jlev(ii), 10)
      write (9, form) ii,jj1,jj2,inlev(ii),(scmat(ii,j), &
                      j = 1, nlevop)
221       continue
    else if(flaghf) then
     do 222 i=1,nj
      ii=jlist(i)
      write (9, form) ii,jlev(ii)+0.5d0,inlev(ii),(scmat(ii,j), &
                      j = 1, nlevop)
222       continue
    else
      do 223 i=1,nj
      ii=jlist(i)
      write (9, form) ii,jlev(ii),inlev(ii),(scmat(ii,j), &
                      j = 1, nlevop)
223       continue
    end if
  else
    if(twomol) then
      do 224 i=1,nj
      ii=jlist(i)
      jj1 = jlev(ii) / 10
      jj2 = mod(jlev(ii), 10)
      write (9, form) ii,jj1,jj2,inlev(ii),(scmat(ii,jlist(j)), &
                      j = 1, nj)
224       continue
    else if(flaghf) then
     do 225 i=1,nj
      ii=jlist(i)
      write (9, form) ii,jlev(ii)+0.5d0,inlev(ii), &
          (scmat(ii,jlist(j)),j = 1, nj)
225       continue
    else
      do 226 i=1,nj
      ii=jlist(i)
      write (9, form) ii,jlev(ii),inlev(ii),(scmat(ii,jlist(j)), &
                      j = 1, nj)
226       continue
    end if
  end if
end if
return
end
! ----------------------------------------------------------------------
subroutine nusum (tsq, tq1, tq2, tq3, &
                 jlev,elev, inlev, jtot, jfirst, &
                 jfinal, jtotd, nu, numin, numax, nud, jlpar, &
                 nerg, ipos, csflag, flaghf, wrpart, prpart, &
                 twomol, nucros, nlevel, nlev, nopen, nmax, tmp_file)
! ---------------------------------------------------------------------------
!
!  subroutine to sum partial cross sections over nu
!
!  latest revision date: 21-mar-1992 by mha
! ---------------------------------------------------------------------------
use constants
use mod_coener, only: energ
use mod_ered, only: ered, rmu
use mod_savfile, only: REC_EN_START, EN_REC_COUNT, EN_REC_PRESENT_INTEGRAL_XS
implicit double precision (a-h,o-z)
logical ipos, csflag, wrpart, prpart, flaghf, twomol, nucros,vrai
dimension tsq(nmax,1),tq1(nmax,1),tq2(nmax,1),tq3(nmax,1), &
          jlev(1),elev(1),inlev(1),nlev(1)
integer :: tmp_file
!
do 100 ien = 1,nerg
  ener = energ(ien)
  ered = ener/econv
  irec=(ien-1)* EN_REC_COUNT + REC_EN_START
  nlevop=nlev(ien)
  call dbin(tmp_file,irec+EN_REC_PRESENT_INTEGRAL_XS,jl,jlp,nn,tsq,nmax,nlevop)
! print partial cross sections, summed over lbar, for given nu
  if(prpart.or.wrpart) then
    vrai=.false.
    call prpartr (tsq,jlev, elev, inlev, jtot, jfirst, &
               jfinal, jtotd, nu, numin, nu, nud, jlpar, ien, &
               ipos, csflag, flaghf, wrpart, prpart, twomol, &
               vrai,nucros,nlevel, nlevop, nopen, nmax)
    endif
! sum up over nu
  irec=(nerg+ien-1)* EN_REC_COUNT + REC_EN_START
  vrai=.true.
  call intpol(irec,nu,numin,numax,nud,jlpar,jlpar, &
              jlev, nmax,nlevop,tq1,tsq,vrai,tmp_file)
100 continue
return
end
! ----------------------------------------------------------------------
subroutine intpol(irec,jl3,j1,j2,jd,jp,jpi,jlev, &
                  nmax,n,q,q3,nucros,tmp_file)
! ----------------------------------------------------------------------
!
! subroutine to sum up and interpolate partial cross sections
!
! ----------------------------------------------------------------------
!
!     this routine stores:
!       the present integral cross sections on record irec
!       the previous integral cross sections on record irec+1
!       the present  partial cross sections on record irec+2
!       the previous partial cross sections on record irec+3 (if available)
!       the second last partial cross sections on record irec+4 (if available)
! ----------------------------------------------------------------------
use mod_savfile, only: EN_REC_PRESENT_INTEGRAL_XS, EN_REC_PREVIOUS_INTEGRAL_XS, EN_REC_PRESENT_PARTIAL_XS, EN_REC_PREVIOUS_PARTIAL_XS, EN_REC_2NDLAST_PARTIAL_XS
implicit none
integer, intent(in) :: irec  ! index of the first record (in the sense of iadr(irec, ifil)) used to store the current energy data in the sav file (EN_REC_COUNT records are used for each energy level) 
integer, intent(in) :: jl3  ! jtot or nu value of present partial cross sections in q3
integer, intent(in) :: j1    ! start value of jtot (nu)
integer, intent(in) :: j2    ! end value of jtot (nu)
integer, intent(in) :: jd    ! increment between jtot (nu) values
integer, intent(in) :: jp    ! present parity
integer, intent(in) :: jpi   ! the first parity
integer, intent(in) :: jlev(n)
integer, intent(in) :: nmax  ! the first dimension of the matrices
integer, intent(in) :: n     ! the actual dimension of the matrices (= nlevop)
real(8), intent(out) :: q(nmax,n)  ! on return contains the present integral cross sections, summed up to jl3
real(8), intent(in) :: q3(nmax,n)
logical, intent(in) :: nucros
integer, intent(in) :: tmp_file
real(8) :: q1(nmax,n) ! scratch array, used for previous partial waves
real(8) :: q2(nmax,n) ! scratch array, used for previous partial waves
integer :: i, j
integer :: jl, num_partial_xs, jl1, jp1, jl2, jp2, nn
integer :: jpl               ! the previous parity (l=last?)
!
!
if(jd <= 1) then
  !
  ! here for jd == 1
  !
  if(jl3 == j1 .and. jp == jpi) then
    ! initialize at first j (or nu)
    do i=1,n
      do j=1,n
       q2(j,i)=0
       q(j,i)=q3(j,i)
      end do
    end do
    call dbout(irec + EN_REC_PRESENT_INTEGRAL_XS,jl3,jp,1,q,nmax,n)
    call dbout(irec + EN_REC_PREVIOUS_INTEGRAL_XS,-1,jp,1,q2,nmax,n)
    call dbout(irec + EN_REC_PRESENT_PARTIAL_XS,jl3,jp,0,q3,nmax,n)
  else
    ! jd == 1 and (jl3,jp) /= (j1,jpi) => simply add
    call dbin(tmp_file, irec + EN_REC_PRESENT_INTEGRAL_XS,jl,jpl,num_partial_xs,q,nmax,n)
    if(jp == jpl .and. jl+1 /= jl3) then
      write(6,15) jp,jl,jl3
15       format(/' ERROR IN INTPOL: JP, JL, JL3:',3i4)
      call exit
    end if
    do i=1,n
      do j=1,n
        q2(j,i)=q(j,i)
        q(j,i)=q(j,i)+q3(j,i)
      end do
    end do
    call dbout(irec + EN_REC_PRESENT_INTEGRAL_XS,jl3,jp,1,q,nmax,n)
    call dbout(irec + EN_REC_PREVIOUS_INTEGRAL_XS,jl,jpl,1,q2,nmax,n)
    call dbout(irec + EN_REC_PRESENT_PARTIAL_XS,jl3,jp,0,q3,nmax,n)
  end if
else
  !
  !.....here for jd /= 1
  !
  if(jl3 == j1 .and. jp == jpi) then
  !.....initialize for jd.ne.1
    do i=1,n
      do j=1,n
        q2(j,i)=0
        q(j,i)=q3(j,i)
      end do
    end do
    jl=-1
    call dbout(irec + EN_REC_PRESENT_INTEGRAL_XS,jl3,jp,1,q,nmax,n)
    call dbout(irec + EN_REC_PREVIOUS_INTEGRAL_XS,jl,jp,1,q2,nmax,n)
    call dbout(irec + EN_REC_PRESENT_PARTIAL_XS,jl3,jp,0,q3,nmax,n)
  else
    !.....jl is value up to which jtot has been summed so far
    call dbin(tmp_file, irec + EN_REC_PRESENT_INTEGRAL_XS,jl,jpl,num_partial_xs,q,nmax,n)
    if(jpl /= jp) then
      do i=1,n
        do j=1,n
          q2(j,i)=q(j,i)
          q(j,i)=q(j,i)+q3(j,i)
        end do
      end do
      call dbout(irec + EN_REC_PRESENT_INTEGRAL_XS,jl3,jp,1,q,nmax,n)
      call dbout(irec + EN_REC_PREVIOUS_INTEGRAL_XS,jl,jpl,1,q2,nmax,n)
      call dbout(irec + EN_REC_PRESENT_PARTIAL_XS,jl3,jp,0,q3,nmax,n)
    else
      nn=min(num_partial_xs+1,3)
      call dbout(irec + EN_REC_PREVIOUS_INTEGRAL_XS,jl,jpl,nn,q,nmax,n)
      call dbin(tmp_file, irec + EN_REC_PRESENT_PARTIAL_XS,jl2,jp2,nn,q2,nmax,n)
      if(num_partial_xs == 1) then
        call intpl2(jl, jl2 ,jl3 , nmax, n, q, q2, q3, jlev, nucros)
        !
        call dbout(irec + EN_REC_PRESENT_INTEGRAL_XS,jl3,jp,2,q,nmax,n)
        call dbout(irec + EN_REC_PRESENT_PARTIAL_XS,jl3,jp,0,q3,nmax,n)
        call dbout(irec + EN_REC_PREVIOUS_PARTIAL_XS,jl2,jp2,0,q2,nmax,n)
      else
        ASSERT (num_partial_xs > 1)
        ! the previous partial cross section is available; use it for interpolation
        call dbin(tmp_file, irec + EN_REC_PREVIOUS_PARTIAL_XS,jl1,jp1,nn,q1,nmax,n)
        call intpl3(jl, jl1, jl2, jl3, nmax, n, q, q1, q2, q3, jlev, nucros)
        !
        call dbout(irec + EN_REC_PRESENT_INTEGRAL_XS,jl3,jp,3,q,nmax,n)
        call dbout(irec + EN_REC_PRESENT_PARTIAL_XS,jl3,jp,0,q3,nmax,n)
        call dbout(irec + EN_REC_PREVIOUS_PARTIAL_XS,jl2,jp2,0,q2,nmax,n)
        call dbout(irec + EN_REC_2NDLAST_PARTIAL_XS,jl1,jp1,0,q1,nmax,n)
      end if
    end if
  end if
end if  ! end of handling for jd /= 1
! q = 0.0d0
q1 = 0.01d0
q2 = 0.01d0

end subroutine

subroutine intpl3(jl,jl1,jl2,jl3,nmax,n,q,q1,q2,q3,jlev,nucros)
use mod_coisc8, only: list => isc8 ! list(20)
implicit none
integer, intent(in) :: jl
integer, intent(in) :: jl1
integer, intent(in) :: jl2
integer, intent(in) :: jl3
integer, intent(in) :: nmax
integer, intent(in) :: n
real(8), intent(out) :: q(nmax, n)
real(8), intent(in) :: q1(nmax, n)
real(8), intent(in) :: q2(nmax, n)
real(8), intent(in) :: q3(nmax, n)
integer, intent(in) :: jlev(n)
logical, intent(in) :: nucros
real(8) :: f0, f1, f2
integer :: i, j, ii, jj
integer :: nleq, nlge
integer :: nu
integer :: jf
!
!.....we will now sum from ji to jf
!
if(jl /= jl2) then
  write(6,32) jl,jl2
32   format(/' ERROR IN INTPL3: JL=',i0,' NOT EQUAL JL2=',i3)
  call exit
end if
jf = jl3
if(.not. nucros) then
  call prefac(jl1,jl2,jl3,jl+1,jf,f0,f1,f2)
  do i=1, n
    do j=1, n
      q(j,i) = q(j,i) + f0 * q2(j,i) &
                      + f1 * (q2(j,i)-q3(j,i)) &
                      + f2 * (q1(j,i)-q2(j,i))
    end do
  end do
else
  ! here for sparse algorithm
  do nu = jl+1, jf
    call nulist(nu,jlev,list,n,nleq,nlge)
    !       write(6,*) 'three-point, jl=',jl,'  jf=',jf,'  nu=',nu,nleq,nlge
    !       do 42 nuu=jl+1,nu
    !       call prefac(jl1,jl2,nu,nuu,nuu,f0,f1,f2)
    call prefac(jl1,jl2,nu,jl+1,nu,f0,f1,f2)
    do ii = 1, nlge
      i = list(ii)
      do jj = 1, nlge
        if(nu < jf .and. jj > nleq .and. ii > nleq) cycle
        j = list(jj)
        q(j,i) = q(j,i) + f0 * q2(j,i) &
                        + f1 * (q2(j,i)-q3(j,i)) &
                        + f2 * (q1(j,i)-q2(j,i))
      end do
    end do
  end do
end if
!
return
end subroutine
!
!.....linear interpolation if only two points available
!
subroutine intpl2(jl,jl2,jl3,nmax,n,q,q2,q3,jlev,nucros)
use mod_coisc8, only: list => isc8 ! list(20)
implicit none
integer, intent(in) :: jl
integer, intent(in) :: jl2
integer, intent(in) :: jl3
integer, intent(in) :: nmax
integer, intent(in) :: n
real(8), intent(out) :: q(nmax, n)
real(8), intent(in) :: q2(nmax, n)
real(8), intent(in) :: q3(nmax, n)
integer, intent(in) :: jlev(n)
logical, intent(in) :: nucros
real(8) :: f0, f1, f2
integer :: i, j, ii, jj
integer :: jf
integer :: nleq, nlge
integer :: nu

if(jl /= jl2) then
  write(6,32) jl, jl2
32   format(/' ERROR IN INTPL2: JL=',i3,' NOT EQUAL JL2=',i3)
  call exit
end if
jf=jl3
if(.not.nucros) then
  call prefac(jl2, jl2, jl3, jl + 1, jf, f0, f1, f2)
  do i = 1, n
    do j = 1, n
      q(j,i) = q(j,i) + f0 * q2(j, i) &
                      + f1 * (q2(j, i)-q3(j, i))
    end do
  end do
else
!......here for sparse algorithm
  do nu = jl + 1, jf
    !.....select contributing levels
    call nulist(nu, jlev, list, n, nleq, nlge)
    !       write(6,*) 'two-point jl=',jl,'  jf=',jf,'  nu=',nu,nleq,nlge
    !       do 62 nuu=jl+1,nu
    !       call prefac(jl2,jl2,nu,nuu,nuu,f0,f1,f2)
    call prefac(jl2, jl2, nu, jl + 1, nu, f0, f1, f2)
    do ii = 1, nlge
      i=list(ii)
      do jj = 1, nlge
        if(nu < jf .and. jj > nleq .and. ii > nleq) cycle
        j=list(jj)
        q(j,i) = q(j,i) + f0 * q2(j,i) &
                        + f1 * (q2(j,i)-q3(j,i))
      end do
    end do
  end do  ! end of loop on nu

!......now all levels for j.ge.nu=jl3
end if
!
return
end subroutine

subroutine nulist(nu,jlev,list,nj,nleq,nlge)
!.....makes a list of all channels with j.ge.nu
implicit none
integer, intent(in) :: nu
integer, intent(in) :: jlev(nj)
integer, intent(out) :: list(nj)
integer, intent(in) :: nj
integer, intent(out) :: nleq
integer, intent(out) :: nlge
integer :: i, nl
nl=0
do 10 i=1,nj
if(jlev(i).eq.nu-1) then
  nl=nl+1
  list(nl)=i
end if
10 continue
nleq=nl
do 20 i=1,nj
if(jlev(i).ge.nu) then
  nl=nl+1
  list(nl)=i
end if
20 continue
nlge=nl
return
end
! ----------------------------------------------------------------------
subroutine prefac(j1,j2,j3,ji,jf,f0,f1,f2)
implicit double precision (a-h,o-z)
!.....compute prefactors for interpolation from ji to jf
!.....values must be known for j1,j2,j3
!.....j2 is taken as expansion point
d1=0
d2=0
do 10 j=ji,jf
d1=d1+(j-j2)
10 d2=d2+(j-j2)**2
x1=j1-j2
x3=j3-j2
f0=jf-ji+1
!     write(6,*) 'prefac:',j1,j2,j3,x1,x3
if(x1.eq.0) then
  f1=-d1/x3
  f2=0.0d0
else if(x3.eq.0) then
  f1=0.0d0
  f2=d1/x1
else
  dd=1.d0/(x1*x3*(x3-x1))
  f1=dd*(d1*x1**2-d2*x1)
  f2=dd*(d1*x3**2-d2*x3)
end if
return
end
! ----------------------------------------------------------------------
subroutine dbout(irec,i1,i2,i3,q,nmax,n)
! ----------------------------------------------------------------------
!
!  subroutine to buffer out cross sections matrices together
!  with their labels
!
! ---------------------------------------------------------------------------
use mod_fileid, only: FILEID_SAV
use mod_savfile, only: REC_LAST_USED
implicit none
integer, intent(in) :: irec
integer, intent(in) :: i1, i2, i3
integer, intent(in) :: nmax, n
real(8), intent(in) :: q(nmax, n)
integer :: ifile, i
ifile = FILEID_SAV
call dbwi(i1,1,ifile,irec)
call dbwi(i2,1,ifile,REC_LAST_USED)
call dbwi(i3,1,ifile,REC_LAST_USED)
do i=1,n
  call dbwr(q(1,i),n,ifile,REC_LAST_USED)
end do
call dbwc(ifile,irec)
return
end

! ----------------------------------------------------------------------
!
!  subroutine to read cross sections matrices together
!  with their labels from the given file unit
!
! ---------------------------------------------------------------------------
subroutine dbin(ifile,irec,i1,i2,i3,q,nmax,n)
use mod_savfile, only: REC_LAST_USED
implicit none
integer, intent(in) :: ifile
integer, intent(in) :: irec
integer, intent(out) :: i1, i2, i3
integer, intent(in) :: nmax, n
real(8), dimension(nmax, n), intent(out) :: q
integer :: i
call dbri(i1,1,ifile,irec)
call dbri(i2,1,ifile,REC_LAST_USED)
call dbri(i3,1,ifile,REC_LAST_USED)
do i=1,n
  call dbrr(q(1,i),n,ifile,REC_LAST_USED)
end do
return
end
! ----------------------------------------------------------------------
subroutine xwrite (zmat, tq3, jlev, elev, inlev, nerg, energ, &
                   jfirst, jfinal, jtotd, csflag, flaghf, &
                   wrxsec, prxsec, ipos, twomol, nucros,nlevel, &
                   nlev, numin, numax, nud, jlpar, nmax, nmx, &
                   ihomo)
! ---------------------------------------------------------------------------
!  subroutine to write out integral cross sections
!  author:  millard alexander
!  latest revision date: 23-feb-2013 by p. dagdigian
!  ------------------------------------------------------------------
!  variables in call list:
!    zmat:    on return:  contains the nlevop x nlevop matrix of integral
!                         cross sections
!    jlev:   rotational angular momenta of each energetically open level
!    elev:   energy (in hartree) of each energetically open level
!    inlev:  additonal quantum index for each energetically open level
!    nerg:    number of total energies
!    energ:   array of total energies (cm-1)
!    jfirst:  initial value of total angular momentum
!    jfinal:  final value of total angular momentum
!    jtotd:   step size for total angular momentum
!    ipos:    if .true., 132 column printer
!             if .false., 80 column printer
!    csflag:  if .true. coupled-states calculation
!             if .false., close-coupled calculation
!    flaghf:  if .true., then system with half-integer spin
!             if .false., then system with integer spin
!    prxsec:  if .true., integral cross sections are printed to file 9
!    wrxsec:  if .true., integral cross sections are save in files 15,16,17 ..
!    twomol:  if .true., then molecule-molecule collision
!             if .false., then atom-molecule or molecule-surface collision
!    nlevel:  number of energetically distinct levels in channel basis
!    nlevop:  number of energetically distinct levels in channel basis which
!             are open asymptotically
!    numin:   initial value of the coupled-states projection index
!    numax:   final value of the coupled-states projection index
!    nud:     step in the coupled-states projection index
!    jlpar:   parity of channels
!               eps * (-1)**(j+l-jtot) = jlpar
!    nmax     on entry:  maximum row dimension of matrices
!    note!!!   if flaghf = .true.( see below), then the true values
!    of the rotational quantum numbers, the total angular momentum,
!    and the coupled-states projection index are equal to the values
!    stored in jlev, jtot, jfirst, jfinal, nu, numin, and numax plus 1/2

! ----------------------------------------------------------------------
use constants
use mod_hibrid2, only: print_integral_cross_sections
use mod_cosysi, only: ispar
use mod_basis, only: basis_get_isa
use mod_par, only: iprint
use funit
use mod_fileid, only: FILEID_SAV
use mod_savfile, only: REC_EN_START, EN_REC_COUNT, EN_REC_PRESENT_INTEGRAL_XS
use mod_parpot, only: potnam=>pot_name, label=>pot_label
use mod_selb, only: ibasty
use mod_ered, only: ered, rmu
use mod_surf, only: flagsu
use mod_hiutil, only: dater
implicit none
real(8), intent(out) :: zmat(nmax, nmax)
real(8), intent(out) :: tq3(nmx, nmx)
integer, intent(in) :: jlev(nmx)
real(8), intent(in) :: elev(nmx)
integer, intent(in) :: inlev(nmx)
integer, intent(in) :: nerg
real(8), intent(in) :: energ(nerg)
integer, intent(in) :: jfirst
integer, intent(in) :: jfinal
integer, intent(in) :: jtotd
logical, intent(in) :: csflag
logical, intent(in) :: flaghf
logical, intent(in) :: wrxsec
logical, intent(in) :: prxsec
logical, intent(in) :: ipos
logical, intent(in) :: twomol
logical, intent(in) :: nucros
integer, intent(in) :: nlevel
integer, intent(in) :: nlev(nerg)
integer, intent(in) :: numin
integer, intent(in) :: numax
integer, intent(in) :: nud
integer, intent(in) :: jlpar
integer, intent(in) :: nmax
integer, intent(in) :: nmx
logical, intent(in) :: ihomo

real(8) :: ener
integer :: i, ien, irec, isa, j, jhold, jj1, jj2, jmin, jphold, nlevmx, nlevop, nn, nxfile
character*20 cdate
integer :: sav_file = FILEID_SAV
!   econv is conversion factor from cm-1 to hartrees
!   xmconv is converson factor from amu to atomic units
nlevmx=0
do 10 ien=1,nerg
10 nlevmx=max(nlevmx,nlev(ien))
if (prxsec) then
  if (.not. twomol) then
    if (.not. flagsu) then
      write (9, 100)
100       format &
       (/' LEVEL LIST FOR INTEGRAL CROSS SECTIONS', &
        /'   N   J  INDEX    EINT(CM-1)  EKIN(CM-1)',/)
    else
      write (9, 105)
105       format &
       (/' LEVEL LIST FOR DEGENERACY AVERAGED', &
         ' TRANSITION PROBABILITIES', &
        /'   N   J  INDEX    EINT(CM-1)  EKIN(CM-1)',/)
    end if
    do 120  i = 1, nlevmx
      if (.not. flaghf .or. ibasty.eq.12) then
        write (9, 110) i, jlev(i), inlev(i), elev(i) * econv, &
                       (energ(ien)-elev(i)*econv,ien=1,nerg)
110         format (i4, i5, i6, f13.3,5f11.3)
      else
        write (9, 115) i, (jlev(i)+0.5), inlev(i), &
         elev(i)*econv, (energ(ien)-elev(i)*econv,ien=1,nerg)
115         format (i4, f5.1, i6, f13.3,5f11.3)
      end if
120     continue
  else
    write (9, 125)
125       format (/' LEVEL LIST FOR INTEGRAL CROSS SECTIONS', &
             /'   N   J1   J2  INDEX    EINT(CM-1)  EKIN(CM-1)'/)
      do 135  i = 1, nlevmx
        jj2 = mod( jlev(i), 10)
        jj1 = jlev(i) / 10
        write (9, 130) i, jj1, jj2, inlev(i), elev(i) * econv, &
          (energ(ien)-elev(i)*econv,ien=1,nerg)
130         format (i4, 2i5, i6, 2f13.3)
135       continue
  end if
  do 210  ien = 1, nerg
    ener = energ(ien)
    if (.not. flagsu) then
      if (.not. flaghf) then
        write (9, 140) ien, rmu * xmconv , ener, jfirst, &
                       jfinal, jtotd
140         format (/' ** INTEGRAL CROSS SECTIONS; IEN=', i2,' **', &
                /'    RMU=', f9.4, '  E=', f10.3,'  JTOT-1=', i3, &
                 '  JTOT-2=', i4,'  JTOT-D=', i3)
      else
        write (9, 145) ien, rmu * xmconv, ener, (jfirst+0.5), &
                       (jfinal+0.5), jtotd
145         format (/' ** INTEGRAL CROSS SECTIONS; IEN=', i2,' **', &
              /'    RMU=', f9.4, '  E=', f10.3,'  JTOT-1=', f5.1, &
                 '  JTOT-2=', f6.1,'  JTOT-D=', i3)
      end if
      if (csflag) then
        if (.not. flaghf) then
          write (9, 160) numin, numax, nud
160           format (29x, &
            '  NUMIN=', i2,'  NUMAX=', i3,'  NU-D=',i2)
        else
          write (9, 165) numin+0.5, numax+0.5, nud
165           format (29x, &
            '  NUMIN=', f5.1,'  NUMAX=', f6.1,'  NU-D=',i2)
        endif
      end if
    else
      if (.not. flaghf) then
        write (9, 180) ien, rmu * xmconv , ener, numin, numax, &
                       nud
180         format (/' ** SUMMED DEGENERACY AVERAGED TRANSITION', &
                 ' PROBABILITIES;  IEN=', i2,' **', &
                /'    RMU=', f9.4, '  E=', f8.3,'  M-MIN=', i3, &
                 '  M-MAX=', i4, '  M-STEP=', i2)
      else
        write (9, 185) ien, rmu * xmconv, ener, (numin+0.5), &
                       (numax+0.5), nud
185         format (/' ** SUMMED DEGENERACY AVERAGED TRANSITION', &
                 ' PROBABILITIES;  IEN=', i2,' **', &
              /'    RMU=', f9.4, '  E=', f8.3,'  M-MIN=', f5.1, &
                 '  M-MAX=', f6.1, '  M-STEP=', i2)
      end if
    end if
    irec=(ien-1)* EN_REC_COUNT + REC_EN_START
    if(nucros) irec = irec + nerg * EN_REC_COUNT
    nlevop=nlev(ien)
    jmin=10000
    do 190 i=1,nlevop
      if (jlev(i) .lt. jmin) jmin=jlev(i)
190     continue
    if ((jmin .le. numax .or. iprint.ge.2) .or. .not.csflag) &
        then
      call dbin(sav_file, irec + EN_REC_PRESENT_INTEGRAL_XS,jhold,jphold,nn,zmat,nmax,nlevop)
      call print_integral_cross_sections(9,zmat,nlevop,nmax,ipos,csflag,flaghf,twomol, &
                 numax,jlev,inlev)
    else
      write (6, 195) jmin, numax
      write (9, 195) jmin, numax
195       format (/'** MIN(J) =',i3,' .GT. NUMAX =',i3, &
          '; NO CS CROSS SECTIONS! ')
    endif
    write (9, 200)
200     format (1h ,30('='))
210   continue
end if
if (wrxsec) then
  do 300  ien = 1, nerg
    nxfile = FUNIT_ICS_START + ien - 1
!         rewind nxfile
    irec=(ien-1)* EN_REC_COUNT + REC_EN_START
    if(nucros) irec = irec + nerg * EN_REC_COUNT
    nlevop=nlev(ien)
    call dbin(sav_file, irec + EN_REC_PRESENT_INTEGRAL_XS,jhold,jphold,nn,zmat,nmax,nlevop)
    call dater (cdate)
    write (nxfile, 215) cdate
215     format (1x,a)
    write (nxfile, 215) label
    write (nxfile, 215) potnam
    write (nxfile, 230) energ(ien), rmu * xmconv
230     format (f10.3, f15.11)
!ABER additional information : IHOMO , ISA
    write (nxfile, 235) csflag, flaghf, flagsu, twomol, ihomo
235     format (5l3)
    isa = basis_get_isa(ibasty, ispar)
    write (nxfile, 240) jfirst, jfinal, jtotd, numin, numax, &
                        nud, jlpar, isa
240     format (24i5)
!ABER
   write (nxfile, 240)  nlevel, nlevop
    write (nxfile, 240) (jlev(i), inlev(i), i = 1, nlevel)
    write (nxfile, 245) (elev(i), i = 1, nlevel)
!245       format (8f16.9)
245     format (8(1pe15.8))
    do 250  i = 1, nlevop
      write (nxfile, 245) (zmat(j,i), j = 1, nlevop)
250     continue
    close (nxfile)
300   continue
end if
return
end
! ---------------------------------------------------------------------------
subroutine readpc (fname, a, scmat, nmax)
!  subroutine to write out selected partial cross sections
!  input:
!    PARTC,JOB,JINI,INDI,IEN,IPRINT
!  from file {fname1}.ics
!  author:  millard alexander
!  current revision date:  3-dec-2007 by mha
!  ------------------------------------------------------------------
!  variables in call list:
!    zmat:    on return:  contains the nlevop x nlevop matrix of integral
!                         cross sections
!    jlev:   rotational angular momenta of each energetically open level
!    elev:   energy (in hartree) of each energetically open level
!    inlev:  additonal quantum index for each energetically open level
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
!    jlpar:   parity of channels
!               eps * (-1)**(j+l-jtot) = jlpar
!             if jlpar = 0, then integral cross sections include contributions
!             of both parities
!    note!!!   if flaghf = .true.( see below), then the true values
!    of the rotational quantum numbers, the total angular momentum,
!    and the coupled-states projection index are equal to the values
!    stored in jlev, jtot, jfirst, jfinal, nu, numin, and numax plus 1/2
!    xmu:       collision reduced mass in (c12) atomic mass units
!    econv:     conversion factor from cm-1 to hartrees
!  ------------------------------------------------------------------
use mod_cosout, only: nnout, jout
use mod_coiout, only: niout, indout
use constants
use mod_cojhld, only: jlev => jhold ! jlev(1)
use mod_coisc1, only: inlev => isc1 ! inlev(1)
use mod_coisc2, only: nj, jlist => isc2 ! nj,jlist(1)
use mod_cosc1, only: elev => sc1 ! elev(1)
use mod_cosc2, only: csum => sc2 ! csum(1)
use mod_cosc3, only: tsum => sc3 ! tsum(1)
use mod_version, only : version
use mod_parpot, only: potnam=>pot_name, label=>pot_label
use mod_selb, only: ibasty
use mod_hiutil, only: gennam
implicit double precision (a-h,o-z)
character*(*) fname
character*20 cdate
character*40 xnam1, xnam2
character*80 line
logical csflag, flaghf, iprint, flagsu, twomol, existf, nucros
dimension  a(4),scmat(nmax,1)

!  input parameters
iprint=.true.
jini=nint(a(1))
indi=nint(a(2))
iene=nint(a(3))
if (a(4) .lt. 0.0) iprint =  .false.
if (iene .le. 0) iene=1
!  open integral cross section file
call gennam (xnam1, fname, iene, 'pcs', lenx)
inquire (file = xnam1, exist = existf)
if (.not. existf) then
  write (6, 10) xnam1(1:lenx)
10   format(/' Partial cross section file ',(a),' not found',/)
  return
end if
open (1, file = xnam1,status = 'old')
!  open output file for partial cross sections
call gennam (xnam2, fname, iene, 'psc', lenx)
call openf(3,xnam2,'sf',0)
call version(3)
read (1, 40) cdate
40 format (1x, a)
read (1, 40) label
read (1, 40) potnam
!  print job information
write (3, 50) xnam1, cdate, label, potnam
if (iprint) write (6, 50) xnam1, cdate, label, potnam
50 format(/' PARTIAL CROSS SECTIONS READ FROM FILE ',(a)/ &
        ' WRITTEN:    ',(a)/ &
        ' LABEL:      ',(a)/ &
        ' POT NAME:   ',(a) )
read (1, *) ered, rmu
read (1, *) csflag, flaghf, flagsu, twomol, nucros
read (1, *) jfirst, jfinal, jtotd,numin, numax, jlpar, nud
read (1, *) nj, nlev
read (1, *) (jlev(i), inlev(i), i = 1, nlev)
read (1, *) (elev(i), i = 1, nlev)
irow=0
if( .not. twomol) then
  write (3, 60) ered
  if(iprint) write (6, 60) ered
60   format (/' LEVEL LIST FOR PARTIAL CROSS SECTIONS AT E=', &
          f8.2,' CM-1'/'   N    J  INDEX  EINT(CM-1)',/)
  do 90  i = 1, nlev
    if(jlev(i).eq.jini.and.inlev(i).eq.indi) irow=i
    if (ibasty.ne.12) then
      if (.not. flaghf) then
        write (3, 70) i, jlev(i), inlev(i), elev(i) * econv
        if(iprint) &
          write (6, 70) i, jlev(i),inlev(i),elev(i)*econv
70         format (i4, i5, i6, f11.3)
      else
        write (3, 80) i, (jlev(i)+0.5d0), inlev(i), &
                     elev(i) * econv
        if(iprint) write (6, 80) i, (jlev(i)+0.5d0), inlev(i), &
                     elev(i) * econv
80         format (i4, f5.1, i6, f11.3)
      endif
    else
        write (3, 81) i, jlev(i), inlev(i)+0.5d0, &
                     elev(i) * econv
        if(iprint) write (6, 81) i,jlev(i),inlev(i)+0.5d0, &
                     elev(i) * econv
81         format (i4, i5,f6.1,f11.3)
    end if
90   continue
else
  write (9, 100) ered
  if(iprint) write (9, 100) ered
100   format (/' LEVEL LIST FOR PARTIAL CROSS SECTIONS AT E=', &
    f8.2,' CM-1'/'   N   J1   J2  INDEX  EINT(CM-1)'/)
  do 120  i = 1, nlev
    if(jlev(i).eq.jini.and.inlev(i).eq.indi) irow=i
    jj1 = jlev(i) / 10
    jj2 = mod(jlev(i), 10)
    write (3, 110) i, jj1, jj2, inlev(i), elev(i) * econv
    if(iprint) write (9, 110) i, jj1, jj2,inlev(i),elev(i)*econv
110     format (i4, 2i5, i6, f11.3)
120   continue
end if
if(irow.eq.0) goto 180
if(flaghf) then
  write(3,130) jini+0.5d0,indi,ered-elev(irow)*econv
  if(iprint) write(6,130) jini+0.5d0,indi,ered-elev(irow)*econv
130   format(/1x,'PARTIAL CROSS SECTIONS FOR INITIAL STATE J=',f5.1, &
             '  INDEX=',i5,'  EKIN=',f10.2)
else
  write(3,131) jini,indi,ered-elev(irow)*econv
  if(iprint) write(6,131) jini,indi,ered-elev(irow)*econv
131   format(/1x,'PARTIAL CROSS SECTIONS FOR INITIAL STATE J=',i3, &
             '  INDEX=',i5,'  EKIN=',f10.2)
end if
njj=0
do 140 j=1,iabs(nnout)
  jj=jout(j)
  in=0
  do 140 i=1,max(1,niout)
    in=indout(i)
    do 140 k=1,nlev
      if(niout.eq.0) in=inlev(k)
      if(jlev(k).eq.jj.and.inlev(k).eq.in) then
        njj=njj+1
        jlist(njj)=k
        csum(njj)=0
      end if
140 continue
write(3,150) (jlist(j),j=1,njj)
if(iprint) then
  if (.not.nucros) then
    write(6,150) (jlist(j),j=1,njj)
150     format(1x,'JTOT',t10,'FINAL STATES'/(t1,10i11))
  else
    write(6,151) (jlist(j),j=1,njj)
151     format(2x,'NU',t10,'FINAL STATES'/(t1,10i11))
  endif
endif
wt=1d0
if (jtotd.gt.1) wt=0.5d0
160 read (1, '(a)',end=210) line

!.....here if RESTART-calculation has been done
if (line(1:14).eq.' ** RESTART **') then
  read (1, 40) cdate
  read (1, 40) label
  read (1, 40) potnam
  read (1, *) ered, rmu
  read (1, *) csflag, flaghf, flagsu, twomol, nucros
  read (1, *) jfirst, jfinal, jtotd,numin, numax, jlpar, nud
  read (1, *) nj, nlev
  read (1, *) (jlev(i), inlev(i), i = 1, nlev)
  read (1, *) (elev(i), i = 1, nlev)
  read (1, '(a)') line
end if
!.....end RESTART handling
read (line, '(i8)') jtot
! loop over initial states as specified in jout
irow=0
do 170 i=1,nj
  read (1,*) ji, in, (scmat(i,j), j = 1, nlev)
  if(ji.eq.jini.and.in.eq.indi) irow=i
170 continue
180 if(irow.eq.0) then
  write(3,190) jini,indi
  if(iprint) write(6,190) jini,indi
190   format(/' INITIAL STATE J=',i3,'  INDEX=',i3,' NOT FOUND')
  goto 210
end if
write(3,200) jtot,(scmat(irow,jlist(j)),j=1,njj)
if(iprint) write(6,200) jtot,(scmat(irow,jlist(j)),j=1,njj)
!     write (6,*) 'jfirst, jfinal:  ', jfirst, jfinal
do 195 j=1,njj
!       wt=1d0
!       if (jtotd.gt.1) then
!          if (jtot.eq.jfirst .or. jtot.eq.jfinal) wt=0.5d0
!       endif
  tsum(j) = scmat(irow, jlist(j))
  csum(j) = csum(j) + wt*scmat(irow,jlist(j))
195 continue
200 format(1x,i3,(t5,10(1pd11.4)))
wt=1d0
goto 160
! if last jtot, and jtotd>1, subtract off 50% of last value, consistent with
! trapezoidal rule interpolation
do 205 j=1,njj
   csum(j)=csum(j)-05d0*tsum(j)
205 continue
210 continue
!      print *, 'jtotd', jtotd
call dscal(njj,dble(jtotd),csum,1)
if(iprint) then
   write(6,220) (csum(j),j=1,njj)
else
   write(6,225) (csum(j),j=1,njj)
endif
write(3,220) (csum(j),j=1,njj)
220 format(/1x,'SUM:',(t5,10(1pd11.4)))
225 format(/1x,'SUM OF PARTIAL CROSS SECTIONS:  ', &
 /5x,(t5,10(1pd11.4)))
close(1)
close(3)
return
end
! ---------------------------------------------------------------------------
subroutine intchk(irec,q,q1,q2,q3,jf,jp,nmax,n,nucros)
! ---------------------------------------------------------------------------
!
!  check consistency of stored integral and partial cross sections
!
!  latest revision date: 18-may-1997 by mha
!
! ---------------------------------------------------------------------------
!
use mod_savfile, only: EN_REC_PRESENT_INTEGRAL_XS, EN_REC_PREVIOUS_INTEGRAL_XS, EN_REC_PRESENT_PARTIAL_XS, EN_REC_PREVIOUS_PARTIAL_XS, EN_REC_2NDLAST_PARTIAL_XS
use mod_fileid, only: FILEID_SAV
implicit none
integer, intent(in) :: irec  ! the base (minimum of the record ranges for an energy) record index to use to read the cross sections in the sav_file.
real(8), intent(out) :: q(nmax,n)
real(8), intent(out) :: q1(nmax,n)
real(8), intent(out) :: q2(nmax,n)
real(8), intent(out) :: q3(nmax,n)
integer, intent(out) :: jf
integer, intent(out) :: jp
integer, intent(in) :: nmax
integer, intent(in) :: n
logical, intent(out) :: nucros
integer :: sav_file = FILEID_SAV
integer :: i, j
integer :: jl, jl1, jl2, jl3
integer :: jpl, jp1, jp2, jp3
integer :: next
integer :: nn
integer :: ierr
integer :: jlev(n)
!
write(6,5)
write(9,5)
5 format(/' CHECKING INTERPOLATION DATA IN RESTART FILE:'/)
call dbin(sav_file, irec + EN_REC_PREVIOUS_INTEGRAL_XS,jl,jpl,next,q,nmax,n)
write(6,6) jl,jpl
write(9,6) jl,jpl
6 format(' READ INTEGRAL CROSS SECTIONS FOR J=',i3,'  JP=',i2)
call dbin(sav_file, irec + EN_REC_PRESENT_PARTIAL_XS,jl3,jp3,nn,q3,nmax,n)
write(6,7) jl3,jp3
write(9,7) jl3,jp3
7 format(' READ PARTIAL  CROSS SECTIONS FOR J=',i3,'  JP=',i2)
if(next.eq.1) then
  do i=1,n
    do j=1,n
      q(j,i) = q(j,i) + q3(j,i)
    end do
  end do
  goto 15
else if(next.eq.2) then
  call dbin(sav_file, irec + EN_REC_PREVIOUS_PARTIAL_XS,jl2,jp2,nn,q2,nmax,n)
  write(6,7) jl2,jp2
  write(9,7) jl2,jp2
  ASSERT( .not. nucros )
  call intpl2(jl,jl2,jl3,nmax,n,q,q2,q3,jlev,nucros)
  goto 15
else if(next.ge.3) then
  call dbin(sav_file, irec + EN_REC_PREVIOUS_PARTIAL_XS,jl2,jp2,nn,q2,nmax,n)
  write(6,7) jl2,jp2
  write(9,7) jl2,jp2
  call dbin(sav_file, irec + EN_REC_2NDLAST_PARTIAL_XS,jl1,jp1,nn,q1,nmax,n)
  write(6,7) jl1,jp1
  write(9,7) jl1,jp1
  ASSERT( .not. nucros )
  call intpl3(jl,jl1,jl2,jl3,nmax,n,q,q1,q2,q3,jlev,nucros)
end if
15 call dbin(sav_file, irec + EN_REC_PRESENT_INTEGRAL_XS,jf,jp,nn,q1,nmax,n)
write(6,6) jf,jp
write(9,6) jf,jp
if(jl3.ne.jf.or.jp3.ne.jp) then
  write(6,20) jl3,jp3,jf,jp
  write(9,20) jl3,jp3,jf,jp
20   format(/' ERROR DETECTED IN INTCHK: JL3,JP3,JF,JP:',4i4)
  call exit
end if
ierr=0
do i=1,n
  do j=1,n
    if(abs(q(j,i)-q1(j,i)) > 1.d-10) ierr=ierr+1
  end do
end do
if(ierr.ne.0) then
  write(6,40) ierr
  write(9,40) ierr
40   format(/1x,i6,' ERRORS DETECTED IN INTCHK. ABORT')
  call exit
end if
write(6,50)
write(9,50)
50 format(/' NO ERRORS DETECTED.')
return
end
! --------------------------------------------------------------------
subroutine restrt (jtot,jtop,jtotd,jlpar,nu,nutop,nud,nerg,nlev, &
               nchmax,rtmn1,rtmx1,dinsid, writs, csflag, nucros)
!  subroutine to read restart information
!  author:  h.-j. werner
!  latest revision:  17-may-1997 by mha
! --------------------------------------------------------------------
use mod_coamat, only: q3 => amat ! q3(1)
use mod_cobmat, only: q => bmat ! q(1)
use mod_coener, only: energ, max_en
use mod_cow, only: q1 => w_as_vec ! q1(1)
use mod_cozmat, only: q2 => zmat_as_vec ! q2(1)
use mod_par, only: wrpart, wrxsec
use mod_parpot, only: label=>pot_label
use mod_fileid, only: FILEID_SAV
use mod_savfile
implicit none
integer, intent(out)   :: jtot
integer, intent(out)   :: jtop
integer, intent(inout) :: jtotd
integer, intent(out)   :: jlpar
integer, intent(out)   :: nu
integer, intent(out)   :: nutop
integer, intent(inout) :: nud
integer, intent(out)   :: nerg
integer, intent(out)   :: nlev(max_en)
integer, intent(out)   :: nchmax
real(8), intent(out)   :: rtmn1
real(8), intent(out)   :: rtmx1
real(8), intent(out)   :: dinsid
logical, intent(out)   :: writs
logical, intent(out)   :: csflag
logical, intent(out)   :: nucros
character*40 oldlab
real(8), allocatable :: word(:)
integer, allocatable :: iword(:)
integer :: i, j, ien, irec, jf, jp, n
allocate(word(get_sav_real_params_count()))
allocate(iword(get_sav_int_params_count()))
if (.not. wrpart .and. .not. wrxsec) then
  write (6, 5)
  write (9, 5)
5   format (' EITHER WRPART OR WRXSEC MUST HAVE BEEN TRUE IN', &
    ' ORIGINAL DATA',/,'    FOR RESTART TO WORK; ABORT')
  call exit
endif
call drest(FILEID_SAV)
call dbri(iword, get_sav_int_params_count(), FILEID_SAV, REC_PARAMS)
j=1
do i = SAV_INTPARAMS_LABEL, SAV_INTPARAMS_LABEL + NUM_INTS_IN_LABEL - 1
  write(oldlab(j:j + NUM_CHARS_IN_INT - 1),'(A4)') iword(i)
  j = j + NUM_CHARS_IN_INT
end do
if (oldlab .ne. label) then
    write (9, 50) oldlab, label
    write (6, 50) oldlab, label
50     format( &
  /' *** LABEL IN FILE TRSTRT DOES NOT MATCH INPUT DATA', &
     '; ABORT ***',/'     TRSTRT:', a, /, &
                        '      INPUT:', a)
    call exit
end if
jtot=iword(SAV_INTPARAMS_JTOT)
if(iword(SAV_INTPARAMS_JTOTD).ne.jtotd) then
    write(6,60) jtotd, iword(SAV_INTPARAMS_JTOTD)
60     format(/' RESTART: JTOTD=',i3,' RESET TO ',i3)
end if
jtotd=iword(SAV_INTPARAMS_JTOTD)
jlpar=iword(SAV_INTPARAMS_JLPAR)
nu=iword(SAV_INTPARAMS_NU)
nutop=iword(SAV_INTPARAMS_NUTOP)
if(iword(SAV_INTPARAMS_NUD).ne.nud) then
    write(6,70) nud, iword(SAV_INTPARAMS_NUD)
70     format(/' RESTART: NUD=',i3,' RESET TO ',i3)
end if
nud=iword(SAV_INTPARAMS_NUD)
nerg=iword(SAV_INTPARAMS_NERG)
nchmax=iword(SAV_INTPARAMS_NCHMAX)
csflag=iword(SAV_INTPARAMS_CSFLAG).ne.0
writs =iword(SAV_INTPARAMS_WRITS).ne.0
nucros=iword(SAV_INTPARAMS_NUCROS).ne.0
jtop=iword(SAV_INTPARAMS_JTOP)
call dbrr(word, get_sav_real_params_count(), FILEID_SAV, REC_LAST_USED)
rtmn1=word(SAV_REALPARAMS_RTMN1)
rtmx1=word(SAV_REALPARAMS_RTMX1)
dinsid=word(SAV_REALPARAMS_DINSID)
write(6,75) jtot,jlpar,nu
write(9,75) jtot,jlpar,nu
75 format(/' RESTART DATA FOUND FOR JTOT=',i3,'  JLPAR=',i2, &
          '  NU=',i3)
do ien=1,nerg
  energ(ien)=word(SAV_REALPARAMS_ENERG_START + ien - 1)
  nlev(ien)=iword(SAV_INTPARAMS_NLEV_START + ien - 1)
  n=nlev(ien)
  irec=(ien-1)* EN_REC_COUNT + REC_EN_START
  if(nucros) irec = irec + EN_REC_COUNT * nerg
  write(6,80) energ(ien),nlev(ien)
  write(9,80) energ(ien),nlev(ien)
80 format(/' ENERGY=',f10.3,' NLEVOP=',i3)
  call intchk(irec,q,q1,q2,q3,jf,jp,n,n,nucros)
  if((.not.nucros.and.(jf.ne.jtot.or.jp.ne.jlpar)) &
     .or. (nucros.and.jf.ne.nu)) then
    write(6,85)
    write(9,85)
85     format(/' INCONSISTENT DATA DETECTED. ABORT')
    call exit
  end if
end do
deallocate(iword)
deallocate(word)
return
end subroutine

subroutine rsave(jtot,jtop,jtotd,jlpar,nu,nutop,nud,nerg,nlev, &
       nchmax,rtmn1,rtmx1,dinsid,writs, csflag, nucros)
use mod_coener, only: energ, max_en
use mod_parpot, only: label=>pot_label
use mod_file, only: savfil
use funit, only: FUNIT_SAV
use mod_fileid, only: FILEID_SAV
use mod_savfile
implicit none
integer, intent(in)   :: jtot
integer, intent(in)   :: jtop
integer, intent(in) :: jtotd
integer, intent(in)   :: jlpar
integer, intent(in)   :: nu
integer, intent(in)   :: nutop
integer, intent(in)   :: nud
integer, intent(in)   :: nerg
integer, intent(in)   :: nlev(max_en)
integer, intent(in)   :: nchmax
real(8), intent(in)   :: rtmn1
real(8), intent(in)   :: rtmx1
real(8), intent(in)   :: dinsid
logical, intent(in)   :: writs
logical, intent(in)   :: csflag
logical, intent(in)   :: nucros
real(8), allocatable :: word(:)
integer, allocatable :: iword(:)
integer :: i, j, ien
allocate(word(get_sav_real_params_count()))
allocate(iword(get_sav_int_params_count()))

j=1
!.....reserve space for restart information
call dres(64, FILEID_SAV, 1)
do i = SAV_INTPARAMS_LABEL, SAV_INTPARAMS_LABEL + NUM_INTS_IN_LABEL - 1
  read(label(j:j + NUM_CHARS_IN_INT - 1),'(A4)') iword(i)
  j = j + NUM_CHARS_IN_INT
end do
iword(SAV_INTPARAMS_JTOT)=jtot
iword(SAV_INTPARAMS_JTOTD)=jtotd
iword(SAV_INTPARAMS_JLPAR)=jlpar
iword(SAV_INTPARAMS_NU)=nu
iword(SAV_INTPARAMS_NUTOP)=nutop
iword(SAV_INTPARAMS_NUD)=nud
iword(SAV_INTPARAMS_NERG)=nerg
iword(SAV_INTPARAMS_NCHMAX)=nchmax
iword(SAV_INTPARAMS_CSFLAG)=0
iword(SAV_INTPARAMS_WRITS)=0
iword(SAV_INTPARAMS_NUCROS)=0
iword(SAV_INTPARAMS_JTOP)=jtop
if(csflag) iword(SAV_INTPARAMS_CSFLAG)=1
if(writs)  iword(SAV_INTPARAMS_WRITS)=1
if(nucros) iword(SAV_INTPARAMS_NUCROS)=1
word(SAV_REALPARAMS_RTMN1)=rtmn1
word(SAV_REALPARAMS_RTMX1)=rtmx1
word(SAV_REALPARAMS_DINSID)=dinsid
do ien = 1, nerg
  word(SAV_REALPARAMS_ENERG_START + ien - 1) = energ(ien)
  iword(SAV_INTPARAMS_NLEV_START + ien - 1) = nlev(ien)
end do
call dbwi(iword, get_sav_int_params_count(), FILEID_SAV, REC_PARAMS)
call dbwr(word, get_sav_real_params_count(), FILEID_SAV, REC_LAST_USED)
call dbwc(FILEID_SAV, 1)
call dclos(FILEID_SAV)
call dopen(FILEID_SAV, FUNIT_SAV, savfil)
deallocate(iword)
deallocate(word)
return
end
! ------------------------------------------------------------------
subroutine intcrs(filnam,a)
!
! driver subroutine to calculate integral cross sections
! from s-matrix elements
!
! author: hjw
! revision date: 22-jan-2008 by mha
! half-integral j's now printed out.  revision date:  28-nov-2011
!    (p.j.dagdigian)
!
! current revision date: 23-feb-2013 by p. dagdigian
!
! WARNING: Due to a revision of the partcr subroutine in may-30-2013,
!     flaghf no longer applies to the second molecule (j2) in
!     molecule-molecule collisions.
! ------------------------------------------------------------------
use mod_codim, only: mmax
use mod_cosout
use constants
use mod_cotq1, only: sc1 => scmat ! sc1(1)
use mod_cotq2, only: sc2 => scmat ! sc2(1)
use mod_cotq3, only: sc3 => scmat ! sc3(1)
use mod_coamat, only: sc4 => toto ! sc4(1)
use mod_coisc2, only: nj, jlist => isc2 ! nj,jlist(1)
use mod_coisc6, only: isc1 => isc6 ! isc1(1)
use mod_coisc7, only: isc2 => isc7 ! isc2(1)
use mod_cosc1, only: elev => sc1 ! elev(1)
use mod_cosc2, only: inlev => sc2int ! inlev(1)
use mod_cosc3, only: jlev => sc3int ! jlev(1)
use mod_coz, only: sreal => z_as_vec ! sreal(1)
use mod_cow, only: simag => w_as_vec ! simag(1)
use mod_cozmat, only: sigma => zmat_as_vec ! sigma(1)
use mod_hibrid2, only: mxoutr
use mod_par, only: batch, ipos
use mod_parpot, only: label=>pot_label
use mod_selb, only: ibasty
use mod_ered, only: ered, rmu
use mod_hiutil, only: gennam
use mod_hiutil, only: mtime, gettim
use mod_hismat, only: rdhead
use mod_fileid, only: FILEID_TMP
implicit double precision (a-h,o-z)
character*(*) filnam
character*40  icsfil, smtfilnam, xname
integer :: smt_unit = 1
integer :: ics_unit = 2
integer :: tmp_unit = 3
integer :: tmp_file = FILEID_TMP
character*20  cdate
character*10  elaps, cpu
character*13  string
logical csflag, flaghf, flagsu, twomol, exstfl, &
        nucros, notequ
!/ nnout, jout(21)
dimension a(3)
!
! initialize timer
call mtime(cpu0,ela0)
! input
eredsv=ered
rmusav=rmu
ien = a(1)
maxjtot = a(2)
if (ien.le. 0) ien = 1
!
! generate filename and check if it is present
!
call gennam(smtfilnam,filnam,ien,'smt',lenfs)
inquire(file = smtfilnam, exist = exstfl)
if (.not. exstfl) then
    write(6,10) smtfilnam(1:lenfs)
10     format(' ** FILE ',(a),' NOT FOUND **')
    return
end if
! open file for interpolation
call dinit
call gennam (xname, filnam, 0, 'tmp', lenx)
call dopen(tmp_file,tmp_unit,xname)
!
! open smatrix-file
!
call openf(smt_unit,smtfilnam,'tu',0)
!
! open file for integral cross sections
!
call gennam(icsfil,filnam,ien,'xxsc',lenft)
call openf(ics_unit,icsfil,'sf',0)
!
! read header of smatrix-file
!
call rdhead(smt_unit,cdate, &
   ered,rmu,csflag,flaghf,flagsu,twomol, &
   nucros,jfirst,jfinal,jtotd,numin,numax,nud,nlevel,nlevop, &
   nnout,jlev,inlev,elev,jout)
maxjt=jfinal
if(maxjtot.gt.0) maxjt=min(maxjtot,jfinal)
spin=0
if(flaghf) spin=0.5d0
!	
write (ics_unit,15) smtfilnam,cdate,label,maxjt
if (.not.batch)write (6,15) smtfilnam,cdate,label,maxjt
15 format( &
/' INTEGRAL CROSS SECTIONS CALCULATED FROM S-MATRICES ON FILE ', &
(a)/ &
        ' WRITTEN:   ',(a)/ &
        ' LABEL:     ',(a)/ &
        ' JTOT_MAX:  ', i3)

! write a header
if (.not.flaghf) then
  write(ics_unit,20) ien,rmu*xmconv,ered*econv,jfirst,maxjt,jtotd
  if(.not.batch) &
  write(6,20) ien,rmu*xmconv,ered*econv,jfirst,maxjt,jtotd
20    format (' IEN=', i2,' RMU=', f9.4,' E=', f9.2, &
          ' JTOT-1=', i3,' JTOT-2=', i4,' JTOT-D=', i3,/)
else
  write(ics_unit,21) &
     ien,rmu*xmconv,ered*econv,jfirst+0.5,maxjt+0.5,jtotd
  if(.not.batch) &
  write(6,21) &
     ien,rmu*xmconv,ered*econv,jfirst+0.5,maxjt+0.5,jtotd
21    format (' IEN=', i2,' RMU=', f9.4,' E=', f9.2, &
          ' JTOT-1=', f5.1,' JTOT-2=', f6.1,' JTOT-D=', i3)
endif
!     stop 'in intcrs'
! make a list of pointers
nj = 0
do 120 i=1, iabs(nnout)
   jo = jout(i)
   do 100 j=1, nlevop
     j1 = jlev(j)
     if(j1.ne.jo) goto 100
     nj = nj + 1
     jlist(nj) = j
100    continue
120 continue
! check if there had been any match
if(nj.eq.0) then
   write(ics_unit,130)
   if(.not. batch) write(6,130)
130    format(' *** NO TRANSITIONS FOUND, ABORT ***')
   goto 300
end if
! now compute cross sections
call intcr(csflag,flaghf,twomol,flagsu,nucros, &
            numin,numax,nud,jfirst,jfinal,jtotd,maxjt, &
            sigma,sreal,simag,sc1,sc2,sc3,sc4,nlevop,mmax,tmp_file)
string=' '
if(nnout.lt.0) string='(COLUMNS)'
write (ics_unit, 210) string
if(.not. batch) write (6, 210) string
210 format (/,' LEVEL LIST FOR INTEGRAL CROSS SECTIONS ',a, &
        /,'   N     J   INDEX  EINT(cm-1)')
do 250 i = 1, nj
  jj = jlist(i)
  if (flaghf .and. ibasty.ne.12) then
    write (ics_unit, 220) i, jlev(jj)+spin, inlev(jj), elev(jj)*econv
    if(.not. batch) &
    write (6, 220) i, jlev(jj)+spin, inlev(jj), elev(jj)*econv
220     format (i4, 1x, f5.1, i6, f11.3)
  else
    write (ics_unit, 221) i, jlev(jj), inlev(jj), elev(jj)*econv
    if(.not. batch) &
    write (6, 221) i, jlev(jj), inlev(jj), elev(jj)*econv
221     format (i4, 1x, i5, i6, f11.3)
  endif
250 continue
if(nnout.lt.0) then
write (ics_unit, 210) '(ROWS)'
  if(.not. batch) write (6, 210) '(ROWS)'
  do 255 i=1,nlevop
    write (ics_unit, 220) i, jlev(i)+spin, inlev(i), elev(i)*econv
    if(.not. batch) &
    write (6, 220) i, jlev(i)+spin, inlev(i), elev(i)*econv
255   continue
end if
if(.not.csflag) then
   write(ics_unit,260)
   if(.not. batch) write(6,260)
260    format(/,' CC INTEGRAL CROSS SECTIONS')
else
   write(ics_unit,270)
   if(.not. batch) write(6,270)
270    format(/,' CS INTEGRAL CROSS SECTIONS')
   if (.not. flaghf) then
     write (ics_unit, 280) numin, numax, nud
     if (.not.batch) write (6, 280) numin, numax, nud
280      format (' ** CS CALCULATION, NUMIN=', i2,', NUMAX=', &
               i2,' NUD=', i2, ' **')
   else
     write (ics_unit, 290) numin + 0.5, numax + 0.5, nud
     if (.not.batch) write (6, 290) numin+0.5, numax+ 0.5, nud
290      format (' ** CS CALCULATION, NUMIN=', f4.1, ', NUMAX=', &
            f4.1,' NUD=', i2, ' **')
   endif
end if
write (ics_unit, 295)
if (.not.batch) write (6, 295)
295   format (/' ** COLUMN HEADINGS ARE FINAL STATES, ROW', &
        ' HEADINGS ARE INITIAL STATES **')
ncol=nj
nrow=nj
if(nnout.lt.0) nrow=nlevop
call mxoutr (ics_unit, sc1, nrow, ncol, mmax, 0, ipos)
! use notequ as scratch variable (ipos = .false. for screen output)
notequ = .false.
if(.not. batch) call mxoutr (6, sc1, nrow, ncol, mmax, 0, notequ)
300 call closf(tmp_file)
close (ics_unit)
!      call dclos(1)
close(unit=smt_unit)
close(unit=tmp_unit,status='delete')
call mtime(cpu1,ela1)
cpu1 = cpu1 - cpu0
ela1 = ela1 - ela0
call gettim(cpu1,cpu)
call gettim(ela1,elaps)
if(.not. batch) write(6,400) elaps, cpu
400 format(/,' ** INTCRS FINISHED. ELAPSED TIME: ',a, &
         '  CPU TIME: ',a,' **')
ered=eredsv
rmu=rmusav
return
end
! ----------------------------------------------------------------------
subroutine intcr(csflag,flaghf,twomol,flagsu,nucros, &
           numin,numax,nud,jfirst,jfinal,jtotd,maxjt, &
           sigma,sreal,simag,sc1,sc2,scmat,tsq,nlevop,nmax,tmp_file)
!
! subroutine to calculate integral cross sections from s-matrix
! elements
!
! author: hjw with revisions by F. de Weerd and mha
! extended to molecule-molecule calculations (p.dagdigian 24-jan-2012)
! current revision date: 8-oct-2012 by q. ma
!
! ----------------------------------------------------------------------
use mod_cosout, only: nnout, jout
use mod_coisc2, only: nj, jlist => isc2 ! nj,jlist(1)
use mod_coisc6, only: isc1 => isc6 ! isc1(1)
use mod_coisc7, only: isc2 => isc7 ! isc2(1)
use mod_cosc1, only: elev => sc1 ! elev(1)
use mod_cosc2, only: inlev => sc2int ! inlev(1)
use mod_cosc3, only: jlev => sc3int ! jlev(1)
use mod_par, only: batch, ipos
use mod_selb, only: ibasty
use mod_hismat, only: sread
use mod_hitypes, only: bqs_type
implicit double precision (a-h,o-z)
logical, intent(in) :: csflag
logical, intent(in) :: flaghf
logical, intent(in) :: twomol
logical, intent(in) :: flagsu
logical, intent(in) :: nucros
integer, intent(in) :: numin
integer, intent(in) :: numax
integer, intent(in) :: nud
integer, intent(in) :: jfirst
integer, intent(in) :: jfinal
integer, intent(in) :: jtotd
integer, intent(in) :: maxjt
real(8), intent(out) :: sigma(nmax, nlevop)
real(8), intent(out) :: sreal(nmax, nlevop)
real(8), intent(out) :: simag(nmax, nlevop)
real(8), intent(out) :: sc1(nmax, nlevop)
real(8), intent(out) :: sc2(nmax, nlevop)
real(8), intent(out) :: scmat(nmax, nlevop)
real(8), intent(out) :: tsq(nmax, nlevop)
integer, intent(in) :: nlevop
integer, intent(in) :: nmax
integer, intent(in) :: tmp_file

type(bqs_type) :: row_bqs
type(bqs_type) :: packed_bqs

! clear sigma array
one=1.0d0
zero=0.0d0

sigma(1:nlevop, 1:nlevop) = 0.d0

jlpold=0
iaddr = 0
!
! read s-matrix for present jtot, nu
!
10 nopen = 0
call sread ( iaddr, sreal, simag, jtot, jlpar, nu, &
  row_bqs, packed_bqs, &
  1, nmax, nopen, ierr)
if(jlpold.eq.0) jlpold=jlpar
if(ierr.eq.-1) goto 100
if (csflag .and. (jtot.gt.maxjt)) goto 100
if (.not.csflag .and. (jtot.gt.maxjt)) then
    if (jlpar.eq.1) goto 10
    if (jlpar.eq.-1) goto 100
endif
if(ierr.lt.-1) then
  write(2,20)
  if(.not.batch) write(6,20)
20   format(' *** READ ERROR, ABORT')
  return
end if
! reset iaddr to 0 to insure sequential read
! print our message only for every 10th jtot
if (jtot .eq. 10*(jtot/10)) then
  if (csflag) then
    write(2,25) jtot,nu
    if (.not.batch) write(6,25) jtot,jlpar,nu
25     format(' READ S-MATRIX FOR JTOT=',i3,' NU=',i3)
  else
    write(2,26) jtot,jlpar
    if (.not. batch)write(6,26) jtot,jlpar
26     format(' READ S-MATRIX FOR JTOT=',i3,'  JLPAR=',i2)
  endif
endif
iaddr = 0
! calculate squared t-matrix
call tsqmat(tsq,sreal,simag,row_bqs, &
   packed_bqs,nopen,nmax)
! calculate partial cross sections
call partcr(tsq,sc1,nopen, &
            row_bqs, packed_bqs, &
            inlev, jlev, elev, jtot, nu, &
            csflag,flaghf,twomol,flagsu, &
            nlevop,nmax)
if(.not.csflag.or.nucros.or.(csflag.and.nu.eq.numin)) then
! sum up partial cross sections over nu
  scmat(1:nlevop, 1:nlevop) = sc1(1:nlevop, 1:nlevop)
else if(nu.gt.numin) then
  scmat(1:nlevop, 1:nlevop) = scmat(1:nlevop, 1:nlevop) + sc1(1:nlevop, 1:nlevop)
end if
! sum up partial cross sections to get integral cross sections
if(csflag.and.nucros) then
  irec=2
  call intpol(irec,jtot,jfirst,jfinal,jtotd,jlpar,jlpold, &
              jlev,nmax,nlevop,tsq,scmat,nucros,tmp_file)
  if(jtot+jtotd.gt.jfinal) then
    irec=7
    call intpol(irec,nu,numin,numax,nud,jlpar,jlpar, &
                jlev,nmax,nlevop,sigma,tsq,nucros, tmp_file)

  end if
else
  if (.not. csflag .or. (csflag .and. (nu .eq. numax) ) ) then
    irec=2
    call intpol(irec,jtot,jfirst,jfinal,jtotd,jlpar,jlpold, &
               jlev,nmax,nlevop,sigma,scmat,nucros,tmp_file)
  end if
end if
! loop back for next partial wave
goto 10
! make a compressed matrix of cross sections
100 continue

if(nnout.gt.0) then
  do i=1,nj
    ii=jlist(i)
    do j=1,nj
      sc1(j,i)=sigma(jlist(j),ii)
    end do
  end do
else
  do i=1,nj
    ii=jlist(i)
    do j=1,nlevop
      sc1(j,i)=sigma(j,ii)
    end do
  end do
end if

return
end
! ----------------------------------------------------------------------
subroutine tsqmat(tsq,sreal,simag,row_bqs, &
          col_bqs,nopen,nmax)
! ----------------------------------------------------------------------
!
!  routine to compute modulus squared t-matrix from given s-matrix
!  current revsion date:  24-jan-2012 by p.dagdigian
!
! ----------------------------------------------------------------------
use mod_hibasis, only: is_j12
use mod_selb, only: ibasty
use mod_hitypes, only: bqs_type
use constants, only: zero, one
implicit none
real(8), intent(out) :: tsq(nmax, nmax)
real(8), intent(in) :: sreal(nmax, nmax)
real(8), intent(in) :: simag(nmax, nmax)
type(bqs_type), intent(in) :: row_bqs
type(bqs_type), intent(in) :: col_bqs
integer, intent(in) :: nopen
integer, intent(in) :: nmax

integer :: icol, irow, in1, j1, l1, in2, j2, l2, j121, j122
complex*8 :: t
real(8) :: t2
logical :: diag
!
do icol = 1, col_bqs%length
   in1 = col_bqs%inq(icol)
   j1 = col_bqs%jq(icol)
   l1 = col_bqs%lq(icol)
   if (is_j12(ibasty)) then
     j121 = col_bqs%j12(icol)
   end if
   do irow = 1, nopen
      in2 = row_bqs%inq(irow)
      j2 = row_bqs%jq(irow)
      l2 = row_bqs%lq(irow)
      if (is_j12(ibasty)) j122 = row_bqs%j12(irow)
      diag = j1.eq.j2 .and. in1.eq.in2 .and. l1.eq.l2
      if (is_j12(ibasty)) diag = diag .and. j121.eq.j122
!
! convert s-matrix to t-matrix: t(j,j1,in1,l1 ; j2,in2,l2) =
!     delta(j1,in1,l1 ; j2,in2,l2) - s(j,j1,in1,l1 ; j2,in2,l2)
!
      t = -cmplx(sreal(irow,icol),simag(irow,icol))
      if (diag) t = t + cmplx(one,zero)
      t2 = real(t * conjg(t))
      tsq(irow,icol) = t2
   end do
end do
return
end
! ----------------------------------------------------------------------
subroutine partcr (tsq,  scmat, nopen, &
                   row_bqs, col_bqs, &
                   inlev, jlev, elev, jtot, nu, &
                   csflag, flaghf,twomol,flagsu, &
                   nlevop,nmax)
! ----------------------------------------------------------------------
!  this routine computes partial cross sections from squared t matrix
!  row_bqs: row indices of t-matrix (nopen values)
!  col_bqs: column indices of t-matrix (ncol values)
!  inlev, jlev: quantum numbers of asymptotic states (nlevop values)
!  elev: energy levels of asymptotic states (nlevop values)
!
!  current revision:  20-oct-2014 by p. dagdigian
!
!  revision:  30-may-2013 by q. ma
!     WARNING: starting from this revision, flaghf no longer applies
!     to the second molecule (j2) if twomol is set true.
!  latest revision:  include correct degeneracy factor, (2*xrow_bqs%jq1+1)*2,
!  [2nd factor is 2*s2+1)] in denominator for ibasty=23 (3P + 2S atom-atom)
!
!  revision:  15-aug-2016 by p.dagdigian
!     corrected degeneracy factor for j2 for ibasty = 12, 13, and 15
!  revision:  30-aug-2016 by p. dagdigian
!     corrected degeneracy factor for ibasty = 23
!
! ----------------------------------------------------------------------
use constants
use mod_hibasis, only: is_j12
use mod_selb, only: ibasty
use mod_ered, only: ered, rmu
use mod_hivector, only: dset
use mod_hitypes, only: bqs_type
implicit double precision (a-h,o-z)
real(8), dimension(nmax,nmax), intent(in) :: tsq
!      real(8), dimension(:,:), intent(in), target :: tototsq
real(8), dimension(nmax,nmax), intent(out) :: scmat
type(bqs_type), intent(in) :: row_bqs
type(bqs_type), intent(in) :: col_bqs
integer, dimension(nlevop), intent(in) :: inlev
integer, dimension(nlevop), intent(in) :: jlev
real(8), dimension(nlevop), intent(in) :: elev
integer, dimension(nmax) :: isc1  ! scratch array
integer, dimension(nmax) :: isc2  ! scratch array
real(8), dimension(nmax) :: sc2   ! scratch array

logical csflag, flaghf, flagsu, twomol
!
ASSERT(row_bqs%length == nopen)
ncol = col_bqs%length

xjtot = jtot
if (flaghf .and. .not. csflag) xjtot = jtot + 0.5d0
!  special fix for ibasty=23 (3P + 2S atom-atom)
if (ibasty.eq.23) xjtot = jtot + 0.5d0
if (flagsu) then
  pj = 1.d0
else
  pj = (pi * (2.d0*xjtot+1.d0)*ang2c)/(2.d0*rmu)
end if
fak=pj
if (csflag.and.(flaghf.or.nu.ne.0)) fak = 2.d0*pj
!  zero out scmat
do 10 icol = 1,nlevop
10 call dset(nlevop, 0.d0, scmat(1, icol), 1)

! 10    call fzero (scmat(1, icol), nlevop)
!  set pointer array for columns (final states)
do 40 i = 1, ncol
  do 20 icol = 1, nlevop
    if (is_j12(ibasty) .and. col_bqs%inq(i) .ne. inlev(icol)) &
          go to 20
    if (.not. twomol .and. col_bqs%inq(i).ne.inlev(icol)) go to 20
    if (col_bqs%jq(i) .ne. jlev(icol)) goto 20
    isc1(i) = icol
    goto 40
20   continue
  write (6, 30) i
  write (9, 30) i
30   format (' *** NO OPEN LEVEL FOUND FOR COLUMN CHANNEL',i3, &
          ' IN PARTCR; ABORT ***')
  call exit
40 continue
! set pointer array and degeneracy factors for rows (initial states)
do 140 j = 1, nopen
  do 120 irow = 1, nlevop
     if (is_j12(ibasty) .and. row_bqs%inq(j).ne.inlev(irow)) &
          go to 120
     if (.not. twomol.and.row_bqs%inq(j).ne.inlev(irow)) go to 120
     if (row_bqs%jq(j) .ne. jlev(irow)) goto 120
     jj = jlev(irow)
     if (.not. twomol) then
       xjrow1 = jj
       if (ibasty.eq.12 .or. ibasty.eq.13 .or. &
            ibasty.eq.15) then
! here for 2P/3P atom + homonuclear molecule in which case
! degeneracy factor is (2jmol+1)*(2*ja+1)
          xj2 = inlev(irow)
          if (flaghf) xj2 = xj2 + 0.5d0
          denrow = (2.d0*xjrow1+1.d0)*(2.d0*xj2+1.d0)
! here for 3P + 2S arom-atom collision
       elseif (ibasty .eq. 23) then
         denrow = (2.d0 * xjrow1 + 1.d0) * 2.d0
       else
         if (flaghf) xjrow1 = xjrow1 + 0.5d0
           denrow = (2.d0 * xjrow1 + 1.d0)
       endif
     else
       jrow1 = jj / 10
       jrow2 = mod (jj, 10)
       xjrow1 = jrow1
       xjrow2 = jrow2
       if (flaghf) then
         xjrow1 = xjrow1 + 0.5d0
! The following statement is removed on May 30, 2013; flaghf now only
! applies to the first molecule (j1).
!$$$               xjrow2 = xjrow2 + 0.5d0
       end if
       denrow = (2.d0 * xjrow1 + 1.d0) * (2.d0 * xjrow2 + 1.d0)
    end if
    if (.not. flagsu) denrow = denrow * (ered - elev(irow) )
    sc2(j) = fak/denrow
    isc2(j) = irow
    goto 140
120   continue
  write (6, 130) j
  write (9, 130) j
130   format (' *** NO OPEN LEVEL FOUND FOR ROW CHANNEL',i3, &
          ' IN PARTCR; ABORT ***')
  call exit
140 continue
!  the array isc1 contains, for each open channel, the index of the
!  corresponding open level
!  the array sc2 contains, for each open level, the normalization factor
!  to convert from t**2 to partial cross sections
do 190  i = 1, ncol
  icol = isc1(i)
  do 180  j = 1, nopen
    irow = isc2(j)
    scmat(irow,icol) = scmat(irow,icol)+tsq(j,i)*sc2(j)
180   continue
190 continue
return
end
end module mod_hibrid5
