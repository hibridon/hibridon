module mod_hypxsc
use omp_lib
implicit none
!*******************************************************************************
!     Module to compute hyperfine-resolved integral cross sections
!
!     reference:  alexander and dagdigian, jcp 83, 2191 (1985)
!     see also corey and mccourt, jpca 87, 2723 (1983) for
!     expression for spin-resolved T-matrix elements
!
!     this subroutine requires close-coupled s-matrix for
!     both parities
!
!     cross sections written to both terminal output and
!     {jobname}n.hfx file
!
!     Three cases can be treated:
!        - atom-molecule collisions with one nuclear spin
!        - atom-molecule collisions with two nuclear spin
!        - molecule-molecule collisions with one nuclear spin
!
!     original author:  p.j. dagdigian
!     this subroutine is a complete rewrite of an earlier
!     subroutine written by j. klos and f. lique
!
!     Refactored as a module and 'parallelized' by b. desrousseaux (jan. 2023)
!
!*******************************************************************************


! DERIVED DATA TYPES USED IN THIS SUBROUTINE
  type hflvl_type
    integer :: n
    real(8), allocatable :: e(:)
    integer, allocatable :: j(:)
    integer, allocatable :: in(:)
    real(8), allocatable :: if(:)
    real(8), allocatable :: f(:)
  end type hflvl_type

  type spin_type
    integer :: i(2) ! nuclear spin I times 2
    real(8) :: nuc(2)
    real(8) :: f = 0d0
    real(8) :: h = 0d0
    logical :: two = .false.
  end type spin_type

contains 

!************************************************************************************************
! Main subroutine that drives hyperine cross sections calculations and printing
!************************************************************************************************
subroutine hypxsc(flname, a)
! Modules used -----------------------------------------------------------------
  use mod_codim, only: mmax           ! Maximum number of channels
  use mod_hismat, only: sread, rdhead, sinqr  ! Read S-matrix from file
  use constants, only: econv                  ! To convert energy units
  use mod_selb, only: ibasty                  ! Basis type
  use mod_hiutil, only: gennam, mtime         ! To generate filenames and print time
  use mod_hitypes, only: bqs_type             ! To store basis set data
  use mod_parpot, only: potnam=>pot_name, label=>pot_label ! To store PES infos
  implicit none

  ! Scratch arrays -------------------------------------------------------------
  integer, allocatable :: jlev(:), inlev(:), jout(:)
  real(8), allocatable :: elev(:)
  real(8), allocatable :: sreal(:), simag(:)
  type(bqs_type)       :: row_bqs
  type(bqs_type)       :: packed_bqs
  ! Arguments ------------------------------------------------------------------
  character*(*), intent(in)         :: flname
  real(8), dimension(4), intent(in) :: a(4)

  ! To store data from argument a(i) -------------------------------------------
  integer :: iener, nucspin, j1min, j2max

  ! For filenames and IO -------------------------------------------------------
  character*40 :: smtfil, hfxfil ! Filenames for S-Matrix and output HFX file
  integer      :: hfxfil_unit = 11
  integer      :: lend
  logical      :: exstfl

  ! Read from S-matrix file ----------------------------------------------------
  integer :: mjtot, mchmx, lngth, ierr, jtot, jlpar, nu, nopen, jlp
  character*20 :: cdate
  integer :: jfrst, jfinl, jtotd, numin, numax, nud, nlevel, nlvop, nnout
  real(8) :: ered, rmu, ee
  logical :: nucrs, csflg, flaghf, flgsu, twmol

  ! To store S-matrix elements and channels info--------------------------------
  real(8), allocatable :: sr(:,:,:), si(:,:,:)
  type(bqs_type), allocatable :: bqs(:,:)

  ! MISC -----------------------------------------------------------------------
  integer :: len2, mchmx2

  ! Data related to spins -----------------------------------------------------
  type(spin_type) :: spins

  ! Data related to hyperfine levels -------------------------------------------
  type(hflvl_type) :: hf

  ! To store T-Matrices and Cross Sections -------------------------------------
  real(8), allocatable :: sigma(:,:)

  ! Define variables from hypxsc command arguments
  iener   = int(a(1))
  nucspin = int(a(2))
  j1min   = int(a(3))
  j2max   = int(a(4))

  ! Generate filename of smt file and check if it is present
  call gennam(smtfil, flname, iener, 'smt', lend)
  inquire(file = smtfil(1:lend), exist = exstfl)
  if (.not. exstfl) then ; write(6,"(3a)") '*** FILE ', smtfil(1:lend), ' NOT FOUND ***' ; return ; endif

  ! Open S-matrix file and read its header
  call openf(1, smtfil, 'tu', 0)
  call sinqr(1, mjtot, mchmx)

  allocate(elev(mmax))
  allocate(jlev(mmax))
  allocate(inlev(mmax))
  allocate(jout(mmax))

  call rdhead(1,cdate,ered,rmu,csflg,flaghf,flgsu, &
     twmol,nucrs,jfrst,jfinl,jtotd,numin,numax,nud, &
     nlevel,nlvop,nnout,jlev,inlev,elev,jout)
  deallocate(jout)

  ! Check if hypxsc supports the basis type and parameters provided in the input
  if (.not. supported(ibasty, flgsu, csflg, jtotd, nnout, nucspin, twmol)) return

  ! Compute spins
  call compute_spins(nucspin, flaghf, spins)

  ! Count and fill hyperfine levels for which XS are to be calculated
  call fill_hf(nlevel, jlev, elev, inlev, j1min, j2max, ered, twmol, spins, hf)
  deallocate(elev, inlev)

  ! Generate output file name and open it
  call gennam(hfxfil,flname,iener,'hfx',lend)
  call openf(hfxfil_unit, hfxfil(1:lend),'sf',0)

  ee = ered*econv ! Convert collision energy

  ! Write info to stdout and output file
  call print_infos(6, label, potnam, smtfil, cdate, jfinl, ee, nucspin)
  call print_infos(hfxfil_unit, label, potnam, smtfil, cdate, jfinl, ee, nucspin)
 
  mchmx2 = mchmx * (mchmx + 1) / 2

  allocate(sr(0:jfinl, 2, mchmx2)) ; sr = 0d0
  allocate(si(0:jfinl, 2, mchmx2)) ; si = 0d0
  allocate(bqs(0:jfinl, 2))                       
  

  allocate(sreal(mchmx2))
  allocate(simag(mchmx2))

  ! Read the S-Matrix elements
  do
    nopen = -1
    call sread (0, sreal, simag, jtot, jlpar, nu, row_bqs, packed_bqs, 1, mchmx2, nopen, ierr)
    if(ierr < -1) then ; write(6,*) '*** READ ERROR IN HYPXSC. ABORT ***' ; return ; endif
    jlp = 1 - (jlpar-1)/2
    len2 = packed_bqs%length*(packed_bqs%length + 1)/2
    
    call bqs(jtot,jlp)%init(packed_bqs%length)
    bqs(jtot,jlp)         = packed_bqs

    sr(jtot,jlp,1:len2)   = sreal(1:len2)
    si(jtot,jlp,1:len2)   = simag(1:len2)

    if(jtot==jfinl .and. jlpar==-1) exit
  enddo

  deallocate(sreal, simag)
  call packed_bqs%deinit()
  call row_bqs%deinit()



  ! Determine the type of collision and call corresponding subroutine to compute squared T-matrix elements
  if (.not. twmol) then 
    if(.not. spins%two) then ! Molecule-Atom with 1 nuclear spin
      write(6,"(a)") "HYPXSC | Molecule-Atom with 1 nuclear spin"
      call molecule_atom_1spin(jfrst, jfinl, nlevel, jlev, bqs, sr, si, spins, hf, sigma)
    else ! Molecule-Atom with 2 nuclear spins
      write(6,"(a)") "HYPXSC | Molecule-Atom with 2 nuclear spins"
      call molecule_atom_2spin(jfrst, jfinl, nlevel, jlev, bqs, sr, si, spins, hf, sigma)
    endif
  else ! Molecule-Molecule with 1 nuclear spin
      write(6,"(a)") "HYPXSC | Molecule-Molecule with 1 nuclear spins"
      call molecule_molecule_1spin(jfrst, jfinl, nlevel, jlev, bqs, sr, si, spins, hf, sigma)
  endif


  deallocate(sr, si, jlev, bqs)

  ! Compute hyperfine XS
  call compute_xs(twmol, rmu, ered, spins, hf, sigma)
  call print_xs(twmol, 6, ered, spins, hf, sigma)

  ! Print hyperfine XS
  call print_xs(twmol, hfxfil_unit, ered, spins, hf, sigma)



  deallocate(sigma)

  close(1)           ! Close S-Matrix file
  close(hfxfil_unit) ! Close output file
end subroutine hypxsc



!************************************************************************************************
! This function checks if hypxsc is implemented for user's case and returns .true. or .false.
!************************************************************************************************
function supported(ibasty, flgsu, csflg, jtotd, nnout, nucspin, twmol)
  implicit none 
  integer, intent(in) :: ibasty, jtotd, nnout, nucspin
  logical, intent(in) :: flgsu, csflg, twmol
  logical :: supported 

  supported = .true.
  if(any((/10,12,13,15,22,23/)==ibasty)) then
    write(6,"(A,I0,A)") '*** HYPERFINE CROSS SECTIONS FOR BASIS TYPE =', ibasty,' NOT IMPLEMENTED ***'
    supported = .false.
  endif
  if(flgsu) then
    write(6,"(A)") '*** HYPERFINE CROSS SECTIONS FOR SURFACE COLLISIONS NOT IMPLEMENTED ***'
    supported = .false.
  endif
  if(csflg) then
    write(6,"(A)") '*** CS HYPERFINE CROSS SECTIONS NOT IMPLEMENTED ***'
    supported = .false.
  endif 
  if (iabs(jtotd).ne.1) then
    write(6,"(A)") '*** DELTA-JTOT MUST BE EQUAL TO ONE ***'
    supported = .false.
  endif
  if (nnout.lt.0) then
    write(6,"(A)") '*** NNOUT < 0, ABORT ***'
    supported = .false.
  endif
  if (twmol .and. nucspin > 100) then
    write(6,"(a)") '*** HYPERFINE CROSS SECTIONS FOR MOLECULE-MOLECULE COLLISIONS NOT IMPLEMENTED FOR TWO NUCLEAR SPINS'
    supported = .false.
  endif
  return
end function supported



!************************************************************************************************
! This subroutine prints calculations info read from the S-matrix
!************************************************************************************************
subroutine print_infos(iunit, label, potnam, smtfil, cdate, jfinl, ered, nucspin)
  implicit none
  ! Arguments
  integer, intent(in) :: iunit, jfinl, nucspin
  character*20, intent(in) :: cdate
  character*48, intent(in) :: label, potnam
  character*40, intent(in) :: smtfil 
  ! Local variables
  real(8), intent(in) :: ered
  write(iunit,*)
  write(iunit,"(a)") '% HYPERFINE-RESOLVED CROSS SECTIONS'
  write(iunit,"(2a)") '%      LABEL:     ', label
  write(iunit,"(2a)") '%      POT NAME:  ', potnam
  write(iunit,"(2a)") '%      S-MATRICES READ FROM FILE: ', smtfil
  write(iunit,"(2a)") '%      WRITTEN:   ', cdate
  write(iunit,"(a,I10)") '%      JTOT2:     ', jfinl
  write(iunit,"(a,f10.3,a)") '%      ETOT:      ', ered, ' (cm-1)'
  if(nucspin<100) write(iunit,"(a,f6.1)") '%      NUCLEAR SPIN:      ', nucspin/2d0
  if(nucspin>100) write(iunit,"(a,2f6.1)") '%      NUCLEAR SPINS:      ', (nucspin/100)/2.d0, mod(nucspin,100)/2.d0
  return
end subroutine print_infos


!************************************************************************************************
! This subroutine compute the different spin values fspin, finuc, f2nuc, fhspin
!************************************************************************************************
subroutine compute_spins(nucspin, flaghf, spins)
  implicit none
  ! Arguments
  integer, intent(in) :: nucspin
  logical, intent(in) :: flaghf
  type(spin_type), intent(out) :: spins
  

  spins%f = 0d0
  if(flaghf) spins%f = 0.5d0 

  spins%h = 0d0
  if (flaghf .and. nucspin.eq.2*(nucspin/2)) spins%h = 0.5d0
  if (.not.flaghf .and. nucspin.ne.2*(nucspin/2)) spins%h = 0.5d0

  if(nucspin<100) then
    spins%two = .false.
    spins%nuc(1) = nucspin/2d0
    spins%i(1) = int(2*spins%nuc(1))
  else
    spins%two = .true.
    spins%nuc(1) = (nucspin/100)/2.d0
    spins%i(1) = int(2*spins%nuc(1))
    spins%nuc(2) = mod(nucspin,100)/2.d0
    spins%i(2) = int(2*spins%nuc(2))
  endif

end subroutine compute_spins



!************************************************************************************************
! This subroutine compute the number of hyperfine levels and fills the hyperfine arrays
!************************************************************************************************
subroutine fill_hf(nlevel, jlev, elev, inlev, j1min, j2max, ered, twmol, spins, hf)
  implicit none
  ! Arguments
  integer, intent(in) :: nlevel
  integer, intent(in) :: jlev(*)
  real(8), intent(in) :: elev(*)
  integer, intent(in) :: inlev(*)
  integer, intent(in) :: j1min, j2max
  real(8), intent(in) :: ered
  logical, intent(in) :: twmol
  type(spin_type), intent(in) :: spins
  type(hflvl_type), intent(out) :: hf
  ! Local variables
  type(hflvl_type) :: hf1, hf2
  integer :: i, ii, ij1, n, n2, nhyp
  real(8) :: ff, fj, ffmin, ffmax, f1, fnmin, fnmax, f2
  logical :: fill

  ! First pass will just count and allocate 
  ! Second pass will fill arrays
  fill = .false.
  do
    n = 0
    do i = 1, nlevel 
      ij1 = jlev(i) ; if(twmol) ij1 = jlev(i)/10
      if(ij1 >= j1min .and. ij1 <= j2max) then
        if(elev(i) <= ered) then
          fj = ij1 + spins%f
          ffmin = abs(fj - spins%nuc(1))
          ffmax = fj + spins%nuc(1)
          nhyp = int(ffmax - ffmin + 1)
          do ii = 1, nhyp
            n = n + 1
            if(fill) then
              hf1%j(n) = jlev(i)
              hf1%in(n) = inlev(i)
              ff = ffmin + (ii - 1)
              hf1%if(n) = int(ff)
              hf1%e(n) = elev(i)
            endif
          enddo
        endif
      endif
    enddo
    if(.not. fill) then
      hf1%n = n 
      allocate(hf1%e(hf1%n), hf1%j(hf1%n), hf1%in(hf1%n), hf1%if(hf1%n))
      fill = .true.
    else
      exit
    endif
  enddo

  hf = hf1


  ! THIS IS FOR TWO SPINS
  ! First pass will just count and allocate 
  ! Second pass will fill arrays
  if(spins%two) then
    fill = .false.
    do
      n2 = 0
      do i = 1, n
        f1 = hf1%if(i) + spins%h
        fnmin = abs(f1 - spins%nuc(2))
        fnmax = f1 + spins%nuc(2)
        nhyp = int(fnmax - fnmin + 1)
        do ii = 1,  nhyp
          n2 = n2 + 1
          if(fill) then
            hf2%j(n2) = hf1%j(i)
            hf2%in(n2) = hf1%in(i)
            hf2%f(n2) = f1
            f2 = fnmin + (ii - 1)
            hf2%if(n2) = f2
            hf2%e(n2) = hf1%e(i)
          endif
        enddo
      enddo
      if(.not. fill) then
        hf2%n = n2
        allocate(hf2%e(hf2%n), hf2%j(hf2%n), hf2%in(hf2%n), hf2%if(hf2%n), hf2%f(hf2%n))
        fill = .true.
      else
        exit
      endif
    enddo
    hf = hf2
  endif

  return
end subroutine fill_hf


!************************************************************************************************
! This subroutine computes XS from T matrix elements 
!************************************************************************************************
subroutine compute_xs(twmol, rmu, ered, spins, hf, sigma)
  use constants, only: ang2 => ang2c
  implicit none
  ! Arguments
  logical, intent(in) :: twmol
  real(8), intent(in) :: rmu, ered
  type(spin_type), intent(in) :: spins
  type(hflvl_type), intent(in) :: hf
  real(8), intent(inout) :: sigma(hf%n, hf%n)
  ! Local variables
  real(8) :: fak, dencol, denrow, ffi, fff
  integer :: ij2, ij2p, i, ii

  fak = acos(-1.d0) * ang2 / (2.0d0 * rmu)

if (.not. spins%two) then
  !$OMP PARALLEL DO PRIVATE(i, ii, ij2, ij2p, ffi, fff, denrow, dencol)
  do i = 1, hf%n
    ij2 = 0.d0
    if (twmol) ij2 = mod(hf%j(i),10)
    ffi = hf%if(i) + spins%h
    denrow = (2.d0 * ffi + 1.d0) * (2.d0 * ij2 + 1.d0) * (ered - hf%e(i))
    do ii = i, hf%n
      ij2p = 0.d0
      if (twmol) ij2p = mod(hf%j(ii),10)
      fff = hf%if(ii) + spins%h
      dencol = (2.d0 * fff + 1.d0) * (2.d0 * ij2p + 1.d0) * (ered - hf%e(ii))
      sigma(i,ii) = sigma(i,ii) * fak / denrow
      if (i.ne.ii) sigma(ii,i) = sigma(ii,i) * fak / dencol
    end do
  end do  
  !$OMP END PARALLEL DO
else
  !$OMP PARALLEL DO PRIVATE(i, ii, ij2, ij2p, ffi, fff, denrow, dencol)
  do i = 1, hf%n
    ij2 = 0.d0
    if (twmol) ij2 = mod(hf%j(i),10)
    ffi = hf%if(i)
    denrow = (2.d0 * ffi + 1.d0) * (2.d0 * ij2 + 1.d0) * (ered - hf%e(i))
    do ii = i, hf%n
      ij2p = 0.d0
      if (twmol) ij2p = mod(hf%j(ii),10)
      fff = hf%if(ii)
      dencol = (2.d0 * fff + 1.d0) * (2.d0 * ij2p + 1.d0) * (ered - hf%e(ii))
      sigma(i,ii) = sigma(i,ii) * fak / denrow
      if (i.ne.ii) sigma(ii,i) = sigma(ii,i) * fak / dencol
    end do
  end do
  !$OMP END PARALLEL DO
end if


end subroutine compute_xs


!************************************************************************************************
! This subroutine prints XS 
!************************************************************************************************
subroutine print_xs(twmol, hfxfil_unit, ered, spins, hf, sigma)
  use constants, only: econv
  implicit none
  ! Arguments
  logical, intent(in) :: twmol
  integer, intent(in) :: hfxfil_unit
  real(8), intent(in) :: ered
  type(spin_type), intent(in) :: spins
  type(hflvl_type), intent(in) :: hf
  real(8), intent(in) :: sigma(hf%n,hf%n)
  ! Local variables
  integer :: i, ii, ij2, ij2p
  real(8) :: xj, xf, xjp, xfp, ee, xf2, xf2p

  ! Write header
  if(.not.spins%two) then
  ! ONE SPIN
    if (.not.twmol) then
    ! ONE MOLECULE
      write(hfxfil_unit,"(a)") '%     E(CM-1)     JI     INI   FI      JF     INF   FF      CROSS SECTION (ANG^2)'
    else
    ! TWO MOLECULES
      write(hfxfil_unit,"(a)") '%    E(CM-1)     JI     INI   FI    J2      JF     INF   FF    J2P      CROSS SECTION (ANG^2)'
    end if
  else 
  ! TWO SPINS
    if (.not.twmol) then
    ! ONE MOLECULE
      write(hfxfil_unit,"(a)") '%     E(CM-1)  JI   INI  F1I   F2I     JF  INF  F1F   F2F  CROSS SECTION (ANG^2)'
    end if
  end if

ee = ered*econv

  if (.not.spins%two) then 
  ! ONE SPIN
    do i = 1, hf%n
      do ii = 1, hf%n
        if (.not. twmol) then
        ! ONE MOLECULE
          xj = hf%j(i) + spins%f
          xf = hf%if(i) + spins%h
          xjp =  hf%j(ii) + spins%f
          xfp = hf%if(ii) + spins%h
          if (sigma(i,ii)>0d0) then
            write(hfxfil_unit,"(f12.3,f8.1,i6,f6.1,3x,f6.1,i6,f6.1,5x,1pe15.4)")&
                  ee,xj,hf%in(i),xf,xjp,hf%in(ii),xfp,sigma(i,ii)
          end if
        else
        ! TWO MOLECULES
          xj = (hf%j(i)/10) + spins%f
          ij2 = mod(hf%j(i),10)
          xf  = hf%if(i) + spins%h
          xjp = (hf%j(ii)/10) + spins%f
          ij2p = mod(hf%j(ii),10)
          xfp  = hf%if(ii) + spins%h
          if (sigma(i,ii)>0d0) then
            write(hfxfil_unit,"(f12.3,f8.1,i6,f6.1,i6,3x,f6.1,i6,f6.1,i6,5x,1pe15.4)")&
                  ee,xj,hf%in(i),xf,ij2,xjp,hf%in(ii),xfp,ij2p,sigma(i,ii)
          endif
        endif
      enddo
    enddo

  else
  ! TWO SPINS
    do i = 1, hf%n
      do ii = 1, hf%n
        xj = hf%j(i) + spins%f
        xf = hf%f(i)
        xf2 = hf%if(i)
        xjp = hf%j(ii) + spins%f
        xfp = hf%f(ii)
        xf2p = hf%if(ii)
        if (sigma(i,ii)>0.d0) then
          write(hfxfil_unit,"(f12.3,f6.1,i4,f6.1,f6.1,f7.1,i4,f6.1,f6.1,2x,1pe15.4)")&
          ee,xj,hf%in(i),xf,xf2,xjp,hf%in(ii),xfp,xf2p,sigma(i,ii)
        endif
      end do
    end do
  end if


end subroutine print_xs




!************************************************************************************************
! This subroutine computes T matrix elements for molecule-atom collisions with one spin
!************************************************************************************************
subroutine molecule_atom_1spin(jfrst, jfinl, nlevel, jlev, bqs, sr, si, spins, hf, sigma)
  use mod_hitypes, only: bqs_type
  use mod_hiutil, only: xf6j
  implicit none
  ! Arguments
  integer, intent(in) :: jfrst, jfinl, nlevel
  integer, intent(in) :: jlev(*)
  type(bqs_type), intent(in) :: bqs(0:jfinl,2)
  real(8), intent(in) :: sr(0:jfinl,2,*) , si(0:jfinl,2,*)
  type(spin_type), intent(in) :: spins
  type(hflvl_type), intent(in) :: hf
  real(8), allocatable, intent(out) :: sigma(:,:)
  ! Local variables
  integer :: jmx, idim
  ! Used within the loops
  real(8), allocatable :: tmatr(:,:), tmati(:,:)
  integer :: iftot, iftmn, iftmx, jlp, jlpar, i, ii, is
  real(8) :: xftot, xj, xjp, xf, xfp, fffp, xjttmn, xjtot, fjtot, xl, xlp, t2sum, phase
  integer :: jttmin, jttmax, jtot, jlparf, jlpt, irow, icol, iflag, iph, ll, lp
  complex(8) :: t


  ! Boundaries for loop over iftot  
  iftmn = max(0,int(jfrst + spins%f - spins%nuc(1) - spins%h))
  iftmx = int(jfinl + spins%f + spins%nuc(1) - spins%h)

  ! Allocate sigma array
  allocate(sigma(hf%n, hf%n))
  sigma = 0d0
  ! Allocate work arrays
  jmx = maxval(jlev(1:nlevel))
  idim = (iftmx + 2*jmx + 3) * (iftmx + jmx + 2)

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(idim,jfrst,jfinl,iftmn,iftmx,spins,bqs,hf,sr,si) REDUCTION(+:sigma)
  allocate(tmatr(idim, idim)) 
  allocate(tmati(idim, idim)) 
!$OMP DO
  do iftot = iftmn, iftmx
    xftot = iftot + spins%h
    do jlp = 1, 2
      jlpar = 1 - (jlp -1)*2
      write(6,"(A,F5.1,A,I4)") ' Computing partial wave J_tot =', xftot, ', jlpar=', jlpar
      do i = 1, hf%n
        do ii = i, hf%n
          xj = hf%j(i) + spins%f
          xjp = hf%j(ii) + spins%f
          xf = hf%if(i) + spins%h
          xfp = hf%if(ii) + spins%h
          fffp = sqrt((2.d0*xf+1.d0)*(2.d0*xfp+1.d0))
          xjttmn = max(spins%h, xftot - spins%nuc(1))
          jttmin = max(jfrst, int(xjttmn-spins%f))
          jttmax = min(jfinl, int(xftot + spins%nuc(1)))
  !     clear T-matrix array
          tmatr = 0.d0
          tmati = 0.d0
  !     sum over jtot consistent with vector addition
  !     ftot = jtot + nucspin
          do jtot=jttmin,jttmax
            xjtot = jtot + spins%f
            fjtot = 2.d0 * xjtot + 1.d0
  !     total parity must be the same for all T-matrix elements
  !     in sum over jtot (for the same parity, jlparf
  !     changes sign for each increase in jtot by 1)
            jlparf = jlpar*(-1)**(iftot - jtot)
            jlpt = 1 + (1 - jlparf)/2
            if (bqs(jtot,jlpt)%length .gt. 0) then
              do irow = 1, bqs(jtot,jlpt)%length
  !     flag to make sure initial level is the bra, final level the ket
                iflag = 0
                if (bqs(jtot,jlpt)%jq(irow) == hf%j(i) .and. bqs(jtot,jlpt)%inq(irow) == hf%in(i)) then
                  iflag = 1
                endif
                do icol = 1, irow
                  if (bqs(jtot,jlpt)%jq(irow) .eq. hf%j(i)  .and. &
                      bqs(jtot,jlpt)%inq(irow).eq. hf%in(i) .and. &
                      bqs(jtot,jlpt)%jq(icol) .eq. hf%j(ii) .and. &
                      bqs(jtot,jlpt)%inq(icol).eq. hf%in(ii) .or. &
                      bqs(jtot,jlpt)%jq(icol) .eq. hf%j(i)  .and. &
                      bqs(jtot,jlpt)%inq(icol).eq. hf%in(i) .and. &
                      bqs(jtot,jlpt)%jq(irow) .eq. hf%j(ii) .and. &
                      bqs(jtot,jlpt)%inq(irow).eq. hf%in(ii)) then
                    is = (irow*(irow - 1))/2 + icol
                    if (iflag.ne.1) then
                      xl = bqs(jtot,jlpt)%lq(icol)
                      xlp = bqs(jtot,jlpt)%lq(irow)
                    else
                      xl = bqs(jtot,jlpt)%lq(irow)
                      xlp = bqs(jtot,jlpt)%lq(icol)
                    end if
  !     convert S-matrix element to T-matrix element
                    t = -cmplx(sr(jtot,jlpt,is), si(jtot,jlpt,is), 8)
  !     next statement for diagonal T-matrix element
                    if (irow.eq.icol) then
                      t = t + cmplx(1.d0, 0.d0, 8)
                    end if
                    iph = int(xfp - xlp - xf + xl)
                    phase = 1.d0
                    if (iph.ne.2*(iph/2)) phase = -1.d0
                    t = t * phase * fffp * fjtot &
                        * xf6j(spins%nuc(1),xj,xf,xl,xftot,xjtot) &
                      * xf6j(spins%nuc(1),xjp,xfp,xlp,xftot,xjtot)
                    ll = int(xl)
                    lp = int(xlp)
                    tmatr(ll+1,lp+1) = tmatr(ll+1,lp+1) + real(t)
                    tmati(ll+1,lp+1) = tmati(ll+1,lp+1) + aimag(t)
  !     for initial level = final level, but l.ne.lp, need to include
  !     both T(l,lp) and T(lp,l)
                    if (hf%j(i) == hf%j(ii) .and. hf%in(i) == hf%in(ii) .and. abs(xl-xlp) > 1d-60) then
  !     note that T(l,lp) = T(lp,l)* (Hermitean matrix)
                      tmatr(lp+1,ll+1) = tmatr(lp+1,ll+1) + real(t)
                      tmati(lp+1,ll+1) = tmati(lp+1,ll+1) - aimag(t)
                    end if
  !     if statement below is end of tests for triangle relations
                  end if
                end do
              end do
  !     end of if statement checking in length
            end if
          end do
          t2sum = sum(tmatr*tmatr + tmati*tmati)
          sigma(i,ii) = sigma(i,ii) + t2sum * (2.d0 * xftot + 1.d0)
          if (i.ne.ii) sigma(ii,i) = sigma(ii,i) + t2sum * (2.d0 * xftot + 1.d0)
        end do
      end do
    end do
  enddo
!$OMP END DO
  deallocate(tmatr)
  deallocate(tmati)
!$OMP END PARALLEL
end subroutine molecule_atom_1spin




!************************************************************************************************
! This subroutine computes T matrix elements for molecule-atom collisions with two spins
!************************************************************************************************
subroutine molecule_atom_2spin(jfrst, jfinl, nlevel, jlev, bqs, sr, si, spins, hf, sigma)
  use mod_hitypes, only: bqs_type
  use mod_hiutil, only: xf6j
  implicit none
  ! Arguments
  integer, intent(in) :: jfrst, jfinl, nlevel
  integer, intent(in) :: jlev(*)
  type(bqs_type), intent(in) :: bqs(0:jfinl,2)
  real(8), intent(in) :: sr(0:jfinl,2,*) , si(0:jfinl,2,*)
  type(spin_type), intent(in) :: spins
  type(hflvl_type), intent(in) :: hf
  real(8), allocatable, intent(out) :: sigma(:,:)
  ! Local variables
  integer :: jmx, idim
  real(8) :: dspin
  ! Used within the loops
  real(8), allocatable :: tmatr(:,:), tmati(:,:)
  integer :: iftot, iftmn, iftmx, jlp, jlpar, i, ii, is, iqmax, iqmin, ir, irmax, irmin
  real(8) :: xftot, xj, xjp, xf, xfp, fffp, xjtot, fjtot, xl, xlp, t2sum, rprod, xf2, xf2p, r, rmin, phase
  integer :: jtot, jlparf, jlpt, irow, icol, iflag, iph, ll, lp
  complex(8) :: t


  ! Boundaries for loop over iftot  
  iftmn = max(0,int(jfrst + spins%f - spins%nuc(1) - spins%nuc(2)))
  iftmx = int(jfinl + spins%f + spins%nuc(1) + spins%nuc(2))

  ! Allocate sigma array
  allocate(sigma(hf%n, hf%n))
  sigma = 0d0

  ! Allocate work arrays
  jmx = maxval(jlev(1:nlevel))
  idim = (iftmx + 2*jmx + 3) * (iftmx + jmx + 2)

  dspin = spins%f + spins%nuc(1)  + spins%nuc(2)
  if (abs(dspin - 2d0*int(dspin/2)) < 1d-60) then
    dspin = 0d0
  else
    dspin = 0.5d0
  endif

  !  NOTE:  xftot is called K in Lara-Moreno et al.
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(idim,jfrst,jfinl,iftmn,iftmx,spins,bqs,hf,sr,si, dspin) REDUCTION(+:sigma)
  allocate(tmatr(idim, idim))
  allocate(tmati(idim, idim))
!$OMP DO
  do iftot = iftmn, iftmx-1
    xftot = float(iftot) + dspin
    do jlp = 1, 2
      jlpar = 1 - (jlp -1)*2
      write(6,"(A,F5.1,A,I4)") ' Computing partial wave K_tot =', xftot, ', jlpar=', jlpar
      do i=1, hf%n
        xj = hf%j(i) + spins%f
        xf = hf%f(i)
        xf2 = hf%if(i)
        do ii=i, hf%n
          xjp = hf%j(ii) + spins%f
          xfp = hf%f(ii)
          xf2p = hf%if(ii)
          fffp = sqrt((2.d0*xf+1.d0)*(2.d0*xfp+1.d0) &
            * (2.d0*xf2+1.d0)*(2.d0*xf2p+1.d0))
          iph = int(xf - xfp - xf2 + xf2p)
          phase = (-1.d0)**iph
  !  clear T-matrix array
          tmatr = 0.d0
          tmati = 0.d0
  !  sum over R = K - I2 consistent with vector addition
          rmin = abs(xftot - spins%nuc(2))
          irmin = int(rmin)
          irmax = int(xftot + spins%nuc(2))
          do ir = irmin, irmax
            r  = rmin + (ir - irmin)
            iqmin = max(jfrst, int(abs(r - spins%nuc(1))))
            iqmax = min(jfinl, int(r + spins%nuc(1)))
            do jtot = iqmin, iqmax
              xjtot = jtot + spins%f
              fjtot = 2.d0 * xjtot + 1.d0
  !  total parity must be the same for all T-matrix elements
  !  in sum over jtot (for the same parity, jlparf
  !  changes sign for each increase in jtot by 1)
              jlparf = jlpar*(-1)**(iftot - jtot)
              jlpt = 1 + (1 - jlparf)/2
              if (bqs(jtot,jlpt)%length .gt. 0) then
              do irow=1,bqs(jtot,jlpt)%length
  !  flag to make sure initial level is the bra, final level the ket
                iflag = 0
                if (bqs(jtot,jlpt)%jq(irow).eq.hf%j(i) .and. &
                    bqs(jtot,jlpt)%inq(irow).eq.hf%in(i)) then
                  iflag = 1
                end if
                do icol=1,irow
                  if (bqs(jtot,jlpt)%jq(irow).eq.hf%j(i) .and. &
                      bqs(jtot,jlpt)%inq(irow).eq.hf%in(i) .and. &
                      bqs(jtot,jlpt)%jq(icol).eq.hf%j(ii) .and. &
                      bqs(jtot,jlpt)%inq(icol).eq.hf%in(ii) .or. &
                      bqs(jtot,jlpt)%jq(icol).eq.hf%j(i) .and. &
                      bqs(jtot,jlpt)%inq(icol).eq.hf%in(i) .and. &
                      bqs(jtot,jlpt)%jq(irow).eq.hf%j(ii) .and. &
                      bqs(jtot,jlpt)%inq(irow).eq.hf%in(ii)) then
                    is = (irow*(irow - 1))/2 + icol
                    if (iflag.ne.1) then
                      xl = bqs(jtot,jlpt)%lq(icol)
                      xlp = bqs(jtot,jlpt)%lq(irow)
                    else
                      xl = bqs(jtot,jlpt)%lq(irow)
                      xlp = bqs(jtot,jlpt)%lq(icol)
                    end if
  !  convert S-matrix element to T-matrix element
                    t = -cmplx(sr(jtot,jlpt,is), si(jtot,jlpt,is), 8)
  !     next statement for diagonal T-matrix element
                    if (irow.eq.icol) then
                      t = t + cmplx(1.d0, 0.d0, 8)
                    end if
                    rprod = (2.d0*r + 1.d0) &
                      * xf6j(spins%nuc(2),xf,xf2,xl,xftot,r) &
                      * xf6j(spins%nuc(2),xfp,xf2p,xlp,xftot,r) &
                      * xf6j(spins%nuc(1),xj,xf,xl,r,xjtot) &
                      * xf6j(spins%nuc(1),xjp,xfp,xlp,r,xjtot)
                    t = t * phase * fffp * (2.d0*xjtot + 1.d0) * rprod
                    ll = int(xl)
                    lp = int(xlp)
                    tmatr(ll+1,lp+1) = tmatr(ll+1,lp+1) + real(t)
                    tmati(ll+1,lp+1) = tmati(ll+1,lp+1) + aimag(t)
  !     for initial level = final level, but l.ne.lp, need to include
  !     both T(l,lp) and T(lp,l)
                    if (hf%j(i) == hf%j(ii) .and. hf%in(i) == hf%in(ii) .and. abs(xl-xlp) > 1d-60) then
  !     note that T(l,lp) = T(lp,l)* (Hermitean matrix)
                      tmatr(lp+1,ll+1) = tmatr(lp+1,ll+1) + real(t)
                      tmati(lp+1,ll+1) = tmati(lp+1,ll+1) - aimag(t)
                    end if
  !     if statement below is end of tests for triangle relations
                  end if
                end do
              end do
  !     end of if statement checking in length
            end if
          end do
          end do
          t2sum = sum(tmatr*tmatr + tmati*tmati)
          sigma(i,ii) = sigma(i,ii) + t2sum * (2.d0 * xftot + 1.d0)
          if (i.ne.ii) then
            sigma(ii,i) = sigma(ii,i) + t2sum * (2.d0 * xftot + 1.d0)
          end if
        end do
      end do
    end do
  end do
!$OMP END DO
  deallocate(tmatr)
  deallocate(tmati)
!$OMP END PARALLEL

end subroutine molecule_atom_2spin



!************************************************************************************************
! This subroutine computes T matrix elements for molecule-molecule collisions with one spin
!************************************************************************************************
subroutine molecule_molecule_1spin(jfrst, jfinl, nlevel, jlev, bqs, sr, si, spins, hf, sigma)
  use mod_hitypes, only: bqs_type
  use mod_hiutil, only: xf6j
  implicit none
  ! Arguments
  integer, intent(in) :: jfrst, jfinl, nlevel
  integer, intent(in) :: jlev(*)
  type(bqs_type), intent(in) :: bqs(0:jfinl,2)
  real(8), intent(in) :: sr(0:jfinl,2,*) , si(0:jfinl,2,*)
  type(spin_type), intent(in) :: spins
  type(hflvl_type), intent(in) :: hf
  real(8), allocatable, intent(out) :: sigma(:,:)
  ! Local variables
  integer :: jmx, jmx2, idim
  ! Used within the loops
  real(8), allocatable :: tmatr(:,:), tmati(:,:)
  integer :: jttmax, jttmin, iftot, iftmn, iftmx, jlp, jlpar, i, ii, is
  real(8) :: xjttmn, xftot, xj, xjp, xf, xfp, fffp, xjtot, fjtot, xl, xlp, t2sum
  real(8) :: fj12, fj12p, fjr, fjrp, ph, xj12, xj12p, xj2, xj2p, xjr, xjrp
  integer :: jtot, jlparf, jlpt, irow, icol, iflag, iph, ll, lp, idimr, isp, jr, jrmax, jrmin, jrp, jrpmax, jrpmin
  complex(8) :: t, tf


  ! Boundaries for loop over iftot  
  iftmn = max(0,int(jfrst + spins%f - spins%nuc(1) - spins%h))
  iftmx = int(jfinl + spins%f + spins%nuc(1) - spins%h)

  ! Allocate sigma array
  allocate(sigma(hf%n, hf%n))
  sigma = 0d0

  ! Allocate work arrays
  jmx = maxval(jlev(1:nlevel)/10)
  jmx2 = maxval(mod(jlev(1:nlevel),10))
  idim = (iftmx + jmx + 2*jmx2 + 3) * (iftmx + jmx + jmx2 + 2)
  idimr = idim / (iftmx + jmx + jmx2 +1)

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(idim,idimr,jfrst,jfinl,iftmn,iftmx,spins,bqs,hf,sr,si) REDUCTION(+:sigma)
  allocate(tmatr(idim, idim))
  allocate(tmati(idim, idim))
!$OMP DO
  do iftot = iftmn, iftmx
    xftot = iftot + spins%h
    do jlp = 1, 2
      jlpar = 1 - (jlp -1) * 2
      write(6,"(A,F5.1,A,I4)") ' Computing partial wave J_tot =', xftot, ', jlpar=', jlpar
      do i = 1, hf%n
        do ii = i, hf%n
          xj = (hf%j(i)/10) + spins%f
          xjp = (hf%j(ii)/10) + spins%f
          xj2 = mod(hf%j(i),10)
          xj2p = mod(hf%j(ii),10)
          xf = hf%if(i) + spins%h
          xfp = hf%if(ii) + spins%h
          fffp = sqrt((2.d0*xf+1.d0)*(2.d0*xfp+1.d0))
          xjttmn = max(spins%h, xftot - spins%nuc(1))
          jttmin = int(xjttmn - spins%f)
          jttmin = max(jfrst, jttmin)
          jttmax = int(xftot + spins%nuc(1))
          jttmax = min(jfinl, jttmax)
  !     clear tmat arrays
          tmatr = 0.d0
          tmati = 0.d0
  !     sum over jtot consistent with vector addition
  !     ftot = jtot + nucspin
          do jtot=jttmin,jttmax
            xjtot = jtot + spins%f
            fjtot = 2.d0 * xjtot + 1.d0
  !     total parity must be the same for all T-matrix elements
  !     in sum over jtot (for the same parity, jlparf
  !     changes sign for each increase in jtot by 1)
            jlparf = jlpar*(-1)**(iftot - jtot)
            jlpt = 1 + (1 - jlparf)/2
            do irow = 1, bqs(jtot,jlpt)%length
  !     flag to make sure initial level is the bra, final level the ket
              iflag = 0
              if (bqs(jtot,jlpt)%jq(irow) == hf%j(i) .and. bqs(jtot,jlpt)%inq(irow) == hf%in(i)) then
                iflag = 1
              end if
              do icol=1,irow
                if (bqs(jtot,jlpt)%jq(irow).eq.hf%j(i) .and. &
                    bqs(jtot,jlpt)%inq(irow).eq.hf%in(i) .and. &
                    bqs(jtot,jlpt)%jq(icol).eq.hf%j(ii) .and. &
                    bqs(jtot,jlpt)%inq(icol).eq.hf%in(ii) .or. &
                    bqs(jtot,jlpt)%jq(icol).eq.hf%j(i) .and. &
                    bqs(jtot,jlpt)%inq(icol).eq.hf%in(i) .and. &
                    bqs(jtot,jlpt)%jq(irow).eq.hf%j(ii) .and. &
                    bqs(jtot,jlpt)%inq(irow).eq.hf%in(ii)) then
                  is = (irow*(irow - 1))/2 + icol
                  if (iflag.ne.1) then
                    xl = bqs(jtot,jlpt)%lq(icol)
                    xlp = bqs(jtot,jlpt)%lq(irow)
                    xj12 = bqs(jtot,jlpt)%j12(icol) + spins%f
                    xj12p = bqs(jtot,jlpt)%j12(irow) + spins%f
                  else
                    xl = bqs(jtot,jlpt)%lq(irow)
                    xlp = bqs(jtot,jlpt)%lq(icol)
                    xj12 = bqs(jtot,jlpt)%j12(irow) + spins%f
                    xj12p = bqs(jtot,jlpt)%j12(icol) + spins%f
                  end if
                  fj12 = sqrt(2.d0 * xj12 + 1.d0)
                  fj12p = sqrt(2.d0 * xj12p + 1.d0)
  !     convert S-matrix element to T-matrix element
                  t = -cmplx(sr(jtot,jlpt,is), si(jtot,jlpt,is), 8)
  !     next statement for diagonal T-matrix element
                  if (irow.eq.icol) then
                    t = t + cmplx(1.d0, 0.d0, 8)
                  end if
  !     sums over jR and jRp
                  jrmin = int(abs(xj2 - xl))
                  jrmax = int(xj2 + xl)
                  jrpmin = int(abs(xj2p - xlp))
                  jrpmax = int(xj2p + xlp)
                  do jr = jrmin, jrmax
                    xjr = jr
                    fjr = sqrt(2.d0 * xjr + 1.d0)
                    do jrp = jrpmin, jrpmax
                      xjrp = jrp
                      fjrp = sqrt(2.d0 * xjrp + 1.d0)
                      iph = int(xjr + xjrp + xl + xlp + xj2 + xj2p)
                      ph = 1.d0
                      if (iph.ne.2*(iph/2)) ph = -1.d0
                      tf = t * ph * fjtot * fffp &
                          * fjr * fjrp * fj12 * fj12p &
                          * xf6j(xj,xj2,xj12,xl,xjtot,xjr) &
                          * xf6j(xjp,xj2p,xj12p,xlp,xjtot,xjrp) &
                          * xf6j(xjr,xj,xjtot,spins%nuc(1),xftot,xf) &
                          * xf6j(xjrp,xjp,xjtot,spins%nuc(1),xftot,xfp)
                      ll = int(xl)
                      lp = int(xlp)
                      is =  (ll + 1) * (idimr - 1) + (jr + 1)
                      isp = (lp + 1) * (idimr - 1) + (jrp + 1)
                      tmatr(is,isp) = tmatr(is,isp) + real(tf)
                      tmati(is,isp) = tmati(is,isp) + aimag(tf)
  !     for initial level = final level, but l.ne.lp or j12.ne.j12p, need to include
  !     both T(l,lp) and T(lp,l)
                      if (hf%j(i) == hf%j(ii) .and. hf%in(i) == hf%in(ii) .and. irow /= icol) then
  !     note that T(l,lp) = T(lp,l)* (Hermitean matrix)
                        tmatr(isp,is) = tmatr(isp,is) + real(tf)
                        tmati(isp,is) = tmati(isp,is) - aimag(tf)
                      end if
                    end do
                  end do
                end if
              end do
            end do
          end do
          t2sum = sum(tmatr*tmatr + tmati*tmati)
  !  do not include contribution from last ftot
          if (iftot.ne.iftmx) then
            sigma(i,ii) = sigma(i,ii) + t2sum * (2.d0 * xftot + 1.d0)
            if (i.ne.ii) then
              sigma(ii,i) = sigma(ii,i) + t2sum * (2.d0 * xftot + 1.d0)
            end if
          end if
        end do
      end do
    end do
  end do
!$OMP END DO
  deallocate(tmatr)
  deallocate(tmati)
!$OMP END PARALLEL
end subroutine molecule_molecule_1spin



end module mod_hypxsc