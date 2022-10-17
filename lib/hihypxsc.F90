module mod_hypxsc
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
!        - atom-molecule collisions with one nuclear spine
!        - atom-molecule collisions with two nuclear spine
!        - molecule-molecule collisions with one nuclear spine
!
!     original author:  p.j. dagdigian
!     this subroutine is a complete rewrite of an earlier
!     subroutine written by j. klos and f. lique
!
!     Refactored as a module and parallelized by 
!              b. desrousseaux (oct. 2022)
!
!*******************************************************************************


! DERIVED DATA TYPES USED IN THIS SUBROUTINE
  type hflvl
    integer :: n
    real(8), allocatable :: e(:)
    integer, allocatable :: j(:)
    integer, allocatable :: in(:)
    integer, allocatable :: if(:)
  end type hflvl

  type spin_data
    integer :: i1 ! nuclear spin I times 2
    integer :: i2 ! nuclear spin I times 2
    real(8) :: nuc1
    real(8) :: nuc2
    real(8) :: f
    real(8) :: h
    logical :: two

  end type spin_data

contains 

!************************************************************************************************
! Main subroutine that drives hyperine cross sections calculations and printing
!************************************************************************************************
subroutine hypxsc(flname, a)
! Modules used -----------------------------------------------------------------
  use mod_coj12, only: j12q => j12
  use mod_codim, only: mmax
  use mod_cojq, only: jq ! jq(1)
  use mod_colq, only: lq ! lq(1)
  use mod_coinq, only: inq ! inq(1)
  use mod_coisc1, only: jlev => isc1 ! jlev(1)
  use mod_coisc3, only: inlev => isc3 ! inlev(1)
  use mod_coisc5, only: jout => isc5 ! jout(1)
  use mod_coisc10, only: ipack => isc10 ! ipack(1)
  use mod_coisc11, only: jpack => isc11 ! jpack(1)
  use mod_coisc12, only: lpack => isc12 ! lpack(1)
  use mod_cosc1, only: elev => sc1 ! elev(1)
  use mod_hismat, only: sread
  use constants, only: econv
  use mod_parpot, only: potnam=>pot_name, label=>pot_label
  use mod_selb, only: ibasty
  use mod_hiutil, only: gennam, mtime
  use mod_hismat, only: sread, rdhead, sinqr
  
  implicit none

  ! Arguments ------------------------------------------------------------------
  character*(*), intent(in) :: flname
  real(8), dimension(4), intent(in) :: a(4)
  
  ! To store data from argument a(i) -------------------------------------------
  integer :: iener, nucspin, j1min, j2max

  ! For filenames and IO -------------------------------------------------------
  character*40 :: smtfil, hfxfil ! Filenames for S-Matrix and output HFX file
  integer :: hfxfil_unit
  integer :: lend
  logical :: exstfl

  ! Read from S-matrix file ----------------------------------------------------
  integer :: mjtot, mchmx, lngth, ierr, jtot, jlpar, nu, nopen, jlp
  character*20 :: cdate
  integer :: jfrst, jfinl, jtotd, numin, numax, nud, nlevel, nlvop, nnout
  real(8) :: ered, rmu, ee
  logical :: nucrs, csflg, flaghf, flgsu, twmol
  real(8), allocatable :: sreal(:), simag(:)
  real(8), allocatable :: sr(:,:,:), si(:,:,:)
  integer, allocatable :: j(:,:,:), in(:,:,:), l(:,:,:), j12(:,:,:), length(:,:)
  logical, allocatable :: exsmtn(:), exsmtp(:)

  ! MISC -----------------------------------------------------------------------
  integer :: len2, mchmx2

  ! Data related to spins -----------------------------------------------------
  type(spin_data) :: spins

  ! Data related to hyperfine levels -------------------------------------------
  type(hflvl) :: hf1, hf2

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
  call rdhead(1,cdate,ered,rmu,csflg,flaghf,flgsu, &
     twmol,nucrs,jfrst,jfinl,jtotd,numin,numax,nud, &
     nlevel,nlvop,nnout,jlev,inlev,elev,jout)

  ! Check if hypxsc supports the basis type and parameters provided in the input
  if (.not. supported(ibasty, flgsu, csflg, jtotd, nnout, nucspin, twmol)) return

  ! Generate output file name and open it
  call gennam(hfxfil,flname,iener,'hfx',lend)
  call openf(hfxfil_unit, hfxfil(1:lend),'sf',0)

  ee = ered*econv ! Convert collision energy
  ! Write info to stdout and output file
  call print_infos(6, label, potnam, smtfil, cdate, jfinl, ee, nucspin)
  call print_infos(hfxfil_unit, label, potnam, smtfil, cdate, jfinl, ee, nucspin)

  ! Allocate and init arrays
  mchmx2 = mchmx * (mchmx + 1) / 2
  allocate(sreal(mchmx2), simag(mchmx2))
  sreal = 0d0 ; simag = 0d0
  allocate(sr(0:jfinl, 2, mchmx2), si(0:jfinl, 2, mchmx2))
  sr = 0d0 ; si = 0d0
  allocate(j(0:jfinl, 2, mchmx), in(0:jfinl, 2, mchmx), l(0:jfinl, 2, mchmx), j12(0:jfinl, 2, mchmx))
  j = 0 ; in = 0 ; l = 0 ; j12 = 0
  allocate(exsmtp(0:jfinl), exsmtn(0:jfinl))
  exsmtp = .false. ; exsmtn = .false.
  allocate(length(0:jfinl, 2))

  ! Read the S-Matrix
  do
    nopen = -1
    call sread (0, sreal, simag, jtot, jlpar, nu, jq, lq, inq, ipack, jpack, lpack, 1, mmax, nopen, lngth, ierr)
    if(ierr < -1) then ; write(6,*) '*** READ ERROR IN HYPXSC. ABORT ***' ; return ; endif
    jlp = 1 - (jlpar-1)/2
    length(jtot, jlp) = lngth 
    len2 = lngth*(lngth + 1)/2
    if(jlpar==1) then ; exsmtp(jtot) = .true. ; else ; exsmtn(jtot) = .true. ; endif
    j(jtot,jlp,1:lngth) = jpack(1:lngth)
    in(jtot,jlp,1:lngth) = ipack(1:lngth)
    l(jtot,jlp,1:lngth) = lpack(1:lngth)
    j12(jtot,jlp,1:lngth) = j12q(1:lngth)
    sr(jtot,jlp,1:len2) = sreal(1:len2)
    si(jtot,jlp,1:len2) = simag(1:len2)
    if(jtot==jfinl .and. jlpar==-1) exit
  enddo

  ! Compute spins
  call compute_spins(nucspin, flaghf, spins)

  ! Count and fill hyperfine levels for which XS are to be calculated
  call fill_hf(nlevel, jlev, elev, inlev, j1min, j2max, ered, twmol, spins, hf1, hf2)

  ! Determine the type of collision and call corresponding subroutine to compute squared T-matrix elements
  if (.not. twmol) then 
    if(nucspin < 100) then ! Molecule-Atom with 1 nuclear spin
      call molecule_atom_1spin(jfrst, jfinl, nlevel, jlev, l, j, in, length, sr, si,  spins, hf1, sigma)
    else ! Molecule-Atom with 2 nuclear spins
      call molecule_atom_2spin()
    endif
  else ! Molecule-Molecule with 1 nuclear spin
      call molecule_molecule_1spin()
  endif

  ! Compute hyperfine XS
  call compute_xs(twmol, rmu, ered, spins, hf1, hf2, sigma)
  ! Print hyperfine cross sections
  call print_xs(twmol, hfxfil_unit, ered, spins, hf1, hf2, sigma)
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
  type(spin_data), intent(out) :: spins
  
 
  if(flaghf) spins%f = 0.5d0 
  if (flaghf .and. nucspin.eq.2*(nucspin/2)) spins%h = 0.5d0
  if (.not.flaghf .and. nucspin.ne.2*(nucspin/2)) spins%h = 0.5d0

  if(nucspin<100) then
    spins%two = .false.
    spins%nuc1 = nucspin/2d0
    spins%i1 = int(2*spins%nuc1)
  else
    spins%two = .true.
    spins%nuc1 = (nucspin/100)/2.d0
    spins%i1 = int(2*spins%nuc1)
    spins%nuc2 = mod(nucspin,100)/2.d0
    spins%i2 = int(2*spins%nuc2)
  endif
  
end subroutine compute_spins



!************************************************************************************************
! This subroutine compute the number of hyperfine levels and fills the hyperfine arrays
!************************************************************************************************
subroutine fill_hf(nlevel, jlev, elev, inlev, j1min, j2max, ered, twmol, spins, hf1, hf2)
  implicit none
  ! Arguments
  integer, intent(in) :: nlevel
  integer, intent(in) :: jlev(*)
  real(8), intent(in) :: elev(*)
  integer, intent(in) :: inlev(*)
  integer, intent(in) :: j1min, j2max
  real(8), intent(in) :: ered
  logical, intent(in) :: twmol
  type(spin_data), intent(in) :: spins
  type(hflvl), intent(out) :: hf1, hf2
  ! Local variables
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
          ffmin = abs(fj - spins%nuc1)
          ffmax = fj + spins%nuc1
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

  ! THIS IS FOR TWO SPINS
  ! First pass will just count and allocate 
  ! Second pass will fill arrays

  if(spins%two) then
    fill = .false.
    do
      n2 = 0
      do i = 1, n
        if (hf1%e(i)<= ered) then
          f1 = hf1%if(i) + spins%f
          fnmin = abs(f1 - spins%nuc2)
          fnmax = f1 + spins%nuc2
          nhyp = int(fnmax - fnmin + 1)
          do ii = 1,  nhyp
            n2 = n2 + 1
            if(fill) then
              hf2%j(n2) = hf1%j(i)
              hf2%in(n2) = hf1%in(i)
              f2 = fnmin + (ii - 1)
              hf2%if(n2) = int(f2)
              hf2%e(n2) = hf1%e(i)
            endif
          enddo
        endif
      enddo
      if(.not. fill) then
        hf2%n = n2
        allocate(hf2%e(hf2%n), hf2%j(hf2%n), hf2%in(hf2%n), hf2%if(hf2%n))
        fill = .true.
      else
        exit
      endif
    enddo
  endif

  return
end subroutine fill_hf





!************************************************************************************************
! This subroutine computes T matrix elements for molecule-atom collisions with one spin
!************************************************************************************************
subroutine molecule_atom_1spin(jfrst, jfinl, nlevel, jlev, l, j, in, length, sr, si,  spins, hf1, sigma)
  use mod_hiutil, only: xf6j
  implicit none
  ! Arguments
  integer, intent(in) :: jfrst, jfinl, nlevel
  integer, intent(in) :: jlev(*), l(0:jfinl,2,*), j(0:jfinl,2,*), in(0:jfinl,2,*), length(0:jfinl, 2)
  real(8), intent(in) :: sr(0:jfinl,2,*) , si(0:jfinl,2,*)
  type(spin_data), intent(in) :: spins
  type(hflvl), intent(in) :: hf1
  real(8), allocatable, intent(out) :: sigma(:,:)
  ! Local variables
  integer :: jmx, idim
  ! Used within the loops
  real(8), allocatable :: tmatr(:,:), tmati(:,:)
  integer :: iftot, iftmn, iftmx, jlp, jlpar, i, ii, is
  real(8) :: xftot, xj, xjp, xf, xfp, fffp, xjttmn, xjtot, fjtot, xl, xlp, t2sum
  integer :: jttmin, jttmax, jtot, jlparf, jlpt, irow, icol, iflag, iph, ll, lp, phase
  complex(8) :: t


  ! Boundaries for loop over iftop   
  iftmn = max(0,int(jfrst + spins%f - spins%nuc1 - spins%h))
  iftmx = int(jfinl + spins%f + spins%nuc1 - spins%h)

  ! Allocate sigma array
  allocate(sigma(hf1%n, hf1%n))
  sigma = 0d0
  ! Allocate work arrays
  jmx = maxval(jlev(1:nlevel))
  idim = (iftmx + 2*jmx + 3) * (iftmx + jmx + 2)
  allocate(tmatr(idim, idim))
  allocate(tmati(idim, idim))

  do iftot = iftmn, iftmx
    xftot = iftot + spins%h
    do jlp = 1, 2
      jlpar = 1 - (jlp -1)*2
      write(6,"(A,F5.1,A,I4)") ' Computing partial wave J_tot =', xftot, ', jlpar=', jlpar
      do i = 1, hf1%n
        do ii = i, hf1%n
          xj = hf1%j(i) + spins%f
          xjp = hf1%j(ii) + spins%f
          xf = hf1%if(i) + spins%h
          xfp = hf1%if(ii) + spins%h
          fffp = sqrt((2.d0*xf+1.d0)*(2.d0*xfp+1.d0))
          xjttmn = max(spins%h, xftot - spins%nuc1)
          jttmin = max(jfrst, int(xjttmn-spins%f))
          jttmax = min(jfinl, int(xftot + spins%nuc1))
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
            if (length(jtot,jlpt) .gt. 0) then
              do irow = 1, length(jtot,jlpt)
  !     flag to make sure initial level is the bra, final level the ket
                iflag = 0
                if (j(jtot,jlpt,irow) == hf1%j(i) .and. in(jtot,jlpt,irow) == hf1%in(i)) then
                  iflag = 1
                endif
                do icol = 1, irow
                  if (j(jtot,jlpt,irow).eq.hf1%j(i) .and. &
                      in(jtot,jlpt,irow).eq.hf1%in(i) .and. &
                      j(jtot,jlpt,icol).eq.hf1%j(ii) .and. &
                      in(jtot,jlpt,icol).eq.hf1%in(ii) .or. &
                      j(jtot,jlpt,icol).eq.hf1%j(i) .and. &
                      in(jtot,jlpt,icol).eq.hf1%in(i) .and. &
                      j(jtot,jlpt,irow).eq.hf1%j(ii) .and. &
                      in(jtot,jlpt,irow).eq.hf1%in(ii)) then
                    is = (irow*(irow - 1))/2 + icol
                    if (iflag.ne.1) then
                      xl = l(jtot,jlpt,icol)
                      xlp = l(jtot,jlpt,irow)
                    else
                      xl = l(jtot,jlpt,irow)
                      xlp = l(jtot,jlpt,icol)
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
                        * xf6j(spins%nuc1,xj,xf,xl,xftot,xjtot) &
                      * xf6j(spins%nuc1,xjp,xfp,xlp,xftot,xjtot)
                    ll = int(xl)
                    lp = int(xlp)
                    tmatr(ll+1,lp+1) = tmatr(ll+1,lp+1) + real(t)
                    tmati(ll+1,lp+1) = tmati(ll+1,lp+1) + aimag(t)
  !     for initial level = final level, but l.ne.lp, need to include
  !     both T(l,lp) and T(lp,l)
                    if (hf1%j(i) == hf1%j(ii) .and. hf1%in(i) == hf1%in(ii) .and. abs(xl-xlp) > 1d-60) then
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
  deallocate(tmatr)
  deallocate(tmati)
end subroutine molecule_atom_1spin




!************************************************************************************************
! This subroutine computes T matrix elements for molecule-atom collisions with two spins
!************************************************************************************************
subroutine molecule_atom_2spin()
  implicit none
end subroutine molecule_atom_2spin



!************************************************************************************************
! This subroutine computes T matrix elements for molecule-molecule collisions with one spin
!************************************************************************************************
subroutine molecule_molecule_1spin()
  implicit none
end subroutine molecule_molecule_1spin


!************************************************************************************************
! This subroutine computes XS from T matrix elements 
!************************************************************************************************
subroutine compute_xs(twmol, rmu, ered, spins, hf1, hf2, sigma)
  use constants, only: ang2 => ang2c
  implicit none
  ! Arguments
  logical, intent(in) :: twmol
  real(8), intent(in) :: rmu, ered
  type(spin_data), intent(in) :: spins
  type(hflvl), intent(in) :: hf1, hf2
  real(8), intent(inout) :: sigma(hf1%n, hf1%n)
  ! Local variables
  real(8) :: fak, dencol, denrow, ffi, fff
  integer :: ij2, ij2p, i, ii

  fak = acos(-1.d0) * ang2 / (2.0d0 * rmu)
  do i = 1, hf1%n
    ij2 = 0.d0
    if (twmol) ij2 = mod(hf1%j(i),10)
    ffi = hf1%if(i) + spins%h
    denrow = (2.d0 * ffi + 1.d0) * (2.d0 * ij2 + 1.d0) * (ered - hf1%e(i))
    do ii = i, hf1%n
      ij2p = 0.d0
      if (twmol) ij2p = mod(hf1%j(ii),10)
      fff = hf1%if(ii) + spins%h
      dencol = (2.d0 * fff + 1.d0) * (2.d0 * ij2p + 1.d0) * (ered - hf1%e(ii))
      sigma(i,ii) = sigma(i,ii) * fak / denrow
      if (i.ne.ii) sigma(ii,i) = sigma(ii,i) * fak / dencol
   end do
end do  
end subroutine compute_xs


!************************************************************************************************
! This subroutine prints XS 
!************************************************************************************************
subroutine print_xs(twmol, hfxfil_unit, ered, spins, hf1, hf2, sigma)
  use constants, only: econv
  implicit none
  ! Arguments
  logical, intent(in) :: twmol
  integer, intent(in) :: hfxfil_unit
  real(8), intent(in) :: ered
  type(spin_data), intent(in) :: spins
  type(hflvl), intent(in) :: hf1, hf2
  real(8), intent(in) :: sigma(hf1%n,hf1%n)
  ! Local variables
  integer :: i, ii, ij2, ij2p
  real(8) :: xj, xf, xjp, xfp, ee

  ! Write header
  write(6,"(a)") '%     E(CM-1)     JI     INI   FI      JF     INF   FF      CROSS SECTION (ANG^2)'
  write(hfxfil_unit,"(a)") '%     E(CM-1)     JI     INI   FI      JF     INF   FF      CROSS SECTION (ANG^2)'

  ee = ered*econv
  do i = 1, hf1%n
    do ii = 1, hf1%n
      if (.not. twmol) then
        xj = hf1%j(i) + spins%f
        xf = hf1%if(i) + spins%h
        xjp =  hf1%j(ii) + spins%f
        xfp = hf1%if(ii) + spins%h
        if (sigma(i,ii)>0d0) then
          write(6,"(f12.3,f8.1,i6,f6.1,3x,f6.1,i6,f6.1,5x,1pe15.4)")&
                ee,xj,hf1%in(i),xf,xjp,hf1%in(ii),xfp,sigma(i,ii)
          write(hfxfil_unit,"(f12.3,f8.1,i6,f6.1,3x,f6.1,i6,f6.1,5x,1pe15.4)")&
                ee,xj,hf1%in(i),xf,xjp,hf1%in(ii),xfp,sigma(i,ii)
        end if
      else
        xj = (hf1%j(i)/10) + spins%f
        ij2 = mod(hf1%j(i),10)
        xf  = hf1%if(i) + spins%h
        xjp = (hf1%j(ii)/10) + spins%f
        ij2p = mod(hf1%j(ii),10)
        xfp  = hf1%if(ii) + spins%h
        if (sigma(i,ii)>0d0) then
          write(6,"(f12.3,f8.1,i6,f6.1,i6,3x,f6.1,i6,f6.1,i6,5x,1pe15.4)")&
                ee,xj,hf1%in(i),xf,ij2,xjp,hf1%in(ii),xfp,ij2p,sigma(i,ii)
          write(hfxfil_unit,"(f12.3,f8.1,i6,f6.1,i6,3x,f6.1,i6,f6.1,i6,5x,1pe15.4)")&
                ee,xj,hf1%in(i),xf,ij2,xjp,hf1%in(ii),xfp,ij2p,sigma(i,ii)
        endif
      endif
    enddo
  enddo
end subroutine print_xs


end module mod_hypxsc