module mod_hypxsc
implicit none



contains 

!************************************************************************************************
! Main subroutine that drives hyperine cross sections calculations and printing
!************************************************************************************************
subroutine hypxsc(flname, a)
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
  ! Arguments
  character*(*), intent(in) :: flname
  real(8), dimension(4), intent(in) :: a(4)
  ! Local variables
  character*40 :: smtfil, hfxfil
  integer :: iener, lend, hfxfil_unit, nucspin, len2, mchmx2, j1min, j2max, nlevelh
  logical :: supported, exstfl
  real(8) :: ee, fspin, fhspin, finuc, f2nuc
  real(8), allocatable :: sigma(:,:)
  ! Read from S-matrix file:
  integer :: mjtot, mchmx, lngth, ierr, jtot, jlpar, nu, nopen, jlp
  character*20 :: cdate
  integer :: jfrst, jfinl, jtotd, numin, numax, nud, nlevel, nlvop, nnout
  real(8) :: ered, rmu
  logical :: nucrs, csflg, flaghf, flgsu, twmol
  real(8), allocatable :: sreal(:), simag(:)
  real(8), allocatable :: sr(:,:,:), si(:,:,:)
  integer, allocatable :: j(:,:,:), in(:,:,:), l(:,:,:), j12(:,:,:), length(:,:)
  logical, allocatable :: exsmtn(:), exsmtp(:)

  !call mtime(cpu0, ela0)! initialize timer 

  ! Generate filename of smt file and check if it is present
  iener = int(a(1))
  call gennam(smtfil, flname, iener, 'smt', lend)
  inquire(file = smtfil(1:lend), exist = exstfl)
  if (.not. exstfl) then
    write(6,"(3a)") '*** FILE ', smtfil(1:lend), ' NOT FOUND ***'
    return
  endif

  ! Open S-matrix file and read its header
  call openf(1, smtfil, 'tu', 0)
  call sinqr(1, mjtot, mchmx)
  call rdhead(1,cdate,ered,rmu,csflg,flaghf,flgsu, &
     twmol,nucrs,jfrst,jfinl,jtotd,numin,numax,nud, &
     nlevel,nlvop,nnout,jlev,inlev,elev,jout)

  ! Define useful variables
  ee = ered*econv ! Convert collision energy
  nucspin = int(a(2)) ! nucspin is nuclear spin I times 2
  j1min = int(a(3))
  j2max = int(a(4))

  ! Check if hypxsc supports the basis type and parameters provided in the input
  call check_if_supported(ibasty, flgsu, csflg, jtotd, nnout, nucspin, twmol, supported)
  if (.not. supported) return

  ! Generate output file name and open it
  call gennam(hfxfil,flname,iener,'hfx',lend)
  call openf(hfxfil_unit, hfxfil(1:lend),'sf',0)

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
    call compute_spins(nucspin, flaghf, fspin, finuc, f2nuc, fhspin)

  ! Determine the type of collision and call corresponding subroutine
  if (.not. twmol) then 
    if(nucspin < 100) then ! Molecule-Atom with 1 nuclear spin
      call molecule_atom_1spin(nucspin, fspin, finuc, fhspin, j1min, j2max, nlevel, jlev, elev, inlev, ered,&
                               jfrst, jfinl, flaghf, rmu, l, in, j, length, sr, si, nlevelh, sigma)
    else ! Molecule-Atom with 2 nuclear spins
      call molecule_atom_2spin()
    endif
  else ! Molecule-Molecule with 1 nuclear spin
      call molecule_molecule_1spin()
  endif
end subroutine hypxsc



!************************************************************************************************
! This subroutine checks if hypxsc is implemented for user's case
!************************************************************************************************
subroutine check_if_supported(ibasty, flgsu, csflg, jtotd, nnout, nucspin, twmol, supported)
  implicit none 
  integer, intent(in) :: ibasty, jtotd, nnout, nucspin
  logical, intent(in) :: flgsu, csflg, twmol
  logical, intent(out):: supported 

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
end subroutine check_if_supported



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
subroutine compute_spins(nucspin, flaghf, fspin, finuc, f2nuc, fhspin)
  implicit none
  ! Arguments
  integer, intent(in) :: nucspin
  logical, intent(in) :: flaghf
  real(8), intent(out) :: fspin, fhspin, finuc, f2nuc
  
  ! fspin
  fspin = 0d0
  if(flaghf) fspin = 0.5d0
  ! finuc and f2nuc
  if(nucspin<100) finuc = nucspin/2d0
  if(nucspin>100) then
    finuc = (nucspin/100)/2.d0
    f2nuc = mod(nucspin,100)/2.d0
  endif
  ! fhspin
  fhspin = 0.d0
  if (flaghf .and. nucspin.eq.2*(nucspin/2)) fhspin = 0.5d0
  if (.not.flaghf .and. nucspin.ne.2*(nucspin/2)) fhspin = 0.5d0
end subroutine compute_spins



!************************************************************************************************
! This subroutine compute the number of hyperfine levels or fills the hyperfine arrays (jlevh, inlevh, iflevh, elevh)
!************************************************************************************************
subroutine fill_levelh(nlevel, jlev, elev, inlev, j1min, j2max, ered, fspin, finuc, twmol, fill,&
                       nlevelh, jlevh, inlevh, iflevh, elevh)
  implicit none
  ! Arguments
  integer, intent(out) :: nlevelh
  real(8), intent(in) :: elev(*), ered, fspin, finuc
  integer, intent(in) :: jlev(*), inlev(*), nlevel, j1min, j2max
  logical, intent(in) :: twmol, fill
  real(8), allocatable, intent(inout) :: elevh(:)
  integer, allocatable, intent(inout) :: jlevh(:), inlevh(:), iflevh(:)
  ! Local variables
  integer :: i, ii, ij1, nhyp
  real(8) :: ff, fj, ffmin, ffmax

  nlevelh = 0
  do i = 1, nlevel 
    if(twmol) then ; ij1 = jlev(i)/10 ; else ; ij1 = jlev(i) ; endif
    if(ij1 >= j1min .and. ij1 <= j2max) then
      if(elev(i) <= ered) then
        fj = ij1 + fspin
        ffmin = abs(fj - finuc)
        ffmax = fj + finuc
        nhyp = int(ffmax - ffmin + 1)
        do ii=1,nhyp
          nlevelh = nlevelh + 1
          if (fill) then
            jlevh(nlevelh) = jlev(i)
            inlevh(nlevelh) = inlev(i)
            ff = ffmin + (ii - 1)
            iflevh(nlevelh) = int(ff)
            elevh(nlevelh) = elev(i)
          endif
        end do
      endif
    endif
  enddo
end subroutine fill_levelh


!************************************************************************************************
! This subroutine computes T matrix elements for molecule-atom collisions with one spin
!************************************************************************************************
subroutine molecule_atom_1spin(nucspin, fspin, finuc, fhspin, j1min, j2max, nlevel, jlev, elev, inlev, ered,&
                               jfrst, jfinl, flaghf, rmu, l, in, j, length, sr, si, nlevelh, sigma)
  use mod_hiutil, only: xf6j
  implicit none
  ! Arguments
  integer, intent(in) :: jlev(*), inlev(*), l(0:jfinl,2,*), j(0:jfinl,2,*), in(0:jfinl,2,*), length(0:jfinl, 2)
  integer, intent(in) :: nucspin, nlevel, jfrst, jfinl, j1min, j2max
  real(8), intent(in) :: elev(*), sr(0:jfinl,2,*) , si(0:jfinl,2,*), ered, rmu, fspin, finuc, fhspin
  logical, intent(in) :: flaghf
  integer, intent(out) :: nlevelh
  real(8), allocatable, intent(inout) :: sigma(:,:)
  ! Local variables
  integer :: jmx, idim
  real(8), allocatable :: tmatr(:,:), tmati(:,:), elevh(:)
  integer, allocatable :: jlevh(:), inlevh(:), iflevh(:)
  ! Used within the loops
  integer :: iftot, iftmn, iftmx, jlp, jlpar, i, ii, is
  real(8) :: xftot, xj, xjp, xf, xfp, fffp, xjttmn, xjtot, fjtot, xl, xlp, t2sum
  integer :: jttmin, jttmax, jtot, jlparf, jlpt, irow, icol, iflag, iph, ll, lp, phase
  complex(8) :: t
  ! For XS
  integer :: ij1, ij2, ij2p
  real(8) :: ffi, fff, fak, dencol, denrow

  ! Boundaries for loop over iftop   
  iftmn = max(0,int(jfrst + fspin - finuc - fhspin))
  iftmx = int(jfinl + fspin + finuc - fhspin)

  ! Count and set up list of hyperfine levels for which XS are to be calculated
  call fill_levelh(nlevel, jlev, elev, inlev, j1min, j2max, ered, fspin, finuc, .false., .false.,&
                   nlevelh, jlevh, inlevh, iflevh, elevh)
  allocate(elevh(nlevelh), jlevh(nlevelh), inlevh(nlevelh), iflevh(nlevelh))
  call fill_levelh(nlevel, jlev, elev, inlev, j1min, j2max, ered, fspin, finuc, .false., .true.,&
                   nlevelh, jlevh, inlevh, iflevh, elevh)

  ! Allocate sigma array
  allocate(sigma(nlevelh,nlevelh))
  sigma = 0d0
  ! Allocate work arrays
  jmx = maxval(jlev(1:nlevel))
  idim = (iftmx + 2*jmx + 3) * (iftmx + jmx + 2)
  allocate(tmatr(idim, idim))
  allocate(tmati(idim, idim))

  do iftot = iftmn, iftmx
    xftot = iftot + fhspin
    do jlp = 1, 2
      jlpar = 1 - (jlp -1)*2
      write(6,"(A,F5.1,A,I4)") ' Computing partial wave J_tot =', xftot, ', jlpar=', jlpar
      do i=1,nlevelh
        do ii=i,nlevelh
          xj = jlevh(i) + fspin
          xjp = jlevh(ii) + fspin
          xf = iflevh(i) + fhspin
          xfp = iflevh(ii) + fhspin
          fffp = sqrt((2.d0*xf+1.d0)*(2.d0*xfp+1.d0))
          xjttmn = max(fhspin, xftot - finuc)
          jttmin = max(jfrst, int(xjttmn-fspin))
          jttmax = min(jfinl, int(xftot + finuc))
  !     clear T-matrix array
          tmatr = 0.d0
          tmati = 0.d0
  !     sum over jtot consistent with vector addition
  !     ftot = jtot + nucspin
          do jtot=jttmin,jttmax
            xjtot = jtot + fspin
            fjtot = 2.d0 * xjtot + 1.d0
  !     total parity must be the same for all T-matrix elements
  !     in sum over jtot (for the same parity, jlparf
  !     changes sign for each increase in jtot by 1)
            jlparf = jlpar*(-1)**(iftot - jtot)
            jlpt = 1 + (1 - jlparf)/2
            if (length(jtot,jlpt) .gt. 0) then
              do irow=1,length(jtot,jlpt)
  !     flag to make sure initial level is the bra, final level the ket
                iflag = 0
                if (j(jtot,jlpt,irow).eq.jlevh(i) .and. &
                    in(jtot,jlpt,irow).eq.inlevh(i)) then
                  iflag = 1
                end if
                do icol=1,irow
                  if (j(jtot,jlpt,irow).eq.jlevh(i) .and. &
                      in(jtot,jlpt,irow).eq.inlevh(i) .and. &
                      j(jtot,jlpt,icol).eq.jlevh(ii) .and. &
                      in(jtot,jlpt,icol).eq.inlevh(ii) .or. &
                      j(jtot,jlpt,icol).eq.jlevh(i) .and. &
                      in(jtot,jlpt,icol).eq.inlevh(i) .and. &
                      j(jtot,jlpt,irow).eq.jlevh(ii) .and. &
                      in(jtot,jlpt,irow).eq.inlevh(ii)) then
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
                        * xf6j(finuc,xj,xf,xl,xftot,xjtot) &
                      * xf6j(finuc,xjp,xfp,xlp,xftot,xjtot)
                    ll = int(xl)
                    lp = int(xlp)
                    tmatr(ll+1,lp+1) = tmatr(ll+1,lp+1) + real(t)
                    tmati(ll+1,lp+1) = tmati(ll+1,lp+1) + aimag(t)
  !     for initial level = final level, but l.ne.lp, need to include
  !     both T(l,lp) and T(lp,l)
                    if (jlevh(i).eq.jlevh(ii) .and. inlevh(i).eq.inlevh(ii) .and. abs(xl-xlp)<1d-60) then
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
          if (i.ne.ii) then
            sigma(ii,i) = sigma(ii,i) + t2sum * (2.d0 * xftot + 1.d0)
          end if
        end do
      end do
    end do
  enddo
  deallocate(tmatr)
  deallocate(tmati)

  ! Compute hyperfine XS
  call compute_xs(nlevelh, rmu, .false., ered, fhspin, iflevh, elevh, jlevh, sigma)
  ! Print hyperfine cross sections
  call print_xs(.false., 3, ered, rmu, nlevelh, elevh, jlevh, iflevh, inlevh, fspin, fhspin, sigma)
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
subroutine compute_xs(nlevelh, rmu, twmol, ered, fhspin, iflevh, elevh, jlevh, sigma)
  use constants, only: ang2 => ang2c
  implicit none
  ! Arguments
  integer, intent(in) :: nlevelh, iflevh(nlevelh), jlevh(nlevelh)
  real(8), intent(in) :: elevh(nlevelh), rmu, fhspin, ered
  logical, intent(in) :: twmol
  real(8), intent(inout) :: sigma(nlevelh,nlevelh)
  ! Local variables
  real(8) :: fak, dencol, denrow, ffi, fff
  integer :: ij12, ij2, ij2p, i, ii

  fak = acos(-1.d0) * ang2 / (2.0d0 * rmu)
  do i=1,nlevelh
    ij2 = 0.d0
    if (twmol) ij2 = mod(jlevh(i),10)
    ffi = iflevh(i) + fhspin
    denrow = (2.d0 * ffi + 1.d0) * (2.d0 * ij2 + 1.d0) * (ered - elevh(i))
    do ii=i,nlevelh
      ij2p = 0.d0
      if (twmol) ij2p = mod(jlevh(ii),10)
      fff = iflevh(ii) + fhspin
      dencol = (2.d0 * fff + 1.d0) * (2.d0 * ij2p + 1.d0) * (ered - elevh(ii))
      sigma(i,ii) = sigma(i,ii) * fak / denrow
      if (i.ne.ii) sigma(ii,i) = sigma(ii,i) * fak / dencol
   end do
end do  
end subroutine compute_xs


!************************************************************************************************
! This subroutine prints XS 
!************************************************************************************************
subroutine print_xs(twmol, hfxfil_unit, ered, rmu, nlevelh, elevh, jlevh, iflevh, inlevh, fspin, fhspin, sigma)
  use constants, only: econv, ang2 => ang2c
  implicit none
  ! Arguments
  integer, intent(in) :: nlevelh
  integer, intent(in) :: iflevh(nlevelh), jlevh(nlevelh), inlevh(nlevelh), hfxfil_unit
  real(8), intent(in) :: rmu, ered, fspin, fhspin, elevh(nlevelh), sigma(nlevelh,nlevelh)
  logical, intent(in) :: twmol
  ! Local variables
  integer :: i, ii, ij2, ij2p
  real(8) :: xj, xf, xjp, xfp, ee

  ee = ered*econv
  do i=1,nlevelh
    do ii=1,nlevelh
      if (.not. twmol) then
        xj = jlevh(i) + fspin
        xf = iflevh(i) + fhspin
        xjp = jlevh(ii) + fspin
        xfp = iflevh(ii) + fhspin
        if (sigma(i,ii)>0d0) then
          write(6,"(f12.3,f8.1,i6,f6.1,3x,f6.1,i6,f6.1,5x,1pe15.4)")&
                ee,xj,inlevh(i),xf,xjp,inlevh(ii),xfp,sigma(i,ii)
          !write(hfxfil_unit,"(f12.3,f8.1,i6,f6.1,3x,f6.1,i6,f6.1,5x,1pe15.4)")&
          !      ee,xj,inlevh(i),xf,xjp,inlevh(ii),xfp,sigma(i,ii)
        end if
      else
        xj = (jlevh(i)/10) + fspin
        ij2 = mod(jlevh(i),10)
        xf  = iflevh(i) + fhspin
        xjp = (jlevh(ii)/10) + fspin
        ij2p = mod(jlevh(ii),10)
        xfp  = iflevh(ii) + fhspin
        if (sigma(i,ii)>0d0) then
          write(6,"(f12.3,f8.1,i6,f6.1,i6,3x,f6.1,i6,f6.1,i6,5x,1pe15.4)")&
                ee,xj,inlevh(i),xf,ij2,xjp,inlevh(ii),xfp,ij2p,sigma(i,ii)
          !write(hfxfil_unit,"(f12.3,f8.1,i6,f6.1,i6,3x,f6.1,i6,f6.1,i6,5x,1pe15.4)")&
          !      ee,xj,inlevh(i),xf,ij2,xjp,inlevh(ii),xfp,ij2p,sigma(i,ii)
        endif
      endif
    enddo
  enddo
end subroutine print_xs


end module mod_hypxsc