!*******************************************************************************
!     This module contains all the procedures related to S-Matrices
!     for hyperfine cross section calculations
!*******************************************************************************
module mod_hyp_smat
    implicit none

    ! Protected variables
    real(8), allocatable, protected :: elev(:)      ! Levels energies 
    integer, allocatable, protected :: jlev(:)      ! Levels rotational angular momenta
    integer, allocatable, protected :: inlev(:)     ! Levels additional quantum numbers
    real(8), allocatable, protected :: sr(:,:,:)    ! Real part of the S-Matrices for each Jtot, paritiy, and channel
    real(8), allocatable, protected :: si(:,:,:)    ! Imaginary part of the S-Matrices for each Jtot, parity, and channel
    integer, allocatable, protected :: j(:,:,:)     ! Channel rotational angular momenta
    integer, allocatable, protected :: l(:,:,:)     ! Channel orbital angular momenta
    integer, allocatable, protected :: in(:,:,:)    ! Channel additional quantum index
    integer, allocatable, protected :: j12(:,:,:)   ! Channel total molecular rotational quantum number 
    integer, allocatable, protected :: length(:,:)  ! Number of channels for each Jtot and parity
    integer,              protected :: jfrst, jfinl ! First and last Jtot values
    integer,              protected :: jtotd        ! Step between Jtot values
    integer,              protected :: nnout        ! Length of JOUT array
    real(8),              protected :: ered         ! Total collision energy
    real(8),              protected :: rmu          ! Reduced mass of the system
    integer,              protected :: nucspin      ! Nuclear spin 
    logical,              protected :: csflg        ! True if CS calculations
    logical,              protected :: flaghf       ! True if total angular momentum is equal to JTOT+1/2
    logical,              protected :: flgsu        ! Ture if surfrace collision
    logical,              protected :: twmol        ! True if two molecules, false if atom-molecule

    ! Private variables 
    character*20,           private :: cdate        ! Date of creation of the S-Matrices file
    character*48,           private :: label        ! Label of the calculation
    character*48,           private :: potnam       ! Potential name used
    character*40,           private :: smtfil       ! Filename of the S-Matrices file

    ! Procedures contained in this module
    public  :: read_s_data
    public  :: print_s_infos
    public :: deallocate_module_arrays
    private :: allocate_module_arrays
    contains

    !************************************************************************************************
    ! This subroutine prints calculations info read from the S-matrices file
    !************************************************************************************************
    subroutine print_s_infos(iunit)
        implicit none
        ! Arguments
        integer, intent(in) :: iunit

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
    end subroutine print_s_infos

    
    !************************************************************************************************
    ! This subroutine reads data from the S-matrices file
    !************************************************************************************************
    subroutine read_s_data(flname, iener)
        use mod_hiutil, only: gennam
        use mod_hismat, only: sread, rdhead, sinqr
        use mod_coj12, only: j12q => j12
        implicit none
        ! Arguments
        character*(*), intent(in) :: flname ! Job name from which smat filename is generated
        integer,       intent(in) :: iener  ! Collision energy index needed
        ! Local variables
        integer :: dum_int
        logical :: dum_log
        integer :: mjtot, mchmx, nlevel, jlpar, jlp, len, len2, ierr
        integer, allocatable :: jout(:)
        real(8), allocatable :: sreal(:), simag(:)
        integer, allocatable :: jq(:), lq(:), inq(:), ipack(:), jpack(:), lpack(:)
        ! Generate filename of smt file and check if it is present
        call gennam(smtfil, flname, iener, 'smt', lend)
        inquire(file = smtfil(1:lend), exist = exstfl)
        if (.not. exstfl) then ; write(6,"(3a)") '*** FILE ', smtfil(1:lend), ' NOT FOUND ***' ; stop ; endif

        ! Open S-matrix file and read headers
        call openf(1, smtfil, 'tu', 0)
        call sinqr(1, mjtot, mchmx)
        ! Allocate local arrays
        allocate(jout(50))
        allocate(jq(mchmx), lq(mchmx), inq(mchmx), ipack(mchmx), jpack(mchmx), lpack(mchmx))
        allocate(sreal(mchmx2), simag(mchmx2))
        ! Init local arrays
        jout = 0 ; sreal = 0d0 ; simag = 0d0
        jq = 0 ; lq = 0 ; inq = 0 ; ipack = 0 ; jpack = 0 ; lpack = 0 ;
        ! Continue reading headers
        call rdhead(1, cdate, ered, rmu, &
                       csflg, flaghf, flgsu, twmol, &
                       dum_log, &
                       jfrst, jfinl, jtotd, &
                       dum_int, dum_int, dum_int, &
                       nlevel, &
                       dum_int, &
                       nnout, &
                       jlev, inlev, elev, jout)

        ! Allocate module arrays
        call allocate_module_arrays(jfinl, mchmx)

        ! Read all S-Matrices
          do
            call sread (0, sreal, simag, jtot, jlpar, dum_int, jq, lq, inq, ipack, jpack, lpack, 1, mchmx, -1, len, ierr)
            if(ierr < -1) then 
                write(6,*) '*** READ ERROR IN HYPXSC. ABORT ***'
                stop 
            endif

            jlp = 1 - (jlpar-1)/2
            length(jtot, jlp) = len 
            len2 = len*(len + 1)/2

            j(jtot,jlp,1:len) = jpack(1:len)
            in(jtot,jlp,1:len) = ipack(1:len)
            l(jtot,jlp,1:len) = lpack(1:len)
            sr(jtot,jlp,1:len2) = sreal(1:len2)
            si(jtot,jlp,1:len2) = simag(1:len2)
            j12(jtot,jlp,1:len) = j12q(1:len)
            
            if(jtot==jfinl .and. jlpar==-1) exit
        enddo
        
        deallocate(jlev,inlev,elev,jout)
        deallocate(jlev,inlev,elev,jout)

    end subroutine read_s_data



    !************************************************************************************************
    ! This subroutine allocates module arrays
    !************************************************************************************************
    subroutine allocate_module_arrays(jfinl, mchmx)
        implicit none
        ! Arguments
        integer, intent(in) :: jfinl, mchmx
        ! Local variavbles
        integer, mchmx2

        mchmx2 = mchmx * (mchmx + 1) / 2

        allocate(  elev(mchmx) )           ; elev = 0d0
        allocate(  jlev(mchmx) )           ; jlev = 0
        allocate( inlev(mchmx) )           ; inlev = 0
        allocate( sr(0:jfinl, 2, mchmx2) ) ; sr = 0d0
        allocate( si(0:jfinl, 2, mchmx2) ) ; si = 0d0
        allocate(  j(0:jfinl, 2, mchmx) )  ; j = 0
        allocate( in(0:jfinl, 2, mchmx) )  ; in = 0
        allocate(  l(0:jfinl, 2, mchmx) )  ; l = 0
        allocate( length(0:jfinl, 2) )     ; length = 0
        allocate( j12(0:jfinl, 2, mchmx) ) ; j12 = 0

    end subroutine allocate_module_arrays

    !************************************************************************************************
    ! This subroutine deallocates module arrays
    !************************************************************************************************
    subroutine deallocate_module_arrays()
        implicit none

        deallocate(elev)
        deallocate(jlev)
        deallocate(inlev)
        deallocate(sr)
        deallocate(si)
        deallocate(j)
        deallocate(in)
        deallocate(l)
        deallocate(length)
        deallocate(j12)

    end subroutine deallocate_module_arrays

end module mod_hyp_smat