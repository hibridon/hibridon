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
!     Refactored, rewritten and parallelized by b. desrousseaux (oct. 2022)
!
!*******************************************************************************
    
    module mod_hypxsc
    implicit none

    !************************************************************************************************
    ! Main subroutine that drives hyperfine cross sections calculations 
    !************************************************************************************************
    subroutine hypxsc(flname, a) 
        use mod_selb, only: ibasty
        use mod_hyp_checks, only: check_if_supported
        use mod_hyp_smat,   only: read_s_data, print_s_infos, flgsu, csflg, jtotd, nnout, twmol
        implicit none
        ! Arguments
        character*(*), intent(in) :: flname
        real(8), dimension(4), intent(in) :: a(4)
        ! Local variables
        integer :: iener, nucspin, j1min, j2max

        ! Define variables from hypxsc command arguments
        iener   = int(a(1))
        nucspin = int(a(2))
        j1min   = int(a(3))
        j2max   = int(a(4))


        ! Generate output file name and open it
        call gennam(hfxfil,flname,iener,'hfx',lend)
        call openf(hfxfil_unit, hfxfil(1:lend),'sf',0)

        ! Read data from the S-Matrices files
        call read_s_data(flname, iener)

        ! Check if hypxsc supports the basis type and parameters read in the S-Matrices file
        call check_if_supported(ibasty, flgsu, csflg, jtotd, nnout, nucspin, twmol)

        ! Write infos to stdout and output file
            
        ! Compute spins

        ! Compute hypfine levels

        ! Compute cross sections

        ! Print hyperfine cross sections




    end subroutine hypxsc