!*******************************************************************************
!     This module contains all the procedures necessary to check 
!     the good execution of hyperfine cross section calculations
!*******************************************************************************
module mod_hyp_checks
    implicit none

    public  :: check_if_supported
    
    contains

    !************************************************************************************************
    ! This function checks if hypxsc is implemented for user's case and returns .true. or .false.
    !************************************************************************************************
    logical function check_if_supported(ibasty, flgsu, csflg, jtotd, nnout, nucspin, twmol)
        implicit none 
        integer, intent(in) :: ibasty, jtotd, nnout, nucspin
        logical, intent(in) :: flgsu, csflg, twmol
        
        check_if_supported = .true.
        ! CHECK IBASTY
        if(any((/10,12,13,15,22,23/)==ibasty)) then
            write(6,"(A,I0,A)") '*** HYPERFINE CROSS SECTIONS FOR BASIS TYPE =', ibasty,' NOT IMPLEMENTED ***'
            check_if_supported = .false.
        endif
        ! CHECK FLGSU
        if(flgsu) then
            write(6,"(A)") '*** HYPERFINE CROSS SECTIONS FOR SURFACE COLLISIONS NOT IMPLEMENTED ***'
            check_if_supported = .false.
        endif
        ! CHECK CSFLG
        if(csflg) then
            write(6,"(A)") '*** CS HYPERFINE CROSS SECTIONS NOT IMPLEMENTED ***'
            check_if_supported = .false.
        endif 
        ! CHECK JTOTD
        if (iabs(jtotd).ne.1) then
            write(6,"(A)") '*** DELTA-JTOT MUST BE EQUAL TO ONE ***'
            check_if_supported = .false.
        endif
        ! CHECK NNOUT
        if (nnout.lt.0) then
            write(6,"(A)") '*** NNOUT < 0, ABORT ***'
            check_if_supported = .false.
        endif
        ! CHECK TWMOL WITH 2 NUC SPIN
        if (twmol .and. nucspin > 100) then
            write(6,"(a)") '*** HYPERFINE CROSS SECTIONS FOR MOLECULE-MOLECULE COLLISIONS NOT IMPLEMENTED FOR TWO NUCLEAR SPINS'
            check_if_supported = .false.
        endif

        return
    end function check_if_supported





end module mod_hyp_checks