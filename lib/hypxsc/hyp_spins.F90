!*******************************************************************************
!     This module contains all the procedures necessary to compute 
!     spins
!*******************************************************************************
module mod_hyp_spins
    implicit none

    ! Protected variables
    integer, protected :: inuc1
    integer, protected :: inuc2
    real(8), protected :: fnuc1
    real(8), protected :: fnuc2
    real(8), protected :: fspin
    real(8), protected :: fhspin

    public  :: 
    
    contains

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
end module mod_hyp_spins