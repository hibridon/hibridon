module mod_hitestptn
    implicit none

    private
    public :: testptn
    

contains 

    subroutine testptn(ihomo)
        implicit none
        logical, intent(in) :: ihomo
        integer :: input

        menuloop: do
            write(6,*) "Enter 0 to quit, 1 for V(R,theta) 2 for V-lambda(R)"
            read(5,*) input
            select case (input)
                case (0) ; exit menuloop
                case (1) ; call v(ihomo)
                case (2) ; call vlam()
            end select       
        enddo menuloop

    end subroutine testptn

    subroutine vlam()
        use mod_hipot, only: driver
        implicit none
        write(6,*)
        write(6,*) "************************************************************"
        write(6,*) "      Entering the DRIVER subroutine of the potential       "
        write(6,*) "                   Press Ctrl+D to go back                  "
        write(6,*) "************************************************************"
        write(6,*)
        call driver
    end subroutine vlam

    subroutine v(ihomo)
        use mod_covvl, only: vvl
        use mod_cosysi, only: ispar
        use mod_parbas, only: lammin, lammax, mproj
        use mod_hipot, only: pot
        implicit none
        logical, intent(in) :: ihomo
        integer :: nterm, nstep
        integer :: n, i, iterm
        real(8) :: vv0, vv, theta, r

        if(.not. pot_loader()) return

        nterm = ispar(1)
        if(ihomo) then 
            nstep = 2
        else
            nstep = 1
        endif

        write(6,*)
        write(6,*) "************************************************************"
        write(6,*) "                         WARNING:                           "
        write(6,*) "    Meaningful values of V(R,) will be given only if        "
        write(6,*) "              the potential used is expanded as:            "
        write(6,*) "      V(R,theta) = V00(R) + V_lam(R) x P_lam(cos(theta))    "
        write(6,*) "************************************************************"
        write(6,*)

        menuloop: do
            write(6,*) "Enter R (bohr), theta (Degrees); R=0 to quit"
            read(5,*) r, theta
            if (r<1d-30 ) then
                return
            else
                call pot(vv0, r)   
                vv = vv0
                n = 1

                do iterm = 1, nterm
                do i=lammin(iterm), lammax(iterm), nstep
                        vv = vv + vvl(n) * dlm0(i,mproj(iterm),theta)
                        n=n+1
                enddo
                enddo

                write(6,"(A4,g13.6,A8)") "V = ", vv , " hartree"
                write(6,"(A4,g13.6,A5)") "V = ", vv*219474.6, " cm-1"
                write(6,"(A4,g13.6,A3)") "V = ", vv*27.211652, " eV"
            endif
        enddo menuloop

    end subroutine v


    logical function pot_loader()
        use mod_file, only: input
        use mod_covvl, only: vvl
        use mod_hipot, only: loapot
        implicit none
        integer ::  ios
        logical :: existf
        character(len=200) :: pot_data_file = ""
        character(len=1000) :: line
    ! Get potdata file by reading last line of input file (if any)
        ! Check that input file is provided
        inquire (file=input, exist=existf) 
        if(.not. existf) then
            write (6, *) "*** INPUT FILE ", trim(input), " DOES NOT EXIST"
            pot_loader = .false.
            return
        endif
        ! Read input file to get potdata filename
        open(unit=1,file=trim(input),status='old')
        do
            read(1,"(a)", iostat=ios) line
            if (ios /= 0) exit ! Exit at last line that should contain potdata filename
        enddo
        ! Check that the last line read is a potdata file. If not, that means this potential does not use a potdata file.
        if(index(line, ".dat") /= 0) then
            pot_data_file = adjustl(trim(line(1:index(line, ".dat")+4)))
        endif
        ! Load potential
        call loapot(1,trim(pot_data_file))
        vvl = 0d0
        pot_loader = .true.
        return
    end function pot_loader

    !function to calculate d_l^m(cos(theta) as defined in "brink and satchler"
    function dlm0(ll,m,theta)
        implicit none
        real(8) :: dlm0
        integer, intent(in) :: ll, m
        real(8), intent(in) :: theta
        
        dlm0=(-1)**m*pm1(ll,m,theta)
        return
    end function dlm0

    !  calculates value of legendre polynomial for l,m,theta
    function pm1(l,m,theta)
        implicit none
        real(8) :: pm1, theta
        real(8) :: zero, one, two, one80, pi, thedeg, x, pp, pm2, rat
        real(8) :: y, al, al2
        integer :: l, m, lmax, imax, ai, i, low
        data pi/3.1415926535897932d0/
        data zero, one, two, one80 /0.d0, 1.d0, 2.d0, 180.d0/
        thedeg=(theta*pi)/one80
        !
        !  if m>l pm1=0 !
        !
        if(m.gt.l) then
        pm1=zero
        return
        end if
        lmax=l
        x=cos(thedeg)
        if (m.ge.0) go to 1
        write (6,100)
        100 format('  NEGATIVE M IN LEGENDRE ROUTINE:  ABORT')
        stop
        !     call exit
        1 if (m.gt.0) go to 5
        !  here for regular legendre polynomials
        pm1=one
        pm2=zero
        do 2 l=1,lmax
        pp=((two*l-one)*x*pm1-(l-one)*pm2)/float(l)
        pm2=pm1
        2 pm1=pp
        return
        !
        !  here for alexander-legendre polynomials
        !
        5 imax=2*m
        rat=1.
        do 6 i=2,imax,2
        ai=i
        6 rat=rat*((ai-one)/ai)
        y=sin(thedeg)
        pm1=sqrt(rat)*(y**m)
        pm2=zero
        low=m+1
        do 10 l=low,lmax
        al=(l+m)*(l-m)
        al=one/al
        al2=((l+m-1)*(l-m-1))*al
        al=sqrt(al)
        al2=sqrt(al2)
        pp=(two*l-one)*x*pm1*al-pm2*al2
        pm2=pm1
        10 pm1=pp
        return
    end function pm1

end module mod_hitestptn