module mod_cosout
   implicit none
   integer, dimension(:), allocatable :: jout
   integer, allocatable               :: nnout
   contains
   subroutine allocate_cosout(n)
      integer, intent(in) :: n
      allocate(jout(n)) ; allocate(nnout)
   end subroutine allocate_cosout
end module mod_cosout

! module mod_coiout
!    implicit none
!    integer, dimension(:), allocatable :: indout
!    integer, allocatable               :: niout
!    contains
!    subroutine allocate_coiout(n)
!       integer, intent(in) :: n
!       allocate(indout(n)) ; allocate(niout)
!    end subroutine allocate_coiout
! end module mod_coiout

! module mod_cov2
!    implicit none
!    real(8), dimension(:), allocatable :: v2
!    integer, allocatable               :: nv2max, ndummy
!    contains
!    subroutine allocate_cov2(n)
!       integer, intent(in) :: n
!       allocate(v2(n)) ; allocate(nv2max) ; allocate(ndummy)
!    end subroutine allocate_cov2
! end module mod_cov2

! module mod_coiv2
!    implicit none
!    real(8), dimension(:), allocatable :: iv2
!    contains
!    subroutine allocate_coiv2(n)
!       integer, intent(in) :: n
!       allocate(iv2(n)) ;
!    end subroutine allocate_coiv2
! end module mod_coiv2

! module mod_cocent
!    implicit none
!    real(8), dimension(:), allocatable :: cent
!    contains
!    subroutine allocate_cocent(n)
!       integer, intent(in) :: n
!       allocate(cent(n)) ;
!    end subroutine allocate_cocent
! end module mod_cocent

! module mod_coeint
!    implicit none
!    real(8), dimension(:), allocatable :: eint
!    contains
!    subroutine allocate_coeint(n)
!       integer, intent(in) :: n
!       allocate(eint(n)) ;
!    end subroutine allocate_coeint
! end module mod_coeint

! module mod_coj12
!    implicit none
!    real(8), dimension(:), allocatable :: j12
!    contains
!    subroutine allocate_coj12(n)
!       integer, intent(in) :: n
!       allocate(j12(n)) ;
!    end subroutine allocate_coj12
! end module mod_coj12

! module mod_coj12pk
!    implicit none
!    real(8), dimension(:), allocatable :: j12pk
!    contains
!    subroutine allocate_coj12pk(n)
!       integer, intent(in) :: n
!       allocate(j12pk(n)) ;
!    end subroutine allocate_coj12pk
! end module mod_coj12pk

! module mod_covvl
!    implicit none
!    real(8), dimension(:), allocatable :: vvl
!    contains
!    subroutine allocate_covvl(n)
!       integer, intent(in) :: n
!       allocate(vvl(n)) ;
!    end subroutine allocate_covvl
! end module mod_covvl

! module mod_cofact
!    implicit none
!    real(8), dimension(:), allocatable :: si
!    contains
!    subroutine allocate_cofact(n)
!       integer, intent(in) :: n
!       allocate(si(n)) ;
!    end subroutine allocate_cofact
! end module mod_cofact

! module mod_coener
!    implicit none
!    real(8), dimension(:), allocatable :: energ
!    contains
!    subroutine allocate_coener(n)
!       integer, intent(in) :: n
!       allocate(energ(n)) ;
!    end subroutine allocate_coener
! end module mod_coener

! module mod_clseg
!    implicit none
!    integer, allocatable :: lseg, intrel, lchar
!    contains
!    subroutine allocate_clseg()
!       allocate(lseg) ; allocate(intrel) ; allocate(lchar) ;
!    end subroutine allocate_clseg
! end module mod_clseg

! module mod_cobuf
!    implicit none
!    integer, dimension(:), allocatable :: ibuf
!    integer, allocatable               :: lbuf
!    contains
!    subroutine allocate_cobuf()
!       allocate(ibuf(1024)) ; allocate(lbuf) ;
!    end subroutine allocate_cobuf
! end module mod_cobuf

! module mod_cofil
!    implicit none
!    integer, dimension(:), allocatable :: maxrec, iofrec
!    integer, allocatable               :: nfl, iofbuf, nwrec
!    contains
!    subroutine allocate_cofil()
!       allocate(maxrec(60)) ; allocate(iofrec(60)) ;
!       allocate(nfl)        ; allocate(iofbuf)      ; allocate(nwrec)
!    end subroutine allocate_cofil
! end module mod_cofil

! module mod_conlam
!    implicit none
!    integer, dimension(:), allocatable :: lamnum
!    integer, allocatable               :: nlam, nlammx
!    contains
!    subroutine allocate_conlam(n)
!       integer, intent(in) :: n
!       allocate(lamnum(n)) ; allocate(nlam) ; allocate(nlammx)
!    end subroutine allocate_conlam
! end module mod_conlam

! module mod_coatpi
!    implicit none
!    integer, dimension(:), allocatable :: isiz
!    integer, allocatable               :: narray
!    contains
!    subroutine allocate_coatpi(n)
!       integer, intent(in) :: n
!       allocate(isiz(n)) ; allocate(narray)
!    end subroutine allocate_coatpi
! end module mod_coatpi

! module mod_coatp3
!    implicit none
!    integer, dimension(:), allocatable :: isizh
!    contains
!    subroutine allocate_coatp3(n)
!       integer, intent(in) :: n
!       allocate(isizh(n))
!    end subroutine allocate_coatp3
! end module mod_coatp3

! module mod_coatpr
!    implicit none
!    real(8), dimension(:), allocatable :: c
!    contains
!    subroutine allocate_coatpr(n)
!       integer, intent(in) :: n
!       allocate(c(n))
!    end subroutine allocate_coatpr
! end module mod_coatpr

! module mod_coatp1
!    implicit none
!    real(8), dimension(:), allocatable :: ctemp
!    contains
!    subroutine allocate_coatp1(n)
!       integer, intent(in) :: n
!       allocate(ctemp(n))
!    end subroutine allocate_coatp1
! end module mod_coatp1

! module mod_coatp2
!    implicit none
!    real(8), dimension(:), allocatable :: chold
!    contains
!    subroutine allocate_coatp2(n)
!       integer, intent(in) :: n
!       allocate(chold(n))
!    end subroutine allocate_coatp2
! end module mod_coatp2

! module mod_coz
!    implicit none
!    real(8), dimension(:,:), allocatable :: z
!    contains
!    subroutine allocate_coz(n)
!       integer, intent(in) :: n
!       allocate(z(n,n))
!    end subroutine allocate_coz
! end module mod_coz

! module mod_cow
!    implicit none
!    real(8), dimension(:,:), allocatable :: w
!    contains
!    subroutine allocate_cow(n)
!       integer, intent(in) :: n
!       allocate(w(n,n))
!    end subroutine allocate_cow
! end module mod_cow

! module mod_cozmat
!    implicit none
!    real(8), dimension(:,:), allocatable :: zmat
!    contains
!    subroutine allocate_cozmat(n)
!       integer, intent(in) :: n
!       allocate(zmat(n,n))
!    end subroutine allocate_cozmat
! end module mod_cozmat

! module mod_coamat
!    implicit none
!    real(8), dimension(:,:), allocatable :: amat
!    contains
!    subroutine allocate_coamat(n)
!       integer, intent(in) :: n
!       allocate(amat(n,n))
!    end subroutine allocate_coamat
! end module mod_coamat

! module mod_cobmat
!    implicit none
!    real(8), dimension(:,:), allocatable :: bmat
!    contains
!    subroutine allocate_cobmat(n)
!       integer, intent(in) :: n
!       allocate(bmat(n,n))
!    end subroutine allocate_cobmat
! end module mod_cobmat

! module mod_cotq1
!    implicit none
!    real(8), dimension(:,:), allocatable :: tq1
!    contains
!    subroutine allocate_cotq1(n)
!       integer, intent(in) :: n
!       allocate(tq1(n,n))
!    end subroutine allocate_cotq1
! end module mod_cotq1

! module mod_cotq2
!    implicit none
!    real(8), dimension(:,:), allocatable :: tq2
!    contains
!    subroutine allocate_cotq2(n)
!       integer, intent(in) :: n
!       allocate(tq2(n,n))
!    end subroutine allocate_cotq2
! end module mod_cotq2

! module mod_cotq3
!    implicit none
!    real(8), dimension(:,:), allocatable :: tq3
!    contains
!    subroutine allocate_cotq3(n)
!       integer, intent(in) :: n
!       allocate(tq3(n,n))
!    end subroutine allocate_cotq3
! end module mod_cotq3

! module mod_cojq
!    implicit none
!    real(8), dimension(:), allocatable :: jq
!    contains
!    subroutine allocate_cojq(n)
!       integer, intent(in) :: n
!       allocate(jq(n))
!    end subroutine allocate_cojq
! end module mod_cojq

! module mod_colq
!    implicit none
!    real(8), dimension(:), allocatable :: lq
!    contains
!    subroutine allocate_colq(n)
!       integer, intent(in) :: n
!       allocate(lq(n))
!    end subroutine allocate_colq
! end module mod_colq

! module mod_coinq
!    implicit none
!    integer, dimension(:), allocatable :: inq
!    contains
!    subroutine allocate_coinq(n)
!       integer, intent(in) :: n
!       allocate(inq(n))
!    end subroutine allocate_coinq
! end module mod_coinq

! module mod_cojhld
!    implicit none
!    integer, dimension(:), allocatable :: jhold
!    contains
!    subroutine allocate_cojhld(n)
!       integer, intent(in) :: n
!       allocate(jhold(n))
!    end subroutine allocate_cojhld
! end module mod_cojhld

! module mod_coehld
!    implicit none
!    real(8), dimension(:), allocatable :: ehold
!    contains
!    subroutine allocate_coehld(n)
!       integer, intent(in) :: n
!       allocate(ehold(n))
!    end subroutine allocate_coehld
! end module mod_coehld

! module mod_coinhl
!    implicit none
!    integer, dimension(:), allocatable :: inhold
!    contains
!    subroutine allocate_coinhl(n)
!       integer, intent(in) :: n
!       allocate(inhold(n))
!    end subroutine allocate_coinhl
! end module mod_coinhl

! module mod_coisc1
!    implicit none
!    integer, dimension(:), allocatable :: isc1
!    contains
!    subroutine allocate_coisc1(n)
!       integer, intent(in) :: n
!       allocate(isc1(n))
!    end subroutine allocate_coisc1
! end module mod_coisc1

! module mod_coisc2
!    implicit none
!    integer, dimension(:), allocatable :: isc2
!    contains
!    subroutine allocate_coisc2(n)
!       integer, intent(in) :: n
!       allocate(isc2(n))
!    end subroutine allocate_coisc2
! end module mod_coisc2

! module mod_coisc3
!    implicit none
!    integer, dimension(:), allocatable :: isc3
!    contains
!    subroutine allocate_coisc3(n)
!       integer, intent(in) :: n
!       allocate(isc3(n))
!    end subroutine allocate_coisc3
! end module mod_coisc3

! module mod_coisc4
!    implicit none
!    integer, dimension(:), allocatable :: isc4
!    contains
!    subroutine allocate_coisc4(n)
!       integer, intent(in) :: n
!       allocate(isc4(n))
!    end subroutine allocate_coisc4
! end module mod_coisc4

! module mod_coisc5
!    implicit none
!    integer, dimension(:), allocatable :: isc5
!    contains
!    subroutine allocate_coisc5(n)
!       integer, intent(in) :: n
!       allocate(isc5(n))
!    end subroutine allocate_coisc5
! end module mod_coisc5

! module mod_coisc6
!    implicit none
!    integer, dimension(:), allocatable :: isc6
!    contains
!    subroutine allocate_coisc6(n)
!       integer, intent(in) :: n
!       allocate(isc6(n))
!    end subroutine allocate_coisc6
! end module mod_coisc6

! module mod_coisc7
!    implicit none
!    integer, dimension(:), allocatable :: isc7
!    contains
!    subroutine allocate_coisc7(n)
!       integer, intent(in) :: n
!       allocate(isc7(n))
!    end subroutine allocate_coisc7
! end module mod_coisc7

! module mod_coisc8
!    implicit none
!    integer, dimension(:), allocatable :: isc8
!    contains
!    subroutine allocate_coisc8(n)
!       integer, intent(in) :: n
!       allocate(isc8(n))
!    end subroutine allocate_coisc8
! end module mod_coisc8

! module mod_coisc9
!    implicit none
!    integer, dimension(:), allocatable :: isc9
!    contains
!    subroutine allocate_coisc9(n)
!       integer, intent(in) :: n
!       allocate(isc9(n))
!    end subroutine allocate_coisc9
! end module mod_coisc9

! module mod_coisc10
!    implicit none
!    integer, dimension(:), allocatable :: isc10
!    contains
!    subroutine allocate_coisc10(n)
!       integer, intent(in) :: n
!       allocate(isc10(n))
!    end subroutine allocate_coisc10
! end module mod_coisc10

! module mod_coisc11
!    implicit none
!    integer, dimension(:), allocatable :: isc11
!    contains
!    subroutine allocate_coisc11(n)
!       integer, intent(in) :: n
!       allocate(isc11(n))
!    end subroutine allocate_coisc11
! end module mod_coisc11

! module mod_coisc12
!    implicit none
!    integer, dimension(:), allocatable :: isc12
!    contains
!    subroutine allocate_coisc12(n)
!       integer, intent(in) :: n
!       allocate(isc12(n))
!    end subroutine allocate_coisc12
! end module mod_coisc12

! module mod_colsc1
!    implicit none
!    integer, dimension(:), allocatable :: lsc1
!    contains
!    subroutine allocate_colsc1(n)
!       integer, intent(in) :: n
!       allocate(lsc1(n))
!    end subroutine allocate_colsc1
! end module mod_colsc1

! module mod_cosc1
!    implicit none
!    real(8), dimension(:), allocatable :: sc1
!    contains
!    subroutine allocate_cosc1(n)
!       integer, intent(in) :: n
!       allocate(sc1(n))
!    end subroutine allocate_cosc1
! end module mod_cosc1

! module mod_cosc2
!    implicit none
!    real(8), dimension(:), allocatable :: sc2
!    contains
!    subroutine allocate_cosc2(n)
!       integer, intent(in) :: n
!       allocate(sc2(n))
!    end subroutine allocate_cosc2
! end module mod_cosc2

! module mod_cosc3
!    implicit none
!    real(8), dimension(:), allocatable :: sc3
!    contains
!    subroutine allocate_cosc3(n)
!       integer, intent(in) :: n
!       allocate(sc3(n))
!    end subroutine allocate_cosc3
! end module mod_cosc3

! module mod_cosc4
!    implicit none
!    real(8), dimension(:), allocatable :: sc4
!    contains
!    subroutine allocate_cosc4(n)
!       integer, intent(in) :: n
!       allocate(sc4(n))
!    end subroutine allocate_cosc4
! end module mod_cosc4

! module mod_cosc5
!    implicit none
!    real(8), dimension(:), allocatable :: sc5
!    contains
!    subroutine allocate_cosc5(n)
!       integer, intent(in) :: n
!       allocate(sc5(n))
!    end subroutine allocate_cosc5
! end module mod_cosc5

! module mod_cosc6
!    implicit none
!    real(8), dimension(:), allocatable :: sc6
!    contains
!    subroutine allocate_cosc6(n)
!       integer, intent(in) :: n
!       allocate(sc6(n))
!    end subroutine allocate_cosc6
! end module mod_cosc6

! module mod_cosc7
!    implicit none
!    real(8), dimension(:), allocatable :: sc7
!    contains
!    subroutine allocate_cosc7(n)
!       integer, intent(in) :: n
!       allocate(sc7(n))
!    end subroutine allocate_cosc7
! end module mod_cosc7

! module mod_cosc8
!    implicit none
!    real(8), dimension(:), allocatable :: sc8
!    contains
!    subroutine allocate_cosc8(n)
!       integer, intent(in) :: n
!       allocate(sc8(n))
!    end subroutine allocate_cosc8
! end module mod_cosc8

! module mod_cosc9
!    implicit none
!    real(8), dimension(:), allocatable :: sc9
!    contains
!    subroutine allocate_cosc9(n)
!       integer, intent(in) :: n
!       allocate(sc9(n))
!    end subroutine allocate_cosc9
! end module mod_cosc9

! module mod_cosc10
!    implicit none
!    real(8), dimension(:), allocatable :: sc10
!    contains
!    subroutine allocate_cosc10(n)
!       integer, intent(in) :: n
!       allocate(sc10(n))
!    end subroutine allocate_cosc10
! end module mod_cosc10

! module mod_cosc11
!    implicit none
!    real(8), dimension(:), allocatable :: sc11
!    contains
!    subroutine allocate_cosc11(n)
!       integer, intent(in) :: n
!       allocate(sc11(n))
!    end subroutine allocate_cosc11
! end module mod_cosc11

module mod_coeig
   ! module for ba3d1p basis routine
   implicit none
   real(8), dimension(:,:), allocatable :: c0, c1, c2
   contains
   subroutine allocate_coeig()
      allocate(c0(4,4)) ; allocate(c1(3,3)) ; allocate(c2(2,2))
   end subroutine allocate_coeig
end module mod_coeig

module mod_coeig2
   ! * module for ba3p2s basis routine
   ! matrices for transformation between atomic and molecular BF functions
   implicit none
   real(8), dimension(:,:), allocatable :: t12, t32
   contains
   subroutine allocate_coeig2()
      allocate(t12(5,5)) ; allocate(t32(3,3))
   end subroutine allocate_coeig2
end module mod_coeig2

module mod_cokaux
   implicit none
   integer, allocatable :: naux
   contains
   subroutine allocate_cokaux(anaux)
      integer, intent(in) :: anaux
      allocate(naux)
      naux = anaux
   end subroutine allocate_cokaux
end module mod_cokaux

module mod_cotble
   ! *     npnt:     max. number of pointer
   ! *     jttble:   array containing pointer to records in s-matrix   
   implicit none
   integer, dimension(:), allocatable :: jttble
   integer, allocatable               :: npnt
   contains
   subroutine allocate_cotble(afact)
      integer, intent(in) :: afact
      allocate(jttble(afact)) ; allocate(npnt)
      npnt = afact
   end subroutine allocate_cotble
end module mod_cotble


module mod_coqvec
   ! *  vectors dimensioned nphoto*nch
   !       q

   ! *  variables in common block /coqvec/
   ! *     mxphot        maximum column dimension of q matrix
   ! *     nphoto        actual collumn dimension of q matrix
   ! *     q             accumulated gamma2 inhomogeneous propagator
   ! *                   only calculated if photof = .true.
   ! *                   this is stored as a column vector for each separate
   ! *                   ground-state wavefunction
   
   ! *     variables in common block /coqvec/
   ! *     mxphot        maximum dimension of q matrix
   ! *     nphoto        actual collumn dimension of q matrix
   ! *     q             accumulated overlap matrix with ground state
   ! *                   only calculated if photof = .true.
   ! *                   this is stored with each wavefunction as a column vector

   ! *     variables in common block /coqvec/
   ! *     mxphot        maximum column dimension of q matrix
   ! *     nphoto        actual collumn dimension of q matrix
   ! *     q             accumulated overlap matrix with ground state
   ! *                   only calculated if photof = .true.
   ! *                   this is stored with each wavefunction as a column vector
   
   ! *     variables in common block /coqvec/
   ! *     mxphot        maximum column dimension of q matrix
   ! *     nphoto        actual column dimension of q matrix (gamma2)
   ! *     q             accumulated overlap matrix with ground state
   ! *                   only calculated if photof = .true.
   ! *                   this is stored with each wavefunction as a column vector

   ! *     variables in common block /coqvec/
   ! *     mxphot        maximum column dimension of q matrix
   ! *     nphoto        actual column dimension of q matrix
   ! *     q             accumulated overlap matrix with ground state
   ! *                   only calculated if photof = .true.
   ! *                   this is stored with each wavefunction as a column vector

   ! *     variables in common block /coqvec/
   ! *     mxphot        maximum dimension of q matrix
   ! *     nphoto        actual collumn dimension of q matrix
   ! *     q             accumulated overlap matrix with ground state
   ! *                   only calculated if photof = .true.
   ! *                   this is stored with each wavefunction as a column vector


   ! mxphot:   maximum number of initial states (maximum columns in q)
   ! nphoto:   actual number of initial states used
   ! q:        accumulated overlap matrix with ground state
   !           only calculated if photof = .true.
   implicit none
   real(8), dimension(:), allocatable :: q
   integer, allocatable               :: mxphot
   integer, allocatable               :: nphoto
   contains
   subroutine allocate_coqvec(aqmax, amxpho, anphot)
      integer, intent(in) :: aqmax
      integer, intent(in) :: amxpho
      integer, intent(in) :: anphot
      allocate(q(aqmax)) ; allocate(mxphot) ; allocate(nphoto)
      mxphot = amxpho
      nphoto = anphot
   end subroutine allocate_coqvec
end module mod_coqvec

module mod_coqvec2
   implicit none
   real(8), dimension(:), allocatable :: q2
   contains
   subroutine allocate_coqvec2(aq2)
      integer, intent(in) :: aq2
      allocate(q2(aq2))
   end subroutine allocate_coqvec2
end module mod_coqvec2

module mod_codim
   implicit none
   save
   integer, allocatable :: mairy, mmax, mbig
   contains
   subroutine allocate_codim(aairy, amax, abig)
      integer, intent(in) :: aairy
      integer, intent(in) :: amax
      integer, intent(in) :: abig
      allocate(mairy) ; allocate(mmax) ; allocate(mbig)
      mairy = aairy
      mmax = amax
      mbig = abig
   end subroutine allocate_codim
end module mod_codim

module mod_comxbs
   implicit none
   save
   integer, allocatable :: maxbas  ! maximum number of allowed basis types 
   ! THIS SHOULD BE INCREASED AS BASIS ROUTINES ARE ADDED
   contains
   subroutine allocate_comxbs(n)
      integer, intent(in) :: n
      allocate(maxbas) ; maxbas = n
   end subroutine allocate_comxbs
end module mod_comxbs

module mod_comxm
   implicit none
   save
   integer, allocatable :: ncache, mxmblk
   contains
   subroutine allocate_comxm()
      allocate(ncache)  ! cache size in words
      allocate(mxmblk)
      !  determine cache and block sizes for matrix multiply
      ncache=4096
#if defined(HIB_UNIX_HP)
      ncache=16000
#endif
      mxmblk=64
   end subroutine allocate_comxm
end module mod_comxm




 ! all the commons blocks from himain.t:
    !!   common /comom/  xmom(3), imom(13)
    !!   common /cosout/ nnout, jout(kout)
    !!   common /coiout/ niout, indout(kout)
    !!   common /cov2/ nv2max, ndummy, v2(kv2max)
    !!   common /coiv2/ iv2(kv2max)
    !!   common /cocent/ cent(kmax)
    !!   common /coeint/ eint(kmax)
    !!   common /coj12/ j12(kmax)
    !!   common /coj12p/ j12pk(kmax)
    !!   common /covvl/  vvl(klammx)
    !!   common /cofact/ si(kfact)
    !!  common /coener/ energ(ken)
    !!   common /clseg/  lseg,intrel,lchar
    !!   common /cobuf/  lbuf,ibuf(1024)
    !!   common /cofil/  nfl,iofbuf,maxrec(60),iofrec(60),nwrec
    !!   common /conlam/ nlam, nlammx, lamnum(klammx)
    !!   common /coatpi/ narray, isiz(krotmx)
    !!   common /coatp3/ isizh(krotmx)
    !!   common /coatpr/ c(krotmx)
    !!   common /coatp1/ ctemp(krotmx)
    !!   common /coatp2/ chold(krotmx)
    !!   common /coz/ z(kmax,kmax)
    !!   common /cow/ w(kmax,kmax)
    !!   common /cozmat/ zmat(kmax,kmax)
    !!   common /coamat/ amat(kmax,kmax)
    !!   common /cobmat/ bmat(kairy,kairy)
    !!   common /cotq1/ tq1(kmax,kmax)
    !!   common /cotq2/ tq2(kmax,kmax)
    !!   common /cotq3/ tq3(kmax,kmax)
    !!   common /cojq/ jq(kmax)
    !!   common /colq/ lq(kmax)
    !!   common /coinq/ inq(kmax)
    !!   common /cojhld/ jhold(kmax)
    !!   common /coehld/ ehold(kmax)
    !!   common /coinhl/ inhold(kmax)
    !!   common /coisc1/ isc1(kmax)
    !!   common /coisc2/ isc2(kmax)
    !!   common /coisc3/ isc3(kmax)
    !!   common /coisc4/ isc4(kmax)
    !!   common /coisc5/ isc5(kmax)
    !!   common /coisc6/ isc6(kmax)
    !!   common /coisc7/ isc7(kmax)
    !!   common /coisc8/ isc8(kmax)
    !!   common /coisc9/ isc9(kmax)
    !!   common /coisc10/ isc10(kmax)
    !!   common /coisc11/ isc11(kmax)
    !!   common /coisc12/ isc12(kmax)
    !!   common /colsc1/ lsc1(kmax)
    !!   common /cosc1/ sc1(kmax)
    !!   common /cosc2/ sc2(kmax)
    !!   common /cosc3/ sc3(kmax)
    !!   common /cosc4/ sc4(kmax)
    !!   common /cosc5/ sc5(kmax)
    !!   common /cosc6/ sc6(kmax)
    !!   common /cosc7/ sc7(kmax)
    !!   common /cosc8/ sc8(kmax)
    !!   common /cosc9/ sc9(kmax)
    !!   common /cosc10/ sc10(kmax)
    !!   common /coeig2/  t12(5,5), t32(3,3)
    !!   common /coeig/  c0(4,4), c1(3,3), c2(2,2)
    !!   common /cosc11/ sc11(kaux3)
    !!   common /cosc11/ sc11(kaux)
    !!   common /cokaux/ naux
    !!   common /cotble/ npnt, jttble(kfact)
    !!   common /coqvec/ mxphot, nphoto, q(kqmax)
    !!   common /coqvec2/ q2(kq2)
    !!   common /codim/ mairy,mmax,mbig
    !!   common /comxbs/ maxbas
    !!   common /comxm/ ncache, mxmblk