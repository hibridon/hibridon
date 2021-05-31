module mod_comom
   implicit none
   real(8), dimension(:), allocatable :: xmom
   integer, dimension(:), allocatable :: imom
   contains
   subroutine allocate_comom() 
      allocate(xmom(3)) ; allocate(imom(13))   
   end subroutine allocate_comom
end module mod_comom

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

module mod_coiout
   implicit none
   integer, dimension(:), allocatable :: indout
   integer, allocatable               :: niout
   contains
   subroutine allocate_coiout(n) 
      integer, intent(in) :: n
      allocate(indout(n)) ; allocate(niout)
   end subroutine allocate_coiout
end module mod_coiout

module mod_cov2
   implicit none
   real(8), dimension(:), allocatable :: v2
   integer, allocatable               :: nv2max, ndummy
   contains
   subroutine allocate_cov2(n) 
      integer, intent(in) :: n
      allocate(v2(n)) ; allocate(nv2max) ; allocate(ndummy)
   end subroutine allocate_cov2
end module mod_cov2

 
 ! All the commons blocks from himain.t:
    ! common /comom/  xmom(3), imom(13)
    !   common /cosout/ nnout, jout(kout)
    !   common /coiout/ niout, indout(kout)
    !   common /cov2/ nv2max, ndummy, v2(kv2max)
    !   common /coiv2/ iv2(kv2max)
    !   common /cocent/ cent(kmax)
    !   common /coeint/ eint(kmax)
    !   common /coj12/ j12(kmax)
    !   common /coj12p/ j12pk(kmax)
    !   common /covvl/  vvl(klammx)
    !   common /cofact/ si(kfact)
    !   common /coener/ energ(ken)
    !   common /clseg/  lseg,intrel,lchar
    !   common /cobuf/  lbuf,ibuf(1024)
    !   common /cofil/  nfl,iofbuf,maxrec(60),iofrec(60),nwrec
    !   common /conlam/ nlam, nlammx, lamnum(klammx)
    !   common /coatpi/ narray, isiz(krotmx)
    !   common /coatp3/ isizh(krotmx)
    !   common /coatpr/ c(krotmx)
    !   common /coatp1/ ctemp(krotmx)
    !   common /coatp2/ chold(krotmx)
    !   common /coz/ z(kmax,kmax)
    !   common /cow/ w(kmax,kmax)
    !   common /cozmat/ zmat(kmax,kmax)
    !   common /coamat/ amat(kmax,kmax)
    !   common /cobmat/ bmat(kairy,kairy)
    !   common /cotq1/ tq1(kmax,kmax)
    !   common /cotq2/ tq2(kmax,kmax)
    !   common /cotq3/ tq3(kmax,kmax)
    !   common /cojq/ jq(kmax)
    !   common /colq/ lq(kmax)
    !   common /coinq/ inq(kmax)
    !   common /cojhld/ jhold(kmax)
    !   common /coehld/ ehold(kmax)
    !   common /coinhl/ inhold(kmax)
    !   common /coisc1/ isc1(kmax)
    !   common /coisc2/ isc2(kmax)
    !   common /coisc3/ isc3(kmax)
    !   common /coisc4/ isc4(kmax)
    !   common /coisc5/ isc5(kmax)
    !   common /coisc6/ isc6(kmax)
    !   common /coisc7/ isc7(kmax)
    !   common /coisc8/ isc8(kmax)
    !   common /coisc9/ isc9(kmax)
    !   common /coisc10/ isc10(kmax)
    !   common /coisc11/ isc11(kmax)
    !   common /coisc12/ isc12(kmax)
    !   common /colsc1/ lsc1(kmax)
    !   common /cosc1/ sc1(kmax)
    !   common /cosc2/ sc2(kmax)
    !   common /cosc3/ sc3(kmax)
    !   common /cosc4/ sc4(kmax)
    !   common /cosc5/ sc5(kmax)
    !   common /cosc6/ sc6(kmax)
    !   common /cosc7/ sc7(kmax)
    !   common /cosc8/ sc8(kmax)
    !   common /cosc9/ sc9(kmax)
    !   common /cosc10/ sc10(kmax)
    !   common /coeig2/  t12(5,5), t32(3,3)
    !   common /coeig/  c0(4,4), c1(3,3), c2(2,2)
    !   common /cosc11/ sc11(kaux3)
    !   common /cosc11/ sc11(kaux)
    !   common /cokaux/ naux
    !   common /cotble/ npnt, jttble(kfact)
    !   common /coqvec/ mxphot, nphoto, q(kqmax)
    !   common /coqvec2/ q2(kq2)
    !   common /codim/ mairy,mmax,mbig
    !   common /comxbs/ maxbas
    !   common /comxm/ ncache, mxmblk