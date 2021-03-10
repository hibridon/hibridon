*                         hibound  library                           *
*                                                                       *
*************************************************************************
*                          routines included:                           *
*                                                                       *
*  1. bound      millard alexander's bound state program                 *
*  2. bound_wavfn
*  3. h_basis

*-------------------------------------------------------------
       subroutine bound(s,t,v,wtemp,scmat,nch,nmax)
*      This routine is called from 'propag in hibrid3'
*      and calling subroutine 'boundwavfn'
*      author: susan k gregurick 3-aug-92
*      modified by moonbong yang 27-Sep-94
*      completely rewritten by mha 24-apr-1997
*      revised by mha 11-may-1997
*      latest revision 23-feb-2004
* definition of variables in call list:
*      nch     number of coupled equations
*      nmax    leading dimension of matrices w
*-----------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical lpar, batch, lparx, wavefl
      character*40 jobnam, input, output, savfil,
     :             evalfil
      character*1 jobz
      integer vmax,nvmax
      include "common/parpot"
*   square matrices
      dimension s(nmax,6)
      dimension t(nmax,6)
      dimension v(nmax,6)
      dimension wtemp(nmax,6)
      dimension scmat(nmax,6)
cstart unix-darwin
* scratch vectors for dsyev
      dimension work(57*nmax),iwork(10*nmax),isuppz(2*nmax)
cend

      common /cotq1/ wr(1)
      common /cotq3/ z(1)
* scratch vectors
      common /cosc1/ sc1(1)
      common /cosc2/ sc2(1)
      common /cosc3/ sc3(15)
      common /cosc4/ sc4(15)
      common /cosc5/ sc5(1)
      common /cosc6/ eval(15)
      common /cosc11/ sc11(1)
      common /coiout/ niout, indout(1)
      common /coisc1/ isc1(1)
      common /corpar/ r1,r2,c,spac,delr,hsimp,eigmin,tolai,xmu
      common /cofile/ input, output, jobnam, savfil

      common /colpar/ lpar(3), batch, lparx(22),wavefl
      common /coconv/ econv, xmconv
      common /coered/ ered, rmu
      data pi,zero,one,half,two /3.14159265358979d0,0d0,1d0,0.5d0,2d0/
      data ione,izero /1,0/
** Boh2; r1=3, r2=20, hsimp=0.025, spac=0.4
** maxch=27 vmax=42 for Jtot=2 B-oH2
** Arhcl: r1=3, r2=20, hsimp=0.025, spac=0.3
** BHAr: r1=4, r2=20, hsimp=0.025, spac=0.4
** ArC2H2_new: r1=5, r2=20, hsimp=0.025, spac=0.4
      naux=3*nmax
*      r1=4.0d0
*      r2=15.d0
*      r2=8.4d0
*       r2=4.4d0
*      rmu=3.105695923666346d+03
* hsimp is integration step size
*      hsimp=0.025d0
* spacing between centers of gaussians
*      spac=0.4d0
*      c=.5d0

      vmax=(r2-r1)/spac+1.d-8
      del=(r2-r1)/vmax
      alph=(c/del)**2
      vmax=vmax+1
      nvmax=nch*vmax
*  open file for eigenvalues
      call gennam (evalfil, jobnam, ienerg, 'evl', lenx)
      call openf(1,evalfil,'sf',0)
      call version(1)
      write (1,1) potnam, label
1     format (' ** POTNAM: ',(a),/,' ** LABEL: ',(a))
      write (6, 5) r1, r2,  spac, del, vmax, c, alph
      write (9, 5) r1, r2,  spac, del, vmax, c, alph
      write (1, 5) r1, r2,  spac, del, vmax, c, alph
5     format (/,' ** BOUND STATE CALCULATION:  R1 =',
     :   f5.2,';  R2 =', f5.2,/,
     :  6x,'DEL-START =',f5.3,', DEL-EXACT =',f8.6,
     :  ', NO. GAUSSIANS =',i3,/,
     :  6x,'C =',f4.2,', ALPH =',g10.5)
      if (nvmax .gt. nmax) then
        write (1, 10) nch, vmax, nmax
        write (9, 10) nch, vmax, nmax
        write (6, 10) nch, vmax, nmax
10      format (' ** PRODUCT OF NCH (',i3,') AND VMAX (',i3,')',
     :    ' .GT. NMAX (',i4,') IN BOUND; ABORT **')
        if (batch) call exit
        return
      endif
* determine vectors containing (ri-rj)^2
* and ri's
      do i=1,vmax
        tt=(i-1)*del
        sc1(i)=tt*tt
        sc2(i)=exp(-half*alph*sc1(i))
        sc3(i)=tt+r1
      enddo
** determining S and T matrices
*      print *, '** determining S and T matrices ...'
      do i=1,vmax
* here for ith diagonal, selt and telt are overlap and kinetic energy
* elements for this diagonal
          selt=sc2(i)
          telt=sc2(i)*alph*half*(one-alph*sc1(i))/rmu
* ndiag is number of elements on ith diagonal
          ndiag=vmax-i+1
          call dset(ndiag,selt,s(i,1),nmax+1)
          call dset(ndiag,selt,s(1,i),nmax+1)
          call dset(ndiag,telt,t(i,1),nmax+1)
          call dset(ndiag,telt,t(1,i),nmax+1)
      enddo
* copy S matrix to scmat
      call matmov(s,scmat,vmax,vmax,nmax,nmax)
** determining eigenvectors of S matrix (returned in sc4)

*      print *, '** determining eigenvalues of S matrix ...'

cstart unix .and. .not. unix-darwin
c;      call rs (nmax, vmax, scmat, sc4, ione, scmat, sc11, sc5, ierr)
cend
cstart unix-darwin
      lwork=57*vmax
      liwork=10*vmax
      abstol=0d0
      lsup=2*vmax

      call dsyevr('N','A','L',vmax,scmat,nmax,vl,vu,il,iu,abstol,mm,
     :   sc4,vecnow,nmax,isuppz,work,lwork,iwork,liwork,ierr)
cend

      evalmin=sc4(idamin(vmax,sc4,1))
      write (1, 19) evalmin
      write (6, 19) evalmin
      write (9, 19) evalmin
19    format (' ** MINIMUM EIGENVALUE OF S MATRIX =',1pe11.4)
      if (evalmin .lt. eigmin) then
        write (6,25) eigmin
        write (9,25) eigmin
25      format ('    THIS IS .LT.',1pe8.1,
     : '; INCREASE PARAMETER C OR DECREASE EIGMIN')
        return
      endif
* move s matrix into larger (nch x vmax) matrix
* each vmax x vmax block of this matrix includes all n channels
* and one distributed gaussian function
* set the overlap and hamiltonian matrices to zero
*      print *,
*     :  '** determining overlap and H matrices in full basis ...'
      ndim=vmax*nch
      do icol=1,ndim
        call dset(ndim,zero,scmat(1,icol),1)
        call dset(ndim,zero,wtemp(1,icol),1)
      enddo
* expand the overlap and kinetic energy matrices into scmat and wmat
      do iv=1,vmax
        irow=(iv-1)*nch+1
        do jv=1,vmax
** determining position of (1,1) element in this block with
* respect to full, expanded matrix
          icol=(jv-1)*nch+1
         call dset(nch,s(iv,jv),scmat(irow,icol),nmax+1)
         call dset(nch,t(iv,jv),wtemp(irow,icol),nmax+1)
        enddo
      enddo
*      write (6, 28)
*      write (9, 28)
*28    format ( ' ** integrating V(R) over R ...  ')
* integrate over r
      nr=(r2-r1+two*delr)/hsimp
      if (nr .lt. 2) nr=2
      mr=nr/2
      if (2*mr .ne. nr) nr=nr-1
      hexac=(r2-r1+two*delr)/nr
      write (1, 30) delr,hsimp, nr+1, hexac
      write (6, 30) delr,hsimp, nr+1, hexac
      write (9, 30) delr,hsimp, nr+1, hexac
30    format
     : (' ** SIMPSON''S INTEGRATION OF POTENTIAL',
     :  /,6x,'DEL-R =',f4.2,', H-START =',f5.3,
     :  ', N-POINTS =',i4,', H-EXACT =',f8.6)
      iu=-1
      do ir=1,nr+1
        rnow=(ir-1)*hexac+r1-delr
* determine gaussian basis functions at this point
        call gbasis(vmax,alph,rnow,sc3,sc4)

        call vmat(nmax,ndim,nch,vmax,v,wr,z,sc4,rnow)
* accumulate simpson's rule sum
        wt=(3+iu)*hexac/3.d0
        if (ir.eq.1) wt=half*wt
        if (ir.eq.nr) wt=half*wt
        do icol=1,ndim
          call daxpy(ndim,wt,v(1,icol),1,wtemp(1,icol),1)
        enddo
        iu=-iu
      enddo

** determining eigenvalues and (optionally) eigenvectors
*        do ip=1,5
*          write (6,40) (wtemp(ip,jk), jk=1,5)
*40        format (5(1pg12.4))
*        enddo
*        call exit
      iflag=izero
      if (wavefl) iflag=ione

cstart unix-ibm
c;      call dsygv(iflag,wtemp,nmax,scmat,nmax,eval,t,nmax,
c;     :           ndim,sc11,naux)
cend
cstart unix-darwin
       if (iflag.eq.1) then
           jobz='v'
       else
           jobz='n'
       endif
       lwork=3*nmax
       call dsygv(1,jobz,'l',ndim,wtemp,nmax,scmat,nmax,eval,
     :            sc11,lwork,ierr)
       if (ierr.ne.0) then
         stop 'dsygv in hibound'
       endif
* copy eigenvectors to matrix t
       if (iflag.eq.1) then
          do i=1,ndim
             call dcopy(ndim,wtemp(1,i),1,t(1,i),1)
          enddo
       endif
cend
cstart unix mac .and. .not. unix-ibm .and. .not. unix-darwin
c;      call rsg(nmax,ndim,wtemp,scmat,eval,iflag,t,sc1,sc11,ierr)
cend
*  print job information
      evalmin=eval(idmin(ndim,eval,1))
      if (evalmin .gt. zero) then
        write (1, 50) evalmin*econv
        write (6, 50) evalmin*econv
        write (9, 50) evalmin*econv
50      format
     :  (' ** MINIMUM EIGENVALUE =',f11.3,' .GT. 0; RETURN')
      else
        if (.not. wavefl) then
          write (6,55) evalfil, label
55        format(' ** ALL EIGENVALUES SAVED IN FILE ',a,
     :      /,6x,'LABEL:',a,/,
     :  ' ** NEGATIVE EIGENVALUES ARE (CM-1):')
        else
          write (6,56) evalfil, label
56        format
     :  (' ** ALL EIGENVALUES AND FIRST 5 EIGENVECTORS SAVED IN FILE ',a,
     :      /,6x,'LABEL:',a,/,
     :  ' ** NEGATIVE EIGENVALUES ARE (CM-1):')
        endif
        write (1,57)
57      format(' ** EIGENVALUES (CM-1)')
        nbound=0
        do k=1,ndim
            write (1,65)  eval(k)*econv
          if (eval(k) .le. 0) then
            nbound=nbound+1
            write (9,65)  eval(k)*econv
            write (6,65)  eval(k)*econv
          end if
        end do
65      format(2x,f11.3)
      endif
      if (wavefl) then
        if (niout .lt. 0 .or. niout .gt. nch) then
          write (6, 70) niout, nch
70        format (' ** NIOUT =', i3, '.LE. 0 .OR. .GT. NCH =',
     :     i3,'; ABORT')
          close(1)
          return
        endif
        write (1, 75) nch, vmax, niout
75      format (' NCH =',i3, '; N-GAUSSIANS = ',i3,'; NSTATE =',i3)
        write (1,80) alph, (sc3(i), i=1,vmax)
80      format (/,' GAUSSIAN FUNCTIONS - ALPH =',1pg12.5,
     :   '; CENTERED AT:',/,0p40f8.4)
        write (1, 85)
85      format (' WAVEFUNCTION EXPANSION COEFFICIENTS:',
     :   ' ROWS ARE CHANNELS, COLUMNS ARE GAUSSIANS')
        if (nbound .le. 0) then
          write (6, 90)
          write (9, 90)
          write (1, 90)
90        format (' ** NO BOUND STATES (E < 0)')
          close(1)
          return
        else
          nstate=min(5,nbound)
          do in = 1, nstate
            write (1, 95) in
95          format (' STATE =',i2)
            do i=1,nch
              write (1, 125)
     :          (t(irow,in), irow=i,(vmax-1)*nch+i,nch)
125           format (40(1pg13.5))
            enddo
          enddo
        endif
      endif
      close(1)
      return
      end
* ---------------------------------
      subroutine vmat(nmax,ndim,nch,vmax,v,wr,z,sc4,r)
      implicit double precision (a-h,o-z)
      integer vmax
      common /coered/ ered, rmu

      dimension v(nmax,30),wr(30),z(nmax,30),sc4(vmax)
      xirmu=0.5d0/rmu
* insert nmax in potmat which is compatible with basis
      ntop=nch
      if (mod(ntop,2) .eq. 0) ntop = ntop + 1

      call potmat(wr,r,nch,ntop)
* move from lower triangular to full matrix
      if (nch .ge. 2) then
        do i=1,nch
          ndiag=nch-i+1
* first multiply by 1/(2*rmu)
          call dscal(ndiag,xirmu,wr(i),ntop+1)
* copy to upper triangle, except if diagonal element
          if (i .ne. 1) then
* ncol is index of (1,i) in packed column form
* with max row dimension of ntop
            ncol=(i-1)*ntop+1
            call dcopy(ndiag,wr(i),ntop+1,wr(ncol),ntop+1)
          endif
        enddo
      else
        wr(1)=wr(1)*xirmu
      endif
* also, add ered back to diagonal elements
      do i=1,nch
        ncol=(i-1)*ntop+i
        wr(ncol)=wr(ncol)+ered
      enddo
* expand into full matrix
      do iv=1,vmax
        do jv=1,vmax
* determine position of (1,1) element in this block with
* respect to full, expanded matrix
          irow=(iv-1)*nch+1
          icol=(jv-1)*nch+1
          do ich=1,nch
* ncol is index of (1,i) in packed column form
* with max row dimension of ntop
            ncol=(ich-1)*ntop+1
            call vsmul(wr(ncol),1,sc4(iv)*sc4(jv),z(1,ich),1,nch)
            call dcopy(nch,z(1,ich),1,v(irow,icol+ich-1),1)
          enddo
        enddo
      enddo
      return
      end
* ---------------------------------------------
      subroutine gbasis(vmax,alph,rnow,sc3,sc4)
* determined vector of gaussian basis functions with exponent
* alpha, centered at points sc3, evaluated at rnow
* these are returned in vector sc4

* e.g. sc4(j)=[(2*alph/pi)**0.25]*exp[-alph*(rnow-sc4(j))**2]
* ---------------------------------------------
      implicit double precision (a-h,o-z)
      integer vmax
      dimension sc3(vmax), sc4(vmax)
      data pi,zero,one,half,two /3.14159265358979d0,0d0,1d0,0.5d0,2d0/
      do i=1,vmax
          sc4(i)=rnow-sc3(i)
      enddo
      do i=1,vmax
          sc4(i)=exp(-alph*sc4(i)*sc4(i))
      enddo
      xnorm=sqrt(sqrt((two*alph/pi)))
      call dscal(vmax,xnorm,sc4,1)
      return
      end
