*system:  F+H2, expansion of Stark-werner PES's r=1.4
*references: M. H. Alexander, J. Chem. Phys. 99, 6014 (1993).
* M.-L. Dubernet and J. M. Hutson, J. Chem. Phys. 101, 1939 (1994).
*
*  Note:  this pot routine requires a data file to be in hibxx/bin/progs/potdata:
*         three.param, two.param
*
      include "common/syusr"
      include "common/ground"
      include "common/bausr"
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      common /coselb/ ibasty
      common /conlam/ nlam, nlammx, lamnum(1)

      potnam='STARK-WERNER F-H2'
      ibasty=1
      lammin(1)=1
      lammax(1)=6
      nlam=lammax(1)-lammin(1)+1
      mproj(1)=0
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* --------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      character *48 potnam
      logical csflag, ljunk, ihomo, lljunk
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(3)
      potnam='STARK-WERNER F-H2'
      print *, potnam
      print *
1      print *, ' r (bohr)'
      read (5, *, end=93) r
      go to 1
      call pot(vv0,r)
      if (.not. csflag .or. (csflag .and. ihomo)) write (6, 100) vv0,vvl
100   format(4(1pe16.8))
      goto 1
93    r=3
      do i=1,100
       call pot(vv0,r)
       write(2,101) r,vv0,vvl
101    format(f8.4,6(1pe16.8))
       r=r+0.2
      enddo

99    end
      subroutine pot (vv0, r)
* ----------------------------------------------------------------------
*  subroutine to calculate the r-dependent coefficients in the
*  f-h2 potential of stark,werner at rh2=1.4 
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*  vv0       totally symmetric potential (v0)
*  variable in common block /covvl/
*    vvl:     vector of length 3 to store r-dependence of each term
*             in potential expansion

* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  9-oct-2006
* ----------------------------------------------------------------------


      implicit double precision (a-h,o-z)
      logical ljunk, ihomo, csflag, lljunk, flagd
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(3)
      dimension vs(4),vzz(4),vd(4),v12(4)
      dimension d0(16),d1(4),d2(9),aa(16)
      dimension xzz(4),xs(4),xd(3),x1(2),kpvt(4),
     :          qraux(4),work(24),rsd(4)
      dimension vvll(13)
      dimension hmat(7),dvmat(3,7)
*    vvll(1-4) expansion coefficients in dl0 (l=0,2,4,6) of vsigma
*    vvll(5-8) expansion coefficients in dl0 (l=0,2,4,6) of vsum
*    vvll(9-11) expansion coefficients in dl2 (l=2,4,6) of vdif
*    vvll(12-13) expansion coefficients in dl1 (l=2,4) of v1
* angles (not needed here, but included for completeness)
*      data beta /90d0,60d0,30d0,0d0/
      data zero, one, half /0.d0,1.d0,0.5d0/
      data two, three, four, five /2.d0,3.d0,4.d0,5.d0/
      data seven,eleven /7.d0,11.d0/
* coefficients for d0 rotation matrices
* stored (by column) for each of 4 angles and for l=0,2,4,6
      data d0/
     : 1d0, 1d0, 1d0, 1d0, -5d-1, -1.25d-1, 6.25d-1, 1d0, 3.75d-1,
     : -2.8906251d-1, 2.3437511d-2, 1d0, -3.125d-1, 3.2324218d-1,
     : -3.7402343d-1, 1d0/

* coefficients for d1 rotation matrices
* stored (by column) for each of 2 angles and for l=2,4
      data d1/
     : -5.3033009d-1, -5.3033009d-1, 3.0257682d-1, -5.4463828d-1/
* coefficients for d2 rotation matrices
* stored (by column) for each of 3 angles and for l=2,4,6
      data d2/
     : 6.1237244d-1, 4.5927933d-1, 1.5309311d-1, -3.9528471d-1,
     : 2.2234765d-1, 4.1999000d-1, 3.2021721d-1, -3.4523418d-1,
     : 4.8532921d-1/
      data ifirst /0/
      if (ifirst.eq.0) then
        call prepot
        ifirst=1
      endif

* determine potentials at angles
      do 200 i=1,4
        theta=(4-i)*30
        rh2=1.4d0
        flagd=.false.
        call vmerge(r,rh2,theta,hmat,dvmat,flagd)
*       hmat(1)=vap
*       hmat(2)=vs
*       hmat(3)=vxz
*       hmat(4)=va2p
*       hmat(5)=vstark
*       hmat(6)=aconst
*       hmat(7)=bconst
* divide vxz by sqrt(2), so that coupling
* is correct in the definite m basis
        v12(i)=hmat(3)/sqrt(two)
        vs(i)=half*(hmat(1)+hmat(4))
        vd(i)=half*(hmat(4)-hmat(1))
        vzz(i)=hmat(2)
*  here for stark-werner potential
        vzz(i)=hmat(5)
200   continue
      sum=0
* determine matrix of d's at angles for least-squares fit
*      call d2matev(9,beta,d0,d1,d2)
* (n.b. these are already in common, subroutine is included
*  only for completeness)
* solve simultaneous equations for solutions
* first for vzz
      tol=1.e-10
      call dcopy(16,d0,1,aa,1)
      call dqrank(aa,4,4,4,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,4,4,4,kr,vzz,xzz,rsd,kpvt,qraux)
      call dcopy(4,xzz,1,vvll,1)
* here is totally symmetric term
      vv0=vvll(1)
      call dcopy(3,vvll(2),1,vvl,1)
* convert to hartree
      econv=1./219474.6
      call dscal(3,econv,vvl,1)
      vv0=vv0*econv
      end
* --------------------------------------------------------------
       subroutine vmerge(r,rh2,theta,hmat,dvmat,flagd)
       implicit double precision (a-h,o-z)
       logical flagd
       dimension rx(11),rrx(8),scr(2,2),diag(2),rbond(3),
     :           dvs(3),dh12(3),dvap(3),dva2p(3),dvsum(3),dvdif(3),
     :           dbvsw(3),dvsw(3),dvmat(3,7),vvec(3),hmat(7),
     :           dswitch(3),dthdq(3),hmatold(4),dx(3),detilde(3,3),
     :           dxmat(3,3),tj(3,3),detilde1(3),detilde2(3),
     :           daconst(3),dbconst(3),dvxz(3)
       data rrx/1.05d0,1.2d0,1.3d0,1.4d0,1.5d0,1.65d0,1.8d0,2.3d0/
       data rx/3,3.3,3.675d0,4d0,4.5d0,5d0,5.5d0,6d0,7d0,8d0,10d0/
       data reh2 /1.4003d0/
       data deh2 /3.836486d+4/
       data zero,fourth,half,one,two,four,five
     :     /0.d0,0.25d0,0.5d0,1d0,2d0,4d0,5d0/
       data evtocm1,hartocm1 /8065.465d0,219474.6d0/
       data pi /3.141592653589793d0/
       call vpi(r,rh2,theta,vap,va2p,vsum,vdif,dvap,dva2p,dvsum,
     :          dvdif,flagd)
       call vsig(r,rh2,theta,vs,dvs,flagd)
       call v12(r,rh2,theta,h12,dh12,flagd)
* multiply v12 and its derivatives by sqrt(2), so that coupling
* is correct in the cartesian basis
* call the result vxz and dvxz
       sq2=sqrt(two)
       vxz=h12*sq2
       do i=1,3
         dvxz(i)=dh12(i)*sq2
       enddo
       call aso(r,rh2,aconst,bconst,daconst,dbconst,flagd)
       poth2=vh2(rh2-reh2)+deh2
* add poth2 to Vap,Vsig, Va2p
       vs=vs+poth2
       vap=vap+poth2
       va2p=va2p+poth2
* save [vap vs h12 va2p] in hmatold (hamiltonian matrix in cartesian
* basis
       hmatold(1)=vap
       hmatold(2)=vs
       hmatold(3)=vxz
       hmatold(4)=va2p
       if (flagd) then
          dpoth2=dvh2(rh2-reh2)
*  add derivatives of h2 potential to Vap, Vsig, Va2p
          dvs(2)=dvs(2)+dpoth2
          dvap(2)=dvap(2)+dpoth2
          dva2p(2)=dva2p(2)+dpoth2
       endif
* determine mixing angle
       thmix=half*atan2(two*vxz,(vs-vap))
       den=sqrt(four*vxz*vxz+(vs-vap)**2)

* cs2 and sn2 are cos(2*thmix) and sin(2*thmix)

       cs2th=(vs-vap)/den
       sn2th=2d0*vxz/den
*  use cos(thmix)=sqrt[0.5*(1+cs2th)] and sin(thmix)=sqrt[0.5*(1-cs2th)]
*      cs=cos(thmix)
*      sn=sin(thmix)

       csth2=half*(one+cs2th)
       snth2=half*(one-cs2th)
* determine adiabatic energies and add hydrogen energy
       av=half*(vs+vap)
       dif=half*(vs-vap)
       diag(1)=av-cs2th*dif-sn2th*vxz
       diag(2)=av+cs2th*dif+sn2th*vxz

* determine stark-werner potential contribution
* convert from jacobi to bond coordinates
       rbond(1)=rh2
       x=cos(theta*pi/180d0)
       rbond(2)=sqrt(r*r+0.25*rh2*rh2-x*r*rh2)
       rbond(3)=sqrt(r*r+0.25*rh2*rh2+x*r*rh2)
       if (flagd) then
* calculate sw surface and derivatives
         call fh2sw_deriv(rbond,vsw,dbvsw)
         vstark=hartocm1*vsw
* derivatives are in bond coordinates, these must to
* changed to jacobi coords
* construct jacobian matrix

         cs=cos(theta*pi/180)
         sn=sin(theta*pi/180)
         tj(1,1)=zero
         tj(1,3)=zero
         tj(1,2)=one
         tj(2,1)=two*r-rh2*cs
         tj(2,2)=half*rh2-r*cs
         tj(2,3)=rh2*r*sn
         call dscal(3,one/(two*rbond(2)),tj(2,1),3)
         tj(3,1)=two*r+rh2*cs
         tj(3,2)=half*rh2+r*cs
         tj(3,3)=-rh2*r*sn
         call dscal(3,one/(two*rbond(3)),tj(3,1),3)
*    calculate derivatives in jacobi frame
         do i=1,3
           dvsw(i)=ddot(3,dbvsw,1,tj(1,i),1)
         enddo
         call dscal(3,hartocm1,dvsw,1)
*   assemble derivative matrix (columns are Vap Vsig Vxz Va2p Vsw aso bso)
*   rows are R, rh2, theta (each element is derivative of
*   column potential by row coordinate)
         call dcopy(3,dvap,1,dvmat(1,1),1)
         call dcopy(3,dvs,1,dvmat(1,2),1)
         call dcopy(3,dvxz,1,dvmat(1,3),1)
         call dcopy(3,dva2p,1,dvmat(1,4),1)
         call dcopy(3,dvsw,1,dvmat(1,5),1)
         call dcopy(3,daconst,1,dvmat(1,6),1)
         call dcopy(3,dbconst,1,dvmat(1,7),1)
*  determine dtheta/dq
*      hmat(1)=vap
*      hmat(2)=vs
*      hmat(3)=vxz
*      hmat(4)=va2p
*      hmat(5)=vstarK
*      return
         vvec(1)=vxz
         vvec(2)=-vxz
         vvec(3)=(vs-vap)
         call dscal(3,one/(den*den),vvec,1)
         do i=1,3
           dthdq(i)=ddot(3,vvec,1,dvmat(i,1),3)
         enddo
       else
          call fh2sw(rbond,vau)
          vstark=hartocm1*vau
       endif

* switching function to ensure that lowest adiabatic surface is original
* stark werner fit outside of entrance channel
* smooth merge occurs at R=3.25 and r=3.25
       dr=r-3.25d0
       drh2=2.4-rh2
       onr=half*(tanh(four*dr)+one)
       onrh2=half*(tanh(five*drh2)+one)
       switch=onr*onrh2
* evaluate derivatives of switching function
       if (flagd) then
          dswitch(1)=onrh2*half*four*(1-tanh(four*dr)**2)
          dswitch(2)=-onr*half*five*(1-tanh(five*drh2)**2)
          dswitch(3)=zero
       endif
* stark-werner merges with lowest potential always
* determine which state is nominally sigma state,
* since first state in hmat is Pi(A') then if |cos(thmix)|>|sin(thmix)|
* then sigma state is 2nd state
       if (csth2.gt.snth2)then
*   here if sigma is second state
         istate=2
       else
         istate=1
       endif
* always lowest state merge
       istate=1
       etilde=diag(istate)
       diag(istate)=switch*etilde+(one-switch)*vstark
* convert back to diabatic energies
*      av=half*(diag(1)+diag(2))
*      dif=half*(diag(1)-diag(2))
*      vxz=sn2th*dif
*      vs=av+cs2th*dif
*      vap=av-cs2th*dif
       vap=csth2*diag(1)+snth2*diag(2)
       vs=csth2*diag(2)+snth2*diag(1)
       vxz=half*sn2th*(diag(2)-diag(1))

* save [vap vs vxz va2p aso bso] in hmat
       hmat(1)=vap
       hmat(2)=vs
       hmat(3)=vxz
       hmat(4)=va2p
       hmat(5)=vstark
       hmat(6)=aconst
       hmat(7)=bconst
* now evaluate derivative matrix in new diabatic basis
       if (flagd) then
* determine first diagonal elements of dEtilde
         do i=1,3
           detilde1(i)=
     :       csth2*dvmat(i,1)+snth2*dvmat(i,2)-sn2th*dvmat(i,3)
           detilde2(i)=
     :       snth2*dvmat(i,1)+csth2*dvmat(i,2)+sn2th*dvmat(i,3)
         enddo
* now determine T*dEtilde*T-transp
*        do i=1,3
*          detilde(i,1)=csth2*detilde1(i)+snth2*detilde2(i)
*          detilde(i,2)=snth2*detilde1(i)+csth2*detilde2(i)
*          detilde(i,3)=half*sn2th*(detilde2(i)-detilde1(i))
*        enddo
* now evaluate (dT/dq)*Etilde*T-transp + T*Etilde*(dT-transp)/dq
*      which equals (dT/dq)*T-transp*H+H*T*(dT-transp)dq
*      =[2vxz h22-h11;h22-h11 -2vxz]*dth/dq
*        do i=1,3
*          term=two*hmatold(3)*dthdq(i)
*          detilde(i,1)=term+detilde(i,1)
*          detilde(i,2)=-term+detilde(i,2)
*          termdiag=hmatold(2)-hmatold(1)
*          detilde(i,3)=termdiag*dthdq(i)+detilde(i,3)
*        enddo
* detild now contains in each row (d/dq)(T*eps-tilde*T-transpose)
* this is all commented out, since detilde=dvmat, which can
* be shown formally and numerically


* determine derivative of T*(1-switch)*(esw-e2tilde)*T-transpose
* first determine x
         x=(one-switch)*(vstark-etilde)
* now determine (dT/dq)*xmatrix*T-transp + T*xmatrix*(dT-transp)/dq
         vvec(1)=-sn2th
         vvec(2)=sn2th
         vvec(3)=-cs2th
         call dscal(3,x,vvec,1)
         if (istate.eq.2) call dscal(3,-one,vvec,1)
         do i=1,3
           do j=1,3
             dxmat(i,j)=dthdq(i)*vvec(j)
           enddo
         enddo
* dxmat now contains (dT/dq)*xmatrix*T-transp + T*xmatrix*(dT-transp)/dq
* now determine switch*(desw/dq-de2tilde/dq)-(dswitch/dq)*(esw-e2tilde)
         if (istate.eq.2) then
            vvec(1)=snth2
            vvec(2)=csth2
            vvec(3)=half*sn2th
         else
            vvec(1)=csth2
            vvec(2)=snth2
            vvec(3)=-half*sn2th
         endif
         do i=1,3
             if (istate.eq.1) then
                x=(one-switch)*(dvmat(i,5)-detilde1(i))
             else
                x=(one-switch)*(dvmat(i,5)-detilde2(i))
             endif
             x=x-dswitch(i)*(vstark-etilde)
             do j=1,3
               dxmat(i,j)=dxmat(i,j)+vvec(j)*x
             enddo
         enddo
* dxmat now contains (d/dq)[T*xmatrix*T-transp]

* finally determine VA', Vsig, and V12 elements of dH/dq
         do i=1,3
           do j=1,3
*             dvmat(i,j)=detilde(i,j)+dxmat(i,j)
* detilde=dvmat!
             dvmat(i,j)=dvmat(i,j)+dxmat(i,j)
           enddo
         enddo
       endif
       return
       end
* ------------------------------------------------------
       subroutine vpi(r,rh2,theta,vap,va2p,vvsum,vvdif,dvap,
     :                dva2p,dvsum,dvdif,flagd)
* --------------------------------------------------------------
* subroutine to determine F-H2 PiA' and PiA" diabatic potentials
* SEC with scale factor of 1.038
* based on Stark-Werner ab initio calculations supplementd by
* further point at r=2.3 bohr
* fitted by m.h. alexander as follows:
* at theta (30 60) the r and R dependence is fitted as

*     v(R,r)=exp(-xlam1*R)*(c1+c2*R+c3*R^2)
*           +(r-re)*exp(-xlam2*R)*(c4+c5*R+c6*R^2)
*           +(r-re)^2*exp(-xlam3*R)*(c7+c8*R+c9*R^2)
*           +(r-re)^3*exp(-xlam4*R)*(c10+c11*R+c12*R^2)
*           -0.5*c13*(tanh(alph*(R-xlam5))+1)*R^(-ipow)
*           -0.5*c14*(r-re)*(tanh(alph*(R-xlam6))+1)*R^(-ipow)
*           -0.5*c15*(r-re)^2*(tanh(alph*(R-xlam7))+1)*R^(-ipow)
* involving 15 linear and 7 non-linear parameters
* where ipow is 5
* for a given value of R=Ri and r=ri, the v(Ri,ri)'s are determined at each
* theta from the above equation, using parameters given below. subseqently
* the average [0.5*(VA'+VA")] is fitted to an expansion in even legendre
* polynomials (l=0,2,4,6) and the difference potential [0.5*(VA"-VA')] is
* fitted to an expansion in even rotation matrix elments of order 2
* (l=2,4,6)
* with this, the potentials vA'(Ri,ri,thetai) and vA" are determined
* errors are
*  mean A' error is 1.71 cm-1 for all points with E < 20,000 cm-1
*  mean A" error is 1.62 cm-1 for all points with E < 20,000 cm-1
*  max A' error is 9.54 cm-1 for all points with E < 20,000 cm-1
*  max A" error is 9.54 cm-1 for all points with E < 20,000 cm-1
*  mean A' relative error is
*     0.16 % for points > 50 cm-1
*     1.18 % for points > 10 cm-1
*  mean A" relative error is
*     0.16 % for points > 50 cm-1
*     1.06 % for points > 10 cm-1
* at small R all terms in vdif and the anisotropies in vsum are
* killed off by the function
*     half*(tanh(4.5d0*(r-2.25d0))+one)
* if flagd=.true., then derivatives are calculated
*   and returned in the vectors dvap,dva2p,dvsum,dvdif with
*   dvap(1)=dvap/dR
*   dvap(2)=dvap/dr
*   dvap(3)=dvap/dtheta
* and similarly for dva2p, dsum, and ddif
* --------------------------------------------------------------
* current revision date:  25-mar-1998
* --------------------------------------------------------------
       dimension xlam(7,7),cc(15,7),xxlam(7),ccc(15),
     :      kpvt(4),cmat(8,4),vsum(4),vdif(4),dvsumr(4),
     :       dvdifr(4),dvsumrh2(4),dvdifrh2(4),dvsum(3),dvdif(3),
     :       xdif(4),xsum(4),pleg(4),dxdifr(4),dxdifrh2(4),
     :       scr1(8),scr2(8),d2mat(9),d2matx(9),xtilsum(4),
     :       xtildif(4),dpleg(4),dd2matx(4),dv(2),
     :       d0mat(16),d0matx(16),dvap(3),dva2p(3),
     :       dxsumrh2(4),dxsumr(4)
       implicit double precision (a-h,o-z)
       logical flagd
       data one,half /1.d0,0.5d0/
       data pi /3.141592653589793d0/
* scale factor for tanh in potential
       data alpha /1.2d0/
* power for long-range tail of potential
       data ipow /5/
* equilibrium h2 distance for potential expansion
       data re /1.4d0/
       data ione /1/
* xlam for pi states rows are lam1,lam2,lam3,lam4...
* columns are A"0 A"30 A"60 A"90 A'30 A'60 A'90
       data xlam /
     : 1.5012129d0, 1.6287682d0, 1.6936217d0, 1.7070310d0, 1.6547615d0,
     : 1.7417981d0, 1.7547183d0, 1.6329877d0, 1.6555711d0, 1.6799885d0,
     : 1.7358484d0, 1.6589411d0, 1.7353634d0, 1.6772495d0, 1.5487196d0,
     : 1.6400137d0, 1.4253960d0, 1.8603522d0, 1.6434217d0, 1.8310540d0,
     : 1.7119064d0, 2.3760952d0, 2.3814315d0, 1.9583216d0, 1.4208832d0,
     : 2.1741882d0, 2.4451020d0, 1.9833472d0, 3.8760026d0, 3.7729244d0,
     : 3.7511013d0, 3.7667134d0, 3.7639745d0, 3.6582959d0, 3.6393376d0,
     : 4.4275460d0, 4.1964861d0, 4.2218247d0, 4.6668232d0, 4.3013701d0,
     : 4.7725499d0, 6.0538442d0, 4.9436300d0, 4.8319987d0, 4.8550601d0,
     : 4.7516793d0, 4.9310403d0, 4.5645319d0, 4.5892214d0/
* linear least-squares coefficients
* rows are A"0 A"30 A"60 A"90 A'30 A'60 A'90
* columns are cc1, cc2 .... cc15
      data cc /
     : 2.8819376d6, -4.3568113d5, -4.0935313d3, -3.4138053d6,
     : 2.1654806d6, -2.5273424d5, -7.3012645d5, 2.0329298d5,
     : 8.9537238d3, 3.2517003d7, -1.9899309d7, 3.1071968d6, 4.5487535d5,
     : 3.1180598d5, 3.4849464d5, 2.2486795d6, 2.2815369d5, -9.2157906d4,
     : -4.0724490d6, 2.3550307d6, -2.6599691d5, -1.1902096d6,
     : 3.2400591d5, 3.2983549d3, 2.0936769d7, -1.2301323d7, 1.8241995d6,
     : 3.9603310d5, 2.2483935d5, 1.8624462d5, 1.7938951d6, 5.0969679d5,
     : -1.2062799d5, -3.6192282d6, 1.8984601d6, -2.0288261d5,
     : -4.6101564d5, 1.0633655d5, -1.0954960d2, 1.2142117d6,
     : -4.9087964d5, 3.7420861d4, 2.6955403d5, 8.3034368d4, 1.6410336d5,
     : 1.7788259d6, 4.8020272d5, -1.1494680d5, -3.3249163d6,
     : 1.6050433d6, -1.5182676d5, 1.2462221d6, -1.1250167d6,
     : 1.7504299d5, 2.3711983d5, -1.0746451d5, 1.0632103d4, 1.9394168d5,
     : 8.1644790d4, 3.3594887d4, 2.1282293d6, 4.0633130d5, -1.1643584d5,
     : -4.2624662d6, 2.4786655d6, -2.8283972d5, -1.2651640d6,
     : 3.3525537d5, 5.3553478d3, 9.6177305d6, -5.4512635d6, 7.6458317d5,
     : 3.6922607d5, 1.6304627d5, 2.1337028d5, 1.5227005d6, 8.9468105d5,
     : -1.7050041d5, -3.9865594d6, 2.1047725d6, -2.1396072d5,
     : 7.8465262d4, -5.6069259d5, 1.2348260d5, -5.7864824d6,
     : 4.4325447d6, -8.6066043d5, 2.0190360d5, 5.4024202d4, 4.6511940d4,
     : 1.6352957d6, 8.0301556d5, -1.5724366d5, -2.6665395d6,
     : 1.3708732d6, -1.3700231d5, 2.2689823d5, -4.1331114d5,
     : 7.3944609d4, -2.5758475d5, 3.2739322d5, -8.4466909d4,
     : 1.1428191d5, -5.0567955d4, 4.4146573d4/


*  this is matrix of dl2(theta) at 30, 60 and 90 degrees
*  l=2 first row, l=4 2nd row, l=6 3rd row
       data d2mat /
     : 1.5309311d-1, 4.5927933d-1, 6.1237244d-1, 4.1999000d-1,
     : 2.2234765d-1, -3.9528471d-1, 4.8532921d-1, -3.4523418d-1,
     : 3.2021721d-1/
*  this is matrix of dl0(theta) at 0, 30, 60 and 90 degrees
*  l=0 first row, l=2 2nd row, l=4 3rd row, l=6 4th row
       data d0mat/
     : 1d0, 1d0, 1d0, 1d0, 1d0, 6.25d-1, -1.25d-1, -5.d-1, 1d0,
     : 2.3437511d-2, -2.8906251d-1, 3.75d-1, 1d0, -3.7402343d-1,
     : 3.2324218d-1, -3.125d-1/

* determine potential at theta=0, 30, 60, and 90 values needed for interpolati
       do ith=1,4
         call dcopy(7,xlam(ith,1),7,xxlam,1)
         call dcopy(15,cc(1,ith),1,ccc,1)
         call dvpiexp(r,rh2,re,alpha,xxlam,ccc,v,dv,ipow,flagd)
         va2p=v
         dva2pr=dv(1)
         dva2prh2=dv(2)
         if (ith .gt. 1) then
           call dcopy(7,xlam(ith+3,1),7,xxlam,1)
           call dcopy(15,cc(1,ith+3),1,ccc,1)
           call dvpiexp(r,rh2,re,alpha,xxlam,ccc,v,dv,ipow,flagd)
           vap=v
           dvapr=dv(1)
           dvaprh2=dv(2)
         else
           vap=va2p
           if (flagd) then
             dvapr=dva2pr
             dvaprh2=dva2prh2
           endif
         endif
         if (flagd) then
            dvsumr(ith)=half*(dvapr+dva2pr)
            dvdifr(ith)=half*(dva2pr-dvapr)
            dvsumrh2(ith)=half*(dvaprh2+dva2prh2)
            dvdifrh2(ith)=half*(dva2prh2-dvaprh2)
         endif
         vdif(ith)=half*(va2p-vap)
         vsum(ith)=half*(vap+va2p)
       enddo

* calculate legendre expansion coefficients for vsum and vdiff
       call dcopy(9,d2mat,1,d2matx,1)
       call dqrank(d2matx,3,3,3,tol,kr,kpvt,scr1,scr2)
       call dqrlss(d2matx,3,3,3,kr,vdif(2),xdif(2),scr2,kpvt,scr1)
       xdif(1)=zero
       if (flagd) then
         call dqrlss(d2matx,3,3,3,kr,dvdifr(2),dxdifr(2),scr2,
     :               kpvt,scr1)
         call dqrlss(d2matx,3,3,3,kr,dvdifrh2(2),dxdifrh2(2),scr2,
     :               kpvt,scr1)
         dxdifr(1)=zero
         dxdifrh2(1)=zero
       endif
       call dcopy(16,d0mat,1,d0matx,1)
       call dqrank(d0matx,4,4,4,tol,kr,kpvt,scr1,scr2)
       call dqrlss(d0matx,4,4,4,kr,vsum,xsum,scr2,kpvt,scr1)
       if (flagd) then
         call dqrlss(d0matx,4,4,4,kr,dvsumr,dxsumr,scr2,
     :               kpvt,scr1)
         call dqrlss(d0matx,4,4,4,kr,dvsumrh2,dxsumrh2,scr2,
     :               kpvt,scr1)

       endif
* determine legendre functions at theta
* x is cos(theta)
       x=cos(theta*pi/180d0)
       x2=x*x
       sn2=one-x2
       x4=x2*x2
       x6=x2*x4
       x8=x4*x4
       pleg(1)=one
       pleg(2)=(3d0*x2-one)*half
       pleg(3)=(3.d0-30d0*x2+35*x4)/8d0
       pleg(4)=(-5d0+105d0*x2-315d0*x4+231d0*x6)/16d0
       d2matx(1)=sqrt(1.5d0)*sn2*half
       d2matx(2)=-sqrt(10d0)*sn2*(one-7d0*x2)/8d0
       d2matx(3)=sqrt(105d0)*sn2*(one-18d0*x2+33d0*x4)/32d0
* kill off higher anisotropies at short range (small R)
       on=half*(tanh(4.5d0*(r-2.25d0))+one)
* save original xsum and xdif if needed for derivatives
       if (flagd) then
         onp=half*4.5d0*(one-(tanh(4.5d0*(r-2.25d0)))**2)
         call dcopy(4,xsum,1,xtilsum,1)
         call dcopy(4,xdif,1,xtildif,1)
* multiply xtilde's by onp
         call dscal(3,onp,xtilsum(2),1)
         call dscal(3,onp,xtildif(2),1)
* first correction is zero, since damping does not involve isotropic
* term in potential
         xtilsum(1)=zero
         xtildif(1)=zero
       endif
* kill off anisotropies
       call dscal(3,on,xsum(2),1)
       call dscal(3,on,xdif(2),1)
* determine potential at theta
       vvsum=ddot(4,xsum,1,pleg,1)
       vvdif=ddot(3,xdif(2),1,d2matx,1)
       va2p=vvsum+vvdif
       vap=vvsum-vvdif
       if (flagd) then
* scale dxsumr, dxsumrh2, dxdifr, and dxdifrh2
         call dscal(3,on,dxsumr(2),1)
         call dscal(3,on,dxsumrh2(2),1)
         call dscal(3,on,dxdifr(2),1)
         call dscal(3,on,dxdifrh2(2),1)
* determine derivatives with respect to rh2
         dvvsumrh2=ddot(4,dxsumrh2,1,pleg,1)
         dvvdifrh2=ddot(3,dxdifrh2(2),1,d2matx,1)
* determine derivatives with respect to R
         dvvsumr=ddot(4,dxsumr,1,pleg,1)+ddot(4,xtilsum,1,pleg,1)
         dvvdifr=ddot(3,dxdifr(2),1,d2matx,1)+
     :           ddot(3,xtildif(2),1,d2matx,1)
* determine derivatives with respect to theta
         sn=sin(theta*pi/180d0)
         dpleg(1)=zero
         dpleg(2)=-3d0*x*sn
         dpleg(3)=(2.5d0)*x*sn*(3d0-7d0*x2)
         dpleg(4)=-(2.625d0)*x*sn*(5d0-30d0*x2+33d0*x4)
         dd2matx(1)=sqrt(1.5d0)*sn*x
         dd2matx(2)=-half*sqrt(10d0)*sn*x*(4d0-7d0*x2)
         dd2matx(3)=sqrt(105d0)*sn*x*(19d0-102d0*x2+99d0*x4)/16d0
         dvvsumth=ddot(4,xsum,1,dpleg,1)
         dvvdifth=ddot(3,xdif(2),1,dd2matx,1)
         dvsum(1)=dvvsumr
         dvdif(1)=dvvdifr
         dvsum(2)=dvvsumrh2
         dvdif(2)=dvvdifrh2
         dvsum(3)=dvvsumth
         dvdif(3)=dvvdifth
         do i=1,3
           dva2p(i)=dvsum(i)+dvdif(i)
           dvap(i)=dvsum(i)-dvdif(i)
         enddo
       endif
       return
       end
* ------------------------------------------------------
       subroutine dvpiexp(r,rh2,re,alph,xlam,c,v,dv,ipow,flagd)

*  to determine v(R,r), dv(R,r)/dr and dv(R,r)/dR, where
*     v(R,r)=exp(-xlam1*R)*(c1+c2*R+c3*R^2)
*           +(r-re)*exp(-xlam2*R)*(c4+c5*R+c6*R^2)
*           +(r-re)^2*exp(-xlam3*R)*(c7+c8*R+c9*R^2)
*           +(r-re)^3*exp(-xlam4*R)*(c10+c11*R+c12*R^2)
*           -0.5*c13*(tanh(alph*(R-xlam5))+1)*R^(-ipow)
*           -0.5*c14*(r-re)*(tanh(alph*(R-xlam6))+1)*R^(-ipow)
*           -0.5*c15*(r-re)^2*(tanh(alph*(R-xlam7))+1)*R^(-ipow)
* on return:  dv(1) is dv/dR and dv(2) is dv/dr
* derivatives determined only if flagd=.true.
       implicit double precision (a-h,o-z)
       logical flagd
       dimension xlam(7),c(15),dv(2),ex(4)
       data one,half,two/1.d0,0.5d0,2d0/
       rmre=rh2-re
       rmre2=rmre*rmre
       rinv=r**(-ipow)
       r2=r*r
       vex=zero
       dvexr=zero
       dvexrh2=zero
       rho=one
       rho2=one
       do i=1,4
         exx=exp(-xlam(i)*r)
         ind=(i-1)*3+1
         vexterm=exx*(c(ind)+c(ind+1)*r+c(ind+2)*r2)
         vex=vex+rho*vexterm
         if (flagd) then
           dvexr=dvexr+rho*exx*(c(ind+1)+two*c(ind+2)*r
     :           -xlam(i)*(c(ind)+c(ind+1)*r+c(ind+2)*r2))
           if (i.gt.1) then
             dvexrh2=dvexrh2+(i-1)*rho2*vexterm
             rho2=rho2*rmre
           endif
         endif
         rho=rho*rmre
       enddo
       vlr=zero
       dvlr=zero
       dvlrrh2=zero
       rho=one
       rho2=one
       do i=1,3
         ind=i+4
         exx=tanh(alph*(r-xlam(ind)))
         vlrterm=-half*(exx+one)*c(i+12)*rinv
         vlr=vlr+rho*vlrterm
         if (flagd) then
           if (i.gt.1) then
             dvlrrh2=dvlrrh2+(i-1)*rho2*vlrterm
             rho2=rho2*rmre
           endif
           dvlr=dvlr-rho*half*c(i+12)*rinv*(alph*(one-exx**2)
     :             -ipow*(1+exx)/r)
         endif
         rho=rho*rmre
       enddo
       v=vex+vlr
       dv(1)=dvexr+dvlr
       dv(2)=dvexrh2+dvlrrh2

       return
       end
* ------------------------------------------------------
       subroutine v12(r,rh2,theta,v,dvv,flagd)
* --------------------------------------------------------------
* subroutine to determine F-H2 diabatic coupling potential
* SEC with scale factor of 1.038
* based on Stark-Werner ab initio calculations supplementd by
* further point at r=2.3 bohr
* fitted by m.h. alexander as follows:
* at each r (1.05 1.2 1.3 1.4 1.5 1.65 1.8 2.3) and theta (30 60)
* the R dependence is fitted as

*     v(R)=exp(-xlam1*R)*(c1+c2*R+c3*R^2+c4*R^3)
*          -0.5*(tanh(alph*(R-xlam2))+1)*R^(-ipow)
* where ipow is 6 here
* for a given value of R=ri, the v(R)'s are determined at each r
* from the above equation, using parameters given below.  with this the
* r dependence for each theta is fitted to a cubic in r by a linear
* least-squares procedure (subroutines dqrank and dqrlss) to obtain
* v(Ri,ri) for theta = 30 and 60.  then, the theta dependence is fit
* to an expansion in even rotational matrix elements of order 1
*    v(Ri,ri,theta)= c2*v2(Ri,ri)*d12(theta) + c4*v4(Ri,ri)*d14(theta)
* with this, the potential v(Ri,ri,thetai) is determined
* errors are
*  mean error is 0.69 cm-1 for all points with R > 3.5
*  max error is 2.55 cm-1 for all points with R > 3.5
*  mean relative error is
*     1.3 % for points > 50 cm-1
*     1.9 % for points > 10 cm-1
*     2.4 % for points >  5 cm-1
* at small R the potential is damped to zero by the function
*     half*(tanh(4d0*(r-3.25d0))+one)
* at large r the potential is damped to zero by the function
*     half*(one-tanh(alphlr*(rh2-rlongr)))
* if flagd=.true., then derivatives are calculated
*   and returned in the vector dvv with
*   dvv(1)=dv12/dR
*   dvv(2)=dv12/dr
*   dvv(3)=dv12/dtheta
* --------------------------------------------------------------
* current revision date:  24-mar-1998
* --------------------------------------------------------------
       dimension xlam130(8),xlam230(8),c130(8),c230(8),c330(8),c430(8),
     :       clr30(8),vang(2),xang(2),dvv(3),dvang(2),dxang(2),
     :       xlam160(8),xlam260(8),c160(8),c260(8),c360(8),c460(8),
     :       clr60(8),xlam(2),cc(5),xtil1230(4),xtil1260(4),
     :       rr(8),v1230(8),v1260(8),kpvt(4),cmat(8,4),
     :       scr1(8),scr2(8),x1230(4),x1260(4),d1mat(4),d1matx(4),
     :       dv1230(8),dv1260(8),dx1230(4),dx1260(4),dvangdr(2),
     :       dxangdr(2)
       implicit double precision (a-h,o-z)
       logical flagd
       data one,half /1.d0,0.5d0/
       data pi /3.141592653589793d0/
* scale factor for tanh in potential
       data alpha /1.2d0/
* power for long-range tail of potential
       data ipow /6/
* scale factor for long-range damping
       data alphlr, rlongr /4d0,3.2d0/
       data ione /1/
       data rr/1.05d0,1.2d0,1.3d0,1.4d0,1.5d0,1.65d0,1.8d0,2.3d0/
       data xlam130/
     : 1.2128797d0, 1.6596660d0, 1.6968504d0, 1.7231682d0, 1.7300095d0,
     : 1.4616455d0, 1.4680763d0, 1.1694348d0/
       data xlam230/
     : 4.5628549d0, 3.4929870d0, 3.5372339d0, 3.5917555d0, 3.6609427d0,
     : 3.5623835d0, 3.6324751d0, 4.1909152d0/
       data c130/
     : -5.6801449d+4, 6.1015109d+4, 1.1892988d+5, 1.9722603d+5,
     : 2.8646479d+5, 1.8877177d+6, 2.7558414d+6, 1.3939615d+6/
       data c230/
     : 1.8401036d+4, -4.3742185d+4, -6.6683004d+4, -9.4959586d+4,
     : -1.2419877d+5, -1.2521877d+6, -1.7931946d+6, -8.0840451d+5/
       data c330/
     : -2.7349040d-1, 1.0440374d0, 9.2507933d-1, -2.2502430d0,
     : 2.3311312d0, 2.5116642d+5, 3.5764668d+5, 1.3939304d+5/
       data c430/
     : 0d0, 0d0, 0d0, 0d0, 0d0, -1.7686714d+4, -2.5062789d+4,
     : -8.2248604d+3/
       data clr30/
     : 1.0510083d+6, -1.2867314d+6, -1.4630225d+6, -1.6557502d+6,
     : -1.8876647d+6, -4.3428791d+6, -5.3785031d+6, -1.1111429d+7/
       data xlam160/
     : 1.6496279d0, 1.6369425d0, 1.7750230d0, 1.8139836d0, 1.8523041d0,
     : 1.9096922d0, 1.9739429d0, 1.0893333d0/
       data xlam260/
     : 3.4441465d0, 3.4788439d0, 3.5036506d0, 3.5333710d0, 3.5635600d0,
     : 3.6098889d0, 3.6531504d0, 4.2765898d0/
       data c160/
     : 2.4854501d+4, 8.3595709d+4, 1.7894607d+5, 2.9426785d+5,
     : 4.5676329d+5, 8.2381803d+5, 1.4378460d+6, 2.4146519d+5/
       data c260/
     : -2.6550043d+4, -5.7488017d+4, -9.1373407d+4, -1.3464132d+5,
     : -1.9324069d+5, -3.2056139d+5, -5.2746651d+5, -1.3955835d+5/
       data c360/
     : -9.2119864d-1, 3.2694307d+3, -1.8042500d0, 2.5393929d-2,
     : 4.2828868d0, 1.3403758d+1, 1.2298473d+1, 2.3004425d+4/
       data c460/
     : 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, -1.2162794d+3/
       data clr60/
     : -9.4327281d+5, -1.1524003d+6, -1.3071537d+6, -1.4584133d+6,
     : -1.6065721d+6, -1.8167485d+6, -1.9970433d+6, -3.7681636d+6/
*  this is matrix of d1m(theta) at 30 and 60 degrees
*  m=2 first row, m=4 2nd row
       data d1mat /
     :  -5.3033009d-1, -5.3033009d-1, -5.4463828d-1, 3.0257682d-1/
* this version uses all values of r, and fits with cubic (least squares)
* determine potential at theta=30 and 60 and r values needed for interpolation
       do ir=1,8
         xlam(1)=xlam130(ir)
         xlam(2)=xlam230(ir)		
         cc(1)=c130(ir)
         cc(2)=c230(ir)
         cc(3)=c330(ir)
         cc(4)=c430(ir)
         cc(5)=clr30(ir)
         v1230(ir)= voneexps(r,alpha,xlam,cc,ipow)
         if (flagd) dv1230(ir)=dvoneexps(r,alpha,xlam,cc,ipow)
         xlam(1)=xlam160(ir)
         xlam(2)=xlam260(ir)		
         cc(1)=c160(ir)
         cc(2)=c260(ir)
         cc(3)=c360(ir)
         cc(4)=c460(ir)
         cc(5)=clr60(ir)
         v1260(ir)= voneexps(r,alpha,xlam,cc,ipow)
         if (flagd) dv1260(ir)=dvoneexps(r,alpha,xlam,cc,ipow)
         cmat(ir,1)=one
         cmat(ir,2)=rr(ir)
         cmat(ir,3)=rr(ir)*cmat(ir,2)
         cmat(ir,4)=rr(ir)*cmat(ir,3)
       enddo
* calculate short R cutoff
       on=half*(tanh(4d0*(r-3.25d0))+one)
* derivative of switching function with respect to R
       if (flagd) then
          onp=2d0*(one-(tanh(4d0*(r-3.25d0)))**2)
          do ir=1,8
            dv1230(ir)=on*dv1230(ir)+onp*v1230(ir)
            dv1260(ir)=on*dv1260(ir)+onp*v1260(ir)
          enddo
       endif
* now scale v1230 and v1260
       call dscal(8,on,v1230,1)
       call dscal(8,on,v1260,1)

* determine r expansion coefficients (cubic fit)
       tol=1.d-10
       call dqrank(cmat,8,8,4,tol,kr,kpvt,scr1,scr2)
       call dqrlss(cmat,8,8,4,kr,v1230,x1230,scr2,kpvt,scr1)
       call dqrlss(cmat,8,8,4,kr,v1260,x1260,scr2,kpvt,scr1)
       if (flagd) then
         call dqrlss(cmat,8,8,4,kr,dv1230,dx1230,scr2,kpvt,scr1)
         call dqrlss(cmat,8,8,4,kr,dv1260,dx1260,scr2,kpvt,scr1)
       endif
* calculate large r damping
       damp=half*(one-tanh(alphlr*(rh2-rlongr)))
* xtil will hold original undamped coefficients
       call dcopy(4,x1230,1,xtil1230,1)
       call dcopy(4,x1260,1,xtil1260,1)
       call dscal(4,damp,x1230,1)
       call dscal(4,damp,x1260,1)
* april 2 correction
       if (flagd) then
          call dscal(4,damp,dx1230,1)
          call dscal(4,damp,dx1260,1)
       endif
* evaluate derivative of damping function and multiply original
* expansion coefficients by this
       if (flagd) then
          ddamp=-half*alphlr*(one-(tanh(alphlr*(rh2-rlongr)))**2)
          call dscal(4,ddamp,xtil1230,1)
          call dscal(4,ddamp,xtil1260,1)
       endif
* interpolate in r (cubic)
       vv1230=x1230(1)
       vv1260=x1260(1)
       if (flagd) then
          dvv1230=dx1230(1)
          dvv1260=dx1260(1)
       endif
       rterm=rh2
       do ir=2,4
         vv1230=vv1230+x1230(ir)*rterm
         vv1260=vv1260+x1260(ir)*rterm
         if (flagd) then
            dvv1230=dvv1230+dx1230(ir)*rterm
            dvv1260=dvv1260+dx1260(ir)*rterm
         endif
         rterm=rterm*rh2
       enddo
* vang is potential at r=rh2
       vang(1)=vv1230
       vang(2)=vv1260
       if (flagd) then
* dvangdR is derivative with respect to R of potential at r=rh2
          dvangdr(1)=dvv1230
          dvangdr(2)=dvv1260
       endif
* calculate derivative in r (if desired)
* f(r)=x(1)+x(2)*r+x(3)*r^2+x(4)*r^3
* df/dr=x(2)+2*x(3)*r+3*x(4)*r^2+
*       x'(1)+x'(2)*r+x'(3)*r^2+x'(4)*r^3
* where x(i)=xtilde(i)*damp
* and   x'(i)=xtilde(i)*ddamp
* NB at this point the array x contains xtilde*damp
*  and the array xtilde contains xtilde*ddamp
       if (flagd) then
          rterm=one
          vv1230=x1230(2)
          vv1260=x1260(2)
          vvtil1230=xtil1230(1)
          vvtil1260=xtil1260(1)
          do ir=2,4
            if (ir.gt.2) then
              vv1230=vv1230+x1230(ir)*rterm*(ir-1)
              vv1260=vv1260+x1260(ir)*rterm*(ir-1)
            endif
            vvtil1230=vvtil1230+xtil1230(ir)*rterm*rh2
            vvtil1260=vvtil1260+xtil1260(ir)*rterm*rh2
            rterm=rterm*rh2
          enddo
          dvang(1)=vvtil1230+vv1230
          dvang(2)=vvtil1260+vv1260
*  dvang now contains derivative of potential with respect to r
*  including damping function
        endif

* calculate  expansion coefficients for d12(theta) and d14(theta)
       call dcopy(4,d1mat,1,d1matx,1)
       call dqrank(d1matx,2,2,2,tol,kr,kpvt,scr1,scr2)
       call dqrlss(d1matx,2,2,2,kr,vang,xang,scr2,kpvt,scr1)
       if (flagd) then
          call dqrlss(d1matx,2,2,2,kr,dvang,dxang,scr2,kpvt,scr1)
          call dqrlss(d1matx,2,2,2,kr,dvangdr,dxangdr,scr2,kpvt,scr1)
       endif
       x=cos(theta*pi/180d0)
       sn=sqrt(one-x*x)
       d2=-sqrt(1.5d0)*sn*x
       d4=0.25d0*sqrt(5d0)*sn*x*(3d0-7d0*x*x)
       v=xang(1)*d2+xang(2)*d4
       if (flagd) then
          dd2=-sqrt(1.5d0)*(2d0*x*x-one)
          dd4=-0.25d0*sqrt(5d0)*(-27d0*x*x+28d0*x**4+3d0)
          dvv(3)=xang(1)*dd2+xang(2)*dd4
          dvv(1)=dxangdr(1)*d2+dxangdr(2)*d4
          dvv(2)=dxang(1)*d2+dxang(2)*d4
       endif
       return
       end
* ------------------------------------------------------
       subroutine vsig(r,rh2,theta,vv,dvv,flagd)
* --------------------------------------------------------------
* subroutine to determine F-H2 diabatic potential for sigma(A')
* SEC with scale factor of 1.038
* based on Stark-Werner ab initio calculations supplementd by
* further point at r=2.3 bohr
* fitted by m.h. alexander as follows:
* at each r (1.05 1.2 1.3 1.4 1.5 1.65 1.8 2.3) and theta (0,30,60,90)
* the R dependence is fitted as

*     v(R)=exp(-xlam1*R)*(c1+c2*R+c3*R^2+c4*R^3)
*          -0.5*(tanh(alph*(R-xlam2))+1)*R^(-ipow)
* where ipow is 5 here
* for a given value of R=ri, the v(R)'s are determined at each r
* from the above equation, using parameters given below.  with this the
* r dependence for each theta is fitted to a quartic in r-re by a linear
* least-squares procedure (subroutines dqrank and dqrlss) to obtain
* v(Ri,ri) for theta = [0:30:90].  then, the theta dependence is fit
* to an expansion in even legendre polynomials (l=0,2,4,6)
* with this, the potential v(Ri,ri,thetai) is determined
* errors are
*  mean error is 1.74 cm-1 for all points with R > 3
*  max error is 11.93 cm-1 for all points with R > 3
*  mean relative error is
*     0.53 % for points > 50 cm-1
*     2.29 % for points > 10 cm-1
* at small R the potential is damped to zero by the function
*     half*(tanh(4d0*(r-3.25d0))+one)
* at large r the quadratic (and higher) in r-re terms in the potential are
* damped to zero by the function
*     half*(one-tanh(alphlr*(rh2-rlongr)))
* if flagd=.true., then derivatives are calculated
*   and returned in the vector dvv with
*   dvv(1)=dvsig/dR
*   dvv(2)=dvsig/dr
*   dvv(3)=dvsig/dtheta
* --------------------------------------------------------------
* current revision date:  24-mar-1998
* --------------------------------------------------------------
       dimension cc(34,5),xlam(34,2),xxlam(2),ccc(5),rbond(3),
     :       rr(8),v(9,4),xth(5),kpvt(5),cmat(8,5),vang(4),
     :       scr1(8),scr2(8),d0mat(16),d0matx(16),dv(9,4),
     :       pleg(4),xang(4),dpleg(4),dvv(3),dvang(4),dxang(4),
     :       dxth(5),dvangdr(4),dxangdr(4),xtilde(5)

       implicit double precision (a-h,o-z)
       logical flagd
       data one,half /1.d0,0.5d0/
       data pi /3.141592653589793d0/
* scale factor for tanh in potential
       data alpha /1.2d0/
* power for long-range tail of potential
       data ipow /5/
* scale factor for long-range damping
       data alphlr, rlongr /4d0,3.2d0/
* conversion from ev to cm-1
       data tocm1 / 8065.465d0/
       data ione /1/
       data rr/1.05d0,1.2d0,1.3d0,1.4d0,1.5d0,1.65d0,1.8d0,2.3d0/

       data cc/
     : 6.2890634d8, 7.6040274d8, 1.0523278d9, 1.4650321d9,
     : 1.8011260d9, 9.9204808d8, -2.4752628e+11, 5.8977005d6,
     : -1.5946205d8, -2.6026113d6, 8.8612483d8, 1.0437225d9,
     : 1.2892826d9, 1.6778238d9, 1.6786030d9, 7.6413524d8, 6.0944475d6,
     : -5.5647822d7, -1.7920106d5, 5.4290184d6, 1.2108633d7,
     : 2.2399485d7, 2.7601340d9, 4.4686953d9, 1.2401323d8,
     : -6.5682107e+10, 4.1219459d5, 5.0969498d6, 4.7134535d6,
     : 7.7425754d6, 1.2554585d7, 1.4732907d7, 4.3265734d7, 1.4117897d8,
     : -4.8225573d8, -5.5991253d8, -7.4804079d8, -1.0056341d9,
     : -1.1972327d9, -8.0486766d8, 2.3775582e+11, -4.0142512d6,
     : 5.7913384d7, 2.0462029d6, -6.5725324d8, -7.5277341d8,
     : -9.0169693d8, -1.1359806d9, -1.0901617d9, -6.2746443d8,
     : -4.3722529d6, 2.1109422d7, 4.6365860d5, -5.9572937d6,
     : -1.2652830d7, -2.2587522d7, -2.4139518d9, -3.8182144d9,
     : -1.1356859d8, 6.4510656e+10, -1.1020606d5, -5.1690783d6,
     : -4.8565749d6, -7.8643370d6, -1.2513404d7, -1.4543110d7,
     : -3.9262269d7, -1.2061158d8, 9.7658037d7, 1.0626178d8,
     : 1.3481868d8, 1.7252768d8, 1.9635525d8, 2.0721959d8,
     : -7.5258311e+10, 6.4548121d5, -5.3488387d6, -2.8338959d5,
     : 1.2603968d8, 1.3826836d8, 1.5822477d8, 1.9025041d8, 1.7147789d8,
     : 1.6426814d8, 7.2920222d5, -2.0488127d6, -8.1704654d4,
     : 2.3450695d6, 4.4831674d6, 7.5399556d6, 6.8223837d8, 1.0566813d9,
     : 3.3594592d7, -2.0768586e+10, 5.7405907d3, 1.9744781d6,
     : 1.8174179d6, 2.7621059d6, 4.2020095d6, 4.7354412d6, 1.1673464d7,
     : 3.2350933d7,0d0,0d0,0d0,0d0,0d0,-1.6888448d7,7.7833319d9,0d0,
     : 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, -1.3859509d7, 0d0, 0d0, 0d0,
     : -2.8745454d5, -5.0746248d5, -8.1294776d5, -6.1239403d7,
     : -9.5031834d7, -3.2712375d6, 2.1530251d9, 0d0, -2.4554310d5,
     : -2.2240816d5, -3.2010603d5, -4.6797756d5, -5.1348391d5,
     : -1.0761377d6, -2.7541840d6, 1.2509443d5, 1.0703694d5,
     : 9.2097995d4, 7.7029394d4, 6.7256907d4, -5.7548141d4,
     : -1.6037808d4, 1.0363696d8, -2.4706938d5, 1.3091573d5,
     : 2.0175792d5, 2.0199954d5, 2.0236672d5, 2.0379734d5, 2.1362737d5,
     : 1.1499859d5, 5.0489318d7, -2.0675897d4, 2.4254857d5, 2.1515005d5,
     : 2.4077851d5, 2.6794892d5, 4.5326857d5, 4.9394601d5, 3.7831220d5,
     : 6.8189939d5, 2.8892658d5, 2.7494329d5, 2.8732767d5, 3.1714555d5,
     : 3.5210836d5, 3.8278795d5, 6.5391444d5, 7.9210598d5/

       data xlam/

     : 3.2563258d0, 3.2072775d0, 3.2106318d0, 3.2149238d0,
     : 3.1882921d0, 2.6473556d0, 4.0056990d0, 9.5292984d-1, 2.0072677d0,
     : 1.9065953d0, 3.3169172d0, 3.3022697d0, 3.2950305d0, 3.2984058d0,
     : 3.1929618d0, 2.5641135d0, 1.1089247d0, 1.8811506d0, 1.7194055d0,
     : 2.0391202d0, 2.1021763d0, 2.1587227d0, 3.4171862d0, 3.5400259d0,
     : 2.3330267d0, 3.9048065d0, 1.3255981d0, 2.0021043d0, 1.9635080d0,
     : 1.9976766d0, 2.0410133d0, 2.0302936d0, 2.3308629d0, 2.4278447d0,
     : 4.8618361d0, 4.8488676d0, 4.8682629d0, 4.8529103d0, 4.5719712d0,
     : 5.0623483d0, 3.1822524d0, 4.9876677d0, 4.0714388d0, 4.1116751d0,
     : 4.7975329d0, 4.8056888d0, 4.8108758d0, 4.8263656d0, 4.5778508d0,
     : 6.3028159d0, 4.8271622d0, 4.1452516d0, 3.9281787d0, 6.3552827d0,
     : 6.2712191d0, 6.1914432d0, 4.7539899d0, 4.7365021d0, 6.0171260d0,
     : 3.5402289d0, 3.9046544d0, 6.3157114d0, 6.4179961d0, 6.3722076d0,
     : 6.3009140d0, 6.3523282d0, 2.0223063d0, 3.3968837d0/
*  this is matrix of dl0(theta) at 0, 30, 60 and 90 degrees
*  l=0 first row, l=2 2nd row, l=4 3rd row, l=6 4th row
       data d0mat/
     : 1d0, 1d0, 1d0, 1d0, 1d0, 6.25d-1, -1.25d-1, -5.d-1, 1d0,
     : 2.3437511d-2, -2.8906251d-1, 3.75d-1, 1d0, -3.7402343d-1,
     : 3.2324218d-1, -3.125d-1/
* this version uses all values of r, and fits with cubic (least squares)
* determine potential at theta=30 and 60 and r values needed for interpolation
* determine potential at theta=0, 30, 60, and 90  and r values needed for inte
* determine potential at theta=30 and 60 and r values needed for interpolation

      do ith=1,4
        if (ith .lt.3) then
           icol0=(ith-1)*9
           itop=9
        else
           icol0=(ith-1)*8+2
           itop=8
        endif
        do ir=1,itop
           icol=icol0+ir
           call dcopy(2,xlam(icol,1),34,xxlam,1)
           call dcopy(5,cc(icol,1),34,ccc,1)
           v(ir,ith)=voneexps(r,alpha,xxlam,ccc,ipow)
           if (flagd) then
             dv(ir,ith)=dvoneexps(r,alpha,xxlam,ccc,ipow)
           endif
        enddo
      enddo
* merge fits for rh2=2.3 and theta=0,30
      on=half*(tanh(5d0*(r-4.8d0))+one)
      if (flagd) then
* derivative of switching function
        onp=2.5d0*(one-(tanh(5d0*(r-4.8d0)))**2)
      endif
      do i=1,2
         v(8,i)=on*v(9,i)+(one-on)*v(8,i)
         if (flagd) then
           dv(8,i)=on*dv(9,i)+(one-on)*dv(8,i)+onp*(v(9,i)-v(8,i))
         endif
      enddo
      re=1.4d0
      do ir=1,8
         cmat(ir,1)=one
         rmre=rr(ir)-re
         cmat(ir,2)=rmre
         cmat(ir,3)=rmre*cmat(ir,2)
         cmat(ir,4)=rmre*cmat(ir,3)
         cmat(ir,5)=rmre*cmat(ir,4)
      enddo
      call dqrank(cmat,8,8,5,tol,kr,kpvt,scr1,scr2)
      do ith=1,4
* determine r expansion coefficients (quartic fit)
        if (ith .lt.3) then
           icol0=(ith-1)*9+1
        else
           icol0=(ith-1)*8+2
        endif
        call dqrlss(cmat,8,8,5,kr,v(1,ith),xth,scr2,kpvt,scr1)
        if (flagd) then
          call dqrlss(cmat,8,8,5,kr,dv(1,ith),dxth,scr2,kpvt,scr1)
        endif

* calculate large rh2 damping
        damp=half*(one-tanh(alphlr*(rh2-rlongr)))
        call dcopy(5,xth,1,xtilde,1)
* xtilde will hold original undamped coefficients
        call dscal(4,damp,xth(2),1)
* new april 2
        if (flagd) call dscal(4,damp,dxth(2),1)
* evaluate derivative of damping function and multiply original
* expansion coefficients by this
        if (flagd) then
          ddamp=-half*alphlr*(one-(tanh(alphlr*(rh2-rlongr)))**2)
          call dscal(4,ddamp,xtilde(2),1)
        endif

* interpolate in r (quartic)
        vv=xth(1)
        if (flagd) dvvdr=dxth(1)
        rmre=rh2-re
        rterm=rmre
        do ir=2,5
          vv=vv+xth(ir)*rterm
          if (flagd) then
            dvvdr=dvvdr+dxth(ir)*rterm
          endif
          rterm=rterm*rmre
        enddo
        vang(ith)=vv
        if (flagd) then
          dvangdr(ith)=dvvdr
        endif
* calculated derivative in r (if desired)
* f(rmre)=x(1)+x(2)*rmre+x(3)*rmre^2+x(4)*rmre^3
* df/drmre=x(2)+2*x(3)*rmre+3*x(4)*rmre^2+
*       x'(1)+x'(2)*rmre+x'(3)*rmre^2+x'(4)*rmre^3
* where x(i)=xtilde(i)*damp
* and   x'(i)=xtilde(i)*ddamp
* NB at this point the array x contains xtilde*damp
*  and the array xtilde contains xtilde*ddamp
       if (flagd) then
          rterm=rmre
          vv=xth(2)
          vvtilde=xtilde(2)*rterm
          do ir=3,5
            vv=vv+xth(ir)*rterm*(ir-1)
            vvtilde=vvtilde+xtilde(ir)*rterm*rmre
            rterm=rterm*rmre
          enddo
* vv now contains derivative of x with respect to r including
* damping function
          dvang(ith)=vv+vvtilde
        endif
      enddo
* calculate  coefficients in legendre expansion
       call dcopy(16,d0mat,1,d0matx,1)
       call dqrank(d0matx,4,4,4,tol,kr,kpvt,scr1,scr2)
       call dqrlss(d0matx,4,4,4,kr,vang,xang,scr2,kpvt,scr1)
       if (flagd) then
         call dqrlss(d0matx,4,4,4,kr,dvang,dxang,scr2,kpvt,scr1)
         call dqrlss(d0matx,4,4,4,kr,dvangdr,dxangdr,scr2,kpvt,scr1)
       endif
       x=cos(theta*pi/180d0)
       x2=x*x
       x4=x2*x2
       x6=x2*x4
       pleg(1)=one
       pleg(2)=(3d0*x2-one)*half
       pleg(3)=(3.d0-30d0*x2+35*x4)/8d0
       pleg(4)=(-5d0+105d0*x2-315d0*x4+231d0*x6)/16d0
       vv=ddot(4,xang,1,pleg,1)
       if (flagd) then
         dvv(1)=ddot(4,dxangdr,1,pleg,1)
         dvv(2)=ddot(4,dxang,1,pleg,1)
       endif
       if (flagd) then
         sn=sin(theta*pi/180d0)
         dpleg(1)=zero
         dpleg(2)=-3d0*x*sn
         dpleg(3)=(2.5d0)*x*sn*(3d0-7d0*x2)
         dpleg(4)=-(2.625d0)*x*sn*(5d0-30d0*x2+33d0*x4)
         dvv(3)=ddot(4,xang,1,dpleg,1)
       endif
       if (.not. flagd) then
         do i=1,3
           dvv(i)=0d0
         enddo
       endif
       return
       end
* ------------------------------------------------------
       double precision function dvoneexps(r,alph,xlam,c,ipow)

*  to determine
*     dv(r)/dr, where
*     v(r)=exp(-xlam1*r)*(c1+c2*r+c3*r^2+c4*r^3)
*          -0.5*c5*(tanh(alph*(r-xlam2))+1)*r^(-ipow)
       implicit double precision (a-h,o-z)
       dimension xlam(2),c(5)
       data zero,one,half/0.d0,1.d0,0.5d0/

       clr=c(5)
*      if (c(4).ne.0d0)  then
*        clr=c(4)
*        c(4)=c(5)
*      endif
       poly=(((c(4)*r+c(3))*r+c(2))*r+c(1))
       ex=exp(-xlam(1)*r)
       poly1=((3d0*c(4)*r+2*c(3))*r+c(2))
       xlr=-half*clr*(tanh(alph*(r-xlam(2)))+one)*r**(-ipow)
       dxlr=ipow*(tanh(alph*(r-xlam(2)))+one)/r
     :   -alph*(one-(tanh(alph*(r-xlam(2))))**2)
       dxlr=dxlr*half*clr*r**(-ipow)
       dvoneexps=ex*(-xlam(1)*poly+poly1)+dxlr
       return
       end
* ------------------------------------------------------
       double precision function voneexps(r,alph,xlam,c,ipow)

*  to determine
*     v(r)=exp(-xlam1*r)*(c1+c2*r+c3*r^2+c4*r^3)
*          -0.5*c5*(tanh(alph*(r-xlam2))+1)*r^(-ipow)
       implicit double precision (a-h,o-z)
       dimension xlam(2),c(5)
       data zero,one,half/0.d0,1.d0,0.5d0/

       clr=c(5)
*      if (c(4).ne.0d0)  then
*        clr=c(4)
*        c(4)=c(5)
*      endif
       poly=(((c(4)*r+c(3))*r+c(2))*r+c(1))
       ex=exp(-xlam(1)*r)
       xlr=-half*clr*(tanh(alph*(r-xlam(2)))+one)*r**(-ipow)
       voneexps=ex*poly+xlr
       return
       end
* ------------------------------------------------------
      double precision function vh2(t)
* h2 potential from scaled ci calculation at R=30 (s=1.038)
      implicit double precision (a-h,o-z)
      dimension xlam(2),cc(6)
      data xlam /2.1192d0, 1.0904d0/
* old data statement
*     data cc /
*    : 6.4178d5, 3.8388d5, 9.5175d4, -6.8014d5, 2.3454d5, -2.1883d4/
* correct so that minimum is -De
*     del=-4.858;
*     x=10.4;
*     cc(1)=cc(1)+x*del;
*     cc(4)=cc(4)+(1-x)*del;
      data cc/
     : 6.417294768d5, 3.8388d5, 9.5175d4, -6.800943348d5,
     : 2.3454d5, -2.1883d4/

*  re= 1.4003;  de= 38364.86
* actual re and de of fitted potential given here is
*  re= 1.3998;  de= 38364.87
*  TITLE  RE(A)    BE     AE      WE     WEXE   WEYE   DE    D0 SIG
*  Exp.  0.7414 60.853  3.0622  4401.21 121.34   .81  4.747 4.47813
*  fit   0.7410 60.923  3.165   4413.4  124.4    2.0  4.757 4.4858


*  t=r-re
      vh2=(cc(1)+cc(2)*t+cc(3)*t*t)*exp(-xlam(1)*t)
     : + (cc(4)+cc(5)*t+cc(6)*t*t)*exp(-xlam(2)*t)
      return
      end
* ------------------------------------------------------
      double precision function dvh2(t)
* derivative with respect to r of h2 potential from scaled
* ci calculation at R=30 (s=1.038)
      implicit double precision (a-h,o-z)
      dimension xlam(2),cc(6)
      data xlam /2.1192d0, 1.0904d0/
* old data statement
*     data cc /
*    : 6.4178d5, 3.8388d5, 9.5175d4, -6.8014d5, 2.3454d5, -2.1883d4/
* correct so that minimum is -De
*     del=-4.858;
*     x=10.4;
*     cc(1)=cc(1)+x*del;
*     cc(4)=cc(4)+(1-x)*del;
      data two /2d0/
      data cc/
     : 6.417294768d5, 3.8388d5, 9.5175d4, -6.800943348d5,
     : 2.3454d5, -2.1883d4/

*  re= 1.4003;  de= 38364.86
* actual re and de of fitted potential given here is
*  re= 1.3998;  de= 38364.87
*  TITLE  RE(A)    BE     AE      WE     WEXE   WEYE   DE    D0 SIG
*  Exp.  0.7414 60.853  3.0622  4401.21 121.34   .81  4.747 4.47813
*  fit   0.7410 60.923  3.165   4413.4  124.4    2.0  4.757 4.4858


*  t=r-re
      dvh2=(-xlam(1)*(cc(1)+cc(2)*t+cc(3)*t*t)+
     :    cc(2)+two*cc(3)*t)*exp(-xlam(1)*t)
     : +(-xlam(2)*(cc(4)+cc(5)*t+cc(6)*t*t)+
     :    cc(5)+two*cc(6)*t)*exp(-xlam(2)*t)
      return
      end
* ------------------------------------------------------
       subroutine aso(r,rh2,a,b,da,db,flagd)
* --------------------------------------------------------------
* subroutine to determine F-H2 spin-orbit couplings from
* SEC with scale factor of 1.038
* fitted by m.h. alexander as follows:
* the r and R dependence is fitted as
* note that a is is pix-piy coupling (what i call A)
* and b is pix-sigma coupling (what i call B)

*     v(R,r)=exp(-xlam1*R)*(c1+c2*R+c3*R^2)
*           +(r-re)*exp(-xlam2*R)*(c4+c5*R+c6*R^2)
*           +(r-re)^2*exp(-xlam3*R)*(c7+c8*R+c9*R^2)
* involving 9 linear and 3 non-linear parameters
* for a given value of R=Ri and r=ri, the v(Ri,ri)'s are determined
* from the above equation, using parameters given below.
* errors are
*  mean error in a is 0.317 cm-1 for all points
*  mean error in b is 0.058 cm-1 for all points
*  max error in a is 1.71 cm-1
*  max error in b is 0.189 cm-1
* at small R and/or large r the pi-sigma spin-orbit coupling is damped
* as follows
*       dampR=0.5*(tanh(3.5*(R-2))+1)
*       damprh2=0.5.*(tanh(4.5*(2.4-rh2))+1)
*       aso=aso.*dampr*damprh2+ainf*dampr
* at small R and/or large r the pi-pi spin-orbit coupling is damped
* as follows
*       dampR=0.5*(tanh(2*(R-2.2))+1)
*       damprh2=0.5.*(tanh(4.5*(2.4-rh2))+1)
*      vso=vso.*dampr*damprh2+(1-dampr)*damprh2*2*v0+ainf
*  where v0 is vso(R=2.7,r)

* if flagd=.true., then derivatives are calculated
*  and returned in the vectors da and db with
*   da(1)=da/dR
*   da(2)=da/dr
*   da(3)=da/dtheta=0
* and similarly for db
* --------------------------------------------------------------
* current revision date:  13-apr-1998
* --------------------------------------------------------------
       implicit double precision (a-h,o-z)
       dimension xlama(6),xlamb(2),cca(3),ccb(4),xxlam(6),ccc(6),
     :           daconst(2),da(3),dbconst(2),db(3),dbconst0(2),
     :           dtemp(3)
       logical flagd
       data one,half,two /1.d0,0.5d0,2.d0/
       data pi /3.141592653589793d0/
* equilibrium h2 distance for potential expansion
       data re /1.4d0/
* asymptotic so constant
       data ainf/ 132.41d0/
       data ione /1/
       data xlama /1.5177d0, 2.8499d0, 1.7330d0, 3.3410d0, 3.0997d0,
     :             3.8876d0/
       data xlamb /1.9942d0, 2.5334e+00/
       data cca /2.9099d1, 3.1968d1, 1.2887d1/
       data ccb /-1.6802d+3, 5.0964d+2, 1.1645d+4, -6.3387d+3/
* determine spin orbit constants
       call dcopy(6,xlama,1,xxlam,1)
       call dcopy(3,cca,1,ccc,1)
       call dasotanh(r,rh2,re,xxlam,ccc,aconst,daconst,nlam,flagd)
       nlam=2
       call dcopy(nlam,xlamb,1,xxlam,1)
       call dcopy(3*nlam,ccb,1,ccc,1)
       call dasoexp(r,rh2,re,xxlam,ccc,bconst,dbconst,nlam,flagd)
* damp pi-sig so constant to zero
       dr=r-2d0
       drh2=2.4-rh2
       dampr=half*(tanh(3.5d0*dr)+one)
       damprh2=half*(tanh(4.5d0*drh2)+one)
       a=aconst*dampr*damprh2+dampr*ainf
       if (flagd) then
          ddampr=half*3.5d0*(one-tanh(3.5d0*dr)**2)
          ddamprh2=-half*4.5d0*(one-tanh(4.5d0*drh2)**2)
          da(1)=aconst*damprh2*ddampr+daconst(1)*dampr*damprh2+ainf*ddampr
          da(2)=aconst*ddamprh2*dampr+daconst(2)*dampr*damprh2
       endif
* damp pi-pi so constant
       dr=r-2.2d0
       dampr=half*(tanh(two*dr)+one)
* evaluate pi-pi constant at 2.7
       call dasoexp(2.7d0,rh2,re,xxlam,ccc,
     :              bconst0,dbconst0,nlam,.false.)
       b=bconst*dampr*damprh2+(one-dampr)*damprh2*two*bconst0+ainf
       if (flagd) then
          ddampr=half*two*(1-tanh(two*dr)**2)
          db(1)=bconst*ddampr*damprh2+dbconst(1)*dampr*damprh2
     :           -ddampr*damprh2*two*bconst0
          db(2)=bconst*dampr*ddamprh2+dbconst(2)*dampr*damprh2
     :           +(one-dampr)*ddamprh2*two*bconst0
       endif
       da(3)=zero
       db(3)=zero

* switch definition of a and b to be consistent with my later
* notation:  A = pix-piy coupling and B = pi-sigma coupling
       temp=a
       call dcopy(3,da,1,dtemp,1)
       a=b
       b=temp
       call dcopy(3,db,1,da,1)
       call dcopy(3,dtemp,1,db,1)

       return
       end
* ------------------------------------------------------
       subroutine dasoexp(r,rh2,re,xlam,c,v,dv,nlam,flagd)
*  to determine aso(R,r), daso(R,r)/dr and daso(R,r)/dR, where
*     asp(R,r)=exp(-xlam1*R)*(c1+c2*R+c3*R^2)
*           +(r-re)*exp(-xlam2*R)*(c4+c5*R+c6*R^2)
* and, if nlam=3
*           +(r-re)^2*exp(-xlam3*R)*(c7+c8*R+c9*R^2)
* on return:  dv(1) is daso/dR and dv(2) is daso/dr
* derivatives determined only if flagd=.true.
       implicit double precision (a-h,o-z)
       logical flagd
       dimension xlam(3),c(9),dv(2),ex(4)
       data one,half,two/1.d0,0.5d0,2d0/
       rmre=rh2-re
       rmre2=rmre*rmre
       r2=r*r
       vex=zero
       dvexr=zero
       dvexrh2=zero
       rho=one
       rho2=one
       do i=1,nlam
         exx=exp(-xlam(i)*r)
         ind=(i-1)*2+1
         vexterm=exx*(c(ind)+c(ind+1)*r)
         vex=vex+rho*vexterm
         if (flagd) then
           dvexrterm=rho*exx*(c(ind+1)
     :           -xlam(i)*(c(ind)+c(ind+1)*r))
           dvexr=dvexr+dvexrterm
           if (i.gt.1) then
             dvexrh2=dvexrh2+(i-1)*rho2*vexterm
             rho2=rho2*rmre
           endif
         endif
         rho=rho*rmre
       enddo
       v=vex
       dv(1)=dvexr
       dv(2)=dvexrh2

       return
       end
* ------------------------------------------------------
       subroutine dasotanh(r,rh2,re,xlam,c,v,dv,nlam,flagd)
*  to determine aso(R,r), daso(R,r)/dr and daso(R,r)/dR, where
*     aso(R,r)=c(1)*(tanh(-xlam1*(R-xlam2))-1)
*           +(r-re)*c(2)*(tanh(-xlam3*(R-xlam4))-1)
*           +(r-re)^2*c(3)*(tanh(-xlam5*(R-xlam6))-1)
* on return:  dv(1) is daso/dR and dv(2) is daso/dr
* derivatives determined only if flagd=.true.
       implicit double precision (a-h,o-z)
       logical flagd
       dimension xlam(6),c(3),dv(2)
       data one,half,two/1.d0,0.5d0,2d0/
* asymptotic spin-orbit constant
       data ainf / 132.41d0/
       rmre=rh2-re
       rmre2=rmre*rmre
       r2=r*r
       vex=zero
       dvexr=zero
       dvexrh2=zero
       rho=one
       rho2=one
       do i=1,3
         ind=2*(i-1)+1
         exx=tanh(xlam(ind)*(r-xlam(ind+1)))
         vexterm=(exx-one)*c(i)
         vex=vex+rho*vexterm
         if (flagd) then
           dvexr=dvexr+rho*c(i)*xlam(ind)*(one-exx*exx)
           if (i.gt.1) then
             dvexrh2=dvexrh2+(i-1)*rho2*vexterm
             rho2=rho2*rmre
           endif
         endif
         rho=rho*rmre
       enddo
       v=vex
       dv(1)=dvexr
       dv(2)=dvexrh2

       return
       end
* --------------------------
      subroutine fh2sw(r,vau)
      implicit double precision (a-h,o-z)
      integer b,c,d
c
c     -----------------------------------------------------------------
c     Stark-Werner F+H2 potential energy surface
c
c     r(1) = H-H distance in bohr
c     r(2) = F-H distance in bohr
c     r(3) = H-F distance in bohr
c
c     vau = potential in hartree from bottom of asymptotic F+H2 valley
c     -----------------------------------------------------------------
c
      dimension r(3)
      parameter (nmx=200)
      dimension a(nmx),b(nmx),c(nmx),d(nmx),p(nmx)
      data icall/0/
      save icall,a,b,c,d,p,ma,np
c
      if (icall .eq. 0) then
         open (unit=1,file='potdata/three.param',
     :     status='old')
         i=1
   2     read (1,*,end=3) nparm,b(i),c(i),d(i),a(i)
            i=i+1
            goto 2
   3     ma=i-1
         open (unit=1,file='potdata/two.param',
    :      status='old')
         i=1
   4     read (1,*,end=5) p(i)
            i=i+1
            goto 4
   5     np=i-1
      endif
      icall = icall+1
c
c     (Modified) Aguado-Paniagua type functions:
c
      x = min(r(2),r(3))
      y = r(1)
      z = max(r(2),r(3))
      b1 = a(ma-5)
      b2 = a(ma-4)
      b3 = a(ma-3)
      x0 = a(ma-2)
      y0 = a(ma-1)
      z0 = a(ma)
      fit = 0.0d0
      do i = 1,ma-6
         expon = b(i)*b1*(x-x0)+c(i)*b2*(y-y0)+d(i)*b3*(z-z0)
         fex = exp(-expon)
         fxy = ((x**b(i))*(y**c(i))*(z**d(i)))*fex
         fit = fit+a(i)*fxy
      enddo
      xr=x-p(3)
      yr=y-p(9)
      zr=z-p(3)
      fx=dexp(-p(2)*xr)
      fy=dexp(-p(8)*yr)
      fz=dexp(-p(2)*zr)
      xval=-p(1)*(1.0d0+p(2)*xr+p(4)*xr**2+p(5)*xr**3)*fx+p(6)
      yval=-p(7)*(1.0d0+p(8)*yr+p(10)*yr**2+p(11)*yr**3)*fy+p(12)
      zval=-p(1)*(1.0d0+p(2)*zr+p(4)*zr**2+p(5)*zr**3)*fz+p(6)
      vagpan=fit+xval+yval+zval
      vau = vagpan
      return
      end
* -------------------------------
      subroutine prepot
c
c     ab initio p.e.s. for f+h2 built by stark and werner
c     modifications by l. banares and v. s ez r banos october 1992
c
      implicit real*8 (a-h,o-z),integer*4(i-n)
      integer*4 b(198),c(198),d(198)
      dimension nparm(198),a(198),p(12)
      dimension r(3),dpe(3),rb(3),dpeb(3)
c
c..... read the three-body-parameters
c..... (computed with the levenberg-marquardt-procedure)
c..... nparm is the actual number of a parameter
c
      open(3,file='potdata/three.param')
      i=1
   10 read(3,*,end=11) nparm(i),b(i),c(i),d(i),a(i)
      i=i+1
      goto 10
c
   11 ma=i-1
      close(3)
c
c.... read the two-body-parameters p
c.... (computed with an extended-rydberg-fit)
c
      open(1,file='potdata/two.param')
      i=1
   12 read(1,*,end=111) p(i)
      i=i+1
      goto 12
  111 np=i-1
      close(1)
c
c.... initialize the non-linear parameters
c
      return
      entry fh2sw_deriv(rb,pe,dpeb)
* convert to same bond ordering convention used in
* this subroutine
      r(1)=rb(2)
      r(2)=rb(1)
      r(3)=rb(3)
c
      b1 = a(ma-5)
      b2 = a(ma-4)
      b3 = a(ma-3)
      x0 = a(ma-2)
      y0 = a(ma-1)
      z0 = a(ma)
      fit = 0.0d0
      dfdx = 0.0d0
      dfdy = 0.0d0
      dfdz = 0.0d0
c
c... the (modified) aguado-paniagua type functions
c... for the three-body-potential
c
      xexpon = b1*(r(1)-x0)
      yexpon = b2*(r(2)-y0)
      zexpon = b3*(r(3)-z0)

           exponx=dexp(-xexpon)
           expony=dexp(-yexpon)
           exponz=dexp(-zexpon)

      fex = r(1)*exponx
      fey = r(2)*expony
      fez = r(3)*exponz


          drhox = exponx*(1-r(1)*b1)
          drhoy = expony*(1-r(2)*b2)
          drhoz = exponz*(1-r(3)*b3)



      do 1010 i=1,ma-6

        pepex = fex
        pepey = fey
        pepez = fez

             if(b(i).eq.0) then
              pepex=1.d0
              go to 1018
             endif
             do 1020 j=1,b(i)-1
             pepex = pepex*fex
 1020        continue
 1018        if(c(i).eq.0) then
              pepey=1.d0
              go to 1019
             endif

             do 1021 j=1,c(i)-1
             pepey = pepey*fey
 1021        continue
 1019        if(d(i).eq.0) then
              pepez=1.d0
              go to 1024
             endif

            do 1022 j=1,d(i)-1
             pepez = pepez*fez
 1022        continue



 1024    fxy = pepex*pepey*pepez
         fit = fit + a(i)*fxy

          dftx = drhox*b(i)*pepex/fex
          dfty = drhoy*c(i)*pepey/fey
          dftz = drhoz*d(i)*pepez/fez

          dftx = dftx*pepey*pepez
          dfty = dfty*pepex*pepez
          dftz = dftz*pepex*pepey

          dfdx = dfdx+a(i)*dftx
          dfdy = dfdy+a(i)*dfty
          dfdz = dfdz+a(i)*dftz
c

 1010 continue
c
c.... two-body-potential : extended-rydberg-functional
c
      xr = r(1)-p(3)
      yr = r(2)-p(9)
      zr = r(3)-p(3)
        xr2=xr*xr
        xr3=xr2*xr
        yr2=yr*yr
        yr3=yr2*yr
        zr2=zr*zr
        zr3=zr2*zr

      fx = dexp(-p(2)*xr)
      fy = dexp(-p(8)*yr)
      fz = dexp(-p(2)*zr)

      ux = -p(1)*(1.0d0+p(2)*xr+p(4)*xr2+p(5)*xr3)
      uy = -p(7)*(1.0d0+p(8)*yr+p(10)*yr2+p(11)*yr3)
      uz = -p(1)*(1.0d0+p(2)*zr+p(4)*zr2+p(5)*zr3)

      xval=ux*fx+p(6)
      yval=uy*fy+p(12)
      zval=uz*fz+p(6)

c

      dux = -p(1)*(p(2)+p(4)*2.d0*xr+p(5)*3.d0*xr2)
      duy = -p(7)*(p(8)+p(10)*2.d0*yr+p(11)*3.d0*yr2)
      duz = -p(1)*(p(2)+p(4)*2.d0*zr+p(5)*3.d0*zr2)
c
      dfx = -p(2)*fx
      dfy = -p(8)*fy
      dfz = -p(2)*fz
c
      dx = dux*fx+ux*dfx
      dy = duy*fy+uy*dfy
      dz = duz*fz+uz*dfz
c
c.... resulting energy in au
c
      dpe(1) = dfdx+dx
      dpe(2) = dfdy+dy
      dpe(3) = dfdz+dz
      dpeb(1)=dpe(2)
      dpeb(2)=dpe(1)
      dpeb(3)=dpe(3)
c
      xe = fit+xval+yval+zval
      pe = fit+xval+yval+zval
c
      return
      end
