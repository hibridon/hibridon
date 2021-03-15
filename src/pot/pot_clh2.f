*system:  Cl(2P)+H2, Dubernet-Hutson expansion of capecchi-werner PES's
*references: G. Capecchi, H.-J. Werner Phys. Chem. Chem. Phys. 6, 4975 (2004)
* M.-L. Dubernet and J. M. Hutson, J. Chem. Phys. 101, 1939 (1994).
*
*  Note:  this pot routine requires data files to be in hibxx/bin/progs/potdata:
*         three.param4, cwfit.dat, BW.3p
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

      potnam='CAPECCHI-WERNER Cl(2P)H2 DUBERNET-HUTSON'
      ibasty=12
      lammin(1)=1
      lammax(1)=9
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
      character *2 frame
      logical csflag, ljunk, ihomo, lljunk
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(9)
      potnam='CAPECCHI-WERNER Cl(2P)H2 DUBERNET-HUTSON'
      print *, potnam
      print *
1      print *, ' r (bohr) frame (sf/bf)'
      read (5, *, end=93) r, frame
      if (frame.eq.'sf') then
        csflag=.false.
        ihomo=.true.
      else if (frame .eq. 'bf') then
        csflag=.true.
        ihomo=.false.
      else
        print *, 'frame must be either "sf" or "bf"'
        go to 1
      endif
      call pot(vv0,r)
      if (.not. csflag .or. (csflag .and. ihomo)) write (6, 100) vv0,vvl
100   format(' V000, V220, V022, V202:  ',4(1pe16.8),/,
     :       ' V222, V224, V404:  ',3(1pe16.8),/,
     :       ' V422, V424, V426:  ',3(1pe16.8))
      if (csflag .and. .not.ihomo) write (6, 110) vv0,vvl
110   format(' v000, v220, v020, v200:  ',4(1pe16.8),/,
     :       ' v222, v221, v400:  ',3(1pe16.8),/,
     :       ' v420, v422, v421:  ',3(1pe16.8))
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
*  f-h2 potentials of alexander,stark,werner at rh2=1.4  using body-frame
*  expansion of dubernet and flower
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*  vv0       totally symmetric potential (v000)
*  variable in common block /covvl/
*    vvl:     vector of length 5 to store r-dependence of each term
*             in potential expansion
*  CS calculations case 1A (csflag=.true. and ihomo=.false)
*  body frame expansion coefficients
*    vvl(1)  v220(R)
*    vvl(2)  v020(R)
*    vvl(3)  v200(R)
*    vvl(4)  v222(R)
*    vvl(5)  v221(R)
*    vvl(6)  v400(R)
*    vvl(7)  v420(R)
*    vvl(8)  v422(R)
*    vvl(9)  v421(R)
*  CC calculations (csflag=.false.)
*    vvl(1)  V220(R)
*    vvl(2)  V022(R)
*    vvl(3)  V202(R)
*    vvl(4)  V222(R)
*    vvl(5)  V224(R)
*    vvl(6)  V404(R)
*    vvl(7)  V422(R)
*    vvl(8)  V424(R)
*    vvl(9)  V426(R)

* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  3-mar-2001
* ----------------------------------------------------------------------


      implicit double precision (a-h,o-z)
      logical ljunk, ihomo, csflag, lljunk, flagd
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(9)
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

* determine potentials at angles
      do 200 i=1,4
        theta=(4-i)*30
        rh2=1.4d0
        flagd=.false.
        call potsub_clh2(r,rh2,theta,hmat)
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
        print *, theta, vs(i), vzz(i), vd(i), v12(i)
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
      call dqrlss(aa,4,4,4,kr,vs,xs,rsd,kpvt,qraux)
      call dcopy(4,xs,1,vvll(5),1)
      call dcopy(9,d2,1,aa,1)
      call dqrank(aa,3,3,3,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,3,3,3,kr,vd,xd,rsd,kpvt,qraux)
      call dcopy(3,xd,1,vvll(9),1)
      call dcopy(4,d1,1,aa,1)
      call dqrank(aa,2,2,2,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,2,2,2,kr,v12(2),x1,rsd,kpvt,qraux)
      call dcopy(2,x1,1,vvll(12),1)
* determine body frame expansion coefficients in dubernet and flower
* expansion
* here is totally symmetric term
      vv0=(two*vvll(5)+vvll(1))/three
      sq6=1.d0/sqrt(6.d0)
* v200
      vvl(1)=(two*vvll(6)+vvll(2))/three
* v400
      vvl(6)=(two*vvll(7)+vvll(3))/three
* v020
      vvl(2)=(vvll(1)-vvll(5))*five/three
* v220
      vvl(3)=(vvll(2)-vvll(6))*five/three
* v420
      vvl(7)=(vvll(3)-vvll(7))*five/three
* v222
      vvl(4)=five*vvll(9)*sq6
* v422
      vvl(9)=five*vvll(10)*sq6
* v221
      vvl(5)=five*vvll(12)*sq6
* v421
      vvl(8)=five*vvll(13)*sq6
* transform to space frame terms if CC or CS case 1C
      if (.not. csflag .or. (csflag .and. ihomo)) then
        call dcopy(9,vvl,1,vvll,1)
        vvl(1)=(two*(vvll(4)-vvll(5))+vvll(3))/sqrt(5.d0)
        vvl(2)=vvll(2)
        vvl(3)=vvll(1)
        vvl(4)=(two*vvll(4)+vvll(5)-vvll(3))*sqrt(two/7.d0)
        vvl(5)=(vvll(4)+four*vvll(5)+3*vvll(3))*sqrt(two/35.d0)
        vvl(6)=vvll(6)
        vvl(7)=sqrt(two/seven)*vvll(7)-sqrt(20.d0/21.d0)*vvl(8)
     :        +sqrt(10d0/21d0)*vvl(9)
        vvl(8)=sqrt(one/77d0)*(-sqrt(20d0)*vvll(7)+sqrt(6d0)*vvl(8)
     :        +sqrt(108d0)*vvl(9))
        vvl(9)=sqrt(five/eleven)*vvll(7)+sqrt(32.d0/33.d0)*vvl(8)
     :        +sqrt(four/33d0)*vvl(9)
      else
* reorder body-frame terms
        vvll(1)=vvl(1)
        vvl(1)=vvl(3)
        vvl(3)=vvll(1)
        vvll(8)=vvl(9)
        vvl(9)=vvl(8)
        vvl(8)=vvll(8)
      endif
* convert to hartree
      econv=1./219474.6
      call dscal(9,econv,vvl,1)
      vv0=vv0*econv
      end
* --------------------------------------------------------------
       subroutine potsub_clh2(rr,r,theta,vcm)
       implicit double precision (a-h,o-z)
* calculates lowest diabatic cl+h2 surfaces and spin-orbit coupling
* from fit to capecchi-bian-werner surfaces including
* spin-orbit coupling
* input variable (dimension 3) (distances in bohr)
*      rbond(1)=rh2
*      rbond(2)=rclh1
*      rbond(3)=rclh2
* subroutine returns (energies in cm-1)
*      vcm(1)=vsig=V_z
*      vcm(2)=vpi=(V_y + V_x)/2
*      vcm(3)=vd=(V_y - V_x)/2
*      vcm(4)=v1=V_xz/sqrt(2)
*      vcm(5)=aconst
*      vcm(6)=bconst
      data toev/27.21139610d0/
      data tocm/219474.63067d0/
      data sq2i/0.707106781d0/

      dimension  vcm(6)
c
c
      v1_cw  =tgcn(rr,r,theta,4)*damp(rr,r,theta,4)*sq2i
c
c     aconst =tgcn(rr,r,theta,5)
c     bconst =tgcn(rr,r,theta,6)
c
      if(r.le.1.d0) then
         aconst=tgcn(rr,1.d0,theta,5)*damp(rr,r,theta,5)
      elseif(r.gt.2.5d0) then
         aconst=tgcn(rr,2.5d0,theta,5)*damp(rr,r,theta,5)
      else
         aconst =tgcn(rr,r,theta,5)*damp(rr,r,theta,5)
      endif
      if(r.le.1.d0) then
         bconst =tgcn(rr,1.d0,theta,6)*damp(rr,r,theta,6)
      elseif(r.gt.2.5d0) then
         bconst=tgcn(rr,2.5d0,theta,6)*damp(rr,r,theta,6)
c     bconst =tgcn(rr,r,theta,6)
      else
         bconst =tgcn(rr,r,theta,6)*damp(rr,r,theta,6)
      endif
c
*  determine bond coordinates

      rh2=r
      cs=cos(theta*3.1415965358979d0/180d0)
      rclh1=sqrt(r*r+rr*rr-2*cs*r*rr)
      rclh2=sqrt(r*r+rr*rr+2*cs*r*rr)

      vsig=dmerge(rclh1,rh2,rclh2,rr,theta)
      vpi=tgcn(rr,r,theta,2)
c     if(r.gt.2.5d0) then
c       vd=tgcn(rr,2.5d0,theta,3)
c     else if(r.lt.1.0d0) then
c       vd=tgcn(rr,1.0d0,theta,3)
c     else
c       vd =tgcn(rr,r,theta,3)
c     end if
c
      vd =tgcn(rr,r,theta,3)*damp(rr,r,theta,3)
c
      vcm(1)=vsig*tocm
      vcm(2)=vpi*tocm
      vcm(3)=vd*tocm
      vcm(4)=v1_cw*tocm
      vcm(5)=aconst
      vcm(6)=bconst
      return
      end
c----------------------------------------------------------------------
      subroutine inifit_bw
c----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      parameter(nmx=600)
      common/cparm/ a(nmx),ma
      common/cint/ ix(nmx),iy(nmx),iz(nmx),mmax
      open(1,file='potdata/three.param4',status='old')
      i=1
10    read(1,*,end=100) nparm,ix(i),iy(i),iz(i),a(i)
      i=i+1
      goto 10
100   ma=i-1
      close(1)
      m=ix(ma-3)
      mmax=m+1
      return
      end
c----------------------------------------------------------------------
      subroutine inifit_cw
c----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      character*32 filnam
      parameter (maxb=15,maxpar=maxb*maxb*maxb)
      common/cdim/ kmax(6),lmax(6),nmax(6)
      common/cpar/ par_rr(10,6),par_r(10,6),par_theta(10,6)
      common/csol/ s(maxpar,6)
      data filnam/'potdata/cwfit.dat'/
      open(1,file=filnam,status='old')
      do i=1,6
        read(1,*) nmax(i),lmax(i),kmax(i)
        read(1,*) npar_rr,(par_rr(k,i),k=1,npar_rr)
        read(1,*) npar_r,(par_r(k,i),k=1,npar_r)
        read(1,*) npar_theta,(par_theta(k,i),k=1,npar_theta)
        npar=kmax(i)*lmax(i)*nmax(i)
        read(1,*) (s(k,i),k=1,npar)
      end do
      return
      end
c----------------------------------------------------------------------
cstart none
c;      subroutine jacobi(r1,r2,r3,rr,r,theta)
c;c----------------------------------------------------------------------
c;      implicit double precision (a-h,o-z)
c;      data pi180/57.29577951308232522583d0/
c;c
c;      r=r2
c;      rr=sqrt(0.5d0*(r3*r3-0.5d0*r2*r2+r1*r1))
c;      argcos=(-r1*r1+rr*rr+0.25d0*r2*r2)/(r2*rr)
c;      if(argcos.gt.1.d0) argcos=1.d0
c;      if(argcos.lt.-1.d0) argcos=-1.d0
c;      if(dabs(argcos).gt.1.d0) write(6,*) 'argcos',argcos
c;      theta=acos(argcos)*pi180
c;      return
c;      end
cend
c----------------------------------------------------------------------
      function tgcn(rr,r,theta,icase)
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (maxb=15,maxpar=maxb*maxb*maxb)
      common/cdim/ kmax(6),lmax(6),nmax(6)
      common/csol/ s(maxpar,6)
      dimension fr(20),frr(20)

      data ifirst/-1/
c
      save ifirst
c
c on first call of this subroutine, read in three-body parameter
c
      if (ifirst.eq.-1) then
         call inifit_cw
         ifirst=0
      end if
c
      nmx=nmax(icase)
      kmx=kmax(icase)
      lmx=lmax(icase)
      do n=1,nmx
        fr(n)=func_r(r,n,icase)
      end do
      do k=1,kmx
        frr(k)=func_rr(rr,k,icase)
      end do
      kk=0
      pot=0
      do l=1,lmx
        ft=func_theta(theta,l,icase)
        do n=1,nmx
          do k=1,kmx
            kk=kk+1
            pot=pot+s(kk,icase)*ft*fr(n)*frr(k)
          end do
        end do
      end do
      if(icase.eq.2) then
c... wall for pi potentials outside valid range
          if(rr.ge.3.0.and.r.lt.2.7) then
            if(pot.gt.0.085) then
             pot=1.00d0
            endif
          else
            pot=1.00d0
          endif
      end if
      tgcn=pot
      return
      end
c--------------------------------------------------------------------
      function func_theta(theta,l,icase)
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      common/cpar/ par_rr(10,6),par_r(10,6),par_theta(10,6)
c
      mm=nint(par_theta(1,icase))
      inc=nint(par_theta(2,icase))
      if(mm.eq.0) then
        ll=inc*(l-1)   !0,2,4...
      else
        ll=inc*l       !2,4,6...
      end if
c
      func_theta=0.d0
      func_theta=(-1)**mm*pm1(ll,mm,theta)
      return
      end
c--------------------------------------------------------------------
      function func_r(r,n,icase)
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      common/cpar/ par_rr(10,6),par_r(10,6),par_theta(10,6)
c
      re=par_r(1,icase)
      func_r=(r-re)**(n-1)
      return
      end
c--------------------------------------------------------------------
      function func_rr(r,k,icase)
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      common/cpar/ par_rr(10,6),par_r(10,6),par_theta(10,6)
c
      re=par_rr(1,icase)
      rdmp=par_rr(2,icase)
      wdmp=par_rr(3,icase)
      kpmx=nint(par_rr(4,icase))
      kimin=nint(par_rr(5,icase))
      rx=par_rr(6,icase)             ! check
c
      if(wdmp.eq.0) wdmp=1.d0
      w=0.5d0*(1.d0+tanh((r-rdmp)/wdmp))
      if(k.eq.1) then
        func_rr=1.d0
        if(icase.gt.4) func_rr=w           ! a and b
      else if(k.le.kpmx) then
        func_rr=(1.d0-w)*(r-re)**(k-1)
      else
        kk=k-kpmx-1+kimin
        func_rr=w*10.d0**kk/r**kk
        if(icase.gt.4) func_rr=w*10.d0**kk/(r-rx)**kk
      end if
      return
      end
cstart none
c;c--------------------------------------------------------------------
c;      function pm1(l,m,theta)
c;c--------------------------------------------------------------------
c;c
c;c  calculates value of legendre polynomial for l,m,theta
c;c
c;      implicit real*8(a-h,o-z)
c;      data pi180/.01745329251994329444d0/
c;c
c;      thedeg=theta*pi180
c;c
c;      if(m.gt.l) then
c;        pm1=0.d0
c;        return
c;      end if
c;      lmax=l
c;      x=cos(thedeg)
c;      if (m.ge.0) go to 1
c;      write (6,100)
c;100   format('  NEGATIVE M IN LEGENDRE ROUTINE:  ABORT')
c;      stop
c;c
c;1     if (m.gt.0) go to 5
c;c  here for regular legendre polynomials
c;      pm1=1.d0
c;      pm2=0.d0
c;      do 2 l=1,lmax
c;      pp=((2*l-1)*x*pm1-(l-1)*pm2)/dble(l)
c;      pm2=pm1
c;2     pm1=pp
c;      return
c;c
c;c  here for alexander-legendre polynomials
c;c
c;5     imax=2*m
c;      rat=1.d0
c;      do 6 i=2,imax,2
c;      ai=i
c;6     rat=rat*((ai-1.d0)/ai)
c;      y=sin(thedeg)
c;      pm1=sqrt(rat)*(y**m)
c;      pm2=0.d0
c;      low=m+1
c;      do 10 l=low,lmax
c;      al=(l+m)*(l-m)
c;      al=1.d0/al
c;      al2=((l+m-1)*(l-m-1))*al
c;      al=sqrt(al)
c;      al2=sqrt(al2)
c;      pp=(2*l-1)*x*pm1*al-pm2*al2
c;      pm2=pm1
c;10    pm1=pp
c;      return
c;      end
cend
c--------------------------------------------------------------------
      function dmerge(rclh1,rh2,rclh2,rr,theta)
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c
      call bw4pot(rclh1,rh2,rclh2,vsig_bw)
      call distjac(rclh1,rclh2,rr,2.5d0,theta)
      call bw4pot(rclh1,2.5d0,rclh2,vbw25)
c
c     vsig_bw=vsig_bw+0.1747736062104906d0
c     vbw25=vbw25+0.1747736062104906d0
      vsig_cw=tgcn(rr,rh2,theta,1)
      vcw25=tgcn(rr,2.5d0,theta,1)
      wr1=0.5d0*(1.d0+tanh((3.50d0-rh2)/0.5d0))
      wr2=0.5d0*(1.d0+tanh((2.500d0-rh2)/0.3d0))
      wrr=0.5d0*(1.d0+tanh((rr-4.00d0)/0.3d0))
c... correct bw potential to approximately agree with cw at rh2=2.5
      vbw_corr=vsig_bw+(vcw25-vbw25)*wrr*wr1
      dmerge=wrr*wr2*vsig_cw+(1.d0-wrr*wr2)*vbw_corr
      if(vsig_bw.gt.1.174773606210491) dmerge=vsig_bw
      return
      end
c
c----------------------------------------------------------------------
      subroutine distjac(r2,r3,rr,r,the)
c----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
c     dimension rbond(3)
      pi=4.d0*datan(1.d0)
c     r1=r
      costhe=dcos(the*pi/180.d0)
      cospithe=dcos(pi-the*pi/180.d0)
c
c     if(dabs(costhe).gt.1.d0) write(6,*) 'costhe', costhe
      arg1=(r/2.d0)**2+rr**2-2.d0*(r/2.d0)*rr*
     >            costhe
      arg2=(r/2.d0)**2+rr**2-2.d0*(r/2.d0)*rr*
     >             cospithe
      if(arg1.lt.0.d0) then
c       write(6,*) 'arg1',arg1
        arg1=0.d0
      endif
      if(arg2.lt.0.d0) then
c       write(6,*) 'arg2',arg2
        arg2=0.d0
      endif
      r2 = dsqrt(arg1)
      r3 = dsqrt(arg2)
c
c     write(7,*) 'rr,r,theta',rr,r,the, 'rbond',rbond
      return
      end
c------------------------------------------------------------------------
      subroutine bw4(x,y,z,v)
c-------------------------------------------------------------------
c
c     System:   ClH2
c     Name:     BW4
c     Author:  Wensheng Bian and Joachim Werner
c     Functional form: Aguado-Paniagua
c     Energy Zero Point: the asymptote Cl+H+H in a.u.
c
c     This subroutine calculates the potential energy for the
c     NON spin-orbit corrected MRCI+Q surface for the system
c     **Scale fact=.948d0**
c
c            Cl + H2
c
c     Input are the three distances x,y,z
c         Cl    H1    H2
c         |_____|
c            x
c               |_____|
c                  y
c         |___________|
c               z
c     This linear geometrie is at 180 Degree.
c
c
c-------------------------------------------------------------------

      implicit real*8(a-h,o-z)
      parameter(nmx=600,np=27,mma=20)
      common/cparm/ a(nmx),ma
      common/cint/ ix(nmx),iy(nmx),iz(nmx),mmax
      dimension xex(0:mma),xey(0:mma),xez(0:mma)
      dimension p(np)
      data ifirst/-1/
      data p/14.81794625, -0.05687046695,1.50963779,
     1       -19.91349307, 58.12148867,-75.88455892,
     1       36.47835698, 1.922975642, 0.7117342021,
     1       1.079946167,-0.02206944094,-7.109456997,
     1        36.79845478,-109.3716794,176.4925683,
     1        -120.4407534,2.351569228,1.082968302,
     1       14.81794625, -0.05687046695,1.50963779,
     1       -19.91349307, 58.12148867,-75.88455892,
     1       36.47835698, 1.922975642, 0.7117342021/

      save ifirst
c
c on first call of this subroutine, read in three-body parameter
c
      if (ifirst.eq.-1) then
         call inifit_bw
         ifirst=0
      end if
c
c.... Three-Body-Potential : Aguado-Paniagua
c
c.... initialize the non-linear parameters
c
      b1 = a(ma-2)
      b2 = a(ma-1)
      b3 = a(ma)
      fit = 0.0d0
      xexpon = b1*x
      yexpon = b2*y
      zexpon = b3*z
      exponx=dexp(-xexpon)
      expony=dexp(-yexpon)
      exponz=dexp(-zexpon)
      fex = x*exponx
      fey = y*expony
      fez = z*exponz
      xex(0)=1
      xey(0)=1
      xez(0)=1
      xex(1)=fex
      xey(1)=fey
      xez(1)=fez
      do m=2,mmax-1
	 xex(m)=xex(m-1)*fex
	 xey(m)=xey(m-1)*fey
	 xez(m)=xez(m-1)*fez
      enddo
c
      do 1010 i=1,ma-3
       fit=fit+xex(ix(i))*xey(iy(i))*xez(iz(i))*a(i)
1010  continue
c
c.... Two-Body-Potential : Aguado-Paniagua
c
c       c0      c1      c2      c3      c4     c5    c6     alpha  beta
c  ClH  p(1)    p(2)    p(3)    p(4)    p(5)   p(6)  p(7)   p(8)   p(9)
c  HH   p(10)   p(11)   p(12)   p(13)   p(14)  p(15) p(16)  p(17)  p(18)
c  ClH  p(19)   p(20)   p(21)   p(22)   p(23)  p(24) p(25)  p(26)  p(27)
c

      rhox=x*dexp(-p(9)*x)
      rhoy=y*dexp(-p(18)*y)
      rhoz=z*dexp(-p(27)*z)
      xval=p(1)*dexp(-p(8)*x)/x
      yval=p(10)*dexp(-p(17)*y)/y
      zval=p(19)*dexp(-p(26)*z)/z
      do 10 i=1,6
        xval=xval+p(i+1)*rhox**i
        yval=yval+p(i+10)*rhoy**i
        zval=zval+p(i+19)*rhoz**i
10    continue
c
c... Total Potential in atomic units
c
      v=fit+xval+yval+zval
cC
      if (x.lt.1.75d0.or.z.lt.1.75d0.or.y.lt.0.8d0) then
         v=1.0d0
      end if
      return
      end
c----------------------------------------------------------------------
      function damp(rr,r,theta,icase)
c----------------------------------------------------------------------
c
      implicit real*8(a-h,o-z)
c
      if(icase.eq.3) then
       w=0.5d0*(1.d0+tanh((rr-3.25d0)/0.5d0))
       w1=0.5d0*(1.d0+tanh((r-2.4d0)/0.2d0))
       w2=0.5d0*(1.d0+tanh((r-1.1d0)/0.2d0))
       damp=w*(1-w1)*w2
      elseif(icase.eq.4) then
       w=0.5d0*(1.d0+tanh((rr-4.00d0)/1.0d0))
       if(r.ge.2.5d0) w=0.5d0*(1.d0+tanh((rr-5.00d0)/0.8d0))
       w1=0.5d0*(1.d0+tanh((r-2.5d0)/0.2d0))
       w2=0.5d0*(1.d0+tanh((r-0.9d0)/0.2d0))
       damp=w*(1.d0-w1)*w2
       elseif(icase.eq.5) then
       call distjac(r2,r3,rr,r,the)
c       damp=w
        w=0.5d0*(1.d0+tanh((r2-3.0d0)/0.2d0))
        w1=0.5d0*(1.d0+tanh((r-2.6d0)/0.2d0))
        w2=0.5d0*(1.d0+tanh((r-0.9d0)/0.2d0))
        damp=w*(1-w1)*w2
      elseif(icase.eq.6) then
       call distjac(r2,r3,rr,r,the)
c      damp=w
        w=0.5d0*(1.d0+tanh((r2-3.0d0)/0.2d0))
        w1=0.5d0*(1.d0+tanh((r-2.6d0)/0.2d0))
        w2=0.5d0*(1.d0+tanh((r-0.9d0)/0.2d0))
        damp=w*(1-w1)*w2
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine bw4pot (x,y,z,v)
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c
c-------------------------------------------------------------------
c
c     System:   ClH2
c     Name:     BW4
c     Author:  Wensheng Bian and Joachim Werner
c     Functional form: Aguado-Paniagua
c     Energy Zero Point: the asymptote Cl+H2(re) in a.u.
c
c     This subroutine calculates the potential energy for the
c     NON spin-orbit corrected MRCI+Q surface for the system
c     **Scale fact=.948d0**
c
c            Cl + H2
c
c     Input are the three distances x,y,z
c         Cl    H1    H2
c         |_____|
c            x
c               |_____|
c                  y
c         |___________|
c               z
c
c-------------------------------------------------------------------
c
      parameter(nmx=600,np=27,mma=20)
      common/cparm/ a(nmx),ma
      common/cint/ ix(nmx),iy(nmx),iz(nmx),mmax
      dimension xex(0:mma),xey(0:mma),xez(0:mma)
      dimension p(np)
      save ifirst,p
      data ifirst/-1/
      data p/14.81794625, -0.05687046695,1.50963779,
     1       -19.91349307, 58.12148867,-75.88455892,
     1       36.47835698, 1.922975642, 0.7117342021,
     1       1.079946167,-0.02206944094,-7.109456997,
     1        36.79845478,-109.3716794,176.4925683,
     1        -120.4407534,2.351569228,1.082968302,
     1       14.81794625, -0.05687046695,1.50963779,
     1       -19.91349307, 58.12148867,-75.88455892,
     1       36.47835698, 1.922975642, 0.7117342021/
c
c on first call of this subroutine, read in three-body parameter
c
      if (ifirst.eq.-1) then
         write (6,61)
         call bw4ini
         ifirst=0
      end if
c
c.... Three-Body-Potential : Aguado-Paniagua
c
c.... initialize the non-linear parameters
c
      b1 = a(ma-2)
      b2 = a(ma-1)
      b3 = a(ma)
      xexpon = b1*x
      yexpon = b2*y
      zexpon = b3*z
      exponx=dexp(-xexpon)
      expony=dexp(-yexpon)
      exponz=dexp(-zexpon)
      fex = x*exponx
      fey = y*expony
      fez = z*exponz
      xex(0)=1
      xey(0)=1
      xez(0)=1
      xex(1)=fex
      xey(1)=fey
      xez(1)=fez
      do m=2,mmax-1
         xex(m)=xex(m-1)*fex
         xey(m)=xey(m-1)*fey
         xez(m)=xez(m-1)*fez
      enddo
      fit = 0.0d0
      do i=1,ma-3
         fit=fit+xex(ix(i))*xey(iy(i))*xez(iz(i))*a(i)
      enddo
c
c.... Two-Body-Potential : Aguado-Paniagua
c
c       c0      c1      c2      c3      c4     c5    c6     alpha  beta
c  ClH  p(1)    p(2)    p(3)    p(4)    p(5)   p(6)  p(7)   p(8)   p(9)
c  HH   p(10)   p(11)   p(12)   p(13)   p(14)  p(15) p(16)  p(17)  p(18)
c  ClH  p(19)   p(20)   p(21)   p(22)   p(23)  p(24) p(25)  p(26)  p(27)
c
      rhox=x*dexp(-p(9)*x)
      rhoy=y*dexp(-p(18)*y)
      rhoz=z*dexp(-p(27)*z)
      xval=p(1)*dexp(-p(8)*x)/x
      yval=p(10)*dexp(-p(17)*y)/y
      zval=p(19)*dexp(-p(26)*z)/z
      do i=1,6
         xval=xval+p(i+1)*rhox**i
         yval=yval+p(i+10)*rhoy**i
         zval=zval+p(i+19)*rhoz**i
      enddo
c
c.... Total Potential in atomic units
c
      v=fit+xval+yval+zval
c
c.... Relative to Cl+H2(re)
c
      v = v+0.1747737310d0
c
c     if (x.lt.1.75d0.or.z.lt.1.75d0.or.y.lt.0.8d0) then
c        v=1.0d0
c     end if
c
      return
  61  format(/1x,
     +'This calculation is using the BW Cl+H2 PES'/1x,
     +'Please cite: ',
     +'W.Bian and H-J.Werner, J. Chem. Phys. 112, 220 (2000)')
      end
c------------------------------------------------------------------------
      subroutine bw4ini
c------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter(nmx=600)
      common/cparm/ a(nmx),ma
      common/cint/ ix(nmx),iy(nmx),iz(nmx),mmax
      open(1,file='potdata/BW.3p',status='old')
      i=1
      rewind 1
10    read(1,*,end=100) nparm,ix(i),iy(i),iz(i),a(i)
      i=i+1
          goto 10
100   ma=i-1
      close(1)
      m=ix(ma-3)
      mmax=m+1
      return
      end
