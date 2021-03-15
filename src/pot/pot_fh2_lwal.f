*system:  F(2P)+H2, Dubernet-Hutson expansion of li-werner-alexander PES's
*references: G. Li, H.-J. Wermer, F. Lique, and M. H. Alexander, J. Chem. Phys.
*            127, 174302 (2007)
* M. H. Alexander, J. Chem. Phys. 99, 6014 (1993).
* M.-L. Dubernet and J. M. Hutson, J. Chem. Phys. 101, 1939 (1994).
*
*  Note:  this pot routine requires a data file to be in hibxx/bin/progs/potdata:
*         fhhfit_1050.dat, fhhfit_1078_new.dat
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

      potnam='LI-WERNER-ALEXANDER F(2P)H2 DUBERNET-HUTSON'
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
      potnam='LI-WERNER-ALEXANDER F(2P)H2 DUBERNET-HUTSON'
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
*  f-h2 potentials of li-werner-alexander at rh2=1.4  using body-frame
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
* latest revision date:  9-oct-2006
* ----------------------------------------------------------------------


      implicit double precision (a-h,o-z)
      logical ljunk, ihomo, csflag, lljunk
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(9)
      dimension vs(4),vzz(4),vd(4),v12(4)
      dimension d0(16),d1(4),d2(9),aa(16)
      dimension xzz(4),xs(4),xd(3),x1(2),kpvt(4),
     :          qraux(4),work(24),rsd(4)
      dimension vvll(13)
      dimension hmat(7),vev(6),rbond(3)
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
      pi=acos(-1d0)

* determine potentials at angles
      do 200 i=1,4
        theta=(4-i)*30
        rh2=1.4d0
* convert from jacobi to bond coordinates
        rbond(1)=rh2
        x=dcos(theta*pi/180d0)
        rbond(2)=dsqrt(r*r+0.25d0*rh2*rh2-x*r*rh2)
        rbond(3)=dsqrt(r*r+0.25d0*rh2*rh2+x*r*rh2)
        iflag=0
        call potsub_fh2(rbond,vev,iflag)
* ---------------------------------------------------------
* calculates lowest diabatic f+h2 surfaces and spin-orbit coupling
* from alexander-li fit to full li-werner surfaces
* spin-orbit coupling taken from earlier alexander-stark-werner fit
* input variables (distances in bohr)
*      rbond(1)=rh2
*      rbond(2)=rfh1
*      rbond(3)=rfh2

* iflag=0:  return 6 diabatic potentials
* iflag=1:  return 3 adiabatic potentials (eigenvalues of H + HSO)
* iflag=2:  return lowest electronically adiabatic potential

* iflag=0
* subroutine returns (energy in eV)
*      vev(1)=vsig
*      vev(2)=vpi
*      vev(3)=vd
*      vev(4)=v1
*      vev(5)=aconst
*      vev(6)=bconst

        vzz(i)=vev(1)
        vs(i)=vev(2)
        vd(i)=vev(3)
        v12(i)=vev(4)
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
* add on damping term
      vv0=vv0+exp(-.5*r)*.5*(tanh(-2*(r-3.5))+1)
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
*     econv=1./219474.6
      econv=8065.465/219474.6
      call dscal(9,econv,vvl,1)
      vv0=vv0*econv
      end
* ---------------------------------------------------
      subroutine potsub_fh2(rbond,vev)
* ---------------------------------------------------------
* calculates lowest diabatic f+h2 surfaces and spin-orbit coupling
* from alexander-li fit to full li-werner surfaces
* spin-orbit coupling taken from earlier alexander-stark-werner fit
* input variables (distances in bohr)
*      rbond(1)=rh2
*      rbond(2)=rfh1
*      rbond(3)=rfh2

* subroutine returns (energy in eV)
*      vev(1)=vsig
*      vev(2)=vpi
*      vev(3)=vd
*      vev(4)=v1
*      vev(5)=aconst
*      vev(6)=bconst
*    NB the spin-orbit constants are defined in Eqs. 20-22 of
*       alexander, manolopoulos, and werner

* ---------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension vev(6),rbond(3),vfin(4),
     :    hmat(6,6),work(20),evalue(6)
c fl
	  dimension vevs78(6),vevs5(6)
c fl	
*  evtocm converts from ev to cm-1
      data evtocm / 8065.465d0/
      data hartoev / 27.21164867d0/
      data zero /0d0/
* vfin contains the diabatic potentials (in cm-1)
* a,b are the two spin-orbit matrix elements (in cm-1)
* determine the coupling potentials required in the scattering calculation

c fl
      call potsub_fh2s5(rbond,vevs5)
      call potsub_fh2s78(rbond,vevs78)

c fct pot=f78*gH2*gHF+(1-gH2*gHF)*f5
c gH2=0.5*(tanh(0.7*(rHH-5.5))+1)	
c gHF=0.5*(tanh(0.7*(rHF>-5.5))+1)		
			 grhh=0.5*(tanh(0.7*(rbond(1)-5.5))+1)
			 if (rbond(2).gt.rbond(3)) then
			 grhf=0.5*(tanh(0.7*(rbond(2)-5.5))+1)
			 else
			 grhf=0.5*(tanh(0.7*(rbond(3)-5.5))+1)
			 endif			
      do i=1,6
       vev(i)=(vevs78(i)*grhh*grhf)+(1-(grhh*grhf))*vevs5(i)
	  enddo
c fl
c      vev(5)=a/evtocm
c      vev(6)=b/evtocm
      vsig=vev(1)
      vpi =vev(2)
      v2= vev(3)
      v1=vev(4)
      a= vev(5)
      b= vev(6)
c... Alexander-Li-Werner diabatic potentials
      return
      end

      subroutine potsub_fh2s5(rbond,vev)
* ---------------------------------------------------------
* calculates lowest diabatic f+h2 surfaces and spin-orbit coupling
* from alexander-li fit to full li-werner surfaces
* spin-orbit coupling taken from earlier alexander-stark-werner fit
* input variables (distances in bohr)
*      rbond(1)=rh2
*      rbond(2)=rfh1
*      rbond(3)=rfh2

* subroutine returns (energy in eV)
*      vev(1)=vsig
*      vev(2)=vpi
*      vev(3)=vd
*      vev(4)=v1
*      vev(5)=aconst
*      vev(6)=bconst
*    NB the spin-orbit constants are defined in Eqs. 20-22 of
*       alexander, manolopoulos, and werner

* ---------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension vev(6),rbond(3),vfin(4),
     :    hmat(6,6),work(20),evalue(6)
*  evtocm converts from ev to cm-1
      data evtocm / 8065.465d0/
      data hartoev / 27.21164867d0/
      data zero /0d0/
      call vfinals5(rbond,vfin,a,b)
* vfin contains the diabatic potentials (in cm-1)
* a,b are the two spin-orbit matrix elements (in cm-1)
* determine the coupling potentials required in the scattering calculation
      call dcopy(4,vfin,1,vev,1)
      call dscal(4,hartoev,vev,1)
c comment since already done in potsub_fh2
      vev(5)=a/evtocm
      vev(6)=b/evtocm
c      vsig=vev(1)
c      vpi =vev(2)
c      v2= vev(3)
c      v1=vev(4)
      a= vev(5)
      b= vev(6)
c... Alexander-Li-Werner diabatic potentials
      return
      end

      subroutine potsub_fh2s78(rbond,vev)
* ---------------------------------------------------------
* calculates lowest diabatic f+h2 surfaces and spin-orbit coupling
* from alexander-li fit to full li-werner surfaces
* spin-orbit coupling taken from earlier alexander-stark-werner fit
* input variables (distances in bohr)
*      rbond(1)=rh2
*      rbond(2)=rfh1
*      rbond(3)=rfh2

* subroutine returns (energy in eV)
*      vev(1)=vsig
*      vev(2)=vpi
*      vev(3)=vd
*      vev(4)=v1
*      vev(5)=aconst
*      vev(6)=bconst
*    NB the spin-orbit constants are defined in Eqs. 20-22 of
*       alexander, manolopoulos, and werner

* ---------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension vev(6),rbond(3),vfin(4),
     :    hmat(6,6),work(20),evalue(6)
*  evtocm converts from ev to cm-1
      data evtocm / 8065.465d0/
      data hartoev / 27.21164867d0/
      data zero /0d0/
      call vfinal(rbond,vfin,a,b)
* vfin contains the diabatic potentials (in cm-1)
* a,b are the two spin-orbit matrix elements (in cm-1)
* determine the coupling potentials required in the scattering calculation
      call dcopy(4,vfin,1,vev,1)
      call dscal(4,hartoev,vev,1)
c comment since already done in potsub_fh2
      vev(5)=a/evtocm
      vev(6)=b/evtocm
c      vsig=vev(1)
c      vpi =vev(2)
c      v2= vev(3)
c      v1=vev(4)
      a= vev(5)
      b= vev(6)
c... Alexander-Li-Werner diabatic potentials
      return
      end

       subroutine vfinals5(rbond,vfin,a,b)
c----------------------------------------------------------------------
* to determine final merged li-werner-alexander fh2 potential
*  on input rbond is a 3-vector containing the HFH cond coordinates

*   rbond(1)=rh
*   rbond(2)=rhf1
*   rbond(3)=rhf2

* on return:
*        vfin contains the  4 diabatic potentials (vsig, vpi, v2, v1)
*        a,b contain the two spin-orbit matrix elements

* energies in cm-1, distances in bohr
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension rbond(3),vdau(4),vau(4),vdnew(4),vfin(4),
     :   evec(3,3),vex(7)
      data zero,half,one,two /0d0,0.5d0,1d0,2d0/
      data hartocm,evtocm /219474.6d0,8065.465d0/
      rh=rbond(1)
      r1=rbond(2)
      r2=rbond(3)
*  theta is the HFH angle
      theta=bondangs5(r1,r2,rh)
*  determine entrance channel (FH2) jacobi coordinates
      call jacobis5(r1,rh,r2,rr,rhh,th_jac)
* determine diabatic entrance-arrangement fh2 potentials
      call pot_fhhs5(rbond,vdau)
* generate damping function for V1 and V2
      dr=rr-3.8d0
      drh2=2.3d0-rh
      xonr=0.5d0*(tanh(4d0*dr)+1d0)
      xonrh2=0.5d0*(tanh(5d0*drh2)+1d0)
      damp=xonr*xonrh2
* damp out Vpi potentials if far outside the reactant region
      dr=rr-3d0
      drh2=3.5-rh
      xonr2=0.5d0*(tanh(4d0*dr)+1d0)
      xonrh2=0.5d0*(tanh(5d0*drh2)+1d0)
      damp2=xonrh2*xonr2
*     vdau(1)=damp2*vdau(1)+(one-damp2)*20000d0/hartocm
      vdau(2)=damp2*vdau(2)+(one-damp2)*40000d0/hartocm
c filter diabatic diagonal potentials, nothing < -12000 and nothing > 40000
      do i=1,2
        if (vdau(i)*hartocm.gt.40d3.or.vdau(i)*hartocm.lt.-12d3) then
           vdau(i)=40d3/hartocm
        elseif (vdau(i)*hartocm.lt.-12d3) then
           vdau(i)=-12d3/hartocm
        endif
      enddo
* damp out V1 and V2
      do i=3,4
        vdau(i)=damp*vdau(i)
      enddo

* determine adiabatic entrance-channel fh2 potentials by diagonalization
      call vadiabs5(vdau,vau,evec)
      rr1=max(r1,r2)
      rr2=min(r1,r2)
* exit channel potential
      call vexits5(rr1,rr2,vex)
* damp and filter out vex at small unphysical values of the coordinates
      rho=3.5d0-rr1
      damp=0.5d0*(tanh(-3d0*(rho-rr2))+1d0)
      do i=1,7
        vex(i)=vex(i)*damp+(1d0-damp)*40000d0
        vex(i)=min(vex(i),40000d0)
      enddo

* interpolate over angle
      call vinterps5(vex,vexintp,theta)
* intermediate region potential
      call vints5(rr1,rr2,rh,vin)
* merge vin/vex
      call vmerges5(r1,r2,rh,theta,vin,vexintp,vmer)
* filter merged potential, with upper limit of 40000
      vmer=min(vmer,40d3)
      xlineent=0.5d0*(tanh(4d0*(rh-2.6d0))+1d0)
      facent=0.5d0*(tanh(-6d0*(rr2-xlineent-2.9d0))+1d0)
* vcomp is obtained by merging vint and vexit and then merging the result with vent
      vcomp=facent*vmer+(1d0-facent)*vau(1)*hartocm
* if still in entrance region, then
* replace vau(1) with composite potential and then redetermine diabatic potentials
      if ((1d0-facent).gt. 0.001d0) then
         vau(1)=vcomp/hartocm
         if (vau(1).gt.vau(2).and.facent.gt.0.4d0) then
             vau(2)=vau(1)
         endif
         call vdiab_adiabs5(evec,vau,vdnew)
* else, just replace lowest (Vsig) potential with composite potential
      else
         vdnew(1)=vcomp/hartocm
         vdnew(2)=max(vau(2),vau(1))
         vdnew(3)=0d0
         vdnew(4)=0d0
      endif
      call dcopy(4,vdnew,1,vfin,1)
* redetermine adiabatic potentials (no longer required here 9/22/06)
*     call vadiab(vdnew,eval,evec)
*     call dcopy(3,eval,1,vau,1)
*     vau(4)=zero
*     call dcopy(4,vau,1,vfin,1)
* determine spin-orbit matrix elements
*  convert to entrance channel jacobi coordinates (for F-H2 arrangment)
*  note, F-H2 jacobi coordinates are used in the spin-orbit fit, and this must
*  be preserved even for F-HD scattering calculations)
      call asos5(rr,rbond(1),a,b)
* damp spin-orbit matrix elements outside of entrance region
      a=a*(1d0-facent)
      b=b*(1d0-facent)
      return
      end
      subroutine vmerges5(r1,r2,rh,ang,vin,vex,vmer)
* merges interaction and exit fh2 potentials
*  vcomp = fi*vin + fex*vex
*     where vin and vex are the computed potentials in the interaction and exit regions
*  and
*     fex=0.5*[tanh(-5*(rl-rho(rl,th))+1)

*     fin=1-fex
*  where
*     rl=min(r1,r2)
*     rg=max(r1,r2)
*     th=the HFH bond angle
*  and
*     rho=1.15+0.1*theta/90+0.5*rg
      implicit double precision (a-h,o-z)
      data zero,half,one,two,three /0d0,0.5d0,1d0,2d0,3d0/
      data xslope / 0.25d0/
      rr1=max(r1,r2)
      rr2=min(r1,r2)
      rho=1.15d0+(0.1d0*ang/90d0)+xslope*rr1
      fex=half*(tanh(-5d0*(rr2-rho))+one)
      fin=one-fex
      vmer=vin*fin+vex*fex
      return
      end
      subroutine vints5(r1,r2,rh,v)
* subroutine to determine interaction region fhh potential (av5z+q, scaling factor 1.078)
*     r1, r2 are the two hf bond distances (r1.ge.r2);  rh is the h2 bond distance
*     on return, v is the potential in cm-1 with respect to f+h2(r=re)

* note that the total potential is equal to vint plus vhf(r1)+vhf(r2)+vhh(rh)
*                        i   j  k     j  k
*	  v =  sum  C   x  [x  x   + x  x  ]  where i+j+k .le. N
*              ijk   ijk 1   2  3     3  2
*     and
*          x = r  exp[-alph  (r  - r )
*           1   HH         HH  HH   1

*          x = r  exp[-alph  (r  - r )
*           2   HF         HF  HF   2
*     and
*          x = r   exp[-alph  (r   - r )
*           3   HF'         HF  HF'   2

*
* 181 points ab initio points with scale=1.05
* chosen by Matlab script interaction_fit5z.m with -7000 .le. E .le. 3000 cm-1
* with E being the total energy (vthree+vtwo) above F+H2(r=re)
* and rhf1+rhf2+rhh .le. 6.5 and min(rhf1,rhf2)<3.5

* least-squares fit with weights chosen as sqrt(1/rb) where rb is the distance
* in bond coordinates from the F-HH barrier, defined by (in bohr)
* rHF=2.916;rH'F=3.783;rHH'=1.457;

* for N = 8 (78 terms)
* weighted rms deviation is 43.2 cm-1 (0.124 kcal/mol)

* fit carried out by matlab script "interaction_fit5z.m"

      implicit double precision (a-h,o-z)
      dimension xlamint(4),ccint(78),ipow(78,3),vecpow(8,3)
* powers i, j and k, stored in column order, first column is i (78 rows), 2nd column is j, etc.
      data ipow /
     : 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,
     : 1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,
     : 2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,4,4,
     : 4,4,4,4,4,4,5,5,5,5,5,6,6,6,7,1,1,1,1,1,1,
     : 1,2,2,2,2,2,3,3,3,4,0,0,0,0,0,0,0,1,1,1,1,
     : 1,1,2,2,2,2,3,3,0,0,0,0,0,0,1,1,1,1,1,2,2,
     : 2,3,0,0,0,0,0,1,1,1,1,2,2,0,0,0,0,1,1,1,2,
     : 0,0,0,1,1,0,0,1,0,1,2,3,4,5,6,7,2,3,4,5,6,
     : 3,4,5,4,1,2,3,4,5,6,7,1,2,3,4,5,6,2,3,4,5,
     : 3,4,1,2,3,4,5,6,1,2,3,4,5,2,3,4,3,1,2,3,4,
     : 5,1,2,3,4,2,3,1,2,3,4,1,2,3,2,1,2,3,1,2,1,
     : 2,1,1/
* nonlinear parameters
      data xlamint /
     : 9.7237749741975121d-1, 1.1775082203719407d+0,
     : 9.2939400486219026d-1, 9.4077298113144014d-1/

* expansion coefficients, column order
      data ccint /
     : 1.0873180316739499d+6, -1.1433503223996919d+7,
     : 2.1213891383117963d+7, -1.0666636552076850d+7,
     : -1.5404858394757487d+7, 2.4173469062621724d+7,
     : -1.0608401613123562d+7, 3.2745517576821446d+7,
     : -1.2006157725571850d+8, 7.6619358688805729d+7,
     : -2.2549900531593929d+6, -4.6765802680971036d+6,
     : 1.0461657219656453d+8, -1.2276937681685026d+8,
     : 1.4723076297063634d+6, 3.5311094088769130d+7,
     : -1.4950632477776124d+6, 1.0444070591387285d+7,
     : -2.9414808693604767d+7, 4.0930816147927091d+7,
     : -2.7898205353006206d+7, 8.5729530781263337d+6,
     : -1.3978326488962572d+6, 4.5490912715259759d+5,
     : 2.4985589689279837d+6, 1.1140181487066777d+7,
     : -2.7988320401775062d+7, 9.4404267997046541d+6,
     : 7.5682434058275139d+6, -3.6741686847277537d+7,
     : 1.1434026549786685d+8, -6.7478267405666724d+7,
     : 1.1777663856016986d+7, -6.3172604659944616d+7,
     : 4.0180811335355595d+7, 5.8126537444099113d+6,
     : -2.9469196842880290d+7, 5.9242780291445784d+7,
     : -6.1611338843251288d+7, 3.1635843008256156d+7,
     : -5.2635443025976196d+6, -6.0802666203220496d+6,
     : 2.0648216690510020d+7, -4.5257703443535469d+7,
     : 6.5607317113019936d+7, -3.8050295968380570d+7,
     : 4.6350662476575896d+6, -2.8235477705365062d+7,
     : 1.9759824183044795d+7, 5.9937994738190342d+6,
     : -9.3477363939526621d+6, 4.2058118745594651d+7,
     : -5.5711782817555688d+7, 2.9988690160162050d+7,
     : -5.8133043224471435d+6, 1.0065660225711560d+7,
     : -2.0012643170973390d+7, 1.3044504849138251d+7,
     : -6.1561582841983403d+5, 1.1254851935326835d+7,
     : -8.0436434619849408d+6, 5.9556088150778320d+6,
     : -3.3125490624115031d+7, 3.1923765377100963d+7,
     : -9.0655994260603711d+6, -9.2400495117019825d+6,
     : 6.8917251313092802d+6, -4.5969271232828693d+6,
     : -8.2153947465614963d+5, 1.4889146034692626d+6,
     : 1.2067568133280948d+7, -6.2923176468872204d+6,
     : 5.4205151299770074d+6, -7.6195999378108629d+5,
     : -3.8827819173716530d+6, -1.5781320013872362d+6,
     : -1.4020551730292330d+6, 1.4619197653594327d+6/
* dissociation energies
      data dehh,dehf /3.849017649305869d+4,
     :  4.964745220485715d+4/

      x2=r1*exp(-xlamint(2)*(r1-xlamint(4)))
      x3=r2*exp(-xlamint(2)*(r2-xlamint(4)))
      x1=rh*exp(-xlamint(1)*(rh-xlamint(3)))
      fact1=1d0
      fact2=1d0
      fact3=1d0
      do npow=0,7
         vecpow(npow+1,1)=fact1
         vecpow(npow+1,2)=fact2
         vecpow(npow+1,3)=fact3
         fact1=fact1*x1
         fact2=fact2*x2
         fact3=fact3*x3
      enddo

      vthree=0d0
      do ind=1,78
         i1=ipow(ind,1)+1
         i2=ipow(ind,2)+1
         i3=ipow(ind,3)+1
         vthree=vthree+ccint(ind)*vecpow(i1,1)*
     :    (vecpow(i2,2)*vecpow(i3,3)+vecpow(i2,3)*vecpow(i3,2))
      enddo
* add on two-body terms
      v=vthree+vhfs5(r1)+vhfs5(r2)+vhhs5(rh)
* shift up energies by de(HH) so that zero of energy is F+H2(re)
      v=v+dehh
      return
      end

      subroutine vexits5(r1,r2,vvec)
* subroutine to determine exit region fhh potential (av5z+q, scaling factor 1.078)
*     r1, r2 are the two hf bond distances (r1.ge.r2);
*     on return, vvec is the potential in cm-1 with respect to f+h2(r=re) at
*        HFH angles 0:30:180

* note that the total potential is equal to vexit plus vhf(r2), where r2 is the shorter
* hf bond distance
      implicit double precision (a-h,o-z)
      dimension xlamexit(7,4),ccexit(36,7),ipow(36,2),ang(7),
     :          vvec(7)
      dimension vecpow(9,2)

* fit to three-body FHH potential in the exit region in aguado-panaguia variables
*                        j  k
*	  v  = sum  C   x  x     where j+k .le.  norder
*                 jk jk  2  3
*     and
*          x = r  exp[-lam(1)(r  - lam(3) )
*           2   HF             HF
*     and
*          x = r   exp[-lam(2)(r   - lam(4) )
*           3   HF'             HF'

*     and r   is the shorter HF bond
*          HF'

*  with j>= 1 (since three-body term goes to zero as rHF goes to infinity

*  fit points (interaction + exit) for each value of theta (HFH angle) = 0:30:180
*  include 791 points with E <20000 above F+H2(re)
*  fit done by matlab script "double_fit_r_exp_exit.m"
* N = 8; 36 coefficients
* theta = 0; 89 points; rms. deviation = 11.9 cm-1
* theta = 30; 110 points; rms. deviation = 17.07 cm-1
* theta = 60; 117 points; rms. deviation = 1.34 cm-1
* theta = 90; 119 points; rms. deviation = 2.51 cm-1
* theta = 120; 118 points; rms. deviation = 2.28 cm-1
* theta = 150; 119 points; rms. deviation = 3.82 cm-1
* theta = 180; 119 points; rms. deviation = 2.45 cm-1

* total deviation is 7.85 cm-1

* dev=sqrt(11.9^2*89+17.07^2*110+1.34^2*117+2.51^2*119+2.28^2*118+3.82^2*119+2.45^2*119)/ ...
*     sqrt(89+110+117+119+118+119+119);


* lambda matrix (stored in column order, 1st colum is lam(1), 2nd column is lam(2), etc. ;
*      rows correspond to theta_bond=0:30:180)
      data xlamexit /
     : 7.4832155231496778d-1, 1.0374830655586664d0,
     : 9.9988070491394310d-1, 1.1562318048817795d0,
     : 1.1404268590404998d0, 1.2017932561853544d0,
     : 1.0796709204834483d0, 8.6177290731423695d-1,
     : 7.8052081435424192d-1, 8.4659239292142008d-1,
     : 8.2123670794651726d-1, 8.8378479718195535d-1,
     : 8.4442108838528629d-1, 9.1467740097755212d-1,
     : 1.0253399980128646d0, 1.1337666631138372d0,
     : 1.1348164807601373d0, 1.0736562123081277d0,
     : 1.1074937974854844d0, 1.1694608246500628d0,
     : 1.2002329851559868d0, 1.0620788859785000d0,
     : 8.1688367906697068d-1, 8.6081728012642378d-1,
     : 8.9785713599839889d-1, 9.2453528650896000d-1,
     : 9.5753020826800173d-1, 9.5561003195740291d-1/

* coefficient matrix (column order. columns, successively, correspond to th_bond=0:30:180
* there are 36 coefficients for each theta
      data (ccexit(j,1),j=1,36)/
     : 2.3736206570603034d8, -1.9122361453458605d9,
     : 6.5688929388139954d9, -1.2480351612881508d10,
     : 1.4170463465823809d10, -9.6191410664089680d9,
     : 3.6157117469332738d9, -5.8070211367844427d8,
     : -5.4937369594956346d7, 3.4007328702802992d8,
     : -8.2353766500739908d8, 9.7660141610105801d8,
     : -5.6779410371578276d8, 1.2897975018574241d8,
     : 6.0589611761345703d5, 5.2309129888335221d7,
     : -4.6156361256902218d8, 1.3813994668187597d9,
     : -1.8841466928868303d9, 1.2110392241722927d9,
     : -2.9893785232231998d8, 1.4923044313937077d8,
     : -6.8925149581672800d8, 1.1415584837441192d9,
     : -8.2266923099382520d8, 2.2069348122726807d8,
     : 5.7689366549483165d7, -1.5386488082755768d8,
     : 1.6204539505310306d8, -6.4445326217074588d7,
     : 4.8404563175906846d6, -4.5900645748539075d7,
     : 3.9009232345378138d7, 1.9616865589563239d7,
     : -1.8454252857085846d7, 3.7696687834629171d4/
      data (ccexit(j,2),j=1,36)/
     : 1.1740715516063435d9, -1.0130583145828569d10,
     : 3.7186659668751762d10, -7.5262692011095886d10,
     : 9.0682852978274033d10, -6.5025223174944656d10,
     : 2.5683272256236912d10, -4.3083342395185719d9,
     : -1.9312730360372117d9, 1.5056394204586435d10,
     : -4.8778785639854561d10, 8.4028194183070541d10,
     : -8.1129278684537674d10, 4.1599937930580215d10,
     : -8.8449614292820034d9, 4.5305696972041380d8,
     : -2.8204977054398303d9, 7.0153572203428917d9,
     : -8.9068011709590015d9, 5.8468625162577124d9,
     : -1.5924975375313072d9, -1.1816375624064630d8,
     : 7.2444600592488527d8, -1.0701222330430840d9,
     : 3.3397335158401090d8, 1.3912242191710919d8,
     : -1.1963484172996163d8, -1.9957003748576382d8,
     : 8.2296131294308996d8, -5.1524074095304173d8,
     : 2.9351507643553299d8, -5.0395027564822358d8,
     : 2.2495594937965760d8, -9.2007840448453322d7,
     : 6.8960791506162420d7, 1.5733431663664030d7/
      data (ccexit(j,3),j=1,36)/
     : -1.2155622892700927d6, 6.5609687132860096d6,
     : -6.0242219569305470d6, -3.3144074116139162d7,
     : 1.0694995590866864d8, -1.3523324559340787d8,
     : 8.1499458097673088d7, -1.9398162456637215d7,
     : 6.6913732165402127d6, -5.0969934059635162d7,
     : 1.5933578410608542d8, -2.6053189915502185d8,
     : 2.3509158307309487d8, -1.1120567464946434d8,
     : 2.1640988834528439d7, -5.1340681810287787d6,
     : 3.4098974097240850d7, -8.9600412163714841d7,
     : 1.1295563017925240d8, -6.5889488823687233d7,
     : 1.3284250369769933d7, -9.9735266882812241d4,
     : -4.0143167440363015d6, 2.1354851595520865d7,
     : -3.5368540822670162d7, 1.9249977175118230d7,
     : 1.8354863485925891d6, -1.1490166837662872d7,
     : 2.2749635278334189d7, -1.5309278860000495d7,
     : 1.3000637181819205d6, -5.1903794921514215d6,
     : 6.5136718741408382d6, 4.3981867225633148d5,
     : -2.1050106806601551d6, 4.3734591568956326d5/
      data (ccexit(j,4),j=1,36)/
     : -1.8721946101258609d7, 1.6163113605798772d8,
     : -5.9892408728954744d8, 1.2379795102153544d9,
     : -1.5460739090203981d9, 1.1697593639688888d9,
     : -4.9751731279622704d8, 9.1869211357912511d7,
     : 3.7857358567199342d7, -2.7085083599253291d8,
     : 7.8297965949657273d8, -1.1539657707671447d9,
     : 8.9551789915059888d8, -3.3226690495917004d8,
     : 4.0829541180287495d7, -3.7887208686317503d7,
     : 2.9767266167574757d8, -9.3320619013842440d8,
     : 1.4343626369985290d9, -1.0761160028807509d9,
     : 3.1464323841777295d8, -5.7737225781269647d7,
     : 3.7114258988228726d8, -8.2675375656569076d8,
     : 7.7820767936674356d8, -2.6233682468779641d8,
     : -6.6284732230212167d7, 2.3538308724771383d8,
     : -2.7298023278557074d8, 9.7694284929761946d7,
     : -4.8898017606926104d6, 5.4577954733096715d6,
     : 7.8688879811068373d6, 1.8538348711176496d6,
     : -7.9851093582122382d6, 1.8885162108907348d6/
      data (ccexit(j,5),j=1,36)/
     : -8.9890227636319294d5, 1.0258883616936402d7,
     : -5.0068094039857790d7, 1.3480709875842977d8,
     : -2.1599983066155091d8, 2.0571602097674966d8,
     : -1.0770363109951158d8, 2.3888666545970522d7,
     : -1.2023019092132917d6, 1.6165457341492368d7,
     : -7.8947084449976966d7, 1.9315477474965444d8,
     : -2.5296635336010414d8, 1.6940746220551521d8,
     : -4.5498731843654007d7, -8.5679339349010736d6,
     : 7.1054527208546788d7, -2.3350066868234554d8,
     : 3.6990942781099671d8, -2.8151314487498808d8,
     : 8.1994982773764357d7, -1.5639050927694717d7,
     : 9.8264436235198230d7, -2.0941162362306046d8,
     : 1.8484587457761016d8, -5.5552220225247175d7,
     : -1.4646173019978052d7, 4.3605474546154603d7,
     : -3.7837319612259060d7, 3.4164465867277249d6,
     : 2.4033458564257682d6, -9.3970537562477421d6,
     : 1.3795977703736229d7, 1.0975127906133854d6,
     : -5.7245140746546565d6, 1.3421818755894545d6/
      data (ccexit(j,6),j=1,36)/
     : -6.5459218645088468d6, 6.2013287818551265d7,
     : -2.5086652504508418d8, 5.6127190167953563d8,
     : -7.4978941519125497d8, 5.9780127849502838d8,
     : -2.6326522821332237d8, 4.9376636572647668d7,
     : 2.0726660551462474d6, -1.3112404174123233d7,
     : 3.5666771362296402d7, -5.4608295293084614d7,
     : 5.3223878605252765d7, -3.2295341637831748d7,
     : 9.1895352086176369d6, -8.3658985638640048d6,
     : 4.6196859950509094d7, -9.2118060265166759d7,
     : 6.7965578797608182d7, 5.1122265159775176d4,
     : -1.4313571589412276d7, -5.6856675937159511d5,
     : -2.0398187904443428d7, 1.0005818710857201d8,
     : -1.4461960619377339d8, 6.7454202660703450d7,
     : 1.2156558225855639d7, -7.3375066492021799d7,
     : 1.2126571008141057d8, -6.3755911532781616d7,
     : 1.3514229618024856d7, -3.3979941812226504d7,
     : 2.4715131228455607d7, 1.5777528383815540d6,
     : -4.3080324400176276d6, 7.5944624153614289d5/
      data (ccexit(j,7),j=1,36)/
     : -1.6738335208459906d5, 3.0309802972347797d6,
     : -1.8756758142484784d7, 5.7632959073840499d7,
     : -9.9199246971820310d7, 9.7608763341773346d7,
     : -5.1404412248845413d7, 1.1253274910337476d7,
     : -1.0945938610250866d6, 1.0742629592215812d7,
     : -4.1734132877416916d7, 8.2848330686146677d7,
     : -8.8434640044789031d7, 4.8083958476547204d7,
     : -1.0372406546021089d7, -2.5928328754342180d6,
     : 1.6711505315170128d7, -3.8705146278779104d7,
     : 3.8902573426537946d7, -1.4631854417543182d7,
     : 2.6414799774344091d5, -8.3756605920417211d5,
     : -2.4123979369338471d6, 1.7637139340368554d7,
     : -2.6676423915811013d7, 1.2504612929741751d7,
     : 2.1593547694505276d6, -9.9907195905217342d6,
     : 1.5597080349832319d7, -8.2583250357681438d6,
     : 6.0829519854375091d5, -2.3515753244663738d6,
     : 2.4380021226242594d6, 2.5316048126764424d5,
     : -7.7196136575518828d5, 1.5501809627226790d5/
* powers j and k, stored in column order, first column is j (36 rows), 2nd column is k
      data ipow /
     : 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3,
     : 4, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 7, 7, 8, 0, 1, 2, 3, 4, 5, 6,
     : 7, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 0, 1,
     : 2, 3, 0, 1, 2, 0, 1, 0/

* dissociation energies
      data dehh /3.849017649305869d4/

      rhf1=max(r1,r2)
      rhf2=min(r1,r2)
* first determine hf two-body potential
      vtwo=vhfs5(r2)

* determine potential at each value of theta_bond
      do ith=1,7
         ang(ith)=(ith-1)*30d0
         x1=rhf1*exp(-xlamexit(ith,1)*(rhf1-xlamexit(ith,3)))
         x2=rhf2*exp(-xlamexit(ith,2)*(rhf2-xlamexit(ith,4)))
         fact1=1d0
         fact2=1d0
         do indpow=0,8
            vecpow(indpow+1,1)=fact1
            vecpow(indpow+1,2)=fact2
            fact1=fact1*x1
            fact2=fact2*x2
         enddo
         vv=0d0
         do ind=1,36
            vv=vv+ccexit(ind,ith)*
     :      vecpow(ipow(ind,1)+1,1)*vecpow(ipow(ind,2)+1,2)
         enddo
*  add on two body potential and hh dissociation energy
         vvec(ith)=vv+vtwo+dehh
      enddo
      return
      end
      double precision function vhfs5(r)
* determine hf potential (in wavenumbers) from fit, zero of energy corresponds to separated atoms
* data for hf potential
      implicit double precision (a-h,o-z)
      dimension cchf(8)
      data alphhf,rminhf,dehf / 0.68d0, 1.733684292230549d0,
     : 4.964745220485715d4/
      data cchf /
     : -3.4823868911924538d-1, 1.5267977622057252d-1,
     : -7.2237267852378628d-1, 1.6357646888480004d0,
     : -2.6815568255367697d0, 2.9634731670446355d0,
     : 2.5008606698426363d-4, -9.9999952824329252d-1/

      xhf=1d0-exp(-alphhf*(r-rminhf))
      vhfs5=polyvals5(cchf,xhf,7)*dehf
      return
      end

      double precision function vhhs5(r)
* determine hh potential (in wavenumbers) from fit, zero of energy corresponds to separated atoms
* data for hh potential
      implicit double precision (a-h,o-z)
      dimension cchh(8)
      data alphhh,rminhh,dehh / 0.593d0, 1.400508623037738d0,
     : 3.849017649305869d4/
      data cchh /
     : -8.4701361177139678d-1, 1.4510683719615993d0,
     : -2.0149691433837162d0, 2.1739804120395574d0,
     : -2.7737259910252416d0, 3.0086702155704170d0,
     : 1.9655538446759183d-3, -9.9997585582414161d-1/

      xhh=1d0-exp(-alphhh*(r-rminhh))
      vhhs5=polyvals5(cchh,xhh,7)*dehh
      return
      end
      double precision function bondangs5(r1,r2,rh)
      implicit double precision (a-h,o-z)
      pi=acos(-1d0)
* to determine angle between bonds r1 and r2
* rh is the third bond distance
      arg=(r1*r1+r2*r2-rh*rh)/(2d0*r1*r2)
      if (arg.gt.1d0) then
         if ((arg-1d0).gt.1.d-12) then
           stop 'arg .gt. 1 in acos in bondang'
         else
           arg=1d0
         endif
      elseif (arg.lt.-1d0) then
         if ((-1d0-arg).gt.1d-12) then
           stop 'arg .lt. -1 in acos in bondang'
         else
           arg=-1d0
         endif
      endif
      bondangs5=acos(arg)*180d0/pi
      return
      end
      double precision function vhh_sigs5(r)
* determine hh potential (in wavenumbers) from fit to guoliang's r=30,theta=0 values,
* for sigma state with scale factor s=1.05
* zero of energy corresponds to separated atoms
* data for guoliang fitted hh sigma potential
      implicit double precision (a-h,o-z)
      dimension cchh(8)
      data alphhh,rminhh,dehh / 0.593d0, 1.399503627960480d0,
     : 5.127448025524163d+4/
      data cchh /
     : 3.881698429370069d0, -5.674190121162882d0,
     : 6.545181723543575d-1, 2.265318660682579d0,
     : -2.393842586647567d0, 2.261637700194363d0,
     : 4.874564738147347d-3, -9.141997200609288d-6/

      xhh=1d0-exp(-alphhh*(r-rminhh))
      vhh_sigs5=polyvals5(cchh,xhh,7)*dehh
      return
      end

      double precision function vhh_pis5(r)
* determine hh potential (in wavenumbers) from fit to guoliang's r=30 value,
* for pi state
* zero of energy corresponds to separated atoms
* data for guoliang fitted hh sigma potential
      implicit double precision (a-h,o-z)
      dimension cchh(8)
      data alphhh,rminhh,dehh / 0.593d0, 1.399451652057266d0,
     : 3.890021378541989d4/
      data cchh /
     : 4.400938417113789d0, -6.076025925696183d0,
     : -3.251786438793819d-1, 3.046186755453157d0,
     : -3.014720762311264d0, 2.964527756937803d0,
     : 4.282510171711792d-3, 1.066373995780180d-4/

      xhh=1d0-exp(-alphhh*(r-rminhh))
      vhh_pis5=polyvals5(cchh,xhh,7)*dehh
      return
      end

      double precision function polyvals5(c,x,n)
* to evaluate a polynomial of degree n in the variable x
* this uses the reverse ordering of matlab
* sum(i,0:n) c(n-i)*x^i
      implicit double precision (a-h,o-z)
      dimension c(1)
      polyvals5=c(n+1)
      if (n.eq.1) return
      xx=x
      do i=2,n+1
         polyvals5=polyvals5+c(n-i+2)*xx
         xx=xx*x
      enddo
      return
      end

      subroutine plegendres5(lmaxs5,theta,p,incp)
c
c  generates regular-legendre polynomials for 0.le.l.le.lmax
c
      implicit double precision (a-h,o-z)
      dimension p(1)
      data zero, one,  two,        rad
     :     /0.d0, 1.d0, 2.d0, 57.29577951308232d0/
      x=cos(theta/rad)
      ll = 1
      p(ll)=one
      ll=ll+incp
      pm1s5=one
      pm2=zero
      do l=1,lmaxs5
         pp=((two*l-one)*x*pm1s5-(l-one)*pm2)/dble(l)
         p(ll)=pp
         ll=ll+incp
         pm2=pm1s5
         pm1s5=pp
      enddo
      return
      end
c----------------------------------------------------------------------
       subroutine pot_fhhs5(rbond,vau)
c----------------------------------------------------------------------
       implicit double precision (a-h,o-z)
*
* calculates lowest diabatic F+H2 surfaces
*
*      rbond(1)=rh2
*      rbond(2)=rfh1
*      rbond(3)=rfh2
*
* return 4 diabatic potentials

* The resulting energy values are in Hartree
*
*      vau(1)=Vsig
*      vau(2)=Vpi= 1/2 [ Vpi(A") + Vpi(A') ]
*      vau(3)=v2 = 1/2 [ Vpi(A") - Vpi(A') ]
*      vau(4)=v1 = < Vsig | H | Vpi(A') > / sqrt(2)
*
c
      parameter (npes=4)
      dimension  vau(npes),rbond(3),h(6,6)
c
      data toev/27.21139610d0/
      data tocm/219474.63067d0/
      data sq2i/0.7071067811865475d0/
      data sq2 /1.4142135623730950d0/
      data dehh /38490.18d0/
      data inifhhs5/0/
      save inifhhs5
c
      rh2=rbond(1)
      rfh1=rbond(2)
      rfh2=rbond(3)
c
      if(inifhhs5.eq.0) then
        call inifit_fhhs5
        inifhhs5=1
      end if
      call jacobis5(rfh1,rh2,rfh2,rr,r,theta)
c
      do icase=1,npes
        vau(icase)=fhhpots5(rr,r,theta,icase)
* ensure that the two-body vhh potential is the same for both states
        if (icase.eq.1) then
           vau(icase)=vau(icase)+(dehh-vhh_sigs5(r)+vhhs5(r))/tocm

        elseif (icase.eq.2) then
           vau(icase)=vau(icase)+(dehh-vhh_pis5(r)+vhhs5(r))/tocm
        endif
      end do
      vau(4)=vau(4)*sq2i
      call dscal(36,0d0,h,1)
      return
      end
c----------------------------------------------------------------------
      subroutine inifit_fhhs5
c----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      character*32 filnam
      parameter (npes=4,maxb=15,maxpar=maxb*maxb*maxb)
      common/cdims5/ kmaxs5(npes),lmaxs5(npes),nmaxs5(npes)
      common/cpars5/ par_rrs5(10,npes),par_rs5(10,npes),par_ts5(10,npes)
      common/csols5/ ss5(maxpar,npes)
      data filnam/'potdata/fhhfit_1050.dat'/
      open(91,file='potdata/fhhfit_1050.dat',
     :  status='old')
      do i=1,npes
        read(91,*) nmaxs5(i),lmaxs5(i),kmaxs5(i)
        read(91,*) npar_rr,(par_rrs5(k,i),k=1,npar_rr)
        read(91,*) npar_r,(par_rs5(k,i),k=1,npar_r)
        read(91,*) npar_theta,(par_ts5(k,i),k=1,npar_theta)
        npar=kmaxs5(i)*lmaxs5(i)*nmaxs5(i)
        read(91,*) (ss5(k,i),k=1,npar)
c        write(6,*) 'kmax(i),lmax(i),nmax(i)',kmax(i),lmax(i),nmax(i)
c        write(6,*) 'npar_R',npar_rr,(par_rr(k,i),k=1,npar_rr)
c        write(6,*) 'npar_r',npar_r,(par_r(k,i),k=1,npar_r)
c        write(6,*) 'npar_theta',npar_theta,(par_theta(k,i),
c     >              k=1,npar_theta)
c        write(6,*) 's',(s(k,i),k=1,npar)
      end do
      return
      end
c--------------------------------------------------------------------
      function fhhpots5(rr,r,theta,icase)
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (npes=4,maxb=15,maxpar=maxb*maxb*maxb)
      common/cdims5/ kmaxs5(npes),lmaxs5(npes),nmaxs5(npes)
      common/csols5/ ss5(maxpar,npes)
      dimension fr(maxb),frr(maxb)
c
      nmx=nmaxs5(icase)
      kmx=kmaxs5(icase)
      lmx=lmaxs5(icase)
      if(nmx.gt.maxb) then
        write(6,*) 'nmaxs5.gt.maxb:',nmx,maxb
        stop
      end if
      if(kmx.gt.maxb) then
        write(6,*) 'kmaxs5.gt.maxb:',kmx,maxb
        stop
      end if
      do n=1,nmx
        fr(n)=func_rs5(r,n,icase)
      end do
      do k=1,kmx
        frr(k)=func_rrs5(rr,k,icase)
      end do
      kk=0
      pot=0
      do l=1,lmx
        ft=func_thetas5(theta,l,icase)
        do n=1,nmx
          do k=1,kmx
            kk=kk+1
            pot=pot+ss5(kk,icase)*ft*fr(n)*frr(k)
          end do
        end do
      end do
      fhhpots5=pot
      return
      end
c--------------------------------------------------------------------
      function func_thetas5(theta,l,icase)
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (npes=4)
      common/cpars5/ par_rrs5(10,npes),par_rs5(10,npes),par_ts5(10,npes)
c
      mm=nint(par_ts5(1,icase))
      inc=nint(par_ts5(2,icase))
      if(mm.eq.0) then
        ll=inc*(l-1)   !0,2,4...
      else
        ll=inc*l       !2,4,6...
      end if
c
      func_thetas5=pm1s5(ll,mm,theta)
      return
      end
c--------------------------------------------------------------------
      function func_rs5(r,n,icase)
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (npes=4)
      common/cpars5/ par_rrs5(10,npes),par_rs5(10,npes),par_ts5(10,npes)
c
      re=par_rs5(1,icase)
      if(icase.le.2)then
       alph=par_rs5(2,icase)
       func_rs5=(1.0d0-exp(-alph*(r-re)))**(n-1)
      else
       func_rs5=(r-re)**(n-1)
      endif
      return
      end
c--------------------------------------------------------------------
      function func_rrs5(r,k,icase)
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (npes=4)
      common/cpars5/ par_rrs5(10,npes),par_rs5(10,npes),par_ts5(10,npes)
      data wthr/1.d-2/
c
      re=par_rrs5(1,icase)
      rdmp=par_rrs5(2,icase)
      wdmp=par_rrs5(3,icase)
      kpmx=nint(par_rrs5(4,icase))
      kimin=nint(par_rrs5(5,icase))
c
      if(wdmp.eq.0) wdmp=1.d0
      w=0.5d0*(1.d0+tanh((r-rdmp)/wdmp))
      if(k.eq.1) then
        func_rrs5=1.d0
      else if(k.le.kpmx) then
        func_rrs5=(1.d0-w)*(r-re)**(k-1)
      else
        kk=k-kpmx-1+kimin
        rr=max(0.2d0,r)
        func_rrs5=w*10.d0**kk/rr**kk
      end if
      return
      end
c--------------------------------------------------------------------
      function pm1s5(l,m,theta)
c--------------------------------------------------------------------
c
c  calculates value of legendre polynomial for l,m,theta
c
      implicit real*8(a-h,o-z)
      data pi180/.01745329251994329444d0/
c
      thedeg=theta*pi180
c
      if(m.gt.l) then
        pm1s5=0.d0
        return
      end if
      lmaxs5=l
      x=cos(thedeg)
      if (m.ge.0) go to 1
      write (6,100)
100   format('  NEGATIVE M IN LEGENDRE ROUTINE:  ABORT')
      stop
c
1     if (m.gt.0) go to 5
c  here for regular legendre polynomials
      pm1s5=1.d0
      pm2=0.d0
      do 2 l=1,lmaxs5
      pp=((2*l-1)*x*pm1s5-(l-1)*pm2)/dble(l)
      pm2=pm1s5
2     pm1s5=pp
      return
c
c  here for alexander-legendre polynomials
c
5     imax=2*m
      rat=1.d0
      do 6 i=2,imax,2
      ai=i
6     rat=rat*((ai-1.d0)/ai)
      y=sin(thedeg)
      pm1s5=sqrt(rat)*(y**m)
      pm2=0.d0
      low=m+1
      do 10 l=low,lmaxs5
      al=(l+m)*(l-m)
      al=1.d0/al
      al2=((l+m-1)*(l-m-1))*al
      al=sqrt(al)
      al2=sqrt(al2)
      pp=(2*l-1)*x*pm1s5*al-pm2*al2
      pm2=pm1s5
10    pm1s5=pp
      return
      end
c----------------------------------------------------------------------
      subroutine distjacs5(r2,r3,rr,r,the)
c----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      data pi/3.141592653589793d0/
      data pi180/.01745329251994329444d0/
c
c... returns bond distances r2,r3 for given rr,r,theta. (r1=r)
c
      costhe=dcos(the*pi180)
      cospithe=dcos(pi-the*pi180)
      arg1=(r/2.d0)**2+rr**2-2.d0*(r/2.d0)*rr*costhe
      arg2=(r/2.d0)**2+rr**2-2.d0*(r/2.d0)*rr*cospithe
      if (arg1.lt.0d0 .or. arg1.lt.0d0) then
           print *, 'r2, r3, rr, r, the ', r2,r3,rr,r,the
           stop 'argument unreasonable in distjac'
      endif
      r2=dsqrt(arg1)
      r3=dsqrt(arg2)
*     r2 = dsqrt((r/2.d0)**2+rr**2-2.d0*(r/2.d0)*rr*costhe)
*     r3 = dsqrt((r/2.d0)**2+rr**2-2.d0*(r/2.d0)*rr*cospithe)
      return
      end
c----------------------------------------------------------------------
      subroutine jacobis5(r1,r2,r3,rr,r,theta)
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      data pi180/57.29577951308232522583d0/
c
      r=r2
      arg=0.5d0*(r3*r3-0.5d0*r2*r2+r1*r1)
      if (arg.lt.0d0) then
         if (arg.gt.-1d-9) then
            arg=0d0
         else
            print *, arg
            stop 'arg.lt.0 in jacobi'
         endif
      endif
      rr=sqrt(arg)
      if(abs(r2*rr).lt.1.d-6) then              ! avoids division by zero
        theta=0
        return
      end if
      argcos=(-r1*r1+rr*rr+0.25d0*r2*r2)/(r2*rr)
      if(argcos.gt.1.d0) argcos=0.999999999999d0
      if(argcos.lt.-1.d0) argcos=-0.999999999999d0
      theta=acos(argcos)*pi180
      return
      end
c----------------------------------------------------------------------
      subroutine vadiabs5(vdau,vau,h)
*  to determine adiabatic fh2 potentials from diabatic potentials
*  on return, the energies are in the vector vau
*             the eigenvectors (column ordered) are in the matrix h
      implicit double precision (a-h,o-z)
      dimension h(3,3),e(3),work(10),vdau(4),vau(3)
      vsig=vdau(1)
      vpi=vdau(2)
      v2=vdau(3)
      v1=vdau(4)

* make sure the two-body vhh potential is the same for both states
c
      h(1,1)=vsig
      h(2,1)=-v1
      h(3,1)=v1
c
      h(1,2)=-v1
      h(2,2)=vpi
      h(3,2)=v2
c
      h(1,3)=v1
      h(2,3)=v2
      h(3,3)=vpi
c
      call dsyev('v','l',3,h,3,e,work,10,ierr)
* eigenvalues come out ordered in terms of increasing energy
      if (ierr.ne.0) stop 'error in dsyev in vadiab'
      call dcopy(3,e,1,vau,1)
      return
      end

      subroutine vdiab_adiabs5(evec,eval,vd)
* to determine diabatic fh2 potentials from adiabatic potentials and
* original matrix of eigenvectors
* on entry:
*     evec contains the matrix of eigenvectors, by columns
*     eval is the vector of adiabatic energies
* on return:  vd contains the diabatic potentials
*     vd(1) = vsig
*     vd(2) = vpi
*     vd(3) = v2
*     vd(4) = v1
      implicit double precision (a-h,o-z)
      dimension evec(3,3),eval(3),vd(4),xmat(3,3)
      do i=1,3
         do j=1,3
            xmat(i,j)=0d0
            do k=1,3
              xmat(i,j)=xmat(i,j)+evec(i,k)*eval(k)*evec(j,k)
            enddo
         enddo
      enddo
      vd(1)=xmat(1,1)
      vd(2)=xmat(2,2)
      vd(4)=xmat(1,3)
      vd(3)=xmat(2,3)
      return
      end

      subroutine vinterps5(v,vintp,ang)
* interpolate angular dependence using double legendre interpolation
      implicit double precision (a-h,o-z)
      dimension v(7), plegen(7,7),tmat(6,6),kpvt(7),vs(7),
     :       ppl(6)

* first 7 legendre polynomials at th=[0:30:180],  stored in column order, so 1st column is Pl=0 at all angles, 2nd
* column is Pl=1 at all angles, etc

      data plegen /
     : 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 8.6602540d-1, 5d-1, 0d0,
     : -5d-1, -8.6602540d-1, -1d0, 1d0, 6.25d-1, -1.25d-1, -5d-1,
     : -1.25d-1, 6.25d-1, 1d0, 1d0, 3.2475953d-1, -4.3750000d-1, 0d0,
     : 4.3750000d-1, -3.2475953d-1, -1d0, 1d0, 2.3437500d-2,
     : -2.8906250d-1, 3.75d-1, -2.8906250d-1, 2.3437500d-2, 1d0, 1d0,
     : -2.2327217d-1, 8.9843750d-2, 0d0, -8.9843750d-2, 2.2327217d-1,
     : -1d0, 1d0, -3.7402344d-1, 3.2324219d-1, -3.1250000d-1,
     : 3.2324219d-1, -3.7402344d-1, 1d0/

* fit first three angles with 3-term legendre expansion
      do i=1,3
         call dcopy(3,plegen(1,i),1,tmat(1,i),1)
      enddo
      call dgetrf(3,3,tmat,6,kpvt,ierr)
      call dcopy(3,v,1,vs,1)
      call dgetrs('N',3,1,tmat,6,kpvt,vs,6,ierr)
      call plegendres5(5,ang,ppl,1)
      vsmallang=ddot(3,vs,1,ppl,1)

* fit angles 2-7 with 6-term legendre expansion
      do i=1,6
         call dcopy(6,plegen(2,i),1,tmat(1,i),1)
      enddo
      call dgetrf(6,6,tmat,6,kpvt,ierr)
      call dcopy(6,v(2),1,vs,1)
      call dgetrs('N',6,1,tmat,6,kpvt,vs,6,ierr)
      vlargang=ddot(6,vs,1,ppl,1)
      fact=0.5d0*(tanh(0.25d0*(ang-15d0))+1d0)
      vintp=(1d0-fact)*vsmallang+fact*vlargang
      return
      end

* ---------------------------------------------------------
      subroutine asos5(r,rh2,a,b)
* --------------------------------------------------------------
* subroutine to determine F-H2 spin-orbit couplings from
* SEC with scale factor of 1.038
* fitted by m.h. alexander as follows:
* the r and R dependence is fitted as
* note that a is is pix-piy coupling (what i call A)
* and b is pix-sigma coupling (what i call B)
*
*     v(R,r)=dexp(-xlam1*R)*(c1+c2*R+c3*R^2)
*           +(r-re)*dexp(-xlam2*R)*(c4+c5*R+c6*R^2)
*           +(r-re)^2*dexp(-xlam3*R)*(c7+c8*R+c9*R^2)
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
* --------------------------------------------------------------
* current revision date:  13-apr-1998
* --------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension xlama(6),xlamb(2),cca(3),ccb(4),xxlam(6),ccc(6)
      data one,half,two /1.d0,0.5d0,2.d0/
      data zero /0d0/
      data pi /3.141592653589793d0/
*  equilibrium h2 distance for potential expansion
      data re /1.4d0/
*  asymptotic so constant
      data ainf/ 132.41d0/
      data ione /1/
      data xlama /1.5177d0,2.8499d0,1.7330d0,3.3410d0,3.0997d0,3.8876d0/
      data xlamb /1.9942d0, 2.5334d0/
      data cca /2.9099d1, 3.1968d1, 1.2887d1/
      data ccb /-1.6802d+3, 5.0964d+2, 1.1645d+4, -6.3387d+3/
*  determine spin orbit constants
      call dcopy(6,xlama,1,xxlam,1)
      call dcopy(3,cca,1,ccc,1)
      call dasotanhs5(r,rh2,re,xxlam,ccc,aconst,nlam)
      nlam=2
      call dcopy(nlam,xlamb,1,xxlam,1)
      call dcopy(3*nlam,ccb,1,ccc,1)
      call dasodexps5(r,rh2,re,xxlam,ccc,bconst,nlam)
*  damp pi-sig so constant to zero
      dr=r-2d0
      drh2=2.4d0-rh2
      dampr=half*(tanh(3.5d0*dr)+one)
      damprh2=half*(tanh(4.5d0*drh2)+one)
      a=aconst*dampr*damprh2+dampr*ainf
*  damp pi-pi so constant
      dr=r-2.2d0
      dampr=half*(tanh(two*dr)+one)
*  evaluate pi-pi constant at 2.7
      call dasodexps5(2.7d0,rh2,re,xxlam,ccc,bconst0,nlam)
      b=bconst*dampr*damprh2+(one-dampr)*damprh2*two*bconst0+ainf

*  switch definition of a and b to be consistent with my later
*  notation:  A = pix-piy coupling and B = pi-sigma coupling
      temp=a
      a=b
      b=temp

* multiply by 1.0174 to scale up to experimental atomic values
      a=a*1.0174d0
      b=b*1.0174d0

      return
      end
* ---------------------------------------------------------

      subroutine dasodexps5(r,rh2,re,xlam,c,v,nlam)
* ------------------------------------------------------
*  to determine aso(R,r), where
*     asp(R,r)=dexp(-xlam1*R)*(c1+c2*R+c3*R^2)
*           +(r-re)*dexp(-xlam2*R)*(c4+c5*R+c6*R^2)
* and, if nlam=3
*           +(r-re)^2*dexp(-xlam3*R)*(c7+c8*R+c9*R^2)
* ------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension xlam(3),c(9)
      data one,half,two/1.d0,0.5d0,2d0/
      data zero /0d0/
      rmre=rh2-re
      vex=zero
      rho=one
      do i=1,nlam
        exx=dexp(-xlam(i)*r)
        ind=(i-1)*2+1
        vexterm=exx*(c(ind)+c(ind+1)*r)
        vex=vex+rho*vexterm
        rho=rho*rmre
      enddo
      v=vex

      return
      end
* ---------------------------------------------------------

      subroutine dasotanhs5(r,rh2,re,xlam,c,v,nlam)
* ----------------------------------------------------
*  to determine aso(R,r), where
*     aso(R,r)=c(1)*(tanh(-xlam1*(R-xlam2))-1)
*           +(r-re)*c(2)*(tanh(-xlam3*(R-xlam4))-1)
*           +(r-re)^2*c(3)*(tanh(-xlam5*(R-xlam6))-1)
* on return:  dv(1) is daso/dR and dv(2) is daso/dr
* ----------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension xlam(6),c(3)
      data one,half,two/1.d0,0.5d0,2d0/
      data zero /0d0/
* asymptotic spin-orbit constant
      data ainf / 132.41d0/
      rmre=rh2-re
      vex=zero
      rho=one
      do i=1,3
        ind=2*(i-1)+1
        exx=tanh(xlam(ind)*(r-xlam(ind+1)))
        vexterm=(exx-one)*c(i)
        vex=vex+rho*vexterm
        rho=rho*rmre
      enddo
      v=vex

      return
      end

c ---------------------------

* ---------------------------------------------------------
      subroutine fh2pot(rbond,vvev)
* ---------------------------------------------------------
* calculates lowest diabatic f+h2 surfaces and spin-orbit coupling
* from alexander-li fit to full li-werner surfaces
* spin-orbit coupling taken from earlier alexander-stark-werner fit
* input variables (distances in bohr)
*      rbond(1)=rh2
*      rbond(2)=rfh1
*      rbond(3)=rfh2

* return lowest adiabatic potentials (eigenvalues of H + HSO)
* iflag=2:  return lowest electronically adiabatic potential

*    NB the spin-orbit constants are defined in Eqs. 20-22 of
*       alexander, manolopoulos, and werner

* ---------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension vev(6),rbond(3),vfin(4),vex(7),
     :   hmat(6,6),work(20),evalue(6)
*  cmtoev converts from cm-1 to ev
      data cmtoev / 8065.465d0/
      data hartoev / 27.21164867d0/
      call vfinal(rbond,vfin,a,b)
* vfin contains the diabatic potentials (in cm-1)
* a,b are the two spin-orbit matrix elements (in cm-1)
* determine the coupling potentials required in the scattering calculation
      call dcopy(4,vfin,1,vev,1)
      call dscal(4,hartoev,vev,1)
      vev(5)=a/cmtoev
      vev(6)=b/cmtoev
c
c... Alexander-Li-Werner adiabatic potentials
c
         vsig=vev(1)
         vpi =vev(2)
         v2= vev(3)
         v1=vev(4)
         a= vev(5)
         b= vev(6)
c
         call dscal(36,0d0,hmat,1)
* construct lower triangle of hamiltonian matrix [Eq.(20) of
* Alexander-Manolopoulos-Werner
         hmat(1,1)=vsig
         hmat(2,2)=vsig
         hmat(3,3)=vpi
         hmat(4,4)=vpi
         hmat(5,5)=vpi
         hmat(6,6)=vpi

         hmat(3,1)=-v1
         hmat(5,1)=v1
         hmat(4,2)=-v1
         hmat(6,2)=v1
         hmat(5,3)=v2
         hmat(6,4)=v2

* add on lower triangle of spin-orbit matrix [Eq. (25) of
* Alexander-Manolopoulos-Werner

         sq2=sqrt(2d0)
         hmat(4,1)=b*sq2
         hmat(5,2)=b*sq2
         hmat(3,3)=hmat(3,3)-a
         hmat(4,4)=hmat(4,4)+a
         hmat(5,5)=hmat(5,5)+a
         hmat(6,6)=hmat(6,6)-a

*        do i=1,6
*           print *, (hmat(i,j),j=1,6)
*        enddo
         call dsyev('v','l',6,hmat,6,evalue,work,20,ierr)
* eigenvalues come out ordered in terms of increasing energy
         if (ierr.ne.0) stop 'error in dsyev in potsub_fh2'
* scale up by deltae-so/3, so that F(2P_{3/2})+H2(re) is zero of energy
         shiftso=134.7d0/cmtoev
         vvev=evalue(1)+shiftso
c
      return
      end
       subroutine vfinal(rbond,vfin,a,b)
c----------------------------------------------------------------------
* to determine final merged li-werner-alexander fh2 potential
*  on input rbond is a 3-vector containing the HFH cond coordinates

*   rbond(1)=rh
*   rbond(2)=rhf1
*   rbond(3)=rhf2

* on return:
*        vfin contains the  4 diabatic potentials (vsig, vpi, v2, v1)
*        a,b contain the two spin-orbit matrix elements

* energies in cm-1, distances in bohr
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension rbond(3),vdau(4),vau(4),vdnew(4),eval(3),vfin(4),
     :   evec(3,3),vex(7)
      data zero,half,one,two /0d0,0.5d0,1d0,2d0/
      data hartocm,evtocm /219474.6d0,8065.465d0/
      pi=acos(-1d0)
      torad=pi/180d0
      rh=rbond(1)
      r1=rbond(2)
      r2=rbond(3)
*  theta is the HFH angle
      theta=bondang(r1,r2,rh)
*  determine entrance channel (FH2) jacobi coordinates
      call jacobi_fh2(r1,rh,r2,rr,rhh,th_jac)
* determine diabatic entrance-arrangement fh2 potentials
      call pot_fhh(rbond,vdau)
* generate damping function for V1 and V2
      dr=rr-3.8d0
      drh2=2.3d0-rh
      xonr=0.5d0*(tanh(4d0*dr)+1d0)
      xonrh2=0.5d0*(tanh(5d0*drh2)+1d0)
      damp=xonr*xonrh2
* damp out Vpi potentials if far outside the reactant region
      dr=rr-3d0
      drh2=3.5-rh
      xonr2=0.5d0*(tanh(4d0*dr)+1d0)
      xonrh2=0.5d0*(tanh(5d0*drh2)+1d0)
      damp2=xonrh2*xonr2
*     vdau(1)=damp2*vdau(1)+(one-damp2)*20000d0/hartocm
      vdau(2)=damp2*vdau(2)+(one-damp2)*40000d0/hartocm
c filter diabatic diagonal potentials, nothing < -12000 and nothing > 40000
      do i=1,2
        if (vdau(i)*hartocm.gt.40d3.or.vdau(i)*hartocm.lt.-12d3) then
           vdau(i)=40d3/hartocm
        elseif (vdau(i)*hartocm.lt.-12d3) then
           vdau(i)=-12d3/hartocm
        endif
      enddo
* damp out V1 and V2
      do i=3,4
        vdau(i)=damp*vdau(i)
      enddo

* determine adiabatic entrance-channel fh2 potentials by diagonalization
      call vadiab(vdau,vau,evec)
      rr1=max(r1,r2)
      rr2=min(r1,r2)
* exit channel potential
      call vexit(rr1,rr2,vex)
* damp and filter out vex at small unphysical values of the coordinates
      rho=3.5d0-rr1
      damp=0.5d0*(tanh(-3d0*(rho-rr2))+1d0)
      do i=1,7
        vex(i)=vex(i)*damp+(1d0-damp)*40000d0
        vex(i)=min(vex(i),40000d0)
      enddo

* interpolate over angle
      call vinterp(vex,vexintp,theta)
* intermediate region potential
      call vint(rr1,rr2,rh,vin)
* merge vin/vex
      call vmerge(r1,r2,rh,theta,vin,vexintp,vmer)
* filter merged potential, with upper limit of 40000
      vmer=min(vmer,40d3)
      xlineent=0.5d0*(tanh(4d0*(rh-2.6d0))+1d0)
      facent=0.5d0*(tanh(-6d0*(rr2-xlineent-2.9d0))+1d0)
* vcomp is obtained by merging vint and vexit and then merging the result with vent
      vcomp=facent*vmer+(1d0-facent)*vau(1)*hartocm
* if still in entrance region, then
* replace vau(1) with composite potential and then redetermine diabatic potentials
      if ((1d0-facent).gt. 0.001d0) then
         vau(1)=vcomp/hartocm
         if (vau(1).gt.vau(2).and.facent.gt.0.4d0) then
             vau(2)=vau(1)
         endif
         call vdiab_adiab(evec,vau,vdnew)
* else, just replace lowest (Vsig) potential with composite potential
      else
         vdnew(1)=vcomp/hartocm
         vdnew(2)=max(vau(2),vau(1))
         vdnew(3)=0d0
         vdnew(4)=0d0
      endif
      call dcopy(4,vdnew,1,vfin,1)
* redetermine adiabatic potentials (no longer required here 9/22/06)
*     call vadiab(vdnew,eval,evec)
*     call dcopy(3,eval,1,vau,1)
*     vau(4)=zero
*     call dcopy(4,vau,1,vfin,1)
* determine spin-orbit matrix elements
*  convert to entrance channel jacobi coordinates (for F-H2 arrangment)
*  note, F-H2 jacobi coordinates are used in the spin-orbit fit, and this must
*  be preserved even for F-HD scattering calculations)
      call aso(rr,rbond(1),a,b)
* damp spin-orbit matrix elements outside of entrance region
      a=a*(1d0-facent)
      b=b*(1d0-facent)
      return
      end
      subroutine vmerge(r1,r2,rh,ang,vin,vex,vmer)
* merges interaction and exit fh2 potentials
*  vcomp = fi*vin + fex*vex
*     where vin and vex are the computed potentials in the interaction and exit regions
*  and
*     fex=0.5*[tanh(-5*(rl-rho(rl,th))+1)

*     fin=1-fex
*  where
*     rl=min(r1,r2)
*     rg=max(r1,r2)
*     th=the HFH bond angle
*  and
*     rho=1.15+0.1*theta/90+0.5*rg
      implicit double precision (a-h,o-z)
      data zero,half,one,two,three /0d0,0.5d0,1d0,2d0,3d0/
      data xslope / 0.25d0/
      rr1=max(r1,r2)
      rr2=min(r1,r2)
      rho=1.15d0+(0.1d0*ang/90d0)+xslope*rr1
      fex=half*(tanh(-5d0*(rr2-rho))+one)
      fin=one-fex
      vmer=vin*fin+vex*fex
      return
      end
      subroutine vint(r1,r2,rh,v)
* subroutine to determine interaction region fhh potential (av5z+q, scaling factor 1.078)
*     r1, r2 are the two hf bond distances (r1.ge.r2);  rh is the h2 bond distance
*     on return, v is the potential in cm-1 with respect to f+h2(r=re)

* note that the total potential is equal to vint plus vhf(r1)+vhf(r2)+vhh(rh)
*                        i   j  k     j  k
*	  v =  sum  C   x  [x  x   + x  x  ]  where i+j+k .le. N
*              ijk   ijk 1   2  3     3  2
*     and
*          x = r  exp[-alph  (r  - r )
*           1   HH         HH  HH   1

*          x = r  exp[-alph  (r  - r )
*           2   HF         HF  HF   2
*     and
*          x = r   exp[-alph  (r   - r )
*           3   HF'         HF  HF'   2

*
* 183 points ab initio points with scale=1.078
* chosen by Matlab script interaction_fit5z.m with -3000 .le. E .le. 5000 cm-1
* with E being the total energy (vthree+vtwo) above F+H2(r=re)
* and rhf1+rhf2+rhh .le. 6 and min(rhf1,rhf2)<3.5

* least-squares fit with weights chosen as sqrt(1/rb) where rb is the distance
* in bond coordinates from the F-H-H barrier, defined by (in bohr)
* rHF=3.783;rH'F=2.916;rHH'=1.457;

* for N = 8 (78 terms)
* weighted rms deviation is 43.75 cm-1 (0.125 kcal/mol)

* fit carried out by matlab script "interaction_fit5z.m"

      implicit double precision (a-h,o-z)
      dimension xlamint(4),ccint(78),ipow(78,3),vecpow(8,3)
* powers i, j and k, stored in column order, first column is i (78 rows), 2nd column is j, etc.
      data ipow /
     : 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,
     : 1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,
     : 2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,4,4,
     : 4,4,4,4,4,4,5,5,5,5,5,6,6,6,7,1,1,1,1,1,1,
     : 1,2,2,2,2,2,3,3,3,4,0,0,0,0,0,0,0,1,1,1,1,
     : 1,1,2,2,2,2,3,3,0,0,0,0,0,0,1,1,1,1,1,2,2,
     : 2,3,0,0,0,0,0,1,1,1,1,2,2,0,0,0,0,1,1,1,2,
     : 0,0,0,1,1,0,0,1,0,1,2,3,4,5,6,7,2,3,4,5,6,
     : 3,4,5,4,1,2,3,4,5,6,7,1,2,3,4,5,6,2,3,4,5,
     : 3,4,1,2,3,4,5,6,1,2,3,4,5,2,3,4,3,1,2,3,4,
     : 5,1,2,3,4,2,3,1,2,3,4,1,2,3,2,1,2,3,1,2,1,
     : 2,1,1/
* nonlinear parameters
      data xlamint /
     : 9.9020116339715725d-1, 1.1790131793220460d0,
     : 7.4819595961496321d-1, 1.0406161571905130d0/
* expansion coefficients, column order
      data ccint /
     : 1.0662357123708795d+6, -1.0646845324696966d+7,
     : 1.9918755864083238d+7, -1.4631486686268322d+7,
     : -1.4575521989598433d+6, 8.5057097107862905d+6,
     : -3.9275016003009258d+6, 2.6780582234398130d+7,
     : -9.0871231162840590d+7, 6.0543830269520380d+7,
     : -1.2057659107850952d+7, 6.1481717136681755d+5,
     : 6.8079370112550795d+7, -7.3054337531326622d+7,
     : 4.3771235566539411d+6, 1.6711752021393748d+7,
     : -9.2470595700706029d+5, 6.0886808779531838d+6,
     : -1.7403366448466156d+7, 2.4822205961477175d+7,
     : -1.6999466035859611d+7, 4.8774748956580469d+6,
     : -5.4620628474739031d+5, -1.0913698065339681d+6,
     : 1.3701249295238266d+7, -1.6440535425101094d+7,
     : 4.9262676654568976d+6, -2.9240719748585513d+6,
     : 4.3799290503262132d+6, -3.7664344474084273d+7,
     : 1.0642329417998898d+8, -6.2447563306067824d+7,
     : 1.2277809558501219d+7, -4.9595492965457246d+7,
     : 2.9075168272162024d+7, 4.7047857401350727d+6,
     : -1.8901837300595589d+7, 3.5572055843661629d+7,
     : -3.9632318078740492d+7, 2.2654691295965727d+7,
     : -4.0401471311823670d+6, -2.7789848426553304d+6,
     : -3.3535602648699754d+6, -7.4089883899556554d+5,
     : 2.4059998091738462d+7, -2.0755119932801228d+7,
     : 1.1135846677257065d+7, -3.3274078146734931d+7,
     : 1.9035455591548670d+7, 5.3780044044083171d+6,
     : -1.0508738910021720d+7, 3.6740382512055844d+7,
     : -3.9393496011739224d+7, 1.8027940451243363d+7,
     : -4.1669606472067372d+6, 9.4993611498711444d+6,
     : -8.0246158756344754d+5, -9.8002328932105824d+6,
     : 9.3401079111355003d+6, 5.8283516734410515d+6,
     : -5.6501088955467688d+6, 9.4671550502373520d+6,
     : -4.0963353259376802d+7, 3.3498293690306894d+7,
     : -7.2942499389726091d+6, -1.5076683554745713d+7,
     : 8.5621966369693971d+5, -3.1075383368317434d+6,
     : 2.3986585261276071d+6, 1.9028844486086518d+6,
     : 1.9436026509143617d+7, -9.2849437378589157d+6,
     : 1.3368032996528506d+7, -3.9200482704314339d+5,
     : -8.6577078437469974d+6, -2.8569976912512463d+6,
     : -4.7074624117781036d+6, 4.1215049573936104d+6/

* dissociation energies
      data dehh,dehf /3.862607006368740d4,4.983413055325700d+4/

      x2=r1*exp(-xlamint(2)*(r1-xlamint(4)))
      x3=r2*exp(-xlamint(2)*(r2-xlamint(4)))
      x1=rh*exp(-xlamint(1)*(rh-xlamint(3)))
      fact1=1d0
      fact2=1d0
      fact3=1d0
      do npow=0,7
         vecpow(npow+1,1)=fact1
         vecpow(npow+1,2)=fact2
         vecpow(npow+1,3)=fact3
         fact1=fact1*x1
         fact2=fact2*x2
         fact3=fact3*x3
      enddo

      vthree=0d0
      do ind=1,78
         i1=ipow(ind,1)+1
         i2=ipow(ind,2)+1
         i3=ipow(ind,3)+1
         vthree=vthree+ccint(ind)*vecpow(i1,1)*
     :    (vecpow(i2,2)*vecpow(i3,3)+vecpow(i2,3)*vecpow(i3,2))
      enddo
* add on two-body terms
      v=vthree+vhf(r1)+vhf(r2)+vhh(rh)
* shift up energies by de(HH) so that zero of energy is F+H2(re)
      v=v+dehh
      return
      end
      subroutine vexit(r1,r2,vvec)
* subroutine to determine exit region fhh potential (av5z+q, scaling factor 1.078)
*     r1, r2 are the two hf bond distances (r1.ge.r2);
*     on return, vvec is the potential in cm-1 with respect to f+h2(r=re) at
*        HFH angles 0:30:180

* note that the total potential is equal to vexit plus vhf(r2), where r2 is the shorter
* hf bond distance
      implicit double precision (a-h,o-z)
      dimension xlamexit(7,4),ccexit(36,7),ipow(36,2),ang(7),
     :          vvec(7)
      dimension coef(7),vecpow(9,2)

* fit to three-body FHH potential in the exit region in aguado-panaguia variables
*                        j  k
*	  v  = sum  C   x  x     where j+k .le.  norder
*                 jk jk  2  3
*     and
*          x = r  exp[-lam(1)(r  - lam(3) )
*           2   HF             HF
*     and
*          x = r   exp[-lam(2)(r   - lam(4) )
*           3   HF'             HF'

*     and r   is the shorter HF bond
*          HF'

*  with j>= 1 (since three-body term goes to zero as rHF goes to infinity

*  fit points (interaction + exit) for each value of theta (HFH angle) = 0:30:180
*  include 791 points with E <20000 above F+H2(re)
*  fit done by matlab script "double_fit_r_exp_exit_3body.m"
* N = 8; 36 coefficients
* theta = 0; 89 points; rms. deviation = 11.8135 cm-1
* theta = 30; 110 points; rms. deviation = 17.0458 cm-1
* theta = 60; 117 points; rms. deviation = 1.5382 cm-1
* theta = 90; 119 points; rms. deviation = 2.741 cm-1
* theta = 120; 118 points; rms. deviation = 2.8811 cm-1
* theta = 150; 119 points; rms. deviation = 4.9039 cm-1
* theta = 180; 119 points; rms. deviation = 3.0921 cm-1

* total deviation is 8.0 cm-1
* dev=sqrt(11.8135^2*89+17.05^2*110+1.5382^2*117+2.741^2*119+2.881^2*118+4.904^2*119+3.092^2*119)/ ...
*     sqrt(89+110+117+119+118+119+119);


* lambda matrix (stored in column order, 1st colum is lam(1), 2nd column is lam(2), etc. ;
*      rows correspond to theta_bond=0:30:180)
      data xlamexit /
     : 7.5766398161223947d-1, 1.0324090724729889d0,
     : 9.8946013389664911d-1, 1.1657241724245146d0,
     : 1.1485879426195940d0, 1.2315573189066709d0, 1.0863686692147301d0,
     : 8.6199160770430350d-1, 7.7962984538214286d-1,
     : 8.4283357201071307d-1, 8.2673912090874102d-1,
     : 8.9963907112857910d-1, 8.4459495932637219d-1,
     : 9.4774758419882910d-1, 1.0323140413221634d0,
     : 7.8000061120413600d-1, 7.9166183442916471d-1,
     : 8.2597129262838753d-1, 7.8049604665520467d-1,
     : 7.6460419189996698d-1, 7.9670397460322828d-1,
     : 1.0528582530640789d0, 9.0947076043328567d-1,
     : 9.5837065611316397d-1, 8.8802504677595395d-1,
     : 9.1464264435010234d-1, 9.9104541737119889d-1,
     : 1.0149604237467238d0/
* coefficient matrix (column order. columns, successively, correspond to th_bond=0:30:180
* there are 36 coefficients for each theta
      data (ccexit(j,1),j=1,36)/
     : 2.3682117378015587d8, -1.9247154166064818d9,
     : 6.6703102758341398d9, -1.2785394679496504d10,
     : 1.4645281475712681d10, -1.0028957130312275d10,
     : 3.8026351165672684d9, -6.1598171366428745d8,
     : -5.3371370041781992d7, 3.2949636914091361d8,
     : -7.9303164531955469d8, 9.2891722770865774d8,
     : -5.2563666557518649d8, 1.0919017009708355d8,
     : 4.4269968697733069d6, 5.6019331515422530d7,
     : -4.9173706850292093d8, 1.4652325458139849d9,
     : -1.9904050547839541d9, 1.2742207223270361d9,
     : -3.1321497079500961d8, 1.6029733854224113d8,
     : -7.2829775720497084d8, 1.1879513433049865d9,
     : -8.4280158322318161d8, 2.2230257398823923d8,
     : 5.2165739157156922d7, -1.3249214781909204d8,
     : 1.3679952754963905d8, -5.4619010326125517d7,
     : 3.6722935315578855d6, -4.4851126351086289d7,
     : 3.8227778565744475d7, 2.0237588866459548d7,
     : -1.8105381471581563d7, -3.8525265296890039d5/
      data (ccexit(j,2),j=1,36)/
     : 1.7560711535637739d9, -1.4073617248176088d10,
     : 4.7978604741864349d10, -9.0174614148685089d10,
     : 1.0088221579201172d11, -6.7154904381196114d10,
     : 2.4617744101855373d10, -3.8315098645974650d9,
     : -4.1335570472241478d9, 2.9955115035800678d10,
     : -9.0203614904174271d10, 1.4442892377858844d11,
     : -1.2961558293670734d11, 6.1780940882909470d10,
     : -1.2211913820440901d10, 1.3680929392253692d9,
     : -7.9370271320080576d9, 1.8381138264522224d10,
     : -2.1694752167203739d10, 1.3217968791497293d10,
     : -3.3381267677964654d9, -4.9098326733843940d8,
     : 2.8361714997850833d9, -3.9621282815320220d9,
     : 1.2219870179386461d9, 4.0847342403660274d8,
     : -7.3912887383132052d8, -1.0130763687869859d9,
     : 4.1848224391883287d9, -2.4855677986976819d9,
     : 2.4920503671938710d9, -4.0535530672677636d9,
     : 1.7334717293719869d9, -1.0724988761262555d9,
     : 7.3393304787747419d8, 2.6869730363999230d8/
      data (ccexit(j,3),j=1,36)/
     : -2.5300830346993688d6, 1.4267515969128178d7,
     : -2.3709090142227259d7, -1.1968418646028960d7,
     : 8.8963364949365512d7, -1.1898737633953933d8,
     : 6.9595063803589836d7, -1.5634712577199891d7,
     : 1.6421981151745604d7, -1.1541351230186166d8,
     : 3.3264091915345806d8, -5.0162447599108857d8,
     : 4.1741997886018002d8, -1.8194626512070718d8,
     : 3.2572862398831893d7, -1.5723112070767459d7,
     : 9.7734912348798037d7, -2.3923103054104349d8,
     : 2.8179036570581084d8, -1.5510677101004231d8,
     : 3.0213404324282579d7, -1.1342490891646121d6,
     : -1.2400594745526424d7, 6.7414895572383568d7,
     : -1.0501881118151481d8, 5.3576963120645992d7,
     : 1.0863672826116296d7, -5.7803876345754802d7,
     : 1.0379859621997039d8, -6.5546863048886217d7,
     : 8.5752323827516679d6, -3.4563933010050319d7,
     : 4.3139797236970201d7, 4.9689061081835246d6,
     : -2.2161124211213950d7, 6.8860131293157600d6/
      data (ccexit(j,4),j=1,36)/
     : -1.8568462266050663d7, 1.6265414631576848d8,
     : -6.1374975736900032d8, 1.2976733471292002d9,
     : -1.6663718323685830d9, 1.3033223602803590d9,
     : -5.7576534859836495d8, 1.1081119315218455d8,
     : 5.1193214008385181d7, -3.5808434884224182d8,
     : 9.9424917573685753d8, -1.3591144226360338d9,
     : 9.0389557071822488d8, -2.2046341650320768d8,
     : -1.1471003684912797d7, -9.1091030183150738d7,
     : 7.2754519280930221d8, -2.3138364454779954d9,
     : 3.5992336704371314d9, -2.7272636993835683d9,
     : 8.0365798713572466d8, -1.8996278334003416d8,
     : 1.2254747049010017d9, -2.7345115419331264d9,
     : 2.5741030176656508d9, -8.6400910933325136d8,
     : -2.8067335041035879d8, 9.8412124755907679d8,
     : -1.1200506538988867d9, 3.8125684112038422d8,
     : -1.2877808114642669d7, -1.2405379280026434d7,
     : 8.7213530758605748d7, 1.8903285351958092d7,
     : -7.6524445409606308d7, 2.2754644430516243d7/
      data (ccexit(j,5),j=1,36)/
     : -1.4270674839446813d5, 3.8533543703966709d6,
     : -2.8513189844582308d7, 9.9558378214779213d7,
     : -1.9134780798080420d8, 2.0870809506162462d8,
     : -1.2141471608937910d8, 2.9301921495752852d7,
     : -4.8224511116927089d6, 5.3414857591200076d7,
     : -2.3601318485510668d8, 5.4472093985094321d8,
     : -6.9047521794679892d8, 4.5493236975471640d8,
     : -1.2148059876112658d8, -2.3475766308063045d7,
     : 2.0194011930385107d8, -6.8794540919547951d8,
     : 1.1247838954986367d9, -8.7989423872865129d8,
     : 2.6197692210517752d8, -6.7346626734612554d7,
     : 4.3383756273987782d8, -9.4240280383012497d8,
     : 8.4374863458773661d8, -2.5241454849881068d8,
     : -9.0976125578014448d7, 2.6949327629343492d8,
     : -2.2237682706835699d8, -4.0495174912630296d6,
     : 2.5783159973446492d7, -1.0524537502770758d8,
     : 1.6321727708382431d8, 1.9190742352617230d7,
     : -9.8494428712488204d7, 3.1957687715119570d7/
      data (ccexit(j,6),j=1,36)/
     : -1.0154173059005873d7, 9.5284037485717386d7,
     : -3.8110684023156524d8, 8.4151548330818856d8,
     : -1.1076361979937983d9, 8.6883588303692567d8,
     : -3.7593730644896460d8, 6.9191914543022946d7,
     : 1.6585412580080036d5, 5.1760777603087770d6,
     : -2.7128275989738338d7, 4.9906201452588744d7,
     : -3.1949416406968839d7, -4.0422659295027452d6,
     : 8.2799682426194260d6, -3.8211312438167907d7,
     : 1.9523559449105552d8, -3.4887320543677205d8,
     : 1.8457048813802388d8, 9.6276660367540479d7,
     : -9.2202957258400485d7, 3.6094481955339159d6,
     : -2.2103068191259220d8, 9.4216699445583773d8,
     : -1.2891691704027736d9, 5.8356817163372540d8,
     : 2.0243258251585332d8, -1.1667936174368122d9,
     : 1.8681687959655535d9, -9.7118357749624467d8,
     : 3.6616967595754004d8, -9.0872847985853803d8,
     : 6.8125861914215362d8, 8.2849902992639467d7,
     : -2.3837150882235098d8, 7.3266531289937183d7/
      data (ccexit(j,7),j=1,36)/
     : 5.1865280760289992d4, 8.8720504444076226d5,
     : -9.7897258049652707d6, 3.6599938876717933d7,
     : -6.9192159118104622d7, 7.1523104507956073d7,
     : -3.8623893525227487d7, 8.5405390076980963d6,
     : -1.8966079033445946d6, 1.9279753763424657d7,
     : -7.5934479922401100d7, 1.5106846934943721d8,
     : -1.6047975733023474d8, 8.6606449048466727d7,
     : -1.8569310931775853d7, -8.8537004409408607d6,
     : 5.5797852791580260d7, -1.2542110259020203d8,
     : 1.2362796671768327d8, -4.8015725939839505d7,
     : 2.8373224239738001d6, -4.1496393159894445d6,
     : -1.5881449279262507d7, 8.8527307851608112d7,
     : -1.1970880168605542d8, 5.1694980417504244d7,
     : 1.8956466173436411d7, -6.8681361663053036d7,
     : 9.4497542616441905d7, -4.7208840265372090d7,
     : -1.0870798380881394d6, -1.3734542828124706d7,
     : 2.1291871075381842d7, 7.3927923794406150d6,
     : -1.5727634315929649d7, 4.1258412535141450d6/

* powers j and k, stored in column order, first column is j (36 rows), 2nd column is k
      data ipow /
     : 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3,
     : 4, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 7, 7, 8, 0, 1, 2, 3, 4, 5, 6,
     : 7, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 0, 1,
     : 2, 3, 0, 1, 2, 0, 1, 0/

* dissociation energies
      data dehh /3.862607006368740d4/

      pi=acos(-1d0)
      rhf1=max(r1,r2)
      rhf2=min(r1,r2)
* first determine hf two-body potential
      vtwo=vhf(r2)

* determine potential at each value of theta_bond
      do ith=1,7
         ang(ith)=(ith-1)*30d0
         x1=rhf1*exp(-xlamexit(ith,1)*(rhf1-xlamexit(ith,3)))
         x2=rhf2*exp(-xlamexit(ith,2)*(rhf2-xlamexit(ith,4)))
         fact1=1d0
         fact2=1d0
         do indpow=0,8
            vecpow(indpow+1,1)=fact1
            vecpow(indpow+1,2)=fact2
            fact1=fact1*x1
            fact2=fact2*x2
         enddo
         vv=0d0
         do ind=1,36
            vv=vv+ccexit(ind,ith)*
     :      vecpow(ipow(ind,1)+1,1)*vecpow(ipow(ind,2)+1,2)
         enddo
*  add on two body potential and hh dissociation energy
         vvec(ith)=vv+vtwo+dehh
      enddo
      return
      end

      double precision function vhf(r)
* determine hf potential (in wavenumbers) from fit, zero of energy corresponds to separated atoms
* data for hf potential
      implicit double precision (a-h,o-z)
      dimension cchf(8)
      data alphhf,rminhf,dehf / 0.68d0, 1.73355299392062d0,
     : 4.98341305532d4/
      data cchf /
     : -0.28938703174552d0,  0.08401953945800d0, -0.73232290770953d0,
     : 1.65168690982660d0, -2.66713872267044d0,  2.95291039409507d0,
     : 0.00023131772196d0, -0.99999950214982d0/

      xhf=1d0-exp(-alphhf*(r-rminhf))
      vhf=polyval(cchf,xhf,7)*dehf
      return
      end

      double precision function vhh(r)
* determine hh potential (in wavenumbers) from fit, zero of energy corresponds to separated atoms
* data for hh potential
      implicit double precision (a-h,o-z)
      dimension cchh(8)
      data alphhh,rminhh,dehh / 0.593d0, 1.399844598435660d0,
     : 3.862607006368740d4/
      data cchh /
     : -8.4944999231924878d-1, 1.4522303609014953d0,
     : -2.0138934690777428d0, 2.1742658991179571d0,
     : -2.7698715219419792d0, 3.0047240341151693d0,
     : 1.9713711176867162d-3, -9.9997673038218271d-1/


      xhh=1d0-exp(-alphhh*(r-rminhh))
      vhh=polyval(cchh,xhh,7)*dehh
      return
      end
      double precision function bondang(r1,r2,rh)
      implicit double precision (a-h,o-z)
      pi=acos(-1d0)
* to determine angle between bonds r1 and r2
* rh is the third bond distance
      arg=(r1*r1+r2*r2-rh*rh)/(2d0*r1*r2)
      if (arg.gt.1d0) then
         if ((arg-1d0).gt.1.d-12) then
           stop 'arg .gt. 1 in acos in bondang'
         else
           arg=1d0
         endif
      elseif (arg.lt.-1d0) then
         if ((-1d0-arg).gt.1d-12) then
           stop 'arg .lt. -1 in acos in bondang'
         else
           arg=-1d0
         endif
      endif
      bondang=acos(arg)*180d0/pi
      return
      end
      double precision function vhh_sig(r)
* determine hh potential (in wavenumbers) from fit to guoliang's r=30,theta=0 values,
* zero of energy corresponds to separated atoms
* data for guoliang fitted hh sigma potential
      implicit double precision (a-h,o-z)
      dimension cchh(8)
      data alphhh,rminhh,dehh / 0.593d0, 1.39884908035893d0,
     : 5.156676206045081d4/
      data cchh /
     : 3.8683968195946541d0, -5.6598766306159716d0,
     : 6.6554016653510495d-1, 2.2526401447621520d0,
     : -2.3858589464424980d0, 2.2542910155027815d0,
     : 4.8841581750617025d-3, -1.5615430201310036d-5/

      xhh=1d0-exp(-alphhh*(r-rminhh))
      vhh_sig=polyval(cchh,xhh,7)*dehh
      return
      end

      double precision function vhh_pi(r)
* determine hh potential (in wavenumbers) from fit to guoliang's r=30 value,
* zero of energy corresponds to separated atoms
* data for guoliang fitted hh sigma potential
      implicit double precision (a-h,o-z)
      dimension cchh(8)
      data alphhh,rminhh,dehh / 0.593d0, 1.398762271644071d0,
     : 3.903611735411907d4/
      data cchh /
     : 4.3983295074177855d0, -6.0679665316530800d0,
     : -3.3001383948384205d-1, 3.0443251427484244d0,
     : -3.0098936168732795d0, 2.9610622385212526d0,
     : 4.1694493553907103d-3, 1.0193872940456477d-4/

      xhh=1d0-exp(-alphhh*(r-rminhh))
      vhh_pi=polyval(cchh,xhh,7)*dehh
      return
      end

      double precision function polyval(c,x,n)
* to evaluate a polynomial of degree n in the variable x
* this uses the reverse ordering of matlab
* sum(i,0:n) c(n-i)*x^i
      implicit double precision (a-h,o-z)
      dimension c(1)
      polyval=c(n+1)
      if (n.eq.1) return
      xx=x
      do i=2,n+1
         polyval=polyval+c(n-i+2)*xx
         xx=xx*x
      enddo
      return
      end

      subroutine plegendre(lmax,theta,p,incp)
c
c  generates regular-legendre polynomials for 0.le.l.le.lmax
c
      implicit double precision (a-h,o-z)
      dimension p(1)
      data zero, one,  two,        rad
     :     /0.d0, 1.d0, 2.d0, 57.29577951308232d0/
      x=cos(theta/rad)
      ll = 1
      p(ll)=one
      ll=ll+incp
      pm1=one
      pm2=zero
      do l=1,lmax
         pp=((two*l-one)*x*pm1-(l-one)*pm2)/dble(l)
         p(ll)=pp
         ll=ll+incp
         pm2=pm1
         pm1=pp
      enddo
      return
      end
c----------------------------------------------------------------------
       subroutine pot_fhh(rbond,vau)
c----------------------------------------------------------------------
       implicit double precision (a-h,o-z)
*
* calculates lowest diabatic F+H2 surfaces
*
*      rbond(1)=rh2
*      rbond(2)=rfh1
*      rbond(3)=rfh2
*
* return 4 diabatic potentials

* The resulting energy values are in Hartree
*
*      vau(1)=Vsig
*      vau(2)=Vpi= 1/2 [ Vpi(A") + Vpi(A') ]
*      vau(3)=v2 = 1/2 [ Vpi(A") - Vpi(A') ]
*      vau(4)=v1 = < Vsig | H | Vpi(A') > / sqrt(2)
*
c
      parameter (npes=4)
      dimension  vau(npes),rbond(3),h(6,6),e(6)
c
      data toev/27.21139610d0/
      data tocm/219474.63067d0/
      data sq2i/0.7071067811865475d0/
      data sq2 /1.4142135623730950d0/
      data dehh /38626.07d0/
      data inifhh/0/
      save inifhh
c
      rh2=rbond(1)
      rfh1=rbond(2)
      rfh2=rbond(3)
c
      if(inifhh.eq.0) then
        call inifit_fhh
        inifhh=1
      end if
      call jacobi_fh2(rfh1,rh2,rfh2,rr,r,theta)
c
      do icase=1,npes
        vau(icase)=fhhpot(rr,r,theta,icase)
* ensure that the two-body vhh potential is the same for both states
        if (icase.eq.1) then
           xsig=vau(icase)*tocm-vhh_sig(r)
           vau(icase)=vau(icase)+(dehh-vhh_sig(r)+vhh(r))/tocm
        elseif (icase.eq.2) then
           xpi=vau(icase)*tocm-vhh_pi(r)
           vau(icase)=vau(icase)+(dehh-vhh_pi(r)+vhh(r))/tocm
        endif
      end do
*     print *,'xsig,xpi', xsig, xpi
      vau(4)=vau(4)*sq2i
      call dscal(36,0d0,h,1)
      return
      end
c----------------------------------------------------------------------
      subroutine inifit_fhh
c----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      character*32 filnam
      parameter (npes=4,maxb=15,maxpar=maxb*maxb*maxb)
      common/cdim/ kmax(npes),lmax(npes),nmax(npes)
      common/cpar/ par_rr(10,npes),par_r(10,npes),par_theta(10,npes)
      common/csol/ s(maxpar,npes)
      data filnam/'potdata/fhhfit_1078_new.dat'/
      open(91,file='potdata/fhhfit_1078_new.dat',
     :  status='old')
      do i=1,npes
        read(91,*) nmax(i),lmax(i),kmax(i)
        read(91,*) npar_rr,(par_rr(k,i),k=1,npar_rr)
        read(91,*) npar_r,(par_r(k,i),k=1,npar_r)
        read(91,*) npar_theta,(par_theta(k,i),k=1,npar_theta)
        npar=kmax(i)*lmax(i)*nmax(i)
        read(91,*) (s(k,i),k=1,npar)
c        write(6,*) 'kmax(i),lmax(i),nmax(i)',kmax(i),lmax(i),nmax(i)
c        write(6,*) 'npar_R',npar_rr,(par_rr(k,i),k=1,npar_rr)
c        write(6,*) 'npar_r',npar_r,(par_r(k,i),k=1,npar_r)
c        write(6,*) 'npar_theta',npar_theta,(par_theta(k,i),
c     >              k=1,npar_theta)
c        write(6,*) 's',(s(k,i),k=1,npar)
      end do
      return
      end
c--------------------------------------------------------------------
      function fhhpot(rr,r,theta,icase)
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (npes=4,maxb=15,maxpar=maxb*maxb*maxb)
      common/cdim/ kmax(npes),lmax(npes),nmax(npes)
      common/csol/ s(maxpar,npes)
      dimension fr(maxb),frr(maxb)
c
      nmx=nmax(icase)
      kmx=kmax(icase)
      lmx=lmax(icase)
      if(nmx.gt.maxb) then
        write(6,*) 'nmax.gt.maxb:',nmx,maxb
        stop
      end if
      if(kmx.gt.maxb) then
        write(6,*) 'kmax.gt.maxb:',kmx,maxb
        stop
      end if
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
      fhhpot=pot
      return
      end
c--------------------------------------------------------------------
      function func_theta(theta,l,icase)
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (npes=4)
      common/cpar/ par_rr(10,npes),par_r(10,npes),par_theta(10,npes)
c
      mm=nint(par_theta(1,icase))
      inc=nint(par_theta(2,icase))
      if(mm.eq.0) then
        ll=inc*(l-1)   !0,2,4...
      else
        ll=inc*l       !2,4,6...
      end if
c
      func_theta=pm1(ll,mm,theta)
      return
      end
c--------------------------------------------------------------------
      function func_r(r,n,icase)
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (npes=4)
      common/cpar/ par_rr(10,npes),par_r(10,npes),par_theta(10,npes)
c
      re=par_r(1,icase)
      if(icase.le.2)then
       alph=par_r(2,icase)
       func_r=(1.0d0-exp(-alph*(r-re)))**(n-1)
      else
       func_r=(r-re)**(n-1)
      endif
      return
      end
c--------------------------------------------------------------------
      function func_rr(r,k,icase)
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (npes=4)
      common/cpar/ par_rr(10,npes),par_r(10,npes),par_theta(10,npes)
      data wthr/1.d-2/
c
      re=par_rr(1,icase)
      rdmp=par_rr(2,icase)
      wdmp=par_rr(3,icase)
      kpmx=nint(par_rr(4,icase))
      kimin=nint(par_rr(5,icase))
c
      if(wdmp.eq.0) wdmp=1.d0
      w=0.5d0*(1.d0+tanh((r-rdmp)/wdmp))
      if(k.eq.1) then
        func_rr=1.d0
      else if(k.le.kpmx) then
        func_rr=(1.d0-w)*(r-re)**(k-1)
      else
        kk=k-kpmx-1+kimin
        rr=max(0.2d0,r)
        func_rr=w*10.d0**kk/rr**kk
      end if
      return
      end
cstart none
c;c pm1 is included in hibrid3.f
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
c----------------------------------------------------------------------
      subroutine distjac(r2,r3,rr,r,the)
c----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      data pi/3.141592653589793d0/
      data pi180/.01745329251994329444d0/
c
c... returns bond distances r2,r3 for given rr,r,theta. (r1=r)
c
      costhe=dcos(the*pi180)
      cospithe=dcos(pi-the*pi180)
      arg1=(r/2.d0)**2+rr**2-2.d0*(r/2.d0)*rr*costhe
      arg2=(r/2.d0)**2+rr**2-2.d0*(r/2.d0)*rr*cospithe
      if (arg1.lt.0d0 .or. arg1.lt.0d0) then
           print *, 'r2, r3, rr, r, the ', r2,r3,rr,r,the
           stop 'argument unreasonable in distjac'
      endif
      r2=dsqrt(arg1)
      r3=dsqrt(arg2)
*     r2 = dsqrt((r/2.d0)**2+rr**2-2.d0*(r/2.d0)*rr*costhe)
*     r3 = dsqrt((r/2.d0)**2+rr**2-2.d0*(r/2.d0)*rr*cospithe)
      return
      end
c----------------------------------------------------------------------
      subroutine jacobi_fh2(r1,r2,r3,rr,r,theta)
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      data pi180/57.29577951308232522583d0/
c
      r=r2
      arg=0.5d0*(r3*r3-0.5d0*r2*r2+r1*r1)
      if (arg.lt.0d0) then
         if (arg.gt.-1d-9) then
            arg=0d0
         else
            print *, arg
            stop 'arg.lt.0 in jacobi'
         endif
      endif
      rr=sqrt(arg)
      if(abs(r2*rr).lt.1.d-6) then              ! avoids division by zero
        theta=0
        return
      end if
      argcos=(-r1*r1+rr*rr+0.25d0*r2*r2)/(r2*rr)
      if(argcos.gt.1.d0) argcos=0.999999999999d0
      if(argcos.lt.-1.d0) argcos=-0.999999999999d0
      theta=acos(argcos)*pi180
      return
      end
c----------------------------------------------------------------------
      subroutine vadiab(vdau,vau,h)
*  to determine adiabatic fh2 potentials from diabatic potentials
*  on return, the energies are in the vector vau
*             the eigenvectors (column ordered) are in the matrix h
      implicit double precision (a-h,o-z)
      dimension h(3,3),e(3),work(10),vdau(4),vau(3)
      vsig=vdau(1)
      vpi=vdau(2)
      v2=vdau(3)
      v1=vdau(4)

* make sure the two-body vhh potential is the same for both states
c
      h(1,1)=vsig
      h(2,1)=-v1
      h(3,1)=v1
c
      h(1,2)=-v1
      h(2,2)=vpi
      h(3,2)=v2
c
      h(1,3)=v1
      h(2,3)=v2
      h(3,3)=vpi
c
      call dsyev('v','l',3,h,3,e,work,10,ierr)
* eigenvalues come out ordered in terms of increasing energy
      if (ierr.ne.0) stop 'error in dsyev in vadiab'
      call dcopy(3,e,1,vau,1)
      return
      end

      subroutine vdiab_adiab(evec,eval,vd)
* to determine diabatic fh2 potentials from adiabatic potentials and
* original matrix of eigenvectors
* on entry:
*     evec contains the matrix of eigenvectors, by columns
*     eval is the vector of adiabatic energies
* on return:  vd contains the diabatic potentials
*     vd(1) = vsig
*     vd(2) = vpi
*     vd(3) = v2
*     vd(4) = v1
      implicit double precision (a-h,o-z)
      dimension evec(3,3),eval(3),vd(4),xmat(3,3)
      do i=1,3
         do j=1,3
            xmat(i,j)=0d0
            do k=1,3
              xmat(i,j)=xmat(i,j)+evec(i,k)*eval(k)*evec(j,k)
            enddo
         enddo
      enddo
      vd(1)=xmat(1,1)
      vd(2)=xmat(2,2)
      vd(4)=xmat(1,3)
      vd(3)=xmat(2,3)
      return
      end

      subroutine vinterp(v,vintp,ang)
* interpolate angular dependence using double legendre interpolation
      implicit double precision (a-h,o-z)
      dimension v(7), plegen(7,7),tmat(6,6),kpvt(7),vs(7),
     :       ppl(6)

* first 7 legendre polynomials at th=[0:30:180],  stored in column order, so 1st column is Pl=0 at all angles, 2nd
* column is Pl=1 at all angles, etc

      data plegen /
     : 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 8.6602540d-1, 5d-1, 0d0,
     : -5d-1, -8.6602540d-1, -1d0, 1d0, 6.25d-1, -1.25d-1, -5d-1,
     : -1.25d-1, 6.25d-1, 1d0, 1d0, 3.2475953d-1, -4.3750000d-1, 0d0,
     : 4.3750000d-1, -3.2475953d-1, -1d0, 1d0, 2.3437500d-2,
     : -2.8906250d-1, 3.75d-1, -2.8906250d-1, 2.3437500d-2, 1d0, 1d0,
     : -2.2327217d-1, 8.9843750d-2, 0d0, -8.9843750d-2, 2.2327217d-1,
     : -1d0, 1d0, -3.7402344d-1, 3.2324219d-1, -3.1250000d-1,
     : 3.2324219d-1, -3.7402344d-1, 1d0/

* fit first three angles with 3-term legendre expansion
      do i=1,3
         call dcopy(3,plegen(1,i),1,tmat(1,i),1)
      enddo
      call dgetrf(3,3,tmat,6,kpvt,ierr)
      call dcopy(3,v,1,vs,1)
      call dgetrs('N',3,1,tmat,6,kpvt,vs,6,ierr)
      call plegendre(5,ang,ppl,1)
      vsmallang=ddot(3,vs,1,ppl,1)

* fit angles 2-7 with 6-term legendre expansion
      do i=1,6
         call dcopy(6,plegen(2,i),1,tmat(1,i),1)
      enddo
      call dgetrf(6,6,tmat,6,kpvt,ierr)
      call dcopy(6,v(2),1,vs,1)
      call dgetrs('N',6,1,tmat,6,kpvt,vs,6,ierr)
      vlargang=ddot(6,vs,1,ppl,1)
      fact=0.5d0*(tanh(0.25d0*(ang-15d0))+1d0)
      vintp=(1d0-fact)*vsmallang+fact*vlargang
      return
      end

* ---------------------------------------------------------
      subroutine aso(r,rh2,a,b)
* --------------------------------------------------------------
* subroutine to determine F-H2 spin-orbit couplings from
* SEC with scale factor of 1.038
* fitted by m.h. alexander as follows:
* the r and R dependence is fitted as
* note that a is is pix-piy coupling (what i call A)
* and b is pix-sigma coupling (what i call B)
*
*     v(R,r)=dexp(-xlam1*R)*(c1+c2*R+c3*R^2)
*           +(r-re)*dexp(-xlam2*R)*(c4+c5*R+c6*R^2)
*           +(r-re)^2*dexp(-xlam3*R)*(c7+c8*R+c9*R^2)
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
* --------------------------------------------------------------
* current revision date:  13-apr-1998
* --------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension xlama(6),xlamb(2),cca(3),ccb(4),xxlam(6),ccc(6)
      data one,half,two /1.d0,0.5d0,2.d0/
      data zero /0d0/
      data pi /3.141592653589793d0/
*  equilibrium h2 distance for potential expansion
      data re /1.4d0/
*  asymptotic so constant
      data ainf/ 132.41d0/
      data ione /1/
      data xlama /1.5177d0,2.8499d0,1.7330d0,3.3410d0,3.0997d0,3.8876d0/
      data xlamb /1.9942d0, 2.5334d0/
      data cca /2.9099d1, 3.1968d1, 1.2887d1/
      data ccb /-1.6802d+3, 5.0964d+2, 1.1645d+4, -6.3387d+3/
*  determine spin orbit constants
      call dcopy(6,xlama,1,xxlam,1)
      call dcopy(3,cca,1,ccc,1)
      call dasotanh(r,rh2,re,xxlam,ccc,aconst,nlam)
      nlam=2
      call dcopy(nlam,xlamb,1,xxlam,1)
      call dcopy(3*nlam,ccb,1,ccc,1)
      call dasodexp(r,rh2,re,xxlam,ccc,bconst,nlam)
*  damp pi-sig so constant to zero
      dr=r-2d0
      drh2=2.4d0-rh2
      dampr=half*(tanh(3.5d0*dr)+one)
      damprh2=half*(tanh(4.5d0*drh2)+one)
      a=aconst*dampr*damprh2+dampr*ainf
*  damp pi-pi so constant
      dr=r-2.2d0
      dampr=half*(tanh(two*dr)+one)
*  evaluate pi-pi constant at 2.7
      call dasodexp(2.7d0,rh2,re,xxlam,ccc,bconst0,nlam)
      b=bconst*dampr*damprh2+(one-dampr)*damprh2*two*bconst0+ainf

*  switch definition of a and b to be consistent with my later
*  notation:  A = pix-piy coupling and B = pi-sigma coupling
      temp=a
      a=b
      b=temp

* multiply by 1.0174 to scale up to experimental atomic values
      a=a*1.0174d0
      b=b*1.0174d0

      return
      end
* ---------------------------------------------------------

      subroutine dasodexp(r,rh2,re,xlam,c,v,nlam)
* ------------------------------------------------------
*  to determine aso(R,r), where
*     asp(R,r)=dexp(-xlam1*R)*(c1+c2*R+c3*R^2)
*           +(r-re)*dexp(-xlam2*R)*(c4+c5*R+c6*R^2)
* and, if nlam=3
*           +(r-re)^2*dexp(-xlam3*R)*(c7+c8*R+c9*R^2)
* ------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension xlam(3),c(9)
      data one,half,two/1.d0,0.5d0,2d0/
      data zero /0d0/
      rmre=rh2-re
      rmre2=rmre*rmre
      r2=r*r
      vex=zero
      rho=one
      rho2=one
      do i=1,nlam
        exx=dexp(-xlam(i)*r)
        ind=(i-1)*2+1
        vexterm=exx*(c(ind)+c(ind+1)*r)
        vex=vex+rho*vexterm
        rho=rho*rmre
      enddo
      v=vex

      return
      end
* ---------------------------------------------------------

      subroutine dasotanh(r,rh2,re,xlam,c,v,nlam)
* ----------------------------------------------------
*  to determine aso(R,r), where
*     aso(R,r)=c(1)*(tanh(-xlam1*(R-xlam2))-1)
*           +(r-re)*c(2)*(tanh(-xlam3*(R-xlam4))-1)
*           +(r-re)^2*c(3)*(tanh(-xlam5*(R-xlam6))-1)
* on return:  dv(1) is daso/dR and dv(2) is daso/dr
* ----------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension xlam(6),c(3)
      data one,half,two/1.d0,0.5d0,2d0/
      data zero /0d0/
* asymptotic spin-orbit constant
      data ainf / 132.41d0/
      rmre=rh2-re
      rmre2=rmre*rmre
      r2=r*r
      vex=zero
      rho=one
      rho2=one
      do i=1,3
        ind=2*(i-1)+1
        exx=tanh(xlam(ind)*(r-xlam(ind+1)))
        vexterm=(exx-one)*c(i)
        vex=vex+rho*vexterm
        rho=rho*rmre
      enddo
      v=vex

      return
      end

