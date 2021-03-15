*system:  F(2P)+H2, Dubernet-Hutson expansion of li-werner-alexander PES's
*references: G. Li, H.-J. Wermer, F. Lique, and M. H. Alexander, J. Chem. Phys.
*            127, 174302 (2007)
* M. H. Alexander, J. Chem. Phys. 99, 6014 (1993).
* M.-L. Dubernet and J. M. Hutson, J. Chem. Phys. 101, 1939 (1994).
*
*  Note:  this pot routine requires a data file to be in hibxx/bin/progs/potdata:
*         fhhfit_1078_new.dat
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
      subroutine potsub_fh2(rbond,vev,iflag)
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
*    NB the spin-orbit constants are defined in Eqs. 20-22 of
*       alexander, manolopoulos, and werner
 
* ---------------------------------------------------------
      
      implicit double precision (a-h,o-z)
      dimension vev(6),rbond(3),vfin(4),
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
      if (iflag.eq.0) then
c... Alexander-Li-Werner diabatic potentials
        return
c
      else if(iflag.eq.1) then
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
         call dcopy(3,evalue,2,vev,1)
c
      else if(iflag.eq.2) then
c
c... lowest electronically-adiabatic li-werner potential surface
c
         call dscal(36,0d0,hmat,1)
         hmat(1,1)=vsig
         hmat(2,1)=-v1
         hmat(3,1)=v1
         hmat(2,2)=vpi
         hmat(3,2)=v2
         hmat(3,3)=vpi
         call dsyev('v','l',3,hmat,6,evalue,work,20,ierr)
         vev(1)=evalue(1)
         vev(2)=zero
         vev(3)=zero
      end if
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
     :   evec(3,3)
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
      call jacobi_n(r1,rh,r2,rr,rhh,th_jac)
* determine diabatic entrance-arrangement fh2 potentials
      call pot_fhh(rbond,vdau)
* generate damping function for V1 and V2
      dr=rr-3.8d0
      drh2=2d0-rh
      xonr=0.5d0*(tanh(4d0*dr)+1d0)
      xonrh2=0.5d0*(tanh(5d0*drh2)+1d0)
      damp=xonr*xonrh2
c filter diabatic diagonal potentials, nothing < -12000 and nothing > 30000
      do i=1,2
        if (vdau(i)*hartocm.gt.30d3.or.vdau(i)*hartocm.lt.-12d3) then
           vdau(i)=30d3/hartocm
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
* intermediate region potential
      call vint(rr1,rr2,rh,vin)
* merge vin/vex
      call vmerge(r1,r2,rh,theta,vin,vex,vmer)
* filter merged potential, with upper limit of 30000
      vmer=min(vmer,30d3)
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
*  convert to entrance channel jacobi_n coordinates (for F-H2 arrangment)
*  note, F-H2 jacobi_n coordinates are used in the spin-orbit fit, and this must
*  be preserved even for F-HD scattering calculations)
      call aso(rr,rbond(1),a,b)
* damp spin-orbit matrix elements outside of entrance region
      a=a*(1d0-facent)
      b=b*(1d0-facent)
      return
      end
      subroutine vmerge(r1,r2,rh,theta,vin,vex,vmer)
* merges interaction and exit fh2 potentials
*  vcomp = fi*vin + fex*vex
*     where vin and vex are the computed potentials in the interaction and exit regions
*  and
*     fin=0.5*[tanh(-alph(theta)*(rg-rho(rl,th))+1)
*     fex=1-fin
*  where
*     rl=min(r1,r2)
*     rg=max(r1,r2)
*     theta=the HFH bond angle in degrees
*  and
*     the value of fin is determined by the subroutine dampint(r1,r2,theta)
      
      implicit double precision (a-h,o-z)
      data zero,half,one,two,three /0d0,0.5d0,1d0,2d0,3d0/
      rr1=max(r1,r2)
      rr2=min(r1,r2)
      fin=dampint(r1,r2,theta)
      vmer=vin*fin+vex*(one-fin)
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
* --------------------------------------------------------------
       subroutine vexit(r1,r2,vex)
* --------------------------------------------------------------
* subroutine to determine FH-H adiabatic potential
* based on av5z+q Li-Werner ab initio calculations scaled by s=1.078
*     on return, vvec is the potential in cm-1 with respect to f+h2(r=re) at
*        HFH angles 0:30:180

* 3body term fitted by m.h. alexander as follows:
* at each r (1.325, 1.4,1.5, 1.6, 1.75, 1.9, 2.1, 2.4, 2.8) and
* and theta = 0:30:180 (colinear HFH corresponds to theta = 180)
* the R dependence is fitted as

* V=  c1*exp(-lam1*R) + (c2+c3*R)exp(-lam2*R)
*     -0.5[tanh{alpha(R-lam3)} +1]*c4/R^^

* here r=min(r1,r2), R=max(r1,r2), and theta=<HFH

* for a given value of R=ri, the v(R)'s are determined at each r and each
* theta from the above equation, using parameters given below.  then, the
* r dependence for each theta is fitted to a 5th order polynomial in
* the morse variable x=1-exp(-alphhf(r-rehf)) by an overdetermined
* least-squares procedure (subroutines dgels) to obtain
* v(rr1,rr2) for each theta = [0:30:180].  at small and large r the
* quadratic (and higher) in (r-re) terms are damped to zero by the functions
*  0.25*(tanh(4*(r-.8))+1)*(tanh(-4*(r-3.3))+1)

* --------------------------------------------------------------
* current revision date:  30-sep-2006
* --------------------------------------------------------------
       implicit double precision (a-h,o-z)
       dimension xlam1(63),xlam2(63),xlam3(63),c1(63),c2(63),
     :      c3(63),c4(63),v(9,7),vander(9,6),scmat(9,6),
     :      cvec(9),aux(390),rhf(9),vex(7)
       data one,half /1.d0,0.5d0/
       data zero /0d0/
* conversion from ev to cm-1
       data tocm1 / 8065.465d0/
       data ione /1/
       data rhf/1.325d0,1.4d0,1.5d0,1.6d0,1.75d0,1.9d0,
     :         2.1d0,2.4d0,2.8d0/
* dissociation energies
      data dehh /3.862607006368740d4/


       data xlam1/
     : 1.9076114d0, 1.5992063d0, 7.5484310d-1, 7.7976308d-1,
     : 1.0525954d0, 1.9445483d0, 7.2321892d-1, 6.5100458d-1,
     : 1.5711723d0, 1.0565989d0, 8.1944590d-1, 1.1716260d0,
     : 8.9896040d-1, 1.4463376d0, 2.2680681d0, 1.5294725d0,
     : 8.4916550d-1, 9.2768905d-1, 7.0084926d-1, 1.2579662d0,
     : 7.1392500d-1, 3.1680495d0, 8.8684289d-1, 8.8253480d-1,
     : 9.6902276d-1, 9.0407506d-1, 1.0278586d0, -2.2567400d-1,
     : 7.6971235d-1, 9.2157485d-1, 4.5322548d0, 1.0803611d0,
     : 8.7217794d-1, 7.6573322d-1, 7.1791023d-1, 7.9363873d-1,
     : 7.4295415d-1, 8.4962304d-1, 8.3068216d-1, 7.9807115d-1,
     : 7.2519798d-1, 7.2181675d-1, 8.1326700d-1, 5.4259385d-1,
     : 8.6581536d-1, 7.9209599d-1, 7.6762514d-1, 7.1703299d-1,
     : 7.5537084d-1, 8.4993151d-1, 8.0691961d-1, 9.0972124d-1,
     : 7.3675237d-1, 7.1974003d-1, 7.0246727d-1, 7.6248453d-1,
     : 6.7853133d-1, 2.0030349d0, 8.4038253d-1, 7.2535139d-1,
     : 7.5186168d-1, 7.5443561d-1, 1.7473501d0/ 
       data xlam2/
     : 9.9131618d-1,8.0409922d-1, 2.8776557d0, 1.9759148d0, 5.8301258d0,
     : 6.6465487d0, 2.2607301d0, 1.8371577d0, 8.2597182d-1, 1.5516773d0,
     : 1.5081550d0, 2.1940102d0, 2.7917569d0, 1.2359453d0, 1.0943465d0,
     : 8.6321516d-1,1.9027327d0, 1.4871332d0, 1.5997524d0, 8.4089986d-1,
     : 2.2182882d0, 1.1586475d0, 1.5345120d0, 1.9247990d0, 1.4686451d0,
     : 1.5405238d0, 1.5079987d0, 1.3233602d0, 2.2503794d0, 1.5044235d0,
     : 1.2327957d0, 1.3736128d0, 1.8080844d0, 2.1337447d0, 2.2258640d0,
     : 2.6893732d0, 1.6419048d0, 1.9442202d0, 1.8665240d0, 1.9734540d0,
     : 2.2339063d0, 2.2417905d0, 3.3924032d0, 1.9101295d0, 1.8615164d0,
     : 1.9607936d0, 2.0611527d0, 2.2736029d0, 2.1030574d0, 2.6634054d0,
     : 2.3224658d0, 1.7618471d0, 2.0803205d0, 2.2002804d0, 2.3811139d0,
     : 2.2210240d0, 1.9793879d0, 9.4850924d-1, 1.8588432d0, 2.1721535d0,
     : 2.2285380d0, 2.2944437d0, 1.7555372d0/
       data xlam3/
     : 5.6606701d0, 5.3929205d0, 2.3215788d0, 3.4936634d0, 4.8458105d0,
     : -7.6454129d0, 4.7783133d0, 4.7031483d0, 5.4340727d0, 3.7243822d0,
     : 4.7353133d0, 4.9112268d0, 3.3849925d0, 4.2269704d0, 4.5064867d0,
     : 5.5293510d0, 4.5272666d0, 4.6980769d0, 4.7617113d0, 5.7486673d0,
     : 4.6264094d0, 5.5884203d0, 4.5105068d0, 5.5537530d0, 4.6634154d0,
     : 4.6685999d0, 4.3321470d0, 4.0706194d0, 3.7156110d0, 4.8504312d0,
     : 4.8887678d0, 4.5931756d0, 4.6196951d0, 4.4114733d0, 4.8290652d0,
     : 3.7163507d0, 5.0011603d0, 4.9667995d0, 4.7415636d0, 4.6910484d0,
     : 4.8501245d0, 4.8009805d0, 3.8979705d0, 2.9601990d0, 4.8857891d0,
     : 4.8338274d0, 4.7311322d0, 4.9077515d0, 4.5839913d0, 4.5128972d0,
     : 4.2313771d0, 4.8733447d0, 4.8611427d0, 4.8773854d0, 5.0289720d0,
     : 3.2522123d0, 5.0509885d0, 4.7641307d0, 5.0952156d0, 3.1736188d0,
     : 3.2085048d0, 3.7885077d0, 4.1905534d0/
       data c1/
     : 3.6457355080663362d+6, 1.4573353289175567d+6,
     : -2.4420975404098743d+4, -1.7009562254799930d+4,
     : 7.1567766669081670d+4, 3.8941157129583880d+6,
     : -1.8342306289282740d+3, -6.1578555800761569d+3,
     : 1.4024024992482997d+6, -1.1035348452735992d+5,
     : -4.1693128846913914d+3, 4.6514181868262909d+4,
     : -5.1146251027461527d+4, 8.6428371791077952d+5,
     : 8.7460752687095180d+6, 1.3100449734468905d+6,
     : -2.0765677606361438d+4, -9.8083045296433338d+3,
     : -4.0324856948029071d+2, 4.2046293657669931d+5,
     : -2.4474189993584064d+3, 7.0401561974936604d+7,
     : -4.0807276724942763d+4, -1.8631379554415547d+4,
     : -1.6039386836978338d+4, -7.1256798878129475d+3,
     : -2.1966186391133757d+4, -1.9087927114946682d-1,
     : -3.9527426920095240d+4, -4.9124867807478659d+4,
     : -2.5665252570612144d+8, -7.5218418851912487d+4,
     : -8.9893217534461528d+3, -4.4969986566601438d+3,
     : -2.9242995792629285d+3, -6.0581254419086450d+4,
     : -6.7491706947863086d+3, -1.2281236110725298d+4,
     : -7.4993771260912717d+3, -5.7442069020237432d+3,
     : -3.3731269012112457d+3, -3.1241613123180632d+3,
     : -7.1113583230460470d+4, 3.5832811837989830d+3,
     : -1.5871389127806366d+4, -6.2982289832543911d+3,
     : -5.0200417084525388d+3, -3.2758834116815524d+3,
     : -4.0940655672390890d+3, -1.3023502805828647d+5,
     : -5.8445919321976624d+4, -2.7185830268766800d+4,
     : -4.2122766392374569d+3, -3.5075386423806481d+3,
     : -3.0469591310436608d+3, 5.9114028977943444d+3,
     : -1.2999224766935712d+5, 8.4619738995280154d+6,
     : -1.8243952815134086d+4, 2.3925255337428421d+3,
     : 2.4828446738476546d+3, -3.8546901393092094d+3,
     : -3.5486715423155026d+9/
       data c2/
     : 1.3869901713303808d+5, -2.5348844396165503d+4,
     : -7.2959553560753993d+5, 1.8502090239612642d+6,
     : -2.2134246026265659d+12, 5.6547156434241391d+13,
     : -4.0275820256380737d+6, 3.8766983062749277d+6,
     : -3.1049104565209582d+4, 3.4014183460467827d+5,
     : 4.7609853914598515d+5, -6.1004352738552541d+5,
     : 2.5348783307846780d+6, -2.5844292861118339d+5,
     : 3.8312894505392917d+5, -4.2842097811494328d+4,
     : -9.5099931252097001d+5, 4.5923458673621091d+5,
     : 5.7496928028322395d+5, -7.7292550516185569d+4,
     : -3.0299690378652373d+6, 9.6193841241419746d+5,
     : 1.2867554781721241d+6, -1.0153570617594868d+6,
     : 4.4686277064389380d+5, 5.0983700189128600d+5,
     : 5.2325390661720658d+5, 5.2468402718793391d+5,
     : 1.7959764423623160d+7, 1.1283809188972958d+6,
     : 5.0903126859089005d+5, 4.5053365021946700d+5,
     : 8.2060936585075004d+4, -2.6176182754719355d+6,
     : -3.2598138503205343d+6, 7.9640722250669390d+7,
     : 1.1670040005370220d+6, -1.4215448470909181d+6,
     : -2.8172861896951293d+5, -8.0473067097785068d+5,
     : -4.5137217095908383d+6, -3.7056058062351798d+6,
     : 9.9506151929210794d+8, -7.8821778193015037d+6,
     : -1.0940663569126553d+6, -9.0593270375349780d+5,
     : -1.7206795839091802d+6, -5.5420736002889853d+6,
     : -1.0953782271962038d+6, 1.2701539800558910d+9,
     : 3.2067978678007495d+7, -1.1754238936514689d+6,
     : -2.2207685213506715d+6, -4.3871833236271245d+6,
     : -1.0651295255202742d+7, -1.8396624247975864d+7,
     : 4.5835001926944768d+8, -5.9712127163287241d+5,
     : -3.5443404166936381d+6, -1.4564943352713587d+7,
     : -1.9020521766527168d+7, -1.4776280533292249d+7,
     : 3.5454855277465487d+9/
       data c3/
     : -3.3785237965957422d+4, 9.5992197672942214d-1,
     : -5.5250659813515255d+6, -3.9707885358504370d+1,
     : 7.3436809255811829d+11, -1.8343874809153504d+13,
     : 2.8581235030270815d+6, -1.6543468888034118d+2,
     : -1.4472401641426309d+0, 2.6595443000862858d+5,
     : 1.0955644669434284d+1, 1.0287497921459285d+6,
     : 6.3179003079547845d+6, -9.4722188833870717d+0,
     : -7.0625115454545085d+4, 2.7415326128105484d+0,
     : 1.0589063296790302d+6, 1.1381369742587799d+1,
     : 6.1334585417188539d+0, 5.0596158189844582d+3,
     : 2.3001356198663251d+6, -1.9039382075198565d+5,
     : 5.3573048535516614d+0, 1.0849803018162050d+6,
     : -1.3244336065724912d+1, 9.7655705272034055d-1,
     : -1.3160269694165722d+1, -9.3057670750093806d+4,
     : 4.8016322442295802d+2, 5.9879710713630885d+1,
     : -9.4644136276377569d+4, -1.2045556312768733d+1,
     : 3.4694072375902324d+5, 1.8197041515593089d+6,
     : 2.3008868199593192d+6, 2.0509315702337713d+3,
     : -5.0931972125955562d+0, 1.1278113959776843d+6,
     : 5.2047682162671973d+5, 8.3478348621548095d+5,
     : 2.6792808169032484d+6, 2.4235406811762964d+6,
     : 3.6886876292652967d+4, 4.7175039702662528d+6,
     : 8.1765651296439313d+5, 8.2526166718227847d+5,
     : 1.2517245450036945d+6, 3.0283524656486679d+6,
     : 1.0879031395537450d+6, -3.2681490071093994d+8,
     : -9.9291375236183349d+6, 6.5011623703442758d+5,
     : 1.3824758865078620d+6, 2.3186602261557719d+6,
     : 4.7722850974694099d+6, 7.3870852472487465d+6,
     : -1.1762796590575413d+8, 1.2504680736615297d+1,
     : 1.1967604890404730d+6, 5.9980774971887916d+6,
     : 7.3963313990726992d+6, 5.1965817819973333d+6,
     : 3.0481892384730451d+7/
       data c4/
     : -2.6006289883899549d+6, -1.6209865037555823d+6,
     : -9.9664810915562809d+6, -5.1987719947070526d+6,
     : 7.6900256287758583d+6, 2.8947291385665499d+6,
     : 1.1384716214299737d+6, 1.6908581361750853d+6,
     : -1.5501362056002012d+6, 1.5328564149576928d+6,
     : 2.8134501895610015d+6, 3.9451054989552922d+6,
     : -9.3037850414925776d+6, 1.0403462426789665d+6,
     : 2.8140730453682044d+6, -1.3282775423038101d+6,
     : -7.0012371530885610d+5, 2.6629529100111416d+6,
     : 2.8805581778235002d+6, -3.9999671790919481d+6,
     : 6.3897494097968296d+5, -3.6183580641298071d+6,
     : 9.3699986791882478d+5, 7.5182821885084512d+5,
     : 2.3376448490897086d+6, 2.0977952793206554d+6,
     : 2.0539362665502750d+6, 8.4497006153091940d+5,
     : -1.0668908841259938d+7, 1.9177704795680770d+6,
     : 6.2241133151068911d+5, 1.5031148034501523d+6,
     : 1.0426381640985068d+6, 4.0123454782332969d+5,
     : 5.0341148663270148d+5, -1.5891699966488680d+7,
     : 4.6964334723130558d+6, 8.6095189295240771d+5,
     : 8.1974878803733177d+5, 6.4012960994996876d+5,
     : 4.0573465802380879d+5, 5.2897534494037437d+5,
     : -1.3647062313072685d+7, 2.3207671129495803d+7,
     : 8.1655723475285002d+5, 5.1900346358843293d+5,
     : 4.5991186945639015d+5, 4.6240925399611820d+5,
     : 7.1180212424474873d+5, -1.8175580971626360d+7,
     : -1.1881197333086058d+7, 8.3675217496155575d+5,
     : 5.4909435392463487d+5, 5.5008194078319601d+5,
     : 5.2933303043585864d+5, 7.0273577812476754d+6,
     : -1.3260788115487826d+8, -7.1986511489418775d+7,
     : 5.7795548355659854d+5, 5.7413435235184180d+6,
     : 5.2652665911315624d+6, 1.1623279713695722d+6,
     : 3.7815257180233458d+6/
* this is 9x6 vandermode matrix, stored in column order. column n is 
* xhf^(6-1) with x=1-exp(-alphhf(r-rehf))
      data vander /

     : -3.3682184d-3, -1.0697538d-3, -1.5107259d-4, -7.7653518d-6,
     : 1.7015602e-10, 1.4034454d-5, 5.2199641d-4, 6.4251714d-3,
     : 3.6497013d-2,  1.0517686d-2, 4.2017198d-3, 8.7770495d-4,
     : 8.1682478d-5, 1.5299522d-8, 1.3114642d-4, 2.3666499d-3,
     : 1.7632252d-2, 7.0762938d-2, -3.2842798d-2, -1.6503282d-2,
     : -5.0993099d-3, -8.5920477d-4, 1.3756514d-6, 1.2255113d-3,
     : 1.0730019d-2, 4.8387238d-2, 1.3720009d-1,  1.0255577d-1,
     : 6.4820674d-2, 2.9626086d-2, 9.0378359d-3, 1.2369123d-4,
     : 1.1451918d-2, 4.8648226d-2, 1.3278649d-1, 2.6601304d-1,
     : -3.2024330d-1, -2.5459905d-1, -1.7212230d-1, -9.5067534d-2,
     : 1.1121656d-2, 1.0701363d-1, 2.2056343d-1, 3.6439880d-1,
     : 5.1576452d-1,  1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0/

      data alphhf / 0.68d0/
      data rehf / 1.736d0/
* this version uses all values of r, and fits with quintic polynomial in r-rej
      rmx=max(r1,r2)
      rmn=min(r1,r2)
* first determine hf two-body potential
      vtwo=vhf(rmn)
* now three-body potential from fit
      ipow=6
      ind=0
      alph=1.2d0
      do ir=1,9
         do iang=1,7
            ind=ind+1
            v(ir,iang)=vtwoexps(rmx,alph,xlam1(ind),xlam2(ind),
     :                 xlam3(ind),
     :           c1(ind),c2(ind),c3(ind),c4(ind),ipow)
          enddo
      enddo
* now interpolate over r
      x=one-exp(-alphhf*(rmn-rehf))
      do iang=1,7
          call dcopy(9,v(1,iang),1,cvec,1)
          do icol=1,6
             call dcopy(9,vander(1,icol),1,scmat(1,icol),1)
          enddo
* solve linear equations for expansion coefficients in x
cstart unix-ibm
c;        call dgells(0,cmat,8,v,9,xth,5,dummy,tol,8,5,4,nk,
c;     :       aux,28)
cend
cstart unix-hp unix-dec unix-aix unix-iris
          lwork=390
          call dgels('N',9,6,1,scmat,9,cvec,9,aux,
     :        lwork,info)
* to check results
*         do i=1,9
*            res=ddot(6,vander(i,1),9,cvec,1)
*            print *, i,res
*         enddo
          if (info.ne.0) then
             write (6,20) info
20           format (
     :       'info = ', i2, ' after return from dgels in vsig; abort')
             stop
          endif
*         print *, (cvec(i),i=1,6)
*         print *, polyval(cvec,x,5)
cend
* damping factor for short and long range
*        xmin = one-exp(-alphhf*(0.7d0-rehf))
         xmin = one-exp(-alphhf*(1.0d0-rehf))
         xmax = one-exp(-alphhf*(3.4d0-rehf))
         damp=0.25d0*(tanh(4d0*(x-xmin))+one)*
     :        (tanh(-4d0*(x-xmax))+one)
         call dscal(4,damp,cvec,1)
* add on two-body potential and hh dissociation energy
*        print *, polyval(cvec,x,5)
         vex(iang)=polyval(cvec,x,5)+vtwo+dehh
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
      if (n.eq.0) return
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
      if(iflag.le.1) then
        if(inifhh.eq.0) then
          call inifit_fhh
          inifhh=1
        end if
        call jacobi_n(rfh1,rh2,rfh2,rr,r,theta)
      end if
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
      data filnam/'fhhfit_1078_new.dat'/
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
      function func_rr(R,k,icase)
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (npes=4)
      common/cpar/ par_rr(10,npes),par_r(10,npes),par_theta(10,npes)
      data wthr/1.d-2/
c
      Re=par_rr(1,icase)
      Rdmp=par_rr(2,icase)
      wdmp=par_rr(3,icase)
      kpmx=nint(par_rr(4,icase))
      kimin=nint(par_rr(5,icase))
c
      if(wdmp.eq.0) wdmp=1.d0
      w=0.5d0*(1.d0+tanh((R-Rdmp)/wdmp))
      if(k.eq.1) then
        func_rr=1.d0
      else if(k.le.kpmx) then
        func_rr=(1.d0-w)*(R-Re)**(k-1)
      else 
        kk=k-kpmx-1+kimin
        rr=max(0.2d0,R)
        func_rr=w*10.d0**kk/rr**kk
      end if  
      return
      end
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
      subroutine jacobi_n(r1,r2,r3,rr,r,theta)
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
            stop 'arg.lt.0 in jacobi_n'
         endif
      endif
      rr=sqrt(arg) 
      if(abs(r2*RR).lt.1.d-6) then              ! avoids division by zero
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

      subroutine potfh2(rbond,vev,iflag)
* calculates diabatic composite F+H2 surfaces
*     
*      rbond(1)=rh2
*      rbond(2)=rfh1
*      rbond(3)=rfh2
*     
* iflag=0:  return 4 diabatic potentials
*      vev(1)=Vsig
*      vev(2)=Vpi= 1/2 [ Vpi(A") + Vpi(A') ]
*      vev(3)=v2 = 1/2 [ Vpi(A") - Vpi(A') ]
*      vev(4)=v1 = < Vsig | H | Vpi(A') > / sqrt(2)
* iflag=1:  return 3 adiabatic potentials
      
* The resulting energy values are in ev
*
      implicit double precision (a-h,o-z)
      dimension vv(91,81),vvsig(91,81),vvpi(91,81),vv2(91,81),
     :   vv1(91,81),rhh(91,81),evec(3,3),xmat(3,3)
      dimension rbond(3),vdau(4),vau(4),vdnew(4),eval(3),vev(4)
      dimension vex(7),vmer(7)
      data zero,half,one,two /0d0,0.5d0,1d0,2d0/
* conversion from hartree and ev to cm-1
      data tocm, toev / 219474.6d0, 8065.465d0/
* conversion from hartree to ev
      data hartoev /27.21165d0/
      bo=0.5d0
      pi=acos(-1d0)
      alpha_merge=2d0
      xline0=1.55d0
      xline1=2.25d0
      xline_alph=0.6d0
      rh=rbond(1)
      r1=rbond(2)
      r2=rbond(3)
* determine bond-angle coordinates
      th_bond=bondang(r1,r2,rh)
* determine diabatic entrance-arrangement fh2 potentials (guoliang li fit)
      call pot_fhh(rbond,vdau)
* construct damping functions for merge and to delineate entrance channel region
      xlineent=0.5d0*(tanh(4d0*(rh-2.6d0))+1d0)
      facent=0.5d0*(tanh(-6d0*(rr2-xlineent-2.9d0))+1d0)
* construct damping function
      call jacobi_n(r1,rh,r2,rrjac,rjac,thjac)
*     xlineent=0.6d0*rrjac-0.33333d0;
*     facent=0.5d0*(tanh(-4d0*(rjac-xlineent))+1d0);
      damp1=half*(tanh(5d0*(rrjac-3.5d0))+1d0)
      damp2=half*(1-tanh(4d0*(rr-3.2d0)))
* damp out coupling potentials if outside of entrance region
      vdau(3)=damp1*damp2*(1d0-facent)*vdau(3)
      vdau(4)=damp1*damp2*(1d0-facent)*vdau(4)
      call vadiab(vdau,vau,evec)
      rr1=max(r1,r2)
      rr2=min(r1,r2)
* determine entrance and exit channel potentials at 7 values of HFH bondangle
* determine lowest A' potential in exit arrangement from mha fit
      call vexit(rr1,rr2,vex)
      do ith=1,7
         th=(ith-1)*30d0
         rh=sqrt(rr1*rr1+rr2*rr2-2*rr1*rr2*cos(th*pi/180d0))
* determine lowest A' potential in interaction arrangement from mha fit
         call vint(rr1,rr2,rh,vin)
         fac=half*(one-cos(th*pi/180d0))
*        ro=6d0-fac
         ro=5-0.5d0*fac
         xline=1.55+2.25*0.5*(tanh(0.6*(rr1-ro))+1d0)
*        xline=1.65d0+1.05d0*0.5d0*(tanh(1.5d0*(rr1-ro))+1d0)+
*       :      0.2d0*0.5d0*(tanh(1d0*(rr1-8d0))+1d0)
         facexit=half*(tanh(-alpha_merge*(rr2-xline))+one)
         facint=half*(tanh(+alpha_merge*(rr2-xline))+one)
         vmer(ith)=vex(ith)*facexit+vin*facint
      enddo
* interpolate over angles
      call vinterp(vmer,vmer_intp,th_bond)
* vcomp2 is obtained by merging vint and vexit and then merging the result with vent
* with switching function facent (determined earlier in this subroutine)
      vcomp2=facent*vmer_intp+(1d0-facent)*vau(1)*tocm
* replace vau(1) with composite potential
*     print *, 'vmerge_interpolated', vmer_intp
      vau(1)=vcomp2/tocm
* redetermine diabatic potentials this is done only if facent < 0.999999
      if (facent.lt.0.999999d0) then
* transform corrected adiabatic potentials into diabatic basis         
         call vdiab_adiab(evec,vau,vdau)
         do i=1,4
           vev(i)=vdau(i)*hartoev
         enddo
      else
* here if were are far outside the entrance region
         vev(1)=vau(1)*hartoev
         do i=2,4
            vev(i)=vdau(i)*hartoev
         enddo
         do i=1,3
            eval(i)=vau(i)
         enddo         
      endif
      if (iflag.eq.1) then
        do i=1,3
          vev(i)=eval(i)*hartoev
        enddo
        vev(4)=0d0
      endif
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
* ------------------------------------------------------
       double precision function vtwoexps(r,alph,xlam1,xlam2,
     :  xlam3, c1,c2,c3,c4,ipow)

*  to determine
*    v(r)= c1*exp(-lam1*R) + (c2+c3*R)exp(-lam2*R)
*     -0.5[tanh{alpha(R-lam3)} +1]*c4/R^ipow

       implicit double precision (a-h,o-z)
       data zero,one,half/0.d0,1.d0,0.5d0/

       clr=c4
       poly=c2+c3*r
       ex=c1*dexp(-xlam1*r)+poly*dexp(-xlam2*r)
       xlr=-half*clr*(tanh(alph*(r-xlam3))+one)*r**(-ipow)
       vtwoexps=ex+xlr
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

c ---------------------------------------------------
      double precision function dampint(r1,r2,theta)
c ---------------------------------------------------

*     to determine damping function between interaction and exit regions
*     for fh2

*     on input:
*       r1, r2 are the two HF distances
*       theta is the HFH angle in degrees

*     on return dampint is the switching function, so that the merged
*     pes is  (1-damp).*vex+(damp).*vin

      implicit double precision (a-h,o-z)
      dimension vmat(2,2),rhs(2),aux(130), ccoef(6), dcoef(6),
     :   tmat(2,2)
      data ccoef /
     : -5.6584362d-10, 2.5042088d-7,-3.8685466d-5,
     : 2.4160354d-3,-4.4042929d-2, 2.0003247d0/
      data dcoef /
     : 4.8010974d-10,-2.2166105d-7,
     : 3.5044893d-5,-1.9823232d-3,-1.0202020d-3, 9.0051948d0/
      data vmat / 3d0, 1.2d0,1d0,1d0/
      cc=polyval(ccoef,theta,5)
      dd=polyval(dcoef,theta,5)
      rhs(1)=dd
      rhs(2)=cc
      lwork=130
      call dcopy(4,vmat,1,tmat,1)
      call dgels('N',2,2,1,tmat,2,rhs,2,aux,
     :        lwork,info)
      rmn=min(r1,r2)
      rmx=max(r1,r2) 
      rho=polyval(rhs,rmn,1)
      alph=1.25d0+0.025d0*(tanh(0.05d0*(theta-105d0))+1d0)
      dampint=0.5d0*(tanh(-alph*(rmx-rho))+1d0)
      return
      end


