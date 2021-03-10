*system:  F(2P)+H2, Dubernet-Hutson expansion of Alexander-stark-werner PES's
*references: M. H. Alexander, J. Chem. Phys. 99, 6014 (1993).
* M.-L. Dubernet and J. M. Hutson, J. Chem. Phys. 101, 1939 (1994).

      include "common/syusr"
      include "common/ground"
      include "common/bausr"
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      common /coselb/ ibasty
      potnam='ALEX-STARK-WERNER F(2P)H2 DUBERNET-HUTSON'
      ibasty=12
      lammin(1)=1
      lammax(1)=9
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
      potnam='ALEX-STARK-WERNER F(2P)H2 DUBERNET-HUTSON'
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
* latest revision date:  23-apr-1998
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
      subroutine pot(v00,r)
      implicit double precision (a-h,o-z)
      dimension b(4)
      common /covvl/ vvl(1)
* parameters for msv F+H2 potential
      data eps /4.12d0/
      data rm /3.34d0/
      data beta /6.3d0/
      data c0 /7.516d3/
      data x3 /1.11d0/
      data x4 /1.5d0/
      data a2 /8.1149d5/
      data alph2 /3.5d0/
      data c2 /9.020d+2/
      data b /-7.500735902794409d-1, 1.627465327650093d0,
-3.903934882318447d0, 
     : -3.723728951541869d0/
      clr=c0/(eps*rm**6)
      rang=r*0.52917715
      xx=rang/rm
      if (x.le.x3) then
         vv0=exp(-2d0*beta*(x-1))-2d0*exp(-beta*(x-1))
      elseif (x.gt. x3 .and. x.lt. x4) then
         vv0=b(1)+(x-x3)*(b(2)+(x-x4)*(b(3)+(x-x3)*b(4)))
      else
         vv0=-clr*x.**(-6)
      endif
      vvl(1)=-a2*exp(-alph2*rr)+c2*rr**(-6)
      
