* system:  Cl(2P)+H2, Dubernet-Hutson expansion of Capecchi-Werner r-dependent
* PES
* References:
*  G. Capecchi, H.-J. Werner Phys. Chem. Chem. Phys. 6, 4975 (2004)
*  M. H. Alexander  to be published
*  M.-L. Dubernet and J. M. Hutson, J. Chem. Phys. 101, 1939 (1994).
*  M. H. Alexander, J. Chem. Phys. 99, 6014 (1993).
*  M. H. Alexander and M. Yang, J. Chem. Phys. 103, 7956 (1995).


      include "common/syusr"
      include "common/ground"
      include "common/bausr"
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      common /coselb/ ibasty
      potnam='CAPECCHI WERNER-<v=0> CL(2P)H2 DUBERNET-HUTSON'
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
*     subroutine driver
      implicit double precision (a-h,o-z)
      character *48 potnam
      character *2 frame
      logical csflag, ljunk, ihomo, lljunk
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(9)
      potnam='CAPECCHI WERNER-<v=0> CL(2P)H2 DUBERNET-HUTSON'
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
      if (r .lt. 0.d0) go to 93
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
93    r=4
      do i=1,100
       call pot(vv0,r)
       write(2,101) r,vv0,vvl
101    format(f8.4,10(1pe16.8))
       r=r+0.2
      enddo

99    end
      subroutine pot (vv0, r)
* ----------------------------------------------------------------------
*  subroutine to calculate the r-dependent coefficients in the
*  cl-h2 potentials of capecchi and werner using body-frame
*  expansion of dubernet and flower (pes is averaged over v=0 vibrational
*  wavefunction of H2
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
* latest revision date:  1-feb-1998
* ----------------------------------------------------------------------


      implicit double precision (a-h,o-z)
      logical ljunk, ihomo, csflag, lljunk
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(9)
      dimension vxxl1(5), vxxl2(5), vxxr0(5), vxxc1(5), vxxc2(5),
     :  vxxc3(5), vxxcl(5)
      dimension vyyl1(5), vyyl2(5), vyyr0(5), vyyc1(5), vyyc2(5),
     :  vyyc3(5), vyycl(5)
      dimension vzzl1(5), vzzl2(5), vzzr0(5), vzzc1(5), vzzc2(5),
     :  vzzc3(5), vzzcl(5)
      dimension vxzlam1(3), vxzlam2(3), vxzc1(3), vxzc2(3)
      dimension vxzc3(3)
*      dimension beta(5)
      dimension vxx(5), vyy(5), vs(5),vzz(5),vd(5),vxz(5)
      dimension d0(25),d1(12),d2(16),aa(25)
      dimension xzz(5),xs(5),xd(4),x1(4),kpvt(5),
     :          qraux(5), work(24),rsd(5)
      dimension vvll(18)
*    vvll(1-5) expansion coefficients in dl0 (l=0,2,4,6,8) of vsigma
*    vvll(6-10) expansion coefficients in dl0 (l=0,2,4,6,8) of vsum
*    vvll(11-14) expansion coefficients in dl2 (l=2,4,6,8) of vdif
*    vvll(15-18) expansion coefficients in dl1 (l=2,4,6,8) of v1
* angles (not needed here, but included for completeness)
*      data beta / 0d0, 22.5d0, 45d0, 67.5d0, 90d0/
      data zero, one, half /0.d0,1.d0,0.5d0/
      data two, three, four, five /2.d0,3.d0,4.d0,5.d0/
      data seven,eleven /7.d0,11.d0/
      data sq2 /1.414213562d0/
* coefficients for d0 rotation matrices
* stored (by column) for each of 5 angles and for l=0,2,4,6,8
      data d0/
     : 1.0000000d0, 1.0000000d0, 1.0000000d0, 1.0000000d0, 1.0000000d0,
     : 1.0000000d0, 7.8033009d-01, 2.5000001d-01, -2.8033008d-01,
     : -5.0000000d-01, 1.0000000d0, 3.6159588d-01, -4.0625000d-01,
     : -8.0345887d-02, 3.7500000d-01, 1.0000000d0, -7.6358299d-02,
     : -1.4843752d-01, 2.7167082d-01, -3.1250000d-01, 1.0000000d0,
     : -3.5735359d-01, 2.9833984d-01, -2.7863272d-01, 2.7343750d-01/
* coefficients for d1 rotation matrices
* stored (by column) for each of 3 angles and for l=2,4,6,8
      data d1/
     : -4.3301270d-01, -6.1237244d-01, -4.3301270d-01,
     : -5.8796105d-01, -1.3975425d-01,  3.9031869d-01,
     : -4.9200540d-01,  3.5441551d-01, -1.8822068d-01,
     : -2.1012777d-01,  1.1186650d-01, -5.0614421d-02/
* coefficients for d2 rotation matrices
* stored (by column) for each of 4 angles and for l=2,4,6,8
      data d2/
     : 8.9679865d-02,  3.0618621d-01,  5.2269256d-01,  6.1237244d-01,
     : 2.8798601d-01,  4.9410589d-01,  8.4775334d-03, -3.9528471d-01,
     : 4.5386125d-01,  4.0027167d-02, -2.5372551d-01,  3.2021721d-01,
     : 4.8368899d-01, -3.2931303d-01,  2.8759680d-01, -2.7731624d-01/
* hyperbolic tangent scaling factor
      data alph /1.2d0/
*  expansion coefficients for vxx
*  by rows: lam1 lam2 r0 c1 c2 c3 cl
      data vxxl1/
     :  4.9209255d-1,4.5552850d-1,5.1674946d-1,5.3207154d-1,
     : 3.4029567d-1/
      data vxxl2/
     :  1.3249044d+0,1.3452039d+0,1.3724598d+0,1.3889515d+0,
     :  1.3799171d+0/
      data vxxr0/
     :  6.5665451d+0,7.6567676d+0,7.3058276d+0,6.2657593d+0,
     :  6.6208885d+0/
      data vxxc1/
     : -1.7594345d+3,-1.0582890d+3,-1.6418101d+3,-9.0280007d+2,
     :  -9.3615624d+1/
      data vxxc2/
     :  4.5969404d+6,4.8824826d+6,4.8568099d+6,4.8206524d+6,
     :  4.5974543d+6/
      data vxxc3/
     : -7.1786723d+5,-7.5962353d+5,-7.1960793d+5,-7.0305420d+5,
     :  -6.7948162d+5/
      data vxxcl/
     : -7.3907130d+5,3.7606042d+5,-6.6716941d+5,5.5407996d+5,
     :  6.5316725d+5/
*  expansion coefficients for vyy
*  by rows: lam1 lam2 r0 c1 c2 c3 cl
      data vyyl1/
     :  4.9208964d-1,4.6713594d-1,4.9691717d-1,4.9527905d-1,
     :  3.7192905d-1/
      data vyyl2/
     :  1.3248988d+0,1.3519813d+0,1.3543451d+0,1.3665110d+0,
     :  1.3648316d+0/
      data vyyr0/
     :  6.5665204d+0,5.7383688d+0,7.5460918d+0,5.9855370d+0,
     :  6.6140186d+0/
       data vyyc1/
     : -1.7593793d+3,-1.2293865d+3,-1.4603445d+3,-9.0373826d+2,
     :  -2.1705396d+2/
      data vyyc2/
     :  4.5968352d+6,4.8712496d+6,4.5221880d+6,4.3679675d+6,
     :  4.2900628d+6/
      data vyyc3/
     : -7.1785238d+5,-7.5007241d+5,-6.8430560d+5,-6.5275684d+5,
     :  -6.4891855d+5/
      data vyycl/
     : -7.3924070d+5,9.0072869d+5,-5.8100154d+5,4.0533927d+5,
     : 5.9305965d+5/
*  expansion coefficients for vzz
*  by rows: lam1 lam2 r0 c1 c2 c3 cl
      data vzzl1/
     :  9.3890115d-1,9.4484583d-1,8.3444589d-1,6.4959055d-1,
     :  7.2862927d-1/
      data vzzl2/
     :  1.9945466d+0,1.8110644d+0,1.9432847d+0,2.1637564d+0,
     :  1.8595869d+0/
      data vzzr0/
     :  8.1307783d+0,7.7509377d+0,7.5271951d+0,7.3819666d+0,
     :  6.6763854d+0/
      data vzzc1/
     : -4.3994319d+4,-6.3811033d+4,-3.1478116d+4,-9.8424586d+3,
     :  -2.0082570d+4/
      data vzzc2/
     : -2.2116574d+7,-6.3843306d+6,-9.2310121d+6,-1.4453759d+7,
     :  -9.9905192d+5/
      data vzzc3/
     :  7.3314878d+6,2.6078408d+6,3.7005850d+6,6.1549428d+6,
     :  1.1144046d+6/
      data vzzcl/
     : -1.9800586d+6,-1.5436109d+6,-5.2671060d+5,-3.1116318d+6,
     : -2.1887290d+6/
*  expansion coefficients for vxz
*  by rows: lam1 lam2 c1 c2 c3
      data vxzlam1/
     :  4.7786427d-1, 4.3287715d-1, 3.7791380d-1/
      data vxzlam2/
     :  1.3680462d+0, 1.3021489d+0, 1.2136737d+0/
      data vxzc1/
     :  4.1957637d+2, 3.5387705d+2, 1.2693966d+2/
      data vxzc2/
     : -7.2028240d+4, -6.8350883d+4, -2.7015677d+4/
      data vxzc3/
     :  2.4288785d+4, 2.4621284d+4,  1.0769822d+4/

    
* determine potentials at angles
      rm1=one/r
      rm2=rm1*rm1
      rm3=one/r**3
      rm5=rm3*rm1*rm1
      rm6=rm5*rm1
* quadal and quadh2 are quadrapole moments of Cl and H2
* averaged over v=0 vibrational wavefunction
* cm and rad are conversion factors from au to cm-1
* and from degrees to radians respectively
* this is currently not implemented here
      cm=2.194746d05
      pi=acos(-1.0d0)
      rad=pi/180.0d0
      quadcl=0*cm
      quadh2=4.8101d-01
      do 200 i=1,5
        if (i.eq.1) then
           ang=0.0d0
        elseif (i.eq.2) then
           ang=22.5d0*rad
        elseif (i.eq.3) then
           ang=45.0d0*rad
        elseif (i.eq.4) then
           ang=67.5d0*rad
        elseif (i.eq.5) then
           ang=90.0d0*rad
        endif
        cs2=(cos(ang))**2.0d0
        qzz=12.0D0*cs2-4.0D0
        qxx=3.0D0-7.0D0*cs2
        qyy=1.0D0-5.0D0*cs2
        alphq=0.6d0
        vxx(i)=vxxc1(i)*dexp(-vxxl1(i)*r)+
     :        (vxxc2(i)+vxxc3(i)*r)*dexp(-vxxl2(i)*r)-
     :        half*(tanh(alph*(r-vxxr0(i)))+one)*vxxcl(i)*rm6
        vyy(i)=vyyc1(i)*dexp(-vyyl1(i)*r)+
     :        (vyyc2(i)+vyyc3(i)*r)*dexp(-vyyl2(i)*r)-
     :        half*(tanh(alph*(r-vyyr0(i)))+one)*vyycl(i)*rm6
        vzz(i)=vzzc1(i)*dexp(-vzzl1(i)*r)+
     :        (vzzc2(i)+vzzc3(i)*r)*dexp(-vzzl2(i)*r)-
     :        half*(tanh(alph*(r-vzzr0(i)))+one)*vzzcl(i)*rm6
        if (i.eq.1 .or. i.eq.5) then
* vxz vanishes at 0 and 90 degrees
           vxz(i)=0.0d0
        else
           vxz(i)=vxzc1(i-1)*dexp(-vxzlam1(i-1)*r)
     :           +(vxzc2(i-1)+vxzc3(i-1)*r)*dexp(-vxzlam2(i-1)*r)
        endif
* include sqrt(2) (table I of mh alexander, jcp 99, 6014)
        vxz(i)=sq2*vxz(i)
        vs(i)=0.5d0*(vxx(i)+vyy(i))
        vd(i)=0.5d0*(vyy(i)-vxx(i))
200   continue
      sum=0
* determine matrix of d's at angles for least-squares fit
*      call d2matev(9,beta,d0,d1,d2)
* (n.b. these are already in common, subroutine is included
*  only for completeness)
* solve simultaneous equations for solutions
* first for vzz
      tol=1.e-10
      call dcopy(25,d0,1,aa,1)
      call dqrank(aa,5,5,5,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,5,5,5,kr,vzz,xzz,rsd,kpvt,qraux)
      call dcopy(5,xzz,1,vvll,1)
      call dqrlss(aa,5,5,5,kr,vs,xs,rsd,kpvt,qraux)
      call dcopy(5,xs,1,vvll(6),1)
*
      call dcopy(16,d2,1,aa,1)
      call dqrank(aa,4,4,4,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,4,4,4,kr,vd(2),xd,rsd,kpvt,qraux)
      call dcopy(4,xd,1,vvll(11),1)
*
      call dcopy(12,d1,1,aa,1)
      call dqrank(aa,3,3,4,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,3,3,4,kr,vxz(2),x1,rsd,kpvt,qraux)
      call dcopy(4,x1,1,vvll(15),1)
* determine body frame expansion coefficients in dubernet and flower
* expansion
* here is totally symmetric term
      vv0=(two*vvll(6)+vvll(1))/three
      sq6=1.d0/sqrt(6.d0)
* v200
      vvl(1)=(two*vvll(7)+vvll(2))/three
* v400
      vvl(6)=(two*vvll(8)+vvll(3))/three
* v020
      vvl(2)=(vvll(1)-vvll(6))*five/three
* v220
      vvl(3)=(vvll(2)-vvll(7))*five/three
* v420
      vvl(7)=(vvll(3)-vvll(8))*five/three
* v222
      vvl(4)=five*vvll(11)*sq6
* v422
      vvl(9)=five*vvll(12)*sq6
* v221
      vvl(5)=five*vvll(15)*sq6
* v421
      vvl(8)=five*vvll(16)*sq6
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
