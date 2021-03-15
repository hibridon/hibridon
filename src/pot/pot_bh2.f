* system:  B(2P)+H2, Dubernet-Hutson expansion of Williams-Alexander PES's
* references: J. Williams and M. H. Alexander, J. Chem. Phys. 112, 5722 (2000).
*   M.-L. Dubernet and J. M. Hutson, J. Chem. Phys. 101, 1939 (1994).
*   M. H. Alexander, J. Chem. Phys. 99, 6014 (1993).
*   M. H. Alexander and M. Yang, J. Chem. Phys. 103, 7956 (1995).


      include "common/syusr"
      include "common/ground"
      include "common/bausr"
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      common /coselb/ ibasty
      potnam='WILLIAMS-ALEXANDER CBS B(2P)H2 DUBERNET-HUTSON'
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
      potnam='WILLIAMS-ALEXANDER CBS B(2P)H2 DUBERNET-HUTSON'
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
*  al-h2 potentials of williams and alexander using body-frame
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
      dimension vxzc1(3), vxzc2(3), vxzc3(3), vxzc4(3)
*      dimension beta(5)
      dimension vxx(5), vyy(5), vs(5),vzz(5),vd(5),vxz(5)
      dimension d0(25),d1(12),d2(16),aa(25)
      dimension fthet(3,5)
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
      data alph /0.6d0/
*  expansion coefficients for vxx
*  by rows: lam1 lam2 r0 c1 c2 c3 cl
      data vxxl1/
     : 4.4670384d-01, 4.3742090d-01, 4.3935411d-01, 5.0692289d-01,
     : 5.1380907d-01/
      data vxxl2/
     : 1.6996970d0, 2.0244423d0, 2.0414867d0, 1.8258795d0,
     : 1.8166005d0/
      data vxxr0/
     : 5.6466938d0, 6.0048944d0, 5.7095172d0, 5.0741013d0,
     : 4.9807616d0/
      data vxxc1/
     : 4.0251766d02, 3.0364868d02, 3.0179821d02, 9.5473025d02,
     : 1.0791195d03/
      data vxxc2/
     : -1.6466810d05,  -7.2475343d06,  -7.5872480d06,  -2.8153405d06,
     : -2.5974526d06/
      data vxxc3/
     : 5.6615293d05, 3.4021097d06, 3.2167559d06, 1.2289871d06,
     : 1.1193141d06/
      data vxxcl/
     : 1.4469377d07, 1.3454264d07, 1.2419143d07, 1.3932886d07,
     : 1.3820838d07/
*  expansion coefficients for vyy
*  by rows: lam1 lam2 r0 c1 c2 c3 cl
      data vyyl1/
     : 4.4668826d-01, 4.6918726d-01, 4.3993516d-01, 4.3808500d-01,
     : 4.4929610d-01/
      data vyyl2/
     : 1.6996795d0, 1.6045513d0, 1.7595547d0, 1.8351076d0,
     : 1.8940367d0/
      data vyyr0/
     : 5.6466669d0, 5.3827504d0, 5.2804372d0, 5.1322688d0,
     : 5.1201533d0/
      data vyyc1/
     : 4.0244101d02, 6.1667676d02, 2.9959718d02, 2.5670936d02,
     : 2.7201522d02/
      data vyyc2/
     : -1.6456278d05, 1.5729801d05, -1.1858770d06,-2.2195718d06,
     : -3.1764182d06/
      data vyyc3/
     : 5.6609358d05, 3.2573659d05, 9.0031874d05, 1.3011085d06,
     : 1.6908730d06/
      data vyycl/
     : 1.4469237d07, 1.5108960d07, 1.2273124d07, 1.1038818d07,
     : 1.0519095d07/
*  expansion coefficients for vzz
*  by rows: lam1 lam2 r0 c1 c2 c3 cl
      data vzzl1/
     : 4.3245145d-01, 4.3176206d-01, 4.3839961d-01, 6.7618076d-01,
     : 6.6245347d-01/
      data vzzl2/
     : 1.3087275d0, 1.2819667d0, 1.4108705d0, 1.6812086d0,
     : 1.6760359d0/
      data vzzr0/
     : 4.9347684d0, 5.0132111d0, 4.8912907d0, 6.3665275d0,
     : 6.2886498d0/
      data vzzc1/
     : 2.1842115d03, 2.1682465d03, 1.4635554d03,  -4.5538411d03,
     : -3.1826428d03/
      data vzzc2/
     : 1.4572167d05, 3.2010430d05,  -1.1521935d05,  -3.3199021d06,
     : -2.8739302d06/
      data vzzc3/
     : 3.3008821d05, 2.5187107d05, 5.9109906d05, 2.2987595d06,
     : 2.1882520d06/
      data vzzcl/
     : 4.7922930d07, 4.8061186d07, 3.3297844d07, 5.7327542d06,
     : 6.6751891d06/
*  expansion coefficients for vxz
*  by rows: c1 c2 c3 c4
      data vxzc1/
     : 5.2063459d02, 5.1853566d02, 5.2608412d02/
      data vxzc2/
     : -3.8463714d02,  -3.8028904d02,  -3.8383943d02/
      data vxzc3/
     : 1.1403683d02, 1.1222302d02, 1.1236680d02/
      data vxzc4/
     : -6.5076944d0,  -6.0216420d0,  -6.3920887d0/
* determine potentials at angles
      rm1=one/r
      rm2=rm1*rm1
      rm3=one/r**3
      rm5=rm3*rm1*rm1
      rm6=rm5*rm1
* quadb and quadh2 are quadrapole moments of B and H2
* cm and rad are conversion factors from au to cm-1
* and from degrees to radians respectively
      cm=2.194746d05
      pi=acos(-1.0d0)
      rad=pi/180.0d0
      quadb=-2.3305d0*cm
      quadh2=4.5455d-01
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
        qzz=12.0d0*cs2-4.0d0
        qxx=3.0d0-7.0d0*cs2
        qyy=1.0d0-5.0d0*cs2
        vxx(i)=vxxc1(i)*dexp(-vxxl1(i)*r)+
     :        (vxxc2(i)+vxxc3(i)*r)*dexp(-vxxl2(i)*r)-
     :        half*(tanh(alph*(r-vxxr0(i)))+one)*vxxcl(i)*rm6
     :        +0.75d0*quadh2*quadb*qxx*half*(tanh(alph*
     :        (r-6.0d0))+one)*rm5
        vyy(i)=vyyc1(i)*dexp(-vyyl1(i)*r)+
     :        (vyyc2(i)+vyyc3(i)*r)*dexp(-vyyl2(i)*r)-
     :        half*(tanh(alph*(r-vyyr0(i)))+one)*vyycl(i)*rm6
     :        +0.75d0*quadh2*quadb*qyy*half*(tanh(alph*
     :        (r-6.0d0))+one)*rm5
        vzz(i)=vzzc1(i)*dexp(-vzzl1(i)*r)+
     :        (vzzc2(i)+vzzc3(i)*r)*dexp(-vzzl2(i)*r)-
     :        half*(tanh(alph*(r-vzzr0(i)))+one)*vzzcl(i)*rm6
     :        +0.75d0*quadh2*quadb*qzz*half*(tanh(alph*
     :        (r-6.0d0))+one)*rm5
        if (i.eq.1 .or. i.eq.5) then
* vxz vanishes at 0 and 90 degrees
           vxz(i)=0.0d0
        else
           vxz(i)=dexp(vxzc4(i-1)+vxzc3(i-1)*rm1+vxzc2(i-1)*rm2
     :            +vxzc1(i-1)*rm3)
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
*----------------Jason's testing
      junk=3.0d0
*---------------------
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
