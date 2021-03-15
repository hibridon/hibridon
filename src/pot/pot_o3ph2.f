*system:  O(3P)+H2, Dubernet-Hutson expansion of scaled Alexander PES's
*references:
*  M. H. Alexander, J. Chem. Phys. 108, 4467 (1998).
*  M.-L. Dubernet and J. M. Hutson, J. Chem. Phys. 101, 1939 (1994).
*  M. H. Alexander, J. Chem. Phys. 99, 6014 (1993).
*  M. H. Alexander and M. Yang, J. Chem. Phys. 103, 7956 (1995).


* NB potential fitted to all ab initio mrci scaled points determined
* with avqz basis with the restrictions:  R .ge. 4 and E .le. 2000 cm-1
* rms fit to these points is 0.0011 (relative error)
      include "common/syusr"
      include "common/ground"
      include "common/bausr"
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      common /coselb/ ibasty
      potnam='ALEXANDER SCALED (s=1.21) O(3P)H2 DUBERNET-HUTSON'
      ibasty=13
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
      econv=219474.6d0
      potnam='ALEXANDER SCALED (s=1.21) O(3P)H2 DUBERNET-HUTSON'
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
      if (.not. csflag .or. (csflag .and. ihomo)) 
     +    write (6, 100) vv0*econv,vvl*econv
100   format(' V000, V220, V022, V202:  ',4(1pe16.8),/,
     :       ' V222, V224, V404:  ',3(1pe16.8),/,
     :       ' V422, V424, V426:  ',3(1pe16.8))
      if (csflag .and. .not.ihomo) write (6, 110) 
     +    vv0*econv,vvl*econv
110   format(' v000, v220, v020, v200:  ',4(1pe16.8),/,
     :       ' v222, v221, v400:  ',3(1pe16.8),/,
     :       ' v420, v422, v421:  ',3(1pe16.8))
      goto 1
93    r=4
      open (unit=2,file='o3ph2_hib_1998_sf_vlm.txt')
      write (2,102)
102   format(' %R/bohr 1V000 1V220 1V022 1V202 1V222 1V224',
     +  ' 1V404 1V422 1V424 1V426')
      do i=1,150
       call pot(vv0,r)
       write(2,101) r,vv0*econv,vvl*econv
101    format(f8.4,10(1pe16.8))
       r=r+0.2
      enddo
      close(2)
99    end
      subroutine pot (vv0, r)
* ----------------------------------------------------------------------
*  subroutine to calculate the r-dependent coefficients in the
*  b-n2 potentials of yang and alexander using body-frame
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
* latest revision date:  2-feb-1998
* potential from fit to matlab oh2pln_2state_lz.m
* coefficients of fit are coefv11s (vxx), coefv22s (vzz), coefvaps (vyy)
*    and coefv12s (vxy)
* ----------------------------------------------------------------------


      implicit double precision (a-h,o-z)
      logical ljunk, ihomo, csflag, lljunk
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(9)
      dimension vxxl1(9), vxxl2(9), vxxr0(9), vxxc1(9), vxxc2(9),
     :  vxxc3(9), vxxcl(9)
      dimension vyyl1(9), vyyl2(9), vyyr0(9), vyyc1(9), vyyc2(9),
     :  vyyc3(9), vyycl(9)
      dimension vzzl1(9), vzzl2(9), vzzr0(9), vzzc1(9), vzzc2(9),
     :  vzzc3(9), vzzcl(9)
      dimension vxzl1(7), vxzl2(7), vxzc1(7), vxzc2(7)
*      dimension beta(9)
      dimension vxx(9), vyy(9), vs(9),vzz(9),vd(9),vxz(9),
     :          v1ap(9),v2ap(9)
      dimension d0(45),d1(36),d2(36),aa(45)
      dimension fthet(3,9)
      dimension xzz(5),xs(5),x1ap(5),x2ap(5),xd(4),x1(4),kpvt(9),
     :          qraux(9), work(24),rsd(9)
      dimension vvll(28)
*    vvll(1-5) expansion coefficients in dl0 (l=0,2,4,6,8) of vsigma
*    vvll(6-10) expansion coefficients in dl0 (l=0,2,4,6,8) of vsum
*    vvll(11-14) expansion coefficients in dl2 (l=2,4,6,8) of vdif
*    vvll(15-18) expansion coefficients in dl1 (l=2,4,6,8) of v1
*    vvll(19-23) expansion coefficients in dl0 (l=0,2,4,6,8) of v1A'
*    vvll(24-28) expansion coefficients in dl0 (l=0,2,4,6,8) of v1A'
* angles (not needed here, but included for completeness)
*      data beta /0d0, 11.25d0, 22.5d0, 33.75d0, 45.d0, 56.25d0, 67.5d0,
*                 78.75d0, 90.d0/

      data zero, one, half /0.d0,1.d0,0.5d0/
      data two, three, four, five /2.d0,3.d0,4.d0,5.d0/
      data seven,eleven /7.d0,11.d0/
      data sq2 /1.414213562d0/
* coefficients for d0 rotation matrices
* stored (by column) for each of 9 angles and for l=0,2,4,6,8
      data d0/
     : 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 9.4290965d-1,
     : 7.8033009d-1, 5.3701258d-1, 2.5d-1, -3.7012563d-2, -2.8033008d-1,
     : -4.4290964d-1, -5.0000d-1, 1d0, 8.1603638d-1, 3.6159588d-1,
     : -1.2648544d-1, -4.0625d-1, -3.6566260d-1,
     :  -8.0345887d-2, 2.3861165d-1,
     : 3.7500d-1, 1d0, 6.3379430d-1, -7.6358299d-2, -4.1470677d-1,
     : -1.4843752d-1, 2.6199014d-1, 2.7167082d-1, -9.0452652d-2,
     :  -3.1250d-1,
     : 1d0, 4.1666541d-1, -3.5735359d-1, -1.7953446d-1, 2.9833984d-1,
     : 8.9800577d-2, -2.7863272d-1, -2.7859278d-2, 2.7343750d-1/
* coefficients for d1 rotation matrices
* stored (by column) for each of 9 angles and for l=2,4,6,8
      data d1/
     : 0d0,-2.3434479d-1,-4.3301270d-1,-5.6575836d-1,-6.1237244d-1,
     : -5.6575836d-1,-4.3301270d-1,-2.3434479d-1,0d0,0d0,
     : -3.9935575d-1,
     : -5.8796105d-1,-4.7499021d-1,-1.3975425d-1,2.1675803d-1,
     :  3.9031869d-1,
     : 2.9239248d-1,0d0,0d0,-5.1753174d-1,-4.9200540d-1,
     : -6.0266556d-3,
     : 3.5441551d-1,2.0878158d-1,-1.8822068d-1,-3.0272350d-1,0d0,0d0,
     : -5.7741337d-1,-2.1012777d-1,3.3709181d-1,1.1186650d-1,
     :  -2.9071290d-1,
     : -5.0614421d-2, 2.7597728d-1, 0d0/
* coefficients for d2 rotation matrices
* stored (by column) for each of 8 angles and for l=2,4,6,8
      data d2/
     : 0d0,2.3307038d-2,8.9679867d-2,1.8901383d-1,3.0618622d-1,
     : 4.2335861d-1,
     : 5.2269257d-1,5.8906540d-1,6.1237244d-1,0d0,8.6259556d-2,
     : 2.8798601d-1,
     : 4.6843615d-1,4.9410588d-1,3.1716725d-1,8.4775167d-3,
     : -2.7893590d-1,
     : -3.9528471d-1,0d0,1.7331796d-1,4.5386126d-1,4.2780021d-1,
     : 4.0027151d-2,
     : -3.1257486d-1,-2.5372550d-1,1.1172821d-1,3.2021721d-1,0d0,
     : 2.7133875d-1,4.8368899d-1,6.1484133d-2,-3.2931303d-1,
     : -4.4640198d-2,
     : 2.8759680d-1, 1.5131952d-2, -2.7731624d-1/
* hyperbolic tangent scaling factor
      data alph /1.2d0/
*  expansion coefficients for vxx
*  by rows: lam1 lam2 r0 c1 c2 c3 cl
      data vxxl1/
     : 1.0009323d0, 9.6435044d-1, 8.6727802d-1, 7.6961256d-1,
     : 7.1100717d-1, 6.7378981d-1, 6.5443595d-1, 6.4131486d-1,
     : 6.3680814d-1/
      data vxxl2/
     : 2.2958728d0, 2.3196989d0, 2.3791880d0, 2.4686627d0, 2.5330163d0,
     : 2.6029127d0, 2.6546652d0, 2.7171115d0, 2.7380774d0/
       data vxxr0/
     : 6.2741230d0, 6.3646212d0, 6.5238134d0, 6.7100984d0, 6.7421066d0,
     : 6.7719166d0, 6.7572671d0, 6.7813085d0, 6.7889758d0/
      data vxxc1/
     : -2.3349853d4, -1.8476771d4, -1.0332806d4, -6.0161871d3,
     : -4.6875659d3, -4.1215917d3, -3.9682107d3, -3.8213805d3,
     : -3.7723711d3/
       data vxxc2/
     : -2.8051620d7, -3.0507290d7, -3.7292260d7, -5.2471062d7,
     : -6.3301872d7, -7.9142856d7, -9.1418388d7, -1.1807198d8,
     : -1.2819986d8/
      data vxxc3/
     : 1.1491992d7, 1.2367283d7, 1.4791489d7, 2.0083937d7, 2.3898246d7,
     : 2.9384205d7, 3.3670273d7, 4.2551504d7, 4.5916662d7/
      data vxxcl/
     : -2.2943544d4, -2.4612627d4, -4.0518782d4, -6.8082608d4,
     : -9.6986294d4, -1.2277418d5, -1.3917957d5, -1.5008960d5,
     : -1.5432418d5/
*  expansion coefficients for vzz
*  by rows: lam1 lam2 r0 c1 c2 c3 cl
      data vzzl1/
     : 5.9872573d-1, 6.0056003d-1, 6.0483487d-1, 6.1548146d-1,
     : 6.3653507d-1, 6.8038073d-1, 7.6913713d-1, 9.6167395d-1,
     : 1.0653241d0/
      data vzzl2/
     : 2.4303626d0, 2.4161642d0, 2.3897642d0, 2.3468610d0, 2.2883743d0,
     : 2.2251604d0, 2.1776059d0, 2.0547793d0, 1.9551856d0/
        data vzzr0/
     : 7.3407810d0, 7.3434237d0, 7.3587963d0, 7.3608507d0, 7.3163311d0,
     : 7.2383750d0, 7.2005983d0, 6.7613598d0, 6.3511642d0/
      data vzzc1/
     : -4.8067992d3, -4.7060428d3, -4.3704401d3, -3.9748707d3,
     : -3.6737434d3, -3.7849873d3, -5.1470033d3, -1.7326603d4,
     : -3.9484695d4/
      data vzzc2/
     : -1.1593870d8, -1.0453941d8, -8.6629358d7, -6.3166580d7,
     : -4.0184668d7, -2.3733405d7, -1.6215832d7, -4.1592267d6,
     : 1.5712474d5/
      data vzzc3/
     : 4.3337521d7, 3.9696623d7, 3.3811939d7, 2.5978123d7, 1.8061750d7,
     : 1.2127428d7, 9.1591077d6, 4.3049584d6, 2.2753008d6/
      data vzzcl/
     : -2.0546085d5, -1.9636213d5, -1.7322281d5, -1.4227678d5,
     : -1.1125370d5, -7.9207108d4, -4.1255144d4, -1.5716731d4,
     :  -2.3294351d4/
*  expansion coefficients for vyy
*  by rows: lam1 lam2 r0 c1 c2 c3 cl
      data vyyl1/
     : 1.0865565d0, 1.0667884d0, 1.0155165d0, 8.1296327d-1,
     : 7.4561122d-1, 7.0528426d-1, 6.8017117d-1, 6.6616900d-1,
     : 6.6277585d-1/
       data vyyl2/
     : 2.2684335d0, 2.2764030d0, 2.2968253d0, 2.4175185d0, 2.4681328d0,
     : 2.5137021d0, 2.5590431d0, 2.5937163d0, 2.6082100d0/
       data vyyr0/
     : 7.5696414d0, 7.7267870d0, 8.1461826d0, 6.6375099d0, 6.6976257d0,
     : 6.7301150d0, 6.7505033d0, 6.7623844d0, 6.7608269d0/
      data vyyc1/
     : -4.0856411d4, -3.6013915d4, -2.6282475d4, -7.1750053d3,
     : -5.0306392d3, -4.1899015d3, -3.7800263d3, -3.5836767d3,
     : -3.5462761d3/
      data vyyc2/
     : -2.5731026d7, -2.5806927d7, -2.5872521d7, -4.0571215d7,
     : -4.5767873d7, -5.0867850d7, -5.8076306d7, -6.5065599d7,
     : -6.8831433d7/
      data vyyc3/
     : 1.0700507d7, 1.0753017d7, 1.0855280d7, 1.6197798d7, 1.8232981d7,
     : 2.0235663d7, 2.2937587d7, 2.5502066d7, 2.6842083d7/
       data vyycl/
     : 7.6329132d3, 1.1270935d4, 2.2870774d4, -5.2028198d4,
     : -7.7160963d4, -9.5677231d4, -1.0771064d5, -1.1526628d5,
     : -1.1585473d5/
*  expansion coefficients for vxz
*  by rows:  lam1 lam2 c1 c2
      data vxzl1/
     : 5.7120349d-1, 5.9112002d-1, 5.9214983d-1, 5.9478707d-1,
     : 5.8086918d-1, 5.7659376d-1, 5.7662435d-1/
      data vxzl2/
     : 1.3174359d0, 1.3513964d0, 1.3601519d0, 1.3628915d0, 1.3198218d0,
     : 1.2943185d0, 1.2958705d0/
      data vxzc1/
     : 3.0966248d2, 6.7852384d2, 8.9040322d2, 9.7177046d2, 7.8236126d2,
     : 5.6599147d2, 3.0510497d2/
      data vxzc2/
     : 1.3281398d4, 2.5481720d4, 3.2043958d4, 3.2161026d4, 2.3927064d4,
     : 1.5607162d4, 8.1097137d3/


* determine potentials at angles
      rm3=one/r**3
      rm1=one/r
      rm5=rm3*rm1*rm1
      rm6=rm3*rm3
      do 200 i=1,9
        vxx(i)=vxxc1(i)*exp(-vxxl1(i)*r)+
     :        (vxxc2(i)+vxxc3(i)*r)*exp(-vxxl2(i)*r)-
     :         half*(tanh(alph*(r-vxxr0(i)))+1)*vxxcl(i)*rm5
        vyy(i)=vyyc1(i)*exp(-vyyl1(i)*r)+
     :        (vyyc2(i)+vyyc3(i)*r)*exp(-vyyl2(i)*r)-
     :         half*(tanh(alph*(r-vyyr0(i)))+1)*vyycl(i)*rm5

* here for double exponential plus tanh for vzz
        vzz(i)=vzzc1(i)*exp(-vzzl1(i)*r)+
     :         (vzzc2(i)+vzzc3(i)*r)*exp(-vzzl2(i)*r)-
     :          half*(tanh(alph*(r-vzzr0(i)))+1)*vzzcl(i)*rm5
        if (i .eq. 1 .or. i .eq. 9) then
* vxz vanishes at 0 and 90 degrees
           vxz(i)=0d0
        else
* double exponential at other angles
           vxz(i)=vxzc1(i-1)*exp(-vxzl1(i-1)*r)+
     :        vxzc2(i-1)*exp(-vxzl2(i-1)*r)
        endif
* include sqrt(2) (table I of mh alexander, jcp 99, 6014)
        vxz(i)=sq2*vxz(i)
        vs(i)=0.5d0*(vxx(i)+vyy(i))
        vd(i)=0.5d0*(vyy(i)-vxx(i))
* determine adiabatic A" potentials
        h11=vxx(i)
        h22=vzz(i)
        if (i .eq. 1) then
          v1ap(i)=min(h11,h22)
          v2ap(i)=max(h11,h22)
        else
          h12=vxz(i)
          ang=0.5d0*atan(2*h12/(h11-h22))
          cs=cos(ang)
          sn=sin(ang)
          e1=h11*cs**2+2*cs*sn*h12+h22*sn**2
          e2=h22*cs**2-2*cs*sn*h12+h11*sn**2
          v1ap(i)=min(e1,e2)
          v2ap(i)=max(e1,e2)
        endif
200   continue
      sum=0
* determine matrix of d's at angles for least-squares fit
*      call d2matev(9,beta,d0,d1,d2)
* (n.b. these are already in common, subroutine is included
*  only for completeness)
* solve simultaneous equations for solutions
* first for vzz
      tol=1.e-10
      call dcopy(45,d0,1,aa,1)
      call dqrank(aa,9,9,5,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,9,9,5,kr,vzz,xzz,rsd,kpvt,qraux)
      call dcopy(5,xzz,1,vvll,1)
      call dqrlss(aa,9,9,5,kr,vs,xs,rsd,kpvt,qraux)
      call dcopy(5,xs,1,vvll(6),1)
      call dqrlss(aa,9,9,5,kr,v1ap,x1ap,rsd,kpvt,qraux)
      call dcopy(5,x1ap,1,vvll(19),1)
      call dqrlss(aa,9,9,5,kr,v2ap,x2ap,rsd,kpvt,qraux)
      call dcopy(5,x2ap,1,vvll(24),1)
      call dcopy(36,d2,1,aa,1)
      call dqrank(aa,9,9,4,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,9,9,4,kr,vd,xd,rsd,kpvt,qraux)
      call dcopy(4,xd,1,vvll(11),1)
      call dcopy(36,d1,1,aa,1)
      call dqrank(aa,9,9,4,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,9,9,4,kr,vxz,x1,rsd,kpvt,qraux)
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
