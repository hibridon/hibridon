*system:  B(2P)+oH2, Dubernet-Hutson expansion of Alexander PES's
*references: M. H. Alexander, J. Chem. Phys. 99, 6014 (1993).
* M.-L. Dubernet and J. M. Hutson, J. Chem. Phys. 101, 1939 (1994).
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
      potnam='ALEXANDER B(2P)H2(J=0,1) DUBERNET-HUTSON'
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
      potnam='ALEXANDER B(2P)H2(J=0,1) DUBERNET-HUTSON'
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
*  b-h2 potentials of alexander, using body-frame
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
* latest revision date:  11-jun-1995
* ----------------------------------------------------------------------


      implicit double precision (a-h,o-z)
      logical ljunk, ihomo, csflag, lljunk
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(9)
      dimension vsl1(9), vsl2(9), vsr0(9), vsc1(9), vsc2(9),
     :  vsc3(9), vscl(9)
      dimension vzzl1(9), vzzl2(9), vzzr0(9), vzzc1(9), vzzc2(9),
     :  vzzc3(9), vzzcl(9)
      dimension vdl1(9), vdl2(9), vdr0(9), vdc1(9), vdc2(9),
     :  vdc3(9), vdcl(9)
      dimension vxzc1(9), vxzc2(9), vxzc3(9), vxzc4(9)
*      dimension beta(9)
      dimension vs(9),vzz(9),vd(9),vxz(9),v1ap(9),v2ap(9)
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
*      data beta / 90d0, 85d0, 75d0, 60d0, 45d0, 30d0, 15d0, 5d0, 1d0/
      data zero, one, half /0.d0,1.d0,0.5d0/
      data two, three, four, five /2.d0,3.d0,4.d0,5.d0/
      data seven,eleven /7.d0,11.d0/
* coefficients for d0 rotation matrices
* stored (by column) for each of 9 angles and for l=0,2,4,6,8
      data d0/1.d+0, 1.d+0, 1.d+0, 1.d+0, 1.d+0, 1.d+0, 1.d+0, 1.d+0,
     : 1.d+0, -5.d-1, -4.8860581d-1, -3.9951904d-1, -1.2499999d-1,
     : 2.5000001d-1, 6.2500001d-1, 8.9951905d-1, 9.8860581d-1,
     : 9.9954312d-1,3.7500000d-1, 3.4676697d-1, 1.4342954d-1,
     : -2.8906251d-1, -4.0625000d-1, 2.3437511d-2, 6.8469544d-1,
     : 9.6227183d-1, 9.9847747d-1,-3.1250000d-1, -2.6378009d-1,
     : 4.3100282d-2, 3.2324218d-1, -1.4843752d-1, -3.7402343d-1,
     : 3.9830600d-1, 9.2159756d-1, 9.9680403d-1,2.7343750d-1,
     : 2.0174615d-1, -1.7021999d-1, -7.3638895d-2, 2.9833984d-1,
     : -3.3877564d-1, 9.6184338d-2, 8.6750724d-1,
     : 9.9452433d-1/
* coefficients for d1 rotation matrices
* stored (by column) for each of 8 angles and for l=2,4,6,8
      data d1/
     :  0.d0, -1.0633736d-1, -3.0618622d-1,
     : -5.3033009d-1, -6.1237244d-1, -5.3033009d-1, -3.0618622d-1,
     : -1.0633736d-1, -2.1371490d-2, 0.d0, 1.4302762d-1,
     : 3.5373043d-1, 3.0257682d-1, -1.3975425d-1, -5.4463828d-1,
     : -4.9348468d-1, -1.9156376d-1, -3.8998025d-2, 0.d0,
     : -1.6789168d-1, -3.1780559d-1, 7.6733208d-2, 3.5441551d-1,
     : -1.8635208d-1, -5.8089087d-1, -2.7179217d-1,
     : -5.6466168d-2, 0.d0, 1.8494685d-1, 2.2351345d-1,
     : -2.8301294d-1, 1.1186650d-1, 2.2022085d-1, -5.5600556d-1,
     : -3.4565694d-1, -7.3847102d-2/
* coefficients for d2 rotation matrices
* stored (by column) for each of 8 angles and for l=2,4,6,8
      data d2/6.1237244d-1, 6.0772078d-1, 5.7135126d-1, 4.5927932d-1,
     : 3.0618621d-1, 1.5309311d-1, 4.1021174d-2, 4.6516566d-3,
     : 1.8652037d-4,-3.9528471d-1, -3.7142331d-1, -1.9586858d-1,
     : 2.2234766d-1, 4.9410589d-1, 4.1999000d-1, 1.4645800d-1,
     : 1.7856130d-2, 7.2213358d-4,3.2021721d-1, 2.7493911d-1,
     : -1.7236033d-2, -3.4523418d-1, 4.0027167d-2, 4.8532921d-1,
     : 2.7741249d-1, 3.8036292d-2, 1.5591157d-3,-2.7731624d-1,
     : -2.0847588d-1, 1.5831807d-1, 1.1374297d-1, -3.2931303d-1,
     : 2.5240112d-1, 3.9848096d-1, 6.4627284d-2, 2.6984111d-3/
* hyperbolic tangent scaling for vsum and vsigma
      data alph /1.2d0/
* quadrupole moments
      data quadb, quadh /-2.3439d0, 0.4615d0/
* f(theta) see Eq.(26) and following paragraph
      data fthet /
     : -4, 2, -1, -3.9088465, 1.9544233, -9.9240388d-1, -3.1961524,
     : 1.5980762, -9.3301270d-1,-1, .5, -.75, 2, -1, -.5, 5, -2.5,
     : -.25,7.1961524, -3.5980762, -6.6987298d-2, 7.9088465,
     : -3.9544233, -7.5961235d-3, 7.9963450, -3.9981725, -3.0458649d-4/
*  coefficients for vsig (alpha = 1.2)
*  by rows: lam1 lam2 r0 c1 c2 c3 cl
      data vzzl1 /
     : 7.3665552d-1, 7.7842593d-1, 7.1656045d-1, 7.1249612d-1,
     : 6.8658552d-1, 6.8374505d-1, 7.1198324d-1, 7.1928441d-1,
     : 7.1774139d-1/
      data vzzl2 /
     : 1.6030247, 1.5866667, 1.6334822, 1.6329619, 1.6560966,
     : 1.6635153, 1.6148457, 1.5929321, 1.6056816/
      data vzzr0 /
     : 7.0697359, 7.0929012, 7.0364426, 5.5981557, 7.0967916,
     : 7.0801266, 7.0187393, 7.0587134, 7.0910992/
      data vzzc1 /
     : -9.1086921d+3, -1.3740979d+4, -7.6380681d+3, -1.4341256d+4,
     : -1.3351106d+4, -1.7308561d+4, -2.8583659d+4, -3.2744364d+4,
     : -3.1383698d+4/
      data vzzc2 /
     : -1.2031329d+6, -9.5166037d+5, -2.0387131d+6, -2.8999479d+6,
     : -3.8692754d+6, -4.6745283d+6, -3.6549076d+6, -3.1204702d+6,
     : -3.4290337d+6/
      data vzzc3 /
     : 1.4616721d+6, 1.3441439d+6, 1.7978931d+6, 2.0285758d+6,
     : 2.3289607d+6, 2.5466954d+6, 2.0824783d+6, 1.8635943d+6,
     : 1.9854653d+6/
      data vzzcl /
     : 4.9130256d+6, 4.7807834d+6, 4.3005085d+6, -1.1848526d+6,
     : -3.1407132d+6, -7.1233496d+6, -1.0703781d+7, -1.2013084d+7,
     : -1.1735345d+7/
*  coefficients for vsum (alpha = 1.2)
*  by rows: lam1 lam2 r0 c1 c2 c3 cl
      data vsl1 /
     : 6.9950648d-1, 7.1750000d-1, 7.0486111d-1, 6.9893504d-1,
     : 7.1361111d-1, 7.2972029d-1, 7.3254302d-1, 7.3500000d-1,
     : 7.2406250d-1/
      data vsl2 /
     : 1.9670949, 1.9350000, 1.9250000, 1.9462084, 1.9325000,
     : 1.8939236, 1.8720154, 1.8675000, 1.8900000/
      data vsr0 /
     : 7.0639772, 7.0000000, 7.0194444, 7.0178014, 7.0777778,
     : 7.0355131, 7.0683865, 7.0729167, 7.0364583/
      data vsc1 /
     : -1.1252720d+4, -1.3363313d+4, -1.2089250d+4, -1.0415273d+4,
     : -9.9565634d+3, -9.6780312d+3, -8.7507134d+3, -8.4568622d+3,
     : -7.4092461d+3/
      data vsc2 /
     : -7.8402201d+6, -6.5298777d+6, -6.0711797d+6, -7.2349938d+6,
     : -7.1548828d+6, -5.9284757d+6, -5.1742241d+6, -5.0252057d+6,
     : -5.8061379d+6/
      data vsc3 /
     : 3.1777619d+6, 2.7140255d+6, 2.5664188d+6, 3.0038986d+6,
     : 3.0329425d+6, 2.6620667d+6, 2.4576221d+6, 2.4260561d+6,
     : 2.7181191d+6/
      data vscl /
     : -3.6393374d+6, -3.3310511d+6, -3.2599576d+6, -1.7632042d+6,
     : 3.3488770d+5, 2.1008190d+6, 3.3128397d+6, 3.7280615d+6,
     : 3.8501318d+6/
*  coefficients for vdif; by rows: lam1 lam2 r0 c1 c2 c3 cl
      data vdl1 /
     : 6.9697773d-1, 7.1133873d-1, 7.1500000d-1, 7.0687500d-1,
     : 7.4375000d-1, 8.8990741d-1, 4.6316793d-1, 4.8169364d-1,
     : 5.5416667d-1/
      data vdl2 /
     : 1.8284255, 1.8348611, 1.8000000, 1.8225000, 1.8675000,
     : 1.8483333, 1.7841860, 1.4386183, 1.4366667/
      data vdr0 /
     : 7.0498797, 7.0361883, 7.0000000, 7.0000000, 7.0000000,
     : 7.3175926, 7.7922568, 7.8082025, 4.6969136/
      data vdc1 /
     : 1.2130992d+3, 1.6624103d+3, 1.6431447d+3, 1.3002096d+3,
     : 1.1760673d+3, 1.7678183d+3, 3.9544378d+1, 1.6051198d+1,
     : 7.6410882/
      data vdc2 /
     : 2.9803558d+5, 4.2851713d+5, 4.0295846d+5, 3.0858331d+5,
     : 2.0027424d+5, 1.8922989d+5, 4.6341570d+3, 5.9327447d+3,
     : 6.5296232d+3/
      data vdc3 /
     : 6.5975773d+3, -6.7231701d+3, -1.7715736d+4, 1.2330285d+3,
     : 1.6416501d+4, -1.5794121d+4, 1.1472723d+4, -4.7066084d+1,
     : -4.2154652d+2/
      data vdcl /
     : 1.3833015d+6, 1.3854401d+6, 1.2659253d+6, 9.8479756d+5,
     : 5.4456294d+5, 1.6657798d+5, 2.6926208d+5, 6.9420392d+4,
     : -2.6340917d+4/
*  coefficients for vxz (includes sqrt(2)); by rows c1 c2 c3 c4
      data vxzc1 /0d0, 5.4753011d+2, 5.4233733d+2, 5.2045768d+2,
     : 4.7823895d+2, 4.4054637d+2, 4.2649081d+2, 5.9573375d+2,
     : 6.0547940d+2/
      data vxzc2 /0d0, -3.7290239d+2, -3.7017322d+2, -3.5845133d+2,
     : -3.3362999d+2, -3.2173985d+2, -3.0712151d+2,-4.0766612d+2,
     : -4.1310290d+2/
      data vxzc3 /0d0, 1.0724078d+2, 1.0679993d+2, 1.0483246d+2,
     : 1.0034229d+2, 9.9434116d+1, 9.5929641d+1, 1.1514076d+2,
     : 1.1610460d+2/
      data vxzc4 /0d0,-7.0112523d+0, -5.9222072d+0, -5.2392203d+0,
     : -4.8147973d+0, -4.9211872d+0, -5.1820009d+0,-7.3938814d+0,
     : -9.0503555d+0/

* determine potentials at angles
      rm3=one/r**3
      rm6=rm3*rm3
      qterm=
     :  0.75d0*quadb*quadh*half*(tanh(alph*(r-7.d0))+1)/r**5
      do 200 i=1,9
        vs(i)=vsc1(i)*exp(-vsl1(i)*r) +
     :        (vsc2(i)+vsc3(i)*r)*exp(-vsl2(i)*r)-
     :         half*(tanh(alph*(r-vsr0(i)))+1)*vscl(i)*rm6+
     :         qterm*fthet(2,i)

        vzz(i)=vzzc1(i)*exp(-vzzl1(i)*r) +
     :        (vzzc2(i)+vzzc3(i)*r)*exp(-vzzl2(i)*r)-
     :         half*(tanh(alph*(r-vzzr0(i)))+1)*vzzcl(i)*rm6+
     :         qterm*fthet(1,i)
        t1=vzzc1(i)*exp(-vzzl1(i)*r)
        t2= (vzzc2(i)+vzzc3(i)*r)*exp(-vzzl2(i)*r)
        t3=- half*(tanh(alph*(r-vzzr0(i)))+1)*vzzcl(i)/r**6
*        print *, t1,t2,t3

        vd(i)=vdc1(i)*exp(-vdl1(i)*r) +
     :        (vdc2(i)+vdc3(i)*r)*exp(-vdl2(i)*r)-
     :         half*(tanh(alph*(r-vdr0(i)))+1)*vdcl(i)*rm6+
     :         qterm*fthet(3,i)
*        vd(i)=vdc1(i)*exp(-vdl1(i)*r) +
*     :        (vdc2(i)+vdc3(i)*r)*exp(-vdl2(i)*r)

        sum=((vxzc1(i)/r+vxzc2(i))/r+vxzc3(i))/r+vxzc4(i)
        vxz(i)=exp(sum)
        if (i .eq. 1) vxz(i)=0.d0
* determine adiabatic A' potentials
        h11=vs(i)-vd(i)
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
