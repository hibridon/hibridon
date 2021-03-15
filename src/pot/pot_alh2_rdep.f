*system:  Al(2P)+H2, Dubernet-Hutson expansion of Williams-Alexander r-dependent
*  PES's references: J. Williams, M. H. Alexander, J. Chem. Phys. 112, 5722 (2000)
*   M.-L. Dubernet and J. M. Hutson, J. Chem. Phys. 101, 1939 (1994).
*   M. H. Alexander, J. Chem. Phys. 99, 6014 (1993).
*   M. H. Alexander and M. Yang, J. Chem. Phys. 103, 7956 (1995).
*   X. Tan, P. J. Dagdigian, J. Williams, and Millard H. Alexander,
*       J. Chem. Phys. 114, 8938 (2001)


      include "common/syusr"
      include "common/ground"
      include "common/bausr"
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      common /coselb/ ibasty
      potnam='WILLIAMS-ALEXANDER CBS-<v=0> AL(2P)H2 DUBERNET-HUTSON'
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
      potnam='WILLIAMS-ALEXANDER CBS-<v=0> AL(2P)H2 DUBERNET-HUTSON'
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

99    return
      end
      subroutine pot (vv0, r)
* ----------------------------------------------------------------------
*  subroutine to calculate the r-dependent coefficients in the
*  al-h2 potentials of williams and alexander using body-frame
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
      data alph /0.6d0/
*  expansion coefficients for vxx
*  by rows: lam1 lam2 r0 c1 c2 c3 cl
      data vxxl1/
     : 4.7635235d-1,4.8405976d-1,5.2041168d-1,9.8360830d-1,8.9640448d-1/
      data vxxl2/
     : 1.4887142d0,1.5059201d0,1.4055146d0,7.6621430d-1,7.3186544d-1/
      data vxxr0/
     : 9.1988983d0,8.8408721d0,8.6744696d0,8.3053343d0,8.0337789d0/
      data vxxc1/
     : -7.8902050d3,-6.8375056d3,-6.9972586d3,1.1607942d6,1.3172020d6/
      data vxxc2/
     : -7.5035838d6,-8.7871061d6,-3.2599856d6,-2.3473725d5,-5.4961680d5/
      data vxxc3/
     : 2.8471078d6,3.2477167d6,1.6325342d6,8.3316365d3,2.8577679d4/
      data vxxcl/
     : -1.0636571d7,-3.8856473d6,-4.2876403d6,-3.0055582d6,-7.0162039d6/
*  expansion coefficients for vyy
*  by rows: lam1 lam2 r0 c1 c2 c3 cl
      data vyyl1/
     : 8.4585977d-1,1.3430527d0,3.9060419d-1,1.4894368d0,1.4249515d0/
      data vyyl2/
     : 1.5868660d0,7.7965334d-1,1.0006896d0,5.9984745d-1,5.8106016d-1/
      data vyyr0/
     : 8.2204911d0,6.4727470d0,7.9537417d0,5.8914182d0,5.7517051d0/
       data vyyc1/
     : -5.4412868d4,1.3329020d6,-3.9900182d2,6.8250412d5,4.0108022d5/
      data vyyc2/
     : -9.4448208d5,-7.3355322d4,4.4370176d5,-7.4820933d3,-7.3807054d3/
      data vyyc3/
     : 1.1314517d6,5.5254227d3,-7.8218285d4,-3.0036811d1,2.9364975d0/
      data vyycl/
     : -4.1379529d6,4.6535974d6,-3.9029090d6,8.2277480d6,8.5409222d6/
*  expansion coefficients for vzz
*  by rows: lam1 lam2 r0 c1 c2 c3 cl
      data vzzl1/
     : 4.3681432d-1,7.7130849d-1,6.6320289d-1,3.3072306d-1,3.1975501d-1/
      data vzzl2/
     : 1.0868973d0,1.6729581d0,1.6280363d0,9.9274825d-1,9.7565900d-1/
      data vzzr0/
     : 8.1166638d0,8.1632089d0,8.1775262d0,7.5862856d0,7.4986816d0/
      data vzzc1/
     : -9.2136607d1,-2.8959449d4,-1.4793822d4,-2.0963876d2,-2.1174175d2/
      data vzzc2/
     : 9.1742747d5,-3.0690691d6,-6.1000126d5,4.8640967d5,4.3654540d5/
      data vzzc3/
     : -1.4209375d5,1.8549317d6,9.2127231d5,-8.3102149d4,-7.6054643d4/
      data vzzcl/
     : -4.6696971d6,-3.9831716d6,-5.0583991d6,-5.3296346d6,-5.9665222d6/
*  expansion coefficients for vxz
*  by rows: lam1 lam2 c1 c2
      data vxzlam1/
     : 4.3627299d-1, 4.0344178d-1, 3.9717913d-1/
      data vxzlam2/
     : 1.0278614d0, 9.8513398d-1, 9.7034323d-1/
      data vxzc1/
     : 8.4049633d2, 7.7159446d2, 4.7195596d2/
      data vxzc2/
     : 4.3351717d4, 5.1675391d4, 3.3274947d4/
* determine potentials at angles
      rm1=one/r
      rm2=rm1*rm1
      rm3=one/r**3
      rm5=rm3*rm1*rm1
      rm6=rm5*rm1
* quadal and quadh2 are quadrapole moments of Al and H2
* averaged over v=0 vibrational wavefunction
* cm and rad are conversion factors from au to cm-1
* and from degrees to radians respectively
      cm=2.194746d05
      pi=acos(-1.0d0)
      rad=pi/180.0d0
      quadal=-5.159d0*cm
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
     :        +0.75d0*quadh2*quadal*qxx*half*(tanh(alphq*
     :        (r-12.0D0))+one)*rm5 
        vyy(i)=vyyc1(i)*dexp(-vyyl1(i)*r)+
     :        (vyyc2(i)+vyyc3(i)*r)*dexp(-vyyl2(i)*r)-
     :        half*(tanh(alph*(r-vyyr0(i)))+one)*vyycl(i)*rm6
     :        +0.75d0*quadh2*quadal*qyy*half*(tanh(alphq*
     :        (r-12.0D0))+one)*rm5
        vzz(i)=vzzc1(i)*dexp(-vzzl1(i)*r)+
     :        (vzzc2(i)+vzzc3(i)*r)*dexp(-vzzl2(i)*r)-
     :        half*(tanh(alph*(r-vzzr0(i)))+one)*vzzcl(i)*rm6
     :        +0.75d0*quadh2*quadal*qzz*half*(tanh(alphq*
     :        (r-12.0D0))+one)*rm5 
        if (i.eq.1 .or. i.eq.5) then
* vxz vanishes at 0 and 90 degrees
           vxz(i)=0.0d0
        else
* NB exponentials are negative, corrected 30 march 2002
           vxz(i)=vxzc1(i-1)*dexp(-vxzlam1(i-1)*r)
     :           +vxzc2(i-1)*dexp(-vxzlam2(i-1)*r)
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
      return
      end
