* System:  PH(A 3Pi)+He, original ab initio CEPA PES's
* Reference: 
* Ch. Kolczewski, K. Fink, V. Staemmler, and L. Neitsch 
* J. Chem. Phys. 106, 7637 (1997)


      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(11)
      common /coconv/ econv
      include "common/parpot"
      potnam='FINK-STAEMMLER HePH(A) CEPA'
      econv=219474.6d0
1      print *, ' r (bohr)'
      rshift=0.5
      xfact=0.8
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vsum',/,7(1pe16.8),/,
     :    '  vdif',/,5e16.8)
      goto 1
99    end
      include "common/syusr"
      include "common/bausr"
      include "common/ground"
* --------------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      potnam='FINK-STAEMMLER HePH(A) CEPA'
      lammin(1)=1
      lammax(1)=6
      lammin(2)=2
      lammax(2)=6
      mproj(1)=0
      mproj(2)=2
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
*  -----------------------------------------------------------------------
      subroutine pot (vv0, r)
*  subroutine to calculate the r-dependent coefficients in the
*  He-PH(A) potentials of fink and staemmler
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*  vv0      spherically symmetric term in V0
*  variable in common block /covvl/
*    vvl:     vector of length 11 to store r-dependence of each term
*             in potential expansion
*    vvl(1-6) expansion coefficients in dl0 (l=1:6) of vsum
*    vvl(7-11) expansion coefficients in dl2 (l=2:6) of vdif

* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  20-mar-1996
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension xlam1(22),xlam2(22),c1(22),c2(22),
     :          clr(22),vsum(11),xsum(11),vdif(11),xdif(11),
     :          ddif(11),vap(11),va2p(11),
     :          d0(77),d2(45),aa(121)
      dimension kpvt(11),qraux(11),work(55),rsd(11)

      common /covvl/ vvl(11)
      common /coconv/ econv

      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.5d0/
* for distances beyond rmax difference potential is damped
      data rmax /11d0/
* coefficicients for d0 rotation matrices
* stored (by column) for each of 11 angles and for l=0:6
* angles are: 0 22.50 45 67.50  90 112.50 135 146.25 157.50 168.75 180
      data d0/
     : 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0,
     : 1.d0, 9.2387953d-1, 7.0710679d-1, 3.8268344d-1, 1.3349125d-8,
     : -3.8268342d-1, -7.0710677d-1, -8.3146960d-1, -9.2387952d-1,
     : -9.8078528d-1, -1.d0, 1.d0, 7.8033009d-1, 2.5000001d-1,
     : -2.8033008d-1, -5.0000000d-1, -2.8033010d-1, 2.4999997d-1,
     : 5.3701254d-1, 7.8033006d-1, 9.4290964d-1, 1.d0, 1.d0,
     : 5.8563198d-1, -1.7677668d-1, -4.3391842d-1, -2.0023687d-8,
     : 4.3391841d-1, 1.7677673d-1, -1.8986961d-1, -5.8563193d-1,
     : -8.8746296d-1, -1.d0, 1.d0, 3.6159588d-1, -4.0625000d-1,
     : -8.0345887d-2, 3.7500000d-1, -8.0345840d-2, -4.0625001d-1,
     : -1.2648549d-1, 3.6159581d-1, 8.1603633d-1, 1.d0, 1.d0,
     : 1.3282228d-1, -3.7565048d-1, 2.9179007d-1, 2.5029609d-8,
     : -2.9179009d-1, 3.7565046d-1, 3.4119961d-1, -1.3282219d-1,
     : -7.3067118d-1, -1.d0, 1.d0, -7.6358299d-2, -1.4843752d-1,
     : 2.7167082d-1, -3.1250000d-1, 2.7167079d-1, -1.4843745d-1,
     : -4.1470677d-1, -7.6358384d-2, 6.3379421d-1, 1.d0/
* coefficicients for d2 rotation matrices
* stored (by column) for each of 9 angles and for l=2:6
      data d2/
     : 8.9679867d-2, 3.0618622d-1, 5.2269257d-1, 6.1237244d-1,
     : 5.2269257d-1, 3.0618622d-1, 1.8901383d-1, 8.9679867d-2,
     : 2.3307038d-2, 1.8526582d-1, 4.8412292d-1, 4.4727126d-1,
     : -0.0000000d+0, -4.4727126d-1, -4.8412292d-1, -3.5141877d-1,
     : -1.8526582d-1, -5.1114725d-2, 2.8798601d-1, 4.9410588d-1,
     : 8.4775167d-3, -3.9528471d-1, 8.4775167d-3, 4.9410588d-1,
     : 4.6843615d-1, 2.8798601d-1, 8.6259556d-2, 3.8249228d-1,
     : 3.2021721d-1, -3.3173380d-1, 0.0000000d+0, 3.3173380d-1,
     : -3.2021721d-1, -4.9929639d-1, -3.8249228d-1, -1.2751613d-1,
     : 4.5386126d-1, 4.0027151d-2, -2.5372550d-1, 3.2021721d-1,
     : -2.5372550d-1, 4.0027151d-2, 4.2780021d-1, 4.5386126d-1,
     : 1.7331796d-1/

* coefficients for expansion of vap(1st ll entries) and
* for va2p (entries 11:22)
      data xlam1/
     : 1.3723853d+0, 1.3511692d+0, 1.3318451d+0, 1.3129986d+0,
     : 1.3068923d+0, 1.3313639d+0, 1.3992370d+0, 1.4436048d+0,
     : 1.4844806d+0, 1.5128079d+0, 1.5257139d+0, 1.3729004d+0,
     : 2.0317691d+0, 1.7829720d+0, 1.7424739d+0, 1.7099923d+0,
     : 1.7145991d+0, 1.9084418d+0, 1.4189241d+0, 1.4760713d+0,
     : 1.5174301d+0, 1.5263291d+0/
      data xlam2/
     : 4.7903887d+0, 4.8890249d+0, 4.9955783d+0, 4.9874717d+0,
     : 4.9213477d+0, 4.8037001d+0, 5.1263281d+0, 5.4059838d+0,
     : 5.6982993d+0, 5.8014434d+0, 5.9185011d+0, 4.7926524d+0,
     : 3.2422307d+0, 3.0625022d+0, 3.8016911d+0, 3.7061485d+0,
     : 3.5730424d+0, 4.5813305d+0, 5.3391111d+0, 5.6733444d+0,
     : 5.8866684d+0, 5.9191158d+0/
      data c1 /
     : 1.3697180d+6, 1.4173103d+6, 1.6869582d+6, 1.8841598d+6,
     : 1.7924871d+6, 1.9880898d+6, 3.9328876d+6, 6.8530491d+6,
     : 1.2738434d+7, 1.9714006d+7, 2.4620807d+7, 1.3719233d+6,
     : -1.0361593d+7, -6.1466047d+5, -1.4620216d+6, -4.1563815d+5,
     : 8.0278270d+4, -1.7155157d+7, 6.9339644d+6, 1.2595208d+7,
     : 2.0951383d+7, 2.4672470d+7/
      data c2/
     : -2.1952907d+5, -2.1060616d+5, -2.2454416d+5, -2.4491838d+5,
     : -2.3081743d+5, -2.6801104d+5, -5.3680377d+5, -8.9835589d+5,
     : -1.6485396d+6, -2.4881234d+6, -3.1611046d+6, -2.1983149d+5,
     : 5.9265121d+6, 2.3292751d+6, 2.8580795d+6, 2.4813716d+6,
     : 2.2783542d+6, 9.3244620d+6, -8.2651884d+5, -1.5548813d+6,
     : -2.6657398d+6, -3.1659782d+6/
      data clr /
     : 1.1904746d+6, 1.2671956d+6, 1.4071969d+6, 1.0619888d+6,
     : 1.1717403d+6, 1.1470695d+6, 2.1376500d+6, 3.5427587d+6,
     : 4.2278950d+6, 5.3101139d+6, 4.5460224d+6, 1.1936689d+6,
     : 2.8434915d+6, 2.9856441d+6, 3.2502732d+6, 3.4711867d+6,
     : 3.4478944d+6, 3.7022722d+6, 2.9042404d+6, 3.9855597d+6,
     : 4.3778701d+6, 4.5626510d+6/
* determine A' and A" potentials at angles
      rm3=one/r**3
      rm6=rm3*rm3
      do 100 i=1,11
        vap(i)=
     :        (c1(i)+c2(i)*r)*exp(-xlam1(i)*r)-
     :         half*(tanh(alph*(r-xlam2(i)))+1)*clr(i)*rm6
        j=i+11
        va2p(i)=
     :        (c1(j)+c2(j)*r)*exp(-xlam1(j)*r)-
     :         half*(tanh(alph*(r-xlam2(j)))+1)*clr(j)*rm6
        vsum(i)=half*(va2p(i)+vap(i))
* don't compute vdif for colinear geometries
        if (i.ne.1 .and. i.ne.11) then
          vdif(i-1)=half*(va2p(i)-vap(i))
* for long range damp out difference potential
          if (r .gt. rmax) then
            damp=-half*(tanh(3.d0*(r-rmax))-one)
            vdif(i-1)=damp*vdif(i-1)
          endif
        endif
100   continue
* solve simultaneous equations for solutions
      lsum=7
      ldif=5
*      lsum=9
*      ldif=7
      tol=1.e-10
      call dscal(11,zero,vvl,1)
      call dcopy(77,d0,1,aa,1)
      call dqrank(aa,11,11,lsum,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,11,11,lsum,kr,vsum,xsum,rsd,kpvt,qraux)
* remove terms less than 0.2 cm-1 in size
      do 110 i=1,11
        if (abs(xsum(i)) .lt. 0.2d0) xsum(i)=zero
110   continue
      vv0=xsum(1)/econv
      call dcopy(lsum-1,xsum(2),1,vvl,1)
      call dcopy(45,d2,1,aa,1)
      call dqrank(aa,9,9,ldif,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,9,9,ldif,kr,vdif,xdif,rsd,kpvt,qraux)
      do 120 i=1,9
        if (abs(xdif(i)) .lt. 0.2d0) xdif(i)=zero
120   continue
      call dcopy(ldif,xdif,1,vvl(7),1)
* convert potential to hartree
      call dscal(11,one/econv,vvl,1)
      end
