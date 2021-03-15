*system:  B(2P)+H2
* reference:  M. H. Alexander, J. Chem. Phys. 99, 6014 (1993).

      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(18)
1      print *, ' r (bohr)'
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vvl
100   format(' vsig',/,5(1pe16.8),/,'  vsum',/,5e16.8,
     :   /,'  vdif',/,4e16.8,/,'  v1  ',/,4e16.8)
      goto 1
99    end
      include "common/syusr"
      include "common/bausr"
      include "common/ground"
      subroutine pot (vv0, r)
*  subroutine to calculate the r-dependent coefficients in the
*  b-h2 potentials of alexander
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*  vv0        (not used)
*  variable in common block /covvl/
*    vvl:     vector of length 18 to store r-dependence of each term
*             in potential expansion
*    vvl(1-5) expansion coefficients in dl0 (l=0,2,4,6,8) of vsigma
*    vvl(6-10) expansion coefficients in dl0 (l=0,2,4,6,8) of vsum
*    vvl(11-14) expansion coefficients in dl2 (l=2,4,6,8) of vdif
*    vvl(15-18) expansion coefficients in dl1 (l=2,4,6,8) of v1

* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  19-dec-1995
* ----------------------------------------------------------------------


      implicit double precision (a-h,o-z)
      dimension vsl1(9), vsl2(9), vsr0(9), vsc1(9), vsc2(9),
     :  vsc3(9), vscl(9)
      dimension vsgl1(9), vsgl2(9), vsgr0(9), vsgc1(9), vsgc2(9),
     :  vsgc3(9), vsgcl(9)
      dimension vdl1(9), vdl2(9), vdc1(9), vdc2(9), vdc3(9)
      dimension v1c1(9), v1c2(9), v1c3(9), v1c4(9)
*      dimension beta(9)
      dimension vs(9),vsg(9),vd(9),v1(9)
      dimension d0(45),d1(36),d2(36),aa(45)
      dimension xsg(5),xs(5),xd(4),x1(4),kpvt(9),qraux(9),
     :          work(24),rsd(9)
      common /covvl/ vvl(18)
* angles (not needed here, but included for completeness)
*      data beta / 90d0, 85d0, 75d0, 60d0, 45d0, 30d0, 15d0, 5d0, 1d0/
      data zero, one, half /0.d0, 1.d0, 0.5d0/
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
*  coefficients for vsum (alpha = 1.2)
*  by rows: lam1 lam2 r0 c1 c2 c3 cl
      data vsl1 /6.5814285d-1,
     : 6.5426580d-1, 6.5513501d-1, 6.7748590d-1,
     : 7.3209427d-1, 8.2916398d-1, 9.7223238d-1, 1.1134409d+0,
     : 1.1517466d+0/
      data vsl2 /2.0021967d0,
     : 2.0327296d+0, 2.0355360d+0, 2.0061479d+0,
     : 1.9432618d+0, 1.8602292d+0, 1.7464777d+0, 1.6555405d+0,
     : 1.6326108d+0/
      data vsr0 /7.045794d0,
     : 7.1324115d+0, 7.1627906d+0, 7.1650301d+0,
     : 7.1421146d+0, 7.2036650d+0, 7.1977878d+0, 6.9948418d+0,
     : 6.9076399d+0/
      data vsc1 /-8.6697587d+03,
     : -8.2423149d+3, -8.0287782d+3, -8.5705405d+3,
     : -1.1118417d+4, -1.9940384d+4, -5.9357430d+4,-2.2069911d+5,
     : -3.3115409d+5/
      data vsc2 /-8.864066d+06,
     : -1.0695870d+7, -1.1145481d+7, -1.0224078d+7,
     : -7.9010281d+6, -5.3974672d+6, -2.8857546d+6,-1.7603403d+6,
     : -1.5018912d+6/
      data vsc3 /3.5250684d+06,
     : 4.1404098d+6, 4.3003393d+6, 4.0245813d+6,
     : 3.2866776d+6, 2.4659017d+6, 1.6280576d+6, 1.3395732d+6,
     : 1.3243259d+6/
      data vscl /-1.3118926d+06,
     : -1.2121730d+6, -1.2022028d+6, -1.0245140d+6,
     : -7.9324370d+5, -6.8554155d+5, -7.5635119d+5,-5.8633283d+5,
     : -5.3313735d+5/
*  coefficients for vsig (alpha = 1.2)
*  by rows: lam1 lam2 r0 c1 c2 c3 cl
      data vsgl1 /1.0601837d0,
     : 1.0273815d+0, 9.0710375d-1, 6.9512329d-1,
     : 6.2438975d-1, 6.0405024d-1, 5.7869165d-1, 5.8559773d-1,
     : 5.8843156d-1/
      data vsgl2 /1.5154820d0,
     : 1.5282215d+0, 1.5714637d+0, 1.6493866d+0,
     : 1.6893976d+0, 1.7081953d+0, 1.7392445d+0, 1.7242208d+0,
     : 1.7129019d+0/
      data vsgr0 /1.0967849d+1,
     : 9.3518916d+0, 8.2459349d+0, 7.9190690d+0,
     : 7.6684286d+0, 7.4769063d+0, 8.0102338d+0, 7.8263202d+0,
     : 7.7633638d+0/
      data vsgc1 / -2.0045503d+05,
     : -1.4018746d+5, -4.5967612d+4, -9.5944255d+3,
     : -7.7909769d+3, -8.9254365d+3, -8.5816805d+3,-9.8265940d+3,
     : -1.0281942d+4/
      data vsgc2 /-2.6310975d+05,
     : -4.2196427d+5, -1.0739839d+6, -3.0085251d+6,
     : -4.9157972d+6, -6.1787693d+6, -8.0834115d+6,-7.3528042d+6,
     : -6.7575738d+6/
      data vsgc3 /1.1432555d+06,
     : 1.1719843d+6, 1.3792184d+6, 2.1024291d+6,
     : 2.7381172d+6, 3.1165946d+6, 3.7303666d+6, 3.4472905d+6,
     : 3.2322326d+6/
      data vsgcl /-2.3713045d+06,
     : -2.6691266d+5, -8.2778277d+4, -5.1635381d+5,
     : -7.1211267d+5, -7.7200390d+5, -1.2257868d+6,-1.5242939d+6,
     : -1.7811467d+6/
*  coefficients for vdif; by rows: lam1 lam2 c1 c2 c3
      data vdl1 /6.8712288d-01,
     : 6.5505249d-1, 6.4624100d-1, 6.2648198d-1,
     : 6.0569878d-1, 6.2705315d-1, 7.1790842d-1, 1.5172412d+0,
     : 9.0468750d-1/
      data vdl2 /1.8188773d0,
     : 1.7643433d+0, 1.7420654d+0, 1.7239882d+0,
     : 1.6832859d+0, 1.7506317d+0, 2.4019188d+0, 4.8845032d-1,
     : 1.5781250d+0/
      data vdc1 /1.3745558d+03,
     : 1.0629404d+3, 9.3273852d+2, 7.0334124d+2,
     : 4.4484644d+2, 3.1528879d+2, 2.6490452d+2, 7.1940640d+3,
     : 2.4906707d+2/
      data vdc2 /3.8406628d5,
     : 3.1439263d+5, 2.8435691d+5, 2.2188838d+5,
     : 1.4663688d+5, 1.1530693d+5, -8.5550494d+5, 4.1531750d+1,
     : 5.8946732d+3/
      data vdc3 /1.9311852d3, 2.4136505d+3, 8.3459388d+2, 2.0206601d+3,
     : -5.8086128d+2, -7.4894720d+2, 3.3064330d+5,-3.5087948d+0,
     : 9.3020299d+1/
*  coefficients for v1 (includes sqrt(2)); by rows c1 c2 c3 c4
      data v1c1 /0d0, 5.4753011d+2, 5.4233733d+2, 5.2045768d+2,
     : 4.7823895d+2, 4.4054637d+2, 4.2649081d+2, 5.9573375d+2,
     : 6.0547940d+2/
      data v1c2 /0d0, -3.7290239d+2, -3.7017322d+2, -3.5845133d+2,
     : -3.3362999d+2, -3.2173985d+2, -3.0712151d+2,-4.0766612d+2,
     : -4.1310290d+2/
      data v1c3 /0d0, 1.0724078d+2, 1.0679993d+2, 1.0483246d+2,
     : 1.0034229d+2, 9.9434116d+1, 9.5929641d+1, 1.1514076d+2,
     : 1.1610460d+2/
      data v1c4 /0d0,-7.0112523d+0, -5.9222072d+0, -5.2392203d+0,
     : -4.8147973d+0, -4.9211872d+0, -5.1820009d+0,-7.3938814d+0,
     : -9.0503555d+0/

* determine potentials at angles
      rm3=one/r**3
      rm6=rm3*rm3
      do 200 i=1,9

        vs(i)=vsc1(i)*exp(-vsl1(i)*r) +
     :        (vsc2(i)+vsc3(i)*r)*exp(-vsl2(i)*r)-
     :         half*(tanh(alph*(r-vsr0(i)))+1)*vscl(i)/r**6

        vsg(i)=vsgc1(i)*exp(-vsgl1(i)*r) +
     :        (vsgc2(i)+vsgc3(i)*r)*exp(-vsgl2(i)*r)-
     :         half*(tanh(alph*(r-vsgr0(i)))+1)*vsgcl(i)*r**(-6)
        t1=vsgc1(i)*exp(-vsgl1(i)*r)
        t2= (vsgc2(i)+vsgc3(i)*r)*exp(-vsgl2(i)*r)
        t3=- half*(tanh(alph*(r-vsgr0(i)))+1)*vsgcl(i)/r**6
*        print *, t1,t2,t3

        vd(i)=vdc1(i)*exp(-vdl1(i)*r) +
     :        (vdc2(i)+vdc3(i)*r)*exp(-vdl2(i)*r)

        sum=((v1c1(i)/r+v1c2(i))/r+v1c3(i))/r+v1c4(i)
100     continue
        v1(i)=exp(sum)
        if (i .eq. 1) v1(i)=0.d0
200   continue
      sum=0
* determine matrix of d's at angles for least-squares fit
*      call d2matev(9,beta,d0,d1,d2)
* (n.b. these are already in common, subroutine is included
*  only for completeness)
* solve simultaneous equations for solutions
* first for vsigma
      tol=1.e-10
      call dcopy(45,d0,1,aa,1)
      call dqrank(aa,9,9,5,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,9,9,5,kr,vsg,xsg,rsd,kpvt,qraux)
      call dcopy(5,xsg,1,vvl,1)
      call dqrlss(aa,9,9,5,kr,vs,xs,rsd,kpvt,qraux)
      call dcopy(5,xs,1,vvl(6),1)
      call dcopy(36,d2,1,aa,1)
      call dqrank(aa,9,9,4,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,9,9,4,kr,vd,xd,rsd,kpvt,qraux)
      call dcopy(4,xd,1,vvl(11),1)
      call dcopy(36,d1,1,aa,1)
      call dqrank(aa,9,9,4,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,9,9,4,kr,v1,x1,rsd,kpvt,qraux)
      call dcopy(4,x1,1,vvl(15),1)
* convert to atomic units
      econv=1.d0/219474.6d0
      call dscal(18,econv,vvl,1)
      end
      subroutine d2matev(nang, theta,vec0,vec1,vec2)
* ---------------------------------------------------
* returns as a vector the first 5 even reduced rotation matrix elements
*  l
* d  (theta)  ;  l=0:2:8
*  00 and the first 4 even reduced
* rotation matrix elements
*  l
* d  (theta)  ;  l=2:2:8
*  10 and the first 4 even reduced rotation matrix elements
*  l
* d (theta)  ;  l=2:2:8
*  20

* variables in call list
*  nang:  number of angles
*  theta: vector of angles (in degrees)
*  vec1:  on return: matrix of d0's, nang rows, 5 columns
*  vec1:  on return: matrix of d1's, nang rows, 4 columns
*  vec2:  on return: matrix of d3's, nang rows, 4 columns
* ---------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension angle(1)
      dimension vec0(8,5),vec1(8,4),vec2(8,4)
      data one, pi /1.d0, 3.141592653589793d0/
      do i=1,nang
        ang=theta(i)*pi/180d0
        x=cos(ang)
        sn=sin(ang)

        x2=x*x
        sn2=one-x2
        x4=x2*x2
        x6=x4*x2
        sn2=1-x2
        vec2(i,1)=sqrt(1.5d0)*sn2/2.d0
        vec2(i,2)=-sqrt(10.d0)*sn2*(one-7.d0*x2)/8.d0
        vec2(i,3)=sqrt(105.d0)*sn2*(one-18.d0*x2+33.d0*x4)/32.d0
        vec2(i,4)=-3d0*sqrt(35.d0)*sn2*(one-33.d0*x2+143.d0*x4-
     :            143.d0*x6)/64.d0
        vec1(i,1)=-sqrt(1.5d0)*sn*x
        vec1(i,2)=sqrt(5.d0)*sn*x*(3d0-7d0*x2)/4d0
        vec1(i,3)=-sqrt(10.5d0)*sn*x*(5d0-30d0*x2+33d0*x4)/8d0
        vec1(i,4)=-3d0*sqrt(0.5d0)*x*sn*(-35d0+385d0*x2
     :            -1001d0*x4+715d0*x6)/32d0
        vec0(i,1)=one
        vec0(i,2)=(-one+3d0*x2)/2d0
        vec0(i,3)=(3d0-30d0*x2+35d0*x4)/8d0
        vec0(i,4)=(-5d0+105d0*x2-315d0*x4+231d0*x6)/16d0
        vec0(i,5)=(35d0-1260d0*x2+6930d0*x4-12012d0*x6
     :            +6435d0*x8)/128d0

      enddo
      return
      end
