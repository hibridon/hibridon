*  system:  H2O-He, SAPT PES from Patkowski et al., J. Molec. Struct.
*  591, 231-243 (2002)
*  subr to compute V vs. (R, theta, phi) provided by B.H.Yang and P.
*  Stancil (U. GA), July-2009
*
*  H2O defined to lie in xz plane, with origin at center of mass
*  x axis is C2 symmetry axis of molecules
*  z axis is a inertial axis of molecule (perpendicular to C2 axis)
*  theta = 90, phi - 0 has He on O side of molecule
*  phi = 0 has all 4 atoms coplanar
*
*  Note:  subr heh2osapt requires 3 data files to be in hibxx/bin/progs/potdata:
*         h2o_coefd.dat, h2o_coefi.dat, h2o_params.dat
*
      include "common/syusr"
      include "common/ground"
      include "common/bausr"
* ------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      dimension xxl(15)
      common /covvl/ vvl(1)
      common /coloapot/ s4pi
      include "common/parpot"
      s4pi = sqrt ( 4.d0 * acos(-1.d0) )
      potnam='Patkowski et al. H2O-He SAPT PES'
      print *, potnam
1     print *, 'R (bohr):'
      read (5, *, end=99) r
      call pot(vv0,r)
*  convert from atomic units
      econv=1.d0/219474.6d0
*  isotropic term multiplying Y00
      xx0 = vv0 / econv * s4pi
      do i=1,15
        xxl(i) = vvl(i) / econv
      end do
      write (6, 100) xx0, (xxl(i), i=1,11)
100   format(' v(lam,0):',4(1pe16.8)/' v(lam,1):',3(1pe16.8)/
     :    ' v(lam,2):'3(1pe16.8)/' v(lam,3):',2(1pe16.8))
      goto 1
99    end
* ------------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
      implicit double precision (a-h,o-z)
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      common /conlam/ nlam, nlammx, lamnum(2)
      common /coloapot/ s4pi
      common /cosysi/ nscode, isicod, nterm
      potnam='Patkowski et al. H2O-He SAPT PES'
*
*  s4pi is factor to normalize isotropic term
      s4pi = sqrt ( 4.d0 * acos(-1.d0) )
      nterm = 4
      mproj(1) = 0
      mproj(2) = 1
      mproj(3) = 2
      mproj(4) = 3
      lammin(1) = 2
      lammin(2) = 1
      lammin(3) = 2
      lammin(4) = 3
      lammax(1) = 6
      lammax(2) = 5
      lammax(3) = 6
      lammax(4) = 5
*
      ipotsy = 2
*
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
*  calculate total number of anisotropic terms
      nlam = 0
      do i=1,nterm
        lmin = lammin(i)
        lmax = lammax(i)
        do lb = lmin, lmax, ipotsy
          nlam = nlam + 1
        end do
      end do
      nlammx = nlam
      return
      end
* ----------------------------------------------------------------------
      subroutine pot (vv0, r)
*  subroutine to calculate the r-dependent coefficients
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:          interparticle distance
*  on return:
*    vv0         contains isotropic term (Y00)
*  variable in common block /covvl/
*    vvl:        vector of length 6 to store r-dependence of each term
*                in potential expansion
*    vvl(1-3):   expansion coefficients in Yl0 (l=2:6:2) of v(lam,0)
*    vvl(4-6):   expansion coefficients in [Yl1 + Y(l,-1)] (l=1:6:2) of v(lam,1)
*    vvl(7-9):  expansion coefficients in [Yl2 + Y(l,-2)] (l=2:6:2) of v(lam,2)
*    vvl(10-11): expansion coefficients in [Yl4 + Y(l,-4)] (l=3:5:2) of v(lam,3)
*  variable in common block /coloapot/
*    s4pi:       normalization factor for isotropic potential
*
*  uses linear least squares routines from lapack
*
* author:  paul dagdigian and millard alexander
* latest revision date:  28-aug-2009
*
      implicit double precision (a-h,o-z)
      dimension iwork(1000)
      dimension ylm(60,16), thetab(60), phib(60), vcalc(60),
     :  aa(960), vfit(16), kpvt(16), qraux(16),swork(16),
     :  work(1812), rsd(16),y00(60),y11(60),y20(60),y22(60),
     :  y31(60),y33(60),y40(60),y42(60),y51(60),y53(60),
     :  y60(60),y62(60),y71(60),y73(60),y80(60),y82(60),
     :  vcalcx(60)
      dimension y11t(60), y22t(60), y31t(60), y33t(60), 
     :  y42t(60), y51t(60), y53t(60), y62t(60), 
     :  y71t(60), y73t(60), y82t(60)
*
      common /covvl/ vvl(1)
      common /coloapot/ s4pi
*
*  coefficients for Ylm terms - all (lambda,mu) combinations up to lambda .le. 8
*  and mu .le. 3, with lambda + mu = even (total of 15 anisotropic terms)
*  determined on a grid of 10 theta and 6 phi values:
*    th = [0:10:90];   %grid of 10 angles
*    ph = [0:36:180];  %grid of  6 angles
*  Note:  The ylm terms for mu .gt. 0 are defined as:
*           [Y(lambda,mu) + (-1)^mu * Y(lambda, -nu)]/2
*
      data y00 /
     : 2.8209479d-1, 2.8209479d-1, 2.8209479d-1, 2.8209479d-1,
     : 2.8209479d-1, 2.8209479d-1, 2.8209479d-1, 2.8209479d-1,
     : 2.8209479d-1, 2.8209479d-1, 2.8209479d-1, 2.8209479d-1,
     : 2.8209479d-1, 2.8209479d-1, 2.8209479d-1, 2.8209479d-1,
     : 2.8209479d-1, 2.8209479d-1, 2.8209479d-1, 2.8209479d-1, 
     : 2.8209479d-1, 2.8209479d-1, 2.8209479d-1, 2.8209479d-1,
     : 2.8209479d-1, 2.8209479d-1, 2.8209479d-1, 2.8209479d-1,
     : 2.8209479d-1, 2.8209479d-1, 2.8209479d-1, 2.8209479d-1,
     : 2.8209479d-1, 2.8209479d-1, 2.8209479d-1, 2.8209479d-1,
     : 2.8209479d-1, 2.8209479d-1, 2.8209479d-1, 2.8209479d-1, 
     : 2.8209479d-1, 2.8209479d-1, 2.8209479d-1, 2.8209479d-1,
     : 2.8209479d-1, 2.8209479d-1, 2.8209479d-1, 2.8209479d-1,
     : 2.8209479d-1, 2.8209479d-1, 2.8209479d-1, 2.8209479d-1,
     : 2.8209479d-1, 2.8209479d-1, 2.8209479d-1, 2.8209479d-1,
     : 2.8209479d-1, 2.8209479d-1, 2.8209479d-1, 2.8209479d-1/
      data y11 /
     : 0d0, -5.9994429d-2, -1.1816596d-1, -1.7274707d-1, -2.2207936d-1,
     : -2.6466387d-1, -2.9920671d-1, -3.2465830d-1, -3.4024532d-1,
     : -3.4549415d-1, 0d0, -4.8536513d-2, -9.5598269d-2, -1.3975532d-1,
     : -1.7966598d-1, -2.1411757d-1, -2.4206331d-1, -2.6265408d-1,
     : -2.7526424d-1, -2.7951064d-1, 0d0, -1.8539298d-2, -3.6515289d-2,
     : -5.3381782d-2, -6.8626296d-2, -8.1785635d-2, -9.2459958d-2,
     : -1.0032493d-1, -1.0514159d-1, -1.0676356d-1, 0d0, 1.8539298d-2,
     : 3.6515289d-2, 5.3381782d-2, 6.8626296d-2, 8.1785635d-2,
     : 9.2459958d-2, 1.0032493d-1, 1.0514159d-1, 1.0676356d-1, 0d0,
     : 4.8536513d-2, 9.5598269d-2, 1.3975532d-1, 1.7966598d-1,
     : 2.1411757d-1, 2.4206331d-1, 2.6265408d-1, 2.7526424d-1,
     : 2.7951064d-1, 0d0, 5.9994429d-2, 1.1816596d-1, 1.7274707d-1,
     : 2.2207936d-1, 2.6466387d-1, 2.9920671d-1, 3.2465830d-1,
     : 3.4024532d-1, 3.4549415d-1/
      data y20 /
     : 6.3078313d-1, 6.0225247d-1, 5.2010172d-1, 3.9423946d-1,
     : 2.3984654d-1, 7.5545027d-2, -7.8847891d-2, -2.0471015d-1,
     : -2.8686091d-1, -3.1539157d-1,  6.3078313d-1, 6.0225247d-1,
     : 5.2010172d-1, 3.9423946d-1, 2.3984654d-1, 7.5545027d-2,
     : -7.8847891d-2, -2.0471015d-1, -2.8686091d-1, -3.1539157d-1, 
     : 6.3078313d-1, 6.0225247d-1, 5.2010172d-1, 3.9423946d-1,
     : 2.3984654d-1, 7.5545027d-2, -7.8847891d-2, -2.0471015d-1,
     : -2.8686091d-1, -3.1539157d-1,  6.3078313d-1, 6.0225247d-1,
     : 5.2010172d-1, 3.9423946d-1, 2.3984654d-1, 7.5545027d-2,
     : -7.8847891d-2, -2.0471015d-1, -2.8686091d-1, -3.1539157d-1, 
     : 6.3078313d-1, 6.0225247d-1, 5.2010172d-1, 3.9423946d-1,
     : 2.3984654d-1, 7.5545027d-2, -7.8847891d-2, -2.0471015d-1,
     : -2.8686091d-1, -3.1539157d-1,  6.3078313d-1, 6.0225247d-1,
     : 5.2010172d-1, 3.9423946d-1, 2.3984654d-1, 7.5545027d-2,
     : -7.8847891d-2, -2.0471015d-1, -2.8686091d-1, -3.1539157d-1/
      data y22 /
     : 0d0, 1.1647592d-2, 4.5185498d-2, 9.6568551d-2, 1.5959920d-1,
     : 2.2667501d-1, 2.8970565d-1, 3.4108870d-1, 3.7462661d-1,
     : 3.8627420d-1, 0d0, 3.5993040d-3, 1.3963087d-2, 2.9841323d-2,
     : 4.9318864d-2, 7.0046429d-2, 8.9523970d-2, 1.0540221d-1,
     : 1.1576599d-1, 1.1936529d-1, 0d0, -9.4231002d-3, -3.6555836d-2,
     : -7.8125598d-2, -1.2911846d-1, -1.8338393d-1, -2.3437680d-1,
     : -2.7594656d-1, -3.0307929d-1, -3.1250239d-1, 0d0, -9.4231002d-3,
     : -3.6555836d-2, -7.8125598d-2, -1.2911846d-1, -1.8338393d-1,
     : -2.3437680d-1, -2.7594656d-1, -3.0307929d-1, -3.1250239d-1, 0d0,
     : 3.5993040d-3, 1.3963087d-2, 2.9841323d-2, 4.9318864d-2,
     : 7.0046429d-2, 8.9523970d-2, 1.0540221d-1, 1.1576599d-1,
     : 1.1936529d-1, 0d0, 1.1647592d-2, 4.5185498d-2, 9.6568551d-2,
     : 1.5959920d-1, 2.2667501d-1, 2.8970565d-1, 3.4108870d-1,
     : 3.7462661d-1, 3.8627420d-1/
      data y31 /
     : 0d0, -2.1601753d-1, -3.7748635d-1, -4.4437275d-1, -4.0178687d-1,
     : -2.6388021d-1, -6.9970562d-2, 1.2606511d-1, 2.7028522d-1,
     : 3.2318018d-1, 0d0, -1.7476185d-1, -3.0539287d-1, -3.5950511d-1,
     : -3.2505240d-1, -2.1348358d-1, -5.6607374d-2, 1.0198881d-1,
     : 2.1866534d-1, 2.6145826d-1, 0d0, -6.6753087d-2, -1.1664970d-1,
     : -1.3731873d-1, -1.2415897d-1, -8.1543470d-2, -2.1622093d-2,
     : 3.8956260d-2, 8.3522728d-2, 9.9868169d-2, 0d0, 6.6753087d-2,
     : 1.1664970d-1, 1.3731873d-1, 1.2415897d-1, 8.1543470d-2,
     : 2.1622093d-2, -3.8956260d-2, -8.3522728d-2, -9.9868169d-2, 0d0,
     : 1.7476185d-1, 3.0539287d-1, 3.5950511d-1, 3.2505240d-1,
     : 2.1348358d-1, 5.6607374d-2, -1.0198881d-1, -2.1866534d-1,
     : -2.6145826d-1, 0d0, 2.1601753d-1, 3.7748635d-1, 4.4437275d-1,
     : 4.0178687d-1, 2.6388021d-1, 6.9970562d-2, -1.2606511d-1,
     : -2.7028522d-1, -3.2318018d-1/
      data y33 /
     : 0d0, -2.1846395d-3, -1.6692606d-2, -5.2152978d-2, -1.1080812d-1,
     : -1.8755602d-1, -2.7099482d-1, -3.4619959d-1, -3.9849555d-1,
     : -4.1722382d-1, 0d0, 6.7509074d-4, 5.1582991d-3, 1.6116156d-2,
     : 3.4241592d-2, 5.7957996d-2, 8.3742006d-2, 1.0698156d-1,
     : 1.2314190d-1, 1.2892925d-1, 0d0, 1.7674105d-3, 1.3504602d-2,
     : 4.2192645d-2, 8.9645653d-2, 1.5173600d-1, 2.1923942d-1,
     : 2.8008135d-1, 3.2238967d-1, 3.3754116d-1, 0d0, -1.7674105d-3,
     : -1.3504602d-2, -4.2192645d-2, -8.9645653d-2, -1.5173600d-1,
     : -2.1923942d-1, -2.8008135d-1, -3.2238967d-1, -3.3754116d-1, 0d0,
     : -6.7509074d-4, -5.1582991d-3, -1.6116156d-2, -3.4241592d-2,
     : -5.7957996d-2, -8.3742006d-2, -1.0698156d-1, -1.2314190d-1,
     : -1.2892925d-1, 0d0, 2.1846395d-3, 1.6692606d-2, 5.2152978d-2,
     : 1.1080812d-1, 1.8755602d-1, 2.7099482d-1, 3.4619959d-1,
     : 3.9849555d-1, 4.1722382d-1/
      data y40 /
     : 8.4628438d-1, 7.2205787d-1, 4.0196624d-1, 1.9834790d-2,
     : -2.6996839d-1, -3.6181573d-1, -2.4462908d-1, -3.2159156d-3,
     : 2.2502838d-1, 3.1735664d-1,  8.4628438d-1, 7.2205787d-1,
     : 4.0196624d-1, 1.9834790d-2, -2.6996839d-1, -3.6181573d-1,
     : -2.4462908d-1, -3.2159156d-3, 2.2502838d-1, 3.1735664d-1, 
     : 8.4628438d-1, 7.2205787d-1, 4.0196624d-1, 1.9834790d-2,
     : -2.6996839d-1, -3.6181573d-1, -2.4462908d-1, -3.2159156d-3,
     : 2.2502838d-1, 3.1735664d-1,  8.4628438d-1, 7.2205787d-1,
     : 4.0196624d-1, 1.9834790d-2, -2.6996839d-1, -3.6181573d-1,
     : -2.4462908d-1, -3.2159156d-3, 2.2502838d-1, 3.1735664d-1, 
     : 8.4628438d-1, 7.2205787d-1, 4.0196624d-1, 1.9834790d-2,
     : -2.6996839d-1, -3.6181573d-1, -2.4462908d-1, -3.2159156d-3,
     : 2.2502838d-1, 3.1735664d-1,  8.4628438d-1, 7.2205787d-1,
     : 4.0196624d-1, 1.9834790d-2, -2.6996839d-1, -3.6181573d-1,
     : -2.4462908d-1, -3.2159156d-3, 2.2502838d-1, 3.1735664d-1/
      data y42 /
     : 0d0, 5.8393520d-2, 2.0274789d-1, 3.5543098d-1, 4.2954632d-1,
     : 3.7145697d-1, 1.8816934d-1, -5.3511807d-2, -2.5595553d-1,
     : -3.3452327d-1, 0d0, 1.8044590d-2, 6.2652543d-2, 1.0983421d-1,
     : 1.3273711d-1, 1.1478652d-1, 5.8147524d-2, -1.6536058d-2,
     : -7.9094609d-2, -1.0337338d-1, 0d0, -4.7241350d-2, -1.6402649d-1,
     : -2.8754970d-1, -3.4751028d-1, -3.0051500d-1, -1.5223219d-1,
     : 4.3291961d-2, 2.0707237d-1, 2.7063501d-1, 0d0, -4.7241350d-2,
     : -1.6402649d-1, -2.8754970d-1, -3.4751028d-1, -3.0051500d-1,
     : -1.5223219d-1, 4.3291961d-2, 2.0707237d-1, 2.7063501d-1, 0d0,
     : 1.8044590d-2, 6.2652543d-2, 1.0983421d-1, 1.3273711d-1,
     : 1.1478652d-1, 5.8147524d-2, -1.6536058d-2, -7.9094609d-2,
     : -1.0337338d-1, 0d0, 5.8393520d-2, 2.0274789d-1, 3.5543098d-1,
     : 4.2954632d-1, 3.7145697d-1, 1.8816934d-1, -5.3511807d-2,
     : -2.5595553d-1, -3.3452327d-1/
      data y51 /
     : 0d0, -3.9903550d-1, -5.4902642d-1, -3.7032566d-1, -3.3095110d-3,
     : 2.9428791d-1, 3.2937930d-1, 1.0543725d-1, -1.8828512d-1,
     : -3.2028165d-1, 0d0, -3.2282650d-1, -4.4417170d-1, -2.9959975d-1,
     : -2.6774507d-3, 2.3808392d-1, 2.6647345d-1, 8.5300526d-2,
     : -1.5232586d-1, -2.5911330d-1, 0d0, -1.2330875d-1, -1.6965849d-1,
     : -1.1443692d-1, -1.0226952d-3, 9.0939966d-2, 1.0178380d-1,
     : 3.2581902d-2, -5.8183303d-2, -9.8972472d-2, 0d0, 1.2330875d-1,
     : 1.6965849d-1, 1.1443692d-1, 1.0226952d-3, -9.0939966d-2,
     : -1.0178380d-1, -3.2581902d-2, 5.8183303d-2, 9.8972472d-2, 0d0,
     : 3.2282650d-1, 4.4417170d-1, 2.9959975d-1, 2.6774507d-3,
     : -2.3808392d-1, -2.6647345d-1, -8.5300526d-2, 1.5232586d-1,
     : 2.5911330d-1, 0d0, 3.9903550d-1, 5.4902642d-1, 3.7032566d-1,
     : 3.3095110d-3, -2.9428791d-1, -3.2937930d-1, -1.0543725d-1,
     : 1.8828512d-1, 3.2028165d-1/
      data y53 /
     : 0d0, -1.3999674d-2, -9.6154653d-2, -2.4864705d-1, -3.9336476d-1,
     : -4.2277566d-1, -2.8087130d-1, -1.5156428d-2, 2.4074596d-1,
     : 3.4594372d-1, 0d0, 4.3261371d-3, 2.9713422d-2, 7.6836163d-2,
     : 1.2155640d-1, 1.3064486d-1, 8.6794004d-2, 4.6835939d-3,
     : -7.4394592d-2, -1.0690249d-1, 0d0, 1.1325974d-2, 7.7790748d-2,
     : 2.0115969d-1, 3.1823877d-1, 3.4203270d-1, 2.2722965d-1,
     : 1.2261808d-2, -1.9476757d-1, -2.7987435d-1, 0d0, -1.1325974d-2,
     : -7.7790748d-2, -2.0115969d-1, -3.1823877d-1, -3.4203270d-1,
     : -2.2722965d-1, -1.2261808d-2, 1.9476757d-1, 2.7987435d-1, 0d0,
     : -4.3261371d-3, -2.9713422d-2, -7.6836163d-2, -1.2155640d-1,
     : -1.3064486d-1, -8.6794004d-2, -4.6835939d-3, 7.4394592d-2,
     : 1.0690249d-1, 0d0, 1.3999674d-2, 9.6154653d-2, 2.4864705d-1,
     : 3.9336476d-1, 4.2277566d-1, 2.8087130d-1, 1.5156428d-2,
     : -2.4074596d-1, -3.4594372d-1/
      data y60 /
     : 1.0171072e+00, 7.1652290d-1, 7.3133063d-2, -3.8042194d-1,
     : -3.2910613d-1, 5.7342740d-2, 3.2877197d-1, 2.1245027d-1,
     : -1.3438157d-1, -3.1784601d-1,  1.0171072e+00, 7.1652290d-1,
     : 7.3133063d-2, -3.8042194d-1, -3.2910613d-1, 5.7342740d-2,
     : 3.2877197d-1, 2.1245027d-1, -1.3438157d-1, -3.1784601d-1, 
     : 1.0171072e+00, 7.1652290d-1, 7.3133063d-2, -3.8042194d-1,
     : -3.2910613d-1, 5.7342740d-2, 3.2877197d-1, 2.1245027d-1,
     : -1.3438157d-1, -3.1784601d-1,  1.0171072e+00, 7.1652290d-1,
     : 7.3133063d-2, -3.8042194d-1, -3.2910613d-1, 5.7342740d-2,
     : 3.2877197d-1, 2.1245027d-1, -1.3438157d-1, -3.1784601d-1, 
     : 1.0171072e+00, 7.1652290d-1, 7.3133063d-2, -3.8042194d-1,
     : -3.2910613d-1, 5.7342740d-2, 3.2877197d-1, 2.1245027d-1,
     : -1.3438157d-1, -3.1784601d-1,  1.0171072e+00, 7.1652290d-1,
     : 7.3133063d-2, -3.8042194d-1, -3.2910613d-1, 5.7342740d-2,
     : 3.2877197d-1, 2.1245027d-1, -1.3438157d-1, -3.1784601d-1/
      data y62 /
     : 0d0, 1.4321472d-1, 4.1286599d-1, 4.9363185d-1, 2.4237694d-1,
     : -1.5358744d-1, -3.5114018d-1, -1.8809782d-1, 1.5390619d-1,
     : 3.2569524d-1, 0d0, 4.4255784d-2, 1.2758261d-1, 1.5254063d-1,
     : 7.4898593d-2, -4.7461130d-2, -1.0850828d-1, -5.8125424d-2,
     : 4.7559628d-2, 1.0064537d-1, 0d0, -1.1586315d-1, -3.3401560d-1,
     : -3.9935656d-1, -1.9608706d-1, 1.2425485d-1, 2.8407838d-1,
     : 1.5217433d-1, -1.2451272d-1, -2.6349299d-1, 0d0, -1.1586315d-1,
     : -3.3401560d-1, -3.9935656d-1, -1.9608706d-1, 1.2425485d-1,
     : 2.8407838d-1, 1.5217433d-1, -1.2451272d-1, -2.6349299d-1, 0d0,
     : 4.4255784d-2, 1.2758261d-1, 1.5254063d-1, 7.4898593d-2,
     : -4.7461130d-2, -1.0850828d-1, -5.8125424d-2, 4.7559628d-2,
     : 1.0064537d-1, 0d0, 1.4321472d-1, 4.1286599d-1, 4.9363185d-1,
     : 2.4237694d-1, -1.5358744d-1, -3.5114018d-1, -1.8809782d-1,
     : 1.5390619d-1, 3.2569524d-1/
      data y71 /
     : 0d0, -5.7323672d-1, -5.1593073d-1, 3.8424259d-2, 3.9192554d-1,
     : 1.6955901d-1, -2.4978896d-1, -2.8241674d-1, 8.6024979d-2,
     : 3.1937046d-1, 0d0, -4.6375824d-1, -4.1739673d-1, 3.1085878d-2,
     : 3.1707443d-1, 1.3717612d-1, -2.0208351d-1, -2.2847994d-1,
     : 6.9595670d-2, 2.5837613d-1, 0d0, -1.7713989d-1, -1.5943136d-1,
     : 1.1873749d-2, 1.2111165d-1, 5.2396615d-2, -7.7189034d-2,
     : -8.7271571d-2, 2.6583180d-2, 9.8690900d-2, 0d0, 1.7713989d-1,
     : 1.5943136d-1, -1.1873749d-2, -1.2111165d-1, -5.2396615d-2,
     : 7.7189034d-2, 8.7271571d-2, -2.6583180d-2, -9.8690900d-2, 0d0,
     : 4.6375824d-1, 4.1739673d-1, -3.1085878d-2, -3.1707443d-1,
     : -1.3717612d-1, 2.0208351d-1, 2.2847994d-1, -6.9595670d-2,
     : -2.5837613d-1, 0d0, 5.7323672d-1, 5.1593073d-1, -3.8424259d-2,
     : -3.9192554d-1, -1.6955901d-1, 2.4978896d-1, 2.8241674d-1,
     : -8.6024979d-2, -3.1937046d-1/
      data y73 /
     : 0d0, -4.2575622d-2, -2.4885375d-1, -4.6932666d-1, -3.9705811d-1,
     : -7.0889627d-3, 3.2785374d-1, 2.5371232d-1, -1.2044765d-1,
     : -3.3189952d-1, 0d0, 1.3156591d-2, 7.6900036d-2, 1.4502992d-1,
     : 1.2269770d-1, 2.1906100d-3, -1.0131238d-1, -7.8401418d-2,
     : 3.7220369d-2, 1.0256259d-1, 0d0, 3.4444402d-2, 2.0132691d-1,
     : 3.7969325d-1, 3.2122676d-1, 5.7350913d-3, -2.6523925d-1,
     : -2.0525758d-1, 9.7444192d-2, 2.6851235d-1, 0d0, -3.4444402d-2,
     : -2.0132691d-1, -3.7969325d-1, -3.2122676d-1, -5.7350913d-3,
     : 2.6523925d-1, 2.0525758d-1, -9.7444192d-2, -2.6851235d-1, 0d0,
     : -1.3156591d-2, -7.6900036d-2, -1.4502992d-1, -1.2269770d-1,
     : -2.1906100d-3, 1.0131238d-1, 7.8401418d-2, -3.7220369d-2,
     : -1.0256259d-1, 0d0, 4.2575622d-2, 2.4885375d-1, 4.6932666d-1,
     : 3.9705811d-1, 7.0889627d-3, -3.2785374d-1, -2.5371232d-1,
     : 1.2044765d-1, 3.3189952d-1/
      data y80 /
     : 1.1631066e+00, 6.0696266d-1, -2.9291611d-1, -3.9403218d-1,
     : 1.6123775d-1, 3.4274692d-1, -8.5649911d-2, -3.2336134d-1,
     : 2.7109515d-2, 3.1803697d-1,  1.1631066e+00, 6.0696266d-1,
     : -2.9291611d-1, -3.9403218d-1, 1.6123775d-1, 3.4274692d-1,
     : -8.5649911d-2, -3.2336134d-1, 2.7109515d-2, 3.1803697d-1, 
     : 1.1631066e+00, 6.0696266d-1, -2.9291611d-1, -3.9403218d-1,
     : 1.6123775d-1, 3.4274692d-1, -8.5649911d-2, -3.2336134d-1,
     : 2.7109515d-2, 3.1803697d-1,  1.1631066e+00, 6.0696266d-1,
     : -2.9291611d-1, -3.9403218d-1, 1.6123775d-1, 3.4274692d-1,
     : -8.5649911d-2, -3.2336134d-1, 2.7109515d-2, 3.1803697d-1, 
     : 1.1631066e+00, 6.0696266d-1, -2.9291611d-1, -3.9403218d-1,
     : 1.6123775d-1, 3.4274692d-1, -8.5649911d-2, -3.2336134d-1,
     : 2.7109515d-2, 3.1803697d-1,  1.1631066e+00, 6.0696266d-1,
     : -2.9291611d-1, -3.9403218d-1, 1.6123775d-1, 3.4274692d-1,
     : -8.5649911d-2, -3.2336134d-1, 2.7109515d-2, 3.1803697d-1/
      data y82 /
     : 0d0, 2.6210731d-1, 5.6960883d-1, 2.9356940d-1, -2.6402938d-1,
     : -3.2003932d-1, 1.3229522d-1, 3.2252245d-1, -4.0988996d-2,
     : -3.2254836d-1, 0d0, 8.0995614d-2, 1.7601881d-1, 9.0717934d-2,
     : -8.1589566d-2, -9.8897588d-2, 4.0881472d-2, 9.9664918d-2,
     : -1.2666296d-2, -9.9672923d-2, 0d0, -2.1204927d-1, -4.6082322d-1,
     : -2.3750263d-1, 2.1360426d-1, 2.5891725d-1, -1.0702908d-1,
     : -2.6092614d-1, 3.3160794d-2, 2.6094710d-1, 0d0, -2.1204927d-1,
     : -4.6082322d-1, -2.3750263d-1, 2.1360426d-1, 2.5891725d-1,
     : -1.0702908d-1, -2.6092614d-1, 3.3160794d-2, 2.6094710d-1, 0d0,
     : 8.0995614d-2, 1.7601881d-1, 9.0717934d-2, -8.1589566d-2,
     : -9.8897588d-2, 4.0881472d-2, 9.9664918d-2, -1.2666296d-2,
     : -9.9672923d-2, 0d0, 2.6210731d-1, 5.6960883d-1, 2.9356940d-1,
     : -2.6402938d-1, -3.2003932d-1, 1.3229522d-1, 3.2252245d-1,
     : -4.0988996d-2, -3.2254836d-1/ 

*  thetab and phib are defined with respect to the C2 symmetry axis, since
*  these are the angles used by the saptpot subr
*
*  grid of theta values w. respect to b inertial axis for input into saptpot subr
      data thetab /
     :   9.0000000d+01, 9.0000000d+01, 9.0000000d+01, 9.0000000d+01
     : , 9.0000000d+01, 9.0000000d+01, 1.0000000d+02, 9.8075873d+01
     : , 9.3075983d+01, 8.6924017d+01, 8.1924127d+01, 8.0000000d+01
     : , 1.1000000d+02, 1.0606336d+02, 9.6066924d+01, 8.3933076d+01
     : , 7.3936645d+01, 7.0000000d+01, 1.2000000d+02, 1.1386033d+02
     : , 9.8888292d+01, 8.1111708d+01, 6.6139669d+01, 6.0000000d+01
     : , 1.3000000d+02, 1.2133400d+02, 1.0145699d+02, 7.8543009d+01
     : , 5.8665998d+01, 5.0000000d+01, 1.4000000d+02, 1.2829737d+02
     : , 1.0369308d+02, 7.6306923d+01, 5.1702633d+01, 4.0000000d+01
     : , 1.5000000d+02, 1.3447751d+02, 1.0552249d+02, 7.4477512d+01
     : , 4.5522488d+01, 3.0000000d+01, 1.6000000d+02, 1.3948424d+02
     : , 1.0688077d+02, 7.3119233d+01, 4.0515760d+01, 2.0000000d+01
     : , 1.7000000d+02, 1.4281861d+02, 1.0771740d+02, 7.2282602d+01
     : , 3.7181394d+01, 1.0000000d+01, 1.8000000d+02, 1.4400000d+02
     : , 1.0800000d+02, 7.2000000d+01, 3.6000000d+01, 0.0000000d+00 /
*  grid of phi values w. respect to b inertial axis for input into saptpot subr
      data phib /
     :   0.0000000d+00, 0.0000000d+00, 0.0000000d+00, 0.0000000d+00
     : , 0.0000000d+00, 0.0000000d+00, 0.0000000d+00, 5.9171456d+00
     : , 9.5197466d+00, 9.5197466d+00, 5.9171456d+00, 1.2372352d-15
     : , 0.0000000d+00, 1.2075617d+01, 1.9093616d+01, 1.9093616d+01
     : , 1.2075617d+01, 2.5538733d-15, 0.0000000d+00, 1.8745053d+01
     : , 2.8770869d+01, 2.8770869d+01, 1.8745053d+01, 4.0510990d-15
     : , 0.0000000d+00, 2.6252994d+01, 3.8590958d+01, 3.8590958d+01
     : , 2.6252994d+01, 5.8877182d-15, 0.0000000d+00, 3.5011057d+01
     : , 4.8578634d+01, 4.8578634d+01, 3.5011057d+01, 8.3621885d-15
     : , 0.0000000d+00, 4.5513129d+01, 5.8739653d+01, 5.8739653d+01
     : , 4.5513129d+01, 1.2153297d-14, 0.0000000d+00, 5.8233250d+01
     : , 6.9058104d+01, 6.9058104d+01, 5.8233250d+01, 1.9278250d-14
     : , 0.0000000d+00, 7.3301524d+01, 7.9496559d+01, 7.9496559d+01
     : , 7.3301524d+01, 3.9793736d-14, 0.0000000d+00, 9.0000000d+01
     : , 9.0000000d+01, 9.0000000d+01, 9.0000000d+01, 6.3434949d+01 /
*
*  use distance units to bohr in saptpot subr
      data iaa / 1 /, zero / 0.d0 /
*  number of vlm coefficients and tolerance of lsq fit
      data nvlm, nthph, nylm / 16, 60, 960 /, tol / 1.d-10 /
*
*  renormalize Y(lam,2) functions
      do i=1,60
        y11t(i) = y11(i) * 2.d0
        y31t(i) = y31(i) * 2.d0
        y51t(i) = y51(i) * 2.d0
        y22t(i) = y22(i) * 2.d0
        y42t(i) = y42(i) * 2.d0
        y62t(i) = y62(i) * 2.d0
        y33t(i) = y33(i) * 2.d0
        y53t(i) = y53(i) * 2.d0
      end do
*
*  determine matrix of ylm's for use in least squares fitting,  use here 
*  only up through lambda = 8. ordering is y00, y20, y40, y60, y11, 
*      y31, y51, y22, y42, y62, y33, y53
      do icol=1,16
         if (icol.eq.1) call dcopy(60,y00,1,ylm(1,icol),1)
         if (icol.eq.2) call dcopy(60,y20,1,ylm(1,icol),1)
         if (icol.eq.3) call dcopy(60,y40,1,ylm(1,icol),1)
         if (icol.eq.4) call dcopy(60,y60,1,ylm(1,icol),1)
         if (icol.eq.5) call dcopy(60,y11t,1,ylm(1,icol),1)
         if (icol.eq.6) call dcopy(60,y31t,1,ylm(1,icol),1)
         if (icol.eq.7) call dcopy(60,y51t,1,ylm(1,icol),1)
         if (icol.eq.8) call dcopy(60,y22t,1,ylm(1,icol),1)
         if (icol.eq.9) call dcopy(60,y42t,1,ylm(1,icol),1)
         if (icol.eq.10) call dcopy(60,y62t,1,ylm(1,icol),1)
         if (icol.eq.11) call dcopy(60,y33t,1,ylm(1,icol),1)
         if (icol.eq.12) call dcopy(60,y53t,1,ylm(1,icol),1)
      enddo
*  compute potential at grid of (theta,phi) values
      do i = 1,60
          th = thetab(i)
          ph = phib(i)
          call heh2osapt(iaa,r,th,ph,energy)
          vcalc(i) = energy
      end do
* reorder vcalc so that all thetas for a given phi occur in column order, e.g.
* first column is all thetas for phi=0, etc
      do iphi=1,6
         indx=(iphi-1)*10+1
         call dcopy(10,vcalc(iphi),6,vcalcx(indx),1)
      end do
*
      lwork=1812
      rcond=1.d-6
      call dgelsd(60,12,1,ylm,60,vcalcx,60,swork,rcond,irank,work,lwork,
     :            iwork,info)
*
*  convert to atomic units and fill vlm coefficients
      econv=1.d0/219474.6d0
*  normalize vv0 so that it multiplies unity instead of 1/sqrt(4*pi)
      vv0 = ( vcalcx(1) * econv ) / s4pi
      call dcopy(11,vcalcx(2),1,vvl,1)
      call dscal(11,econv,vvl,1)
      return
      end
* ----------------------------------------------------------------------
c  saptpot.f
c
c  program obtained from phillip stancil and benyui yang (uga)
c  compute H2O-He SAPT potential [see Patkowski et al., J. Molec.
c  Struct. Theochem 591, 231-241 (2002)]
c
c  original name of subr was arh2osapt (written for Ar-H2O PES, but
c  never published)
c
c...    actual calculation of potential
        subroutine heh2osapt(iaa,r,theta,phi,energy)
c       subroutine arh2osapt(iaa,r,theta,phi,energy,tima,timb)
        implicit real*8 (a-h,o-z)
        logical first
        parameter (mxl=100)
        parameter (maxb=4*mxl,maxp=800)
        parameter (maxc=5000)
        dimension rpt(6)
        common/damp/ ldum(5,mxl),cdum(maxb),ndum
        common/exch/ lex(5,mxl),cex(maxb),nex,irpowex,iwex,iexback,igrx,
     1               idonex,exscale,lmaxex
        common/spher/lsp(5,mxl),csp(maxb),nsp,irpowsp,iwsp,ispback,igrt,
     $               lmaxli
        common/dind/ cd(maxc),ci(maxc),ld(6,maxc),li(6,maxc),
     1               idisp,iind,lmaxi,lmaxd
        data a0 /0.529177249d0/
        data first /.true./
        data xkcal /627.51d0/
        data cm    /219474.63d0/
        data icount /0/
c
c       Modified version of subroutine computing the Ar-H2O potential.
c       The calls to Wormer's almre were replaced by a calculation of
c       spherical harmonics -- subroutine ssh.
c
c       INPUT:
c       ======
c       dimer geometry (R, theta, phi)   -- via header or read in
c       long-range coefficients          -- in files h2o_coefi.dat and h2o_coefd.dat
c       short-range optimized parameters -- in params
c
c       irpowex = 1  now read in rdexch
c       irpowsp = 3  now read in rdlin
c       imode   = i  controls terms computed by potentot
c                 3  only asymptotics
c                 4  both short-range and asymptotics
c
c       Lines deactivated by c1 compute derviatives.
c
c       icount=icount+1
c       if((icount/100)*100 .eq. icount) write(6,*) 'icount=', icount
        pi = dacos(-1.d0)
        d2rad = pi/180.d0
c...    activate next two lines for stand alone program
c       write(6,*)'enter: i=1 or 2 for bohr/ang, R, theta, phi in deg'
c       read(5,*) iaa,r,theta,phi
        do j=1,6
          rpt(j)=0.0d0
        enddo
c...    program uses R in angstroms
        if(iaa.eq.1) rpt(1) = a0*r
        if(iaa.eq.2) rpt(1) = r
        rpt(2) = d2rad*theta
        rpt(3) = d2rad*phi
        if(first) then
*          call rdmom    ! now superfluous call
          call rdcoef   ! read induction and dispersion coeffs
          call rdexch   ! read short-range exponential parms
          call rdlin    ! read short-range linear parms
          first=.false.
        endif
        ndum = nex
        ipow = 2
        do i=1,nex
         cdum(i) = cex(ipow*i-ipow+2)
         do k=1,5
          ldum(k,i) = lex(k,i)
         end do
        end do
c
c       compute table of A_lm functions.  There functions are
c       for this particular case identical to those produced by almre.
c
        lmax = max0(lmaxex,lmaxli,lmaxi,lmaxd)
c       write(6,*) 'lmax', lmax
c       call gclock(tim1)
        call ssh(lmax,rpt(2),rpt(3))
c       call gclock(tim2)
c       tima=tima+tim2-tim1
c
        imode = 4
c       call potentot(rpt,vex,vsp,valas,der,imode,timb)
        call potentot(rpt,vex,vsp,valas,der,imode)
        tot = vsp + valas
        energy=tot*cm/xkcal
c       write(6,*)'Total energy:',energy
c       energy=tot
c       stop
        return
        end

c
c Read the multipole moment components
*  commented subr rdmom out, since listed in calling subr as superfluous
*  (paul dagdigian, aug-2009)
c
*        subroutine rdmom
*        implicit real*8 (a-h,o-z)
*       parameter (mxl=100)
*       parameter (maxb=4*mxl,maxp=800)
*        common/moments/ qa(100,3), qb(100,3),la(maxb),ka(maxb),lb(maxb),
*     1 kb(maxb),nmoma, nmomb
c
c---- read the moments for molecule A
*        open(unit=3,file='moments.a',form='formatted')
c
*        do i=1,1000
*         read(3,*,END=100) la(i),ka(i),(qa(i,j),j=1,3)
*        end do
* 100    continue
*        nmoma = i - 1
*        close(3,STATUS='delete')
c---- read the moments for molecule B
*        open(unit=3,file='moments.b',form='formatted')
c
*        do i=1,1000
*         read(3,*,END=101) lb(i),kb(i),(qb(i,j),j=1,3)
*        end do
* 101    continue
*        nmomb = i - 1
*        close(3,STATUS='delete')
c
*        return
*        end
c
c Read the dispersion ond induction coefficients.
c It is assumed that the induction coefficients are
c complete, i.e., that both A->B and B->A are included,
c with the proper phases after the induct calculation...
c
        subroutine rdcoef
        implicit real*8 (a-h,o-z)
        parameter (maxc=5000)
        common/dind/ cd(maxc),ci(maxc),ld(6,maxc),li(6,maxc),
     1 idisp,iind,lmaxi,lmaxd
c
        open(unit=1,file=
     :    'potdata/h2o_coefd.dat',
     :    form='formatted')
        open(unit=2,file=
     :    'potdata/h2o_coefi.dat',
     :    form='formatted')
c
        lmaxd=0
        nmin=1000
        nmax=0
        do i=1,100000
         read(1,*,end=100)(ld(j,i),j=1,6), c0,c1,ctot
c !!!! update !!!!
c        cd(i) = c0    ! KSz
         cd(i) = ctot  ! KSz
         lmaxd=max0(lmaxd,ld(1,i))
         nmin=min0(nmin,ld(6,i))
         nmax=max0(nmax,ld(6,i))
        end do
 100    continue
        if(nmin.lt.6 .or. nmax.gt.12) then
          write(6,*) 'n out of range in rdcoef-disp',nmin,nmax
          stop
        endif
        idisp = i-1
c
        lmaxi=0
        nmin=1000
        nmax=0
        do i=1,100000
         read(2,*,END=200)(li(j,i),j=1,6),c0,c1,ctot
c !!!! update !!!!
c        ci(i) = c0    ! KSz
         ci(i) = ctot  ! KSz
         lmaxi=max0(lmaxi,li(1,i))
         nmin=min0(nmin,li(6,i))
         nmax=max0(nmax,li(6,i))
        end do
 200    continue
        if(nmin.lt.6 .or. nmax.gt.12) then
          write(6,*) 'n out of range in rdcoef-ind',nmin,nmax
          stop
        endif
        iind = i-1
c
        close(1)
        close(2)
        return
        end
c
c Subroutine for the calculation of the multipole
c part of the electrostatic energy for a given
c dimer conformation. The angles are expressed
c in radians. 
c
        subroutine elst(rpt,el)
        implicit real*8 (a-h,o-z)
        complex*16 a, alm
       parameter (mxl=100)
       parameter (maxb=4*mxl,maxp=800)
        dimension rrev(20), oa(3), ob(3), oc(2), rpt(6), el(20)
        common/factorial/ fac(0:40)
        common/moments/ qa(100,3), qb(100,3),la(maxb),ka(maxb),lb(maxb),
     1 kb(maxb),nmoma, nmomb
        common/realm/ almre(0:50,0:50)
        data Pi /3.1415926535897932D0/
        data efact /627.51d0/
        data a0 /0.529177249d0/
        data oa(1) /0.d0/, oc /0.d0,0.d0/
c
        r = rpt(1)
        oa(2) = rpt(2)
        oa(3) = rpt(3)
        ob(1) = rpt(4)
        ob(2) = rpt(5)
        ob(3) = rpt(6)
        rrev(1) = 1.d0/(R/a0)
        do i=2,20
         rrev(i) = rrev(i-1)*rrev(1)
        end do
c
         els = 0.d0
        do i=1,20
         el(i) = 0.d0
        end do
c
        do ia=1,nmoma
        do ib=1,nmomb
         phase = fac(2*la(ia)+2*lb(ib)+1)
         phase = phase/(fac(2*la(ia))*fac(2*lb(ib)))
         phase = (-1.d0)**la(ia) * dsqrt(phase)
         ll = la(ia) + lb(ib)
         rrr = rrev(ll+1)
c...     replace call to almre by substitution of precomputed value
c        aa  = almre(la(ia),ka(ia),lb(ib),kb(ib),ll,OA,OB,OC)
         if(ka(ia) .lt. 0) then
           write(6,*) 'ka .lt. 0 in elst'
           stop
         endif
         aa  = almre(la(ia),ka(ia))
c calculate the product of multipole moments at the proper level
c (currently -- MP3-resp)
         qpr = qa(ia,1)*qb(ib,1) + qa(ia,1)*(qb(ib,2)+qb(ib,3))
         qpr = qpr + (qa(ia,2)+qa(ia,3))*qb(ib,1)
c limitation to R^-11 in electrostatics....
         if((ll+1).le.11)
c     1   el(ll+1) = el(ll+1) + phase*rrr*real(a)*qpr
     1   el(ll+1) = el(ll+1) + phase*rrr*aa*qpr
        end do
        end do
c
        do i=1,20
         el(i) = efact*el(i)
        end do
c
c That's all
c
        return
        end
c
c
c   Calculate the total damped asymptotics from elst and dispind..
c   th1, th2, phi in radians...., R in Angstroms.
c
        subroutine asymp(rpt,value,valder)
        implicit real*8 (a-h,o-z)
       parameter (mxl=100)
       parameter (maxb=4*mxl,maxp=800)
        dimension rpt(6),oa(3), ob(3), oc(2)
        dimension en(20), el(20), ei(20), ed(20)
      common/damp/ ldum(5,mxl),cdum(maxb),ndum
      common/realm/ almre(0:50,0:50)
      data oa(1) /0.d0/, oc /0.d0,0.d0/
c
        r = rpt(1)
        oa(2) = rpt(2)
        oa(3) = rpt(3)
        ob(1) = rpt(4)
        ob(2) = rpt(5)
        ob(3) = rpt(6)
c
c---- call the asymptotics procedures to be changed....
c
        call dispind(rpt,ed,ei)
c       call elst(rpt,el) ! KSz deacive for ArH2O
c        call shrtas(R,th1,th2,phi,srval)
c
c---- Compute the damping constant for a given geometry
c
        dump = 0.d0
        do i=1,ndum
         la = ldum(1,i)
         ka = ldum(2,i)
         lb = ldum(3,i)
         kb = ldum(4,i)
         l  = ldum(5,i)
c glamc changed into glam to test the real alm version...
c...     replace calls to almre by substitution of precomputed value
c        glam = almre(la,ka,lb,kb,l,oa,ob,oc)
         if(ka .lt. 0) then
           write(6,*) 'ka .lt. 0 in asympt'
           stop
         endif
         glam = almre(la,ka) ! KSz
c -- symmetrize iif molecules identical...
         if(iperm.eq.1) then
          iphase = (-1)**(la+lb)
c         glam = glam + iphase*almre(lb,kb,la,ka,l,oa,ob,oc)
          if(kb .lt. 0) then
            write(6,*) 'kb .lt. 0 in asympt'
            stop
          endif
          glam = glam + iphase*almre(lb,kb)
         endif
c         glam = real(glamc)
         dump = dump - cdum(i)*glam
        end do
c
c--- Damping factor ready, proceed to the asymptotics...
c
c       do i=1,20
        do i=6,12
c        en(i) = el(i) + ed(i) + ei(i)
         en(i) = ed(i) + ei(i)
        end do
        value = 0.d0
        valder = 0.d0
        do i=6,12 ! KSz changed from 1,20
         ddd = d(i,dump,r)
         value = value + ddd*en(i)
c1     valder = valder + en(i)*(dd(i,dump,r)-i*ddd/r)
        end do
       return
       end
c
      function d(n,beta,r)
c
c     calculate the damping factor (small R correct)
c
      implicit real*8 (a-h,o-z)
      br=beta*r
      sum=1.0d0
      term=1.0d0
      ncn=n
      do i=1,ncn
        term=term*br/i
        sum=sum+term
      enddo
      d=1.0d0 - dexp(-br)*sum
c     in case of d --> 0 use
c     d=1.0d0 - dexp(-br)*sum = sum_m=ncn+1^\infty br^m/m!
      if(dabs(d).lt.1.0d-8) then
        d=0.0d0
        do i=ncn+1,1000
          term=term*br/i
          d=d+term
          if(term/d .lt. 1.0d-8) go to 111
        enddo
  111 continue
      d=d*dexp(-br)
      endif
c     write(6,'(i4,2f10.5,e20.10)') n,beta,r,d
      return
      end
c
      function dd(n,b,r)
      implicit real*8 (a-h,o-z)
      common/factorial/ f(0:40)
      br = b*r
      dd = b*dexp(-br)*(br)**n/f(n)
      return
      end
c
c Subroutine for the calculation of the multipole
c part of the dispersion and induction energies for a given
c dimer conformation. 
c
c       Ar-H2O version:
c       *** NOTICE ***
c       For Ar-H2O C_n^lm = C_n^l,-m,  however,  Polcor suite
c       includes both types of coefficients in output.
c       To avoid calculation of Alm for negative m we use
c       A_lm = A_l,-m for this symmetry.
c
        subroutine dispind(rpt,eld,eli)
        implicit real*8 (a-h,o-z)
        parameter (maxc=5000)
        dimension rr(20),oa(3),ob(3),oc(2), rpt(6), eld(2), eli(20)
        complex*16 a, alm
c--- common/dind/ contains dispersion and induction coefs.
        common/dind/ cd(maxc),ci(maxc),ld(6,maxc),li(6,maxc),
     1               idisp,iind
        common/realm/ almre(0:50,0:50)
        data Pi /3.1415926535897932D0/
        data a0 /0.529177249d0/
        data efact /627.51d0/        
        data oa(1) /0.d0/, oc /0.d0,0.d0/
c
        r = rpt(1)
        oa(2) = rpt(2)
        oa(3) = rpt(3)
        ob(1) = rpt(4)
        ob(2) = rpt(5)
        ob(3) = rpt(6)
        rrev = 1.d0/(R/a0)
        rr(1) = rrev
        do i=2,12
         rr(i) = rr(i-1)*rrev
        end do
        do i=1,20
         eld(i) = 0.d0
         eli(i) = 0.d0
        end do
c
        do i=1,idisp
         la = ld(1,i)
         ka = ld(2,i)
         lb = ld(3,i)
         kb = ld(4,i)
         l  = ld(5,i)
         n  = ld(6,i)
c...     replace calls to almre by substitution of precomputed value
c        aa = almre(la,ka,lb,kb,l,oa,ob,oc)
         aa = almre(la,iabs(ka)) ! KSz: notice iabs.
c         eld(n) = eld(n) + cd(i)*rr(n)*real(a)
         eld(n) = eld(n) + cd(i)*rr(n)*aa
        end do
c
        do i=1,iind
         la = li(1,i)
         ka = li(2,i)
         lb = li(3,i)
         kb = li(4,i)
         l  = li(5,i)
         n  = li(6,i)
c        aa = almre(la,ka,lb,kb,l,oa,ob,oc)
         aa = almre(la,iabs(ka)) ! KSz: notice iabs.
c         eli(n) = eli(n) + ci(i)*rr(n)*real(a)
         eli(n) = eli(n) + ci(i)*rr(n)*aa
        end do
c       do i=1,20
        do i=6,12
         eli(i) = -efact*eli(i)
         eld(i) = -efact*eld(i)
        end do
c
c That's all
c
        return
        end
c
c  A short subroutine to calculate the total interaction
c  potential from the exchange, "spherical", and
c  asymptotic parts...
c
        subroutine potentot(rpt,vex,vsp,valas,der,imode)
c       subroutine potentot(rpt,vex,vsp,valas,der,imode,timb)
        implicit real*8 (a-h,o-z)
        complex*16 glamc, alm
       parameter (mxl=100)
       parameter (maxb=4*mxl,maxp=800)
        dimension rr(10), rpt(6), oa(3), ob(3), oc(2)
      common/exch/lex(5,mxl),cex(maxb),nex,irpowex,iwex,iexback,igrx,
     1 idonex,exscale
      common/spher/lsp(5,mxl),csp(maxb),nsp,irpowsp,iwsp,ispback,igrt
      common/realm/ almre(0:50,0:50)
      data xkcal /627.51d0/
      data cm    /219474.63d0/
      data oa(1) /0.d0/, oc /0.d0,0.d0/
c
        r = rpt(1)
        oa(2) = rpt(2)
        oa(3) = rpt(3)
        ob(1) = rpt(4)
        ob(2) = rpt(5)
        ob(3) = rpt(6)
        vex = 0.d0
        vsp = 0.d0
        valas = 0.d0
        if(imode.ne.3) then
        rr(1) = 1.d0
        do i=2,max(irpowex,irpowsp)+1
         rr(i) = rr(i-1)*r
        end do
c---- Compute the exchange part
        ipow = irpowex+1
c       write(6,*) 'ipow:', ipow
c       write(6,*) 'iexback:', iexback
        rrev = r**(-iexback)
        damp = 0.d0
        do i=1,nex
         la = lex(1,i)
         ka = lex(2,i)
         lb = lex(3,i)
         kb = lex(4,i)
         l  = lex(5,i)
c glamc changed into glam to test the real alm version...
c...     replace calls to almre by substitution of precomputed value
c        glam = almre(la,ka,lb,kb,l,oa,ob,oc)
         if(ka .lt. 0) then
           write(6,*) 'ka .lt. 0 in potentot'
           stop
         endif
         glam = almre(la,ka)
c --- symmetrize if molecules identical...
c        write(6,*) 'iperm:', iperm
         if(iperm.eq.1) then
          iphase = (-1)**(la+lb)
c         glam = glam + iphase*almre(lb,kb,la,ka,l,oa,ob,oc)
          if(kb .lt. 0) then
            write(6,*) 'kb .lt. 0 in potentot'
            stop
          endif
          glam = glam + iphase*almre(lb,kb)
         endif
c         glam = real(glamc)
c1       damp = damp + glam*cex(2*i)
         do k=1,ipow
         vex = vex + glam*cex(ipow*i-ipow+k)*rr(k)*rrev
c        write(6,'(2i2,4f16.8)') i,k,glam,cex(ipow*i-ipow+k),
c    $                 rr(k),vex
         end do
        end do
        vex = dexp(vex)
        if(imode.gt.1) then
c---- Compute the "spherical" part
        ipow = irpowsp + 1
c       write(6,*) 'ipow:', ipow
        rrev = r**(-ispback)
        rrev1 = rrev/r
        vsp = 0.d0
        vsp1 = 0.d0
        do i=1,nsp
         la = lsp(1,i)
         ka = lsp(2,i)
         lb = lsp(3,i)
         kb = lsp(4,i)
         l  = lsp(5,i)
c glamc changed into glam to test the real alm version...
c        glam = almre(la,ka,lb,kb,l,oa,ob,oc)
         if(ka .lt. 0) then
           write(6,*) 'ka .lt. 0 in potentot'
           stop
         endif
         glam = almre(la,ka)
c --- symmetrize if molecules identical....
         if(iperm.eq.1) then
          iphase = (-1)**(la+lb)
c         glam = glam + iphase*almre(lb,kb,la,ka,l,oa,ob,oc)
          if(kb .lt. 0) then
            write(6,*) 'kb .lt. 0 in potentot'
            stop
          endif
          glam = glam + iphase*almre(lb,kb)
         endif
c         glam = real(glamc)
         do k=1,ipow
c---- skip dividing by r when done with experiments, 
c---- update gradient if proved successful
         vsp = vsp + glam*csp(ipow*i-ipow+k)*rr(k)*rrev
c        write(6,'(2i2,4f16.8)') i,k,glam,csp(ipow*i-ipow+k)*cm/xkcal,
c    $                 rr(k),vsp*cm/xkcal
c1       vsp1 = vsp1+ (k-ispback-1)*glam*csp(ipow*i-ipow+k)*rr(k)*rrev1
         end do
        end do
        vsp = vex*vsp
c1      vsp1 = vsp1*vex + damp*vsp
        endif
        endif
c----- Compute the asymptotics...
c       call gclock(tim1)
        if(imode.ge.3)
     1  call asymp(rpt,valas,valdr)
c       call gclock(tim2)
c       timb=timb+tim2-tim1
c1      der = vsp1 + valdr
        return
        end
c
      subroutine rdexch
c
c     reads in the parameters of the expansion of exponential
c      
      implicit real*8 (a-h,o-z)
      parameter (mxl=100)
      parameter (maxb=4*mxl,maxp=800)
      common/exch/lex(5,mxl),cex(maxb),nex,irpowex,iwex,iexback,igrx,
     1           idonex,exscale,lmaxex
c
        open(unit=7,file=
     :    'potdata/h2o_params.dat',
     :    form='formatted')
c
c       read nex     - number of angular functions in expansion of exponent
c            irpowex - largest power of R in exponent
c
        read(7,*) nex,irpowex
        lmaxex=0
        do i=1,nex
         read(7,*)num,(lex(j,i),j=1,5)
         lmaxex=max0(lmaxex,lex(1,i))
        end do
        do i=1,(irpowex+1)*nex
         read(7,*) num, cex(i)
        end do
 10     format(6i4)
       return
       end
c
      subroutine rdlin
c
c     reads in the parameters of the expansion of linear factor
c
      implicit real*8 (a-h,o-z)
      parameter (mxl=100)
      parameter (maxb=4*mxl,maxp=800)
      common/spher/lsp(5,mxl),csp(maxb),nsp,irpowsp,iwsp,ispback,igrt,
     $             lmaxli
c
      read(7,*) nsp,irpowsp
      lmaxli=0
      do k=1,irpowsp+1
      do i=1,nsp
       read(7,*) (lsp(j,i),j=1,5),csp((irpowsp+1)*i-irpowsp-1+k)
       lmaxli=max0(lmaxli,lsp(1,i))
      end do
      end do
      return
      end
      subroutine ssh(lmax,theta,phi)
c     computes symmetrized spherical harmonics
c     which for Ar-H2O are equivalent to real part of
c     the function Alm.
c
c     theta and phi in radians
c
      implicit real*8 (a-h,o-z)
      common/realm/ almre(0:50,0:50)
      common/factorial/ fac(0:40)
      dimension cosm(0:8)
      if(lmax.gt.50) then
        write(6,*) 'lmax gt 50'
        stop
      endif
      nmax=2*lmax
      if(nmax.gt.54) then
        write(6,*) 'nmax gt 40'
        stop
      endif
      call fct(nmax)
      x=dcos(theta)
c...  Plm is Robert's program computing associated Legendre polynomials.
c...  Normalization was changed to standard,  i.e., P_l0(1) = 1.
      call xplm(almre,x,lmax)
      do m=0,lmax
        cosm(m)=dcos(m*phi)
      enddo
      do l=0,lmax
c...    now we need only even m but if odd needed,  change the line below
        do m=0,l,2
c         write(6,*) l,m, almre(l,m)
          almre(l,m) = almre(l,m)*(-1)**(l+m)
     $                 *dsqrt(fac(l-m)/fac(l+m)/(2*l+1))
     $                 *cosm(m)
c         write(6,*) 'TESST',l,m, almre(l,m)
        enddo
      enddo
      return
      end 
c 
c Compute the set of associated Legendre polynomials P_lm
c for l=0,1,...,lmax, and m=0,1,...,l. First the standard
c polynomials
c
c   P^m_l(x) = (1/2^l l!)(1-x^2)^(m/2) (d^(l+m)/d x^(l+m))(x^2 -1)^l
c
c are computed, and then multiplied by
c
c  (-1)^m sqrt[(2l+1)(l-m)!/2(l+m)!]/sqrt(2Pi)
c
c to get the P_lm polynomials....
c
        subroutine xplm(p,x,lmax)
        implicit real*8 (a-h,o-z)
        dimension p(0:50,0:50)
        common/factorial/ fact(0:40)
c inverse of dsqrt(2Pi)
        data twopinv /0.3989422804014d0/
c
c starting value
c
        p(0,0) = 1.d0
        u = dsqrt(1-x*x)
c
c compute the diagonal elements
c
        do l=1,lmax
         p(l,l) = (2*l-1)*p(l-1,l-1)*u
        end do
c
c compute P_lm along the columns with fixed m
c
        do m = 0,lmax-1
        do l = m,lmax-1
         if((l-1).lt.m) then
           pp = 0
         else
           pp = p(l-1,m)
         endif
         p(l+1,m) = ((2*l+1)*x*p(l,m)-(l+m)*pp)/(l-m+1)
        end do 
        end do
c
c Renormalize values...
c
c       do l=0,lmax
c       mm = 1
c       do m=0,l
c        dnorm = fact(l-m)*(2*l+1)/(2*fact(l+m))
c        p(l,m) = mm*twopinv*dsqrt(dnorm)*p(l,m)
c        mm = -mm
c       end do
c       end do
c
        return
        end
c
c compute the matrix of N!
c
        subroutine fct(nmax)
        implicit real*8 (a-h,o-z)
        common/factorial/ f(0:40)
c       
        f(0) = 1.d0
        do i=1,nmax
         f(i) = f(i-1)*i
        end do
        return
        end


