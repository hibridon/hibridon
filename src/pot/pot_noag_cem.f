* System:  NO(X 2Pi)+Ag(111), CEM PES's
* Reference: A. E. DePristo and M. H. Alexander, J. Chem. Phys. 94, 8454 (1991).

      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(15)
1      print *, ' r (bohr)'
      rshift=0.5
      xfact=0.8
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vsum',/,9(1pe16.8),/,
     :    '  vdif',/,7e16.8)
      goto 1
99    end
      include "common/syusr"
      include "common/bausr"
      include "common/ground"
* --------------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parpot"
      include "common/parbas"
      common /coselb/ ibasty
      potnam='DEPRISTO-ALEXANDER NO-Ag(111)'
      ibasty=3
      lammin(1)=1
      lammax(1)=8
      lammin(2)=2
      lammax(2)=8
      mproj(1)=0
      mproj(2)=2
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
*  -----------------------------------------------------------------------
      subroutine pot (vv0, z)
* ----------------------------------------------------
*   subroutine to calculate the z-dependent coefficients in the 
*  NO(X 2PI)-Ag(111) potential of DePristo and Alexander

*  vv0 contains the isotropic term [V  (Z)]
*                                    00
*  vvl(i = 1,8) contains the potential expansion coefficients 
*  for the anistropic components in the sum potential for lam = 1:8
*  that is vvl(i)=V  (Z)
*                  i0
*  vvl(i = 9:15) contains the expansion coefficients for the
*  difference potential for lam = 2:8
*  that is vvl(i)=V  (Z), where l=i-7
*                  l2
*  variable in common block /conlam/
*    nlam:      the total number of angular coupling terms used
*               this should be 15
*    nlammx:    the maximum number of anisotropic terms allowed

*  variable in common block /covvl/
*    vvl        array to store z-dependence of each angular term in the 
*               potential                               

*  written by:  m. h. alexander
*  latest revision date:  16-jun-1990
* ------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension coefa(16,8), coefb(16,8), coefc(16,8)
      dimension pp(16), vpot(16,3), vsum(9,3), vdif(7,3)
      dimension v0exp(9,9), v2exp(7,7), vleg(16,3)
      common /covvl/ vvl(15) 
      common /coconv/ econv, xmconv, ang2
 
*  iselect    if iselect=1; return just the atop potential
*             if iselect=2; return just the bridge potential
*             if iselect=3; return just the center potential
*             if iselect=4; return the site-averaged potential

      data iselect / 4/
      data zero, half, quart, three/0.d0, 0.5d0, 0.25d0, 3.d0/
      
* coefficients for the atop site, first 10 are for V-, last 6 are for V+

      data coefa /
     : 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.3d0, 1.2d0, 1.d0, 1.d0, 1.d0,
     : 1.2d0, 1.2d0, 1.2d0, 1.2d0, 1.1d0, 1.d0,

     : 2.8585448d0, 2.3558802d0, 2.4857090d0, 2.4013132d0, 2.6477683d0,
     : 3.2229382d0, 3.3618559d0, 3.0784431d0, 2.7951324d0, 2.2456580d0,
     : 2.2879875d0, 2.4008073d0, 2.3275859d0, 2.6473418d0, 2.7286295d0,
     : 2.8382770d0,

     : 9.6606505d-1, 8.7249799d-1, 8.6113452d-1, 8.5513346d-1,
     : 8.2890123d-1, 8.7747474d-1, 8.8565759d-1, 8.8566387d-1,
     : 9.2170319d-1, 7.7629755d-1, 8.5703213d-1, 9.0630810d-1,
     : 8.3738965d-1, 8.9991307d-1, 8.7783600d-1, 8.8993423d-1,

     : 5.6714953d0, 4.9792302d0, 4.9576375d0, 5.9534608d0, 7.1339279d0,
     : 5.0515218d0, 5.1614905d0, 5.5015663d0, 5.1403008d0, 6.1702585d0,
     : 6.3113187d0, 6.4468154d0, 7.3691746d0, 6.1496206d0, 4.9624900d0,
     : 5.4481485d0,

     : 1.1500878d+9, 1.8336955d+8, 1.3960941d+8, 4.2128734d+7,
     : 4.5090415d+7, 3.5252033d+8, 2.1630840d+9, 2.6797992d+9,
     : 1.3959588d+9, 1.4031610d+8, 7.6494519d+7, 4.0950075d+7,
     : 2.1385072d+7, 7.6663328d+7, 3.2237821d+8, 1.2257388d+9,

     : 2.3137271d+6, 8.5247014d+5, 6.5781463d+5, 4.2380000d+5,
     : 2.4860562d+5, 3.8223607d+5, 6.6716363d+5, 1.1979052d+6,
     : 1.9769398d+6, 3.2300580d+5, 5.5073560d+5, 6.5794840d+5,
     : 4.1714362d+5, 7.7489303d+5, 9.7128577d+5, 1.3818845d+6,

     : -4.0855719d+5, -1.6133334d+5, -1.2636688d+5, -9.5101201d+4,
     : -6.9875990d+4, -1.1605453d+5, -1.6337062d+5, -2.2218718d+5,
     : -3.5119917d+5, -5.7534765d+4, -1.0470818d+5, -1.3305373d+5,
     : -8.9343956d+4, -1.5680153d+5, -1.8020472d+5, -2.4058045d+5,

     : -1.9292672d+7, -1.3488460d+7, -9.6979146d+5, -1.0613173d+6,
     : -7.2683882d+6, -1.1970106d+7, -1.2763912d+7, 3.9049753d+6,
     : -1.0230041d+7, 2.2023477d+7, 9.9465648d+6, 9.3150852d+6,
     : -8.3937860d+6, 3.4814822d+6, 4.4660156d+6, 9.2157637d+6/

* coefficients for the two fold bridge site, first 10 are for V-, last 6 are for V+
      data coefb /
     : 1.2d0, 1.2d0, 1.2d0, 1.2d0, 1.2d0, 1.2d0, 1.2d0, 1.2d0, 1.2d0,
     : 1.2d0, 1.2d0, 1.2d0, 1.2d0, 1.2d0, 1.2d0, 1.2d0,

     : 4.9293258d-1, 5.2067901d-1, 4.5402182d-1, 4.9639912d-1,
     : 5.2304463d-1, 5.4585953d-1, 5.2637647d-1, 5.7871158d-1,
     : 5.3785825d-1, 4.7985816d-1, 5.0671779d-1, 5.3547361d-1,
     : 5.6454285d-1, 3.8341013d-1, 3.5999828d-1, 5.7590346d-1,

     : 1.6565818d0, 1.5913580d0, 1.7307263d0, 1.6555184d0, 2.0500315d0,
     : 1.8673908d0, 1.8524646d0, 1.5346307d0, 1.5870437d0, 1.6921836d0,
     : 1.0199074d0, 1.3504087d0, 1.4176228d0, 1.7364820d0, 1.0243073d0,
     : 1.1037072d0,

     : 6.3625289d0, 4.9777778d0, 5.9796519d0, 5.6569029d0, 5.3253428d0,
     : 5.5130101d0, 5.8116543d0, 4.7871794d0, 4.7109165d0, 4.6995052d0,
     : 3.4868988d0, 4.3058477d0, 4.1564190d0, 5.5547420d0, 3.3140588d0,
     : 4.9073984d0,

     : -1.5879509d+4, -2.7749708d+4, -1.2151247d+4, -1.6679088d+4,
     : -2.2642811d+4, -3.0328465d+4, -2.8867976d+4, -5.9862634d+4,
     : -4.3640235d+4, -1.9402591d+4, -1.1124109d+4, -2.7149459d+4,
     : -3.5690858d+4, -6.5648815d+3, -2.7213703d+3, -5.4888111d+4,

     : -1.9828967d+6, -3.4939374d+5, 1.6534318d+6, 3.6043987d+6,
     : 6.0261815d+6, 8.1274765d+6, 3.2308463d+6, 8.3630735d+5,
     : -9.9464302d+5, -1.5364288d+6, 1.4240107d+6, 2.7685458d+6,
     : 3.0224888d+6, 5.1658611d+6, 2.0360124d+6, 1.7188263d+6,

     : 2.1212221d+6, 1.4294300d+6, 1.1845632d+6, -1.8293481d+5,
     : 6.0028900d+5, -3.0004617d+5, 2.2951153d+6, 9.1476972d+5,
     : 1.5642637d+6, 2.2858109d+6, -2.1606582d+5, -3.7145189d+5,
     : -3.7731494d+5, 7.6227558d+4, -3.3435201d+5, -1.7953098d+5,

     : 1.7485930d+7, -7.9611413d+6, 1.2026072d+7, 1.2429503d+7,
     : 1.2572140d+7, 1.4571608d+7, 1.5024689d+7, -1.1928631d+7,
     : -1.3564250d+7, -7.9665827d+6, 1.5096066d+7, -1.0986529d+7,
     : -1.2244710d+7, 2.1651350d+7, 2.5744286d+7, -6.0853738d+6/

* coefficients for the center site, first 10 are for V-, last 6 are for V+
      data coefc /
     : 1.3d0, 1.5d0, 1.3d0, 1.d0, 1.d0, 1.5d0, 8.0000000d-1, 1.5d0,
     : 1.5d0, 1.d0, 8.0000000d-1, 1.d0, 1.5d0, 1.5d0, 1.3d0, 1.3d0,

     : 5.9206854d-1, 4.8211355d-1, 5.0083443d-1, 4.3162945d-1,
     : 4.6083676d-1, 2.5162946d-1, 1.0148686d0, 6.0993908d-1,
     : 6.8381173d-1, 5.1655333d-1, 1.0055724d0, 4.9566139d-1,
     : 5.1757318d-1, 5.0916757d-1, 5.7708333d-1, 5.1291120d-1,

     : 1.2909739d0, 1.4251563d0, 1.4093689d0, 1.0325265d0, 1.1386868d0,
     : 1.0603548d0, 8.1579456d-1, 1.3322932d0, 1.4243109d0,
     : 1.4768801d0, 4.3979766d-1, 9.6542541d-1, 1.3669254d0,
     : 1.3318004d0, 1.2152778d0, 1.5044678d0,

     : 4.0219577d0, 3.8637067d0, 3.7086225d0, 3.4333523d0, 4.5131442d0,
     : 4.4123470d0, 3.7895744d0, 3.6871304d0, 6.1343128d0, 5.2366869d0,
     : 3.8633199d0, 4.9151075d0, 4.3070811d0, 4.1664637d0, 4.8777778d0,
     : 4.7004567d0,

     : -5.2939967d+4, -1.9809274d+4, -2.2362215d+4, -9.3354592d+3,
     : -1.4326361d+4, -1.9210918d+3, 1.0061153d+6, -7.7678721d+4,
     : -8.7588970d+4, -2.7921653d+4, 7.0503426d+5, -1.3315369d+4,
     : -2.5879416d+4, -2.6159040d+4, -5.6782287d+4, -3.3776902d+4,

     : -1.9641389d+5, -3.7854512d+5, 1.7404331d+5, 7.3691134d+5,
     : 1.0039207d+6, 1.1244527d+6, -9.0691549d+4, -1.2596253d+5,
     : -4.3675957d+5, -6.4259432d+5, -4.4213572d+4, 6.7090224d+5,
     : 8.5367757d+5, 9.8684780d+5, 7.7706394d+5, -5.1927850d+5,

     : 5.5400296d+5, 7.2823893d+5, 4.7028454d+5, -1.2312854d+5,
     : -2.1837767d+5, -2.6959662d+5, -4.6702160d+4, 6.2191440d+5,
     : 8.0546863d+5, 9.4401782d+5, 3.4369411d+3, -1.2059807d+5,
     : 6.1710921d+4, 8.2825055d+3, 7.2977200d+4, 9.0268158d+5,

     : 9.0350398d+6, 4.2489671d+6, 5.9321616d+6, 6.3568979d+6,
     : -6.5880019d+6, -2.6214274d+6, 3.0065304d+7, 5.6211583d+6,
     : 3.0642294d+7, -8.1481983d+6, 2.5813038d+7, -8.5521428d+6,
     : -2.0424441d+6, -2.5152705d+6, -5.2242232d+6, -7.2568474d+6/

*  v0exp converts the potential at the 9 orientation angles
*  (0:pi/8:pi)  into the first 9 regular legendre expansion terms (v0 ... v8)
*  v0exp is the inverse of the matrix with columns corresponding to
*  the 9 angles, and rows corresponding to the first 9 legendre
*  polynomials
*  in other words v0exp=inv(v0mat), where v0mat(i,j)=p (ang(j))
*                                                     i
      data v0exp /
     : 7.93650780852878d-03, 7.31093233835851d-02, 1.39682537577978d-01,
     : 1.80858927057298d-01, 1.96825395169334d-01, 1.80858929060267d-01,
     : 1.39682541273573d-01, 7.31093284353578d-02, 7.93651023407854d-03,
     : 2.38095233409424d-02, 2.02632622978926d-01, 2.96311410412662d-01,
     : 2.07635150128442d-01, 7.71305340277810d-09,
     : -2.07635138706197d-01, -2.96311410679779d-01,
     : -2.02632634401169d-01, -2.38095307868796d-02,
     : 4.18470411350253d-02, 2.80918020064780d-01, 1.78932183226667d-01,
     : -2.57829987513341d-01, -4.87734483521231d-01,
     : -2.57830016011649d-01, 1.78932159945856d-01,
     : 2.80918029636257d-01, 4.18470530376361d-02, 6.06060592575147d-02,
     : 2.90373987044323d-01, -1.65705819195380d-01,
     : -5.53211633416018d-01, -2.80221719265120d-08,
     : 5.53211627920203d-01, 1.65705864239898d-01,
     : -2.90373981548508d-01, -6.06060762798620d-02,
     : 8.79120861001281d-02, 2.08285569842524d-01,
     : -4.89110876663031d-01, -1.44349529644506d-01,
     : 6.74525468594754d-01, -1.44349454843105d-01,
     : -4.89110911826977d-01, 2.08285542517218d-01,
     : 8.79121059229947d-02, 1.17216113601702d-01, 5.83061242880968d-02,
     : -5.52562079546045d-01, 5.73937671524467d-01,
     : 5.23905479895759d-08, -5.73937730130473d-01,
     : 5.52562052013123d-01, -5.83060656820945d-02,
     : -1.17216138459325d-01, 2.03174592250444d-01,
     : -2.44053368724218d-01, -1.47763386103429d-01,
     : 5.39580127463835d-01, -7.01875910518908d-01,
     : 5.39580061463744d-01, -1.47763292446648d-01,
     : -2.44053429149969d-01, 2.03174605765149d-01,
     : 2.98368303799842d-01, -5.51312734311346d-01,
     : 4.21956488328764d-01, -2.28361188236892d-01,
     : -3.20814302152483d-08, 2.28361240916468d-01,
     : -4.21956505573245d-01, 5.51312681631777d-01,
     : -2.98368254473938d-01, 1.59129772705873d-01,
     : -3.18259544566671d-01, 3.18259541961814d-01,
     : -3.18259537363286d-01, 3.18259530276053d-01,
     : -3.18259519669259d-01, 3.18259503054198d-01,
     : -3.18259471438869d-01, 1.59129725040146d-01/

*  v2exp converts the potential at the 7 orientation angles
*  (pi/8:pi/8:7pi/8) into the first 7 reduced rotation matrix element d2l
*  expansion terms (v22...v82)
*  v2exp is the inverse of the matrix with columns corresponding to
*  the 7 angles, and rows corresponding to the first 7 reduced
*  rotation matrix elements
*                                                     i
*  in other words v2exp=inv(v2mat), where v2mat(i,j)=d  (ang(j))
*                                                     20  

      data v2exp /
     : 3.45494799896523d-02, 2.12077027958067d-01, 4.74435390563992d-01,
     : 6.00884924586591d-01, 4.74435406286922d-01, 2.12077045006018d-01,
     : 3.45494863940010d-02, 9.99240203311050d-02, 4.69452516651283d-01,
     : 5.68368220094310d-01, 2.52058069214459d-08,
     : -5.68368202486280d-01, -4.69452541664563d-01,
     : -9.99240379391353d-02, 2.04128466892967d-01,
     : 6.11605343863769d-01, 1.82734990822429d-02,
     : -7.02587957480906d-01, 1.82734266379968d-02,
     : 6.11605347523234d-01, 2.04128498047839d-01, 3.33992522045737d-01,
     : 4.80443122885568d-01, -6.58372092487126d-01,
     : -5.15689978057847d-08, 6.58372137427699d-01,
     : -4.80443069896770d-01, -3.33992566986309d-01,
     : 5.35323563771289d-01, 3.73473257874748d-08,
     : -5.35323621330841d-01, 7.57061878958283d-01,
     : -5.35323542381495d-01, -6.59317295605231d-08,
     : 5.35323599941047d-01, 7.21837810449948d-01,
     : -5.52470727926944d-01, 2.98994980253032d-01,
     : 4.20044534693153d-08, -2.98995049226790d-01,
     : 5.52470750505272d-01, -7.21837741476202d-01,
     : 4.03467612297724d-01, -4.03467608995464d-01,
     : 4.03467603165766d-01, -4.03467594181057d-01,
     : 4.03467580734490d-01, -4.03467559671056d-01,
     : 4.03467519591310d-01/

*  initialize all potentials to zero
      
      do  20 ipot=1, 3
        call dscal(16, zero, vpot(1,ipot), 1)
20    continue
      if (iselect .eq. 1 .or. iselect .eq. 4) then
        call pot_noag(z, pp, coefa)
        call dcopy(16, pp, 1, vpot(1,1), 1)
      
      end if  
      if (iselect .eq. 2 .or. iselect .eq. 4) then

        call pot_noag(z, pp, coefb)
        call dcopy(16, pp, 1, vpot(1,2), 1)

      endif
      if (iselect .eq. 3) then

        call pot_noag(z, pp, coefc)

        call dcopy(16, pp, 1, vpot(1,3), 1)
                          
      endif
* convert to vav and vdif

      do 30 i= 1,3

          vsum(1,i)=vpot(1,i)
          vsum(9,i)=vpot(9,i)
        do 30 j=1,7
          vdif(j,i)=half*(vpot(j+9,i)-vpot(j+1,i))
          vsum(j+1,i)=half*(vpot(j+9,i)+vpot(j+1,i))
30      continue
      
* convert potential to legendre coefficients

      if (iselect .eq. 1 .or. iselect .eq. 4) then
        call mxva(v0exp,9,1,vsum(1,1),1,vleg(1,1),1,9,9)
        call mxva(v2exp,7,1,vdif(1,1),1,vleg(10,1),1,7,7)
      end if  
      if (iselect .eq. 2 .or. iselect .eq. 4) then
        call mxva(v0exp,9,1,vsum(1,2),1,vleg(1,2),1,9,9)
        call mxva(v2exp,7,1,vdif(1,2),1,vleg(10,2),1,7,7)
      end if  
      if (iselect .eq. 3) then
        call mxva(v0exp,9,1,vsum(1,3),1,vleg(1,3),1,9,9)
        call mxva(v2exp,7,1,vdif(1,3),1,vleg(10,3),1,7,7)
      end if  

* damp out higher order anisotropies at large distances
      damp1=-half*(tanh(z-9.d0)-1.d0)
      damp2=-half*(tanh(z-7.d0)-1.d0)
      do 40 i=4,9
        if (i.le.5) then
          damp=damp1
        else
          damp=damp2
        endif
        call dscal(3,damp,vleg(i,1),16)
40    continue
      do 50 i=10,16
        if (i.le.11) then
          damp=damp1
        else
          damp=damp2
        endif
        call dscal(3,damp,vleg(i,1),16)
50    continue
                       
* store desired coefficients in vvl

      if (iselect .le. 3) then
* here if atop, bridge, or center site desired
        vv0=vleg(1,iselect)/econv
        do 60 i=2,16
          vvl(i-1)=vleg(i,iselect)/econv
60      continue
      else
* here if site-averaged potential desired
        vv0=quart*(vleg(1,1)+three*vleg(1,2))/econv
        do 70 i=2,16
          vvl(i-1)=quart*(vleg(i,1)+3*vleg(i,2))/econv
70      continue
      endif
      return
      end

                 
      subroutine pot_noag(z, vpot, c)

* to return fit to depristo no-ag(111) surface
* variables in call list

*   z     value of distance from surface
*   pot   on return: a vector of length 16 containing, respectively
*         vminus(0:22.5:180) then vplus(22.5:22.5:180), where an angle of 0 
*         corresponds to the o end down
*         form of fit is
*         c(5)*exp(-c(2)*z) + (c(6)+c(7)*z)exp(-c(3)*z)
*             -0.5[tanh{c(1)(z-c(4)))} +1]*c(8)/z^6
*   c     array of length 16*8 containing the coefficients used in the fit

      implicit double precision (a-h,o-z)
      dimension vpot(16), c(16,8)
      data half, zero, one /0.5d0, 0.d0, 1.d0/
      do 40 i = 1, 16
        vpot(i)=c(i,5)*exp(-c(i,2)*z)
     :         + (c(i,6)+c(i,7)*z)*exp(-c(i,3)*z)
     :         -half*(tanh(c(i,1)*(z-c(i,4))) +one)
     :               *c(i,8)/(z**6)
40    continue

      return

      end             
