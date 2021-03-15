* System:  OH(X 2Pi)+Ar, UMP4 ab initio PES's
* Reference: 
* J. Klos, G. Chalasinski, M. T. Berry, R. A. Kendall, 
* R. Burcl, and M. M. Szczesniak, J. Chem. Phys. 112, 4952 (2000)


      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(13)
      common /coconv/ econv
      include "common/parpot"
      potnam='KLOS-CHALASINSKI Ar-OH UMP4'
      econv=219474.6d0
1      print *, ' r (bohr)'
      rshift=0.5
      xfact=0.8
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vsum',/,6(1pe16.8),/,
     :    '  vdif',/,8e16.8)
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
*  variable in common block /conlam/
*    nlam:      the number of angular coupling terms actually used
*    nlammx:    the maximum number of angular coupling terms
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx
      common /conlam/ nlam, nlammx, lamnum(2)
      potnam='KLOS-CHALASINSKI Ar-OH UMP4'
      lammin(1)=1
      lammax(1)=5
      lammin(2)=2
      lammax(2)=9
      lamnum(1)=lammax(1)-lammin(1)+1
      lamnum(2)=lammax(2)-lammin(2)+1
      nlam=lamnum(1)+lamnum(2)
      nlammx=nlam
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
*  Ar-OH potentials of klos and chalasinski
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*  vv0      spherically symmetric term in V0
*  variable in common block /covvl/
*    vvl:     vector of length 13 to store r-dependence of each term
*             in potential expansion
*    vvl(1-8) expansion coefficients in dl0 (l=1:8) of vsum
*    vvl(9-15) expansion coefficients in dl2 (l=2:8) of vdif

* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  11-mar-1999
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension xlam1(17),xlam2(17),vsum(10),xsum(10),vdif(10),xdif(10),
     :          ddif(10),
     :          d0(100),d2(64),aa(121)
      dimension kpvt(10),qraux(10),work(55),rsd(10), costh(10)

      common /covvl/ vvl(13)

      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
* coefficicients for d0 rotation matrices
* stored (by column) for each of 10 angles [0:20:180]
* and for l=0:9
      data d0/
     : 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0,
     : 9.3969262d-1, 7.6604445d-1, 5.0000001d-1, 1.7364819d-1,
     : -1.7364816d-1, -4.9999998d-1, -7.6604443d-1, -9.3969261d-1,
     : -1d0,1d0, 8.2453334d-1, 3.8023614d-1, -1.2499999d-1,
     : -4.5476946d-1, -4.5476947d-1, -1.2500002d-1, 3.8023610d-1,
     : 8.2453331d-1, 1d0,1d0, 6.6488474d-1, -2.5233323d-2,
     : -4.3750000d-1, -2.4738195d-1, 2.4738191d-1, 4.3750001d-1,
     : 2.5233373d-2, -6.6488469d-1, -1d0,1d0, 4.7497774d-1,
     : -3.1900434d-1, -2.8906251d-1, 2.6590160d-1, 2.6590163d-1,
     : -2.8906248d-1, -3.1900437d-1, 4.7497767d-1, 1d0,1d0,
     : 2.7149176d-1, -4.1968205d-1, 8.9843733d-2, 2.8101755d-1,
     : -2.8101752d-1, -8.9843784d-2, 4.1968205d-1, -2.7149167d-1,
     : -1d0,1d0, 7.1903012d-2, -3.2357074d-1, 3.2324218d-1,
     : -1.3212132d-1, -1.3212137d-1, 3.2324220d-1, -3.2357069d-1,
     : 7.1902917d-2, 1d0,1d0, -1.0722615d-1, -1.0060172d-1,
     : 2.2314455d-1, -2.8347993d-1, 2.8347991d-1, -2.2314450d-1,
     : 1.0060165d-1, 1.0722624d-1, -1d0,1d0, -2.5183942d-1,
     : 1.3862678d-1, -7.3638895d-2, 2.3307822d-2, 2.3307885d-2,
     : -7.3638959d-2, 1.3862685d-1, -2.5183950d-1, 1d0,1d0,
     : -3.5169654d-1, 2.9001294d-1, -2.6789855d-1, 2.5962717d-1,
     : -2.5962718d-1, 2.6789857d-1, -2.9001298d-1, 3.5169659d-1, -1d0/
* coefficicients for d2 rotation matrices
* stored (by column) for each of 8 angles [20:20:160]
* and for l=2:9
      data d2/
     : 7.1633967d-2, 2.5301754d-1, 4.5927933d-1, 5.9390715d-1,
     : 5.9390715d-1, 4.5927933d-1, 2.5301754d-1, 7.1633967d-2,
     : 1.5051848d-1, 4.3340069d-1, 5.1348990d-1, 2.3060769d-1,
     : -2.3060769d-1, -5.1348990d-1, -4.3340069d-1, -1.5051848d-1,
     : 2.3957418d-1, 5.0756736d-1, 2.2234765d-1, -3.0244624d-1,
     : -3.0244624d-1, 2.2234765d-1, 5.0756736d-1, 2.3957418d-1,
     : 3.2835759d-1, 4.3600553d-1, -1.6982082d-1, -2.7746877d-1,
     : 2.7746877d-1, 1.6982082d-1, -4.3600553d-1, -3.2835759d-1,
     : 4.0592179d-1, 2.3830028d-1, -3.4523418d-1, 1.5131756d-1,
     : 1.5131756d-1, -3.4523418d-1, 2.3830028d-1, 4.0592179d-1,
     : 4.6231022d-1, -1.3906540d-2, -1.9131358d-1, 2.8490318d-1,
     : -2.8490318d-1, 1.9131358d-1, 1.3906540d-2, -4.6231022d-1,
     : 4.8973053d-1, -2.2700359d-1, 1.1374299d-1, -3.5240961d-2,
     : -3.5240961d-2, 1.1374299d-1, -2.2700359d-1, 4.8973053d-1,
     : 4.8345441d-1, -3.2461587d-1, 2.7905801d-1, -2.6334950d-1,
     : 2.6334950d-1, -2.7905801d-1, 3.2461587d-1, -4.8345441d-1/
* cos of angle
      data costh /
     : 1d0, 9.3969262d-1, 7.6604444d-1, 5.d-1, 1.7364818d-1,
     : -1.7364818d-1, -5.d-1, -7.6604444d-1, -9.3969262d-1,-1d0/
* conversion from kcal/mol to hartree
      data econv /1.593563611216621d-3/


* determine A' and A" potentials at angles
* NOTE, that function vvdif works if angle is defined as pi-theta
* this checks with paper of klos et al
      rm3=one/r**3
      rm6=rm3*rm3
      do 100 i=1,10
*        vsum(i)=vvsum(r,costh(10-i+1))
        vsum(i)=vvsum(r,costh(i))
        if (i.ne.1 .or. i.ne.10) then
*          vdif(i)=vvdif(r,costh(i))
           vdif(i)=vvdif(r,costh(10-i+1))
        endif
100   continue
* diagonostic
*     do i=1,10
*        write (6,101) acos(costh(i))*180/3.1415927,
*    :             econv*1.e+6*(vsum(i)-vdif(i)),
*    :             econv*1.e+6*(vsum(i)+vdif(i))
*101      format(f5.1,2f8.1)
*     enddo
* solve simultaneous equations for solutions
      lsum=10
      ldif=8
      tol=1.d-10
      call dscal(13,zero,vvl,1)
      call dcopy(100,d0,1,aa,1)
      call dqrank(aa,10,10,lsum,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,10,10,lsum,kr,vsum,xsum,rsd,kpvt,qraux)
* remove terms less than 0.2 cm-1 in size
*     do 110 i=1,10
*       if (abs(xsum(i)) .lt. 0.2d0) xsum(i)=zero
*110   continue
      vv0=xsum(1)*econv
      call dcopy(5,xsum(2),1,vvl,1)
      call dcopy(64,d2,1,aa,1)
      call dqrank(aa,8,8,ldif,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,8,8,ldif,kr,vdif(2),xdif,rsd,kpvt,qraux)
*     do 120 i=1,8
*       if (abs(xdif(i)) .lt. 0.2d0) xdif(i)=zero
*120   continue
      call dcopy(ldif,xdif,1,vvl(6),1)
* convert potential to hartree
      call dscal(13,econv,vvl,1)
      end
      function vvsum(rr,y)
cCode for Vsum at MP4 level for  ArOH: Vsum=(A'+A")*1/2
c input: rr in a.u., y = cos(theta) theta=0 for Ar---HO
c output:interaction  energy in kcal/mol
c       ***********************************
c       *       JACEK KLOS                *
c       * Department of Chemistry         *
c       * Warsaw University               *
c       * Pasteura 1 Warsaw               *
c       * e-mail:jakl@tiger.chem.uw.edu.pl*
c       ***********************************
      implicit double precision (a-h,o-z)
      dimension d0(6),b0(6),g01(6),g02(6),g03(6),g04(6)
      dimension a0(10)
cPARAMETERS OF FIT
      data d0/14.3734874987198609d0,
     * 3.89630740092608452d0,
     * 1.55112894938056978d0,
     * 0.616459963478951933d0,
     * -0.534630164003549346d0,
     * -0.372685013274676014d0/
      data b0/1.22677731376697441d0,
     * -0.460276797308268315d-01,
     * -0.371931914390388665d-02,
     * -0.395580590551371306d-01,
     * -0.141398788676112303d0,
     * -0.704892824091889331d-01/
      data g01/0.332676554686949649d-01,
     * -0.923096247874362896d-01,
     *  0.685104110461919324d-01,
     * -0.319982385622900184d-01,
     * 0.116997210082613878d-01,
     * -0.261302012157420011d-02/
      data g02/-0.117243697167210160d-01,
     * 0.328498945044529700d-01,
     * -0.249797806757429339d-01,
     * 0.121336924362730227d-01,
     * -0.455258250360721702d-02,
     * 0.101416757246625047d-02/
      data g03/0.140732203577235634d-02,
     * -0.398594174329758907d-02,
     * 0.310302167970479686d-02,
     * -0.154389143153907719d-02,
     * 0.561374714121383130d-03,
     * -0.114377468377337450d-03/
      data g04/-0.602345630749664135d-04,
     * 0.172026827567311780d-03,
     * -0.136076045671759211d-03,
     * 0.678913697015956356d-04,
     * -0.221022353460932092d-04,
     * 0.338057804086432707d-05/
      data c60/-25450.5667812686115d0/
c Associated Legendre polynomial factors (calculated from Legendre polynomials
c     in this case)


      y2 = y*y
      y3 = y2*y
      y4 = y3*y
      y5 = y4*y
      y6 = y5*y
      y7 = y6*y
      y8 = y7*y
      y9 = y8*y
      y10 = y9*y
      a0(1) = 1.d0
      a0(2) = y
      a0(3) = 0.5d0*(3*y2-1)
      a0(4) = 0.5d0*(-3*y + 5*y3)
      a0(5) = 0.125d0*(35*y4-30*y2+3)
      a0(6) = 0.125d0*(15*y - 70*y3 + 63*y5)
      a0(7) = 0.0625d0*(231*y6-315*y4+105*y2-5)
      a0(8) = 0.0625d0*(429*y7-693*y5+315*y3-35*y)
      a0(9) = 0.0078125d0*(6435*y8-12012*y6+6930*y4-1260*y2+35)
      a0(10) = 0.0078125d0*(12155*y9-25740*y7+18018*y5-4620*y3+
     1 315*y)
      do i = 2,10
         fac = 2.d0*dfloat(i-1) + 1.d0
         a0(i) = a0(i)/dsqrt(fac)
      end do

c Vsh(R, theta) short range part:

      dth = 0.
      bth = 0.
      do i = 1,6
         dth = dth + d0(i)*a0(i)
         bth = bth + b0(i)*a0(i)
      end do

      r2 = rr*rr
      r3 = r2*rr
      grth = 0.
      do i = 1, 6
         grth = grth + (g01(i)+g02(i)*rr+g03(i)*r2+g04(i)*r3)*a0(i)
      end do
      vsh = grth*exp(dth - bth*rr)

c Vas(R, theta) long range part:

      bthr = bth*rr
      sum = 1.
      term = 1.
      do k = 1, 6
         term = term*bthr/k
         sum  = sum + term
      end do

      f6 = 1. - exp(-bthr)*sum

       bthr = bth*rr
      sum7 = 1.
      term7 = 1.
      do k = 1, 7
         term7 = term7*bthr/k
         sum7  = sum7 + term7
      end do


c     in case of f6 --> 0 use
c     d=1.0d0 - dexp(-br)*sum = sum_m=ncn+1^\infty br^m/m!
c (This section copied from slow version; see nhar12.f)

      if(abs(f6).lt.1.0e-8) then
         f6=0.
         do i = 7, 1000
            term=term*bthr/i
            f6 = f6 + term
            if((term/f6) .lt. 1.0d-8) go to 111
         end do
         write(12,*) 'No convergence in d'
 111     continue
         f6 = f6*exp(-bthr)
      endif

      r6 = r3*r3

      vas6 = (f6/r6)*(c60*a0(1))

       vvsum = vsh + vas6

       return
      end
cBELOW IS FUNCTION FOR vVdif (V_{2})
      function vvdif(rr,y)
c Code for Vdif at MP4 level for ArOH: Vdif= (A"-A')*1/2
c  y = cos(theta) theta=0 for Ar--HO
c units:
c input: rr in bohr, y = cos(theta) theta=0 for Ar--HO
c output in kcal/mol
c       ***********************************
c       *       JACEK KLOS                *
c       * Department of Chemistry         *
c       * Warsaw University               *
c       * Pasteura 1 Warsaw               *
c       * e-mail:jakl@tiger.chem.uw.edu.pl*
c       ***********************************
      implicit double precision (a-h,o-z)
      dimension d0(2),b0(2),g01(8),g02(8),g03(8),g04(8)
      dimension a0(10),a2(10)
      data d0/5.47501692013515484d0,
     * -2.25684943779012004d0/
      data b0/0.905145465036452990d0,
     * 0.570197159104062282d-01/
      data g01/1.96616606083388867d0,
     * 0.564599191229976705d0,
     * 0.136815038952975809d0,
     * 0.361642302149937853d-01,
     * 0.104998766848347462d-01,
     * 0.327212351055252131d-02,
     * -0.800107556727183223d-03,
     * -0.101415467561405165d-02/
      data g02/-0.816145146897791718d0,
     * -0.245918950088674093d0,
     * -0.592896168329055950d-01,
     * -0.146272786729374337d-01,
     * -0.395178028779387489d-02,
     * -0.121836640095142031d-02,
     * 0.338436660100453095d-03,
     * 0.409085314004364398d-03/
      data g03/0.104518822590993363d0,
     * 0.325376676369335677d-01,
     * 0.789019652727181338d-02,
     * 0.188888069834253691d-02,
     * 0.490244247245910517d-03,
     * 0.150552242193093601d-03,
     * -0.465675838715012402d-04,
     * -0.541733098832178591d-04/
      data g04/-0.491521516545659429d-02,
     * -0.157945560064353687d-02,
     * -0.385869053625031500d-03,
     * -0.885520741508413528d-04,
     * -0.214321944144399716d-04,
     * -0.630635761051296448d-05,
     * 0.209259707830972255d-05,
     * 0.236213077118763488d-05/
      data c62/13909.3430548761644d0/
      data c73/2956.56967099912572d0/
      data c84/488.446028923574204d0/
      data c95/-4269.83079783630910d0/
c Associated Legendre polynomial factors (calculated from Legendre polynomials
c     in this case)


      y2 = y*y
      y3 = y2*y
      y4 = y3*y
      y5 = y4*y
      y6 = y5*y
      y7 = y6*y
      y8 = y7*y
      y9 = y8*y
      y10 = y9*y
      y11 = y10*y
      a0(1) = 1.
      a0(2) = y
      a0(3) = 0.5d0*(3*y2-1)
      a0(4) = 0.5d0*(-3*y + 5*y3)
      a0(5) = 0.125d0*(35*y4-30*y2+3)
      a0(6) = 0.125d0*(15*y - 70*y3 + 63*y5)
      a0(7) = 0.0625d0*(231*y6-315*y4+105*y2-5)
      a0(8) = 0.0625d0*(429*y7-693*y5+315*y3-35*y)
      a0(9) = 0.0078125d0*(6435*y8-12012*y6+6930*y4-1260*y2+35)
      a0(10) = 0.0078125d0*(12155*y9-25740*y7+18018*y5-4620*y3+
     1 315*y)
      do i = 2, 10
         fac = 2.d0*dfloat(i-1) + 1.d0
         a0(i) = a0(i)/dsqrt(fac)
      end do
        a2(1) = 3.d0*(1-y2)
        a2(2) = 15.d0*(y-y3)
        a2(3) = 15.d0*(-1+8*y2-7*y4)/2.0
        a2(4) = 105.d0*(-y+4*y3-3*y5)/2.0
        a2(5) = 105.d0*(1-19*y2+51*y4-33*y6)/8.0
        a2(6) = 63.d0*(15*y-125*y3+253*y5-143*y7)/8.0
        a2(7) = 315.d0*(-1+34*y2-176*y4+286*y6-
     1 143*y8)/16.0
        a2(8)= 495.d0*(-7*y+98*y3-364*y5+494*y7-
     1 221*y9)/16.0
        a2(9) = 495.d0*(7-371*y2+3094*y4-8918*y6+10387*y8-
     1 4199*y10)/128.0
        a2(10) = 2145.d0*(21*y-441*y3+2562*y5-6018*y7+6137*y9-
     1 2261*y11)/128.0
c Vsh(R, theta) part:

      dth = 0.
      bth = 0.
      do i = 1, 2
         dth = dth + d0(i)*a0(i)
         bth = bth + b0(i)*a0(i)
      end do

      r2 = rr*rr
      r3 = r2*rr
      grth = 0.
      do i = 1, 8
      grth = grth + (g01(i)+g02(i)*rr+g03(i)*r2+g04(i)*r3)*a2(i)
      end do
      vsh = grth*dexp(dth - bth*rr)

c Vas(R, theta) part:

      r6 = r3*r3
      r7 = r6*rr
      r8 = r7*rr
      r9 = r6*r3

      vas6 = (1.d0/r6)*c62*a2(1)
      vas7 = (1.d0/r7)*c73*a2(2)
      vas8 = (1.d0/r8)*c84*a2(3)
      vas9 = (1.d0/r9)*c95*a2(4)

      vvdif = vsh+vas6+vas7+vas8+vas9

       return
      end
