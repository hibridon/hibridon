*reference:G. Capecchi, H.-J. Werner Phys. Chem. Chem. Phys. 6, 4975 (2004)
*system: Cl-H2 using CW-ADiabatic
*
      include "common/syusr"
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(2)
      include "common/parpot"
      potnam='CAPECCHI-WERNER CLH2 ADIABATIC POT # 2' 
      print *, potnam 
1     print *, ' r (bohr)'
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vlambda',/,7(1pe16.8))
      goto 1
99    end
      include "common/bausr"
      include "common/ground"
* --------------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parpot"
      include "common/parbas"
      potnam='CAPECCHI-WERNER CLH2 ADIABATIC POT # 2' 
      lammin(1)=2
      lammax(1)=4
      mproj(1)=0
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end

      subroutine pot (vv0, r)
*  -----------------------------------------------------------------------

*  subroutine to calculate the r-dependent coefficients in the
*  collision of a homonuclear diatomic with a structureless target
*  in units of hartree for energy and bohr for distance

*  on return:
*  vv0 contains the isotropic term (n=0) in the potential
*  the coefficients for each angular term in the coupling potential
*  [ vvl(i) for i = 1, nlam ] are returned in common block /covvl/
*  vvl(1:2) contains the anisotropic (n=2,4) terms in the potential

*  variable in common block /conlam/
*    nlammx:    the maximum number of anisotropic terms
*    nlam:      the total number of angular coupling terms
*  variable in common block /covvl/
*    vvl        array to store r-dependence of each angular term in the
*               potential

*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      common /conlam/ nlam,nlammx

      include "common/parbas"
      common /covvl/ vvl(2)

      dimension rr(31),c10(31),c12(31),c14(31),c20(31),c22(31),c24(31),
     :          c30(31),c32(31),c34(31),xmat(31,3),cspline(31,3),
     :          vlam(3)
      data rr/
     : 4d0, 4.25d0, 4.5d0, 4.75d0, 5d0, 5.25d0, 5.5d0, 5.75d0, 6d0,
     : 6.25d0, 6.5d0, 6.75d0, 7d0, 7.25d0, 7.5d0, 7.75d0, 8d0, 8.5d0,
     : 9d0, 9.5d0, 1d1, 1.05d1, 1.1d1, 1.15d1, 1.2d1, 1.4d1, 1.6d1,
     : 1.8d1, 2d1, 2.2d1, 2.4d1/
      data c10/
     : 1.2880353d3, 7.4434706d2, 3.6305993d2, 9.8221511d1, -8.2960212d1,
     : -2.0496318d2, -2.8512954d2, -3.3487224d2, -3.6269188d2,
     : -3.7584040d2, -3.7943740d2, -3.7682757d2, -3.7065354d2,
     : -3.6277668d2, -3.5439867d2, -3.4625110d2, -3.3873865d2,
     : -3.2625291d2, -3.1717209d2, -3.1091109d2, -3.0673258d2,
     : -3.0401069d2, -3.0227805d2, -3.0121857d2, -3.0064037d2,
     : -2.9945716d2, -2.9620100d2, -2.9437477d2, -2.9419992d2,
     : -2.9414053d2, -2.9411898d2/
      data c12/
     : 5.3934429d2, 4.9066272d2, 4.0658021d2, 3.1660624d2, 2.3442022d2,
     : 1.6611400d2, 1.1309135d2, 7.4123885d1, 4.6340554d1, 2.3825721d1,
     : 6.2915360d0, -4.6265952d0, -1.0686955d1, -1.3500996d1,
     : -1.4294577d1, -1.3905791d1, -1.2880154d1, -1.0194586d1,
     : -7.6060201d0, -5.4537966d0, -3.7415335d0, -2.3708553d0,
     : -1.2768711d0, -4.3994296d-1, -1.0063081d-1, -4.2458530d-1,
     : -8.3210634d-2, 1.1573419d-1, 5.8679008d-2, 2.5562538d-2,
     : 3.5140445d-3/
      data c14/
     : 2.8458721d1, 2.3320584d1, 2.6086159d1, 2.5164172d1, 2.1667791d1,
     : 1.7214985d1, 1.3912817d1, 1.2214307d1, 1.1592148d1, 5.9230729d0,
     : -6.5411340d-1, -3.1039665d0, -3.4699424d0, -3.0261285d0,
     : -2.3944105d0, -1.8287867d0, -1.4051834d0, -9.4716734d-1,
     : -7.6745680d-1, -6.7732841d-1, -5.7217863d-1, -4.1267871d-1,
     : -1.6578158d-1, 1.8653842d-1, 1.7370196d-1, -2.7202381d-1,
     : -1.8817712d-1, 7.4861171d-1, 6.5711776d-1, 1.9172763d-1,
     : 3.4139535d-2/
      data c20/
     : 7.0761536d3, 4.4759085d3, 2.7371484d3, 1.5904665d3, 8.4536163d2,
     : 3.6902114d2, 7.02748d1, -1.1280513d2, -2.2160411d2, -2.8303573d2,
     : -3.1514340d2, -3.2998089d2, -3.3487878d2, -3.3431901d2,
     : -3.3100995d2, -3.2655361d2, -3.2187072d2, -3.1355097d2,
     : -3.0748342d2, -3.0350199d2, -3.0102426d2, -2.9950349d2,
     : -2.9852361d2, -2.9781891d2, -2.9720121d2, -2.9531040d2,
     : -2.9453073d2, -2.9384191d2, -2.9360762d2, -2.9395491d2,
     : -2.9408872d2/
      data c22/
     : 6.8345084d2, 3.9740870d2, 2.1084908d2, 9.281171d1, 2.154166d1,
     : -1.9138261d1, -4.0211225d1, -4.8814728d1, -4.9435542d1,
     : -4.23267d1,-3.1886893d1,-2.300382d1, -1.6089502d1, -1.0902321d1,
     : -7.0929328d0, -4.3409268d0, -2.3787983d0, -1.2835370d-1,
     : 7.1642992d-1, 7.4309475d-1, 3.2435477d-1, -3.0742363d-1,
     : -9.9843313d-1, -1.5962818d0, -1.7675166d0, -8.0411098d-1,
     : -3.0930376d-1, 2.4709062d-1, 4.8111619d-1, 1.6206025d-1,
     : 2.9849248d-2/
      data c24/
     : 4.5402856d1,1.7544993d1,3.6652224d0,-2.8458707d0,-5.9075048d0,
     : -7.3012176d0, -8.3719701d0, -9.4342559d0, -1.0318606d1,
     : -5.3203054d0, 1.0741660d0, 3.5920251d0, 4.1063031d0, 3.7886689d0,
     : 3.2331504d0, 2.6708040d0, 2.1756442d0, 1.4488819d0, 9.3105340d-1,
     : 5.2411015d-1, 1.6823964d-1, -1.6250418d-1, -5.0226084d-1,
     : -8.9145154d-1, -8.5430773d-1, -1.3655977d-3, 5.2520259d-1,
     : 5.3530789d-2, 2.1303068d-2, 1.9391255d-2, 7.9203720d-3/
      data c30/
     : 7.7367199d3, 5.1252496d3, 3.3940126d3, 2.2676404d3, 1.5518785d3,
     : 1.1098330d3, 8.4526187d2, 6.9178402d2, 6.0584035d2, 5.6019536d2,
     : 5.3827226d2, 5.3007197d2, 5.2957865d2, 5.3316397d2, 5.3867254d2,
     : 5.4483041d2, 5.5091754d2, 5.6160650d2, 5.6974323d2, 5.7561215d2,
     : 5.7979195d2, 5.8282305d2, 5.8510839d2, 5.8689472d2, 5.8832931d2,
     : 5.9133659d2, 5.9013357d2, 5.8891081d2, 5.8855328d2, 5.8831407d2,
     : 5.8824371d2/
      data c32/
     : 6.2667222d2, 3.7086655d2, 2.0140967d2, 9.578425d1, 3.5720381d1,
     : 5.3492307d0, -8.2769110d0, -1.3808292d1, -1.5914712d1,
     : -1.6455459d1, -1.6080556d1, -1.5092878d1, -1.3685133d1,
     : -1.2068853d1, -1.0394458d1, -8.7792205d0, -7.3074156d0,
     : -4.9163317d0, -3.2770495d0, -2.2447142d0, -1.6377902d0,
     : -1.2934064d0, -1.1050135d0, -9.9678404d-1, -9.1938023d-1,
     : -6.1155525d-1, -1.9331977d-1, 1.8413163d-1, 2.7035279d-1,
     : 9.8311258d-2, 1.7023041d-2/
      data c34/
     : 4.9592459d1, 1.9637064d1, 5.0717711d0, -9.6409958d-1,
     : -2.2777894d0, -1.3819683d0, -2.3761123d-1, 4.2223766d-1,
     : 6.8224434d-1, 7.5871689d-1, 7.6226740d-1, 7.7128024d-1,
     : 7.6445234d-1, 7.5188900d-1, 7.1938818d-1, 6.5684456d-1,
     : 5.7600216d-1, 3.7198433d-1, 1.5321821d-1, -3.9601875d-2,
     : -1.8353415d-1, -2.8049042d-1, -3.3484056d-1, -3.5505116d-1,
     : -3.4385340d-1, -1.3737749d-1, 1.6468911d-1, 3.9820353d-1,
     : 3.3811795d-1, 1.0077991d-1, 2.0210604d-2/
      data one,half,alph /1d0,0.5d0,1.2d0/
      data cmtohar / 219474.6d0/
      data ifirst /0/
      ipot =2
* ipot = 1 for lowest adiabatic pes
* ipot = 2 for next adiabatic pes
* ipot = 3 for highest adiabatic pes

*  determine spline coefficients on first pass
      if (ifirst.eq.0) then
         if (ipot .eq. 1) then
             call dcopy(31,c10,1,xmat(1,1),1)
             call dcopy(31,c12,1,xmat(1,2),1)
             call dcopy(31,c14,1,xmat(1,3),1)
         elseif (ipot.eq.2) then
             call dcopy(31,c20,1,xmat(1,1),1)
             call dcopy(31,c22,1,xmat(1,2),1)
             call dcopy(31,c24,1,xmat(1,3),1)
         else if (ipot.eq.3) then
             call dcopy(31,c30,1,xmat(1,1),1)
             call dcopy(31,c32,1,xmat(1,2),1)
             call dcopy(31,c34,1,xmat(1,3),1)
         endif
         do i=1,3
*  estimate first derivatives
            d1=(xmat(2,i)-xmat(1,i))/(rr(2)-rr(1))
            dn=(xmat(31,i)-xmat(30,i))/(rr(31)-rr(30))
* determine spline coefficient
            call dspline(rr,xmat(1,i),31,d1,dn,cspline(1,i))
         enddo
         ifirst=1
      endif
      do i=1,3
        vlam(i)=0d0
      enddo
* spline interpolation      
      do i=1,3
          call dsplint(rr,xmat(1,i),cspline(1,i),31,r,vlam(i))
      enddo
* check to see if r is outside range of spline knots
      if (r.lt.rr(1)) then
         do i=2,3
            vlam(i)=xmat(1,i)
         enddo
* exponential extrapolation for spherically symmetric term
         bexp=log(xmat(1,1)/xmat(2,1))/(rr(2)-rr(1))
         a=xmat(1,1)*exp(bexp*rr(1))
         vlam(1)=a*exp(-bexp*r)
      elseif(r.gt.rr(31)) then
         do i=2,3
            vlam(i)=0d0
         enddo
         vlam(1)=xmat(31,1)
      endif
* convert to hartree
      do i=1,3
         vlam(i)=vlam(i)/cmtohar
      enddo
      vv0=vlam(1)
      do i=2,3
         vvl(i-1)=vlam(i)
      enddo
      end
