* System:  NO(X 2Pi)+Ne, new avqz basis set CCSDT ab initio PES's
* Reference: m. h. alexander, p. soldan, t. g. wright, y. kim,
* h. meyer, p. j. dagdigian, and e. p. f. lee,  j. chem. phys. (in press)
* Vlam components correspond to theta=0 for NeNO


      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(15)
      common /coconv/ econv
      include "common/parpot"
      potnam='WRIGHT Ne-NO CCSDT'
      econv=219474.6d0
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
      include "common/parbas"
      include "common/parpot"
      potnam='WRIGHT Ne-NO CCSDT'
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
      subroutine pot (vv0, r)
*  subroutine to calculate the r-dependent coefficients in the
*  Ne-NO potentials of wright
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*  vv0      spherically symmetric term in V0
*  variable in common block /covvl/
*    vvl:     vector of length 15 to store r-dependence of each term
*             in potential expansion
*    vvl(1-8) expansion coefficients in dl0 (l=1:8) of vsum
*    vvl(9-15) expansion coefficients in dl2 (l=2:8) of vdif

* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  12-apr-2000
* NB Vlam components correspond to theta=0 for ArON!!!,
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      common /covvl/ vvl(15)
      common /coconv/ econv

      call dcopy(ldif,xdif,1,vvl(9),1)
* convert distance to angstroms
      rang=r*0.52917715d0
      do ilam=1,9
         if (ilam.eq.1) vv0=vlam0(ilam-1,rang)
         if (ilam.gt.1) then
             vvl(ilam-1)=vlam0(ilam-1,rang)
             if (ilam.gt.2) then
* change sign of difference potential so that it is defined at 0.5*(VA"-VA')
*              vvl(ilam+6)=-vlam2(ilam-1,rang)
* NO don't change it
               vvl(ilam+6)=vlam2(ilam-1,rang)
             endif
         endif
         x=2
      enddo
      end
*Date: Fri, 31 Mar 2000 10:23:53 +0100 (BST)
*From: Pavel Soldan <Pavel.Soldan@durham.ac.uk>
*To: Timothy Wright <T.G.Wright@sussex.ac.uk>
*Subject: NeNO
*
*And here are the modified subroutines. Regards, Pavel
       function vlam2(lam,rvdw)

c
c     an RKHS-AP subroutine for the NeNO Vdif potential
c
c     CCSD(T)/CAVQZ level of ab initio
c     written by Pavel Soldan, 27/03/2000
c     version 2 by PS 31/03/2000
c     Department of Chemistry, University of Durham
c     pavel.soldan@durham.ac.uk
c
c     coordinates: Jacobi with fixed diatom at re(NO) = 1.15258893 Ang
c                  theta = 0  for Ne.NO
c                  theta = pi for Ne.ON
c                  rvdw = R
c
c     units: R - Angstroem
c            E - Hartree
c
c     radial interpolation wrt R^2 using RKHS with m=2 and n=2
c     the function generates the angular projections V_lam2(R)
c     to d^lam_20, lam=2,3,...,8
c
c     WARNING! The user should check the consistency of definitions of
c     d^lam_20. The corresponding function rotmat(x,l)=d^l_20(x),
c     x=cos(theta), is attached at the end of the file.

       implicit none

       integer nlp, nr, k, lam
       parameter(nlp=7,nr=12)
       double precision vlam2, rvdw
       double precision rk
       double precision r(nr), al(nlp,nr), sum

        r( 1)= .23d+01
        r( 2)= .28d+01
        r( 3)= .31d+01
        r( 4)= .33d+01
        r( 5)= .35d+01
        r( 6)= .36d+01
        r( 7)= .37d+01
        r( 8)= .38d+01
        r( 9)= .40d+01
        r(10)= .45d+01
        r(11)= .50d+01
        r(12)= .70d+01
        al( 1, 1)=  -.1063895717637725d+02
        al( 1, 2)=   .1926318551397530d+02
        al( 1, 3)=  -.1326278593140001d+02
        al( 1, 4)=   .6663895152080782d+01
        al( 1, 5)=   .1558495990408475d+00
        al( 1, 6)=  -.3401960818800028d+01
        al( 1, 7)=   .4275159313334511d+01
        al( 1, 8)=  -.2167842339649448d+01
        al( 1, 9)=   .4243265847327822d+00
        al( 1,10)=  -.5467875935828594d+00
        al( 1,11)=  -.4683833808359152d+00
        al( 1,12)=  -.4115378983868183d+00
        al( 2, 1)=  -.7817991683399577d+01
        al( 2, 2)=   .1308953263065585d+02
        al( 2, 3)=  -.8626730888102136d+01
        al( 2, 4)=   .3983933831772284d+01
        al( 2, 5)=   .3543953190868268d+01
        al( 2, 6)=  -.9084406719423859d+01
        al( 2, 7)=   .9819286974632400d+01
        al( 2, 8)=  -.4632415315362341d+01
        al( 2, 9)=   .7332585713326709d+00
        al( 2,10)=  -.8994560785204848d-01
        al( 2,11)=  -.7275226148740468d+00
        al( 2,12)=  -.1961761720246335d+00
        al( 3, 1)=  -.1743432063872205d+02
        al( 3, 2)=   .3184768816412170d+02
        al( 3, 3)=  -.2281349620889748d+02
        al( 3, 4)=   .1317423300752003d+02
        al( 3, 5)=  -.5635213544736379d+01
        al( 3, 6)=   .3391964380020958d+01
        al( 3, 7)=  -.2981245939068647d+00
        al( 3, 8)=  -.3324115566485215d+00
        al( 3, 9)=  -.3565350520216878d-01
        al( 3,10)=  -.6621086248738322d+00
        al( 3,11)=  -.8475435209079638d+00
        al( 3,12)=  -.4942238890710146d+00
        al( 4, 1)=  -.4262898412923615d+01
        al( 4, 2)=   .8409058739956844d+01
        al( 4, 3)=  -.6220937260309124d+01
        al( 4, 4)=   .2757932254330065d+01
        al( 4, 5)=   .9479193299355574d+00
        al( 4, 6)=  -.3128931361853506d+01
        al( 4, 7)=   .2811968118014744d+01
        al( 4, 8)=  -.8013255076640041d+00
        al( 4, 9)=  -.3402523774771581d+00
        al( 4,10)=   .3523052180509561d+00
        al( 4,11)=  -.5469113163615668d+00
        al( 4,12)=   .1000663753319285d+00
        al( 5, 1)=  -.5211289278665402d+01
        al( 5, 2)=   .1117439844221021d+02
        al( 5, 3)=  -.8779428520631548d+01
        al( 5, 4)=   .4367888323544024d+01
        al( 5, 5)=  -.2347529670374828d+01
        al( 5, 6)=   .2639177059022312d+01
        al( 5, 7)=  -.2351299475349796d+01
        al( 5, 8)=   .1128980397957800d+01
        al( 5, 9)=  -.3359151801967284d+00
        al( 5,10)=  -.5159316047173260d-01
        al( 5,11)=  -.1749232191662927d+00
        al( 5,12)=  -.1288764968500672d+00
        al( 6, 1)=  -.1132932943276275d+01
        al( 6, 2)=   .2625342103182053d+01
        al( 6, 3)=  -.2234238148785230d+01
        al( 6, 4)=   .1142420007662185d+01
        al( 6, 5)=   .4499838071145829d-01
        al( 6, 6)=  -.1154144131378224d+01
        al( 6, 7)=   .9346566078146841d+00
        al( 6, 8)=   .8248363002545162d-01
        al( 6, 9)=  -.4063623223322678d+00
        al( 6,10)=   .2955454550138065d+00
        al( 6,11)=  -.2226690640949870d+00
        al( 6,12)=   .4692473466773672d-02
        al( 7, 1)=  -.1454500180611368d+01
        al( 7, 2)=   .3575946599691577d+01
        al( 7, 3)=  -.3272899190010239d+01
        al( 7, 4)=   .2051612123477050d+01
        al( 7, 5)=  -.2728745092971328d+01
        al( 7, 6)=   .3783299379174364d+01
        al( 7, 7)=  -.2150034447255586d+01
        al( 7, 8)=  -.6653785834382193d-01
        al( 7, 9)=   .3606216358716023d+00
        al( 7,10)=  -.6728749941300431d-01
        al( 7,11)=   .3355624794796382d-01
        al( 7,12)=  -.2595864233579441d+00

        sum=0.d0
        do k=1,nr
          sum=sum+al(lam-1,k)*rk(r(k),rvdw)
        enddo

        vlam2=sum

        return

        end



c      function rotmat(x,l)
c
c     implicit none
c     integer l
c     double precision rotmat,x
c     integer ll
c     double precision pll,pmm,pmmp1,fact
c
c     fact=dsqrt(1.d0/dble((l+2)*(l+1)*l*(l-1)))
c     pmm=3.d0*(1.d0-x**2)
c
c      if(l.eq.2) then
c        rotmat=fact*pmm
c        return
c      endif
c      pmmp1=5*x*pmm
c
c      if(l.eq.3) then
c        rotmat=fact*pmmp1
c        return
c      endif
c
c      do 12 ll=4,l
c        pll=(x*dble(2*ll-1)*pmmp1-dble(ll+1)*pmm)/dble(ll-2)
c        pmm=pmmp1
c        pmmp1=pll
c12    continue
c
c      rotmat=fact*pll
c
c      return
c      end
       function vlam0(lam,rvdw)

c
c     an RKHS-AP subroutine for the NeNO Vsum potential
c
c     CCSD(T)/CAVQZ level of ab initio
c     written by Pavel Soldan, 27/03/2000
c     version 2 by PS 31/03/2000
c     Department of Chemistry, University of Durham
c     pavel.soldan@durham.ac.uk
c
c     coordinates: Jacobi with fixed diatom at re(NO) = 1.15258893 Ang
c                  theta = 0  for Ne.NO
c                  theta = pi for Ne.ON
c                  rvdw = R
c
c     units: R - Angstroem
c            E - Hartree
c
c     radial interpolation wrt R^2 using RKHS with m=2 and n=2
c     the function generates the angular projection V_lam(R)
c     to P_lam, lam=0,1,...,8
c
c     WARNING! The user should check the consistency of definitions
c     of P_lam. The corresponding subroutine pleg(x,n)=P_n(x) is attached at
c     the end of this file.



       implicit none

       integer nlp, nr, k, lam
       parameter(nlp=9,nr=12)
       double precision vlam0, rvdw
       double precision rk
       double precision r(nr), al(nlp,nr), sum


        r( 1)= .23d+01
        r( 2)= .28d+01
        r( 3)= .31d+01
        r( 4)= .33d+01
        r( 5)= .35d+01
        r( 6)= .36d+01
        r( 7)= .37d+01
        r( 8)= .38d+01
        r( 9)= .40d+01
        r(10)= .45d+01
        r(11)= .50d+01
        r(12)= .70d+01
        al( 1, 1)=   .7903636094748597d+02
        al( 1, 2)=  -.1681262955615253d+03
        al( 1, 3)=   .1204207072936247d+03
        al( 1, 4)=  -.6465982121069091d+02
        al( 1, 5)=   .3290517294649015d+02
        al( 1, 6)=  -.1786457534340750d+02
        al( 1, 7)=   .7839582059649057d+01
        al( 1, 8)=  -.1243141994008235d+01
        al( 1, 9)=   .4639777592502777d+01
        al( 1,10)=   .3790009167979761d+01
        al( 1,11)=   .1731533675076544d+01
        al( 1,12)=   .3348108567212993d+00
        al( 2, 1)=   .2560124787505748d+02
        al( 2, 2)=  -.3958840221905596d+02
        al( 2, 3)=   .2008824423954857d+02
        al( 2, 4)=  -.1370670616515334d+02
        al( 2, 5)=   .2658787662620160d+01
        al( 2, 6)=   .1897164282868059d+01
        al( 2, 7)=  -.4181636676949754d+01
        al( 2, 8)=   .2465007037236447d+01
        al( 2, 9)=   .7803233753303294d+00
        al( 2,10)=   .2184025337736138d+01
        al( 2,11)=   .1459555803672008d+01
        al( 2,12)=   .3786058619457134d+00
        al( 3, 1)=   .1246024492976184d+03
        al( 3, 2)=  -.2601185738897672d+03
        al( 3, 3)=   .1870404398652330d+03
        al( 3, 4)=  -.9950575956446720d+02
        al( 3, 5)=   .4922402542296579d+02
        al( 3, 6)=  -.2339067867187471d+02
        al( 3, 7)=   .6635963378453869d+01
        al( 3, 8)=   .1468477042942264d+01
        al( 3, 9)=   .5708905434321164d+01
        al( 3,10)=   .5653669506857648d+01
        al( 3,11)=   .2273101778096229d+01
        al( 3,12)=   .4628978711135691d+00
        al( 4, 1)=   .2022647224512621d+02
        al( 4, 2)=  -.3721410788373363d+02
        al( 4, 3)=   .2315146142477892d+02
        al( 4, 4)=  -.1335137743000569d+02
        al( 4, 5)=   .6416858700215443d+01
        al( 4, 6)=  -.4770150986319783d+01
        al( 4, 7)=   .3879574257395865d+01
        al( 4, 8)=  -.1892953967347905d+01
        al( 4, 9)=   .1877193096548985d+01
        al( 4,10)=   .1095063971546565d+01
        al( 4,11)=   .6066502711582280d+00
        al( 4,12)=   .6036736631340286d-02
        al( 5, 1)=   .3426949621278042d+02
        al( 5, 2)=  -.7732838909525007d+02
        al( 5, 3)=   .5884476801582800d+02
        al( 5, 4)=  -.2992025583878842d+02
        al( 5, 5)=   .1652478497608870d+02
        al( 5, 6)=  -.9864685595233302d+01
        al( 5, 7)=   .6240414407056997d+01
        al( 5, 8)=  -.2079304053643933d+01
        al( 5, 9)=   .2294510873549914d+01
        al( 5,10)=   .9593106812238902d+00
        al( 5,11)=   .2442622741203566d+00
        al( 5,12)=  -.3339310718700594d+00
        al( 6, 1)=   .5979518866988759d+01
        al( 6, 2)=  -.1309442162723346d+02
        al( 6, 3)=   .9759669061575920d+01
        al( 6, 4)=  -.5150747939009410d+01
        al( 6, 5)=   .2888367224918531d+01
        al( 6, 6)=  -.2099694022165240d+01
        al( 6, 7)=   .2230822583389663d+01
        al( 6, 8)=  -.1440330053884334d+01
        al( 6, 9)=   .7483946801297372d+00
        al( 6,10)=   .8356497278152557d-01
        al( 6,11)=   .1585091200211814d+00
        al( 6,12)=  -.1669126784504235d+00
        al( 7, 1)=   .6197700122998614d+01
        al( 7, 2)=  -.1554948349896701d+02
        al( 7, 3)=   .1307363524313145d+02
        al( 7, 4)=  -.6214499642704915d+01
        al( 7, 5)=   .3486070131044158d+01
        al( 7, 6)=  -.3134450147862741d+01
        al( 7, 7)=   .3264214769381234d+01
        al( 7, 8)=  -.1647957866396965d+01
        al( 7, 9)=   .4345431416164996d+00
        al( 7,10)=   .1821392147796374d+00
        al( 7,11)=  -.9897198605301573d-01
        al( 7,12)=   .4037539445787783d-01
        al( 8, 1)=   .1899375359358360d+01
        al( 8, 2)=  -.4928039338087592d+01
        al( 8, 3)=   .4330654488070556d+01
        al( 8, 4)=  -.2059162122815359d+01
        al( 8, 5)=   .2085587886839157d+01
        al( 8, 6)=  -.4311678967363362d+01
        al( 8, 7)=   .5351510235828156d+01
        al( 8, 8)=  -.2821014653247558d+01
        al( 8, 9)=   .4714116107506663d+00
        al( 8,10)=  -.2271053444503596d-01
        al( 8,11)=  -.1397690087759670d-01
        al( 8,12)=   .4234882358971209d-01
        al( 9, 1)=   .2594516274949892d+01
        al( 9, 2)=  -.7044239890835624d+01
        al( 9, 3)=   .6534875899382561d+01
        al( 9, 4)=  -.3111914990212040d+01
        al( 9, 5)=   .1022350551567566d+01
        al( 9, 6)=  -.5054651232081595d+00
        al( 9, 7)=   .1447079919112991d+01
        al( 9, 8)=  -.1317718504807079d+01
        al( 9, 9)=   .4549160600041151d+00
        al( 9,10)=  -.9508422633464764d-01
        al( 9,11)=   .1888883695394340d-01
        al( 9,12)=   .7578194755683785d-02

        sum=0.d0
        do k=1,nr
          sum=sum+al(lam+1,k)*rk(r(k),rvdw)
        enddo

        vlam0=sum

        return

        end


         function rk(y,z)

         implicit none
         double precision rk,y,z
         double precision x_l,x_s,h

         x_l=z**2
         x_s=y**2

         if (x_l.lt.x_s) then
           x_s=z**2
           x_l=y**2
         endif

         h=x_s/x_l

         rk=1.d0*(1.d0-3.d0*h/5.d0)/3.d0/x_l**3

         return
         end


c      function pleg(x,n)
c
c      implicit none
c      integer j,n
c      double precision y(0:n+1),x,pleg
c
c      y(0)=1.d0
c      y(1)=x
c
c      if (n.ge.2) then
c        do j=1,n-1
c          y(j+1)=(dble(2*j+1)*x*y(j)-dble(j)*y(j-1))/dble(j+1)
c        enddo
c      endif
c
c      pleg=y(n)
c
c      return
c      end
