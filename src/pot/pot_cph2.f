*system:  C+(2P)+H2, Dubernet-Hutson expansion of Halvick PES's
*authors francois Lique 21/02/2013 
*
*  F. Lique, G. Werfelli, P. Halvick, T. Stoecklin, A. Faure, L. Wiesenfeld,
*  and P. J. Dagdigian, J. Chem. Phys. 138, 204314 (2013).
*
      include "common/syusr"
      include "common/ground"
      include "common/bausr"
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      common /coselb/ ibasty
      potnam='HALVICK Cp(2P)-H2 PES'
      ibasty=12
      lammin(1)=1
      lammax(1)=4
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
      potnam='HALVICK C+(2P)-H2 PES'
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
93    r=1.
      do i=1,500
       call pot(vv0,r)
       write(2,101) r,vv0,vvl
101    format(f8.4,10(1pe16.8))
       r=r+0.2
      enddo

99    end
      subroutine pot (vv0, r)
* ----------------------------------------------------------------------
*  subroutine to calculate the r-dependent coefficients in the
*  al-h2 potentials of williams and alexander using body-frame
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
*  CC calculations (csflag=.false.)
*    vvl(1)  V220(R)
*    vvl(2)  V022(R)
*    vvl(3)  V202(R)
*    vvl(4)  V222(R)
* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  1-feb-1998
* ----------------------------------------------------------------------


      implicit double precision (a-h,o-z)
      logical ljunk, ihomo, csflag, lljunk
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(4)
      dimension evec(5)
      dimension vvll(18)
      data zero, one, half /0.d0,1.d0,0.5d0/
      data two, three, four, five /2.d0,3.d0,4.d0,5.d0/
      data seven,eleven /7.d0,11.d0/
      data sq2 /1.414213562d0/

      evec(1)=  fpot_per1(r)
      evec(2)=  fpot_per2(r)
      evec(3)=  fpot_per3(r)
      evec(4) = fpot_lin1(r)
      evec(5) = fpot_lin2(r)
      rac5 = sqrt(5.d0)

* here is totally symmetric term
      vv0=(2.d0*(evec(1)+evec(2)+evec(3)+evec(5))+evec(4))/9.d0
* v200
      vvl(1)=(2.d0/(9.d0*rac5))*(-evec(1)-evec(2)-evec(3)+evec(4)
     *       +2.d0*evec(5))
* v020
      vvl(2)=(rac5/9.d0)*(2.d0*evec(1)-evec(2)-evec(3)
     *+evec(4)-evec(5))
* v220
      vvl(3)=(2.d0*(evec(4)-evec(1)-evec(5))+evec(2)+evec(3))/9.d0
* v222
      vvl(4)=(evec(3)-evec(2))/3.d0
* transform to space frame terms if CC or CS case 1C
      if (.not. csflag .or. (csflag .and. ihomo)) then
        call dcopy(4,vvl,1,vvll,1)
        vvl(1)=5d0*(two*(vvll(4))+vvll(3))/sqrt(5.d0)
c        vvl(1)=5d0*(vvll(4)+vvll(3))/sqrt(5.d0)
        vvl(2)=dsqrt(5d0)*vvll(2)
        vvl(3)=dsqrt(5d0)*vvll(1)
        vvl(4)=5d0*(two*vvll(4)-vvll(3))*sqrt(two/7.d0)
c        vvl(4)=5d0*(vvll(4)-vvll(3))*sqrt(two/7.d0)
      else
* reorder body-frame terms
        vvll(1)=vvl(1)
        vvl(1)=vvl(3)
        vvl(3)=vvll(1)
      endif
* convert to hartree
      econv=1./219474.6
      call dscal(4,econv,vvl,1)
      vv0=vv0*econv
c      write(*,*) vv0,vvl(1),vvl(2),vvl(3),vvl(4)

      end

c---------------------------------------------------------
      double precision function fpol20(x)
c---------------------------------------------------------
      real(8) :: x
      fpol20=0.5d0*(3.d0*x**2-1.d0)
      return
      end
c---------------------------------------------------------
      double precision function fpol22(x)
c---------------------------------------------------------
      real(8) :: x
      fpol22=3.d0*(1.d0-x**2)
      return
      end
c-----------------------------------------------------------------------
      double precision function fpot_lin1(y)
c-----------------------------------------------------------------------
c
c     potentiel RKHS pour l'etat Sigma de H2 + C+ en gÈomÈtrie linÈaire
c     distance en bohr, Ènergie en cm-1
c     18/02/13 nouveau fit.
c
      implicit double precision(a-h,o-z)
      integer, parameter :: nx=44,mq2m=2
      dimension xgr(nx),coef(nx)
c
      data xgr/ 40.00d0,35.00d0,30.00d0,28.00d0,26.00d0,24.00d0,
     >          22.00d0,20.00d0,16.00d0,14.00d0,12.00d0,11.00d0,
     >          10.00d0, 9.00d0, 8.00d0, 7.50d0, 7.00d0, 6.50d0,
     >           6.00d0, 5.60d0, 5.20d0, 5.00d0, 4.80d0, 4.60d0,
     >           4.40d0, 4.20d0, 4.00d0, 3.80d0, 3.60d0, 3.40d0,
     >           3.30d0, 3.20d0, 3.10d0, 3.00d0, 2.90d0, 2.80d0,
     >           2.70d0, 2.60d0, 2.50d0, 2.40d0, 2.30d0, 2.20d0,
     >           2.10d0, 2.00d0/
       data coef/ 0.491913754828200D+06, -0.389614553294182D+06,
     >            0.105477283936191D+06, -0.710176850920096D+05,
     >            0.343896718753763D+05, -0.589590962989149D+05,
     >            0.291534767976850D+05, -0.894194817649573D+04,
     >            0.496735828223487D+05, -0.128072282970141D+06,
     >           -0.649536954822289D+05,  0.113064790922805D+06,
     >           -0.137219507293418D+06,  0.475462665887163D+06,
     >            0.107027739037149D+07,  0.737243033847262D+06,
     >            0.789969183787059D+06,  0.122419871220695D+07,
     >            0.720862074381618D+06, -0.760041942592218D+06,
     >           -0.121326689221023D+07, -0.805685614645418D+06,
     >           -0.109601364248240D+07, -0.938972562577993D+06,
     >           -0.680704186234914D+06, -0.357739765567921D+06,
     >            0.138219294092268D+06,  0.978656518573068D+06,
     >            0.190668686731875D+07,  0.161950027155279D+07,
     >            0.676512333576357D+06,  0.864117788340936D+06,
     >            0.618104120347438D+06,  0.425780308388997D+06,
     >            0.192075163116185D+06,  0.663079310771450D+05,
     >           -0.344729798261138D+06,  0.214460849006149D+06,
     >           -0.242630687156751D+07,  0.563631264524441D+07,
     >           -0.213864639054639D+08,  0.669736784265141D+08,
     >           -0.219633576376215D+09,  0.164839840328613D+09/
c
      fpot_lin1 = frecip(y,mq2m,nx,xgr,coef)
c
      return
      end
c-----------------------------------------------------------------------
      double precision function fpot_lin2(y)
c-----------------------------------------------------------------------
c
c     potentiel RKHS pour l'Ètat Pi de H2 + C+ en gÈometrie LinÈaire
c     distance en bohr, Ènergie en cm-1
c     18/02/13 nouveau fit.
c
      implicit double precision(a-h,o-z)
      integer, parameter :: nx=42,mq2m=2
      dimension xgr(nx),coef(nx)
c
      data xgr/ 40.00d0,35.00d0,30.00d0,28.00d0,26.00d0,24.00d0,
     >          22.00d0,20.00d0,16.00d0,14.00d0,12.00d0,11.00d0,
     >          10.00d0, 9.00d0, 8.00d0, 7.50d0, 7.00d0, 6.50d0,
     >           6.00d0, 5.60d0, 5.20d0, 5.00d0, 4.80d0, 4.60d0,
     >           4.40d0, 4.20d0, 4.00d0, 3.80d0, 3.60d0, 3.40d0,
     >           3.30d0, 3.20d0, 3.10d0, 3.00d0, 2.90d0, 2.80d0,
     >           2.70d0, 2.60d0, 2.50d0, 2.40d0, 2.30d0, 2.20d0/
      data coef/ -0.102423380113364D+07,  0.800281550858557D+06,
     >           -0.210902191201732D+06,  0.863368091983516D+05,
     >           -0.591350328831654D+04, -0.228764266510822D+05,
     >           -0.472812773696892D+04, -0.746357733673602D+04,
     >           -0.384465890498897D+05,  0.549937448861031D+04,
     >           -0.202074486933749D+06, -0.228440181518286D+05,
     >           -0.342270644904602D+06, -0.416877900891068D+06,
     >           -0.358906670194745D+06, -0.739937246622927D+05,
     >           -0.109617860144651D+06,  0.486711399498861D+06,
     >            0.108820759246288D+07,  0.178306653538519D+07,
     >            0.174507428033313D+07,  0.100859271143262D+07,
     >            0.138625953888559D+07,  0.129401091884049D+07,
     >            0.105678884309736D+07,  0.537594868398659D+06,
     >            0.969042512216419D+04, -0.236024577648222D+06,
     >           -0.400352190889070D+06, -0.368112191223879D+06,
     >           -0.212698792762753D+06, -0.270190460321086D+06,
     >           -0.302761247769324D+06, -0.262214721955683D+06,
     >           -0.476234063878375D+06,  0.145287176651437D+06,
     >           -0.209022132005335D+07,  0.511935958585044D+07,
     >           -0.191701353821519D+08,  0.608200376268663D+08,
     >           -0.202729474816972D+09,  0.151955318041259D+09/
c
       fpot_lin2 = frecip(y,mq2m,nx,xgr,coef)
c
      return
      end
c-----------------------------------------------------------------------
      double precision function fpot_per1(y)
c-----------------------------------------------------------------------
c
c     potentiel RKHS pour l'Ètat A1 de H2 + C+ en gÈometrie perpendiculaire
c     distance en bohr, Ènergie en cm-1
c     18/02/13 nouveau fit.
c
      implicit double precision(a-h,o-z)
      integer, parameter :: nx=38,mq2m=2
      dimension xgr(nx),coef(nx)
c
      data xgr/ 40.00d0,35.00d0,30.00d0,28.00d0,26.00d0,24.00d0,
     >          22.00d0,20.00d0,16.00d0,14.00d0,12.00d0,11.00d0,
     >          10.00d0, 9.00d0, 8.00d0, 7.50d0, 7.00d0, 6.50d0,
     >           6.00d0, 5.60d0, 5.20d0, 5.00d0, 4.80d0, 4.60d0,
     >           4.40d0, 4.20d0, 4.00d0, 3.80d0, 3.60d0, 3.40d0,
     >           3.30d0, 3.20d0, 3.10d0, 3.00d0, 2.90d0, 2.80d0,
     >           2.70d0, 2.60d0/
      data coef/  0.497312498192842D+06, -0.389362449521100D+06,
     >            0.103527946816084D+06, -0.487796806147103D+05,
     >            0.948099111235072D+04, -0.323846887910337D+04,
     >           -0.119593097424195D+05, -0.462132253629155D+03,
     >            0.173214954850661D+05,  0.715301757539238D+04,
     >           -0.551887732025340D+05, -0.683460760152683D+05,
     >            0.111277239582690D+06,  0.166811802097189D+06,
     >            0.912018066104369D+06,  0.248396644216687D+06,
     >            0.124156166502989D+07,  0.849323184364350D+06,
     >            0.633960584960461D+06,  0.405718187857017D+06,
     >            0.270418216828844D+05, -0.245877475104116D+06,
     >           -0.391139137461271D+06, -0.540462910746460D+06,
     >           -0.677405964641019D+06, -0.794357681383695D+06,
     >           -0.950939768885401D+06, -0.101399178450056D+07,
     >           -0.124621897509947D+07, -0.850889109742134D+06,
     >           -0.857459451295006D+06,  0.152146669844492D+07,
     >           -0.521701400048766D+07,  0.231422493371648D+07,
     >           -0.193193542712688D+07, -0.628330518445135D+07,
     >           -0.359474777877695D+08,  0.484690367485688D+08/
c
      fpot_per1 = frecip(y,mq2m,nx,xgr,coef)
c
      return
      end
c-----------------------------------------------------------------------
      double precision function fpot_per2(y)
c-----------------------------------------------------------------------
c
c     potentiel RKHS pour l'etat B1 de H2 + C+ en geometrie perpendiculaire
c     distance en bohr, energie en cm-1
c     18/02/13 nouveau fit.
c
      implicit double precision(a-h,o-z)
      integer, parameter :: nx=46,mq2m=2
      dimension xgr(nx),coef(nx)
c
      data xgr/   40.00d0,35.00d0,30.00d0,28.00d0,26.00d0,24.00d0,
     >           22.00d0,20.00d0,16.00d0,14.00d0,12.00d0,11.00d0,
     >           10.00d0, 9.00d0, 8.00d0, 7.50d0, 7.00d0, 6.50d0,
     >           6.00d0, 5.60d0, 5.20d0, 5.00d0, 4.80d0, 4.60d0,
     >           4.40d0, 4.20d0, 4.00d0, 3.80d0, 3.60d0, 3.40d0,
     >           3.30d0, 3.20d0, 3.10d0, 3.00d0, 2.90d0, 2.80d0,
     >           2.70d0, 2.60d0, 2.50d0, 2.40d0, 2.30d0, 2.20d0,
     >           2.10d0, 2.00d0, 1.90d0, 1.80d0/
      data coef/  0.492828716719614D+06, -0.388389788109108D+06,
     >            0.922857767024105D+05, -0.151407584199640D+05,
     >           -0.987405714565953D+05,  0.309575955309148D+06,
     >           -0.916834130193340D+06,  0.135897518996735D+07,
     >           -0.447312016654355D+07,  0.855990586794297D+07,
     >           -0.958094884891660D+07,  0.553905178308592D+07,
     >           -0.124892426463643D+07, -0.401101886716017D+06,
     >            0.397907664877789D+06, -0.625063970834337D+06,
     >            0.422460583100225D+06,  0.439115653485784D+06,
     >            0.671367515088769D+06,  0.133514005011755D+07,
     >            0.101042025279355D+07,  0.706268730222982D+06,
     >            0.736190600392109D+06,  0.599533462817998D+06,
     >            0.397832851620759D+06,  0.953969101499394D+05,
     >           -0.271390055025436D+06, -0.626772742787768D+06,
     >           -0.100559680155128D+07, -0.921325599133609D+06,
     >           -0.574766572889371D+06, -0.762471715664800D+06,
     >           -0.534991366598009D+06, -0.370684167599710D+06,
     >           -0.298157237760607D+06, -0.232247607609851D+06,
     >           -0.167696163083141D+06, -0.958041555478132D+05,
     >           -0.898316897797301D+05,  0.162004099637610D+05,
     >           -0.332661050583040D+06,  0.647610928779313D+06,
     >           -0.285000603524927D+07,  0.808473509656521D+07,
     >           -0.273876072213895D+08,  0.223685816358945D+08/
c
      fpot_per2 = frecip(y,mq2m,nx,xgr,coef)
c
      return
      end
c-----------------------------------------------------------------------
      double precision function fpot_per3(y)
c-----------------------------------------------------------------------
c
c     potentiel RKHS pour l'etat B2 de H2 + C+ en geometrie perpendiculaire
c     distance en bohr, energie en cm-1
c     18/02/13 nouveau fit.
c
      implicit double precision(a-h,o-z)
      integer, parameter :: nx=49,mq2m=2
      dimension xgr(nx),coef(nx)
c
      data xgr/ 40.00d0,35.00d0,30.00d0,28.00d0,26.00d0,24.00d0,
     >          22.00d0,20.00d0,16.00d0,14.00d0,12.00d0,11.00d0,
     >          10.00d0, 9.00d0, 8.00d0, 7.50d0, 7.00d0, 6.50d0,
     >           6.00d0, 5.60d0, 5.20d0, 5.00d0, 4.80d0, 4.60d0,
     >           4.40d0, 4.20d0, 4.00d0, 3.80d0, 3.60d0, 3.40d0,
     >           3.30d0, 3.20d0, 3.10d0, 3.00d0, 2.90d0, 2.80d0,
     >           2.70d0, 2.60d0, 2.50d0, 2.40d0, 2.30d0, 2.20d0,
     >           2.10d0, 2.00d0, 1.90d0, 1.80d0, 1.70d0, 1.60d0,
     >           1.50d0/
      data coef/  0.492003994527753D+06, -0.388577738953952D+06,
     >            0.945873730147467D+05, -0.261159300683036D+05,
     >           -0.656199845401742D+05,  0.210629055998609D+06,
     >           -0.634113202008125D+06,  0.928955241437864D+06,
     >           -0.308328715436063D+07,  0.576984818169066D+07,
     >           -0.624019322117409D+07,  0.353794552179810D+07,
     >           -0.990133508000862D+06, -0.354656076639005D+06,
     >            0.306173800930716D+06, -0.584050303537930D+06,
     >            0.424344787881094D+06,  0.349213312823969D+06,
     >            0.726029556404268D+06,  0.108508541799525D+07,
     >            0.114103230787236D+07,  0.573499431998308D+06,
     >            0.696911256176775D+06,  0.543284714920591D+06,
     >            0.327504846587345D+06,  0.360353708783535D+05,
     >           -0.324452632926787D+06, -0.664048177438341D+06,
     >           -0.953945125796022D+06, -0.839245773077455D+06,
     >           -0.530709598515846D+06, -0.570091145896127D+06,
     >           -0.306411409557792D+06, -0.182144403914101D+06,
     >           -0.812582772255084D+05,  0.149579782034294D+05,
     >            0.946326535285505D+05,  0.135472425050693D+06,
     >            0.119689417229309D+06,  0.427923179727981D+05,
     >           -0.691132014144378D+05, -0.227080768508594D+06,
     >           -0.271443559997739D+06, -0.592772701292631D+06,
     >            0.159779879777001D+06, -0.229182622438273D+07,
     >            0.559574463188051D+07, -0.187156775390270D+08,
     >            0.155922298373337D+08/
c
      fpot_per3 = frecip(y,mq2m,nx,xgr,coef)
c
      return
      end
c-----------------------------------------------------------------------
      double precision function frecip(y,mq2m,nx,xgr,coef)
c-----------------------------------------------------------------------
c
c     Calcul de la fonction "reciprocal power" au point "y"
c
      implicit double precision(a-h,o-z)
      dimension xgr(nx),coef(nx)
      logical first
      save fac1,fac2
      data first/.true./
c
      if (first) then
	first = .false.
        fac1 = 4.d0 / ((mq2m+1.d0)*(mq2m+2.d0))
        fac2 = (mq2m+1.d0) / (mq2m+3.d0)
      end if
c
      sum = 0.d0
      do 30 i=1,nx
        if (xgr(i).gt.y) then
          xsup = xgr(i)
          xinf = y
        else
	  xsup = y
	  xinf = xgr(i)
        end if
        q2m = fac1*xsup**(-mq2m-1) * (1.d0 - fac2*xinf/xsup)  
        sum = sum + coef(i)*q2m
   30 continue
      frecip = sum
c
      return
      end
