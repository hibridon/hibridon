*  system:  NH3-He, scaled perturbation theory PES from M. P. Hodges
*  and R. J. Wheatley, J. Chem. Phys. 14 , 8836 (2001).
*  subr to compute V vs. (R, theta, phi) provided by P. Stancil 
*  (U. GA), June-2011
*
*  C3 axis of NH3 defined to lie along z axis.  One N-H bond lies
*  in xz plane.
c
c Wrapper inherited from pot_ch3he_ccsdt.f  (q. ma, 20-jun-2011)
c
*  Note:  this pot routine requires a data file to be in hibxx/bin/progs/potdata:
*         

      include "common/syusr"
      include "common/ground"
      include "common/bausr"
c
* ------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(21)
      include "common/parpot"
      integer i
      double precision PI, DEG, S4PI
      parameter (PI=3.14159265358979323846d0, DEG=PI/180d0,
     +S4PI=2d0*dsqrt(PI))
      potnam='Wheatley/Hodges NH3-He PES'
      print *, potnam
1     print *, 'R (bohr):'
      read (5, *, end=99) r
*
* if r < 0, then output list of vlm's into a file (after 1st call)
      if (r .le. 0.d0) goto 99
      call pot(vv0,r)
      write (6, 100) vv0 * S4PI, vvl
100   format(' v',/,22(1pe12.4))
      goto 1
99    rr=3.5d0
      dr=0.5d0
      open (unit=12,file='nh3he_vlms.dat')
      write(12,109)
109   format(' R/bohr v00  v10  v20  v30  v40  v50  v60  v70',
     :   '  v80  v90  v33  v43  v53  v63  v73  v83  v93',
     :   '  v66  v76  v86  v96  v99')
      do i=1,40
        call pot(vv0,rr)
        write (12,110) rr, vv0 * S4PI, (vvl(j),j=1,22)
110     format(f7.2,23(1pe16.8))
        rr = rr + dr
      enddo
      close(12)
      end
! 99    end
c
* ------------------------------------------------------------------------
c
      subroutine loapot(iunit,filnam)
      implicit double precision (a-h,o-z)
      character*(*) filnam
      include "common/parbas"
      include "common/parpot" 
      common /conlam/ nlam, nlammx, lamnum(2)
      common /cosysi/ nscode, isicod, nterm
      potnam='Wheatley/Hodges NH3-He PES'
*
      nterm = 4
      mproj(1) = 0
      mproj(2) = 3
      mproj(3) = 6
      mproj(4) = 9
      lammin(1) = 1
      lammin(2) = 3
      lammin(3) = 6
      lammin(4) = 9
      lammax(1) = 9
      lammax(2) = 9
      lammax(3) = 9
      lammax(4) = 9
*
      ipotsy = 3
*
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
*  calculate total number of anisotropic terms
      nlam = 0
      do i=1,nterm
        lmin = lammin(i)
        lmax = lammax(i)
        do lb = lmin, lmax
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
*    vvl:        vector of length 21 to store r-dependence of each term
*                in potential expansion
*  variable in common block /coloapot/
*    s4pi:       normalization factor for isotropic potential
*
*  uses linear least squares routines from lapack
*
      implicit double precision (a-h,o-z)
c Define the sizes of grids
c NLMTRM: number of lambda-mu tuples for the fit
c NTHETA, NPHI: number of theta and pi
c NANGLE: number of theta-phi tuples
c LWORK, LIWORK: size of working (temp) arrays
      integer NLMTRM, NTHETA, NPHI, NANGLE, LWORK, LIWORK
      parameter (NLMTRM=22, NTHETA=18, NPHI=7, NANGLE=NTHETA*NPHI)
      parameter (LWORK=NANGLE*NLMTRM, LIWORK=LWORK/3)
c
      dimension iwork(LIWORK),ylm(NANGLE,NLMTRM)
      dimension swork(NLMTRM), work(LWORK)
      integer info
*
      common /covvl/ vvl(NLMTRM-1)
      double precision vsp(NANGLE)
      integer i, j, k, ind
      integer llist(NLMTRM), mlist(NLMTRM)
      data llist/0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 3, 4, 5, 6, 7, 8, 9, 6,
     +7, 8, 9, 9/
      data mlist/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 6,
     +6, 6, 6, 9/
      double precision theta, phi
      integer lambda, mu
      logical isfst
      data isfst/.TRUE./
*
      double precision stylmc
c This is a C/gsl function that calculate Ylm terms
c
      double precision PI, DEG, S4PI
      parameter (PI=3.14159265358979323846d0, DEG=PI/180d0,
     +S4PI=2d0*dsqrt(PI))
      double precision YLMCAL(NANGLE,NLMTRM)
c
      if (isfst) then
       open (unit=10, file=
     +   'potdata/pot_nh3he_wheatley_ylm.txt')
       read (10, *) YLMCAL
  97   close (10)
       isfst = .false.
      endif
c
c Coefficients for Ylm terms - all (lambda,mu) combinations of lambda<=9,
c   with mu a multiple of 3.
*    th = [5:10:175];   %grid of 19 angles
*    ph = [0:10:60];    %grid of  7 angles
*  All values of Ylm have been normalized
*  Note:  The ylm terms for mu .gt. 0 are defined as:
*           [Y(lambda,mu) + (-1)^mu * Y(lambda, -nu)]/2
!       do i = 1, NTHETA
!        do j = 1, NPHI
!         ind = NTHETA * (j - 1) + i
!         theta = (i - 1) * 1d1 + 5d0
!         phi = (j - 1) * 1d1
!         do k = 1, NLMTRM
!          lambda = llist(k)
!          mu = mlist(k)
!          ylm(ind, k) = stylmc(lambda, mu, theta, phi)
! c The stylmc function is implemented by C/gsl under gsl_legendre.c
!         enddo
!        enddo
!       enddo
!       write (6, 99) ylm
!   99  format (693(4(d17.10, 1x), /))
      ylm = YLMCAL
c Get potential values at distance r
      do i = 1, NTHETA
       do j = 1, NPHI
        ind = NTHETA * (j - 1) + i
        theta = ((i - 1) * 1d1 + 5d0) * DEG
        phi = (j - 1) * 1d1 * DEG
        vsp(ind) = potcm(r, theta, phi)
       enddo
      enddo
c Linear least-square fit
      rcond = 1.d-6
      call dgelsd(NANGLE, NLMTRM, 1, ylm, NANGLE, vsp, NANGLE,
     +swork, rcond, irank, work, LWORK, iwork, info)
*  save terms
      call dcopy(NLMTRM - 1, vsp(2), 1, vvl, 1)
      vv0 = vsp(1) / S4PI
c Long-range interactions have been taken care of in potcm function
      return
      end
c
* ---------------------------------------------------
*  system:  NH3-He, scaled perturbation theory PES from M. P. Hodges
*  and R. J. Wheatley, J. Chem. Phys. 14 , 8836 (2001).
*  subr to compute V vs. (R, theta, phi) provided by P. Stancil 
*  (U. GA), June-2011
*
*  original name of subroutine:  potnh3.f

cc potential subroutine for NH3-He, JCP, 114, 8836(2001)
cc the pes of JCP is in Cartesian coordinates
cc do the coordinate transformation in Body Fixed frame, origin is
cc on the center of mass of NH3
cc Aug. 07, 2007
cc -----------------Atomic Units------------------
cc inversion angle is fixed at experimental equilibrium geometry
c alpha_e=112.14, r_e=1.9132 bohr
cc
      double precision function potcm(r, theta, phi)
      implicit real*8 (a-h,o-z)
      parameter(pi=3.1415926d0, alphae=111.6d0)
         zn=0.128d0
         rex=r*r+zn*zn-2.0d0*r*zn*cos(theta)
         rn=dsqrt(rex)
         athe=acos((r*cos(theta)-zn)/rn)
         aphi=phi
         alpha=alphae/180.0d0*pi
         x=rn*sin(athe)*cos(aphi)
         y=rn*sin(athe)*sin(aphi)
         z=rn*cos(athe)
         potcm=pot4(x,y,z,alpha)
      return
      end
ccc-------------------------------------------------------------
ccc JCP, 114, 8836(2001)

      double precision function pot4(x,y,z,theta)
c
c Calculates He-NH3 intermolecular potential in terms of Cartesian
c coordinates x, y, z (Bohr) of He relative to N as centre, with the NH3
c symmetry axis being the z axis, and one H atom lying in the zx plane at
c negative z and positive x.  The angle theta is the NH3 umbrella angle
c (Radians), and theta=pi/2 when the NH3 is planar.  The intramolecular
c distortion energy of NH3 is not included.
c 
      implicit double precision(a-h,o-z)
      data ifirst,thlast/0,0d0/
      save
      if (ifirst.eq.0.or.theta.ne.thlast) then
      ifirst=1
      thlast=theta
      call potini(theta)
      endif
      pot4=potcal(x,y,z)
      end

      subroutine potini(t)
      implicit double precision(a-h,o-z)
      common/cnh3he/x1,y1,z1,x2,y2,z2,x3,y3,z3,
     1 p1,p2,p3,p4,p5,p6,p7,p8,b,x4,y4,z4,x5,y5,z5,x6,y6,z6,
     1 qn1,qn2,qn3,qn4,qn5,qn6,qn7,qn8,qn9,qn10,qh1,qh2,qh3,
     1 qh4,qh5,qh6,qc1,qc2,alpha,r
      save
      data p11,p21,p31,p41,p51,p61,p71,p81,p12,p22,p32,p42,p52,p62,p72
     1,p82,p13,p23,p33,p43,p53,p63,p73,p83,h1,h2,h3,h4/
     a 11.6523488375684d0,-2.98158118204410d0,1.49393389265970d0,
     b 17.2811199030927d0,379.968752899844d0,-17.3581436591624d0,
     c 0.261720480745098d0,29.5523802147692d0,-7.12983123895760d0,
     d 6.60000388701621d0,-0.404967288493762d0,78.7087154209477d0,
     e -534.987689713811d0,-62.3083854882989d0,1.93535982924854d0,    
     f -10.1837463629195d0,-1.09198084573485d0,1.75627013552562d0,
     g -1.40704330357101d0,286.638844403380d0,-223.587329773642d0,
     h 118.110532011749d0,0.668417880283956d0,-6.22354215723589d0,
     i 1.8726d0,0.224d0,-0.18d0,0.45d0/
      data qn1a,qn2a,qn3a,qn4a,qn5a,qn6a,qn7a,qn8a,qn9a,qn10a,qh1a,qh2a,
     1 qh3a,qh4a,qh5a,qh6a,qc1a,qc2a,qn1b,qn2b,qn3b,qn4b,qn5b,qn6b/
     a 105.358582685727d0,-93.5137650766694d0,-20.0556731797762d0,
     b 53.9601454589709d0,-13.3141673258754d0,-26.7620176122042d0,
     c 26.6988260179006d0,4.76344131933886d0,-13.4722471409789d0,
     d 3.55461655537341d0,77.7307451245117d0,83.4124804317485d0,
     e 10.8174097511536d0,-22.9172044461988d0,-22.8113636404374d0,
     f -2.74141173554461d0,-103.846480504968d0,31.2016079201890d0,
     g -124.151847166268d0,417.282168854078d0,22.9881057315202d0,
     h -254.773867022760d0,-18.5021336999981d0,34.8299600687467d0/
      data qn7b,qn8b,qn9b,qn10b,qh1b,qh2b,qh3b,qh4b,qh5b,qh6b,qc1b,qc2b,
     1 qn1c,qn2c,qn3c,qn4c,qn5c,qn6c,qn7c,qn8c,qn9c,qn10c,qh1c,qh2c/
     a -96.5301259714376d0,-5.06116647156490d0,62.5367249014338d0,
     b -5.03407376699773d0,-137.948336849227d0,-162.002983441311d0,
     c -22.7412630993916d0,33.5154303928627d0,35.3285466381656d0,
     d 4.55272224376798d0,180.648860035816d0,-45.1230870070718d0,
     e 331.449373618478d0,-597.207250158190d0,-43.2352558277175d0,
     f 274.524585476137d0,203.911154564190d0,-72.3569784341317d0,
     g 120.361515733561d0,7.26811312052723d0,-64.9599619089990d0,
     h -25.1531925783774d0,345.667422849100d0,411.949407885026d0/
      data qh3c,qh4c,qh5c,qh6c,qc1c,qc2c/
     a 60.0114831548612d0,-72.5666681773184d0,-80.7346698141146d0,
     b -11.0424451043878d0,-441.365662544351d0,94.6162278421202d0/
      b=2.7d0
      alpha=1.89d0
      c=cos(t)
      c2=c*c
      p1=p11+c2*(p12+p13*c2)
      p2=p21+c2*(p22+p23*c2)
      p3=p31+c2*(p32+p33*c2)
      p4=p41+c2*(p42+p43*c2)
      p5=p51+c2*(p52+p53*c2)
      p6=p61+c2*(p62+p63*c2)
      p7=p71+c2*(p72+p73*c2)
      p8=p81+c2*(p82+p83*c2)
      qn1=qn1a+c2*(qn1b+c2*qn1c)
      qn2=qn2a+c2*(qn2b+c2*qn2c)
      qn3=qn3a+c2*(qn3b+c2*qn3c)
      qn4=qn4a+c2*(qn4b+c2*qn4c)
      qn5=qn5a+c2*(qn5b+c2*qn5c)
      qn6=qn6a+c2*(qn6b+c2*qn6c)
      qn7=qn7a+c2*(qn7b+c2*qn7c)
      qn8=qn8a+c2*(qn8b+c2*qn8c)
      qn9=qn9a+c2*(qn9b+c2*qn9c)
      qn10=qn10a+c2*(qn10b+c2*qn10c)
      qh1=qh1a+c2*(qh1b+c2*qh1c)
      qh2=qh2a+c2*(qh2b+c2*qh2c)
      qh3=qh3a+c2*(qh3b+c2*qh3c)
      qh4=qh4a+c2*(qh4b+c2*qh4c)
      qh5=qh5a+c2*(qh5b+c2*qh5c)
      qh6=qh6a+c2*(qh6b+c2*qh6c)
      qc1=qc1a+c2*(qc1b+c2*qc1c)
      qc2=qc2a+c2*(qc2b+c2*qc2c)
      r=h1+c2*(h2+c2*(h3+c2*h4))
      x1=r*sin(t)
      y1=0
      z1=r*c
      x2=-x1/2
      y2=x1*sqrt(0.75d0)
      z2=z1
      x3=x2
      y3=-y2
      z3=z2
      x4=x1/2
      y4=y1/2
      z4=z1/2
      x5=x2/2
      y5=y2/2
      z5=z2/2
      x6=x3/2
      y6=y3/2
      z6=z3/2
      end
      double precision function potcal(x,y,z)
      implicit double precision(a-h,o-z)
      common/cnh3he/x1,y1,z1,x2,y2,z2,x3,y3,z3,
     1 p1,p2,p3,p4,p5,p6,p7,p8,b,x4,y4,z4,x5,y5,z5,x6,y6,z6,
     1 qn1,qn2,qn3,qn4,qn5,qn6,qn7,qn8,qn9,qn10,qh1,qh2,qh3,
     1 qh4,qh5,qh6,qc1,qc2,alpha,rnh
      s=(x-x1)**2+(y-y1)**2+(z-z1)**2
      r=sqrt(s)
      br=b*r
      ebr=exp(-br)
      damp6=1+br*(1+br*(1+br*(1+br*(1+br*(1+br/6)/5)/4)/3)/2)
      damp8=damp6+br**7*(1+br/8)/5040
      potcal=(-p7*(1-damp6*ebr)-p8*(1-damp8*ebr)/s)/s**3
      y10=((x-x1)*x1+(y-y1)*y1+(z-z1)*z1)/(rnh*r)
      y20=3*y10*y10-1
      c00=qh1+qh4*r
      c10=qh2+qh5*r
      c20=qh3+qh6*r
      potrep=(c00-c10*y10+c20*y20)*exp(-alpha*r)
      s=(x-x2)**2+(y-y2)**2+(z-z2)**2
      r=sqrt(s)
      br=b*r
      ebr=exp(-br)
      damp6=1+br*(1+br*(1+br*(1+br*(1+br*(1+br/6)/5)/4)/3)/2)
      damp8=damp6+br**7*(1+br/8)/5040
      potcal=potcal-(p7*(1-damp6*ebr)+p8*(1-damp8*ebr)/s)/s**3
      y10=((x-x2)*x2+(y-y2)*y2+(z-z2)*z2)/(rnh*r)
      y20=3*y10*y10-1
      c00=qh1+qh4*r
      c10=qh2+qh5*r
      c20=qh3+qh6*r
      potrep=potrep+(c00-c10*y10+c20*y20)*exp(-alpha*r)
      s=(x-x3)**2+(y-y3)**2+(z-z3)**2
      r=sqrt(s)
      br=b*r
      ebr=exp(-br)
      damp6=1+br*(1+br*(1+br*(1+br*(1+br*(1+br/6)/5)/4)/3)/2)
      damp8=damp6+br**7*(1+br/8)/5040
      potcal=potcal-(p7*(1-damp6*ebr)+p8*(1-damp8*ebr)/s)/s**3
      y10=((x-x3)*x3+(y-y3)*y3+(z-z3)*z3)/(rnh*r)
      y20=3*y10*y10-1
      c00=qh1+qh4*r
      c10=qh2+qh5*r
      c20=qh3+qh6*r
      potrep=potrep+(c00-c10*y10+c20*y20)*exp(-alpha*r)
      s=(x-x4)**2+(y-y4)**2+(z-z4)**2
      r=sqrt(s)
      c00=qc1+qc2*r
      potrep=potrep+c00*exp(-alpha*r)
      s=(x-x5)**2+(y-y5)**2+(z-z5)**2
      r=sqrt(s)
      c00=qc1+qc2*r
      potrep=potrep+c00*exp(-alpha*r)
      s=(x-x6)**2+(y-y6)**2+(z-z6)**2
      r=sqrt(s)
      c00=qc1+qc2*r
      potrep=potrep+c00*exp(-alpha*r)
      s=x**2+y**2+z**2
      r=sqrt(s)
      br=b*r
      ebr=exp(-br)
      damp6=1+br*(1+br*(1+br*(1+br*(1+br*(1+br/6)/5)/4)/3)/2)
      damp8=damp6+br**7*(1+br/8)/5040
      y20=z*z/s
      y33=x*(x*x-3*y*y)/(r*s)
      y40=y20*y20
      y53=y20*y33
      c6=p1+p2*y20+p3*y33
      c8=p4+p5*y20+p6*y33
      potcal=potcal-(c6*(1-damp6*ebr)+c8*(1-damp8*ebr)/s)/s**3
      c00=qn1+qn6*r
      c20=qn2+qn7*r
      c33=qn3+qn8*r
      c40=qn4+qn9*r
      c53=qn5+qn10*r
      potrep=potrep+(c00+c20*y20+c33*y33+c40*y40+c53*y53)*exp(-alpha*r)
      potcal=potrep+potcal
      end
