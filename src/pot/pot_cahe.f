*  System:  Ca(4s5p) + He
*  References
*   B. Pouilly, J.-M. Robbe, and M. H. Alexander, J. Chem. Phys. 91, 1658 (1989).
*   T. Duhoo and B. Pouilly, J. Chem. Phys. 101, 7554 (1994).
      include "common/syusr"
      include "common/bausr"
      include "common/ground"
* --------------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      common /coselb/ ibasty
      include "common/parbas"
      include "common/parpot"
      potnam='Ca(4s5p)-He'
      ibasty=7
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
      common /cosysi/ junk(5), npot
      common /cosysr/ isrcod, junkr, en(4), de(4), re(4), be(4),
     :                        rl(4), cl(4), cmix, alphg, rgaus,
     :                        agaus, demor, remor, bemor, dissmor
      common /covvl/ vvl(4)
      include "common/parpot"
      potnam='Ca(4s5p)-He'
      print *, potnam
      print *, 
     + ' npot = 1:  Pouilly-Robbe-Alexander, JCP 91, 1658 (1989)'
      print *, ' npot = 2:  Pouilly JCP '
      print *, 
     + ' npot = 3:  Czuchaj, Chem. Phys. Lett. 182, 191 (1991)'

1      print *, ' npot, r (bohr)?  '
      read (5, *, end=94) npot, r
      call pot(vv0,r)
      write (6, 100) vvl
100   format(' v1pi =',1pe16.8,'; v1sig =',1pe16.8,/,
     :       ' v3pi =',1pe16.8,'; v3sig =',1pe16.8)
      goto 1
94    end
* --------------------------------------------------------------------------
      block data cahepot
      implicit double precision (a-h,o-z)
      common /cosysr/ isrcod, junkr, en(4), de(4), re(4), be(4),
     :                        rl(4), cl(4), cmix, alphg, rgaus,
     :                        agaus, demor, remor, bemor, dissmor
      data en /0,7.076,27.44,183.9 /
      data de /1.d-4,8.d-6,1.d-4,8.d-6/
      data re /9,13.8,9,13.8/
      data be /.4,.47,.4,.45/
      data rl /13,18,13,18/
      data cl /40,50,40,50/
      data cmix, rgaus, agaus, alphg, demor, remor, bemor, dissmor
     :  / .9747,8.1,1E-4,7.5E-2,3.5E-3,6.5,.43,4.3E-3/
      end
      
* --------------------------------------------------------------------------
      subroutine pot (vv0, r)
*  subroutine to calculate the rz-dependent coefficients in the
*  ca-he msv potentials of pouilly et al.
*  in units of hartree for energy and bohr for distance
* ----------------------------------------------------------------------
*  on entry:
*    rz:      interparticle distance
*  on return:
*  vv0        (not used)
*  vvl(1,2)   contains the 3pi and 3sig potentials
*  vvl(3,4)   contains the 1pi and 1sig potentials
*  variable in common block /conlam/ used here
*    nlam:    the number of angular coupling terms actually used
*  variable in common block /covvl/
*    vvl:     array to store r-dependence of each
*             angular term in the potential
*  variable in common block /coconv/
*   econv:     conversion factor from cm-1 to hartrees
*  variables in common /cosysr/
*    isrcod:  number of real parameters
*    en(1):    asymptotic energy of j=0 fine-structure level of triplet
*    en(2):    asymptotic energy of j=1 fine-structure level of triplet
*    en(3):    asymptotic energy of j=2 fine-structure level of triplet
*    en(4):    asymptotic energy of singlet state (cm-1)
*    cmix:     mixing coefficient of j=1 singlet and triplet levels
*    de, re, be, rl, c:  msv parameters for 3pi, 3sig, 1pi, 1sig
 
* author:  millard alexander and brigitte pouilly
* latest revision date:  17-may-1994 by bp
* ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      common /conlam/ nlam
      common /covvl/ vvl(4)
      common /coconv/ econv
      common /cosysi/ nscode, isicod, nterm, nstate, ipol, npot
      common /cosysr/ isrcod, junkr, en(4), de(4), re(4), be(4),
     :                        rl(4), cl(4), cmix, alphg, rgaus,
     :                        agaus, demor, remor, bemor, dissmor
      data onethr, fivthr /0.3333333333333333d0, 1.66666666666667d0/
      if (npot .lt. 1 .or. npot .gt. 3) then
        write (6, 5) npot
5       format (' *** NPOT = ',i2, ' NOT ALLOWED IN POT_13P; ABORT')
        call exit
      endif
      if (npot. eq. 1) then
* here for original msv potentials
        call potmsv1(r,vp,vs,re(1),de(1),be(1),rl(1),cl(1))
        vv0=0.d0;
c Commented ou by J. Klos 2012.11.22
c        vvl(1)=vp
c        vvl(2)=vs
c The coupling formula which is coded is identical to
c Aquilanti's formula Eq. 2.10 in Aquilanti, Grossi
c J. Chem. Phys. 73,1165 (1980)
c where coupling elements multiply v_mu, not v_pi and v_sigma
c vvl(1) is v_0, isotropic part
c vvl(2) is v_2, anisotropic part
        vvl(1)=(vs+2.d0*vp)/3.d0
        vvl(2)=5.d0*(vs-vp)/3.d0
        vvl(3)=(vs+2.d0*vp)/3.d0
        vvl(4)=5.d0*(vs-vp)/3.d0
      else if (npot. eq. 2) then
* here for new msv potentials + coupling for the 1sigma state
        call pot1(r,vp1,vs1,vp3,vs3)
      else if (npot. eq. 3) then
* here for czuchaj potentials
        call pot2(r,vp1,vs1,vp3,vs3)
      endif
      if (npot .gt. 1) then
        vvl(1)=(vs1+2.d0*vp1)/3.d0
        vvl(2)=5.d0*(vs1-vp1)/3.d0
        vvl(3)=(vs3+2.d0*vp3)/3.d0
        vvl(4)=5.d0*(vs3-vp3)/3.d0
      endif
      return
      end
      subroutine potmsv1 (r, vp, vs, re, de, be, rl, c)
 
c -----------------------------------------------------------------------
c
c  subroutine to calculate the r-dependent coefficients in the
c  ca-he inelastic problem
c  in units of hartree for energy and bohr for distance
c
c  on return:
c  vp contains the pi potential
c  vs contains the sigma potential
c
c -----------------------------------------------------------------------
c
c  morse-spline-van der waals potential used for all curves
c
c      if r .le. re,  v(r) = morse curve with parameters re, de, and be
c      if re .lt. r .le. rc,  v(r) = spline function
c      if r. ge. rl, v(r) = - c r**(-nn)
c
c  in smp notation we have:
c
c   vmorse[r,re,de,be] : de*exp[-(be*(r + -re))]*(-2 + exp[-(be*(r\
c                        + -re))])
c   vspline[r,re,de,a,b] : -de + (a + b*(r + -re))*(r + -re)^2
c       with
c        a[re,rl,de,c,n] : (3((de + -(c*rl^(-n)))/(-re + rl)) + -(c\
c                          *n*rl^(-1 + -n)))/(-re + rl)
c       and
c        b[re,de,rl,c,n] : (-2((de + -(c*rl^(-n)))/(-re + rl)) + c\
c                          *n*rl^(-1 + -n))/(-re + rl)^2
c
c   the simplified expression for vspline with the explicit values of
c   the coefficients a and b (preceding lines) incorporated is
c
c   vspline[r,re,de,rl,c,n] : (-(de*(-rl + r)^2*(-rl + -2r + 3re))\
c                            + c*rl^(-1 + -n)*(r + -re)^2*(-2(rl*r) + -(rl*re)\
c                            + n*(-rl + r)*(-rl + re) + 3(rl^2)))/(-rcl + re)^3
c
c
c   vlr[r,c,n] : -(c*r^(-n))
c
c
c   note that all parameters with energy should
c   be in hartree, all parameters with distance units [re, rl] should
c   be in bohr, and all parameters with inverse distance units [be]
c   should be in inverse bohr
c
c -----------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
      dimension re(2), de(2), be(2), rl(2), c(2)
      dimension v(2)
c
c -----------------------------------------------------------------------
c now determine potential
c
      do 30  i = 1, 2
c
c
c  i = 1 for pi; i=2 for sigma
 
      if (r .le. re(i) ) then
c
c  here for morse potential
c
       v(i) = de(i) * exp(- be(i) * (r - re(i))) * (-2. +
     :         exp(- be(i) * (r - re(i))))
c
      else if (r .gt. rl(i) ) then
c
c  here for long-range potential
c
       v(i) = -c(i) * r ** (-6)
c
      else if (r. ge. re(i) .and. r .le. rl(i)) then
c
         v(i) = (- de(i) * (r  - rl(i)) ** 2 * (-2 * r + 3 * re(i)
     :       - rl(i)) + c(i) * rl(i) ** (-1 - 6) * (r - re(i)) ** 2
     :       * (-2. * r * rl(i)  - re(i) * rl(i) + 6 * (r - rl(i))
     :       * (re(i) - rl(i))+ 3 * rl(i)** 2)) / (re(i) - rl(i)) ** 3
      end if
30    continue
c
c -----------------------------------------------------------------------
c  now set up values to be returned to calling program
c
       vp = v(1)
       vs = v(2)
       return
       end
* ---------------------------------------------------------------------
* pot1.f
      subroutine pot1(r,vp1,vs1,vp3,vs3)
      implicit double precision (a-h,o-z)
      common /cosysr/ isrcod, junkr, en(4), de(4), re(4), be(4),
     :                        rl(4), cl(4), cmix, alphg, rgaus,
     :                        agaus, demor, remor, bemor, dissmor
      call potmsv(vp3,vs3,vp1,vs1,r)
      call potmorse (vs2, r)
 
*  here vs1 contains the msv part of the excited singlet sigma state
*  and vs2 contains the morse part of the excited singlet sigma state
 
*  here to obtain a gaussian form for the coupling
*  coupsi(1) = h0 * exp(- alpha (r - rmax)**2)
*  to obtain a maximum of coupsi at r = 8.1 bohr
*  alpha = 0.075
*      coupling = coupsi(1) * dexp( -(0.075) * ( r - 8.1)**2)
      coupling = agaus*dexp(-alphg*(r-rgaus)**2)
      vpert = (0.5d0 *(vs1 + vs2) - 0.5d0 * dsqrt((vs1 - vs2)**2 +
     #          4.d0 * coupling**2))
      vs1=vpert
      return
      end
 
      subroutine potmsv(vp,vs,vp1,vs1,r)
      implicit double precision (a-h,o-z)
      common /cosysr/ isrcod, junkr, en(4), de(4), re(4), be(4),
     :                        rl(4), c(4)
      dimension v(4)
      do 30 j=1,4
c
c     j=1 for triplet pi
c     j=2 for triplet sigma
c     j=3 for singlet pi
c     j=4 for singlet sigma
c
      if(r.le.re(j)) then
      v(j)=de(j)*dexp(-be(j)*(r-re(j)))*(-2.+dexp(-be(j)*(r-re(j))))
      else if (r.gt.rl(j)) then
      v(j)=-c(j)*r**(-6)
      else if (r.gt.re(j).and.r.le.rl(j)) then
      v(j)=(-de(j)*(r-rl(j))**2*(-2*r+3*re(j)
     :       -rl(j))+c(j)*rl(j)**(-7)*(r-re(j))**2
     :       *(-2.*r*rl(j)-re(j)*rl(j)+6*(r-rl(j))
     :       *(re(j)-rl(j))+3*rl(j)**2))/(re(j)-rl(j))**3
      end if
 30   continue
c
      vp=v(1)
      vs=v(2)
      vp1=v(3)
      vs1=v(4)
c
      return
      end
      subroutine potmorse ( vs, r)
c -----------------------------------------------------------------------
c
c  subroutine to calculate the r-dependent coefficients for the morse
c  part of the excited singlet sigma state of the ca - he system.
c  units are in hartree or cm-1 for energy and bohr for distance
 
 
c  on return:
c  vs contains the morse potential
c -----------------------------------------------------------------------
 
c  in smp notation we have:
c
c   vmorse[r,re,de,be] : de*exp[-(be*(r + -re))]*(-2 + exp[-(be*(r\
c                        + -re))])
c
c   de should be in hartree, re in bohr and de in inverse bohr
 
c   note that the morse potential is  shifted by an amount of
c   energy equal to delrsi
c   e.g.
c       v(r) = v(r) - delrsi
c
c -----------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
      common /cosysr/ isrcod, junkr, en(4), de(4), re(4), be(4),
     :                        rl(4), cl(4), cmix, alphg, rgaus,
     :                        agaus, demor, remor, bemor, dissmor
c
c -----------------------------------------------------------------------
c now determine potential
c
 
*  i = 1 for the morse potential
*  coup is the coupling between the  morse and the msv parts
*  of the potential
 
c
       v = demor * dexp(- bemor * (r - remor)) * (-2.d0 +
     :         dexp(- bemor * (r - remor)))
 
       v = v + dissmor
c
c -----------------------------------------------------------------------
c  now set up values to be returned to calling program
c
       vs = v
       return
       end
c -----------------------------------------------------------------------
c pot2
      block data blk2
c -----------------------------------------------------------------------
* nouveau fit des potentiels
 
      implicit double precision (a-h,o-z)
      parameter (n1si = 32,n3si =36)
      parameter (n1pi = 34,n3pi =32)
      common /ppot1pi/ x1pi(n1pi),vl1pi(n1pi),dvlf1pi,dvli1pi,
     #                dvl1pi(n1pi)
      common /ppot3pi/ x3pi(n3pi),vl3pi(n3pi),dvlf3pi,dvli3pi,
     #                dvl3pi(n3pi)
      common /ppot1si/ x1si(n1si),vl1si(n1si),dvlf1si,dvli1si,
     #                dvl1si(n1si)
      common /ppot3si/ x3si(n3si),vl3si(n3si),dvlf3si,dvli3si,
     #                dvl3si(n3si)
 
       data x1si /
     # 4.00,   4.13,   4.24,   4.40,
     # 4.50,   5.00,   5.50,   6.00,
     # 6.50,   7.00,   7.50,   8.00,   8.50,
     # 9.00,  10.00,  11.00,  12.00,  14.00,
     #16.00,  18.00,  20.00,  22.00,  25.00,
     #30.00,  30.62,  31.66,    34.72,
     #35.75,  36.75,  37.74,  38.78,  39.81/
 
       data x3si /
     # 4.00, 4.15,  4.24,
     # 4.50, 4.73,  5.00,   5.50,   6.00,
     # 6.50,   7.00,   7.50,   8.00,   8.50,
     # 9.00,  10.00,  11.00,  12.00,  14.00,
     #16.00,  18.00,  20.00,  22.00,  25.00,
     #30.00,  30.26,  30.82,  31.50,  32.41,
     #33.45,  34.84,  35.83,  36.51,  36.87,
     #37.46,  38.34,  39.14/
 
       data x1pi /
     # 4.00,   4.12,   4.19,   4.33,   4.42,
     # 4.50,   4.61,   5.00,   5.50,   6.00,
     # 6.50,   7.00,   7.50,   8.00,   8.50,
     # 9.00,  10.00,  11.00,  12.00,  14.00,
     #16.00,  18.00,  20.00,  22.00,
     #30.00,  30.62,  31.66,  32.69,  33.68,
     #34.72,  35.75,  36.75,  37.74,  38.78/
 
       data x3pi /
     # 4.00,   4.12,   4.19,   4.32,   4.40,
     # 4.50,   4.58,   5.00,   5.50,   6.00,
     # 6.50,   7.00,   7.50,   8.00,   8.50,
     # 9.00,  10.00,  11.00,  12.00,  14.00,
     #16.00,  18.00,  20.00,
     #30.00,  32.69,  33.68,  34.72,  35.75,
     #36.75,  37.62,  38.50,  39.33/
 
* 1sigma potential
       data (vl1si(i), i =1, n1si) /
     # .01241803,  .01028470,  .00908363,  .00681495,
     # .00588107,  .00309300,  .00213408,  .00139403,
     # .00081205,  .00050616,  .00055909,  .00070214,  .00084519,
     # .00097418,  .00118613,  .00128102,  .00126314,  .00100422,
     # .00061917,  .00032711,  .00014615,  .00005722,  .00001121,
     #-.00000095,
     #-1.02278e-06, -1.07616e-06,  -9.26690e-07,
     #-7.87900e-07,  -6.27758e-07,  -4.46263e-07, -2.43416e-07,
     #-2.98933e-08/
 
* 1pi potential
       data (vl1pi(i), i =1, n1pi) /
     # .00929594,  .00828205,  .00770513,  .00661538,  .00603846,
     # .00555086,  .00507692,  .00405383,  .00288582,  .00184798,
     # .00104594,  .00052500,  .00022697,  .00007987,  .00001597,
     #-.00000906, -.00002718, -.00002909, -.00003099, -.00003004,
     #-.00001717, -.00000715, -.00000501, -.00000215,
     #-.00000119,
     #-1.14615e-06, -1.05641e-06, -9.60256e-07, -8.44872e-07,
     #-7.23077e-07, -5.94872e-07, -4.60256e-07, -3.19231e-07,
     #-1.71795e-07/
 
*3sigma potential
       data (vl3si(i), i =1, n3si) /
     # .01611304,  .01402140,  .01268680,
     # .00977206,   .00828292, .00722814,  .00580907,  .00455213,
     # .00325894,  .00205517,  .00129008,  .00111008,  .00121117,
     # .00135994,  .00158715,  .00165010,  .00156498,  .00113606,
     # .00065804,  .00032902,  .00013900,  .00006294,  .00008416,
     #-.00001192,
     #-1.78803e-05, -2.76060e-05, -3.73317e-05, -4.63092e-05,
     #-5.07980e-05, -4.96758e-05, -4.44389e-05, -3.95761e-05,
     #-3.62095e-05, -2.98504e-05, -2.04988e-05, -1.07731e-05/
 
*3pi potential
       data (vl3pi(i), i =1, n3pi) /
     # .00855804,  .00751282,  .00693590,  .00584615,  .00526923,
     # .00466991,  .00430769,  .00311494,  .00211501,  .00147080,
     # .00100183,  .00060296,  .00024700,  .00005603, -.00002217,
     #-.00004506, -.00005817, -.00005412, -.00004911, -.00003719,
     #-.00002813, -.00001597, -.00001216,
     #-.00000906,
     #-9.22420e-06, -8.63701e-06, -7.67616e-06, -6.50178e-06,
     #-5.16726e-06, -3.83274e-06, -2.44484e-06, -1.05694e-06/
      end 
* ----------------------------------------------------------------------
      subroutine pot2(rz,vp1,vs1,vp3,vs3)
 
* ----------------------------------------------------------------------
*  on entry:
*    rz:      interparticle distance
 
*  distance is in bohr
*  energy is in hartree
 
*  on return:
 
*  ---------------------------------------------------------------------
*  subroutines called:
 
*   spline  (from numerical recipes)
*      to evaluate the second derivatives of the interpolating function
*      this will be best done if the two first (two last) values of r
*      are chosen close from each other in order to obtain a correct
*      first derivative for the lower (higher) value of r
 
*   splint  (from numerical recipes)
*      to return a cubic-spline interpolated value for a given value rz
 
*  variables in common block / copot/
*    contains the arrays of the singlet and triplet potentials
*  ---------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (n1si = 32,n3si =36)
      parameter (n1pi = 34,n3pi =32)
      common /coenerg/ en, rmu
      common /ppot1pi/ x1pi(n1pi),vl1pi(n1pi),dvlf1pi,dvli1pi,
     #                dvl1pi(n1pi)
      common /ppot3pi/ x3pi(n3pi),vl3pi(n3pi),dvlf3pi,dvli3pi,
     #                dvl3pi(n3pi)
      common /ppot1si/ x1si(n1si),vl1si(n1si),dvlf1si,dvli1si,
     #                dvl1si(n1si)
      common /ppot3si/ x3si(n3si),vl3si(n3si),dvlf3si,dvli3si,
     #                dvl3si(n3si)
      dimension vvl(4)
      dimension c6(3), rl2(3)

      data interp / 1/
 
* here are the parameters for extrapolation at long range ( -c *r^-6)
* for r > rl2 in order 1sigma,  1pi, 3pi
 
      data c6 / 851.44, 320.0,  556.32 /
      data rl2 / 30.4, 20.0, 18.3 /
 
* here are the parameters for the extrapolation at long range for
* the 3sigma state. ( alpha * exp(beta r))
 
      data alpha / 1.74322d0/
      data beta  /-0.4719d0/
      data rextr / 20.0/
 
 
*     save ddvl
*     save vl
 
      rmu2=2.*rmu
 
      if (interp .eq. 1) then
 
*  here for first pass through pot subroutine
 
*  interpolate the  coefficients
 
      ns1 = n1si -1
        dvlf1si =( vl1si(n1si) - vl1si(ns1)) / (x1si(n1si) - x1si(ns1))
        dvli1si =( vl1si(2) - vl1si(1)) / (x1si(2) - x1si(1))
        call dspline1 ( x1si, vl1si(1), n1si, dvlf1si, dvli1si, dvl1si)
      ns3 = n3si -1
        dvlf3si = (vl3si(n3si) - vl3si(ns3)) / (x3si(n3si) - x3si(ns3))
        dvli3si = (vl3si(2) - vl3si(1)) / (x3si(2) - x3si(1))
        call dspline1 ( x3si, vl3si(1), n3si, dvlf3si, dvli3si, dvl3si)
      nt1 = n1pi -1
        dvlf1pi =( vl1pi(n1pi) - vl1pi(nt1)) / (x1pi(n1pi) - x1pi(nt1))
        dvli1pi =( vl1pi(2) - vl1pi(1)) / (x1pi(2) - x1pi(1))
        call dspline1 ( x1pi, vl1pi(1), n1pi, dvlf1pi, dvli1pi, dvl1pi)
      nt3 = n3pi -1
        dvlf3pi = (vl3pi(n3pi) - vl3pi(nt3)) / (x3pi(n3pi) - x3pi(nt3))
        dvli3pi = (vl3pi(2) - vl3pi(1)) / (x3pi(2) - x3pi(1))
        call dspline1 ( x3pi, vl3pi(1), n3pi, dvlf3pi, dvli3pi, dvl3pi)
 
        interp = 0
      end if
 
* here is the entry point for subsequent entries into pot subroutine
 
      do 40 i = 1, 4
        vvl(i) = 0.
40    continue
* here for long range extrapolation
      if (rz .gt. rl2(1)) then
        vvl(1) = - c6(1) * rz**(-6)
        do i = 3, 4
        vvl(i) = - c6(i-1) * rz**(-6)
        end do
        vvl(2) = alpha * dexp( beta * rz)
      end if
 
* here for extrapolation in the intermediate region
      rl1 = 20.0
      if (rz .gt. rl1 .and. rz .le. rl2(1) ) then
        vvl(2) = alpha * dexp( beta * rz)
        do i = 3, 4
        vvl(i) = - c6(i-1) * rz**(-6)
        end do
        call dsplint1 (x1si, vl1si(1), dvl1si(1), n1si, rz, vvl(1))
      end if
 
      if (rz .gt. rl2(3) .and. rz .le. rl1 ) then
        vvl(4) = - c6(3) * rz**(-6)
        call dsplint1 (x1si, vl1si(1), dvl1si(1), n1si, rz, vvl(1))
        call dsplint1 (x3si, vl3si(1), dvl3si(1), n3si, rz, vvl(2))
        call dsplint1 (x1pi, vl1pi(1), dvl1pi(1), n1pi, rz, vvl(3))
      end if
 
* here for spline interpolation
      if ( rz .le. rl2(3) ) then
      call dsplint1 (x1si, vl1si(1), dvl1si(1), n1si, rz, vvl(1))
      call dsplint1 (x3si, vl3si(1), dvl3si(1), n3si, rz, vvl(2))
      call dsplint1 (x1pi, vl1pi(1), dvl1pi(1), n1pi, rz, vvl(3))
      call dsplint1 (x3pi, vl3pi(1), dvl3pi(1), n3pi, rz, vvl(4))
      end if
 
*     write (6, 1000) rz, vvl(1), vvl(2), vvl(3), vvl(4)
*1000  format (1x,f5.2,4(1x,e15.4))
 
      vs1 = vvl(1)
      vp1 = vvl(3)
      vs3 = vvl(2)
      vp3 = vvl(4)
      return
      end
* ------------------------------------
      subroutine dspline1(x,y,n,yp1,ypn,y2)
      implicit double precision(a-h,o-z)
      parameter (nmax=300)
      dimension x(n),y(n),y2(n),u(nmax)
      save
      if (n. gt. nmax) then
        write (6,5) n
        write (9,5) n
5       format (' *** N = ',i3, ' .GT. NMAX IN DSPLINE; ABORT')
        call exit
      endif
      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     *      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      end
      subroutine dsplint1(xa,ya,y2a,n,x,y)
      implicit double precision(a-h,o-z)
      dimension xa(n),ya(n),y2a(n)
      save
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
*      if (h.eq.0.d0) pause 'bad xa input.'
      if (h.eq.0.d0) then
        write (6, 15) h
        write (9, 15) h
15      format (' *** bad xa input in splint; abort')
        call exit
      endif
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     *      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
      return
      end
