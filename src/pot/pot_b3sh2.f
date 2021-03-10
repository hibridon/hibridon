*system: B(2s^2 3s)H2(J=0,1) Using Alexander ab initio potential and 
*  expansion of dubernet and hutson
* reference:  m.h. alexander and m. yang, j. chem. phys. 103, 7956 (1995)
*             m.-l. dubernet et al., J. Chem. Phys. 94, 7602 (1991).
      include "common/syusr"
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(5)
      include "common/parpot"
      potnam='ALEXANDER B(2S)H2(J=0,1) DUBERNET-HUTSON'
      print *, potnam
1      print *, ' r (bohr)'
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl(3)
100   format(' V000, V202',/,2(1pe16.8))
      goto 1
99    end
      include "common/ground"
      include "common/bausr"
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      common /coselb/ ibasty
      potnam='ALEXANDER B(2S)H2(J=0,1) DUBERNET-HUTSON'
      ibasty=12
      lammin(1)=1
      lammax(1)=3
      mproj(1)=0
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* --------------------------------------------------------------------------
      subroutine pot (vv0, rz)
* ----------------------------------------------------------------------
*  subroutine to calculate the r-dependent coefficients in the
*  b(2s)-h2 potentials of alexander, using body-frame 
*  expansion of dubernet and flower
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    rz:      interparticle distance
*  on return:
*  vv0       totally symmetric potential (v000)
*  variable in common block /covvl/
*    vvl:     vector of length 1 to store r-dependence of each term
*             in potential expansion
*  CS calculations case 1A (csflag=.true. and ihomo=.false)
*  body frame expansion coefficients
*    vvl(1)  v200(R)
*  CC calculations (csflag=.false.)
*    vvl(1)  V220(R)

*  subroutines used:
*   bint4,bvalu  (spline routines from cmlib bspline package)

* author:  millard alexander
* latest revision date:  20-dec-1995
* ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (nn=44, nnp6=nn+6, nw=5*nnp6)
      dimension rr(nn), v0(nn), v2(nn)
      dimension work(nw)
      dimension t0(nnp6),t2(nnp6),b0(nnp6),b2(nnp6)
      dimension a(2), alph(2), clr(2), vx(2)
      common /coered/ ered, rmu
      common /covvl/ vvl(5)

      data rr/
     : 3.4d0, 3.6d0, 3.8d0, 4.0d0, 4.2d0, 4.4d0, 4.6d0, 4.8d0, 5.0d0,
     : 5.2d0, 5.4d0, 5.6d0, 5.8d0, 6.0d0, 6.2d0, 6.4d0, 6.6d0, 6.8d0,
     : 7.0d0, 7.2d0, 7.4d0, 7.6d0, 7.8d0, 8.0d0, 8.2d0, 8.4d0, 8.6d0,
     : 8.8d0, 9.0d0, 9.2d0, 9.4d0, 9.6d0, 9.8d0, 1.0d1, 1.02d1, 1.04d1,
     : 1.06d1, 1.08d1, 1.1d1, 1.12d1, 1.14d1, 1.16d1, 1.18d1, 1.2d1/

      data v0 /
     : 1.0782976d-2, 8.1342629d-3, 6.4392151d-3, 5.3633356d-3,
     : 4.6225462d-3, 4.1134057d-3, 3.7648734d-3, 3.5267130d-3,
     : 3.3607782d-3, 3.2367904d-3, 3.1332747d-3, 3.0348723d-3,
     : 2.9316906d-3, 2.8186109d-3, 2.6924480d-3, 2.5538590d-3,
     : 2.4044841d-3, 2.2467111d-3, 2.0836807d-3, 1.9184882d-3,
     : 1.7539274d-3, 1.5926922d-3, 1.4369416d-3, 1.2883000d-3,
     : 1.1481967d-3, 1.0173795d-3, 8.9641116d-4, 7.8551862d-4,
     : 6.8459289d-4, 5.9344371d-4, 5.1161460d-4, 4.3858038d-4,
     : 3.7376541d-4, 3.1654371d-4, 2.6629327d-4, 2.2241700d-4,
     : 1.8432401d-4, 1.5142342d-4, 1.2312434d-4, 9.8868099d-5,
     : 7.8225021d-5, 6.0797670d-5, 4.6188619d-5, 3.4000447d-5/
      data v2 /
     : 6.7186575d-3, 5.2651291d-3, 4.4971122d-3, 4.0643579d-3,
     : 3.6850711d-3, 3.3216755d-3, 2.9699547d-3, 2.6380393d-3,
     : 2.3334848d-3, 2.0586940d-3, 1.8133300d-3, 1.5953844d-3,
     : 1.4021624d-3, 1.2307918d-3, 1.0785600d-3, 9.4323188d-4,
     : 8.2271719d-4, 7.1530775d-4, 6.1967405d-4, 5.3457572d-4,
     : 4.5905041d-4, 3.9221026d-4, 3.3325497d-4, 2.8147191d-4,
     : 2.3617769d-4, 1.9678487d-4, 1.6272983d-4, 1.3344517d-4,
     : 1.0835983d-4, 8.6938441d-5, 6.8788569d-5, 5.3552636d-5,
     : 4.0853536d-5, 3.0294628d-5, 2.1517769d-5, 1.4322106d-5,
     : 8.5461184d-6, 4.0282901d-6, 6.0711290d-7, -1.8778951d-6,
     : -3.5831261d-6, -4.6639418d-6, -5.2756986d-6, -5.5737553d-6/
      data a / 1.29982844255d0, 0.42377137115d0/
      data alph /1.40941734378d0, 1.21891345267d0/
      data clr / 1.01524789975d+2, -1.664314420225d+1/
      data interp /1/
      data zero,half,one,two /0.d0,0.5d0,1.d0,2.d0/
      if (interp .eq. 1) then
*
*  here for first pass through pot subroutine
*  define spline initial conditions
        ibcl = 1
        ibcr = 2
        fbcl = zero
        fbcr = zero
        kntopt = 2
        inbv=1
        dr=rr(2)-rr(1)
*  compute initial derivatives
        d0=(v0(2)-v0(1))/dr
        d2=(v2(2)-v2(1))/dr

*  interpolate the  coefficients for the energies
         call bint4 (rr, v0, nn, ibcl,ibcr,d0,fbcr,kntopt,
     :               t0, b0,nc0,kord, work)
         call bint4 (rr, v2, nn, ibcl,ibcr,d2,fbcr,kntopt,
     :               t2, b2,nc2,kord, work)
         interp=0
       endif
* here is  the entry point for subsequent entries into pot subroutine
      vvl(1)=zero
      vv0 = zero
* here to determine potential
      kord=4
      izero=0
      if (rz .gt. rr(nn-2)) then
*  here for inverse power extrapolation at large r
        r6=rz**6
        do 50 ii=1,2
          vx(ii)=clr(ii)/r6
50      continue
      else if (rz. gt. rr(3)) then
* here for spline range
        kord=4
        vx(1)= bvalu (t0, b0, nc0, kord,izero, rz, inb0, work)
        vx(2)= bvalu (t2, b2, nc2, kord,izero, rz, inb2, work)
* exponential extrapolation at short r
      else
        vx(1)=a(1)*exp(-alph(1)*rz)
        vx(2)=a(2)*exp(-alph(2)*rz)
      endif
      vv0=vx(1)
      call dscal(5,zero,vvl,1)
      vvl(3)=vx(2)
      return
      end
