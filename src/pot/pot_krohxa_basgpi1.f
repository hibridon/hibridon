*system:  test routine for OH(A/X)-Kr system
*         to debug basgpi1 basis routine
*
*   written by P.j.dagdigian (dec-2011)
*   modified by J. Klos (jan-2012)
*
*   divided V1 PES by sqrt(2).  see mha and corey, jcp 84, 100 (1986)
*
      subroutine driver
      implicit double precision (a-h,o-z)
      character *48 potnam
      logical csflag, ljunk, ihomo, lljunk
      include "common/parbas"
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(69)
      common /covpot/ numvib,ivibpi(5)
      potnam='OH(A/X)-Kr Quenching Klos PESs'
      print *, potnam
*  consider only v=0 and 1 vib levels of 2pi state
      numvib = 2
      ivibpi(1)=0
      ivibpi(2)=1
      nterm=4 + 3*(numvib - 1)
      write (6,89) numvib, (ivibpi(i),i=1,numvpi)
89    format(' Number of 2pi vibrational levels =',i3/
     :  5x,'v =',10i3)
      print *
1      print *, ' r (bohr) '
      read (5, *, end=93) r
      if (r .lt. 0.d0) go to 93
      call pot(vv0,r)
      write (6, 100) vvl
100   format(' vsig  ',11(pe16.8)/
     :  4(' vpisum'11(1pe16.8),/,
     :  ' vpi2  ',9(1pe16.8),/,
     :  ' vpi1  ',9(1pe16.8)/))
      goto 1
93    r=4

*      do i=1,100
*       call pot(vv0,r)
*       write(2,101) r,vvl
*101    format(f8.4,8(1pe16.8))
*       r=r+0.2
*      enddo

99    end
* ----------------------------------------------------------------------
      include "common/syusr"
      include "common/bausr"
      include "common/ground"
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      common /coselb/ ibasty
      common /covpot/ numvib,ivibpi(5)
      potnam='OH(A/X)-Kr Quenching Klos PESs'
*      ibasty=19
      numvib=2
*  consider only v=0 and 1 vib levels of 2pi state
      ivibpi(1)=0
      ivibpi(2)=1
      nterm=4 + 3*(numvib - 1)
*  for 2sigma state
      lammin(1)=0
      lammax(1)=10
      mproj(1)=0
*  for first 2pi level
      lammin(2)=0
      lammax(2)=10
      mproj(2)=0
      lammin(3)=2
      lammax(3)=10
      mproj(3)=2
      lammin(4)=1
      lammax(4)=9
      mproj(4)=1
*  for second 2pi level
      lammin(5)=0
      lammax(5)=10
      mproj(5)=0
      lammin(6)=2
      lammax(6)=10
      mproj(6)=2
      lammin(7)=1
      lammax(7)=9
      mproj(7)=1
      return
      end
* ----------------------------------------------------------------------
      subroutine pot (vv0, r)
* ----------------------------------------------------------------------
*  subroutine to calculate the r-dependent coefficients 
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*    vv0      dummy term here
*  variable in common block /covvl/
*    vvl:     vector to store r-dependence of each term
*             in potential expansion
*    vvl(1:x) expansion of vsigma in d(l,0) terms, for l=0,...
*    vvl(x+1:y) expansion of vsum for pi(v=0) in d(l,0) terms
*    vvl(y+1:z) expansion of vdif for pi(v=0) in d(l,2) terms
*    vvl(z+1:a) expansion of v1 coupling of sigma to pi(v=0) in d(l,1) terms
*    and more terms if numvib > 0
*
* author:  paul dagdigian
* latest revision date:  24-jan-2012
* latest revision by jacek klos 24-jan-2012
* ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical ljunk, ihomo, csflag, lljunk
      include "common/parbas"
      dimension vsig(11), vpi(11), v2(9),v1(10)
      dimension d0(121),d1(81),d2(81),aa(121)
      dimension xsig(11),xpi(11),x1(10),x2(9),kpvt(11),
     :          qraux(11), work(121),rsd(11)
      dimension thta(11)
      dimension FCF(2) ! stores Franck-Condon factors
      dimension vvll(40)
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(69)
      common /covpot/ numvib,ivibpi(5)

      data zero, one, half /0.d0,1.d0,0.5d0/
* hyperbolic tangent scaling factor
      data alph /1.2d0/
      data rmax /13d0/
      aa=zero
* Here are Franck-Condon factors (as squares of the <Av=0|Xv=v'>
* overlap integrals) between transitions from vibrational
* state v=0 of A state
* to vibrational states v' of the X state:
* FCF_0_0 = 0.917664288481494
* FCF_0_1 = 0.0800257653502253
* FCF_0_2 =0.00227074305148058
* FCF_0_3 =3.83580435449925e-05
*Obtained from MRCISD+Q/AVQZ diatomic X and A potentials.
*J. Klos 2011 November
       FCF(1)=0.917664288481494D0
       FCF(2)=0.0800257653502253D0
* below 2 lines for debug only
*       FCF(1)=1.D0
*       FCF(2)=1.D0

      nangle=11
      thta(1)=0.D0
      thta(2)=20.D0
      thta(3)=40.D0
      thta(4)=60.D0
      thta(5)=80.D0
      thta(6)=90.D0
      thta(7)=100.D0
      thta(8)=120.D0
      thta(9)=140.D0
      thta(10)=160.D0
      thta(11)=180.D0

      lmax=10
      index=0
      index1=0
      index2=0
      do l=0,lmax
      do i=1,nangle
      index=index+1
      d0(index)=DLM0(l,0,thta(i))
      if (l.ge.1.and.l.lt.lmax.and.i.gt.1.and.i.lt.nangle) then
      index1=index1+1
      d1(index1)=DLM0(l,1,thta(i))
c      write(6,*) index1,thta(i),d1(index1)
      endif
      if (l.ge.2.and.i.gt.1.and.i.lt.nangle) then
      index2=index2+1
      d2(index2)=DLM0(l,2,thta(i))
      endif
      enddo
      enddo
C    Setup values of Vsigma, VPi, V1 and V2 PESs at angular grid
      do i=1,nangle
        vsig(i)=V_KROH_A(r,thta(i))
C   Vpi and V2 are Vsum and Vdif of  A' and A" of the X state
        vpi(i)=Vsum_rre(r,thta(i))
        if (i.gt.1.and.i.lt.nangle) then
*  divide V1 by sqrt(2)
          v1(i-1)=H12_KROH(r,thta(i)) / sqrt(2.d0)
          v2(i-1)=Vdif_rre(r,thta(i))
        endif
        if (r .gt. rmax) then
          damp=-half*(tanh(3.d0*(r-rmax))-one)
          v2(i)=v2(i)*damp
          v1(i)=v1(i)*damp
        endif
      enddo
*
*     vv0 term is a dummy
      vv0 = 0.d0
* solve simultaneous equations for solutions
* first for vsig
      tol=1.e-10
      call dcopy(121,d0,1,aa,1)
      call dqrank(aa,11,11,11,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,11,11,11,kr,vsig,xsig,rsd,kpvt,qraux)
      call dcopy(11,xsig,1,vvll,1)
* then for vpi
      call dcopy(121,d0,1,aa,1)
      call dqrank(aa,11,11,11,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,11,11,11,kr,vpi,xpi,rsd,kpvt,qraux)
      call dcopy(11,xpi,1,vvll(12),1)
* then for v2
      call dcopy(81,d2,1,aa,1)
      call dqrank(aa,9,9,9,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,9,9,9,kr,v2,x2,rsd,kpvt,qraux)
      call dcopy(9,x2,1,vvll(23),1)
* then for v1
      call dcopy(81,d1,1,aa,1)
      call dqrank(aa,9,9,9,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,9,9,9,kr,v1,x1,rsd,kpvt,qraux)
      call dcopy(9,x1,1,vvll(32),1)
* convert to hartree
      econv=1./219474.6
      call dscal(40,econv,vvll,1)
      vv0=vv0*econv
      il=lmax+1
      do i=1,lmax+1
      vvl(i)=vvll(i)
      enddo
* scale only V1 term with FCF factors
      do i=1,numvib
      do ilam=lmax+2,40
      il=il+1
      if (ilam.lt.(40-(lmax-2))) then
       vvl(il)=vvll(ilam)
      else
       vvl(il)=FCF(i)*vvll(ilam)
      endif
      enddo
      enddo
      return
      end


C******************************************************

      FUNCTION PLGNDR(L,M,x)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.)pause
     *'bad arguments in plgndr'
      pmm=1.
      if(m.gt.0) then
        somx2=dsqrt((1.-x)*(1.+x))
        fact=1.
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.
11      continue
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
         do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndr=pll
        endif
      endif
      return
      END


      FUNCTION H11_KROH(RR,TT) 
C Kr-OH(X-A) H11 DIABAT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      TE=32803.75D0 
      X_APRIM=Vsum_rre(RR,TT)-Vdif_rre(RR,TT)
      A_APRIM=V_KROH_A(RR,TT)+TE
      SIN2A=SIN2ALFA(RR,TT)
      SINSQ=0.5D0*(1.D0-DSQRT(1.D0-SIN2A**2))
      COSSQ=0.5D0*(1.D0+DSQRT(1.D0-SIN2A**2))
      H11_KROH=X_APRIM*SINSQ+A_APRIM*COSSQ
      RETURN
      END   


      FUNCTION H22_KROH(RR,TT)
C Kr-OH(X-A) H22 DIABAT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      TE=32803.75D0
      X_APRIM=Vsum_rre(RR,TT)-Vdif_rre(RR,TT)
      A_APRIM=V_KROH_A(RR,TT)+TE
      SIN2A=SIN2ALFA(RR,TT)
      SINSQ=0.5D0*(1.D0-DSQRT(1.D0-SIN2A**2))
      COSSQ=0.5D0*(1.D0+DSQRT(1.D0-SIN2A**2))
      H22_KROH=X_APRIM*COSSQ+A_APRIM*SINSQ
      RETURN
      END

      FUNCTION H12_KROH(RR,TT)
C Kr-OH(X-A) H12 DIABAT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      TE=32803.75D0
       TE=32440.6D0 ! Value of the asymptotic E of the A state from Paul
C      TE=0.D0 ! Not sure  if sHift is zero for Hibridon Pot
      X_APRIM=Vsum_rre(RR,TT)-Vdif_rre(RR,TT)
      A_APRIM=V_KROH_A(RR,TT)+TE
      SIN2A=SIN2ALFA(RR,TT)
      H12_KROH=(A_APRIM-X_APRIM)*0.5D0*SIN2A
      RETURN
      END


      FUNCTION X1APRIM_KROH(RR,TT)
C Kr-OH(X) 1A' ADIABAT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      X1APRIM_KROH=Vsum_rre(RR,TT)-Vdif_rre(RR,TT)
      RETURN
      END

      FUNCTION A2APRIM_KROH(RR,TT)
C Kr-OH(A) 2A'  ADIABAT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      TE=32803.75D0
      A2APRIM_KROH=V_KROH_A(RR,TT)+TE
      RETURN
      END

      FUNCTION X1ABIS_KROH(RR,TT)
C Kr-OH(X) 1A'' ADIABAT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      X1ABIS_KROH=Vsum_rre(RR,TT)+Vdif_rre(RR,TT)
      RETURN
      END










C PROGRAM TO GENERATE 2-D Kr-OH(r=re)(A^2Sigma) 
C 1A'-2A' sin(2*alfa) WHERE alfa IS MIXING ANGLE
C BETWEEN X(1A') and A(2A') ADIABATS OBTAINED
C FROM MRCISD/AVTZ-DK CALCULATIONS
C Units: input in (bohrs,degrees)
C Reference: J. Klos, to be published
      FUNCTION SIN2ALFA(RR,TT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=962)
      DIMENSION VAB(NPT),VCOEF(NPT),VCOEF1(481),VCOEF2(481)
      DIMENSION R1(481),R2(481),THETA1(481),THETA2(481) 
      DIMENSION R(NPT),THETA(NPT)
      DATA R1/
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.2D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .              3.3D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .             3.35D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .              3.4D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .             3.45D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .              3.5D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .             3.55D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .              3.6D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .             3.85D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .              3.9D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .             3.95D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .                4D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0,
     .              4.1D0/
      DATA R2/
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .                4.2D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .               4.25D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.3D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.4D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .                4.5D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .               4.75D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                4.9D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .                5.1D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .               5.25D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.4D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.5D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.6D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0,
     .                5.7D0/
      DATA THETA1/
     .                 0D0,
     .               0.5D0,
     .                 1D0,
     .                 2D0,
     .                 3D0,
     .                 4D0,
     .                 5D0,
     .                10D0,
     .                15D0,
     .                20D0,
     .                25D0,
     .                30D0,
     .                35D0,
     .                40D0,
     .                45D0,
     .                50D0,
     .                60D0,
     .                70D0,
     .                80D0,
     .                90D0,
     .               100D0,
     .               110D0,
     .               120D0,
     .               130D0,
     .               135D0,
     .               150D0,
     .               155D0,
     .               160D0,
     .               165D0,
     .               170D0,
     .               175D0,
     .               176D0,
     .               177D0,
     .               178D0,
     .               179D0,
     .             179.5D0,
     .               180D0,
     .                 0D0,
     .               0.5D0,
     .                 1D0,
     .                 2D0,
     .                 3D0,
     .                 4D0,
     .                 5D0,
     .                10D0,
     .                15D0,
     .                20D0,
     .                25D0,
     .                30D0,
     .                35D0,
     .                40D0,
     .                45D0,
     .                50D0,
     .                60D0,
     .                70D0,
     .                80D0,
     .                90D0,
     .               100D0,
     .               110D0,
     .               120D0,
     .               130D0,
     .               135D0,
     .               150D0,
     .               155D0,
     .               160D0,
     .               165D0,
     .               170D0,
     .               175D0,
     .               176D0,
     .               177D0,
     .               178D0,
     .               179D0,
     .             179.5D0,
     .               180D0,
     .                 0D0,
     .               0.5D0,
     .                 1D0,
     .                 2D0,
     .                 3D0,
     .                 4D0,
     .                 5D0,
     .                10D0,
     .                15D0,
     .                20D0,
     .                25D0,
     .                30D0,
     .                35D0,
     .                40D0,
     .                45D0,
     .                50D0,
     .                60D0,
     .                70D0,
     .                80D0,
     .                90D0,
     .               100D0,
     .               110D0,
     .               120D0,
     .               130D0,
     .               135D0,
     .               150D0,
     .               155D0,
     .               160D0,
     .               165D0,
     .               170D0,
     .               175D0,
     .               176D0,
     .               177D0,
     .               178D0,
     .               179D0,
     .             179.5D0,
     .               180D0,
     .                 0D0,
     .               0.5D0,
     .                 1D0,
     .                 2D0,
     .                 3D0,
     .                 4D0,
     .                 5D0,
     .                10D0,
     .                15D0,
     .                20D0,
     .                25D0,
     .                30D0,
     .                35D0,
     .                40D0,
     .                45D0,
     .                50D0,
     .                60D0,
     .                70D0,
     .                80D0,
     .                90D0,
     .               100D0,
     .               110D0,
     .               120D0,
     .               130D0,
     .               135D0,
     .               150D0,
     .               155D0,
     .               160D0,
     .               165D0,
     .               170D0,
     .               175D0,
     .               176D0,
     .               177D0,
     .               178D0,
     .               179D0,
     .             179.5D0,
     .               180D0,
     .                 0D0,
     .               0.5D0,
     .                 1D0,
     .                 2D0,
     .                 3D0,
     .                 4D0,
     .                 5D0,
     .                10D0,
     .                15D0,
     .                20D0,
     .                25D0,
     .                30D0,
     .                35D0,
     .                40D0,
     .                45D0,
     .                50D0,
     .                60D0,
     .                70D0,
     .                80D0,
     .                90D0,
     .               100D0,
     .               110D0,
     .               120D0,
     .               130D0,
     .               135D0,
     .               150D0,
     .               155D0,
     .               160D0,
     .               165D0,
     .               170D0,
     .               175D0,
     .               176D0,
     .               177D0,
     .               178D0,
     .               179D0,
     .             179.5D0,
     .               180D0,
     .                 0D0,
     .               0.5D0,
     .                 1D0,
     .                 2D0,
     .                 3D0,
     .                 4D0,
     .                 5D0,
     .                10D0,
     .                15D0,
     .                20D0,
     .                25D0,
     .                30D0,
     .                35D0,
     .                40D0,
     .                45D0,
     .                50D0,
     .                60D0,
     .                70D0,
     .                80D0,
     .                90D0,
     .               100D0,
     .               110D0,
     .               120D0,
     .               130D0,
     .               135D0,
     .               150D0,
     .               155D0,
     .               160D0,
     .               165D0,
     .               170D0,
     .               175D0,
     .               176D0,
     .               177D0,
     .               178D0,
     .               179D0,
     .             179.5D0,
     .               180D0,
     .                 0D0,
     .               0.5D0,
     .                 1D0,
     .                 2D0,
     .                 3D0,
     .                 4D0,
     .                 5D0,
     .                10D0,
     .                15D0,
     .                20D0,
     .                25D0,
     .                30D0,
     .                35D0,
     .                40D0,
     .                45D0,
     .                50D0,
     .                60D0,
     .                70D0,
     .                80D0,
     .                90D0,
     .               100D0,
     .               110D0,
     .               120D0,
     .               130D0,
     .               135D0,
     .               150D0,
     .               155D0,
     .               160D0,
     .               165D0,
     .               170D0,
     .               175D0,
     .               176D0,
     .               177D0,
     .               178D0,
     .               179D0,
     .             179.5D0,
     .               180D0,
     .                 0D0,
     .               0.5D0,
     .                 1D0,
     .                 2D0,
     .                 3D0,
     .                 4D0,
     .                 5D0,
     .                10D0,
     .                15D0,
     .                20D0,
     .                25D0,
     .                30D0,
     .                35D0,
     .                40D0,
     .                45D0,
     .                50D0,
     .                60D0,
     .                70D0,
     .                80D0,
     .                90D0,
     .               100D0,
     .               110D0,
     .               120D0,
     .               130D0,
     .               135D0,
     .               150D0,
     .               155D0,
     .               160D0,
     .               165D0,
     .               170D0,
     .               175D0,
     .               176D0,
     .               177D0,
     .               178D0,
     .               179D0,
     .             179.5D0,
     .               180D0,
     .                 0D0,
     .               0.5D0,
     .                 1D0,
     .                 2D0,
     .                 3D0,
     .                 4D0,
     .                 5D0,
     .                10D0,
     .                15D0,
     .                20D0,
     .                25D0,
     .                30D0,
     .                35D0,
     .                40D0,
     .                45D0,
     .                50D0,
     .                60D0,
     .                70D0,
     .                80D0,
     .                90D0,
     .               100D0,
     .               110D0,
     .               120D0,
     .               130D0,
     .               135D0,
     .               150D0,
     .               155D0,
     .               160D0,
     .               165D0,
     .               170D0,
     .               175D0,
     .               176D0,
     .               177D0,
     .               178D0,
     .               179D0,
     .             179.5D0,
     .               180D0,
     .                 0D0,
     .               0.5D0,
     .                 1D0,
     .                 2D0,
     .                 3D0,
     .                 4D0,
     .                 5D0,
     .                10D0,
     .                15D0,
     .                20D0,
     .                25D0,
     .                30D0,
     .                35D0,
     .                40D0,
     .                45D0,
     .                50D0,
     .                60D0,
     .                70D0,
     .                80D0,
     .                90D0,
     .               100D0,
     .               110D0,
     .               120D0,
     .               130D0,
     .               135D0,
     .               150D0,
     .               155D0,
     .               160D0,
     .               165D0,
     .               170D0,
     .               175D0,
     .               176D0,
     .               177D0,
     .               178D0,
     .               179D0,
     .             179.5D0,
     .               180D0,
     .                 0D0,
     .               0.5D0,
     .                 1D0,
     .                 2D0,
     .                 3D0,
     .                 4D0,
     .                 5D0,
     .                10D0,
     .                15D0,
     .                20D0,
     .                25D0,
     .                30D0,
     .                35D0,
     .                40D0,
     .                45D0,
     .                50D0,
     .                60D0,
     .                70D0,
     .                80D0,
     .                90D0,
     .               100D0,
     .               110D0,
     .               120D0,
     .               130D0,
     .               135D0,
     .               150D0,
     .               155D0,
     .               160D0,
     .               165D0,
     .               170D0,
     .               175D0,
     .               176D0,
     .               177D0,
     .               178D0,
     .               179D0,
     .             179.5D0,
     .               180D0,
     .                 0D0,
     .               0.5D0,
     .                 1D0,
     .                 2D0,
     .                 3D0,
     .                 4D0,
     .                 5D0,
     .                10D0,
     .                15D0,
     .                20D0,
     .                25D0,
     .                30D0,
     .                35D0,
     .                40D0,
     .                45D0,
     .                50D0,
     .                60D0,
     .                70D0,
     .                80D0,
     .                90D0,
     .               100D0,
     .               110D0,
     .               120D0,
     .               130D0,
     .               135D0,
     .               150D0,
     .               155D0,
     .               160D0,
     .               165D0,
     .               170D0,
     .               175D0,
     .               176D0,
     .               177D0,
     .               178D0,
     .               179D0,
     .             179.5D0,
     .               180D0,
     .                 0D0,
     .               0.5D0,
     .                 1D0,
     .                 2D0,
     .                 3D0,
     .                 4D0,
     .                 5D0,
     .                10D0,
     .                15D0,
     .                20D0,
     .                25D0,
     .                30D0,
     .                35D0,
     .                40D0,
     .                45D0,
     .                50D0,
     .                60D0,
     .                70D0,
     .                80D0,
     .                90D0,
     .               100D0,
     .               110D0,
     .               120D0,
     .               130D0,
     .               135D0,
     .               150D0,
     .               155D0,
     .               160D0,
     .               165D0,
     .               170D0,
     .               175D0,
     .               176D0,
     .               177D0,
     .               178D0,
     .               179D0,
     .             179.5D0,
     .               180D0/
      DATA THETA2 /
     .                  0D0,
     .                0.5D0,
     .                  1D0,
     .                  2D0,
     .                  3D0,
     .                  4D0,
     .                  5D0,
     .                 10D0,
     .                 15D0,
     .                 20D0,
     .                 25D0,
     .                 30D0,
     .                 35D0,
     .                 40D0,
     .                 45D0,
     .                 50D0,
     .                 60D0,
     .                 70D0,
     .                 80D0,
     .                 90D0,
     .                100D0,
     .                110D0,
     .                120D0,
     .                130D0,
     .                135D0,
     .                150D0,
     .                155D0,
     .                160D0,
     .                165D0,
     .                170D0,
     .                175D0,
     .                176D0,
     .                177D0,
     .                178D0,
     .                179D0,
     .              179.5D0,
     .                180D0,
     .                  0D0,
     .                0.5D0,
     .                  1D0,
     .                  2D0,
     .                  3D0,
     .                  4D0,
     .                  5D0,
     .                 10D0,
     .                 15D0,
     .                 20D0,
     .                 25D0,
     .                 30D0,
     .                 35D0,
     .                 40D0,
     .                 45D0,
     .                 50D0,
     .                 60D0,
     .                 70D0,
     .                 80D0,
     .                 90D0,
     .                100D0,
     .                110D0,
     .                120D0,
     .                130D0,
     .                135D0,
     .                150D0,
     .                155D0,
     .                160D0,
     .                165D0,
     .                170D0,
     .                175D0,
     .                176D0,
     .                177D0,
     .                178D0,
     .                179D0,
     .              179.5D0,
     .                180D0,
     .                  0D0,
     .                0.5D0,
     .                  1D0,
     .                  2D0,
     .                  3D0,
     .                  4D0,
     .                  5D0,
     .                 10D0,
     .                 15D0,
     .                 20D0,
     .                 25D0,
     .                 30D0,
     .                 35D0,
     .                 40D0,
     .                 45D0,
     .                 50D0,
     .                 60D0,
     .                 70D0,
     .                 80D0,
     .                 90D0,
     .                100D0,
     .                110D0,
     .                120D0,
     .                130D0,
     .                135D0,
     .                150D0,
     .                155D0,
     .                160D0,
     .                165D0,
     .                170D0,
     .                175D0,
     .                176D0,
     .                177D0,
     .                178D0,
     .                179D0,
     .              179.5D0,
     .                180D0,
     .                  0D0,
     .                0.5D0,
     .                  1D0,
     .                  2D0,
     .                  3D0,
     .                  4D0,
     .                  5D0,
     .                 10D0,
     .                 15D0,
     .                 20D0,
     .                 25D0,
     .                 30D0,
     .                 35D0,
     .                 40D0,
     .                 45D0,
     .                 50D0,
     .                 60D0,
     .                 70D0,
     .                 80D0,
     .                 90D0,
     .                100D0,
     .                110D0,
     .                120D0,
     .                130D0,
     .                135D0,
     .                150D0,
     .                155D0,
     .                160D0,
     .                165D0,
     .                170D0,
     .                175D0,
     .                176D0,
     .                177D0,
     .                178D0,
     .                179D0,
     .              179.5D0,
     .                180D0,
     .                  0D0,
     .                0.5D0,
     .                  1D0,
     .                  2D0,
     .                  3D0,
     .                  4D0,
     .                  5D0,
     .                 10D0,
     .                 15D0,
     .                 20D0,
     .                 25D0,
     .                 30D0,
     .                 35D0,
     .                 40D0,
     .                 45D0,
     .                 50D0,
     .                 60D0,
     .                 70D0,
     .                 80D0,
     .                 90D0,
     .                100D0,
     .                110D0,
     .                120D0,
     .                130D0,
     .                135D0,
     .                150D0,
     .                155D0,
     .                160D0,
     .                165D0,
     .                170D0,
     .                175D0,
     .                176D0,
     .                177D0,
     .                178D0,
     .                179D0,
     .              179.5D0,
     .                180D0,
     .                  0D0,
     .                0.5D0,
     .                  1D0,
     .                  2D0,
     .                  3D0,
     .                  4D0,
     .                  5D0,
     .                 10D0,
     .                 15D0,
     .                 20D0,
     .                 25D0,
     .                 30D0,
     .                 35D0,
     .                 40D0,
     .                 45D0,
     .                 50D0,
     .                 60D0,
     .                 70D0,
     .                 80D0,
     .                 90D0,
     .                100D0,
     .                110D0,
     .                120D0,
     .                130D0,
     .                135D0,
     .                150D0,
     .                155D0,
     .                160D0,
     .                165D0,
     .                170D0,
     .                175D0,
     .                176D0,
     .                177D0,
     .                178D0,
     .                179D0,
     .              179.5D0,
     .                180D0,
     .                  0D0,
     .                0.5D0,
     .                  1D0,
     .                  2D0,
     .                  3D0,
     .                  4D0,
     .                  5D0,
     .                 10D0,
     .                 15D0,
     .                 20D0,
     .                 25D0,
     .                 30D0,
     .                 35D0,
     .                 40D0,
     .                 45D0,
     .                 50D0,
     .                 60D0,
     .                 70D0,
     .                 80D0,
     .                 90D0,
     .                100D0,
     .                110D0,
     .                120D0,
     .                130D0,
     .                135D0,
     .                150D0,
     .                155D0,
     .                160D0,
     .                165D0,
     .                170D0,
     .                175D0,
     .                176D0,
     .                177D0,
     .                178D0,
     .                179D0,
     .              179.5D0,
     .                180D0,
     .                  0D0,
     .                0.5D0,
     .                  1D0,
     .                  2D0,
     .                  3D0,
     .                  4D0,
     .                  5D0,
     .                 10D0,
     .                 15D0,
     .                 20D0,
     .                 25D0,
     .                 30D0,
     .                 35D0,
     .                 40D0,
     .                 45D0,
     .                 50D0,
     .                 60D0,
     .                 70D0,
     .                 80D0,
     .                 90D0,
     .                100D0,
     .                110D0,
     .                120D0,
     .                130D0,
     .                135D0,
     .                150D0,
     .                155D0,
     .                160D0,
     .                165D0,
     .                170D0,
     .                175D0,
     .                176D0,
     .                177D0,
     .                178D0,
     .                179D0,
     .              179.5D0,
     .                180D0,
     .                  0D0,
     .                0.5D0,
     .                  1D0,
     .                  2D0,
     .                  3D0,
     .                  4D0,
     .                  5D0,
     .                 10D0,
     .                 15D0,
     .                 20D0,
     .                 25D0,
     .                 30D0,
     .                 35D0,
     .                 40D0,
     .                 45D0,
     .                 50D0,
     .                 60D0,
     .                 70D0,
     .                 80D0,
     .                 90D0,
     .                100D0,
     .                110D0,
     .                120D0,
     .                130D0,
     .                135D0,
     .                150D0,
     .                155D0,
     .                160D0,
     .                165D0,
     .                170D0,
     .                175D0,
     .                176D0,
     .                177D0,
     .                178D0,
     .                179D0,
     .              179.5D0,
     .                180D0,
     .                  0D0,
     .                0.5D0,
     .                  1D0,
     .                  2D0,
     .                  3D0,
     .                  4D0,
     .                  5D0,
     .                 10D0,
     .                 15D0,
     .                 20D0,
     .                 25D0,
     .                 30D0,
     .                 35D0,
     .                 40D0,
     .                 45D0,
     .                 50D0,
     .                 60D0,
     .                 70D0,
     .                 80D0,
     .                 90D0,
     .                100D0,
     .                110D0,
     .                120D0,
     .                130D0,
     .                135D0,
     .                150D0,
     .                155D0,
     .                160D0,
     .                165D0,
     .                170D0,
     .                175D0,
     .                176D0,
     .                177D0,
     .                178D0,
     .                179D0,
     .              179.5D0,
     .                180D0,
     .                  0D0,
     .                0.5D0,
     .                  1D0,
     .                  2D0,
     .                  3D0,
     .                  4D0,
     .                  5D0,
     .                 10D0,
     .                 15D0,
     .                 20D0,
     .                 25D0,
     .                 30D0,
     .                 35D0,
     .                 40D0,
     .                 45D0,
     .                 50D0,
     .                 60D0,
     .                 70D0,
     .                 80D0,
     .                 90D0,
     .                100D0,
     .                110D0,
     .                120D0,
     .                130D0,
     .                135D0,
     .                150D0,
     .                155D0,
     .                160D0,
     .                165D0,
     .                170D0,
     .                175D0,
     .                176D0,
     .                177D0,
     .                178D0,
     .                179D0,
     .              179.5D0,
     .                180D0,
     .                  0D0,
     .                0.5D0,
     .                  1D0,
     .                  2D0,
     .                  3D0,
     .                  4D0,
     .                  5D0,
     .                 10D0,
     .                 15D0,
     .                 20D0,
     .                 25D0,
     .                 30D0,
     .                 35D0,
     .                 40D0,
     .                 45D0,
     .                 50D0,
     .                 60D0,
     .                 70D0,
     .                 80D0,
     .                 90D0,
     .                100D0,
     .                110D0,
     .                120D0,
     .                130D0,
     .                135D0,
     .                150D0,
     .                155D0,
     .                160D0,
     .                165D0,
     .                170D0,
     .                175D0,
     .                176D0,
     .                177D0,
     .                178D0,
     .                179D0,
     .              179.5D0,
     .                180D0,
     .                  0D0,
     .                0.5D0,
     .                  1D0,
     .                  2D0,
     .                  3D0,
     .                  4D0,
     .                  5D0,
     .                 10D0,
     .                 15D0,
     .                 20D0,
     .                 25D0,
     .                 30D0,
     .                 35D0,
     .                 40D0,
     .                 45D0,
     .                 50D0,
     .                 60D0,
     .                 70D0,
     .                 80D0,
     .                 90D0,
     .                100D0,
     .                110D0,
     .                120D0,
     .                130D0,
     .                135D0,
     .                150D0,
     .                155D0,
     .                160D0,
     .                165D0,
     .                170D0,
     .                175D0,
     .                176D0,
     .                177D0,
     .                178D0,
     .                179D0,
     .              179.5D0,
     .                180D0/
      DATA VCOEF1/
     .                 -3176894078705.7753906250000000D0,
     .                  -483862676545.8272705078125000D0,
     .                  1476117171348.3085937500000000D0,
     .                -12039920758258.9492187500000000D0,
     .                 24337902329084.2773437500000000D0,
     .                 -6966364838079.2832031250000000D0,
     .                 -3465007159330.9560546875000000D0,
     .                   375907538596.8083496093750000D0,
     .                   -66863369888.2847290039062500D0,
     .                    10616795899.9810409545898438D0,
     .                    -2329806645.5189266204833984D0,
     .                      867870000.2212023735046387D0,
     .                     -169375829.6556873321533203D0,
     .                      -53183862.3686604499816895D0,
     .                      129885050.5654751062393188D0,
     .                      -58731887.4158499538898468D0,
     .                     -163161165.0691993534564972D0,
     .                      374156790.9578830003738403D0,
     .                     -412366856.3683130741119385D0,
     .                      272152146.7076357603073120D0,
     .                     -128793981.3090383559465408D0,
     .                       72592856.6273121982812881D0,
     .                      -65502185.9983831942081451D0,
     .                       76928940.4729885458946228D0,
     .                      -35767002.3828404545783997D0,
     .                    -2004927237.0508267879486084D0,
     .                    10327804024.3898162841796875D0,
     .                   -36951085191.3492584228515625D0,
     .                   170227477842.6778259277343750D0,
     .                 -1008160296057.0817871093750000D0,
     .                 17206844852012.0937500000000000D0,
     .                 -9632880130037.3867187500000000D0,
     .                -32790365124500.7304687500000000D0,
     .                 22278002376721.6875000000000000D0,
     .                  6153517603679.8603515625000000D0,
     .                -14177917794420.3437500000000000D0,
     .                 11829351393284.7167968750000000D0,
     .                 -5129768011192.8212890625000000D0,
     .                 22110137563897.0898437500000000D0,
     .                 -8723298314513.0507812500000000D0,
     .                -20764974001038.7929687500000000D0,
     .                  7064767710095.9902343750000000D0,
     .                  5834263297418.5986328125000000D0,
     .                   734663471390.2509765625000000D0,
     .                 -1371796674315.8452148437500000D0,
     .                   294316672388.6135253906250000D0,
     .                   -60202146337.7350463867187500D0,
     .                    18170478381.3372192382812500D0,
     .                    -8419583143.8525733947753906D0,
     .                     2103794280.9416284561157227D0,
     .                      194412757.5608372688293457D0,
     .                     -205919778.2813071012496948D0,
     .                     -144804757.9165096282958984D0,
     .                      827563007.4411693811416626D0,
     .                    -1610551307.2103805541992187D0,
     .                     1818457878.9050478935241699D0,
     .                    -1293307165.8772618770599365D0,
     .                      666591907.4666564464569092D0,
     .                     -372767697.4322705268859863D0,
     .                      299161483.2983213067054749D0,
     .                     -286193791.0897908210754395D0,
     .                      -30093540.5606286525726318D0,
     .                    10920919964.9936027526855469D0,
     .                   -53822531854.4219970703125000D0,
     .                   171668393043.8169860839843750D0,
     .                  -631362875337.5561523437500000D0,
     .                  2957467327560.4907226562500000D0,
     .                -15788831935965.6914062500000000D0,
     .                  5298272642898.7226562500000000D0,
     .                  5181522712125.8916015625000000D0,
     .                 17844764138634.1562500000000000D0,
     .                 -1677624377583.1958007812500000D0,
     .                -17898628923077.6328125000000000D0,
     .                  4585827696308.9453125000000000D0,
     .                 -2690851815734.5854492187500000D0,
     .                 21336133295334.8554687500000000D0,
     .                 -2976620499691.9399414062500000D0,
     .                -11369098065853.6542968750000000D0,
     .                 -3093044978062.3964843750000000D0,
     .                 -1371684860548.7653808593750000D0,
     .                 -1510642240695.7519531250000000D0,
     .                  2077855590042.0869140625000000D0,
     .                  -495722492946.6790771484375000D0,
     .                   121110102511.4160766601562500D0,
     .                   -41498838064.4053649902343750D0,
     .                    18451508041.0705413818359375D0,
     .                    -3685352504.1193442344665527D0,
     .                    -1136954509.2055070400238037D0,
     .                      118994315.0893921852111816D0,
     .                      583052203.6283757686614990D0,
     .                     -817644557.5357816219329834D0,
     .                     1432948948.9815449714660645D0,
     .                    -1712618524.4018895626068115D0,
     .                     1327852010.8455779552459717D0,
     .                     -774195625.2654943466186523D0,
     .                      463546924.8205135464668274D0,
     .                     -265861433.6653360128402710D0,
     .                      -80632406.6090581417083740D0,
     .                      586991977.0597629547119141D0,
     .                   -19045410214.0929794311523438D0,
     .                    96651992730.7175903320312500D0,
     .                  -309383861852.1696777343750000D0,
     .                   980685556769.0246582031250000D0,
     .                 -3644595933561.3217773437500000D0,
     .                 -6356724675306.2402343750000000D0,
     .                 14834550944027.4726562500000000D0,
     .                 12266853607041.7246093750000000D0,
     .                 11007701811013.1250000000000000D0,
     .                 -3452348473305.2822265625000000D0,
     .                -24802923995787.5117187500000000D0,
     .                  -601848386151.5905761718750000D0,
     .                 -8104752457803.8535156250000000D0,
     .                 21924535067660.1875000000000000D0,
     .                 -3867470092578.5844726562500000D0,
     .                 -8792596427710.9746093750000000D0,
     .                 -7807597629086.9277343750000000D0,
     .                  1382851115454.0717773437500000D0,
     .                  6153214686605.7910156250000000D0,
     .                 -1079125392742.9331054687500000D0,
     .                   252194535644.0242919921875000D0,
     .                   -88613760577.6320495605468750D0,
     .                    41320960748.1824340820312500D0,
     .                   -17247064020.1139144897460938D0,
     .                      594337893.3306903839111328D0,
     .                     4605243258.8116407394409180D0,
     .                    -2159862370.6984310150146484D0,
     .                      398090625.1802301406860352D0,
     .                     -425915109.5084649920463562D0,
     .                      464447627.5434807538986206D0,
     .                     -229438210.7527505159378052D0,
     .                       -6864643.9075227454304695D0,
     .                      139739086.7293491363525391D0,
     .                     -149383361.4026186466217041D0,
     .                     -175266754.9316487312316895D0,
     .                     1150557485.4510178565979004D0,
     .                    -1310964716.6948556900024414D0,
     .                    15040207735.6601943969726563D0,
     .                   -88425552661.2119445800781250D0,
     .                   290566349558.0510253906250000D0,
     .                  -761552158398.6149902343750000D0,
     .                  2463187487423.8417968750000000D0,
     .                -13472938514502.6601562500000000D0,
     .                  4650739073516.0791015625000000D0,
     .                  9591097563969.3554687500000000D0,
     .                  6928643423648.4589843750000000D0,
     .                  5092614770388.7744140625000000D0,
     .                -22029285355124.4570312500000000D0,
     .                  7320704467646.1621093750000000D0,
     .                -19018224719560.4453125000000000D0,
     .                 20736545179648.5546875000000000D0,
     .                    45575760731.0229492187500000D0,
     .                  1404220821848.7114257812500000D0,
     .                -13162237464860.0976562500000000D0,
     .                  2107617349082.8759765625000000D0,
     .                 10175855897573.5605468750000000D0,
     .                 -2771629933830.9096679687500000D0,
     .                   549688371437.5336914062500000D0,
     .                   -64271910589.5263061523437500D0,
     .                   -10660858337.0392761230468750D0,
     .                     7890504376.9015808105468750D0,
     .                     5745453137.1455249786376953D0,
     .                   -11920291302.2521953582763672D0,
     .                     8689775174.9033164978027344D0,
     .                    -3656128093.0107297897338867D0,
     .                     1473717468.4635763168334961D0,
     .                    -1290400509.8141212463378906D0,
     .                      945890792.5267701148986816D0,
     .                     -522005244.1226327419281006D0,
     .                      203549022.9501446485519409D0,
     .                      -88670012.6921858787536621D0,
     .                      542725759.0758690834045410D0,
     .                    -1816917383.1172590255737305D0,
     .                     1191772302.4877214431762695D0,
     .                     3498171912.1315460205078125D0,
     .                    -2776715893.2432250976562500D0,
     .                    72550169930.6799163818359375D0,
     .                  -932396851095.2076416015625000D0,
     .                  5962817049809.6738281250000000D0,
     .                -14927914347124.6230468750000000D0,
     .                 -8884046996150.8710937500000000D0,
     .                 11180881972513.0214843750000000D0,
     .                  -405980619065.3969726562500000D0,
     .                  1928526402505.0063476562500000D0,
     .                 -6863626668719.6435546875000000D0,
     .                 12868600912294.6699218750000000D0,
     .                -17959471585023.8632812500000000D0,
     .                 19420859804988.0156250000000000D0,
     .                  1581883428820.3286132812500000D0,
     .                  9049213953957.0410156250000000D0,
     .                -10893870543022.5000000000000000D0,
     .                 -1433496875858.2958984375000000D0,
     .                 -1317729124992.6665039062500000D0,
     .                  2028994265025.1738281250000000D0,
     .                  -570901964186.4980468750000000D0,
     .                    95932064560.5099639892578125D0,
     .                     5106543187.0550727844238281D0,
     .                    -5450010221.9946842193603516D0,
     .                   -10536122427.1817951202392578D0,
     .                    19881900975.5229644775390625D0,
     .                   -16234566235.3094406127929688D0,
     .                     6896470531.5141143798828125D0,
     .                    -1632446207.1983919143676758D0,
     .                      902508669.3792951107025146D0,
     .                     -544540144.1505579948425293D0,
     .                      330520856.1528505682945251D0,
     .                     -228399374.9624173939228058D0,
     .                      276110286.8139804005622864D0,
     .                     -944686531.8060768842697144D0,
     .                     2471856886.4207267761230469D0,
     .                    -1146718331.8682775497436523D0,
     .                   -21632212825.9034881591796875D0,
     .                    95676277197.9009094238281250D0,
     .                  -382415281153.8804931640625000D0,
     .                  1739133900890.4111328125000000D0,
     .                 -8136834555973.8300781250000000D0,
     .                 19409869117386.7421875000000000D0,
     .                  6933273095136.4101562500000000D0,
     .                  1435919612078.4326171875000000D0,
     .                -17937841064113.7500000000000000D0,
     .                  -775707071984.9682617187500000D0,
     .                -13318410826441.7910156250000000D0,
     .                 10958407265382.7832031250000000D0,
     .                -11353579960275.6796875000000000D0,
     .                 14145982889452.9531250000000000D0,
     .                  1331802841423.3686523437500000D0,
     .                 15119351812460.5000000000000000D0,
     .                 -9395053673271.6464843750000000D0,
     .                 -4704215318618.1005859375000000D0,
     .                 -8433591123373.0839843750000000D0,
     .                  3475854496457.5585937500000000D0,
     .                  -190926569912.8056640625000000D0,
     .                    24458589112.5879516601562500D0,
     .                   -28096802570.5057067871093750D0,
     .                     7688376234.3148689270019531D0,
     .                     8705294530.2878131866455078D0,
     .                   -18331004154.9448051452636719D0,
     .                    15497353969.7095222473144531D0,
     .                    -6342287583.2994184494018555D0,
     .                     1003916652.7681008577346802D0,
     .                     -254824383.5225356221199036D0,
     .                       72604332.2462242841720581D0,
     .                      -99127128.4176373481750488D0,
     .                      175167972.1380694210529327D0,
     .                     -374329637.6171214580535889D0,
     .                     1227771272.0728406906127930D0,
     .                    -3316308473.2840037345886230D0,
     .                     2165351526.9557838439941406D0,
     .                    18960303598.3812484741210938D0,
     .                   -65126082361.4609069824218750D0,
     .                   131334856895.5074462890625000D0,
     .                   463538719465.5562744140625000D0,
     .                 -7117609208656.5957031250000000D0,
     .                 24665704430019.3671875000000000D0,
     .                  4134892936030.3032226562500000D0,
     .                     6707212146.5260620117187500D0,
     .                -24497749252015.5312500000000000D0,
     .                 -3396233692419.5615234375000000D0,
     .                -10931576018964.3964843750000000D0,
     .                 16587350521153.2304687500000000D0,
     .                -28128361078897.8085937500000000D0,
     .                  5623150655412.3769531250000000D0,
     .                 -3756526785424.7568359375000000D0,
     .                 21073620515914.9492187500000000D0,
     .                  1285325479283.9355468750000000D0,
     .                  6875856671231.8701171875000000D0,
     .                 -1130838111437.3691406250000000D0,
     .                 -1927021918946.1982421875000000D0,
     .                    98021276417.9443359375000000D0,
     .                   -25638842623.7084960937500000D0,
     .                    15579892375.9526290893554688D0,
     .                    -3510071114.6323623657226562D0,
     .                    -2589344753.7663478851318359D0,
     .                     6754196139.5595245361328125D0,
     .                    -5961946779.0138931274414062D0,
     .                     2395204500.6470098495483398D0,
     .                     -279353927.6462996006011963D0,
     .                      -20134827.2151447534561157D0,
     .                       75302062.6168452501296997D0,
     .                      -23730077.0261340141296387D0,
     .                      -36143214.2646696120500565D0,
     .                      166192651.2817730903625488D0,
     .                     -691153180.8949660062789917D0,
     .                     2100999645.0239455699920654D0,
     .                    -1680892931.9520685672760010D0,
     .                    -6215650560.2353248596191406D0,
     .                     6800026282.1806945800781250D0,
     .                    79341269710.9949798583984375D0,
     .                 -1180754416293.3886718750000000D0,
     .                  9478760905555.6796875000000000D0,
     .                 -7140778658078.1611328125000000D0,
     .                -22892240115566.2343750000000000D0,
     .                -11362670624894.4960937500000000D0,
     .                -19449449070501.6679687500000000D0,
     .                 11853546306889.4218750000000000D0,
     .                  5042024152844.0146484375000000D0,
     .                 35571769027340.3125000000000000D0,
     .                -39529801314211.2812500000000000D0,
     .                 -9918441830814.1640625000000000D0,
     .                -14384873200195.5000000000000000D0,
     .                 24553773009508.4179687500000000D0,
     .                 20220118965802.3984375000000000D0,
     .                 21001288234964.0937500000000000D0,
     .                  7202940680549.8242187500000000D0,
     .                -10336178341945.7187500000000000D0,
     .                  1260348670316.7395019531250000D0,
     .                   -85648209547.9070739746093750D0,
     .                    18794666970.4974822998046875D0,
     .                     1975666666.8475265502929687D0,
     .                   -10029265186.3187522888183594D0,
     .                     8752751870.9984760284423828D0,
     .                    -3744696142.4763402938842773D0,
     .                      642052566.5320791006088257D0,
     .                      243727811.0157482624053955D0,
     .                     -217165686.7454786300659180D0,
     .                        5828609.7747744321823120D0,
     .                      153947658.6124812960624695D0,
     .                     -300528150.9194294214248657D0,
     .                      230839139.9471551179885864D0,
     .                      789711885.8516000509262085D0,
     .                    -3482949095.6228971481323242D0,
     .                     2618973302.6808042526245117D0,
     .                    10489599514.6521968841552734D0,
     .                   -25827817391.9594001770019531D0,
     .                   -20708933596.1514816284179688D0,
     .                   130867888756.2688293457031250D0,
     .                  1381820997896.9057617187500000D0,
     .                  2239332970997.7177734375000000D0,
     .                 -7339172797103.2480468750000000D0,
     .                 -1512559485394.0097656250000000D0,
     .                -20416563869331.5156250000000000D0,
     .                  5801778249341.4941406250000000D0,
     .                 -3272640991176.0688476562500000D0,
     .                 23023224459305.7617187500000000D0,
     .                -14279684469952.0312500000000000D0,
     .                  9878655406729.2597656250000000D0,
     .                   590302760200.1134033203125000D0,
     .                 24745339329742.5507812500000000D0,
     .                  6279085915066.8613281250000000D0,
     .                 -7878234502244.2158203125000000D0,
     .                -26419747684044.6796875000000000D0,
     .                  7574364832961.1171875000000000D0,
     .                  -338005458114.5653686523437500D0,
     .                  -169396000886.1477050781250000D0,
     .                    15586172609.2707824707031250D0,
     .                    -7952530500.4545593261718750D0,
     .                    24305429754.9551811218261719D0,
     .                   -23257676711.7902526855468750D0,
     .                    10748101383.3869400024414063D0,
     .                    -1823507855.0064978599548340D0,
     .                     -707262503.8284575939178467D0,
     .                      565598611.2406927347183228D0,
     .                      -56357209.8536190986633301D0,
     .                     -283222515.8532115221023560D0,
     .                      556274160.4966069459915161D0,
     .                     -324513589.2625260949134827D0,
     .                    -1939735433.7637796401977539D0,
     .                     7392954155.5727138519287109D0,
     .                    -4982684893.1262416839599609D0,
     .                   -24374231575.9767875671386719D0,
     .                    41448108770.8028259277343750D0,
     .                   103710646332.6472930908203125D0,
     .                  -337477024937.7340087890625000D0,
     .                 -2276621954637.5639648437500000D0,
     .                  7401393036979.4287109375000000D0,
     .                  4927379018325.9042968750000000D0,
     .                  8516967326662.8398437500000000D0,
     .                -20629225236270.2421875000000000D0,
     .                   721646358223.1562500000000000D0,
     .                -11654319852045.3789062500000000D0,
     .                 13208965804385.1230468750000000D0,
     .                -12212030214642.7187500000000000D0,
     .                 17183237715586.6523437500000000D0,
     .                  4025882021968.5434570312500000D0,
     .                 18500630474209.2460937500000000D0,
     .                 -5025772727722.9912109375000000D0,
     .                -15083751399167.7304687500000000D0,
     .                -14548681695730.2167968750000000D0,
     .                  8850266457231.1933593750000000D0,
     .                 -2108943297599.0261230468750000D0,
     .                   471076757520.5272216796875000D0,
     .                   -42779441580.0658569335937500D0,
     .                    -5765806117.6573944091796875D0,
     .                   -16756666033.7788391113281250D0,
     .                    22603979167.1909561157226563D0,
     .                   -11241598864.0654296875000000D0,
     .                     1534830023.4662623405456543D0,
     .                      937038526.5187036991119385D0,
     .                     -551791719.5220906734466553D0,
     .                       40277264.4967642128467560D0,
     .                      182820037.1963351368904114D0,
     .                     -291183896.8827682733535767D0,
     .                     -129555331.1161661148071289D0,
     .                     2222932574.4976754188537598D0,
     .                    -6359710174.7713127136230469D0,
     .                     3350095298.1675925254821777D0,
     .                    21789839235.7062454223632813D0,
     .                    -1575552960.4203400611877441D0,
     .                  -260937464237.6711425781250000D0,
     .                  1148715993808.2294921875000000D0,
     .                 -3859934395891.7832031250000000D0,
     .                  1087893130559.9753417968750000D0,
     .                 10684373019950.2363281250000000D0,
     .                  8033190209448.5927734375000000D0,
     .                -12519369807017.7890625000000000D0,
     .                 -1066858008264.0490722656250000D0,
     .                -13865578884546.9667968750000000D0,
     .                 10599382063606.5273437500000000D0,
     .                -18591732460053.2382812500000000D0,
     .                 12351560573662.4472656250000000D0,
     .                 -1899540964016.7661132812500000D0,
     .                  7237621376423.4775390625000000D0,
     .                 -8982701107344.6171875000000000D0,
     .                 -4178244936622.8925781250000000D0,
     .                 18210607559331.7148437500000000D0,
     .                 -4963671042587.4042968750000000D0,
     .                   968564127850.5185546875000000D0,
     .                  -154061560913.6757507324218750D0,
     .                   -16392586733.2933959960937500D0,
     .                    23272383229.0142364501953125D0,
     .                     -797985909.0328063964843750D0,
     .                    -8690602970.8906517028808594D0,
     .                     4871515241.0221519470214844D0,
     .                     -265788707.6617106199264526D0,
     .                     -636673300.3734767436981201D0,
     .                      279403702.6274839639663696D0,
     .                      -46718272.4635147452354431D0,
     .                       -1344925.1677712649106979D0,
     .                      -17442815.9592985957860947D0,
     .                      324195823.6357161998748779D0,
     .                    -1158776838.0914447307586670D0,
     .                     2272965964.1495513916015625D0,
     .                     -733713786.0286161899566650D0,
     .                    -5820149318.2469282150268555D0,
     .                   -22089674528.2214202880859375D0,
     .                   183222844747.4926757812500000D0,
     .                  -808311335583.6230468750000000D0,
     .                  3571533896406.7416992187500000D0,
     .                -14572718869862.7871093750000000D0,
     .                  4346483475016.1572265625000000D0,
     .                  8899699958203.0371093750000000D0,
     .                 -4173671455862.1914062500000000D0,
     .                  3506817924737.5380859375000000D0,
     .                -11776987755741.6367187500000000D0,
     .                 10851160638301.1328125000000000D0,
     .                 -1604217426731.5117187500000000D0,
     .                 19248735432167.0546875000000000D0,
     .                 -1107932701153.6071777343750000D0,
     .                 -9116794856694.5605468750000000D0,
     .                -20869808638644.9140625000000000D0,
     .                 -8663215502320.9609375000000000D0,
     .                 25637094843819.1210937500000000D0,
     .                 -4261818303766.5634765625000000D0,
     .                   890808217896.7390136718750000D0,
     .                  -188388884633.0145874023437500D0,
     .                    48816140355.2877960205078125D0,
     .                   -18883844111.9981842041015625D0,
     .                     6352108094.0888786315917969D0,
     .                     -644699289.5473518371582031D0,
     .                      520219257.3358759880065918D0,
     .                     -947946447.9473376274108887D0,
     .                      483123229.9674471616744995D0,
     .                     -268817731.2039058208465576D0,
     .                      218443303.4020178318023682D0,
     .                     -115945066.4690704941749573D0,
     .                      -52927034.7453688383102417D0,
     .                      132531179.8238504230976105D0,
     .                     -266949652.8533353209495544D0,
     .                      222655788.1611878275871277D0,
     .                      721322582.1887493133544922D0,
     .                   -14882193061.4029617309570313D0,
     .                    46653971538.2760620117187500D0,
     .                   -59146039590.5998535156250000D0,
     .                   -61566811290.6383056640625000D0,
     .                   936934011912.0330810546875000D0,
     .                -11168086093819.7773437500000000D0,
     .                  4368559794426.4912109375000000D0,
     .                 13722264982609.2031250000000000D0,
     .                  3601724899869.4707031250000000D0,
     .                  1546268537336.3879394531250000D0,
     .                -13775801987507.3164062500000000D0,
     .                   856329330257.9182128906250000D0/
      DATA VCOEF2/               
     .               8220216028235.2568359375000000D0,
     .              26447452726942.8515625000000000D0,
     .                551949529565.3088378906250000D0,
     .             -11987385475578.8515625000000000D0,
     .             -19228423172860.7812500000000000D0,
     .             -10543738386122.5605468750000000D0,
     .               2309117437765.6181640625000000D0,
     .               5307887712316.2851562500000000D0,
     .              -1306299325370.2673339843750000D0,
     .                252851577950.7065429687500000D0,
     .                -26256349491.8215942382812500D0,
     .                  7761527003.2064208984375000D0,
     .                 -7670577887.6835937500000000D0,
     .                  4620598656.8596324920654297D0,
     .                 -4676284632.3404445648193359D0,
     .                  3482470911.3160858154296875D0,
     .                 -1380703943.7851724624633789D0,
     .                  1046856219.8273091316223145D0,
     .                 -1076105847.0909764766693115D0,
     .                   456300667.0229959487915039D0,
     .                   790843655.9011635780334473D0,
     .                 -1967486576.4446916580200195D0,
     .                  3148504871.8771948814392090D0,
     .                 -1707636790.7984004020690918D0,
     .                 -5086521099.4443778991699219D0,
     .                 77926123779.2750244140625000D0,
     .               -225526036019.5408630371093750D0,
     .                279893998771.3645019531250000D0,
     .                 11770460275.7504882812500000D0,
     .              -1171227163159.7922363281250000D0,
     .               -464714824014.8496093750000000D0,
     .               1911511498742.6157226562500000D0,
     .               2559414099766.9941406250000000D0,
     .              10322441018759.4765625000000000D0,
     .               2185049386607.2338867187500000D0,
     .             -14519947254135.1093750000000000D0,
     .               -961705366136.2264404296875000D0,
     .              -4958591855461.0302734375000000D0,
     .              21287239040177.5078125000000000D0,
     .              -3689367153055.0024414062500000D0,
     .             -15025495198452.3750000000000000D0,
     .              -5898260086726.1835937500000000D0,
     .               3887558331445.5478515625000000D0,
     .               5666087350048.8300781250000000D0,
     .              -1672493945047.9833984375000000D0,
     .                480658026961.5245361328125000D0,
     .                -49820609910.1936187744140625D0,
     .                -43627145481.6951293945312500D0,
     .                 17095551278.9188537597656250D0,
     .                 -1248053634.7038269042968750D0,
     .                 -1107051217.1276702880859375D0,
     .                  4821037073.3707714080810547D0,
     .                 -4744803455.0864048004150391D0,
     .                  2000755377.0037732124328613D0,
     .                 -1656895660.0225710868835449D0,
     .                  1943436526.4409995079040527D0,
     .                  -937022410.1881250143051147D0,
     .                 -1408608812.2289006710052490D0,
     .                  3705858811.5636377334594727D0,
     .                 -5731580479.8251380920410156D0,
     .                  2771771306.2281322479248047D0,
     .                  8938916015.1393356323242187D0,
     .               -129399404252.4094238281250000D0,
     .                380005720332.2591552734375000D0,
     .               -517273160409.0620117187500000D0,
     .                247002329907.1469726562500000D0,
     .                720024932088.9042968750000000D0,
     .               -995212556035.7690429687500000D0,
     .              -5561246809293.6630859375000000D0,
     .                322895031894.5234375000000000D0,
     .              13711909601765.4941406250000000D0,
     .               3673752846602.9746093750000000D0,
     .             -15283968829436.8242187500000000D0,
     .               3423180287047.2802734375000000D0,
     .              -7293060084983.3388671875000000D0,
     .              14714280931155.5957031250000000D0,
     .              -9919442606611.4472656250000000D0,
     .             -16646296192142.7539062500000000D0,
     .               2195322765273.2231445312500000D0,
     .              15160746257399.5839843750000000D0,
     .               5956339940684.0507812500000000D0,
     .              -4968151966137.5312500000000000D0,
     .                966029647003.2592773437500000D0,
     .               -226879592636.1197509765625000D0,
     .                 84290240837.0091857910156250D0,
     .                -30924852973.8851928710937500D0,
     .                 12250981391.4344482421875000D0,
     .                 -5637215817.1196823120117187D0,
     .                  -664265288.3207979202270508D0,
     .                  2564744608.2758798599243164D0,
     .                 -1116111500.2637720108032227D0,
     .                   942927181.0760573148727417D0,
     .                 -1310673733.9121353626251221D0,
     .                   754212730.2803201675415039D0,
     .                   832276401.3179782629013062D0,
     .                 -2447430105.6201095581054687D0,
     .                  3836920543.2272787094116211D0,
     .                 -2070348355.6171541213989258D0,
     .                 -5374375751.8770542144775391D0,
     .                 81086391398.6670989990234375D0,
     .               -245830418189.0079650878906250D0,
     .                380918821131.2060546875000000D0,
     .               -424139121734.2021484375000000D0,
     .                743000657801.3041992187500000D0,
     .               4710627107113.8291015625000000D0,
     .             -12038420777900.2343750000000000D0,
     .              -2401804960901.9628906250000000D0,
     .              13089403103463.6523437500000000D0,
     .               6052822574424.1757812500000000D0,
     .             -13519758034427.3320312500000000D0,
     .               3577278421658.1083984375000000D0,
     .              -5037408112710.8193359375000000D0,
     .              15844673652854.8906250000000000D0,
     .              -5828891524938.2197265625000000D0,
     .             -12053627275681.4179687500000000D0,
     .              16353112359553.0898437500000000D0,
     .              12176561403799.3320312500000000D0,
     .             -25799888400215.1562500000000000D0,
     .               5179970247014.4296875000000000D0,
     .               -991700742390.5556640625000000D0,
     .                194052768037.0940551757812500D0,
     .                -50204955914.2410583496093750D0,
     .                 20128319174.9123268127441406D0,
     .                -10970709069.6577606201171875D0,
     .                  5536779874.8109855651855469D0,
     .                  -723551277.1913290023803711D0,
     .                  -846169605.2305299043655396D0,
     .                   243691223.1203035712242126D0,
     .                  -151724254.0223848819732666D0,
     .                   354106351.0228716135025024D0,
     .                  -272175223.2592694759368896D0,
     .                  -152595651.7846420407295227D0,
     .                   648834175.1033505201339722D0,
     .                 -1191870985.6247887611389160D0,
     .                  1135546120.3060653209686279D0,
     .                   908123551.5413062572479248D0,
     .                -22682835740.2995109558105469D0,
     .                 75264173834.1635894775390625D0,
     .               -157410230009.0435791015625000D0,
     .                388314815296.2020263671875000D0,
     .              -1538005659555.2685546875000000D0,
     .              12380075212602.0351562500000000D0,
     .              -9991029308577.0351562500000000D0,
     .             -13140186678000.0117187500000000D0,
     .              15874351819446.3476562500000000D0,
     .               8739932696242.3828125000000000D0,
     .             -12929357112787.7851562500000000D0,
     .                319437152844.7934570312500000D0,
     .             -11249104392900.8984375000000000D0,
     .               8084300269338.1464843750000000D0,
     .             -12177217025469.6640625000000000D0,
     .             -11462735708671.3945312500000000D0,
     .              25031059640259.0898437500000000D0,
     .              18285563965238.3789062500000000D0,
     .             -15559313390838.1562500000000000D0,
     .              -1186698442056.1110839843750000D0,
     .                275662483061.0477294921875000D0,
     .                -50157530705.0065612792968750D0,
     .                 13352975741.4750823974609375D0,
     .                 -7904614583.1340980529785156D0,
     .                  5118789605.3659610748291016D0,
     .                 -2412832677.4355220794677734D0,
     .                   135156301.7290229797363281D0,
     .                   429711482.2719335556030273D0,
     .                   -60027379.2486447095870972D0,
     .                    13715847.7013690769672394D0,
     .                  -102471166.1177999526262283D0,
     .                    93664468.2787855565547943D0,
     .                    26552367.8656969368457794D0,
     .                  -190399730.6645207405090332D0,
     .                   420538886.8050022721290588D0,
     .                  -513415970.9205474853515625D0,
     .                  -175821017.5376217365264893D0,
     .                  6898923872.4506559371948242D0,
     .                -22428120031.6482086181640625D0,
     .                 47692059074.1427612304687500D0,
     .               -110131198845.0323486328125000D0,
     .                225430147206.2656250000000000D0,
     .               8114425283362.9326171875000000D0,
     .              -8146543430932.3476562500000000D0,
     .             -16850595221917.1523437500000000D0,
     .              16276105851111.3808593750000000D0,
     .               7968276754929.9794921875000000D0,
     .              -8684282233203.8789062500000000D0,
     .               1175559759653.9208984375000000D0,
     .               1520710241440.9863281250000000D0,
     .              15862220834187.4804687500000000D0,
     .              -6950808302544.2216796875000000D0,
     .             -13486684467802.4785156250000000D0,
     .               7833636940684.9316406250000000D0,
     .             -11237185908634.2812500000000000D0,
     .               6332300442518.1757812500000000D0,
     .                173767218830.6147460937500000D0,
     .                -54265258544.4855346679687500D0,
     .                  7766973695.8122100830078125D0,
     .                 -3895029150.1319007873535156D0,
     .                  3744087115.8826088905334473D0,
     .                 -1709669593.1113615036010742D0,
     .                   287528296.1183981895446777D0,
     .                   374229305.3185383677482605D0,
     .                  -280762509.7946823239326477D0,
     .                    -7900812.9999456778168678D0,
     .                    33451754.2139889001846313D0,
     .                    15702151.9227334856987000D0,
     .                   -33428098.7244412153959274D0,
     .                     3233465.0523898601531982D0,
     .                    45785918.8055144548416138D0,
     .                   -99623032.4652772545814514D0,
     .                   -92408736.1022633314132690D0,
     .                   424100650.9856808185577393D0,
     .                 -1493701916.2100114822387695D0,
     .                  2323098063.6566238403320312D0,
     .                 -2433294993.7510604858398437D0,
     .                  -753812808.2877273559570313D0,
     .                116863928348.5753784179687500D0,
     .              -7002198022770.1093750000000000D0,
     .              12321578405280.2265625000000000D0,
     .              -9328541283824.9921875000000000D0,
     .              12162337926539.6796875000000000D0,
     .               4761077290252.0175781250000000D0,
     .             -10568208305422.8554687500000000D0,
     .              -2460820226963.3120117187500000D0,
     .               2215232185816.5898437500000000D0,
     .              13855621321575.2988281250000000D0,
     .              -4618533930423.3447265625000000D0,
     .              -9162845083599.4804687500000000D0,
     .               2498153993863.7041015625000000D0,
     .             -19293681468581.9335937500000000D0,
     .              15197428243438.3632812500000000D0,
     .               -824891356442.2226562500000000D0,
     .                158222183174.9830322265625000D0,
     .                -30811867593.6039123535156250D0,
     .                  9645096410.2400169372558594D0,
     .                 -4682770315.6128053665161133D0,
     .                  1219987900.6359682083129883D0,
     .                   289174421.6496243476867676D0,
     .                  -748089184.7992402315139771D0,
     .                   399920926.7348154783248901D0,
     .                    42350949.5833997577428818D0,
     .                  -104426677.4717843532562256D0,
     .                    47619151.3822000920772552D0,
     .                     7255412.8735608216375113D0,
     .                    -5705737.5837154388427734D0,
     .                   -29915669.4960079193115234D0,
     .                    71992009.6962184309959412D0,
     .                   190356446.9679546356201172D0,
     .                  -581846327.8431129455566406D0,
     .                  2473058190.2687635421752930D0,
     .                 -6551802711.0572462081909180D0,
     .                 17760139329.1109695434570313D0,
     .                -67396414843.6660766601562500D0,
     .                379455434820.6376953125000000D0,
     .             -10132224626390.5468750000000000D0,
     .              16720872098573.7265625000000000D0,
     .              -6369876382381.8798828125000000D0,
     .               3068972035466.5600585937500000D0,
     .               6291774834841.7304687500000000D0,
     .             -11780505243367.2734375000000000D0,
     .               1875591645022.5205078125000000D0,
     .              -6474152555739.8320312500000000D0,
     .              10537740296272.5468750000000000D0,
     .              -6770227323532.7568359375000000D0,
     .               2402421806985.8867187500000000D0,
     .               7289506507788.5605468750000000D0,
     .             -13159946600383.4941406250000000D0,
     .               6321846139169.0800781250000000D0,
     .                -95180616923.8870849609375000D0,
     .                -72452962352.1259460449218750D0,
     .                 26138727963.7827072143554688D0,
     .                 -7145391744.7685031890869141D0,
     .                   764114756.9109725952148438D0,
     .                  1653245365.0422716140747070D0,
     .                 -2140695146.6031336784362793D0,
     .                  1955881396.4934880733489990D0,
     .                  -815369944.7058086395263672D0,
     .                   -52946271.8839163407683372D0,
     .                   176618551.7260352373123169D0,
     .                  -120472063.2182369232177734D0,
     .                    20772317.4681120514869690D0,
     .                    14231599.4074413776397705D0,
     .                    22338060.1259117871522903D0,
     .                  -114147866.1301196366548538D0,
     .                   101442080.8475131988525391D0,
     .                   284793698.7790446281433105D0,
     .                 -4826622249.7606544494628906D0,
     .                 15376237678.6992721557617188D0,
     .                -25348657368.2009429931640625D0,
     .                  8100734487.4126586914062500D0,
     .                152313578979.7239990234375000D0,
     .              -2914063443589.0346679687500000D0,
     .               9605305511712.4062500000000000D0,
     .             -10605541849400.8613281250000000D0,
     .               -467704889622.3820800781250000D0,
     .               5715938512443.6279296875000000D0,
     .              -7717295854613.3378906250000000D0,
     .               6237448909475.8808593750000000D0,
     .              -8651063127005.9941406250000000D0,
     .               5750225344335.8066406250000000D0,
     .              -7032000838804.8828125000000000D0,
     .               6285192380692.5878906250000000D0,
     .               9995149312854.8261718750000000D0,
     .              -2677529729307.4174804687500000D0,
     .              -3907816824954.7480468750000000D0,
     .                139212355048.5507812500000000D0,
     .                137249523842.3317565917968750D0,
     .                -45801943811.9998321533203125D0,
     .                  4942200024.0634193420410156D0,
     .                  5895536882.5649576187133789D0,
     .                 -6277277244.2456703186035156D0,
     .                  5284083758.4443092346191406D0,
     .                 -4273516480.6779966354370117D0,
     .                  1720637090.6336352825164795D0,
     .                   -29116458.5205169320106506D0,
     .                  -163142588.8644876480102539D0,
     .                   111362467.5122568905353546D0,
     .                   -10625941.7149281203746796D0,
     .                   -12735542.5379977226257324D0,
     .                   -64195372.4043535590171814D0,
     .                   236321454.2009912729263306D0,
     .                  -547458534.2822141647338867D0,
     .                   161673986.7235527038574219D0,
     .                  7157183685.7939910888671875D0,
     .                -24736827468.1243438720703125D0,
     .                 40234211236.1340637207031250D0,
     .                -14488111996.9454345703125000D0,
     .               -201779534832.3078613281250000D0,
     .                947918452917.1933593750000000D0,
     .               4425020599166.8544921875000000D0,
     .              -6552528355754.3583984375000000D0,
     .              -5874810657664.4160156250000000D0,
     .               5925699494370.4667968750000000D0,
     .              -4505867851633.1113281250000000D0,
     .               5828391196114.2597656250000000D0,
     .             -10465794406824.6835937500000000D0,
     .               5127153448242.4833984375000000D0,
     .              -8170919124206.2070312500000000D0,
     .               9414393030070.9492187500000000D0,
     .               6167978218491.9101562500000000D0,
     .               3853181110803.7827148437500000D0,
     .              -6301628779913.2500000000000000D0,
     .                694802762131.4688720703125000D0,
     .               -403891898502.1836547851562500D0,
     .                 89862265123.2917785644531250D0,
     .                  7272744485.6880569458007812D0,
     .                -23533509952.2748565673828125D0,
     .                 18049806279.0352478027343750D0,
     .                -12640203948.6520729064941406D0,
     .                  8933626350.2122344970703125D0,
     .                 -3540164411.7511091232299805D0,
     .                   292365523.3059643507003784D0,
     .                   122662135.6453586220741272D0,
     .                  -188320449.3270677030086517D0,
     .                   162383614.0140390992164612D0,
     .                  -151622924.0571578443050385D0,
     .                   234236425.5373077392578125D0,
     .                  -443977999.8625612258911133D0,
     .                  1221152533.8946495056152344D0,
     .                 -1015384952.0651664733886719D0,
     .                 -7815351382.8600463867187500D0,
     .                 30537095552.7946777343750000D0,
     .                -59111526812.6416320800781250D0,
     .                113060233793.2938232421875000D0,
     .               -425986691662.6858520507812500D0,
     .               4282847616928.8808593750000000D0,
     .              -1371476539924.7998046875000000D0,
     .              -3375029104790.0053710937500000D0,
     .              -7282576035286.3642578125000000D0,
     .               3799513099801.1054687500000000D0,
     .              -2886091375101.2583007812500000D0,
     .               7182216151036.9355468750000000D0,
     .              -7251995331441.7636718750000000D0,
     .               5897344583599.9384765625000000D0,
     .              -5704581413350.4296875000000000D0,
     .               7800659088591.6894531250000000D0,
     .              -1861930971035.3278808593750000D0,
     .               6284992802555.4882812500000000D0,
     .              -4922718307619.1884765625000000D0,
     .               -661531221993.2785644531250000D0,
     .                522199805831.6440429687500000D0,
     .                -96805684155.6965332031250000D0,
     .                -29406975140.3800735473632813D0,
     .                 44059404153.9812164306640625D0,
     .                -33372467212.8697471618652344D0,
     .                 21496183480.3173141479492188D0,
     .                -12268070936.4581127166748047D0,
     .                  4238355427.6266021728515625D0,
     .                  -335303348.1418012380599976D0,
     .                  -257246152.4487919509410858D0,
     .                   491273283.7460340261459351D0,
     .                  -492143004.6952762603759766D0,
     .                   357106955.3916598558425903D0,
     .                  -278820602.6009342670440674D0,
     .                   336952845.0472679138183594D0,
     .                 -1211128615.2075653076171875D0,
     .                  1363446673.4822311401367187D0,
     .                  5410784585.7003173828125000D0,
     .                -22005680772.6704711914062500D0,
     .                 33252696564.4508819580078125D0,
     .                -48314480560.5352020263671875D0,
     .                213325656687.1235961914062500D0,
     .               2144957869800.7919921875000000D0,
     .              -5511413291017.0371093750000000D0,
     .               4674740955202.7226562500000000D0,
     .              -5910503852593.8886718750000000D0,
     .               1593672457265.0187988281250000D0,
     .              -2785132043007.9140625000000000D0,
     .               5611655019547.7226562500000000D0,
     .              -3356749954160.7441406250000000D0,
     .               7523474279207.3378906250000000D0,
     .              -4130477232861.9804687500000000D0,
     .               6941432509876.8027343750000000D0,
     .             -12403715248484.1464843750000000D0,
     .               5937935259901.8115234375000000D0,
     .               -990985456143.5292968750000000D0,
     .                844448536803.8811035156250000D0,
     .               -449569258762.5197753906250000D0,
     .                 78228830720.7163085937500000D0,
     .                 24987558152.6415023803710938D0,
     .                -36884048342.9413223266601563D0,
     .                 30314102280.4933319091796875D0,
     .                -18981748930.4731674194335938D0,
     .                  8685950392.0432720184326172D0,
     .                 -2264155115.3290820121765137D0,
     .                    22728953.1041066646575928D0,
     .                   346735804.9413455724716187D0,
     .                  -553403214.0774860382080078D0,
     .                   502688612.9758108854293823D0,
     .                  -267092390.1482289135456085D0,
     .                    51107766.3393158912658691D0,
     .                    81719014.9087758064270020D0,
     .                   251635956.1480102539062500D0,
     .                  -559129809.5795364379882813D0,
     .                 -1732764011.6990051269531250D0,
     .                  5975260251.2357177734375000D0,
     .                  4387824302.5941543579101562D0,
     .                -44250317470.7684326171875000D0,
     .                 87165286532.6243896484375000D0,
     .               1170284412357.0468750000000000D0,
     .              -7466354435266.9824218750000000D0,
     .              12788554295410.5058593750000000D0,
     .              -6863711816382.2421875000000000D0,
     .               1619122581711.6601562500000000D0,
     .              -5765865847339.2246093750000000D0,
     .               4466668604559.3759765625000000D0,
     .              -4387308921012.9409179687500000D0,
     .              11270263957040.3886718750000000D0,
     .               -675826025815.3356933593750000D0,
     .               5414161624847.9208984375000000D0,
     .             -20980023178992.8398437500000000D0,
     .               5564974323239.5351562500000000D0,
     .               4238398740535.0581054687500000D0,
     .               -604931731175.8751220703125000D0,
     .                196058101700.6907348632812500D0,
     .                -35709301636.5499801635742188D0,
     .                 -5552898794.7671737670898437D0,
     .                 11568881335.2105350494384766D0,
     .                -10612060291.2894058227539063D0,
     .                  6655627276.8213119506835937D0,
     .                 -2581457582.8096332550048828D0,
     .                   444608608.1544632911682129D0,
     .                    74011356.9334512650966644D0,
     .                  -155223583.6997753977775574D0,
     .                   213459407.2186323106288910D0,
     .                  -173742191.5832203626632690D0,
     .                    62268542.1875815540552139D0,
     .                    43656880.8476594686508179D0,
     .                  -126669453.0354866981506348D0,
     .                   125231374.0341720581054688D0,
     .                    -2426552.5528583526611328D0,
     .                   272356835.9782714843750000D0,
     .                  -304929447.5205993652343750D0,
     .                 -5474794614.6634063720703125D0,
     .                 21796986106.1268310546875000D0,
     .                -11804646432.7533798217773438D0,
     .               -544479723845.8056640625000000D0,
     .              -7727591469700.6582031250000000D0,
     .              18298597138319.2617187500000000D0,
     .              -6581645555475.3554687500000000D0,
     .              -1070008013507.2830810546875000D0,
     .              -6797407662036.6054687500000000D0,
     .               4418009457006.3457031250000000D0/
C Kernel parameters
      M=2
C
C      IF (RR.LE.3.2D0) THEN
C      RR=3.2D0
C      ENDIF
      DO I=1,NPT 
      IF (I.LE.481) THEN 
       VCOEF(I)=VCOEF1(I)
       R(I)=R1(I)
       THETA(I)=THETA1(I)
      ELSE
       VCOEF(I)=VCOEF2(I-481)
       R(I)=R2(I-481)
       THETA(I)=THETA2(I-481)
      ENDIF
      ENDDO
      PI=DACOS(-1.D0)
      FACT=PI/180.D0
C      EVALUATE PES    
       SUMA=ZERO
       XTT=(1.d0-DCOS(TT*FACT))/2.d0
       DO I=1,NPT
        TI=(1.d0-DCOS(THETA(I)*FACT))/2.d0
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)*RKHS_AKN2(TI,XTT)
       END DO
       IF (SUMA.GT.1D0) THEN
       SUMA=1.D0
       ENDIF
       IF (SUMA.LT.-1D0) THEN
       SUMA=-1.D0
       ENDIF


       IF(TT.EQ.0.D0.OR.TT.EQ.180.D0) THEN 
       SIN2ALFA=0.D0
       ELSE
       SIN2ALFA=SUMA
       ENDIF
       RETURN
       END



      FUNCTION V_KROH_A(RR,TT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=405)
      DIMENSION VAB(NPT),VCOEF(NPT),VCOEF1(198),VCOEF2(207)
      DIMENSION R1(198),R2(207),THETA1(198),THETA2(207) 
      DIMENSION R(NPT),THETA(NPT)
      DATA R1/
     .              3.000D0,
     .              3.100D0,
     .              3.250D0,
     .              3.300D0,
     .              3.400D0,
     .              3.500D0,
     .              3.750D0,
     .              3.800D0,
     .              3.900D0,
     .              4.000D0,
     .              4.100D0,
     .              4.200D0,
     .              4.250D0,
     .              4.500D0,
     .              4.750D0,
     .              4.800D0,
     .              4.900D0,
     .              5.000D0,
     .              5.150D0,
     .              5.250D0,
     .              5.350D0,
     .              5.500D0,
     .              5.750D0,
     .              6.000D0,
     .              6.250D0,
     .              6.500D0,
     .              6.750D0,
     .              7.000D0,
     .              7.250D0,
     .              7.500D0,
     .              7.750D0,
     .              8.000D0,
     .              8.250D0,
     .              8.500D0,
     .              9.000D0,
     .              9.500D0,
     .             10.000D0,
     .             12.000D0,
     .             14.000D0,
     .             16.000D0,
     .             18.000D0,
     .             20.000D0,
     .             22.000D0,
     .             24.000D0,
     .              2.500D0,
     .              2.750D0,
     .              3.000D0,
     .              3.100D0,
     .              3.250D0,
     .              3.300D0,
     .              3.400D0,
     .              3.500D0,
     .              3.750D0,
     .              3.800D0,
     .              3.900D0,
     .              4.000D0,
     .              4.100D0,
     .              4.200D0,
     .              4.250D0,
     .              4.500D0,
     .              4.750D0,
     .              4.800D0,
     .              4.900D0,
     .              5.000D0,
     .              5.150D0,
     .              5.250D0,
     .              5.350D0,
     .              5.500D0,
     .              5.750D0,
     .              6.000D0,
     .              6.250D0,
     .              6.500D0,
     .              6.750D0,
     .              7.000D0,
     .              7.250D0,
     .              7.500D0,
     .              7.750D0,
     .              8.000D0,
     .              8.250D0,
     .              8.500D0,
     .              9.000D0,
     .              9.500D0,
     .             10.000D0,
     .             12.000D0,
     .             14.000D0,
     .             16.000D0,
     .             18.000D0,
     .             20.000D0,
     .             22.000D0,
     .             24.000D0,
     .              2.500D0,
     .              2.750D0,
     .              3.000D0,
     .              3.100D0,
     .              3.250D0,
     .              3.300D0,
     .              3.400D0,
     .              3.500D0,
     .              3.750D0,
     .              3.800D0,
     .              3.900D0,
     .              4.000D0,
     .              4.100D0,
     .              4.200D0,
     .              4.250D0,
     .              4.500D0,
     .              4.750D0,
     .              4.800D0,
     .              4.900D0,
     .              5.000D0,
     .              5.150D0,
     .              5.250D0,
     .              5.350D0,
     .              5.500D0,
     .              5.750D0,
     .              6.000D0,
     .              6.250D0,
     .              6.500D0,
     .              6.750D0,
     .              7.000D0,
     .              7.250D0,
     .              7.500D0,
     .              7.750D0,
     .              8.000D0,
     .              8.250D0,
     .              8.500D0,
     .              9.000D0,
     .              9.500D0,
     .             10.000D0,
     .             12.000D0,
     .             14.000D0,
     .             16.000D0,
     .             18.000D0,
     .             20.000D0,
     .             22.000D0,
     .             24.000D0,
     .              3.750D0,
     .              3.800D0,
     .              3.900D0,
     .              4.000D0,
     .              4.100D0,
     .              4.200D0,
     .              4.250D0,
     .              4.500D0,
     .              4.750D0,
     .              4.800D0,
     .              4.900D0,
     .              5.000D0,
     .              5.150D0,
     .              5.250D0,
     .              5.350D0,
     .              5.500D0,
     .              5.750D0,
     .              6.000D0,
     .              6.250D0,
     .              6.500D0,
     .              6.750D0,
     .              7.000D0,
     .              7.250D0,
     .              7.500D0,
     .              7.750D0,
     .              8.000D0,
     .              8.250D0,
     .              8.500D0,
     .              9.000D0,
     .              9.500D0,
     .             10.000D0,
     .             12.000D0,
     .             14.000D0,
     .             16.000D0,
     .             18.000D0,
     .             20.000D0,
     .             22.000D0,
     .             24.000D0,
     .              4.750D0,
     .              4.800D0,
     .              4.900D0,
     .              5.000D0,
     .              5.150D0,
     .              5.250D0,
     .              5.350D0,
     .              5.500D0,
     .              5.750D0,
     .              6.000D0,
     .              6.250D0,
     .              6.500D0,
     .              6.750D0,
     .              7.000D0,
     .              7.250D0,
     .              7.500D0,
     .              7.750D0,
     .              8.000D0,
     .              8.250D0,
     .              8.500D0,
     .              9.000D0,
     .              9.500D0,
     .             10.000D0,
     .             12.000D0/
      DATA R2/
     .           14.000D0,
     .           16.000D0,
     .           18.000D0,
     .           20.000D0,
     .           22.000D0,
     .           24.000D0,
     .            5.000D0,
     .            5.150D0,
     .            5.250D0,
     .            5.350D0,
     .            5.500D0,
     .            5.750D0,
     .            6.000D0,
     .            6.250D0,
     .            6.500D0,
     .            6.750D0,
     .            7.000D0,
     .            7.250D0,
     .            7.500D0,
     .            7.750D0,
     .            8.000D0,
     .            8.250D0,
     .            8.500D0,
     .            9.000D0,
     .            9.500D0,
     .           10.000D0,
     .           12.000D0,
     .           14.000D0,
     .           16.000D0,
     .           18.000D0,
     .           20.000D0,
     .           22.000D0,
     .           24.000D0,
     .            5.150D0,
     .            5.250D0,
     .            5.350D0,
     .            5.500D0,
     .            5.750D0,
     .            6.000D0,
     .            6.250D0,
     .            6.500D0,
     .            6.750D0,
     .            7.000D0,
     .            7.250D0,
     .            7.500D0,
     .            7.750D0,
     .            8.000D0,
     .            8.250D0,
     .            8.500D0,
     .            9.000D0,
     .            9.500D0,
     .           10.000D0,
     .           12.000D0,
     .           14.000D0,
     .           16.000D0,
     .           18.000D0,
     .           20.000D0,
     .           22.000D0,
     .           24.000D0,
     .            4.500D0,
     .            4.750D0,
     .            4.800D0,
     .            4.900D0,
     .            5.000D0,
     .            5.150D0,
     .            5.250D0,
     .            5.350D0,
     .            5.500D0,
     .            5.750D0,
     .            6.000D0,
     .            6.250D0,
     .            6.500D0,
     .            6.750D0,
     .            7.000D0,
     .            7.250D0,
     .            7.500D0,
     .            7.750D0,
     .            8.000D0,
     .            8.250D0,
     .            8.500D0,
     .            9.000D0,
     .            9.500D0,
     .           10.000D0,
     .           12.000D0,
     .           14.000D0,
     .           16.000D0,
     .           18.000D0,
     .           20.000D0,
     .           22.000D0,
     .           24.000D0,
     .            4.100D0,
     .            4.200D0,
     .            4.250D0,
     .            4.500D0,
     .            4.750D0,
     .            4.800D0,
     .            4.900D0,
     .            5.000D0,
     .            5.150D0,
     .            5.250D0,
     .            5.350D0,
     .            5.500D0,
     .            5.750D0,
     .            6.000D0,
     .            6.250D0,
     .            6.500D0,
     .            6.750D0,
     .            7.000D0,
     .            7.250D0,
     .            7.500D0,
     .            7.750D0,
     .            8.000D0,
     .            8.250D0,
     .            8.500D0,
     .            9.000D0,
     .            9.500D0,
     .           10.000D0,
     .           12.000D0,
     .           14.000D0,
     .           16.000D0,
     .           18.000D0,
     .           20.000D0,
     .           22.000D0,
     .           24.000D0,
     .            3.800D0,
     .            3.900D0,
     .            4.000D0,
     .            4.100D0,
     .            4.200D0,
     .            4.250D0,
     .            4.500D0,
     .            4.750D0,
     .            4.800D0,
     .            4.900D0,
     .            5.000D0,
     .            5.150D0,
     .            5.250D0,
     .            5.350D0,
     .            5.500D0,
     .            5.750D0,
     .            6.000D0,
     .            6.250D0,
     .            6.500D0,
     .            6.750D0,
     .            7.000D0,
     .            7.250D0,
     .            7.500D0,
     .            7.750D0,
     .            8.000D0,
     .            8.250D0,
     .            8.500D0,
     .            9.000D0,
     .            9.500D0,
     .           10.000D0,
     .           12.000D0,
     .           14.000D0,
     .           16.000D0,
     .           18.000D0,
     .           20.000D0,
     .           22.000D0,
     .           24.000D0,
     .            2.500D0,
     .            2.750D0,
     .            3.000D0,
     .            3.100D0,
     .            3.250D0,
     .            3.300D0,
     .            3.400D0,
     .            3.500D0,
     .            3.750D0,
     .            3.800D0,
     .            3.900D0,
     .            4.000D0,
     .            4.100D0,
     .            4.200D0,
     .            4.250D0,
     .            4.500D0,
     .            4.750D0,
     .            4.800D0,
     .            4.900D0,
     .            5.000D0,
     .            5.150D0,
     .            5.250D0,
     .            5.350D0,
     .            5.500D0,
     .            5.750D0,
     .            6.000D0,
     .            6.250D0,
     .            6.500D0,
     .            6.750D0,
     .            7.000D0,
     .            7.250D0,
     .            7.500D0,
     .            7.750D0,
     .            8.000D0,
     .            8.250D0,
     .            8.500D0,
     .            9.000D0,
     .            9.500D0,
     .           10.000D0,
     .           12.000D0,
     .           14.000D0,
     .           16.000D0,
     .           18.000D0,
     .           20.000D0,
     .           22.000D0,
     .           24.000D0/

      DATA THETA1/
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .              0.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             20.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             40.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             60.000D0,
     .             80.000D0,
     .             80.000D0,
     .             80.000D0,
     .             80.000D0,
     .             80.000D0,
     .             80.000D0,
     .             80.000D0,
     .             80.000D0,
     .             80.000D0,
     .             80.000D0,
     .             80.000D0,
     .             80.000D0,
     .             80.000D0,
     .             80.000D0,
     .             80.000D0,
     .             80.000D0,
     .             80.000D0,
     .             80.000D0,
     .             80.000D0,
     .             80.000D0,
     .             80.000D0,
     .             80.000D0,
     .             80.000D0,
     .             80.000D0/
       DATA THETA2 /
     .              80.000D0,
     .              80.000D0,
     .              80.000D0,
     .              80.000D0,
     .              80.000D0,
     .              80.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .              90.000D0,
     .             100.000D0,
     .             100.000D0,
     .             100.000D0,
     .             100.000D0,
     .             100.000D0,
     .             100.000D0,
     .             100.000D0,
     .             100.000D0,
     .             100.000D0,
     .             100.000D0,
     .             100.000D0,
     .             100.000D0,
     .             100.000D0,
     .             100.000D0,
     .             100.000D0,
     .             100.000D0,
     .             100.000D0,
     .             100.000D0,
     .             100.000D0,
     .             100.000D0,
     .             100.000D0,
     .             100.000D0,
     .             100.000D0,
     .             100.000D0,
     .             100.000D0,
     .             100.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             120.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             140.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             160.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0,
     .             180.000D0/
     
      DATA VCOEF1/
     .            23092795740.3262100219726562D0, 
     .           -34479574574.5429916381835938D0,
     .            23276661953.8395195007324219D0,
     .           -13736153575.2621154785156250D0,
     .             1508607254.0388567447662354D0,
     .             -444120761.3152050971984863D0,
     .             -324955884.9531007409095764D0,
     .             1212777166.1337935924530029D0,
     .             -675499532.9111161231994629D0,
     .              364433410.7020100951194763D0,
     .             -427015570.9131940007209778D0,
     .              625260098.7054095268249512D0,
     .             -218178753.4012194275856018D0,
     .               13787062.2106401063501835D0,
     .            -1935232618.7997238636016846D0,
     .             2013125491.8292107582092285D0,
     .             -250444274.4069581329822540D0,
     .              674658153.3941773176193237D0,
     .            -1412061503.5253696441650391D0,
     .              950937116.5246207714080811D0,
     .             -260351104.1752514839172363D0,
     .               22862925.6166418865323067D0,
     .               13696462.1397831458598375D0,
     .               61170955.7017747685313225D0,
     .               72046581.6635076999664307D0,
     .               77757319.8328657597303391D0,
     .               65482144.6390266194939613D0,
     .               54770742.3451209142804146D0,
     .               40058173.1051778048276901D0,
     .               32861160.9170640893280506D0,
     .               17788555.7530217282474041D0,
     .               17698519.6325375847518444D0,
     .                 604847.6256132097914815D0,
     .                6417071.1112622125074267D0,
     .               -3905861.8200243311002851D0,
     .               -5377361.0215423218905926D0,
     .              -12405864.3904476072639227D0,
     .              -27250150.5743684619665146D0,
     .                -377385.1440169084235094D0,
     .               -2766733.1514719841070473D0,
     .                1030630.8534696570131928D0,
     .               -4122592.4137536738999188D0,
     .                3981928.7843842366710305D0,
     .               -3342121.8424896649084985D0,
     .            12830653865.3734359741210938D0,
     .           -43772494645.9513473510742188D0,
     .           169174948632.3132019042968750D0,
     .          -214852020202.3634338378906250D0,
     .           164227864476.4625854492187500D0,
     .           -95579488574.7215118408203125D0,
     .            12404825205.8885478973388672D0,
     .            -3337410641.6981334686279297D0,
     .             2133783455.8042824268341064D0,
     .            -4810105038.3495063781738281D0,
     .             1500959092.5701594352722168D0,
     .            -1219761927.5457611083984375D0,
     .              779944733.5794268846511841D0,
     .             -907414835.5075492858886719D0,
     .              502063222.4599015116691589D0,
     .             -244437047.8075387775897980D0,
     .             2991576656.3008041381835938D0,
     .            -3400029815.7168316841125488D0,
     .              226104950.6644559502601624D0,
     .            -1780111387.0827076435089111D0,
     .             2843773529.9348206520080566D0,
     .            -2006576843.6969065666198730D0,
     .              449356748.5986661911010742D0,
     .               34529348.7354581430554390D0,
     .              300723942.4589837789535522D0,
     .              348395708.3899065852165222D0,
     .              349558262.3656947612762451D0,
     .              295615728.8401453495025635D0,
     .              228206860.3118911087512970D0,
     .              177049070.8814584910869598D0,
     .              127865248.6336607486009598D0,
     .               83973279.3382417261600494D0,
     .               71496870.4241472780704498D0,
     .               12238947.3583563435822725D0,
     .               38581891.0445498675107956D0,
     .               -5284212.1503330040723085D0,
     .              -16847319.7345648705959320D0,
     .               20131284.0186607688665390D0,
     .              -76055208.1426191776990891D0,
     .             -138541825.5153657197952271D0,
     .                5089279.8050776114687324D0,
     .               -5717114.1660041064023972D0,
     .              -20304556.0414973571896553D0,
     .                9780071.0784453563392162D0,
     .               -6852475.0399842550978065D0,
     .               -2481930.6819703797809780D0,
     .           -11026125920.3238010406494141D0,
     .            65369621091.4787216186523438D0,
     .          -387100540238.2198486328125000D0,
     .           540406315612.0474243164062500D0,
     .          -435382820750.4566040039062500D0,
     .           253806552376.5062866210937500D0,
     .           -31864172569.4651260375976562D0,
     .             7937582738.4396982192993164D0,
     .             1286171701.8283634185791016D0,
     .            -2268571711.9191851615905762D0,
     .              131379111.9767262339591980D0,
     .             -181666626.7206954658031464D0,
     .            -1454838454.6774349212646484D0,
     .             1745606164.4670231342315674D0,
     .            -1268133346.9400966167449951D0,
     .             -936596806.0129894018173218D0,
     .            -6484792180.2632122039794922D0,
     .             6874076398.2152957916259766D0,
     .             -849778536.1268777847290039D0,
     .             1170521460.7846522331237793D0,
     .            -1520383787.1477541923522949D0,
     .             1176925150.7966823577880859D0,
     .             -114400629.9582824856042862D0,
     .              298105728.1545095443725586D0,
     .              208321175.9000206887722015D0,
     .               94509688.1326137483119965D0,
     .                9256915.5948263090103865D0,
     .              -23509620.1657352410256863D0,
     .              -32327430.9978984370827675D0,
     .              -35427044.7403171584010124D0,
     .              -16830085.0931801237165928D0,
     .              -15678972.1798456851392984D0,
     .               -6008439.2606058884412050D0,
     .                 491177.2796594904502854D0,
     .                2064493.7053639260120690D0,
     .               31199647.0750845409929752D0,
     .               29932187.8241446502506733D0,
     .               75940764.3524071425199509D0,
     .              -53754567.7434320673346519D0,
     .              -30100836.2695675119757652D0,
     .                5446053.1920506739988923D0,
     .               -6066938.3752986053004861D0,
     .               -7566128.7713700924068689D0,
     .                2752261.6420178408734500D0,
     .               -3954871.4784936788491905D0,
     .               -1968156.5354250418022275D0,
     .            -8993126865.7792949676513672D0,
     .             7585567458.0669231414794922D0,
     .             3363668215.3686456680297852D0,
     .             -219499767.0817124545574188D0,
     .            -1573405613.2189488410949707D0,
     .             2927392611.1513614654541016D0,
     .             -944606536.4629153013229370D0,
     .             -835273639.4814721345901489D0,
     .           -15051710691.6903305053710938D0,
     .            16363198044.6101474761962891D0,
     .            -1926083480.2537236213684082D0,
     .              870487938.2682449817657471D0,
     .            -3898263354.5669608116149902D0,
     .             4280197128.6601500511169434D0,
     .            -1100728618.6890413761138916D0,
     .              -13259080.7656563222408295D0,
     .             -315534614.7156488299369812D0,
     .             -268334201.4993233084678650D0,
     .             -241486093.1114495396614075D0,
     .             -193161095.4701724052429199D0,
     .             -137039664.9401140213012695D0,
     .              -98722828.5748409479856491D0,
     .              -59188702.8740815669298172D0,
     .              -11256015.7726005613803864D0,
     .              -22157481.4263973608613014D0,
     .               72801318.6246849894523621D0,
     .              -17188275.9656277038156986D0,
     .              100291680.0607105642557144D0,
     .              105036512.2739842087030411D0,
     .               15474130.8297282271087170D0,
     .              116034465.8642234057188034D0,
     .              106374151.7831316590309143D0,
     .              -10314240.4977144077420235D0,
     .               -2494292.9811759446747601D0,
     .               16800599.3138881400227547D0,
     .              -32558099.3035462312400341D0,
     .               38205810.8568465411663055D0,
     .              -21597090.3387611731886864D0,
     .            45675415515.8980560302734375D0,
     .           -47821422060.9256896972656250D0,
     .             5442723623.6706295013427734D0,
     .           -14607340407.9325275421142578D0,
     .            29311872065.3623580932617188D0,
     .           -22070684868.6398239135742188D0,
     .             5111766155.2340393066406250D0,
     .            -1222688290.2185878753662109D0,
     .             -152070035.6863631308078766D0,
     .             -248061377.3075999617576599D0,
     .              -99643642.0464322715997696D0,
     .              -69347783.0662887990474701D0,
     .              -27762341.8582338765263557D0,
     .               -7254915.1896373881027102D0,
     .               48035452.7850887700915337D0,
     .              -12826790.1786715872585773D0,
     .              133568596.9037150144577026D0,
     .              -25037342.4826090596616268D0,
     .              135724300.3480714261531830D0,
     .               40429672.1492394730448723D0,
     .              138972768.1975903809070587D0,
     .               20628811.2504454441368580D0,
     .              207232665.3169269561767578D0,
     .               99768424.5603480786085129D0/
         DATA VCOEF2/               
     .               -28122246.5028472542762756D0,
     .                22094082.3755243979394436D0,
     .               -37504647.0343720316886902D0,
     .                64147215.5668610632419586D0,
     .               -91979458.7732433229684830D0,
     .                58888456.4456582441926003D0,
     .             15926967599.7907409667968750D0,
     .            -47533584524.6642837524414062D0,
     .             40293969410.8920135498046875D0,
     .             -9832482274.2044658660888672D0,
     .              1558610065.7148866653442383D0,
     .              -309999129.9916262626647949D0,
     .                53637025.5591742545366287D0,
     .               -62077704.8743233159184456D0,
     .               -48211471.7366920337080956D0,
     .               -28544338.0699616856873035D0,
     .                -3371983.8063026177696884D0,
     .               -59504171.2188396975398064D0,
     .                79046014.6358001977205276D0,
     .              -133656103.0452174842357635D0,
     .               123745106.2477076500654221D0,
     .              -106550555.1770142763853073D0,
     .                63622814.1936504691839218D0,
     .                -8235207.7664829222485423D0,
     .                10095885.0227325037121773D0,
     .                 1596889.8445444190874696D0,
     .                12360144.9405959118157625D0,
     .                 2611409.8344933209009469D0,
     .                -6255267.6728565199300647D0,
     .                27887777.4587011262774467D0,
     .               -67916761.8314860165119171D0,
     .               107738477.1377761363983154D0,
     .               -70897171.0120177865028381D0,
     .             23943223663.5134086608886719D0,
     .            -29867603908.0105705261230469D0,
     .              7052852758.0172090530395508D0,
     .             -1482244610.2648630142211914D0,
     .               -58888468.4489063769578934D0,
     .              -202548492.5052129328250885D0,
     .               -88633913.6942991316318512D0,
     .                 -763584.1842180109815672D0,
     .               -31792823.7653017304837704D0,
     .                18969314.7199526689946651D0,
     .                57801054.7767284289002419D0,
     .               -16545550.9542314503341913D0,
     .               139200445.2366498410701752D0,
     .               -27143465.0227306522428989D0,
     .               110544854.7392248809337616D0,
     .                51304615.5349221676588058D0,
     .               121287534.9524920582771301D0,
     .                44931058.5313494056463242D0,
     .               163131871.0476610362529755D0,
     .                75407369.3545038700103760D0,
     .               -26520028.7249931953847408D0,
     .                12979371.7818287573754787D0,
     .               -23362594.9152855053544044D0,
     .                56277539.6633210182189941D0,
     .               -83986330.1559239327907562D0,
     .                48227801.3744531348347664D0,
     .              6171122228.6857881546020508D0,
     .            -25799785113.2354927062988281D0,
     .             22379606317.5750007629394531D0,
     .             -2607634646.2776207923889160D0,
     .              -819482456.4573781490325928D0,
     .             -3254439831.0805683135986328D0,
     .              5382684325.0818786621093750D0,
     .             -1426070788.6073603630065918D0,
     .               132643859.3174381405115128D0,
     .              -170445744.5331101119518280D0,
     .               -82514668.4185730665922165D0,
     .               -93585007.9128130674362183D0,
     .               -33388676.3793441578745842D0,
     .               -51455287.1408639997243881D0,
     .                 1185059.7128611248917878D0,
     .               -20099897.7685120590031147D0,
     .                35239166.4964914470911026D0,
     .               -22537303.4857313148677349D0,
     .                49673138.1289434283971786D0,
     .                -3852866.5897127431817353D0,
     .                44298938.3351686000823975D0,
     .                52994182.3860580101609230D0,
     .                21829134.5077456384897232D0,
     .                68557600.8195753693580627D0,
     .                32975710.0010868944227695D0,
     .               -10616880.5258261431008577D0,
     .                -3713519.0621770191937685D0,
     .                15625966.0655471347272396D0,
     .               -21664548.6310537979006767D0,
     .                14485429.1578396130353212D0,
     .                -4930490.1387307457625866D0,
     .             10899578246.6567840576171875D0,
     .            -17463245187.0746803283691406D0,
     .              5608144215.4342842102050781D0,
     .             -2604283487.4475183486938477D0,
     .              5167863654.8778171539306641D0,
     .             -3181440427.7678885459899902D0,
     .               643413899.9423127174377441D0,
     .              -551936215.0400342941284180D0,
     .              -557560978.7214215993881226D0,
     .              2288474489.6445078849792480D0,
     .              -279114591.3941955566406250D0,
     .               280557143.7190266251564026D0,
     .                93631189.9949005097150803D0,
     .                 5748305.8516647890210152D0,
     .                 2921576.2965152724646032D0,
     .              -109315895.9553678929805756D0,
     .                27883098.8301900476217270D0,
     .               -71708277.1131147146224976D0,
     .               -18775230.7446656264364719D0,
     .               -22090165.4178576730191708D0,
     .               -35024272.2549113482236862D0,
     .                -6046235.9355522012338042D0,
     .               -21849033.2432330399751663D0,
     .               -16860282.7925511896610260D0,
     .               -30748153.2653949446976185D0,
     .               -12506467.5255542546510696D0,
     .               -36694374.9924012944102287D0,
     .               -21036775.8774046674370766D0,
     .                 2056733.6930611152201891D0,
     .                -7378590.8412494836375117D0,
     .                11908016.1380871646106243D0,
     .               -26656157.2819432280957699D0,
     .                39407840.8962592259049416D0,
     .               -23239514.9259936474263668D0,
     .             14936750061.2865066528320312D0,
     .            -20660642822.5058288574218750D0,
     .              4570267493.4560670852661133D0,
     .             -5224348431.4024085998535156D0,
     .              6382674813.0300350189208984D0,
     .             -1756545387.2290563583374023D0,
     .              -488063553.1159021854400635D0,
     .              1058085032.1983376741409302D0,
     .                 7866809.2449813345447183D0,
     .               191590794.3435797393321991D0,
     .               161112371.2858384549617767D0,
     .              2734717886.5814695358276367D0,
     .             -2650623742.7555212974548340D0,
     .               870096175.9684103727340698D0,
     .               111839190.2819045782089233D0,
     .               147560645.0797961354255676D0,
     .                11514712.1891295593231916D0,
     .               -21019319.8598006106913090D0,
     .               -16121700.5322087816894054D0,
     .               -48693205.0565529540181160D0,
     .               -22377942.0051368847489357D0,
     .               -30578202.2145823724567890D0,
     .               -34303961.9715307652950287D0,
     .                -5826288.4207230899482965D0,
     .               -40732498.9296812117099762D0,
     .                -8363792.0733832828700542D0,
     .               -36247933.6974939182400703D0,
     .               -47224182.2863222658634186D0,
     .               -22690556.0834364481270313D0,
     .               -59909945.6751262098550797D0,
     .               -25204902.9004962779581547D0,
     .                 3126066.9368129004724324D0,
     .                -1431896.7226848003920168D0,
     .                -4785134.3649913575500250D0,
     .                 9107246.0375414658337831D0,
     .                -6932985.6672314051538706D0,
     .                 1599356.8020151318050921D0,
     .               668985755.4271870851516724D0,
     .             -3170916805.7788109779357910D0,
     .             16932064861.2310619354248047D0,
     .            -23738632831.5366020202636719D0,
     .             19706505094.3647232055664062D0,
     .            -11479644905.1028327941894531D0,
     .              1449141958.7311894893646240D0,
     .              -378790602.9345162510871887D0,
     .               284102298.4686219692230225D0,
     .             -4460027588.3497619628906250D0,
     .              5688867095.1325941085815430D0,
     .             -1402131353.2420535087585449D0,
     .              -156078333.1812246739864349D0,
     .               535836658.8011304140090942D0,
     .              -322678567.1989965438842773D0,
     .               167540272.2295531630516052D0,
     .              -466928748.5284626483917236D0,
     .               197580812.0104706585407257D0,
     .               -35430611.8552500680088997D0,
     .                -3637772.5080013275146484D0,
     .             -1102593525.7106561660766602D0,
     .              1372035958.2584323883056641D0,
     .              -321082650.0508201122283936D0,
     .                75227938.6422338932752609D0,
     .                 5467206.3604799322783947D0,
     .                 6745426.9343924149870872D0,
     .                -1499194.7421312453225255D0,
     .                  152126.6809985425206833D0,
     .                -5745869.0065905349329114D0,
     .                 -850743.4340431911405176D0,
     .                -4834545.4471013033762574D0,
     .                  -52278.8193659601165564D0,
     .                -7404126.0374316461384296D0,
     .                 1048836.0329411532729864D0,
     .                -4287993.1948530143126845D0,
     .                -5485548.1378772668540478D0,
     .                -8028056.5779850482940674D0,
     .                -3916460.7844366407953203D0,
     .                -8340875.7449064617976546D0,
     .                -4167433.7130199009552598D0,
     .                  655226.0731320534832776D0,
     .                 -522298.6189026608481072D0,
     .                  296907.5643561522592790D0,
     .                  131241.9507320474367589D0,
     .                -1479353.0029812389984727D0,
     .                 1266227.1406982007902116D0/
         
C Kernel parameters
      N1=2
      M=5
      N2=10
      M2=0
C
      DO I=1,NPT 
      IF (I.LE.198) THEN 
       VCOEF(I)=VCOEF1(I)
       R(I)=R1(I)
       THETA(I)=THETA1(I)
      ELSE
       VCOEF(I)=VCOEF2(I-198)
       R(I)=R2(I-198)
       THETA(I)=THETA2(I-198)
      ENDIF
      ENDDO
      PI=DACOS(-1.D0)
      FACT=PI/180.D0
C      EVALUATE PES    
       SUMA=ZERO
       XTT=DCOS(TT*FACT)
       DO I=1,NPT
        TI=DCOS(THETA(I)*FACT)
        SUMA=SUMA+VCOEF(I)*RKHS_DK(R(I),RR,N1,M)*RKHS_OP(TI,XTT,N2,M2)
       END DO
       V_KROH_A=SUMA
       RETURN
       END


      FUNCTION FACT(N)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      FACT=1.D0
      DO  J=2,N
      FACT=FACT*DFLOAT(J)
      ENDDO
      RETURN
      END


      FUNCTION Vsum_rre(R,theta)
C*********************************
C System: Kr-OH(X2Pi)
C Method:RCCSD(T)
C PES: Vsum=1/2(A''+A')
C r=re OH distance
C Variables:
C R in Bohr
C Theta n Degrees
C Output in cm-1
C Theta=0 for Kr---HO geom.
C********************************
C Authors:
C Dr Jacek Klos
C Department of Chemistry and Biochemistry
C University of Maryland
C College Park, MD 20742
C email:jklos@mail.umd.edu
C**********************************
      implicit double precision(a-h, o-z)
      dimension V0(11)
      dimension T(11)
      pi=dacos(-1.d0)
      conv=1.D0 ! IN CM-1
      V0(1)=VL0_0(r)*conv    
      V0(2)=VL0_1(r)*conv
      V0(3)=VL0_2(r)*conv
      V0(4)=VL0_3(r)*conv
      V0(5)=VL0_4(r)*conv
      V0(6)=VL0_5(r)*conv
      V0(7)=VL0_6(r)*conv
      V0(8)=VL0_7(r)*conv
      V0(9)=VL0_8(r)*conv
      V0(10)=VL0_9(r)*conv
      V0(11)=VL0_10(r)*conv
       do j=1,11
       T(j)=PLGNDR((j-1),0,DCOS(theta*pi/180.D0))
       enddo
       s=0.D0
       do i=1,11
       s=s+V0(i)*T(i)
       enddo
       Vsum_rre=s
       return
       end

      FUNCTION Vdif_rre(R,theta)
C*********************************
C System: Kr-OH(X2Pi)
C Method:RCCSD(T)
C PES: Vdif=1/2(A''-A')
C r=re r OH distance
C Variables:
C R in Bohr
C Theta n Degrees
C Output in cm-1
C Theta=0 for Kr---HO geom.
C********************************
C Authors:
C Dr Jacek Klos
C Department of Chemistry and Biochemistry
C University of Maryland
C College Park, MD 20742
C email:jklos@mail.umd.edu
C**********************************
      implicit double precision(a-h, o-z)
      dimension V2(9)
      dimension T2(9)
      pi=dacos(-1.d0)
      conv=1.D0 ! IN CM-1
      V2(1)=VL2_2(r)*conv 
      V2(2)=VL2_3(r)*conv
      V2(3)=VL2_4(r)*conv
      V2(4)=VL2_5(r)*conv
      V2(5)=VL2_6(r)*conv
      V2(6)=VL2_7(r)*conv
      V2(7)=VL2_8(r)*conv
      V2(8)=VL2_9(r)*conv
      V2(9)=VL2_10(r)*conv
       do j=1,9
       T2(j)=DSQRT(FACT(j+1-2)/FACT(j+1+2))*
     .    PLGNDR((j+1),2,DCOS(theta*pi/180.D0))
       enddo
       s=0.D0
       do i=1,9
       s=s+V2(i)*T2(i)
       enddo
       Vdif_rre=s
       return
       end
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN LEGENDRE
C ANGULAR EXPANSION OF Kr-OH(X,r=re) VSUM RCCSD(T) PES
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_0(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=39)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .  3.000D0,
     .  3.250D0,
     .  3.500D0,
     .  3.750D0,
     .  4.000D0,
     .  4.250D0,
     .  4.500D0,
     .  4.750D0,
     .  5.000D0,
     .  5.250D0,
     .  5.500D0,
     .  5.750D0,
     .  6.000D0,
     .  6.250D0,
     .  6.500D0,
     .  6.750D0,
     .  7.000D0,
     .  7.250D0,
     .  7.500D0,
     .  7.750D0,
     .  8.000D0,
     .  8.250D0,
     .  8.500D0,
     .  8.750D0,
     .  9.000D0,
     .  9.250D0,
     .  9.500D0,
     .  9.750D0,
     . 10.000D0,
     . 11.000D0,
     . 12.000D0,
     . 13.000D0,
     . 14.000D0,
     . 15.000D0,
     . 16.000D0,
     . 18.000D0,
     . 20.000D0,
     . 22.000D0,
     . 24.000D0/
      DATA VCOEF/
     *    0.650154996248683357D+10,
     *   -0.538653518640884209D+10,
     *    0.147313802359313440D+10,
     *   -0.454877060238906503D+09,
     *    0.104986145549243677D+10,
     *   -0.123910688400931668D+10,
     *   -0.359160348042816281D+09,
     *   -0.803008585911992908D+09,
     *   -0.702502758386389017D+09,
     *   -0.708718805878646374D+09,
     *   -0.611315128743668556D+09,
     *   -0.499841295542309701D+09,
     *   -0.370335469596802950D+09,
     *   -0.245754723744827390D+09,
     *   -0.126828139676472068D+09,
     *   -0.181504328287209272D+08,
     *    0.732536913714449704D+08,
     *    0.144317449045113742D+09,
     *    0.195583782929069519D+09,
     *    0.224118770170559764D+09,
     *    0.232105078774194121D+09,
     *    0.226596344619236231D+09,
     *    0.215851861692186952D+09,
     *    0.197832908679841161D+09,
     *    0.180824340703056216D+09,
     *    0.148340621732523203D+09,
     *    0.148921336538859129D+09,
     *    0.129828053809721470D+08,
     *    0.262309126399464369D+09,
     *    0.158501765367674351D+09,
     *    0.164200619856357574D+06,
     *    0.650664366078782082D+07,
     *   -0.136737285040930510D+08,
     *   -0.452015258186113834D+07,
     *   -0.180136860422211885D+08,
     *   -0.908826303136044741D+07,
     *   -0.685832386917352676D+06,
     *   -0.257736517491457462D+08,
     *    0.904576673472821712D+05/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_0=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN LEGENDRE
C ANGULAR EXPANSION OF Kr-OH(X,r=re) VSUM RCCSD(T) PES
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_10(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=39)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .  3.000D0,
     .  3.250D0,
     .  3.500D0,
     .  3.750D0,
     .  4.000D0,
     .  4.250D0,
     .  4.500D0,
     .  4.750D0,
     .  5.000D0,
     .  5.250D0,
     .  5.500D0,
     .  5.750D0,
     .  6.000D0,
     .  6.250D0,
     .  6.500D0,
     .  6.750D0,
     .  7.000D0,
     .  7.250D0,
     .  7.500D0,
     .  7.750D0,
     .  8.000D0,
     .  8.250D0,
     .  8.500D0,
     .  8.750D0,
     .  9.000D0,
     .  9.250D0,
     .  9.500D0,
     .  9.750D0,
     . 10.000D0,
     . 11.000D0,
     . 12.000D0,
     . 13.000D0,
     . 14.000D0,
     . 15.000D0,
     . 16.000D0,
     . 18.000D0,
     . 20.000D0,
     . 22.000D0,
     . 24.000D0/
      DATA VCOEF/
     *   -0.625806467471389055D+09,
     *    0.456149117940544844D+09,
     *    0.363472242455147362D+10,
     *   -0.102141850607749214D+11,
     *    0.128059242527312946D+11,
     *   -0.844925619050274563D+10,
     *    0.309319697348318195D+10,
     *   -0.992254902941572428D+09,
     *    0.361639812059536874D+09,
     *   -0.110724416509756744D+09,
     *    0.401027726161134243D+08,
     *   -0.111086456277090907D+08,
     *    0.868500256131374836D+07,
     *    0.524336180788278580D+06,
     *   -0.288939718078970909D+06,
     *    0.211768046420812607D+07,
     *    0.181816120910370350D+07,
     *   -0.422387894880586863D+07,
     *    0.209869672628563643D+07,
     *    0.130574609131479263D+07,
     *   -0.220289547643192112D+07,
     *    0.309583565291321278D+07,
     *   -0.212994013590168953D+07,
     *    0.212978683960556984D+07,
     *   -0.648777725316405296D+07,
     *    0.157627177353174686D+08,
     *   -0.228374644308629036D+08,
     *    0.174323891517794132D+08,
     *   -0.498897370571255684D+07,
     *   -0.689724810494184494D+05,
     *   -0.891602497327327728D+06,
     *    0.152730455986070633D+07,
     *   -0.105895729438066483D+07,
     *    0.375997099634647369D+06,
     *   -0.108406684974193573D+06,
     *    0.206886174726486206D+05,
     *   -0.872024910211563110D+04,
     *    0.335501787281036377D+04,
     *   -0.873019350051879883D+03/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_10=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN LEGENDRE
C ANGULAR EXPANSION OF Kr-OH(X,r=re) VSUM RCCSD(T) PES
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_1(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=39)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .  3.000D0,
     .  3.250D0,
     .  3.500D0,
     .  3.750D0,
     .  4.000D0,
     .  4.250D0,
     .  4.500D0,
     .  4.750D0,
     .  5.000D0,
     .  5.250D0,
     .  5.500D0,
     .  5.750D0,
     .  6.000D0,
     .  6.250D0,
     .  6.500D0,
     .  6.750D0,
     .  7.000D0,
     .  7.250D0,
     .  7.500D0,
     .  7.750D0,
     .  8.000D0,
     .  8.250D0,
     .  8.500D0,
     .  8.750D0,
     .  9.000D0,
     .  9.250D0,
     .  9.500D0,
     .  9.750D0,
     . 10.000D0,
     . 11.000D0,
     . 12.000D0,
     . 13.000D0,
     . 14.000D0,
     . 15.000D0,
     . 16.000D0,
     . 18.000D0,
     . 20.000D0,
     . 22.000D0,
     . 24.000D0/
      DATA VCOEF/
     *    0.790456848902496815D+10,
     *   -0.108711464382417641D+11,
     *    0.591628826792322540D+10,
     *   -0.464124585421898079D+10,
     *    0.342595716552711821D+10,
     *   -0.221981375188533688D+10,
     *    0.827523980086966515D+09,
     *   -0.316329622176288366D+09,
     *   -0.308938889089569747D+08,
     *   -0.195506331817254305D+09,
     *   -0.167330331021493137D+09,
     *   -0.177677356851020157D+09,
     *   -0.156968311134513974D+09,
     *   -0.127084768843011320D+09,
     *   -0.778797804117343128D+08,
     *   -0.251017932370996550D+08,
     *    0.142399033007270545D+08,
     *    0.411271055448095351D+08,
     *    0.669467329929246604D+08,
     *    0.834515586525505185D+08,
     *    0.959325723696631193D+08,
     *    0.986405307259423733D+08,
     *    0.965931036518476009D+08,
     *    0.932149832659511566D+08,
     *    0.817321512292857170D+08,
     *    0.699854339923188686D+08,
     *    0.639025993918507099D+08,
     *    0.138183152639424801D+08,
     *    0.998617790556564331D+08,
     *    0.642591207184023857D+08,
     *   -0.108211018935341835D+08,
     *   -0.636658587855052948D+07,
     *   -0.907289298352813721D+07,
     *   -0.324523997056162357D+08,
     *    0.398057408762412071D+08,
     *    0.349873322143146992D+08,
     *   -0.113348285562472582D+09,
     *   -0.915380662447822094D+07,
     *    0.679594993086765409D+08/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_1=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN LEGENDRE
C ANGULAR EXPANSION OF Kr-OH(X,r=re) VSUM RCCSD(T) PES
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_2(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=39)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .  3.000D0,
     .  3.250D0,
     .  3.500D0,
     .  3.750D0,
     .  4.000D0,
     .  4.250D0,
     .  4.500D0,
     .  4.750D0,
     .  5.000D0,
     .  5.250D0,
     .  5.500D0,
     .  5.750D0,
     .  6.000D0,
     .  6.250D0,
     .  6.500D0,
     .  6.750D0,
     .  7.000D0,
     .  7.250D0,
     .  7.500D0,
     .  7.750D0,
     .  8.000D0,
     .  8.250D0,
     .  8.500D0,
     .  8.750D0,
     .  9.000D0,
     .  9.250D0,
     .  9.500D0,
     .  9.750D0,
     . 10.000D0,
     . 11.000D0,
     . 12.000D0,
     . 13.000D0,
     . 14.000D0,
     . 15.000D0,
     . 16.000D0,
     . 18.000D0,
     . 20.000D0,
     . 22.000D0,
     . 24.000D0/
      DATA VCOEF/
     *    0.107248616232689743D+11,
     *   -0.151265488774406090D+11,
     *    0.554393198131965065D+10,
     *   -0.254096323252676010D+10,
     *    0.599113896929283237D+10,
     *   -0.437682982148953724D+10,
     *    0.938084542288794041D+09,
     *   -0.882719321144740224D+09,
     *   -0.255421231587736726D+09,
     *   -0.433957613197937131D+09,
     *   -0.314612507450466037D+09,
     *   -0.279648523522281528D+09,
     *   -0.204317749131337285D+09,
     *   -0.129803369299804091D+09,
     *   -0.492359249163770080D+08,
     *    0.232679882303901613D+08,
     *    0.738602284505146593D+08,
     *    0.110542353359060466D+09,
     *    0.130827628675404668D+09,
     *    0.144026009177417278D+09,
     *    0.146117157966813326D+09,
     *    0.139911806130347729D+09,
     *    0.127444399539057970D+09,
     *    0.116533526609477043D+09,
     *    0.101267218349838257D+09,
     *    0.763210883155441284D+08,
     *    0.789265001109118462D+08,
     *    0.130536266559424400D+08,
     *    0.939625543636751175D+08,
     *    0.571484442999901772D+08,
     *   -0.174761614689388275D+08,
     *   -0.769233069229841232D+07,
     *   -0.919684397630071640D+07,
     *   -0.145359583360872269D+08,
     *    0.493616831508922577D+07,
     *    0.105739155455636978D+08,
     *   -0.446885001644976139D+08,
     *    0.129884371708338261D+08,
     *    0.721092130824542046D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_2=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN LEGENDRE
C ANGULAR EXPANSION OF Kr-OH(X,r=re) VSUM RCCSD(T) PES
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_3(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=39)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .  3.000D0,
     .  3.250D0,
     .  3.500D0,
     .  3.750D0,
     .  4.000D0,
     .  4.250D0,
     .  4.500D0,
     .  4.750D0,
     .  5.000D0,
     .  5.250D0,
     .  5.500D0,
     .  5.750D0,
     .  6.000D0,
     .  6.250D0,
     .  6.500D0,
     .  6.750D0,
     .  7.000D0,
     .  7.250D0,
     .  7.500D0,
     .  7.750D0,
     .  8.000D0,
     .  8.250D0,
     .  8.500D0,
     .  8.750D0,
     .  9.000D0,
     .  9.250D0,
     .  9.500D0,
     .  9.750D0,
     . 10.000D0,
     . 11.000D0,
     . 12.000D0,
     . 13.000D0,
     . 14.000D0,
     . 15.000D0,
     . 16.000D0,
     . 18.000D0,
     . 20.000D0,
     . 22.000D0,
     . 24.000D0/
      DATA VCOEF/
     *    0.747340001091994858D+10,
     *   -0.906788615501340294D+10,
     *    0.403229945860146666D+10,
     *   -0.700927250751128769D+10,
     *    0.800193912476531792D+10,
     *   -0.492893505859770012D+10,
     *    0.203273167924659586D+10,
     *   -0.594544118656662703D+09,
     *    0.141771414099769026D+09,
     *   -0.193068130343837589D+09,
     *   -0.119164822517496258D+09,
     *   -0.150350360126012057D+09,
     *   -0.124814156767583311D+09,
     *   -0.104868646792289495D+09,
     *   -0.663103076369976997D+08,
     *   -0.169964318346613050D+08,
     *    0.141727491604927778D+08,
     *    0.291611809366096258D+08,
     *    0.477237506743898392D+08,
     *    0.640872731078896523D+08,
     *    0.767548427884361744D+08,
     *    0.751696422735466957D+08,
     *    0.766909827815840244D+08,
     *    0.592605607167608738D+08,
     *    0.657616346662329435D+08,
     *    0.414274715107479095D+08,
     *    0.524965259351714849D+08,
     *   -0.136584309700131416D+07,
     *    0.806168255410342216D+08,
     *    0.464196389842672348D+08,
     *   -0.145093818512091637D+08,
     *    0.108393249806747437D+08,
     *   -0.146321587884035110D+08,
     *   -0.326054832617044449D+07,
     *   -0.235684829232215881D+07,
     *   -0.485475408069992065D+07,
     *    0.565222289511036873D+07,
     *   -0.250208376553955078D+08,
     *    0.175404615570994616D+08/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_3=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN LEGENDRE
C ANGULAR EXPANSION OF Kr-OH(X,r=re) VSUM RCCSD(T) PES
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_4(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=39)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .  3.000D0,
     .  3.250D0,
     .  3.500D0,
     .  3.750D0,
     .  4.000D0,
     .  4.250D0,
     .  4.500D0,
     .  4.750D0,
     .  5.000D0,
     .  5.250D0,
     .  5.500D0,
     .  5.750D0,
     .  6.000D0,
     .  6.250D0,
     .  6.500D0,
     .  6.750D0,
     .  7.000D0,
     .  7.250D0,
     .  7.500D0,
     .  7.750D0,
     .  8.000D0,
     .  8.250D0,
     .  8.500D0,
     .  8.750D0,
     .  9.000D0,
     .  9.250D0,
     .  9.500D0,
     .  9.750D0,
     . 10.000D0,
     . 11.000D0,
     . 12.000D0,
     . 13.000D0,
     . 14.000D0,
     . 15.000D0,
     . 16.000D0,
     . 18.000D0,
     . 20.000D0,
     . 22.000D0,
     . 24.000D0/
      DATA VCOEF/
     *    0.790948467278751945D+10,
     *   -0.110747161433036842D+11,
     *    0.494257303834481812D+10,
     *   -0.668695157548302078D+10,
     *    0.104094688864126396D+11,
     *   -0.807780904218128300D+10,
     *    0.262746739000534821D+10,
     *   -0.801437049974289179D+09,
     *    0.395746134043645978D+09,
     *   -0.377339063407847285D+08,
     *    0.739278678168736696D+08,
     *    0.216972166859626770D+07,
     *    0.269987883487930894D+07,
     *   -0.711072983071462810D+07,
     *    0.275552285212099552D+06,
     *    0.162949794473858476D+08,
     *    0.208446113580332398D+08,
     *    0.196563718403483629D+08,
     *    0.230353203017060757D+08,
     *    0.280315382713446617D+08,
     *    0.330677084368937016D+08,
     *    0.345803872117404938D+08,
     *    0.378378256641170979D+08,
     *    0.206659178498990536D+08,
     *    0.280078751822465658D+08,
     *    0.218637214182413816D+08,
     *    0.426347851343953609D+07,
     *    0.127464475397220850D+08,
     *    0.219563858710591793D+08,
     *    0.229215320688014030D+08,
     *   -0.982757721016359329D+07,
     *   -0.245463397744560242D+07,
     *   -0.153398414276695251D+07,
     *   -0.715553458886051178D+07,
     *    0.281788770740461349D+07,
     *    0.190332698167848587D+07,
     *   -0.918741874832797050D+07,
     *    0.767778018621706963D+07,
     *   -0.543716424407398701D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_4=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN LEGENDRE
C ANGULAR EXPANSION OF Kr-OH(X,r=re) VSUM RCCSD(T) PES
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_5(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=39)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .  3.000D0,
     .  3.250D0,
     .  3.500D0,
     .  3.750D0,
     .  4.000D0,
     .  4.250D0,
     .  4.500D0,
     .  4.750D0,
     .  5.000D0,
     .  5.250D0,
     .  5.500D0,
     .  5.750D0,
     .  6.000D0,
     .  6.250D0,
     .  6.500D0,
     .  6.750D0,
     .  7.000D0,
     .  7.250D0,
     .  7.500D0,
     .  7.750D0,
     .  8.000D0,
     .  8.250D0,
     .  8.500D0,
     .  8.750D0,
     .  9.000D0,
     .  9.250D0,
     .  9.500D0,
     .  9.750D0,
     . 10.000D0,
     . 11.000D0,
     . 12.000D0,
     . 13.000D0,
     . 14.000D0,
     . 15.000D0,
     . 16.000D0,
     . 18.000D0,
     . 20.000D0,
     . 22.000D0,
     . 24.000D0/
      DATA VCOEF/
     *    0.634862155591963005D+10,
     *   -0.101606720132869797D+11,
     *    0.731332854276379585D+10,
     *   -0.116038672657824326D+11,
     *    0.149212394213599052D+11,
     *   -0.977523481196464920D+10,
     *    0.352090067702706528D+10,
     *   -0.113352118507166314D+10,
     *    0.425523286447757721D+09,
     *   -0.106074124364168286D+09,
     *    0.686437231487681866D+08,
     *    0.233760721902692318D+07,
     *    0.186089739275537953D+08,
     *    0.356705531473678350D+07,
     *    0.581991629367935658D+07,
     *    0.206910965367119312D+08,
     *    0.199613628168695569D+08,
     *    0.106169757638845444D+08,
     *    0.756596152483141422D+07,
     *    0.911996797011852264D+07,
     *    0.133602872204751968D+08,
     *    0.191395810804677010D+08,
     *    0.922440818934190273D+07,
     *    0.197336684892741442D+08,
     *    0.144946624238115549D+07,
     *    0.115116947406084910D+08,
     *    0.408197793913128972D+07,
     *   -0.954611120195919275D+07,
     *    0.193771868726526499D+08,
     *    0.312932079893636703D+07,
     *   -0.799431869456291199D+07,
     *    0.894359325200366974D+07,
     *   -0.496325500638389587D+07,
     *   -0.364133218446540833D+07,
     *   -0.368906735531806946D+06,
     *    0.513249838176345825D+07,
     *   -0.123198399321360588D+08,
     *    0.119410431372797489D+08,
     *   -0.534182341556429863D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_5=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN LEGENDRE
C ANGULAR EXPANSION OF Kr-OH(X,r=re) VSUM RCCSD(T) PES
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_6(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=39)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .  3.000D0,
     .  3.250D0,
     .  3.500D0,
     .  3.750D0,
     .  4.000D0,
     .  4.250D0,
     .  4.500D0,
     .  4.750D0,
     .  5.000D0,
     .  5.250D0,
     .  5.500D0,
     .  5.750D0,
     .  6.000D0,
     .  6.250D0,
     .  6.500D0,
     .  6.750D0,
     .  7.000D0,
     .  7.250D0,
     .  7.500D0,
     .  7.750D0,
     .  8.000D0,
     .  8.250D0,
     .  8.500D0,
     .  8.750D0,
     .  9.000D0,
     .  9.250D0,
     .  9.500D0,
     .  9.750D0,
     . 10.000D0,
     . 11.000D0,
     . 12.000D0,
     . 13.000D0,
     . 14.000D0,
     . 15.000D0,
     . 16.000D0,
     . 18.000D0,
     . 20.000D0,
     . 22.000D0,
     . 24.000D0/
      DATA VCOEF/
     *    0.211459165138363695D+10,
     *   -0.277667784015477371D+10,
     *    0.567390147123888779D+10,
     *   -0.162963424282780418D+11,
     *    0.218107113095778275D+11,
     *   -0.152523046118422394D+11,
     *    0.556233713737663841D+10,
     *   -0.164105841099759245D+10,
     *    0.715084228244249225D+09,
     *   -0.141301516626129210D+09,
     *    0.110794094239706278D+09,
     *    0.536013305357372761D+07,
     *    0.273684002277901694D+08,
     *    0.809750701315307617D+07,
     *    0.841480321764039993D+07,
     *    0.186199227437106371D+08,
     *    0.114774545782145262D+08,
     *    0.344989970779275894D+07,
     *    0.254354619241905212D+07,
     *    0.323102589848589897D+07,
     *    0.451512502846920490D+07,
     *    0.110846934804904461D+08,
     *    0.402589237529695034D+06,
     *    0.133839895076580048D+08,
     *   -0.332883927162909508D+07,
     *    0.232075844002056122D+07,
     *    0.709380135325622559D+07,
     *   -0.610590898336791992D+07,
     *    0.237564336565232277D+07,
     *    0.636640621435117722D+07,
     *   -0.652782204353523254D+07,
     *    0.338269531751632690D+06,
     *    0.321533226963043213D+06,
     *    0.122451898520374298D+07,
     *    0.278705855759620667D+06,
     *   -0.348987968663406372D+07,
     *    0.443365649380970001D+07,
     *   -0.753888166499280930D+07,
     *    0.522519271144866943D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_6=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN LEGENDRE
C ANGULAR EXPANSION OF Kr-OH(X,r=re) VSUM RCCSD(T) PES
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_7(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=39)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .  3.000D0,
     .  3.250D0,
     .  3.500D0,
     .  3.750D0,
     .  4.000D0,
     .  4.250D0,
     .  4.500D0,
     .  4.750D0,
     .  5.000D0,
     .  5.250D0,
     .  5.500D0,
     .  5.750D0,
     .  6.000D0,
     .  6.250D0,
     .  6.500D0,
     .  6.750D0,
     .  7.000D0,
     .  7.250D0,
     .  7.500D0,
     .  7.750D0,
     .  8.000D0,
     .  8.250D0,
     .  8.500D0,
     .  8.750D0,
     .  9.000D0,
     .  9.250D0,
     .  9.500D0,
     .  9.750D0,
     . 10.000D0,
     . 11.000D0,
     . 12.000D0,
     . 13.000D0,
     . 14.000D0,
     . 15.000D0,
     . 16.000D0,
     . 18.000D0,
     . 20.000D0,
     . 22.000D0,
     . 24.000D0/
      DATA VCOEF/
     *    0.595818783992917389D+08,
     *   -0.757049339965861797D+09,
     *    0.554581323014209652D+10,
     *   -0.196399049937673454D+11,
     *    0.307429031577237625D+11,
     *   -0.217367435793850861D+11,
     *    0.745118342227116585D+10,
     *   -0.248924471881365299D+10,
     *    0.900105542937775850D+09,
     *   -0.254051345048080683D+09,
     *    0.119017275648574114D+09,
     *   -0.121792034964179993D+08,
     *    0.287053694049052000D+08,
     *    0.419901966822576523D+07,
     *    0.202314451997256279D+07,
     *    0.131170985101206303D+08,
     *    0.804399884903240204D+07,
     *    0.141072739501690865D+07,
     *    0.922321982784271240D+06,
     *    0.100748971696567535D+07,
     *    0.109254268729138374D+07,
     *    0.198702963756448030D+07,
     *    0.628717153744673729D+07,
     *   -0.747529705280542374D+06,
     *    0.175903761447048187D+07,
     *   -0.304385715439081192D+07,
     *    0.665881251525497437D+07,
     *   -0.463833267470312119D+07,
     *    0.414458773277616501D+07,
     *   -0.303554805325520039D+07,
     *    0.208211842515659332D+07,
     *   -0.269528382381820679D+07,
     *    0.588800446224594116D+07,
     *   -0.112066152963256836D+08,
     *    0.979975830324840546D+07,
     *   -0.455968967960643768D+07,
     *    0.192190851966094971D+07,
     *   -0.739432658108234406D+06,
     *    0.192409933963775635D+06/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_7=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN LEGENDRE
C ANGULAR EXPANSION OF Kr-OH(X,r=re) VSUM RCCSD(T) PES
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_8(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=39)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .  3.000D0,
     .  3.250D0,
     .  3.500D0,
     .  3.750D0,
     .  4.000D0,
     .  4.250D0,
     .  4.500D0,
     .  4.750D0,
     .  5.000D0,
     .  5.250D0,
     .  5.500D0,
     .  5.750D0,
     .  6.000D0,
     .  6.250D0,
     .  6.500D0,
     .  6.750D0,
     .  7.000D0,
     .  7.250D0,
     .  7.500D0,
     .  7.750D0,
     .  8.000D0,
     .  8.250D0,
     .  8.500D0,
     .  8.750D0,
     .  9.000D0,
     .  9.250D0,
     .  9.500D0,
     .  9.750D0,
     . 10.000D0,
     . 11.000D0,
     . 12.000D0,
     . 13.000D0,
     . 14.000D0,
     . 15.000D0,
     . 16.000D0,
     . 18.000D0,
     . 20.000D0,
     . 22.000D0,
     . 24.000D0/
      DATA VCOEF/
     *   -0.735737155481506348D+09,
     *    0.365498076925222993D+09,
     *    0.774575735424332619D+10,
     *   -0.221631425477595215D+11,
     *    0.281567971064589233D+11,
     *   -0.188113495253790779D+11,
     *    0.690400189389943409D+10,
     *   -0.216942075707533932D+10,
     *    0.816197507118257523D+09,
     *   -0.232153234567023516D+09,
     *    0.978275868398387432D+08,
     *   -0.153465073651915789D+08,
     *    0.235006509776092619D+08,
     *    0.256866745871543884D+06,
     *    0.875281513861656189D+06,
     *    0.100782181370038986D+08,
     *    0.832286400806427002D+06,
     *    0.118549461744451523D+07,
     *   -0.864278815328836441D+06,
     *    0.599738151105642319D+05,
     *    0.418097155445015430D+07,
     *   -0.504603484314996004D+07,
     *    0.516617191177153587D+07,
     *   -0.156360344460964203D+07,
     *    0.840511159457445145D+07,
     *   -0.222879726779298782D+08,
     *    0.339736071936345100D+08,
     *   -0.257900580814843178D+08,
     *    0.990517060397386551D+07,
     *   -0.121757752096700668D+07,
     *   -0.333432123307132721D+07,
     *    0.521524146025371552D+07,
     *   -0.367352287057209015D+07,
     *    0.247296227717494965D+07,
     *   -0.145482403345489502D+07,
     *    0.277642441337585449D+06,
     *   -0.117026248213291168D+06,
     *    0.450245314059257507D+05,
     *   -0.117159652369022369D+05/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_8=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C COEFFICIENTS IN LEGENDRE
C ANGULAR EXPANSION OF Kr-OH(X,r=re) VSUM RCCSD(T) PES
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without 
C contacting me before.  
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_9(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=39)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .  3.000D0,
     .  3.250D0,
     .  3.500D0,
     .  3.750D0,
     .  4.000D0,
     .  4.250D0,
     .  4.500D0,
     .  4.750D0,
     .  5.000D0,
     .  5.250D0,
     .  5.500D0,
     .  5.750D0,
     .  6.000D0,
     .  6.250D0,
     .  6.500D0,
     .  6.750D0,
     .  7.000D0,
     .  7.250D0,
     .  7.500D0,
     .  7.750D0,
     .  8.000D0,
     .  8.250D0,
     .  8.500D0,
     .  8.750D0,
     .  9.000D0,
     .  9.250D0,
     .  9.500D0,
     .  9.750D0,
     . 10.000D0,
     . 11.000D0,
     . 12.000D0,
     . 13.000D0,
     . 14.000D0,
     . 15.000D0,
     . 16.000D0,
     . 18.000D0,
     . 20.000D0,
     . 22.000D0,
     . 24.000D0/
      DATA VCOEF/
     *   -0.108575860545568800D+10,
     *    0.772296987666342020D+09,
     *    0.394585919458055973D+10,
     *   -0.146667628108104019D+11,
     *    0.234271943850725479D+11,
     *   -0.167528431153319359D+11,
     *    0.571658083453688526D+10,
     *   -0.191633821395803261D+10,
     *    0.682624663659066439D+09,
     *   -0.207895511212280989D+09,
     *    0.780367356473298073D+08,
     *   -0.184135795539050102D+08,
     *    0.173776406149927378D+08,
     *    0.125901462086534500D+07,
     *   -0.869333591265320778D+06,
     *    0.476010804544663429D+07,
     *    0.105027950483489037D+07,
     *   -0.102438519457435608D+07,
     *    0.162689474926543236D+07,
     *    0.864935831467390060D+06,
     *   -0.449743886288166046D+06,
     *   -0.723862235954403877D+05,
     *   -0.434064651839733124D+06,
     *    0.152996886946702003D+07,
     *   -0.454287486869335175D+06,
     *   -0.337242590674495697D+07,
     *    0.655415799069023132D+07,
     *   -0.447825728284454346D+07,
     *    0.240913883853673935D+07,
     *   -0.417810386234462261D+06,
     *   -0.120459433232259750D+07,
     *   -0.481273151935577393D+06,
     *    0.330200406666564941D+07,
     *   -0.274130883965682983D+07,
     *    0.790368320182800293D+06,
     *   -0.150835964381217957D+06,
     *    0.635773363690376282D+05,
     *   -0.244606644482612610D+05,
     *    0.636498074150085449D+04/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_9=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2   
C COEFFICIENTS IN LEGENDRE PL2(x)
C ANGULAR EXPANSION OF Kr-OH(X,r=re) VDIF RCCSD(T) PES
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without
C contacting me before.
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS

      FUNCTION VL2_10(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=39)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .  3.000D0,
     .  3.250D0,
     .  3.500D0,
     .  3.750D0,
     .  4.000D0,
     .  4.250D0,
     .  4.500D0,
     .  4.750D0,
     .  5.000D0,
     .  5.250D0,
     .  5.500D0,
     .  5.750D0,
     .  6.000D0,
     .  6.250D0,
     .  6.500D0,
     .  6.750D0,
     .  7.000D0,
     .  7.250D0,
     .  7.500D0,
     .  7.750D0,
     .  8.000D0,
     .  8.250D0,
     .  8.500D0,
     .  8.750D0,
     .  9.000D0,
     .  9.250D0,
     .  9.500D0,
     .  9.750D0,
     . 10.000D0,
     . 11.000D0,
     . 12.000D0,
     . 13.000D0,
     . 14.000D0,
     . 15.000D0,
     . 16.000D0,
     . 18.000D0,
     . 20.000D0,
     . 22.000D0,
     . 24.000D0/
      DATA VCOEF/
     *   -0.148233706161373943D+09,
     *    0.524558053023571849D+09,
     *    0.826357135361385942D+09,
     *   -0.204812328062248182D+10,
     *    0.403377359623183727D+09,
     *    0.335551641672020793D+09,
     *    0.784234646639250815D+08,
     *    0.219527043773946241D+08,
     *    0.437700396152780950D+07,
     *    0.716866979242026806D+06,
     *    0.150607086232636124D+07,
     *   -0.576661662327349186D+06,
     *   -0.736304186484217644D+05,
     *    0.271306331577599049D+06,
     *   -0.196024997001886368D+06,
     *   -0.411172920289903879D+06,
     *    0.129667437186788023D+07,
     *   -0.124757220953269303D+07,
     *    0.948351362592935562D+06,
     *   -0.621224460508257151D+06,
     *    0.191753425403892994D+06,
     *   -0.589340072635412216D+05,
     *    0.180396303763985634D+05,
     *   -0.550083478543162346D+04,
     *    0.167130362761020660D+04,
     *   -0.505961287498474121D+03,
     *    0.152277967885136604D+03,
     *   -0.439992869794368744D+02,
     *    0.943406204134225845D+01,
     *   -0.568265765905380249D+00,
     *    0.222721278667449951D+00,
     *   -0.849287509918212891D-01,
     *    0.317519903182983398D-01,
     *   -0.107471048831939697D-01,
     *    0.253668427467346191D-02,
     *   -0.615492463111877441D-03,
     *    0.144541263580322266D-03,
     *    0.558570027351379395D-04,
     *    0.190660357475280762D-04/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_10=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2   
C COEFFICIENTS IN LEGENDRE PL2(x)
C ANGULAR EXPANSION OF Kr-OH(X,r=re) VDIF RCCSD(T) PES
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without
C contacting me before.
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS

      FUNCTION VL2_2(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=39)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .  3.000D0,
     .  3.250D0,
     .  3.500D0,
     .  3.750D0,
     .  4.000D0,
     .  4.250D0,
     .  4.500D0,
     .  4.750D0,
     .  5.000D0,
     .  5.250D0,
     .  5.500D0,
     .  5.750D0,
     .  6.000D0,
     .  6.250D0,
     .  6.500D0,
     .  6.750D0,
     .  7.000D0,
     .  7.250D0,
     .  7.500D0,
     .  7.750D0,
     .  8.000D0,
     .  8.250D0,
     .  8.500D0,
     .  8.750D0,
     .  9.000D0,
     .  9.250D0,
     .  9.500D0,
     .  9.750D0,
     . 10.000D0,
     . 11.000D0,
     . 12.000D0,
     . 13.000D0,
     . 14.000D0,
     . 15.000D0,
     . 16.000D0,
     . 18.000D0,
     . 20.000D0,
     . 22.000D0,
     . 24.000D0/
      DATA VCOEF/
     *   -0.238442658022454929D+10,
     *    0.734788013428048515D+10,
     *   -0.309172174349280834D+10,
     *    0.601736457669748783D+09,
     *   -0.109670078630409288D+10,
     *   -0.291165905530556917D+09,
     *   -0.383148337411976755D+09,
     *   -0.232047219098335713D+09,
     *   -0.227027672576747149D+09,
     *   -0.202192871234510422D+09,
     *   -0.193875344813789606D+09,
     *   -0.177686350698948711D+09,
     *   -0.161100002361591160D+09,
     *   -0.141828891429747194D+09,
     *   -0.118152714462816775D+09,
     *   -0.884219918993771076D+08,
     *   -0.543440370585142970D+08,
     *   -0.192015624138791859D+08,
     *    0.131735629921401441D+08,
     *    0.366130494041590542D+08,
     *    0.563885897666253448D+08,
     *    0.645943928121656179D+08,
     *    0.733087015698451996D+08,
     *    0.797106566291652918D+08,
     *    0.793421081920348406D+08,
     *    0.754541743505449295D+08,
     *    0.853884029569187164D+08,
     *   -0.882085556148052216D+06,
     *    0.172085728750956059D+09,
     *    0.141951966992311954D+09,
     *    0.343861107772841454D+08,
     *    0.116382065095996857D+08,
     *    0.913989955182719231D+07,
     *    0.598729677242231369D+07,
     *   -0.274024321573605537D+08,
     *    0.235744008489114046D+07,
     *    0.133428062872217000D+08,
     *    0.134109317120519280D+07,
     *   -0.197877845514203459D+08/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_2=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2   
C COEFFICIENTS IN LEGENDRE PL2(x)
C ANGULAR EXPANSION OF Kr-OH(X,r=re) VDIF RCCSD(T) PES
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without
C contacting me before.
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS

      FUNCTION VL2_3(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=39)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .  3.000D0,
     .  3.250D0,
     .  3.500D0,
     .  3.750D0,
     .  4.000D0,
     .  4.250D0,
     .  4.500D0,
     .  4.750D0,
     .  5.000D0,
     .  5.250D0,
     .  5.500D0,
     .  5.750D0,
     .  6.000D0,
     .  6.250D0,
     .  6.500D0,
     .  6.750D0,
     .  7.000D0,
     .  7.250D0,
     .  7.500D0,
     .  7.750D0,
     .  8.000D0,
     .  8.250D0,
     .  8.500D0,
     .  8.750D0,
     .  9.000D0,
     .  9.250D0,
     .  9.500D0,
     .  9.750D0,
     . 10.000D0,
     . 11.000D0,
     . 12.000D0,
     . 13.000D0,
     . 14.000D0,
     . 15.000D0,
     . 16.000D0,
     . 18.000D0,
     . 20.000D0,
     . 22.000D0,
     . 24.000D0/
      DATA VCOEF/
     *    0.227706289167344666D+10,
     *   -0.492168545351502228D+10,
     *    0.231549816311410713D+10,
     *   -0.528303832765132129D+09,
     *    0.897512944141575336D+09,
     *    0.704145149454316050D+08,
     *    0.124347863082717091D+09,
     *   -0.479271772531816363D+07,
     *    0.129530405228869617D+07,
     *   -0.123949733941587657D+08,
     *   -0.124100039159778059D+08,
     *   -0.155111747648834586D+08,
     *   -0.211767119147433490D+08,
     *   -0.254671202550031394D+08,
     *   -0.280176211744493544D+08,
     *   -0.276601203910383880D+08,
     *   -0.267825009223000109D+08,
     *   -0.242464666381783038D+08,
     *   -0.218476838756548502D+08,
     *   -0.160952802759158164D+08,
     *   -0.145216166052054465D+08,
     *   -0.937393046739017963D+07,
     *   -0.399845034557414055D+07,
     *   -0.391510680922961235D+07,
     *    0.126290371943724155D+07,
     *   -0.149895284766542912D+07,
     *    0.664750963974857330D+07,
     *   -0.763895619446098804D+07,
     *    0.161339692421791553D+08,
     *    0.112158463832218647D+08,
     *    0.200258587914943695D+06,
     *    0.700629662154722214D+07,
     *    0.456370910092735291D+07,
     *   -0.869250252006399632D+07,
     *    0.150110611888504028D+07,
     *    0.648049040597647429D+07,
     *   -0.103942359132092297D+08,
     *    0.583313851878318191D+07,
     *   -0.135094610610172153D+06/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_3=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2   
C COEFFICIENTS IN LEGENDRE PL2(x)
C ANGULAR EXPANSION OF Kr-OH(X,r=re) VDIF RCCSD(T) PES
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without
C contacting me before.
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS

      FUNCTION VL2_4(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=39)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .  3.000D0,
     .  3.250D0,
     .  3.500D0,
     .  3.750D0,
     .  4.000D0,
     .  4.250D0,
     .  4.500D0,
     .  4.750D0,
     .  5.000D0,
     .  5.250D0,
     .  5.500D0,
     .  5.750D0,
     .  6.000D0,
     .  6.250D0,
     .  6.500D0,
     .  6.750D0,
     .  7.000D0,
     .  7.250D0,
     .  7.500D0,
     .  7.750D0,
     .  8.000D0,
     .  8.250D0,
     .  8.500D0,
     .  8.750D0,
     .  9.000D0,
     .  9.250D0,
     .  9.500D0,
     .  9.750D0,
     . 10.000D0,
     . 11.000D0,
     . 12.000D0,
     . 13.000D0,
     . 14.000D0,
     . 15.000D0,
     . 16.000D0,
     . 18.000D0,
     . 20.000D0,
     . 22.000D0,
     . 24.000D0/
      DATA VCOEF/
     *    0.311064685794629231D+08,
     *    0.117388265636479950D+10,
     *    0.121840301026444769D+10,
     *   -0.247574687936656857D+10,
     *   -0.733191295395989656D+09,
     *   -0.264128101818752401D+08,
     *    0.153217255388659418D+09,
     *    0.190457693750714600D+09,
     *    0.161197048584909260D+09,
     *    0.123210204992129624D+09,
     *    0.885363743259236217D+08,
     *    0.581353117859297991D+08,
     *    0.347605971138305664D+08,
     *    0.165644641523060203D+08,
     *    0.452446048090869188D+07,
     *   -0.252921170420813560D+07,
     *   -0.646390618471749499D+07,
     *   -0.852230756052464247D+07,
     *   -0.712024786020284891D+07,
     *   -0.761668612328529358D+07,
     *   -0.455931315277683735D+07,
     *   -0.443433121139097214D+07,
     *    0.964403379449307919D+06,
     *   -0.413607755568051338D+07,
     *    0.380415234390106797D+07,
     *    0.130716456111666560D+07,
     *    0.956647260039553046D+06,
     *    0.286747679755660892D+07,
     *    0.502420183310452662D+07,
     *    0.746382095096322894D+07,
     *   -0.852775237643957138D+06,
     *    0.721639288433274627D+07,
     *   -0.469611572901654243D+07,
     *    0.714261202100783587D+07,
     *   -0.826190911101830006D+07,
     *    0.770028547871826589D+07,
     *   -0.628678971365910769D+07,
     *   -0.371330830430652946D+07,
     *    0.737284777487060800D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_4=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2   
C COEFFICIENTS IN LEGENDRE PL2(x)
C ANGULAR EXPANSION OF Kr-OH(X,r=re) VDIF RCCSD(T) PES
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without
C contacting me before.
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS

      FUNCTION VL2_5(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=39)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .  3.000D0,
     .  3.250D0,
     .  3.500D0,
     .  3.750D0,
     .  4.000D0,
     .  4.250D0,
     .  4.500D0,
     .  4.750D0,
     .  5.000D0,
     .  5.250D0,
     .  5.500D0,
     .  5.750D0,
     .  6.000D0,
     .  6.250D0,
     .  6.500D0,
     .  6.750D0,
     .  7.000D0,
     .  7.250D0,
     .  7.500D0,
     .  7.750D0,
     .  8.000D0,
     .  8.250D0,
     .  8.500D0,
     .  8.750D0,
     .  9.000D0,
     .  9.250D0,
     .  9.500D0,
     .  9.750D0,
     . 10.000D0,
     . 11.000D0,
     . 12.000D0,
     . 13.000D0,
     . 14.000D0,
     . 15.000D0,
     . 16.000D0,
     . 18.000D0,
     . 20.000D0,
     . 22.000D0,
     . 24.000D0/
      DATA VCOEF/
     *   -0.950560052971527696D+09,
     *    0.246175609542858601D+10,
     *   -0.335897863766864014D+10,
     *    0.217625590367156601D+10,
     *    0.537346527591993064D+08,
     *   -0.490266843262237459D+08,
     *   -0.172856147039380729D+09,
     *   -0.863041290536779016D+08,
     *   -0.513931202600138336D+08,
     *   -0.185914860249467939D+08,
     *   -0.445168233025008440D+07,
     *    0.144591708701276779D+07,
     *    0.462875038025569916D+07,
     *    0.229199681891128421D+07,
     *    0.128103163063663384D+07,
     *    0.252547655400007963D+06,
     *   -0.242415178226613998D+07,
     *   -0.472973615728676319D+06,
     *   -0.405538038472408056D+07,
     *   -0.897624480593204498D+06,
     *   -0.343860001917552948D+07,
     *    0.107258573143434525D+07,
     *   -0.563191802089100704D+07,
     *    0.327940752903138101D+07,
     *   -0.278654247641196847D+07,
     *    0.253761373237431049D+06,
     *    0.299761126169377565D+07,
     *   -0.336000628909349442D+07,
     *    0.306753348415863514D+07,
     *    0.198293934589147568D+06,
     *    0.213413456926155090D+07,
     *    0.266914328844904900D+06,
     *   -0.728077796970963478D+06,
     *    0.150444188519847393D+07,
     *   -0.158738044400352240D+07,
     *    0.623831111377179623D+07,
     *   -0.161530941655363441D+08,
     *    0.204400933133132756D+08,
     *   -0.103677548785188645D+08/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_5=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2   
C COEFFICIENTS IN LEGENDRE PL2(x)
C ANGULAR EXPANSION OF Kr-OH(X,r=re) VDIF RCCSD(T) PES
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without
C contacting me before.
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS

      FUNCTION VL2_6(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=39)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .  3.000D0,
     .  3.250D0,
     .  3.500D0,
     .  3.750D0,
     .  4.000D0,
     .  4.250D0,
     .  4.500D0,
     .  4.750D0,
     .  5.000D0,
     .  5.250D0,
     .  5.500D0,
     .  5.750D0,
     .  6.000D0,
     .  6.250D0,
     .  6.500D0,
     .  6.750D0,
     .  7.000D0,
     .  7.250D0,
     .  7.500D0,
     .  7.750D0,
     .  8.000D0,
     .  8.250D0,
     .  8.500D0,
     .  8.750D0,
     .  9.000D0,
     .  9.250D0,
     .  9.500D0,
     .  9.750D0,
     . 10.000D0,
     . 11.000D0,
     . 12.000D0,
     . 13.000D0,
     . 14.000D0,
     . 15.000D0,
     . 16.000D0,
     . 18.000D0,
     . 20.000D0,
     . 22.000D0,
     . 24.000D0/
      DATA VCOEF/
     *    0.597834064734796286D+09,
     *   -0.574404904531731129D+09,
     *    0.347093583646252441D+10,
     *   -0.512456233753853226D+10,
     *    0.356737639903022408D+09,
     *    0.606855252772457600D+09,
     *    0.362867807385929167D+09,
     *    0.164994393042727351D+09,
     *    0.766002623069429398D+08,
     *    0.341562703329478353D+08,
     *    0.159623862470363379D+08,
     *    0.839859165975689888D+07,
     *    0.415783437650048733D+07,
     *    0.235814253926801682D+07,
     *    0.642818657255768776D+06,
     *    0.294513842133134604D+06,
     *    0.115061818026527762D+06,
     *   -0.638248799156546593D+06,
     *   -0.122011628325676918D+07,
     *   -0.109877149844312668D+07,
     *   -0.701680111161470413D+06,
     *    0.413394505396842957D+06,
     *   -0.450979057692337036D+07,
     *    0.689273847724455595D+07,
     *   -0.688372644974449277D+07,
     *    0.157869783029861562D+07,
     *    0.152513190719872713D+07,
     *   -0.182463815380281210D+07,
     *    0.263410369994658232D+07,
     *   -0.231371977934134007D+07,
     *    0.303128584358513355D+07,
     *   -0.233232065290617943D+07,
     *    0.823594227469003201D+07,
     *   -0.182639204802335501D+08,
     *    0.146742532000879049D+08,
     *    0.260186656977236271D+06,
     *   -0.147705997190236449D+08,
     *    0.199422016477705836D+08,
     *   -0.862763226192814112D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_6=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2   
C COEFFICIENTS IN LEGENDRE PL2(x)
C ANGULAR EXPANSION OF Kr-OH(X,r=re) VDIF RCCSD(T) PES
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without
C contacting me before.
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS

      FUNCTION VL2_7(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=39)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .  3.000D0,
     .  3.250D0,
     .  3.500D0,
     .  3.750D0,
     .  4.000D0,
     .  4.250D0,
     .  4.500D0,
     .  4.750D0,
     .  5.000D0,
     .  5.250D0,
     .  5.500D0,
     .  5.750D0,
     .  6.000D0,
     .  6.250D0,
     .  6.500D0,
     .  6.750D0,
     .  7.000D0,
     .  7.250D0,
     .  7.500D0,
     .  7.750D0,
     .  8.000D0,
     .  8.250D0,
     .  8.500D0,
     .  8.750D0,
     .  9.000D0,
     .  9.250D0,
     .  9.500D0,
     .  9.750D0,
     . 10.000D0,
     . 11.000D0,
     . 12.000D0,
     . 13.000D0,
     . 14.000D0,
     . 15.000D0,
     . 16.000D0,
     . 18.000D0,
     . 20.000D0,
     . 22.000D0,
     . 24.000D0/
      DATA VCOEF/
     *    0.507709960668299675D+09,
     *   -0.808986698711983442D+09,
     *   -0.123319598600527358D+10,
     *    0.215605397105296946D+10,
     *   -0.683600200191270411D+08,
     *   -0.420310824348253667D+09,
     *   -0.995106620528929532D+08,
     *   -0.370689866248770654D+08,
     *   -0.203680070004735887D+07,
     *    0.627528526598460972D+06,
     *    0.236302818598511815D+07,
     *    0.210418506050437689D+07,
     *    0.980474210268020630D+06,
     *    0.869095672365367413D+06,
     *   -0.168615256461381912D+06,
     *   -0.778703651215136051D+05,
     *    0.146273585824665800D+07,
     *   -0.228195864045292139D+07,
     *    0.201380608371037245D+07,
     *   -0.175996965328538418D+07,
     *   -0.375356657347500324D+06,
     *   -0.896291793411016464D+06,
     *    0.256407167122411728D+07,
     *   -0.433587718188565969D+07,
     *    0.312079481139218807D+07,
     *    0.387530632578134537D+06,
     *   -0.198000499693757296D+07,
     *    0.231816601789489388D+06,
     *    0.898183672901822254D+06,
     *   -0.483966561332523823D+06,
     *   -0.239336277479201555D+06,
     *    0.230146304767632484D+07,
     *   -0.290005728531929851D+07,
     *    0.160572027060836554D+07,
     *   -0.375217756875216961D+06,
     *    0.716075421975255013D+05,
     *   -0.301825682303458452D+05,
     *    0.116124031085446477D+05,
     *   -0.302169741436094046D+04/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_7=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2   
C COEFFICIENTS IN LEGENDRE PL2(x)
C ANGULAR EXPANSION OF Kr-OH(X,r=re) VDIF RCCSD(T) PES
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without
C contacting me before.
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS

      FUNCTION VL2_8(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=39)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .  3.000D0,
     .  3.250D0,
     .  3.500D0,
     .  3.750D0,
     .  4.000D0,
     .  4.250D0,
     .  4.500D0,
     .  4.750D0,
     .  5.000D0,
     .  5.250D0,
     .  5.500D0,
     .  5.750D0,
     .  6.000D0,
     .  6.250D0,
     .  6.500D0,
     .  6.750D0,
     .  7.000D0,
     .  7.250D0,
     .  7.500D0,
     .  7.750D0,
     .  8.000D0,
     .  8.250D0,
     .  8.500D0,
     .  8.750D0,
     .  9.000D0,
     .  9.250D0,
     .  9.500D0,
     .  9.750D0,
     . 10.000D0,
     . 11.000D0,
     . 12.000D0,
     . 13.000D0,
     . 14.000D0,
     . 15.000D0,
     . 16.000D0,
     . 18.000D0,
     . 20.000D0,
     . 22.000D0,
     . 24.000D0/
      DATA VCOEF/
     *   -0.114963814468320549D+09,
     *    0.795268078393739104D+09,
     *    0.241505953010832644D+10,
     *   -0.507587151291350937D+10,
     *    0.825808383566712260D+09,
     *    0.793372230981934071D+09,
     *    0.247102599541449308D+09,
     *    0.807251139179649949D+08,
     *    0.212558928937915266D+08,
     *    0.850109112017114460D+07,
     *    0.178669960082298517D+07,
     *    0.147172356397247314D+07,
     *    0.724283558368682861D+06,
     *   -0.220401596828937531D+06,
     *   -0.693016594386100769D+05,
     *    0.504137438805818558D+06,
     *   -0.192661010573804379D+04,
     *    0.463616058325201273D+06,
     *   -0.647747317159175873D+05,
     *   -0.155567702466785908D+07,
     *    0.222481107727611065D+07,
     *   -0.327637248184370995D+07,
     *    0.305335235439455509D+07,
     *   -0.261817519723975658D+07,
     *    0.226549141055756807D+07,
     *   -0.209584579915183783D+07,
     *    0.298790260708409548D+07,
     *   -0.267623539531166852D+07,
     *    0.769581747289836407D+06,
     *   -0.898138696950674057D+05,
     *    0.220385977011680603D+06,
     *   -0.838340502458810806D+05,
     *    0.310509440327882767D+05,
     *   -0.110250568944215775D+05,
     *    0.317871780931949615D+04,
     *   -0.606634947180747986D+03,
     *    0.255696807056665421D+03,
     *   -0.983756519556045532D+02,
     *    0.255984668731689453D+02/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_8=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2   
C COEFFICIENTS IN LEGENDRE PL2(x)
C ANGULAR EXPANSION OF Kr-OH(X,r=re) VDIF RCCSD(T) PES
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD
C Contact: Jacek Klos jklos@umd.edu
C Reference: J. Klos, to be published
C Please do not pass this PES further without
C contacting me before.
C
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS

      FUNCTION VL2_9(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=39)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .  3.000D0,
     .  3.250D0,
     .  3.500D0,
     .  3.750D0,
     .  4.000D0,
     .  4.250D0,
     .  4.500D0,
     .  4.750D0,
     .  5.000D0,
     .  5.250D0,
     .  5.500D0,
     .  5.750D0,
     .  6.000D0,
     .  6.250D0,
     .  6.500D0,
     .  6.750D0,
     .  7.000D0,
     .  7.250D0,
     .  7.500D0,
     .  7.750D0,
     .  8.000D0,
     .  8.250D0,
     .  8.500D0,
     .  8.750D0,
     .  9.000D0,
     .  9.250D0,
     .  9.500D0,
     .  9.750D0,
     . 10.000D0,
     . 11.000D0,
     . 12.000D0,
     . 13.000D0,
     . 14.000D0,
     . 15.000D0,
     . 16.000D0,
     . 18.000D0,
     . 20.000D0,
     . 22.000D0,
     . 24.000D0/
      DATA VCOEF/
     *    0.341592293849776834D+08,
     *    0.177918577087092519D+09,
     *   -0.131120024291889524D+10,
     *    0.163115899307917690D+10,
     *   -0.239828448845629156D+09,
     *   -0.252494358383107871D+09,
     *   -0.409723457105881721D+08,
     *   -0.462409288543490320D+07,
     *    0.327115060558655858D+07,
     *    0.132386580238880962D+07,
     *    0.107459173024300486D+07,
     *   -0.134292678983330727D+06,
     *    0.315810237276136875D+06,
     *    0.734378943487703800D+05,
     *   -0.223966466982364655D+06,
     *    0.204041204646252096D+06,
     *    0.436654086911678314D+05,
     *    0.736787073677957058D+06,
     *   -0.105149265162914991D+07,
     *    0.326057291964858770D+06,
     *   -0.100644142509818077D+06,
     *    0.309322372418045998D+05,
     *   -0.946831993350386620D+04,
     *    0.288718045619130135D+04,
     *   -0.877206655859947205D+03,
     *    0.265563629593700171D+03,
     *   -0.799295763373374939D+02,
     *    0.230994731187820435D+02,
     *   -0.495518079400062561D+01,
     *    0.299023032188415527D+00,
     *   -0.117313742637634277D+00,
     *    0.447872281074523926D-01,
     *   -0.168235301971435547D-01,
     *    0.556522607803344727D-02,
     *   -0.100460648536682129D-02,
     *    0.281333923339843750D-04,
     *    0.199913978576660156D-03,
     *   -0.392660498619079590D-03,
     *    0.172972679138183594D-03/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_9=SUMA
       RETURN
       END

C------RKHS ROUTINES

C------------------------RKHS---------------------------------
C Fortran  Deck with distance-like and angular-like 
C kernels for RKHS interpolation method
C****************************************************
C* By J. KLOS  12/05/2007   University of Maryland  *
C* REVISON: 1.0                                     * 
C* jasztol@gmail.com                                *
C****************************************************
C-------------------------------------------------------------
C FUNCTIONS:
C
C RKHS_DK(X,Y,N,M): Reproducing kernel for  distance-like (DK) 
C                   variables, see Ho, Rabitz Eq(17) of
C                   J. Chem. Phys. 104, 2584 (1996) 
C                   Range = [0, inf]
C
C  INPUT: X,Y-DOUBLE PRECISION
C          N:INTEGER:ORDER OF SMOOTHNESS OF THE FUNCTION
C          M>=0:INTEGER:RECIPROCAL POWER OF THE WEIGHT W(X)=X**(-M)
C  OUTPUT: DOUBLE PRECISION
C -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -            
C RKHS_AK(X,Y,N): Reproducing kernel for  angular-like (AK) 
C                   variables, see Ho, Rabitz Eq(23) of
C                   J. Chem. Phys. 104, 2584 (1996) 
C                    Range = [0, 1]
C    
C  INPUT: X,Y-DOUBLE PRECISION
C          N:INTEGER:ORDER OF SMOOTHNESS OF THE FUNCTION
C  OUPUT: DOUBLE PRECISION
C  HINT: SCALE ANGLES TO [0,1] INTEGRAL FOR ANGLES, FOR EXAMPLE: 
C          X=(1.d0-cosd(ANGLE))/2.d0
C -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
C RKHS_OP(X,Y,N,M):Reproducing kernel for angular like 
C                    variables with a use of orthogonal 
C                    polynomials Plm(L,M,X) where X is 
C                    K=SUM_L=M^N PLM(X)PLM(X')
C                    variable transformed from angles as
C                    x=cos(theta)
C-------------------------------------------------------------
       DOUBLE PRECISION FUNCTION RKHS_DKN2(X,Y,M)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IF (X .LT. Y) THEN
           XB=X
           XU=Y
       ELSE
           XB=Y
           XU=X
       ENDIF

C       M=ABS(M)! FORCE M>=0
       NSQR=4
       XUPOWM=XU**(-(M+1))

C      TERM WITH BETA FUNCTION
C       CALL BETA(DFLOAT(M+1),DFLOAT(N),BETATERM)
C      N=2 case direct expression
       BETATERM=1.D0/DFLOAT((M+1)*(M+2))

C      GET TERM WITH 2F1 HYPERGEOMETRIC  GAUSS FUNCTION
C       CALL HYGFX(DFLOAT(-N+1),DFLOAT(M+1),DFLOAT(N+M+1),XB/XU,HF)
C      DIRECT EXPRESSION FOR 2F1 WITH N=2 case
       HF=1.D0-(XB/XU)*(DFLOAT(M+1)/DFLOAT(M+3))
       RKHS_DKN2=NSQR*XUPOWM*BETATERM*HF
       RETURN
       END

       DOUBLE PRECISION FUNCTION RKHS_DK(X,Y,N,M)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IF (X .LT. Y) THEN
           XB=X
           XU=Y
       ELSE
           XB=Y
           XU=X
       ENDIF

C       M=ABS(M)! FORCE M>=0
       NSQR=N*N
       XUPOWM=XU**(-(M+1))

C      TERM WITH BETA FUNCTION
       CALL BETA(DFLOAT(M+1),DFLOAT(N),BETATERM)

C      GET TERM WITH 2F1 HYPERGEOMETRIC  GAUSS FUNCTION
       CALL HYGFX(DFLOAT(-N+1),DFLOAT(M+1),DFLOAT(N+M+1),XB/XU,HF)
       RKHS_DK=NSQR*XUPOWM*BETATERM*HF
       RETURN
       END

C-------------------------------------------------------------
       DOUBLE PRECISION FUNCTION RKHS_EXPR(X,Y,XE,N)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER (ZERO=0.D0,ONE=1.D0)
       
       SUMA=ZERO
       Z1=(X-XE)/XE
       Z2=(Y-XE)/XE
       DO I=0,N
       SUMA=SUMA+(Z1**I)*(Z2**I)
       ENDDO

       RKHS_EXPR=SUMA
       RETURN
       END


C-------------------------------------------------------------
       DOUBLE PRECISION FUNCTION RKHS_AK(X,Y,N)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER (ZERO=0.D0,ONE=1.D0)
       
       IF ((X.EQ.ZERO).OR.(Y.EQ.ZERO)) THEN
          RKHS_AK=ONE
          RETURN
       ENDIF
       IF (X .LT. Y) THEN
           XB=X
           XU=Y
       ELSE
           XB=Y
           XU=X
       ENDIF

       SUMA=ZERO
C      GET POLYNOMIAL TERM
       DO I=0,N-1
       SUMA=SUMA+(XU**I)*(XB**I)
       END DO

C      GET TERM WITH 2F1 HYPERGEOMETRIC  GAUSS FUNCTION
       CALL HYGFX(DFLOAT(1),DFLOAT(-N+1),DFLOAT(N+1),XB/XU,HF)
       
       RKHS_AK=SUMA+N*(XB**N)*(XU**(N-1))*HF
       RETURN
       END

       DOUBLE PRECISION FUNCTION RKHS_AKN2(X,Y)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER (ZERO=0.D0,ONE=1.D0)
       
       IF ((X.EQ.ZERO).OR.(Y.EQ.ZERO)) THEN
          RKHS_AKN2=ONE
          RETURN
       ENDIF 
       IF (X .LT. Y) THEN
           XB=X
           XU=Y
       ELSE
           XB=Y
           XU=X
       ENDIF

C      GET POLYNOMIAL TERM
       SUMA=1.D0+XU*XB

C      GET TERM WITH 2F1 HYPERGEOMETRIC  GAUSS FUNCTION
C       CALL HYGFX(DFLOAT(1),DFLOAT(-N+1),DFLOAT(N+1),XB/XU,HF)
       HF=1.D0-(XB/XU)/3.D0
       
       RKHS_AKN2=SUMA+2.D0*(XB**2)*XU*HF
       RETURN
       END


       DOUBLE PRECISION FUNCTION RKHS_OP(X,Y,N,M)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER (ZERO=0.D0,ONE=1.D0)

       SUMA=ZERO
C      GET ORTHOGONAL POLYNOMIAL TERM
       DO I=ABS(M),N
       SUMA=SUMA+PLGNDR(I,M,X)*PLGNDR(I,M,Y)
       END DO
       RKHS_OP=SUMA
       RETURN
       END

       DOUBLE PRECISION FUNCTION RKHS_OP2L(X,Y,N,M)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER (ZERO=0.D0,ONE=1.D0)
 
       SUMA=ZERO
C      GET ORTHOGONAL POLYNOMIAL TERM
       DO I=ABS(M),N,2
       SUMA=SUMA+PLGNDR(I,M,X)*PLGNDR(I,M,Y)
       END DO
       RKHS_OP2L=SUMA
       RETURN
       END


            



        SUBROUTINE BETA(P,Q,BT)
C
C       ==========================================
C       Purpose: Compute the beta function B(p,q)
C       Input :  p  --- Parameter  ( p > 0 )
C                q  --- Parameter  ( q > 0 )
C       Output:  BT --- B(p,q)
C       Routine called: GAMMA for computing â(x)
C       ==========================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        CALL GAMMA(P,GP)
        CALL GAMMA(Q,GQ)
        PPQ=P+Q
        CALL GAMMA(PPQ,GPQ)
        BT=GP*GQ/GPQ
        RETURN
        END


        SUBROUTINE GAMMA(X,GA)
C
C       ==================================================
C       Purpose: Compute gamma function (x)
C       Input :  x  --- Argument of ( x )
C                       ( x is not equal to 0,-1,-2,...)
C       Output:  GA --- (x)
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,
     &          -0.6558780715202538D0, -0.420026350340952D-1,
     &          0.1665386113822915D0,-.421977345555443D-1,
     &          -.96219715278770D-2, .72189432466630D-2,
     &          -.11651675918591D-2, -.2152416741149D-3,
     &          .1280502823882D-3, -.201348547807D-4,
     &          -.12504934821D-5, .11330272320D-5,
     &          -.2056338417D-6, .61160950D-8,
     &          .50020075D-8, -.11812746D-8,
     &          .1043427D-9, .77823D-11,
     &          -.36968D-11, .51D-12,
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END



        SUBROUTINE HYGFX(A,B,C,X,HF)
C
C       ====================================================
C       Purpose: Compute hypergeometric function F(a,b,c,x)
C       Input :  a --- Parameter
C                b --- Parameter
C                c --- Parameter, c <> 0,-1,-2,...
C                x --- Argument   ( x < 1 )
C       Output:  HF --- F(a,b,c,x)
C       Routines called:
C            (1) GAMMA for computing gamma function
C            (2) PSI_AAA for computing psi function
C       ====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        LOGICAL L0,L1,L2,L3,L4,L5
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        L0=C.EQ.INT(C).AND.C.LT.0.0
        L1=1.0D0-X.LT.1.0D-15.AND.C-A-B.LE.0.0
        L2=A.EQ.INT(A).AND.A.LT.0.0
        L3=B.EQ.INT(B).AND.B.LT.0.0
        L4=C-A.EQ.INT(C-A).AND.C-A.LE.0.0
        L5=C-B.EQ.INT(C-B).AND.C-B.LE.0.0
        IF (L0.OR.L1) THEN
           WRITE(*,*)'The hypergeometric series is divergent'
           RETURN
        ENDIF
        EPS=1.0D-15
        IF (X.GT.0.95) EPS=1.0D-8
        IF (X.EQ.0.0.OR.A.EQ.0.0.OR.B.EQ.0.0) THEN
           HF=1.0D0
           RETURN
        ELSE IF (1.0D0-X.EQ.EPS.AND.C-A-B.GT.0.0) THEN
           CALL GAMMA(C,GC)
           CALL GAMMA(C-A-B,GCAB)
           CALL GAMMA(C-A,GCA)
           CALL GAMMA(C-B,GCB)
           HF=GC*GCAB/(GCA*GCB)
           RETURN
        ELSE IF (1.0D0+X.LE.EPS.AND.DABS(C-A+B-1.0).LE.EPS) THEN
           G0=DSQRT(PI)*2.0D0**(-A)
           CALL GAMMA(C,G1)
           CALL GAMMA(1.0D0+A/2.0-B,G2)
           CALL GAMMA(0.5D0+0.5*A,G3)
           HF=G0*G1/(G2*G3)
           RETURN
        ELSE IF (L2.OR.L3) THEN
           IF (L2) NM=INT(ABS(A))
           IF (L3) NM=INT(ABS(B))
           HF=1.0D0
           R=1.0D0
           DO 10 K=1,NM
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X
10            HF=HF+R
           RETURN
        ELSE IF (L4.OR.L5) THEN
           IF (L4) NM=INT(ABS(C-A))
           IF (L5) NM=INT(ABS(C-B))
           HF=1.0D0
           R=1.0D0
           DO 15 K=1,NM
              R=R*(C-A+K-1.0D0)*(C-B+K-1.0D0)/(K*(C+K-1.0D0))*X
15            HF=HF+R
           HF=(1.0D0-X)**(C-A-B)*HF
           RETURN
        ENDIF
        AA=A
        BB=B
        X1=X
        IF (X.LT.0.0D0) THEN
           X=X/(X-1.0D0)
           IF (C.GT.A.AND.B.LT.A.AND.B.GT.0.0) THEN
              A=BB
              B=AA
           ENDIF
           B=C-B
        ENDIF
        IF (X.GE.0.75D0) THEN
           GM=0.0D0
           IF (DABS(C-A-B-INT(C-A-B)).LT.1.0D-15) THEN
              M=INT(C-A-B)
              CALL GAMMA(A,GA)
              CALL GAMMA(B,GB)
              CALL GAMMA(C,GC)
              CALL GAMMA(A+M,GAM)
              CALL GAMMA(B+M,GBM)
              CALL PSI_AAA(A,PA)
              CALL PSI_AAA(B,PB)
              IF (M.NE.0) GM=1.0D0
              DO 30 J=1,ABS(M)-1
30               GM=GM*J
              RM=1.0D0
              DO 35 J=1,ABS(M)
35               RM=RM*J
              F0=1.0D0
              R0=1.0D0
              R1=1.0D0
              SP0=0.D0
              SP=0.0D0
              IF (M.GE.0) THEN
                 C0=GM*GC/(GAM*GBM)
                 C1=-GC*(X-1.0D0)**M/(GA*GB*RM)
                 DO 40 K=1,M-1
                    R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(K-M))*(1.0-X)
40                  F0=F0+R0
                 DO 45 K=1,M
45                  SP0=SP0+1.0D0/(A+K-1.0)+1.0/(B+K-1.0)-1.0/K
                 F1=PA+PB+SP0+2.0D0*EL+DLOG(1.0D0-X)
                 DO 55 K=1,250
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0))
                    SM=0.0D0
                    DO 50 J=1,M
50                     SM=SM+(1.0D0-A)/((J+K)*(A+J+K-1.0))+1.0/
     &                    (B+J+K-1.0)
                    RP=PA+PB+2.0D0*EL+SP+SM+DLOG(1.0D0-X)
                    R1=R1*(A+M+K-1.0D0)*(B+M+K-1.0)/(K*(M+K))*(1.0-X)
                    F1=F1+R1*RP
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 60
55                  HW=F1
60               HF=F0*C0+F1*C1
              ELSE IF (M.LT.0) THEN
                 M=-M
                 C0=GM*GC/(GA*GB*(1.0D0-X)**M)
                 C1=-(-1)**M*GC/(GAM*GBM*RM)
                 DO 65 K=1,M-1
                    R0=R0*(A-M+K-1.0D0)*(B-M+K-1.0)/(K*(K-M))*(1.0-X)
65                  F0=F0+R0
                 DO 70 K=1,M
70                  SP0=SP0+1.0D0/K
                 F1=PA+PB-SP0+2.0D0*EL+DLOG(1.0D0-X)
                 DO 80 K=1,250
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0))
                    SM=0.0D0
                    DO 75 J=1,M
75                     SM=SM+1.0D0/(J+K)
                    RP=PA+PB+2.0D0*EL+SP-SM+DLOG(1.0D0-X)
                    R1=R1*(A+K-1.0D0)*(B+K-1.0)/(K*(M+K))*(1.0-X)
                    F1=F1+R1*RP
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 85
80                  HW=F1
85               HF=F0*C0+F1*C1
              ENDIF
           ELSE
              CALL GAMMA(A,GA)
              CALL GAMMA(B,GB)
              CALL GAMMA(C,GC)
              CALL GAMMA(C-A,GCA)
              CALL GAMMA(C-B,GCB)
              CALL GAMMA(C-A-B,GCAB)
              CALL GAMMA(A+B-C,GABC)
              C0=GC*GCAB/(GCA*GCB)
              C1=GC*GABC/(GA*GB)*(1.0D0-X)**(C-A-B)
              HF=0.0D0
              R0=C0
              R1=C1
              DO 90 K=1,250
                 R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(A+B-C+K))*(1.0-X)
                 R1=R1*(C-A+K-1.0D0)*(C-B+K-1.0)/(K*(C-A-B+K))
     &              *(1.0-X)
                 HF=HF+R0+R1
                 IF (DABS(HF-HW).LT.DABS(HF)*EPS) GO TO 95
90               HW=HF
95            HF=HF+C0+C1
           ENDIF
        ELSE
           A0=1.0D0
           IF (C.GT.A.AND.C.LT.2.0D0*A.AND.
     &         C.GT.B.AND.C.LT.2.0D0*B) THEN
              A0=(1.0D0-X)**(C-A-B)
              A=C-A
              B=C-B
           ENDIF
           HF=1.0D0
           R=1.0D0
           DO 100 K=1,250
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X
              HF=HF+R
              IF (DABS(HF-HW).LE.DABS(HF)*EPS) GO TO 105
100           HW=HF
105        HF=A0*HF
        ENDIF
        IF (X1.LT.0.0D0) THEN
           X=X1
           C0=1.0D0/(1.0D0-X)**AA
           HF=C0*HF
        ENDIF
        A=AA
        B=BB
        IF (K.GT.120) WRITE(*,115)
115     FORMAT(1X,'Warning! You should check the accuracy')
        RETURN
        END



        SUBROUTINE PSI_AAA(X,PS)
C
C       ======================================
C       Purpose: Compute Psi function
C       Input :  x  --- Argument of psi(x)
C       Output:  PS --- psi(x)
C       ======================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XA=DABS(X)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        S=0.0D0
        IF (X.EQ.INT(X).AND.X.LE.0.0) THEN
           PS=1.0D+300
           RETURN
        ELSE IF (XA.EQ.INT(XA)) THEN
           N=XA
           DO 10 K=1 ,N-1
10            S=S+1.0D0/K
           PS=-EL+S
        ELSE IF (XA+.5.EQ.INT(XA+.5)) THEN
           N=XA-.5
           DO 20 K=1,N
20            S=S+1.0/(2.0D0*K-1.0D0)
           PS=-EL+2.0D0*S-1.386294361119891D0
        ELSE
           IF (XA.LT.10.0) THEN
              N=10-INT(XA)
              DO 30 K=0,N-1
30               S=S+1.0D0/(XA+K)
              XA=XA+N
           ENDIF
           X2=1.0D0/(XA*XA)
           A1=-.8333333333333D-01
           A2=.83333333333333333D-02
           A3=-.39682539682539683D-02
           A4=.41666666666666667D-02
           A5=-.75757575757575758D-02
           A6=.21092796092796093D-01
           A7=-.83333333333333333D-01
           A8=.4432598039215686D0
           PS=DLOG(XA)-.5D0/XA+X2*(((((((A8*X2+A7)*X2+
     &        A6)*X2+A5)*X2+A4)*X2+A3)*X2+A2)*X2+A1)
           PS=PS-S
        ENDIF
        IF (X.LT.0.0) PS=PS-PI*DCOS(PI*X)/DSIN(PI*X)-1.0D0/X
        RETURN
        END
