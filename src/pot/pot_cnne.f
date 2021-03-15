*  routines to calculate the potentials for CN-Ne interaction
*  of sig/pi + atom scattering
*  This use ab initio points directly and spline fit them.
*  author: Moonbong Yang 27-sep-1995
*  revision: Moonbong Yang 23-oct-1996
*  latest revision:  mha 4-apr-1997
*  v1 pot is scaled by 0.17921 of vibrational matrix element of
*  7/3  then
*  divided by sqrt(2) since V1(R,q)=sq(2) * Sum(V1_l(R) * d_10(q))
*
* reference: 
* Moonbong Yang and M.  H. Alexander J. Chem. Phys. 107, 7148 (1997)

      include "common/syusr"
      include "common/bausr"
      include "common/ground"
*--------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
*--------------------------------------------------------------------

      character*(*) filnam
      logical ifirst
      common /coselb/ ibasty
      include "common/parbas"
      include "common/parpot"
      data ifirst/.true./

*     first call cnne_setup once to evaluate the necessary
*     long and short-range coefficients and spline polynomials

      if(ifirst) then
         call cnne_setup
         ifirst=.false.
      endif

      potnam='Yang CNNe (A,X) mrci+q avtz'

* Sig (iterm=1)
      iterm=1
      lammin(iterm)=0
      lammax(iterm)=6
      mproj(iterm)=0
      ntv(iterm)=1
      ivrow(ntv(iterm),iterm)=7
      ivcol(ntv(iterm),iterm)=7
* Pi(sum)
      iterm=2
      lammin(iterm)=0
      lammax(iterm)=6
      mproj(iterm)=0
      ntv(iterm)=1
      ivrow(ntv(iterm),iterm)=3
      ivcol(ntv(iterm),iterm)=3
* V1(Sig/Pi mixing term)
      iterm=3
      lammin(iterm)=1
      lammax(iterm)=5
      mproj(iterm)=1
      ntv(iterm)=1
      ivrow(ntv(iterm),iterm)=7
      ivcol(ntv(iterm),iterm)=3
* V2(Pi_dif)
      iterm=4
      lammin(iterm)=2
      lammax(iterm)=6
      mproj(iterm)=2
      ntv(iterm)=1
      ivrow(ntv(iterm),iterm)=3
      ivcol(ntv(iterm),iterm)=3
      if (ibasty .eq. 3) then
* Pi(sum)
         iterm=1
         lammin(iterm)=0
         lammax(iterm)=6
         mproj(iterm)=0
         ntv(iterm)=1
         ivrow(ntv(iterm),iterm)=3
         ivcol(ntv(iterm),iterm)=3
* V2(Pi_dif)
         iterm=2
         lammin(iterm)=2
         lammax(iterm)=6
         mproj(iterm)=2
         ntv(iterm)=1
         ivrow(ntv(iterm),iterm)=3
         ivcol(ntv(iterm),iterm)=3
      endif
      return
      end

*23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

* ----------------------------------------------------------------------
      subroutine cnne_setup
* ----------------------------------------------------------------------

*  this routine is to interpolate ab initio data
*  using piecewise cubic spline routine
*  For the range shorter than ab initio points,
*  Use single expontial for extrapolate using
*  two of the innermost r points.
*  For the range longer than ab initio points,
*  Use leading long range r^(-n) term to extrapolate,
*  while c_coeffincent were obtained from the
*  outermost ab initio points.
*
*     calculate the short range potential
*     V = A * exp (-alpha * r)
*     evaluate A and alpha from the innermost two ab initio points
*                  1           v1
*     alpha =  ---------  log -----
*               r2 - r1        v2
*
*       A   =  v1 * exp (alpha * r1)
*     using alpha and A, evaluate the potential shorter than
*     calculated ab initio points
*
*     calculate long range potential from the outermost ab initio value
*     for dispersion and induction, the leading long range term may vary
*     if r^(-n) term is leading
*     coef = - v1 * r1 ^ n
*     then potential energy can be evaluated as
*             coef
*     v = - --------
*            r1 ^ n
*
*     For both dispersion (c6) and induction(c8) exist,
*     use two of the last points
*     c6 = (r1^8 * v1 - r2^8 * v2)/(r2^2-r1^2)
*     c8 = -0.5 [c6*(r1^2 + r2^2) + r1^8*v1 + r2^8*v2]

* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      parameter (nang=7,na12=5,nsp=7,n12=5)
      parameter (nr=13,npsp=nr*nsp,np12=nr*n12)

      dimension x(nr),vs(npsp),vp(npsp),v1(np12),v2(np12),
     :          aas(nang),aap(nang),aa1(na12),aa2(na12),
     :          alps(nang),alpp(nang),alp1(na12),alp2(na12),
     :          cfs6(nang),cfp6(nang),cfs8(nang),cfp8(nang),
     :          cf16(na12),cf26(na12),cf18(na12),cf28(na12),
     :          ccs(nang,4,nr),ccp(nang,4,nr),cc1(na12,4,nr),
     :          cc2(na12,4,nr),vv1(nr),c1(4,nr),vv2(nr),c2(4,nr)

      common /setup/ aas,aap,aa1,aa2,alps,alpp,alp1,alp2,
     :               cfs6,cfs8,cfp6,cfp8,ccs,ccp,cc1,cc2

      data half, zero, one, two /0.5d0, 0.d0, 1.d0, 2.d0/

c     x where ab initio data available
      data x/   4.0,     5.0,     5.5,     6.0,     6.5,    7.0,
     : 7.5,     8.0,     8.5,     9.0,    10.0,    11.0,   12.0/


**    ab initio data points Here
**    Vsig, Vpi for angles 0 30 60 90 120 150 180
**    V1 and V2 for angles  30 60 90 120 150

      data vs/
     : 29734.39, 4178.69, 1356.34, 380.72,  69.73, -16.05, -31.34,
     :   -27.67,  -20.91,  -15.08,  -7.70,  -4.10,  -2.33,
     : 22263.00, 2793.49,  876.67, 226.79,  25.57, -25.57, -31.22,
     :   -25.63,  -18.98,  -13.61,  -6.97,  -3.73,  -2.15,
     :  9301.08, 1087.96,  306.48,  51.98, -19.07, -31.31, -27.45,
     :   -20.89,  -15.18,  -10.88,  -5.67,  -3.09,  -1.77,
     :  5728.99,  653.81,  171.05,  15.49, -25.56, -30.11, -25.04,
     :   -18.81,  -13.63,   -9.77,  -5.13,  -2.86,  -1.65,
     :  9690.05, 1304.22,  404.23,  89.39,  -8.14, -30.30, -29.17,
     :   -22.80,  -16.64,  -11.91,  -6.12,  -3.31,  -1.91,
     : 26382.11, 3702.98, 1309.36, 405.60,  88.86,  -8.75, -30.29,
     :   -28.70,  -22.16,  -16.04,  -8.14,  -4.26,  -2.38,
     : 34917.18, 5792.35, 2148.33, 716.06, 194.02,  22.10, -24.26,
     :   -29.75,  -24.49,  -18.10,  -9.26,  -4.82,  -2.65/

      data vp/
     : 41831.63, 5396.26, 1799.09, 547.60, 133.94,   9.16, -21.08,
     :   -23.23,  -18.78,  -13.90,  -7.17,  -3.79,  -2.13,
     : 23770.45, 3241.52, 1065.04, 304.34,  57.58, -12.02, -25.19,
     :   -22.68,  -17.36,  -12.62,  -6.48,  -3.46,  -1.97,
     :  8183.45, 1004.57,  287.68,  49.96, -17.41, -29.26, -25.71,
     :   -19.54,  -14.15,  -10.12,  -5.26,  -2.86,  -1.64,
     :  4402.45,  492.68,  116.84,  -1.70, -30.50, -30.78, -24.72,
     :   -18.25,  -13.12,   -9.37,  -4.92,  -2.72,  -1.59,
     :  8651.19, 1379.38,  468.63, 125.47,   9.30, -22.54, -26.04,
     :   -21.68,  -16.32,  -11.86,  -6.16,  -3.34,  -1.91,
     : 28895.40, 5692.97, 2184.27, 772.20, 238.51,  50.98,  -7.36,
     :   -20.63,  -19.83,  -15.71,  -8.46,  -4.49,  -2.48,
     : 58224.86, 10333.74, 3953.63, 1429.59, 474.52, 131.39, 17.34,
     :   -14.91,  -20.00,  -17.28,  -9.72,  -5.11,  -2.78/

      data v1/
     : 4517.88,  589.18, 218.14,  81.04, 30.02, 11.06, 3.89, 1.18,
     :    0.24,    0.05,   0.11,   0.06,  0.03,
     : 1417.73,  283.16, 116.27,  46.58, 18.15,  6.74, 2.31, 0.67,
     :    0.11,    0.06,   0.11,   0.08,  0.05,
     :  780.19,  197.14,  88.21,  37.97, 15.74,  6.21, 2.29, 0.76,
     :    0.20,    0.01,   0.03,   0.02,  0.01,
     : 2457.08,  557.51, 239.54, 100.17, 41.07, 16.50, 6.41, 2.31,
     :    0.68,    0.08,   0.11,   0.05,  0.01,
     : 8218.96, 1123.24, 422.65, 162.72, 63.21, 24.36, 9.26, 3.49,
     :    1.24,    0.35,   0.07,   0.07,  0.03/

      data v2/
     : 2217.05,  193.36,  64.89,  22.40,  7.73,  2.60, 0.81, 0.22,
     :    0.03,    0.0,    0.0,    0.0,   0.0,
     : 1862.09,  257.71,  95.28,  34.70, 12.26,  4.10, 1.18, 0.21,
     :   -0.05,   -0.11,  -0.08,  -0.04, -0.01,
     : 1466.67,  231.68,  89.35,  33.38, 12.14,  3.88, 1.17, 0.19,
     :   -0.09,   -0.13,  -0.09,  -0.04, -0.01,
     : 1937.99,  269.41,  98.96,  36.12, 12.96,  4.45, 1.39, 0.32,
     :   -0.01,   -0.07,  -0.07,  -0.04, -0.02,
     : 3724.55,  213.70,  65.23,  20.85,  6.78,  2.18, 0.65, 0.17,
     :    0.02,    0.0,    0.0,    0.0,   0.0/


c------------------------------------------------------------------
      ibcbeg=0
      ibcend=0

** for Vsig
      do i = 1, nang
         do j = 1, nr
            i1=(i-1)*nr
            k = i1 + j
            vv1(j)=vs(k)
            c1(1,j)=vs(k)
            vv2(j)=vp(k)
            c2(1,j)=vp(k)
         enddo
         call extpol (nr,x,vv1,aas(i),alps(i),cfs6(i),cfs8(i),1)
         call extpol (nr,x,vv2,aap(i),alpp(i),cfp6(i),cfp8(i),1)
         call cubspl (x,c1,nr,ibcbeg,ibcend)
         call cubspl (x,c2,nr,ibcbeg,ibcend)
         do l = 1, 4
           do j = 1, nr
             ccs(i,l,j) = c1(l,j)
             ccp(i,l,j) = c2(l,j)
           enddo
         enddo
      enddo

** for V1 and V2
      do i = 1, na12
         do j = 1, nr
            i1=(i-1)*nr
            k = i1 + j
            vv1(j)=v1(k)
            c1(1,j)=v1(k)
            vv2(j)=v2(k)
            c2(1,j)=v2(k)
         enddo
         call extpol (nr,x,vv1,aa1(i),alp1(i),cf16(i),cf18(i),0)
         call extpol (nr,x,vv2,aa2(i),alp2(i),cf26(i),cf28(i),0)
         call cubspl (x,c1,nr,ibcbeg,ibcend)
         call cubspl (x,c2,nr,ibcbeg,ibcend)
         do l = 1, 4
           do j = 1, nr
             cc1(i,l,j) = c1(l,j)
             cc2(i,l,j) = c2(l,j)
           enddo
         enddo
      enddo

      return
      end


*23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

* ----------------------------------------------------------------------
      subroutine pot (vv0, r)
* ----------------------------------------------------------------------

*     subroutine to calculate the r-dependent coefficients for
*     CN-Ne potentials from ab initio MRCI+Q potentials
*     in atomic units (distance and energy)

*     on entry:
*       r:  interparticle distance
*     on return:
*       vv0 contains isotropic term l=0 in V_sig (iterm=1)
*       vvl variable in common block /covvl/
*           contains all oter terms
*           iterm=1 : l=1 to 6 vvl(1:6)
*           iterm=2 : l=0 to 6 vvl(7:13)
*           iterm=3 : l=1 to 5 vvl(14:18)
*           iterm=4 : l=2 to 6 vvl(19:23)

*      uses linear least squares routines from cmlib

*      author: moonbong yang
*      latest revision date:  4-apr-1997
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      parameter (nang=7,na12=5,nsp=7,n12=5)
      parameter (nr=13,npsp=nr*nsp,np12=nr*n12)
      dimension
     :  x(nr),d0(nang*nsp),d1(na12*n12),d2(na12*n12),
     :  vs(nang),vp(nang),v1(na12),v2(na12),
     :  vls(nsp),vlp(nsp),vl1(n12),vl2(n12),
     :  aaa(nang*nsp),kpvt(nsp),qraux(nsp),work(nsp),rsd(nsp)

      dimension
     :          aas(nang),aap(nang),aa1(na12),aa2(na12),
     :          alps(nang),alpp(nang),alp1(na12),alp2(na12),
     :          cfs6(nang),cfp6(nang),cfs8(nang),cfp8(nang),
     :          cf16(na12),cf26(na12),cf18(na12),cf28(na12),
     :          ccs(nang,4,nr),ccp(nang,4,nr),cc1(na12,4,nr),
     :          cc2(na12,4,nr),c1(4,nr),c2(4,nr)

      common /covvl/ vvl(2*nsp+2*n12-1)
      common /coselb/ ibasty
      common /setup/ aas,aap,aa1,aa2,alps,alpp,alp1,alp2,
     :               cfs6,cfs8,cfp6,cfp8,ccs,ccp,cc1,cc2
* common scale is transmitted from hisystem_cn/sysgpi
*      common /scale/ scalev1
* use scale in this pot only
      data scalev1 /0.17921d0/

      data half, zero, one, two /0.5d0, 0.d0, 1.d0, 2.d0/

c     x where ab initio data available
      data x/   4.0,     5.0,     5.5,     6.0,     6.5,    7.0,
     : 7.5,     8.0,     8.5,     9.0,    10.0,    11.0,   12.0/


****---- check for very small numbers ----****
* coefficicients for d0 rotation matrices
* stored (by column) for each of 7 angles and for l=0:6(row)
* number of l expansion is same as number of angles
* angles are 0 30 60 90 120 150 180

      data d0/
     : 1d0,  1d0,  1d0,  1d0,  1d0,  1d0,  1d0,
     : 1d0,  8.6602541d-1,  5.d-1,         0d0,     -5.d-1,
     :      -8.6602539d-1, -1d0,
     : 1d0,  6.25d-1,      -1.25d-1,      -5.d-1,   -1.25d-1,
     :       6.25d-1,       1d0,
     : 1d0,  3.2475954d-1, -4.375d-1,      0d0,      4.3751d-1,
     :      -3.2475948d-1, -1d0,
     : 1d0,  2.3437511d-2, -2.8906251d-1,  3.75d-1, -2.8906248d-1,
     :       2.3437446d-2,  1d0,
     : 1d0, -2.2327216d-1,  8.9843733d-2,  0d0,     -8.9843784d-2,
     :       2.2327222d-1, -1d0,
     : 1d0, -3.7402343d-1,  3.2324218d-1, -3.125d-1, 3.2324220d-1,
     :      -3.7402346d-1,  1d0/

* coefficicients for d1 rotation matrices
* stored (by column) for each of 5 angles and for l=1:5
* angles are 30 60 90 120 150
      data d1/
     :-3.5355339d-1, -6.1237243d-1, -7.0710678d-1, -6.1237243d-1,
     :-3.5355339d-1,
     :-5.3033008d-1, -5.3033008d-1,  0.d0,          5.3033008d-1,
     : 5.3033008d-1,
     :-5.9539246d-1, -9.3750000d-2,  4.3301270d-1, -9.3749999d-2,
     :-5.9539246d-1,
     :-5.4463828d-1,  3.0257682d-1,  0.d0,         -3.0257682d-1,
     : 5.4463828d-1,
     :-3.9581512d-1,  3.5205044d-1, -3.4232659d-1,  3.5205044d-1,
     :-3.9581512d-1/

* coefficicients for d2 rotation matrices
* stored (by column) for each of 5 angles and for l=2:6
* angles are 30 60 90 120 150
      data d2/
     : 1.5309311d-1,  4.5927932d-1,  6.1237244d-1,  4.5927934d-1,
     : 1.5309312d-1,
     : 2.9646353d-1,  5.1348990d-1,  0.d0,         -5.1348989d-1,
     :-2.9646355d-1,
     : 4.1999000d-1,  2.2234766d-1, -3.9528471d-1,  2.2234762d-1,
     : 4.1999002d-1,
     : 4.9023048d-1, -1.6982081d-1,  0.d0,          1.6982085d-1,
     :-4.9023049d-1,
     : 4.8532921d-1, -3.4523418d-1,  3.2021721d-1, -3.4523418d-1,
     : 4.8532920d-1/


*    determine Vs, Vp, V1, V2 potentials at angles
      xmin=x(1)
      xmax=x(nr)
      rr=r*r
      rrr=rr*r
      r6=rrr*rrr
      r8=r6*rr

* for Vsig and Vpi
      do i = 1, nang
         if (r .lt. x(2)) then
           vs(i) = aas(i) * dexp (-alps(i) * r)
           vp(i) = aap(i) * dexp (-alpp(i) * r)
         else if (r .gt. x(nr-1)) then
           vs(i) = - cfs6(i) / r6  - cfs8(i) / r8
           vp(i) = - cfp6(i) / r6  - cfp8(i) / r8
         else
           do l = 1 , 4
            do j = 1, nr
             c1(l,j) = ccs(i,l,j)
             c2(l,j) = ccp(i,l,j)
            enddo
           enddo
          call splval (r,vs(i),x,c1,nr)
          call splval (r,vp(i),x,c2,nr)
         endif
       enddo

      sq2 = dsqrt(two)

* for V1 and V2
      do i = 1, na12
         if (r .lt. xmin) then
           v1(i) = aa1(i) * dexp (-alp1(i) * r) * scalev1 / sq2
           v2(i) = aa2(i) * dexp (-alp2(i) * r)
         else if (r .gt. xmax) then
           v1(i) = zero
           v2(i) = zero
         else
           do l = 1 , 4
            do j = 1, nr
             c1(l,j) = cc1(i,l,j)
             c2(l,j) = cc2(i,l,j)
            enddo
           enddo
          call splval (r,v1i,x,c1,nr)
          v1(i) = v1i * scalev1 / sq2
          call splval (r,v2(i),x,c2,nr)
         endif
       enddo

* convert to hartree
      conv=1.d0/219474.6
      tol=1.e-10

* solve simultaneous equations for solutions
* first for Vsigma
      call dcopy(nsp*nsp,d0,1,aaa,1)
      call dqrank(aaa,nsp,nsp,nsp,tol,kr,kpvt,qraux,work)
      call dqrlss(aaa,nsp,nsp,nsp,kr,vs,vls,rsd,kpvt,qraux)
      call dscal(nsp,conv,vls,1)
      vv0=0.d0
      call dcopy(nsp,vls(1),1,vvl,1)

* next for Vpi
      call dcopy(nsp*nsp,d0,1,aaa,1)
      call dqrank(aaa,nsp,nsp,nsp,tol,kr,kpvt,qraux,work)
      call dqrlss(aaa,nsp,nsp,nsp,kr,vp,vlp,rsd,kpvt,qraux)
      call dscal(nsp,conv,vlp,1)
      call dcopy(nsp,vlp,1,vvl(nsp+1),1)

* next for V1
      call dcopy(n12*n12,d1,1,aaa,1)
      call dqrank(aaa,n12,n12,n12,tol,kr,kpvt,qraux,work)
      call dqrlss(aaa,n12,n12,n12,kr,v1,vl1,rsd,kpvt,qraux)
      call dscal(n12,conv,vl1,1)
      call dcopy(n12,vl1,1,vvl(2*nsp+1),1)

* next for V2
      call dcopy(n12*n12,d2,1,aaa,1)
      call dqrank(aaa,n12,n12,n12,tol,kr,kpvt,qraux,work)
      call dqrlss(aaa,n12,n12,n12,kr,v2,vl2,rsd,kpvt,qraux)
      call dscal(n12,conv,vl2,1)
      call dcopy(n12,vl2,1,vvl(2*nsp+n12+1),1)

* if just pi state, reduce number of vvl matrix elements
      if (ibasty .eq. 3) then
        call dscal(2*nsp+n12+2*n12,0.d0,vvl,1)
        vv0=vlp(1)
        call dcopy(nsp-1,vlp(2),1,vvl,1)
        call dcopy(n12,vl2,1,vvl(nsp),1)
      endif

      return
      end


*23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

* ----------------------------------------------------------------------
      subroutine extpol(n,x,v,aa,alpha,coef6,coef8,il)
* ----------------------------------------------------------------------

*     obtain parameters for short and longe extrapolation
*     input :
*       n : number of data point
*       x(n) : r points of ab initio data
*       v(n) : ab initio potential at aan angle
*       il : if 1, then include long range
*            if 0, then long range pot is zero
*     output :
*       aa    : pre-expontial coeff
*       alpha : expontial coefficient
*       coef6 : long range coefficient

      implicit double precision (a-h,o-z)
      dimension x(n),v(n)
      data zero / 0.d0 /

      x1 = x(1)
      x2 = x(2)
      xn1 = x(n)
      xn2 = x(n-1)
      dx = x2 - x1

      xn12 = xn1 * xn1
      xn14 = xn12 * xn12
      xn16 = xn12 * xn14
      xn18 = xn14 * xn14

      xn22 = xn2 * xn2
      xn24 = xn22 * xn22
      xn26 = xn22 * xn24
      xn28 = xn24 * xn24


      v1 = v(1)
      v2 = v(2)
      vn1 = v(n)
      vn2 = v(n-1)
      alpha = dlog(v1/v2) / dx
      aa = v1 * dexp(alpha * x1)

      if (il.eq.0) then
         coef6=zero
         coef8=zero
      else

*     c6 = (r1^8 * v1 - r2^8 * v2)/(r2^2-r1^2)
*     c8 = -0.5 [c6*(r1^2 + r2^2) + r1^8*v1 + r2^8*v2]
*        = -c6*r1^2 - r1^8*v1


         coef6 = (xn18*vn1 - xn28*vn2)/(xn22-xn12)
         coef8 = - coef6*xn12 - xn18*vn1
      endif

      return
      end


*23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

* ----------------------------------------------------------------------
      subroutine cubspl (tau,c,n,ibcbeg,ibcend)
* ----------------------------------------------------------------------

c The following two subroutines calculte
c the piecewise cubic spline interpolation.
c First one cubspl caslculate the polynomial coefficients
c Second polval evaluate the spline for each r from coefficients
c This is same as the spline in MATLAB
c Only nmax>3 and "not-a-knot" condition is used
c For reference see  page 57 in
c         "a practical guide to splines" by c. de boor, applied
c         mathematics series 27, springer, 1979.

      implicit double precision (a-h,o-z)
      dimension c(4,n),tau(n)

c     for ibcbeg and ibcend=0 not-a-knot condition

      if (ibcbeg .ne. 0 .or. ibcend .ne. 0) then
         print*,'use only for not-a-knot condition'
         stop
      endif

      l=n-1

      do m=2,n
         c(3,m)=tau(m)-tau(m-1)
         c(4,m)=(c(1,m)-c(1,m-1))/c(3,m)
      enddo
c12
      c(4,1)=c(3,3)
      c(3,1)=c(3,2)+c(3,3)
      c(2,1)=((c(3,2)+2.*c(3,1))*c(4,2)*c(3,3)+c(3,2)**2*c(4,3))/c(3,1)
c19
      do m=2,l
         g=-c(3,m+1)/c(4,m-1)
         c(2,m)=g*c(2,m-1)+3.*(c(3,m)*c(4,m+1)+c(3,m+1)*c(4,m))
         c(4,m)=g*c(3,m-1)+2.*(c(3,m)+c(3,m+1))
      enddo
c21
      g=c(3,n-1)+c(3,n)
      c(2,n)=((c(3,n)+2.*g)*c(4,n)*c(3,n-1)
     +          +c(3,n)**2*2*(c(1,n-1)-c(1,n-2))/c(3,n-1))/g
      g=-g/c(4,n-1)
      c(4,n)=c(3,n-1)

c29
      c(4,n)=g*c(3,n-1)+c(4,n)
      c(2,n)=(g*c(2,n-1)+c(2,n))/c(4,n)
c30
      do j=l,1,-1
         c(2,j)=(c(2,j)-c(3,j)*c(2,j+1))/c(4,j)
      enddo
c50
      do i=2,n
         dtau=c(3,i)
         divdf1=(c(1,i)-c(1,i-1))/dtau
         divdf3=c(2,i-1)+c(2,i)-2.*divdf1
         c(3,i-1)=2.*(divdf1-c(2,i-1)-divdf3)/dtau
         c(4,i-1)=(divdf3/dtau)*(6./dtau)
      enddo

      return
      end

*23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

* ----------------------------------------------------------------------
      subroutine splval(r,f,tau,c,n)
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension c(4,n),tau(n)

      do i = 1 , n-1
        if ((r .ge. tau(i)) .and. (r .le. tau(i+1))) then
          h = r-tau(i)
          f = c(1,i)+h*(c(2,i)+h*(c(3,i)+h*c(4,i)/3.)/2.)
        endif
      enddo
      return
      end
