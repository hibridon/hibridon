* System:  NH(a 1Delta)+Ar, original ab initio MRCI PES's
* Reference: M. Yang et al. J. Chem. Phys. 103, 905 (1995)

      include "common/syusr"
      include "common/bausr"
      include "common/ground"
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(6)
1      print *, ' r (bohr)'
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vsum',/,6(1pe16.8),/,
     :    '  vdif',/,1e16.8)
      goto 1
99    end
      subroutine loapot(iunit,filnam)
*-----------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      potnam='YANG-ALEXANDER Ar-NH(a) mod MRCI'
      lammin(1)=1
      lammax(1)=5
      lammin(2)=4
      lammax(2)=4
      mproj(1)=0
      mproj(2)=4
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* ----------------------------------------------------------------------
      subroutine pot (vv0, rr)
*  subroutine to calculate the r-dependent coefficients in the
*  Ar-NH(a) potentials of jansen and hess
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*  vv0        contains isotropic term (d00 term in vsum)
*  variable in common block /covvl/
*    vvl:     vector of length 6 to store r-dependence of each term
*             in potential expansion
*    vvl(1-5) expansion coefficients in dl0 (l=1:5) of vsum
*    vvl(6) expansion coefficients in dl4 (l=4) of vdif

* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  20-sep-1993
* ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension xlam1(12),xlam2(12),r0(12),c1(12),c2(12),c3(12),
     :          clr(12),vsum(6),xsum(6),vdif(6),xdif(6),
     :          ddif(6),vap(6),va2p(6),
     :          d0(36),d4(4),aa(64)
      dimension kpvt(8),qraux(6),work(55),rsd(6), re(12)
      dimension depth(12),csa(6)

      common /covvl/ vvl(6)
      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/

** following cos of angles are for 0 45 67.5 90 135 180 degrees
** and needed to make A' and A" degenerate for linear geometry
** cos of angle
      data csa /
     : 1.0000,   0.7071,   0.3827,   0.0000,  -0.7071,  -1.0000/   
      data depth /
     : -79.6625, -75.7607, -75.2627, -70.3909, -64.4370, -67.1671,
     : -80.1856, -74.9602, -70.2489, -63.4581, -63.6724, -67.4099/ 
      data re /
     : 7.6568, 7.2626, 6.9822, 6.9669, 7.2185, 7.2508,
     : 7.6537, 7.2803, 7.0692, 7.0946, 7.2368, 7.2501/

* for distances beyond rmax difference potential is damped
      data rmax /12.5d0/
* coefficicients for d0 rotation matrices
* stored (by column) for each of 6 angles and for l=0:5
      data d0/
     : 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 7.0710679d-1, 3.8268344d-1,
     : 1.3349125d-8, -7.0710677d-1, -1, 1d0, 2.5000001d-1,
     : -2.8033008d-1, -5d-1, 2.4999997d-1, 1, 1d0,
     : -1.7677668d-1, -4.3391842d-1, -2.0023687d-8, 1.7677673d-1,
     : -1, 1d0, -4.0625000d-1, -8.0345887d-2, 3.7500000d-1,
     : -4.0625001d-1, 1, 1d0, -3.7565048d-1, 2.9179007d-1,
     : 2.5029609d-8, 3.7565046d-1, -1d0/
* coefficicients for d4 rotation matrices
* stored (by column) for each of 4 angles and for l=4
      data d4/
     : 1.3072813d-1, 3.8096968d-1, 5.2291252d-1, 1.3072814d-1/

* coefficients for expansion of vap(1st 6 entries) and
* for va2p (entries 7:12)
      data xlam1/
     : 6.7687817d-1, 6.7782801d-1, 7.0250000d-1, 7.0000000d-1,
     : 7.2432966d-1, 7.5897013d-1, 6.7737338d-1, 6.5263268d-1,
     : 6.1838165d-1, 6.3915468d-1, 6.4147363d-1, 7.1688272d-1/
      data xlam2/
     : 2.3162397, 2.2950485, 2.3541667, 2.3458333, 2.2753600, 2.2271122,
     : 2.3149360, 2.2896971, 2.3195567, 2.1828759, 2.2866119, 2.3210905/
      data r0 /
     : 8.2285238, 7.8451395, 8.5888889, 8.8111111, 1.0444211d+1,
     : 1.0048512d+1, 8.2293732, 7.9877038, 8.1420384, 7.9247038,
     : 7.9561804, 9.2837449/
      data c1 /
     : -2.2215425d+4, -1.6335887d+4, -1.5420109d+4, -1.4122589d+4,
     : -1.8703332d+4, -2.6374488d+4, -2.2402907d+4, -1.3457105d+4,
     : -8.1654233d+3, -9.0975250d+3, -1.0096571d+4, -1.8511021d+4/
      data c2/
     : -1.5142465d+9, -5.1879156d+8, -5.7780884d+8, -5.9418437d+8,
     : -3.9158528d+8, -2.0471668d+8, -1.4955193d+9, -4.3488891d+8,
     : -2.8883592d+8, -9.5212596d+7, -3.3482669d+8, -4.6368317d+8/
      data c3/
     : 4.7397810d+8, 1.6826932d+8, 1.5963645d+8, 1.5230062d+8,
     : 1.2172435d+8, 8.5546362d+7, 4.6912112d+8, 1.5143026d+8,
     : 9.9883598d+7, 3.7937376d+7, 1.1335367d+8, 1.6220293d+8/
       data clr/
     : -3.0307095d+6, -1.9003203d+6, 1.1295889d+5, 3.7027712d+5,
     : 3.9585314d+6, 4.3698954d+6, -3.0385391d+6, -2.7458469d+6,
     : -2.3721513d+6, -1.5594725d+6, -1.9296503d+6, 1.0122546d+6/

* determine A' and A" potentials at angles
* first shift radius out
      rrshift=0.5
      xfact=0.35
      rshift=0.5
      r=rr+rrshift
      rm3=one/r**3
      rm6=rm3*rm3
      do 100 i=1,6
        vap(i)=c1(i)*exp(-xlam1(i)*r) +
     :        (c2(i)+c3(i)*r)*exp(-xlam2(i)*r)-
     :         half*(tanh(alph*(r-r0(i)))+1)*clr(i)*rm6
        fact=xfact*(tanh(alph*(r-re(i)+rshift))+1)
        vap(i)=vap(i)+fact*vap(i)
        j=i+6
        va2p(i)=c1(j)*exp(-xlam1(j)*r) +
     :        (c2(j)+c3(j)*r)*exp(-xlam2(j)*r)-
     :         half*(tanh(alph*(r-r0(j)))+1)*clr(j)*rm6
        fact=xfact*(tanh(alph*(r-re(j)+rshift))+1)
        va2p(i)=va2p(i)+fact*va2p(i)
100   continue

* correct for differences in linear geometries
      dif0=vap(1)-va2p(1)
      dif180=vap(6)-va2p(6)
      ac=0.25*(dif0+dif180)
      bc=0.25*(dif0-dif180)
      do 110 i=1,6
        vap(i)=vap(i)-ac-bc*csa(i)
        va2p(i)=va2p(i)+ac+bc*csa(i)
        vsum(i)=half*(va2p(i)+vap(i))
* don't compute vdif for colinear geometries
        if (i.ne.1 .and. i.ne.6) then
          vdif(i-1)=half*(vap(i)-va2p(i))
* for long range damp out difference potential
          if (r .gt. rmax) then
            damp=-half*(tanh(3.d0*(r-rmax))-one)
            vdif(i-1)=vdif(i-1)*damp
          endif
        endif
110   continue


* solve simultaneous equations for solutions
* first for vsigma
      tol=1.e-10
      call dcopy(36,d0,1,aa,1)
      call dqrank(aa,6,6,6,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,6,6,6,kr,vsum,xsum,rsd,kpvt,qraux)
* convert to hartree
      conv=1.d0/219474.6
      call dscal(6,conv,xsum,1)
      vv0=xsum(1)
      call dcopy(5,xsum(2),1,vvl,1)
      call dcopy(4,d4,1,aa,1)
      call dqrank(aa,4,4,1,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,4,4,1,kr,vdif,xdif,rsd,kpvt,qraux)
      vvl(6)=conv*xdif(1)
      end
