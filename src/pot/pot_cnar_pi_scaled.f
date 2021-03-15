* System:  CN(A 2Pi)+Ar, from MRCI-PES
* Reference: A. Berning, Thesis, Universitaet Stuttgart, 1995.
*
* scaled potential, following J. Han, M. C. Heaven, and U. Schnupf
* JCP 128, 224309 (2008)
* p. dagdigian - jul 2010


      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(11)
      common /coconv/ econv
      include "common/parpot"
      dimension xxl(11)
      potnam='BERNING CN(A)-AR MRCI NEWFIT'
      econv=219474.6d0
1      print *, ' r (bohr)'
      rshift=0.5
      xfact=0.8
      read (5, *, end=99) r
      call pot(vv0,r)
*
* addition  by pjd (16-jun-2010)
* multiply vvl by econv to convert to cm^-1
      xx0 = vv0*econv
      do i=1,11
        xxl(i) = vvl(i)*econv
      end do
      write (6, 100) xx0, xxl
*
100   format(' vsum',/,7(1pe16.8),/,
     :    '  vdif',/,5e16.8)
      goto 1
99    r=4
      write (2,*) 'Vsum CN(A)-Ar'
      do i=1,40
         call pot(vv0,r)
         write (2,101) r,vv0,(vvl(j), j=1,6)
101      format(f8.4,7(pe16.8))
         r=r+0.2d0
      enddo
      r=4
      write (2,*) 'Vdif CN(A)-Ar'
      do i=1,40
         call pot(vv0,r)
         write (2,102) r,(vvl(j), j=7,11)
102      format(7(pe16.8))
         r=r+0.2d0
      enddo
      end
      include "common/syusr"
      include "common/bausr"
      include "common/ground"
* --------------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      potnam='BERNING CN(A)-AR MRCI NEWFIT'
      lammin(1)=1
      lammax(1)=6
      lammin(2)=2
      lammax(2)=6
      mproj(1)=0
      mproj(2)=2
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* ----------------------------------------------------------------------
      subroutine pot (vv0, r)
*  subroutine to calculate the r-dependent coefficients
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*  vv0        contains isotropic term (d00 term in vsum)
*  variable in common block /covvl/
*    vvl:     vector of length 6 to store r-dependence of each term
*             in potential expansion
*    vvl(1-6) expansion coefficients in dl0 (l=1:6) of vsum
*    vvl(7-11) expansion coefficients in dl2 (l=2:6) of vdif

* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  21-Aug-1998
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension xlam1(13),xlam2(13),r0(13),c1(13),c2(13),c3(13),
     :          clr(13),vsum(7),xsum(7),vdif(7),xdif(7),
     :          ddif(7),vap(7),va2p(7),
     :          d0(49),d2(25),aa(64)
      dimension kpvt(8),qraux(7),work(55),rsd(7),re(14)

      common /covvl/ vvl(11)
      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
* for distances beyond rmax difference potential is damped
      data rmax /13d0/
* coefficicients for d0 rotation matrices
* stored (by column) for each of 7 angles and for l=0:6
* angles are 0 30 60 90 120 150 180
      data d0/
     : 1d0,  1d0,  1d0,  1d0,  1d0,  1d0,  1d0,
     : 1d0,  8.6602541d-1,  5.0000001d-1,  1.3349125d-8, -4.9999998d-1,
     : -8.6602539d-1, -1d0,
     : 1d0,  6.2500001d-1, -1.2499999d-1, -5.0000000d-1, -1.2500002d-1,
     : 6.2499997d-1,  1d0,
     : 1d0,  3.2475954d-1, -4.3750000d-1, -2.0023687d-8,  4.3750001d-1,
     : -3.2475948d-1, -1d0,
     : 1d0,  2.3437511d-2, -2.8906251d-1,  3.7500000d-1, -2.8906248d-1,
     : 2.3437446d-2,  1d0,
     : 1d0, -2.2327216d-1,  8.9843733d-2,  2.5029609d-8, -8.9843784d-2,
     : 2.2327222d-1, -1d0,
     : 1d0, -3.7402343d-1,  3.2324218d-1, -3.1250000d-1,  3.2324220d-1,
     : -3.7402346d-1,  1d0/
* coefficicients for d2 rotation matrices
* stored (by column) for each of 5 angles and for l=2:6
* angles are 30 60 90 120 150
      data d2/
     : 1.5309311d-1,  4.5927932d-1,  6.1237244d-1,  4.5927934d-1,
     : 1.5309312d-1,
     : 2.9646353d-1,  5.1348990d-1,  1.8279042d-8, -5.1348989d-1,
     : -2.9646355d-1,
     : 4.1999000d-1,  2.2234766d-1, -3.9528471d-1,  2.2234762d-1,
     : 4.1999002d-1,
     : 4.9023048d-1, -1.6982081d-1, -2.4180900d-8,  1.6982085d-1,
     : -4.9023049d-1,
     : 4.8532921d-1, -3.4523418d-1,  3.2021721d-1, -3.4523418d-1,
     : 4.8532920d-1/

* coefficients for expansion of vap(1st 7 entries) and
* for va2p (entries 8:13)
      data xlam1/
     : 6.1447233d-1, 1.0609287d0, 5.4049025d-1, 6.1964594d-1,
     : 5.4153220d-1, 6.9260741d-1, 6.0317842d-1, 0d0, 6.6452456d-1,
     : 6.6271397d-1, 4.4375934d-1, 5.2046076d-1, 4.8299745d-1/

      data xlam2/
     : 2.3220336d0, 1.6486351d0, 1.6089740d0, 2.2603255d0, 1.5664448d0,
     : 2.1793710d0, 1.4905442d0, 0d0, 2.1933894d0, 2.1112711d0,
     : 1.3779602d0, 1.4255272d0, 1.3569892d0/

      data r0 /
     : 8.9448682d0, 6.1578261d0, 5.2608382d0, 8.0144937d0, 5.2718489d0,
     : 7.5318135d0, 5.9521662d0, 0d0, 8.3358418d0, 7.0524468d0,
     : 5.8726669d0, 6.1215945d0, 6.7123478d0/

      data c1 /
     : -1.4782126d+4, -4.4180810d+5, -4.5678934d+4, -1.0669847d+4,
     : -4.3052470d+4, -2.8747702d+4, -2.2278078d+4, 0d0, -2.1779772d+4,
     : -9.3813580d+3, -9.6992364d+2, -5.5701698d+3, -1.4054028d+3/

      data c2/
     : -5.8907918d+9, 4.6546445d+7, 1.5250938d+8, -1.5569607d+8,
     : 1.0499715d+8, -3.6262657d+9, 1.3984915d+8, 0d0, -1.1268121d+9,
     : 2.7175445d+8, 1.1669959d+7, 2.9558818d+7, 4.1716535d+7/

      data c3/
     : 1.5633281d+9, 1.0565527d+5, -2.9665089d+7, 5.5735619d+7,
     : -1.9890282d+7, 8.7614934d+8, -1.5617455d+7, 0d0, 3.7639462d+8,
     : -5.1404904d+4, -1.6341526d+6, -4.0322853d+6, -5.0150429d+6/

       data clr/
     : -5.7097820d+6, 1.5676693d+7, -1.9918062d+8, -3.2082715d+6,
     : -1.8191348d+8, 4.1180496d+6, -2.4927192d+7, 0d0, -3.8845716d+6,
     : 7.0885170d+6, 1.7350606d+6, -1.4746174d+7, 8.3559868d+6/

* determine A' and A" potentials at angles
      rm3=one/r**3
      rm6=rm3*rm3
      do 100 i=1,7
        vap(i)=c1(i)*exp(-xlam1(i)*r) +
     :        (c2(i)+c3(i)*r)*exp(-xlam2(i)*r)-
     :         half*(tanh(alph*(r-r0(i)))+1)*clr(i)*rm6
        j=i+7
* don't compute va2p and vdif for colinear geometries
        if (i.ne. 1. and. i .ne. 7 ) then
           va2p(i)=c1(j)*exp(-xlam1(j)*r) +
     :        (c2(j)+c3(j)*r)*exp(-xlam2(j)*r)-
     :         half*(tanh(alph*(r-r0(j)))+1)*clr(j)*rm6
           vsum(i)=half*(va2p(i)+vap(i))
           vdif(i-1)=half*(va2p(i)-vap(i))
* for long range damp out difference potential

           if (r .gt. rmax) then
             damp=-half*(tanh(3.d0*(r-rmax))-one)
             vdif(i-1)=vdif(i-1)*damp
           endif
        else
           vsum(i)=vap(i)
        endif
100   continue
* solve simultaneous equations for solutions
* first for vsigma
      tol=1.e-10
      call dcopy(49,d0,1,aa,1)
      call dqrank(aa,7,7,7,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,7,7,7,kr,vsum,xsum,rsd,kpvt,qraux)
* convert to hartree
      conv=1.d0/219474.6
      call dscal(7,conv,xsum,1)
*
* addition by p. dagdigian - jul 2010
* scale potential
      scfact = 1.76
      call dscal(7,scfact,xsum,1)
      vv0=xsum(1)
      call dcopy(6,xsum(2),1,vvl,1)
      call dcopy(25,d2,1,aa,1)
      call dqrank(aa,5,5,5,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,5,5,5,kr,vdif,xdif,rsd,kpvt,qraux)
      call dscal(5,conv,xdif,1)
* scale potential
      call dscal(5,scfact,xdif,1)
      call dcopy(5,xdif,1,vvl(7),1)
      end
