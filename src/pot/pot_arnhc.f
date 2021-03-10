* System:  NH(c 1Pi)+Ar, original ab initio MRCI+Q PES's
* Reference: M. Yang et al., J. Chem. Phys. 102, 4069 (1995)

      include "common/syusr"
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(13)
1      print *, ' r (bohr)'
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vsum',/,8(1pe16.8),/,
     :    '  vdif',/,6e16.8)
      goto 1
99    end
      include "common/bausr"
      include "common/ground"
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      potnam='ALEXANDER-WERNER Ar-NH(c) MRCI'
      lammin(1)=1
      lammax(1)=7
      lammin(2)=2
      lammax(2)=7
      mproj(1)=0
      mproj(2)=2
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* ----------------------------------------------------------------------
      subroutine pot (vv0, r)
*  subroutine to calculate the r-dependent coefficients in the
*  Ar-NH(c) potentials of alexander and werner
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*  vv0        isotropic term in V0 (d00 term in vsum)
*  variable in common block /covvl/
*    vvl:     vector of length 13 to store r-dependence of each term
*             in potential expansion
*    vvl(1-7) expansion coefficients in dl0 (l=1:7) of vsum
*    vvl(8-13) expansion coefficients in dl2 (l=2:7) of vdif

* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  20-sep-1993
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension xlam1(16),xlam2(16),r0(16),c1(16),c2(16),c3(16),
     :          clr(16),vsum(8),xsum(8),vdif(8),xdif(8),
     :          ddif(8),vap(8),va2p(8),
     :          d0(64),d2(36),aa(64)
      dimension kpvt(8),qraux(8),work(55),rsd(8)

      common /covvl/ vvl(13)
      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
* for distances beyond rmax difference potential is damped
      data rmax /11.5d0/
* coefficicients for d0 rotation matrices
* stored (by column) for each of 8 angles and for l=0:7
      data d0/
     : 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0,
     : 1d0, 9.6592583d-1, 8.6602541d-1, 5.0000001d-1, 1.3349125d-8,
     : -4.9999998d-1, -8.6602539d-1, -1d0,
     : 1d0, 8.9951905d-1, 6.2500001d-1, -1.2499999d-1, -5.0000000d-1,
     : -1.2500002d-1, 6.2499997d-1, 1d0,
     : 1d0, 8.0416393d-1, 3.2475954d-1, -4.3750000d-1, -2.0023687d-8,
     : 4.3750001d-1, -3.2475948d-1, -1d0,
     : 1d0, 6.8469544d-1, 2.3437511d-2, -2.8906251d-1, 3.7500000d-1,
     : -2.8906248d-1, 2.3437446d-2, 1d0,
     : 1d0, 5.4712588d-1, -2.2327216d-1, 8.9843733d-2, 2.5029609d-8,
     : -8.9843784d-2, 2.2327222d-1, -1d0,
     : 1d0, 3.9830600d-1, -3.7402343d-1, 3.2324218d-1, -3.1250000d-1,
     : 3.2324220d-1, -3.7402346d-1, 1d0,
     : 1d0, 2.4554105d-1, -4.1017805d-1, 2.2314455d-1, -2.9201210d-8,
     : -2.2314450d-1, 4.1017804d-1, -1d0/
* coefficicients for d2 rotation matrices
* stored (by column) for each of 6 angles and for l=2:7
      data d2/
     : 4.1021174d-2, 1.5309311d-1, 4.5927932d-1, 6.1237244d-1,
     : 4.5927934d-1, 1.5309312d-1, 8.8600642d-2, 2.9646353d-1,
     : 5.1348990d-1, 1.8279042d-8, -5.1348989d-1, -2.9646355d-1,
     : 1.4645800d-1, 4.1999000d-1, 2.2234766d-1, -3.9528471d-1,
     : 2.2234762d-1, 4.1999002d-1, 2.1086100d-1, 4.9023048d-1,
     : -1.6982081d-1, -2.4180900d-8, 1.6982085d-1, -4.9023049d-1,
     : 2.7741249d-1, 4.8532921d-1, -3.4523418d-1, 3.2021721d-1,
     : -3.4523418d-1, 4.8532920d-1, 3.4147272d-1, 4.0112588d-1,
     : -1.9131360d-1, 2.8675019d-8, 1.9131355d-1, -4.0112583d-1/

* coefficients for expansion of vap(1st 8 entries) and
* for va2p (entries 9:16)
      data xlam1/
     : 6.8251080d-1, 7.0563344d-1, 6.9646178d-1, 6.8106362d-1,
     : 6.6555984d-1, 6.6737903d-1, 6.7688238d-1, 7.0982932d-1,
     : 6.9166667d-1, 6.9098549d-1, 6.6827021d-1, 6.6654194d-1,
     : 6.4811548d-1, 6.2681388d-1, 6.8624117d-1, 7.0591451d-1/
      data xlam2/
     : 2.1788007, 2.5908122, 2.2148505, 2.1804703, 2.1100794, 2.1897645,
     : 1.7884042, 1.7023734, 2.1175000, 2.3200120, 2.0293552, 1.9423573,
     : 1.9999588, 1.9894804, 1.6173039, 1.7038259/
      data r0 /
     : 7.8133471, 7.6255886, 7.5355740, 7.5119147, 7.5896642, 7.5566877,
     : 7.4887842, 7.2944857, 7.7106944, 7.6534275, 7.8709667, 7.7710341,
     : 7.8059462, 7.8574895, 7.5236112, 7.3331908/
      data c1 /
     : -3.1811460d+4, -3.6206406d+4, -2.8380902d+4, -1.5668995d+4,
     : -1.1705039d+4, -1.2133446d+4, -1.4425077d+4, -2.0619980d+4,
     : -3.2567816d+4, -3.1034294d+4, -2.0798106d+4, -1.4911249d+4,
     : -1.1753416d+4, -1.0226346d+4, -1.7252103d+4, -1.9710314d+4/
      data c2/
     : 4.4270072d+8, 1.0721135d+8, 7.8568951d+7, -7.0160749d+7,
     : -4.2032021d+7, -4.5669775d+7, 5.7814359d+6, -7.1972605d+5,
     : 3.5370732d+8, 3.6791082d+8, 5.8967489d+7, -1.0897849d+7,
     : -4.0199947d+7, -2.4104713d+7, 5.2709064d+6, -6.0423782d+5/
      data c3/
     : -5.9760239d+7, 1.8894869d+8, 1.1304768d+7, 3.4083391d+7,
     : 2.1066307d+7, 2.3071235d+7, 3.2445967d+4, 8.4614357d+5,
     : -4.9473198d+7, -1.7292996d+7, 9.1794356d+5, 1.2994616d+7,
     : 2.3818167d+7, 1.6818622d+7, -7.7024882d+4, 8.2721962d+5/
       data clr/
     : -6.4938083d+6, -4.8046045d+6, -3.9463805d+6, -1.6003944d+6,
     : -1.6314313d+6, -1.8690112d+6, -3.1793799d+6, -3.8647226d+6,
     : -4.5310210d+6, -4.6327361d+6, -3.1430370d+6, -2.1421713d+6,
     : -2.5134226d+6, -4.0509923d+6, -3.9890903d+6, -3.8187158d+6/

* determine A' and A" potentials at angles
      rm3=one/r**3
      rm6=rm3*rm3
      do 100 i=1,8
        vap(i)=c1(i)*exp(-xlam1(i)*r) +
     :        (c2(i)+c3(i)*r)*exp(-xlam2(i)*r)-
     :         half*(tanh(alph*(r-r0(i)))+1)*clr(i)*rm6
        j=i+8
        va2p(i)=c1(j)*exp(-xlam1(j)*r) +
     :        (c2(j)+c3(j)*r)*exp(-xlam2(j)*r)-
     :         half*(tanh(alph*(r-r0(j)))+1)*clr(j)*rm6
        vsum(i)=half*(va2p(i)+vap(i))
* don't compute vdif for colinear geometries
        if (i.ne.1 .and. i.ne.8) then
          vdif(i-1)=half*(va2p(i)-vap(i))
* for long range damp out difference potential

          if (r .gt. rmax) then
            damp=-half*(tanh(3.d0*(r-11.5d0))-one)
            vdif(i-1)=vdif(i-1)*damp
          endif
        endif
100   continue
* solve simultaneous equations for solutions
* first for vsigma
      tol=1.e-10
      call dcopy(64,d0,1,aa,1)
      call dqrank(aa,8,8,8,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,8,8,8,kr,vsum,xsum,rsd,kpvt,qraux)
* convert to hartree 
      conv=1.d0/219474.6
      call dscal(8,conv,xsum,1)
      vv0=xsum(1)
      call dcopy(7,xsum(2),1,vvl,1)
      call dcopy(36,d2,1,aa,1)
      call dqrank(aa,6,6,6,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,6,6,6,kr,vdif,xdif,rsd,kpvt,qraux)
      call dscal(6,conv,xdif,1)
      call dcopy(6,xdif,1,vvl(8),1)
      end
