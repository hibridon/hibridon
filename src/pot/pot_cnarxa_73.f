*system:  CN(Av=3,Xv=7)-Ar, refit of Berning's original PES's
*references: Berning, A. (1995). 
*   Energieübertragungsprozesse in Atom-Molekül Stößen, Ph. D. thesis,Universität Stuttgart.
*   Alexander, M. H., X. Yang, et al. (2000). ³CN-Ar PES and low-J scattering.² J. Chem. Phys. 112: 781-791.
      include "common/syusr"
      include "common/ground"
      include "common/bausr"
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      common /coselb/ ibasty
      potnam='REFIT BERNING CN[X,v=7,A,v=3)-AR PESs'
      ibasty=4
      lammin(1)=0
      lammax(1)=6
      mproj(1)=0
      lammin(2)=0
      lammax(2)=6
      mproj(2)=0
      lammin(3)=1
      lammax(3)=5
      mproj(3)=1
      lammin(4)=2
      lammax(4)=6
      mproj(4)=2
      ntv(1)=1
      ntv(2)=1
      ntv(3)=1
      ntv(4)=1
      ivrow(1,1)=7
      ivrow(1,2)=3
      ivrow(1,3)=7
      ivrow(1,4)=3
      ivcol(1,1)=7
      ivcol(1,2)=3
      ivcol(1,3)=3
      ivcol(1,4)=3
      return
      end
* --------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      character *48 potnam
      logical csflag, ljunk, ihomo, lljunk
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(24)
      potnam='REFIT BERNING CN[X,v=7,A,v=3)-AR PESs'
      print *, potnam
      print *
1      print *, ' r (bohr) '
      read (5, *, end=93) r
      if (r .lt. 0.d0) go to 93
      call pot(vv0,r)
      write (6, 100) vvl
100   format(' vvl(0:6) vsigma',7(1pe16.8),/,
     :       ' vvl(7:14)   vpi',7(1pe16.8),/,
     :       ' vvl(15:19)   v1',5(1pe16.8),/,
     :       ' vvl(20:24)   v2',5(1pe16.8))
      goto 1
93    r=4
      do i=1,100
       call pot(vv0,r)
       write(2,101) r,vvl
101    format(f8.4,8(1pe16.8))
       r=r+0.2
      enddo

99    end
      subroutine pot (vv0, r)
* ----------------------------------------------------------------------
*  subroutine to calculate the r-dependent coefficients in the
*  cn(A,X) v=3,7 refit to Berning's potential
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*  vv0       dummy here
*  variable in common block /covvl/
*    vvl:     vector of length 24 to store r-dependence of each term
*             in potential expansion
*    vvl(1:7) expansion of vsig in legendre polynomials of order 0
*    vvl(8:14) expansion of vpi in legendre polynomials of order 0
*    vvl(15:19) expansion of v1 in legendre polynomials of order 1
*    vvl(20:24) expansion of v2 (vdif) in legendre polynomials of order 2
* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  24-sept-2001
* ----------------------------------------------------------------------


      implicit double precision (a-h,o-z)
      logical ljunk, ihomo, csflag, lljunk
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(24)
      dimension vsigl1(7), vsigl2(7), vsigr0(7), vsigc1(7), vsigc2(7),
     :  vsigc3(7), vsigcl(7)
      dimension xlam1(13),xlam2(13),r0(13),c1(13),c2(13),c3(13),
     :          clr(13)

      dimension v1l1(5), v1l2(5), v1c1(5), v1c2(5), v1c3(5), v1c4(5)
      dimension vsig(7), vap(7), va2p(7), vpi(7), v2(5),v1(5)
      dimension d0(49),d1(25),d2(25),aa(49)
      dimension xsig(7),xpi(7),x1(5),x2(5),kpvt(7),
     :          qraux(7), work(50),rsd(7)
*    vvl(1-7) expansion coefficients in dl0 (l=0:6) of vsigma
*    vvl(8:14) expansion coefficients in dl0 (l=0:6) of vpi
*    vvl(15:19) expansion coefficients in dl1 (l=1:5) of v1
*    vvl(20:24) expansion coefficients in dl2 (l=2:6) of v2
* angles are 0 30 60 90 120 150 180
      data zero, one, half /0.d0,1.d0,0.5d0/
* coefficients for d0 rotation matrices
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
* coefficients for d1 rotation matrices
* stored (by column) for each of 5 angles and for l=1:5
      data d1/
     : -3.5355339d-1, -6.1237244d-1, -7.0710678d-1, -6.1237244d-1,
     : -3.5355339d-1, -5.3033009d-1, -5.3033009d-1, 0d0, 5.3033009d-1,
     : 5.3033009d-1, -5.9539247d-1, -9.3750000d-2, 4.3301270d-1, 
     : -9.3750000d-2,
     : -5.9539247d-1, -5.4463828d-1, 3.0257682d-1, 0d0, -3.0257682d-1,
     : 5.4463828d-1, -3.9581513d-1, 3.5205044d-1, -3.4232660d-1, 
     :  3.5205044d-1,
     : -3.9581513d-1/
* coefficients for d2 rotation matrices
* stored (by column) for each of 5 angles and for l=2:6
      data d2/
     : 1.5309311d-1,4.5927933d-1,6.1237244d-1,4.5927933d-1,1.5309311d-1,
     : 2.9646353d-1, 5.1348990d-1, -0d0, -5.1348990d-1, -2.9646353d-1,
     : 4.1999000d-1, 2.2234765d-1, -3.9528471d-1, 2.2234765d-1, 
     :  4.1999000d-1,
     : 4.9023048d-1, -1.6982082d-1, 0d0, 1.6982082d-1, -4.9023048d-1,
     : 4.8532921d-1, -3.4523418d-1, 3.2021721d-1, -3.4523418d-1, 
     : 4.8532921d-1/
* hyperbolic tangent scaling factor
      data alph /1.2d0/
*  expansion coefficients for vsigma
*  by rows: lam1 lam2 r0 c1 c2 c3 cl
      data vsigl1/
     : 5.9952171d-1,6.9875525d-1,4.8851642d-1,4.7308802d-1,4.9950844d-1,
     : 5.5639041d-1, 3.1156829d-1/
      data vsigl2/
     : 1.4022638d0, 1.4446592d0, 1.5660267d0, 1.4188301d0, 1.3817310d0,
     : 1.3794546d0, 1.3023817d0/
      data vsigr0/
     : 6.0737174d0, 5.7846066d0, 5.4523692d0, 5.4828829d0, 4.4811091d0,
     : 6.1865680d0, 6.6083840d0/
      data vsigc1/
     : 3.5273093d3,-3.1710258d3,8.7879793d3,-1.4667112d3,-4.0320752d3,
     : -7.5859233d3, -1.6775191d2/
      data vsigc2/
     : 3.3980592d7, 3.3816379d7, -1.6003967d7, 1.3650140d7, 1.6349294d7,
     : 2.9615955d7, 2.2847126d7/
      data vsigc3/
     : -4.1434101d6,-4.1646824d6,5.8581572d6,-1.9412723d6,-2.3894461d6,
     : -3.8516599d6, -2.8684311d6/
      data vsigcl/
     : 3.0804905d7, 1.8297311d7, 9.3572419d7, 1.8364610d6, -1.1988674d7,
     : -7.6733061d6, 1.4427108d7/

* coefficients for expansion of vpiap(1st 7 entries) and
* for vpia2p (entries 8:13)
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

*  expansion coefficients for v1
*  by rows: lam1 lam2 c1 c2 c3
      data v1l1/
     : 1.6686230d0, 1.7586914d0, 1.8883447d0, 1.7168457d0, 3.1427425d0/
      data v1l2/
     : 0d0, 0d0, 0d0, 0d0, 1.5463528d0/
      data v1c1/
     : 0d0, 0d0, 0d0, 0d0, 6.0548759e+09/
      data v1c2/
     : 2.1684873d7, -1.9732962d7, -3.4766127d7, -2.4953277d7, 
     : 1.9205265d7/ 
      data v1c3/
     : -2.0784264d6, 9.1259265d6, 1.2176257d7, 1.1599075d7, 
     : -1.8469688d6/
      data v1c4/
     : 0d0, -7.1671651d5, -6.1166906d5, -8.5531639d5, 0d0/

      data rmax /13d0/
* v=3-v=7 franck-condon factor
      data fcfact / 0.17921d0/
* determine potentials at angles
      rm1=one/r
      rm2=rm1*rm1
      rm3=one/r**3
      rm6=rm3*rm3

* cm and rad are conversion factors from au to cm-1
* and from degrees to radians respectively
      cm=2.194746d05
      pi=acos(-1.0d0)
      rad=pi/180.0d0
      do 200 i=1,7
        ang=(i-1)*30*rad
        vsig(i)=vsigc1(i)*dexp(-vsigl1(i)*r)+
     :        (vsigc2(i)+vsigc3(i)*r)*dexp(-vsigl2(i)*r)-
     :        half*(tanh(alph*(r-vsigr0(i)))+one)*vsigcl(i)*rm6
        vap(i)=c1(i)*exp(-xlam1(i)*r) +
     :        (c2(i)+c3(i)*r)*exp(-xlam2(i)*r)-
     :         half*(tanh(alph*(r-r0(i)))+1)*clr(i)*rm6
        j=i+7
* don't compute va2p and vdif for colinear geometries
        if (i.ne. 1. and. i .ne. 7 ) then
           va2p(i)=c1(j)*exp(-xlam1(j)*r) +
     :        (c2(j)+c3(j)*r)*exp(-xlam2(j)*r)-
     :         half*(tanh(alph*(r-r0(j)))+1)*clr(j)*rm6
           vpi(i)=half*(va2p(i)+vap(i))
           v2(i-1)=half*(va2p(i)-vap(i))

        else
           vpi(i)=vap(i)
        endif

        if (i.gt.1 .and. i.lt.7) then
* here for non-colinear angles
           j=i-1
* fit to v1, if lam(2)=0, then single exponential
           if (v1l2(j).eq.0d0) then
             v1(j)=exp(-v1l1(j)*r)*(v1c2(j)+v1c3(j)*r+v1c4(j)*r*r)
* otherwise double exponential
           else
               v1(j)=v1c1(j)*dexp(-v1l1(j)*r)+
     :        (v1c2(j)+v1c3(j)*r)*dexp(-v1l2(j)*r)
           endif
* for long range damp out difference potential and v1 potential
           if (r .gt. rmax) then
              damp=-half*(tanh(3.d0*(r-rmax))-one)
              v2(j)=v2(j)*damp
              v1(j)=v1(j)*damp
           endif
        endif
200   continue
      sum=0
* determine matrix of d's at angles for least-squares fit
*      call d2matev(9,beta,d0,d1,d2)
* (n.b. these are already in common, subroutine is included
*  only for completeness)
      vv0=0d0
* solve simultaneous equations for solutions
* first for vsig
      tol=1.e-10
      call dcopy(49,d0,1,aa,1)
      call dqrank(aa,7,7,7,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,7,7,7,kr,vsig,xsig,rsd,kpvt,qraux)
      call dcopy(7,xsig,1,vvl,1)
* then for vpi
      call dqrlss(aa,7,7,7,kr,vpi,xpi,rsd,kpvt,qraux)
      call dcopy(7,xpi,1,vvl(8),1)
* then for v1
      call dcopy(25,d1,1,aa,1)
      call dqrank(aa,5,5,5,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,5,5,5,kr,v1,x1,rsd,kpvt,qraux)
* multiply by franck-condon factor
      call dscal(5,fcfact,x1,1)
      call dcopy(5,x1,1,vvl(15),1)
* then for v2
      call dcopy(25,d2,1,aa,1)
      call dqrank(aa,5,5,5,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,5,5,5,kr,v2,x2,rsd,kpvt,qraux)
      call dcopy(5,x2,1,vvl(20),1)
* convert to hartree
      econv=1./219474.6
      call dscal(24,econv,vvl,1)
      vv0=vv0*econv

      end
