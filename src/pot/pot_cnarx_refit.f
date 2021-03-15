*system:  CN(Av=3,Xv=7)-Ar, refit of Berning's original PES's
*  this just calculates Vsigma (CN(X)+Ar) PES
*  references: Berning, A. (1995). 
*  Energieübertragungsprozesse in Atom-Molekule Stossen, Ph. D. thesis,
*  Universitaet Stuttgart.
*  Alexander, M. H., X. Yang, et al., "CN-Ar PES and low-J scattering,"
*  J. Chem. Phys. 112: 781-791 (2000)
*
      subroutine driver
      implicit double precision (a-h,o-z)
      character *48 potnam
      logical csflag, ljunk, ihomo, lljunk
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(7)
      potnam='ALEXANDER REFIT BERNING CN(X)-AR PESs'
      print *, potnam
      print *
1      print *, ' r (bohr) '
      read (5, *, end=93) r
      if (r .lt. 0.d0) go to 93
      call pot(vv0,r)
      write (6, 100) vvl
100   format('vsigma',7(1pe16.8))
      goto 1
93    r=4
      do i=1,25
       call pot(vv0,r)
       write(2,101) r,vvl
101    format(f8.4,7(1pe16.8))
       r=r+0.2
      enddo

99    end
* -----------------------------
*comdeck syusr
      subroutine syusr (irpot, readp, iread)
*  dummy syusr subroutine
      logical readp
      character*(*)fname
      if(.not.readp.or.iread.eq.0) then
        call loapot(1,' ')
        return
      endif
      entry ptrusr (fname)
      entry savusr (readp)
      entry chkusr
      return
      end
*comdeck ground
      subroutine ground(wf, r, nch, nphoto, mxphot)
      implicit double precision (a-h,o-z)
      entry wfintern(wf,yymin,nnvib,nny)
*  dummy ground subroutine
      return
      end
*comdeck bausr
* --------------------------------------------------------------------
      subroutine bausr (j, l, is, jhold, ehold, ishold, nlevel, nlevop,
     :                  sc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu,
     :                  csflag, clist, bastst, ihomo, nu, numin, jlpar,
     :                  n, nmax, ntop)
      return
      end
* --------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
*comdeck parbas
*  revised march 1992, c.r. 13-may-1997 by mha
      parameter (maxtrm=12,maxvib=10,maxvb2=maxvib**2)
*  variables in common block /cobspt/
*    lammin:   array containing minimum value of lambda for each term
*    lammax:   array containing maximum value of lambda for each term
*    mproj:    array containing the order of the reduced rotation matrix
*              elements for each term.  lammin can not be less than mproj.
*              for homonuclear molecules, the allowed values of lambda for
*              each term range from lammin to lammax in steps of 2
      common /cobsp2/ ntv(maxtrm),ivcol(maxvb2,maxtrm),
     :                ivrow(maxvb2,maxtrm)
      common /cobspt/ lammin(maxtrm), lammax(maxtrm), mproj(maxtrm)
      common /cobsptln/ lam2(maxtrm), m2proj(maxtrm)
*comdeck parpot
      character*48 potnam, label
      common /coptnm/ potnam, label
      common /coselb/ ibasty
      potnam='ALEXANDER REFIT BERNING CN(X)-AR PESs'
      ibasty=1
      lammin(1)=1
      lammax(1)=6
      mproj(1)=0
*     lammin(2)=0
*     lammax(2)=6
*     mproj(2)=0
*     lammin(3)=1
*     lammax(3)=5
*     mproj(3)=1
*     lammin(4)=2
*     lammax(4)=6
*     mproj(4)=2
*     ntv(1)=1
*     ntv(2)=1
*     ntv(3)=1
*     ntv(4)=1
*     ivrow(1,1)=7
*     ivrow(1,2)=3
*     ivrow(1,3)=7
*     ivrow(1,4)=3
*     ivcol(1,1)=7
*     ivcol(1,2)=3
*     ivcol(1,3)=3
*     ivcol(1,4)=3
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* --------------------------------------------------------------------------
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
*    vvl(1:7) expansion of vsig in legendre polynomials of order 0
* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  22-aug-2003
* ----------------------------------------------------------------------


      implicit double precision (a-h,o-z)
      logical ljunk, ihomo, csflag, lljunk
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(7)
      dimension vsigl1(7), vsigl2(7), vsigr0(7), vsigc1(7), vsigc2(7),
     :  vsigc3(7), vsigcl(7)

      dimension vsig(7)
      dimension d0(49),aa(49)
      dimension xsig(7),kpvt(7),
     :          qraux(7), work(50),rsd(7)
*    vvl(1-7) expansion coefficients in dl0 (l=0:6) of vsigma
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


      data rmax /13d0/
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
* convert to hartree
      econv=1./219474.6
      call dscal(24,econv,vvl,1)
      vv0=vv0*econv

      end
