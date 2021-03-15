*system:  OH(X2PI)+Ne CCSD(T) PES's
*reference:  M. Yang and M. H. Alexander, J. Chem. Phys. 103, 3400 (1995).


      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(9)
      include "common/parpot"
      potnam='YANG-ALEXANDER Ne-OH(X) CCSD[T] PES'      
      print *, potnam
1      print *, ' r (bohr)'
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vsum',/,6(1pe16.8),/,
     :    '  vdif',/,5e16.8)
      goto 1
99    end

      include "common/syusr"
      include "common/bausr"
      subroutine loapot(iunit,filnam)
* ------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      potnam='YANG-ALEXANDER Ne-OH(X) CCSD[T] PES'      
      lammin(1)=1
      lammax(1)=5
      lammin(2)=2
      lammax(2)=5
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
*  Ne-OH(X) potentials of Yang and Alexander
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
*    vvl(6-9) expansion coefficients in dl2 (l=2:5) of vdif

* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  20-dec-1995
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension xlam1(12),xlam2(12),r0(12),c1(12),c2(12),c3(12),
     :          clr(12),vsum(6),xsum(6),vdif(6),xdif(6),
     :          ddif(6),vap(6),va2p(6),
     :          d0(36),d2(16),aa(64)
      dimension kpvt(8),qraux(6),work(55),rsd(6)       

      common /covvl/ vvl(9)
      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
* for distances beyond rmax difference potential is damped
      data rmax /12d0/
* coefficicients for d0 rotation matrices
* stored (by column) for each of 6 angles and for l=0:5
* angles are 0 30 60 90 135 180
      data d0/
     : 1d0,  1d0,  1d0,  1d0,  1d0,  1d0,
     : 1d0, 8.6602541d-1, 5.0000001d-1, 1.3349125d-8,-7.071067d-1,-1d0,
     : 1d0, 6.2500001d-1,-1.2499999d-1,-5.0000000d-1, 2.499999d-1, 1d0,
     : 1d0, 3.2475954d-1,-4.3750000d-1,-2.0023687d-8, 1.767767d-1,-1d0,
     : 1d0, 2.3437511d-2,-2.8906251d-1, 3.7500000d-1,-4.062500d-1, 1d0,
     : 1d0,-2.2327216d-1, 8.9843733d-2, 2.5029609d-8, 3.756505d-1,-1d0/
* coefficicients for d2 rotation matrices
* stored (by column) for each of 4 angles and for l=2:5
* angles are 30 60 90 135
      data d2/
     : 1.5309311d-1,  4.5927932d-1,  6.1237244d-1,  3.06186217d-1,
     : 2.9646353d-1,  5.1348990d-1,  1.8279042d-8, -4.84122918d-1,
     : 4.1999000d-1,  2.2234766d-1, -3.9528471d-1,  4.94105884d-1,
     : 4.9023048d-1, -1.6982081d-1, -2.4180900d-8, -3.20217211d-1/

* coefficients for expansion of vap(1st 6 entries) and
* for va2p (entries 7:12)
      data xlam1/
     : 0.6376533, 0.6285311, 0.6360937, 0.4134702,    
     : 0.6540505, 0.6401758,                 
     : 0.6376533, 0.6295150, 0.6277373, 0.6359025,    
     : 0.6952167, 0.6401758/
      data xlam2/
     : 1.7967713, 1.7514748, 1.8014062, 1.4708627,                      
     : 1.7508831, 1.8531631,                                             
     : 1.7967713, 1.7523494, 1.7829534, 1.8110027,
     : 1.8462222, 1.8531631/            
      data r0 /
     : 8.0713904, 8.2406903, 8.6112500, 5.4169603,                      
     : 6.7743249, 10.5344656,                                           
     : 8.0713904, 8.2396051, 7.6037129, 7.3092139,
     : 10.057997, 10.5344656/         
      data c1 /
     : -3.53077731d+3, -2.56764176d+3, -2.10634765d+3, -0.29690723d+3,
     : -0.13363108d+3, -1.97361795d+3,                            
     : -3.53077731d+3, -2.72322312d+3, -2.29710833d+3, -2.32339971d+3,
     : -3.27460606d+3, -1.97361795d+3/
      data c2/
     :  3.5960481d+7,  1.9759046d+7,  1.1215005d+7,  0.2828463d+7,
     :  1.0820751d+7,  1.7167923d+7,                              
     :  3.5960481d+7,  2.0568998d+7,  1.3017215d+7,  1.2747974d+7,
     :  1.5925241d+7,  1.7167923d+7/                    
      data c3/
     : -5.42806586d+6, -3.13085373d+6, -1.89225782d+6, -0.54053520d+6,
     : -1.94539294d+6, -2.78974765d+6,                                
     : -5.42806586d+6, -3.12582261d+6, -1.84685057d+6, -1.70659135d+6,
     : -2.30383120d+6, -2.78974765d+6/                
       data clr/
     : -4.37496405d+5, -1.05204438d+5, -1.60800199d+5, -2.71842141d+6,
     :  2.25695878d+6,  6.86244259d+5,                               
     : -4.37496405d+5, -2.88043669d+5, -4.26491186d+5, -4.50596926d+5,
     :  9.20593447d+5,  6.86244259d+5/                

* determine A' and A" potentials at angles
      rm3=one/r**3
      rm6=rm3*rm3
      do 100 i=1,6
        vap(i)=c1(i)*exp(-xlam1(i)*r) +
     :        (c2(i)+c3(i)*r)*exp(-xlam2(i)*r)-
     :         half*(tanh(alph*(r-r0(i)))+1)*clr(i)*rm6
        j=i+6
        va2p(i)=c1(j)*exp(-xlam1(j)*r) +
     :        (c2(j)+c3(j)*r)*exp(-xlam2(j)*r)-
     :         half*(tanh(alph*(r-r0(j)))+1)*clr(j)*rm6
        vsum(i)=half*(va2p(i)+vap(i))
* don't compute vdif for colinear geometries
        if (i.ne.1 .and. i.ne.6) then
          vdif(i-1)=half*(-vap(i)+va2p(i))
* for long range damp out difference potential

          if (r .gt. rmax) then
            damp=-half*(tanh(3.d0*(r-rmax))-one)
            vdif(i-1)=vdif(i-1)*damp
          endif
        endif
100   continue
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
      call dcopy(16,d2,1,aa,1)
      call dqrank(aa,4,4,4,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,4,4,4,kr,vdif,xdif,rsd,kpvt,qraux)
      call dscal(4,conv,xdif,1)
      call dcopy(4,xdif,1,vvl(6),1)
      end
c  -----------------------------------------------------------------
      subroutine ground(wf, r, nch, nphoto, mxphot)
c  -----------------------------------------------------------------
*  driven equation source vector for Ne-OH(X) ccsd
*  author:  millard alexander
*  current revision date:  14-nov-1994
c  -----------------------------------------------------------------
*     variables in call list:

*     wf        array of dimension nch*nphoto, containing, on return,
*               ground state wavefunction in each of nch components
*               nphoto is number of difference ground state wavefunctions
*     r         value of separation coordinate
*     nch       total number of channels (row dimension of q)
*     nphoto    number of different wavefunctions calculated
*               column index of q vector
*     mxphot    maximum size of q vector (mxphot .ge. nch*nphoto)
*  -------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (nl =  52)
      dimension wf(36)
      common /cotq1/tmat(36)
      common /coipar/ jtot
      common /cosysi/ junk(5), isa
*  isa = 1, for v=0; isa = 2 for v=1 in lower lambda component
*  isa = -1, for v=0; isa = -2 for v=1 in lower lambda component
      common /coered/ ered, rmu
      dimension rl(nl), psi0(nl), psi1(nl), dpsi0(nl), dpsi1(nl)
* these are the spline knots and values for the v=0 and 1 adiabatic
* bender wavefunctions in the upper spin-orbit manifold
* these correspond to n=5 (diabatized)
      data rl /
     : 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2,
     : 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8,
     : 1.00d+1, 1.02d+1, 1.04d+1, 1.06d+1, 1.08d+1, 1.10d+1, 1.12d+1,
     : 1.14d+1, 1.16d+1, 1.18d+1, 1.20d+1, 1.22d+1, 1.24d+1, 1.26d+1,
     : 1.28d+1, 1.30d+1, 1.32d+1, 1.34d+1, 1.36d+1, 1.38d+1, 1.40d+1,
     : 1.42d+1, 1.44d+1, 1.46d+1, 1.48d+1, 1.50d+1 /
      data psi0 /
     : 9.3184411d-4, 5.1243555d-3, 2.0133373d-2, 6.0433611d-2,
     : 1.4289714d-1, 2.7447390d-1, 4.4231926d-1, 6.1490664d-1,
     : 7.5567855d-1, 8.3867174d-1, 8.5603754d-1, 8.1598009d-1,
     : 7.3550614d-1, 6.3337372d-1, 5.2544056d-1, 4.2277767d-1,
     : 3.3175409d-1, 2.5501958d-1, 1.9273372d-1, 1.4362754d-1,
     : 1.0578984d-1, 7.7163215d-2, 5.5823048d-2, 4.0106461d-2,
     : 2.8646830d-2, 2.0360119d-2, 1.4408751d-2, 1.0159601d-2,
     : 7.1410456d-3, 5.0055841d-3, 3.5003419d-3, 2.4426764d-3,
     : 1.7015179d-3, 1.1833832d-3, 8.2190189d-4, 5.7015848d-4,
     : 3.9510827d-4, 2.7355070d-4, 1.8923726d-4, 1.3081646d-4,
     : 9.0372226d-5, 6.2395395d-5, 4.3053954d-5, 2.9688195d-5,
     : 2.0453769d-5, 1.4072395d-5, 9.6578068d-6, 6.5956154d-6,
     : 4.4582717d-6, 2.9469800d-6, 1.8500216d-6, 1.0134445d-6/
      data psi1 /
     : 1.0869120d-3, 5.8071814d-3, 2.2021198d-2, 6.3274410d-2,
     : 1.4168630d-1, 2.5389092d-1, 3.7369727d-1, 4.5988954d-1,
     : 4.7622306d-1, 4.0834204d-1, 2.6722318d-1, 8.0800808d-2,
     : -1.1842014d-1, -3.0285207d-1, -4.5438977d-1, -5.6472963d-1,
     : -6.3333321d-1, -6.6463544d-1, -6.6554846d-1, -6.4361685d-1,
     : -6.0590393d-1, -5.5842063d-1, -5.0594229d-1, -4.5205430d-1,
     : -3.9929147d-1, -3.4932359d-1, -3.0314706d-1, -2.6127936d-1,
     : -2.2388873d-1, -1.9089016d-1, -1.6205483d-1, -1.3706492d-1,
     : -1.1555667d-1, -9.7153264d-2, -8.1484195d-2, -6.8199356d-2,
     : -5.6976546d-2, -4.7525179d-2, -3.9586939d-2, -3.2935308d-2,
     : -2.7373164d-2, -2.2730946d-2, -1.8862469d-2, -1.5643386d-2,
     : -1.2968443d-2, -1.0748771d-2, -8.9093535d-3, -7.3874558d-3,
     : -6.1302211d-3, -5.0938285d-3, -4.2417981d-3, -3.5439297d-3/
      data one, two /1.d0,2.d0/
      data interp / 1/
      if (interp .eq. 1) then
        if (isa .eq. 1 .or. isa .eq. 2) then
          write (6,100) isa-1
          write (9,100) isa-1
100       format  (' ** V = ',i2,' LOWER LAMBDA FOR NE-OH; ADIABATIC')
        elseif (isa .eq. -1 .or. isa .eq. -2) then
          write (6,101) -isa-1
          write (9,101) -isa-1
101       format  (' ** V = ',i2,' UPPER LAMBDA FOR NE-OH; ADIABATIC')
        else        
          write (6,102) isa
          write (9,102) isa
102       format  (' ** ISA = ',i2,' NOT ALLOWED FOR NE-OH; ABORT')
          call exit
        endif
*  here for first pass through source subroutine
*  interpolate the coefficients for spline expansion of adiabatic
*  bender wavefunctions
        nlr=nl-15
        nlm1 = nlr -1
        dpsi0f =( psi0(nlr) - psi0(nlm1)) / (rl(nlr) - rl(nlm1))
        dpsi0i =( psi0(2) - psi0(1)) / (rl(2) - rl(1))
        call dspline( rl, psi0, nlr, dpsi0f, dpsi0i, dpsi0)
        nlr=nl
        nlm1 = nlr -1
        dpsi1f =( psi1(nlr) - psi1(nlm1)) / (rl(nlr) - rl(nlm1))
        dpsi1i =( psi1(2) - psi1(1)) / (rl(2) - rl(1))
        call dspline( rl, psi1, nlr, dpsi1f, dpsi1i, dpsi1)
        interp=0
      endif
* determine wavefunction
      psifunc=0.d0
      if (iabs(isa) .eq. 1) then
         if (r .lt. 4.8d0 .or. r .gt. 12d0) then
           psifun=0.d0
         else
           call dsplint (rl,psi0(1),dpsi0(1),nl-15,r,psifunc)
         endif
      endif
      if (iabs(isa) .eq. 2) then
         if (r .lt. 4.8d0 .or. r .gt. 15d0) then
           psifunc=0.d0
         else
           call dsplint (rl,psi1(1),dpsi1(1),nl,r,psifunc)
         endif
      endif
* initial source vector is  (nb jmax=6)
* these are read down from maximum energy (for jmax=6)
*     jtot=0:  21st or 22nd eigenvector
*     jtot=1:  41st or 42nd eigenvector
*     jtot=2:  59th or 60th eigenvector
      if (jtot .eq. 0) then
        nvec0=22
        nvec1=21
      else if (jtot .eq. 1) then
        nvec0=42
        nvec1=41
      else if (jtot .eq. 2) then
        nvec0=60
        nvec1=59
      else
        write (6, 110) jtot
110     format (' *** JTOT = ',i3,' > 2 IN GROUND; ABORT ***')
        call exit
      endif
      nvec=nvec0
      if (isa .lt. 0) nvec=nvec1
*      if (r .le. 6.0882d0) then
*        nvec=nvec1
*        if (isa .lt. 0) nvec=nvec0
*      endif
      call dcopy(nch, tmat(nvec), nch, wf, 1)
      call dscal(nch,two*rmu*psifunc,wf,1)
      entry wfintern(wf,yymin,nnvib,nny)
      return
      end
