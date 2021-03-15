*  System:  O(3P) + H(2S)
*
*   calculation of potential energy curves and spin-orbit matrix
*   elements, taken from m. alexander's notes (2018)
*
*   MRCI+Q calculation of richard dawes, nov 2018
*
*   two spin-orbit matrix parameters are included:  apar and aperp
*   the atomic spin-orbit splittings are computed from the large-R
*   values of these parameters.
*   here, the atomic values (apar = aperp) are used for all R
*  
*   written by p. dagdigian
*   current revision date:  20-oct-2018
*
      include "common/syusr"
      include "common/bausr"
* --------------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      common /coselb/ ibasty
      include "common/parbas"
      include "common/parpot"
      common /conlam/ nlam, nlammx, lamnum(1)
      potnam='O(3P)-H(2S) DAWES 2018'
      npot=1
      ibasty=23
      lammin(1)=1
      lammax(1)=6
      mproj(1)=0
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
*  count number of anisotropic terms
      nlam = 0
      do il = lammin(1), lammax(1)
        nlam = nlam + 1
      enddo
      return
      end
* --------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysi/ junk(5), npot
      common /covvl/ vvl(8)
      include "common/parpot"
      dimension wf(9)
      potnam='O(3P)-H(2S) DAWES 2018'
      econv=219474.6d0
      print *, potnam
1     print *, 'R (bohr):'
      read (5, *, end=99) r
      if (r.le.0.d0) goto 99
      call pot(vv0,r)
      write (6, 100) (econv*vvl(i), i=1,6)
      nphoto=1
      mxphot=9
      call ground(wf, r, nch, nphoto, mxphot)
*  here we print out the vector (dimension 9) wf
      write (6,100) (wf(i), i=1,9)
      
100   format('  ',9(1pe16.8))
      goto 1
99    rr= 1.075d0
      dr=0.2d0
      open (unit=12,file='sh_vlms_dawes.txt')
      write(12,109)
109   format(' %R/bohr X2Pi  4Pi  2Sig-  4Sig-  apar   aperp')
      do i=1,426
        call pot(vv0,rr)
        write (12,110) rr,(econv*vvl(j),j=1,6)
110     format(f11.5,6(1pe16.8))
        rr = rr + dr
      enddo
      close(12)
      end
* --------------------------------------------------------------------------
      subroutine pot (vv0, r)
*  subroutine to calculate the r-dependent coefficients
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*    vv0        (not used)
*    vvl(1,4)   contains the X2Pi, 4Pi, 2Sig-, 4Sig- energies
*    vvl(5,6)   contains the spin-orbit matrix parameters apar and aperp at R
*  variable in common block /conlam/ used here
*    nlam:      the number of angular coupling terms actually used
*  variable in common block /covvl/
*    vvl:       array to store r-dependence of each
*               angular term in the potential
*  variable in common block /coconv/
*    econv:     conversion factor from cm-1 to hartrees
* 
* author:  paul dagdigian and millard alexander
* latest revision date:  30-dec-2018
* ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical ljunk, chlist
      include "common/parbas"
      common /covvl/ vvl(6)
      common /colpar/ ljunk(4), chlist
      parameter (npsi = 157)
      dimension rr(51),v2pi(51),v4pi(51),v2sigm(51),v4sigm(51),
     :          cs2p(51), cs4p(51), cs2s(51), cs4s(51)
      dimension asr(5), bsr(5), clr(5), dV_sr(5), dv_lr(5)
      dimension vl(51), vec(51), v(6)

      dimension d1(4), d2(4), vlr(4)
      data npts / 51/
*
c  number of r values= 51
      data rr/
     + 9.0000d-1, 9.9109d-1, 1.0822d0, 1.1733d0, 1.2644d0, 1.3555d0,
     + 1.4465d0, 1.5376d0, 1.6287d0, 1.7198d0, 1.8109d0, 1.9020d0,
     + 1.9931d0, 2.0842d0, 2.1753d0, 2.2664d0, 2.3575d0, 2.4485d0,
     + 2.5396d0, 2.6307d0, 2.7218d0, 2.8129d0, 2.9040d0, 2.9951d0,
     + 3.0862d0, 3.1773d0, 3.2684d0, 3.3595d0, 3.4506d0, 3.5416d0,
     + 3.6327d0, 3.7238d0, 3.9060d0, 4.0882d0, 4.2704d0, 4.4526d0,
     + 4.6347d0, 4.8169d0, 5.2724d0, 5.6367d0, 6.0011d0, 6.3655d0,
     + 6.7298d0, 7.0942d0, 7.4586d0, 7.8229d0, 8.1873d0, 8.5517d0,
     + 8.9160d0, 9.2804d0, 9.6447d0/
c X2Pi
       data v2pi/ 
     + 1.6737356d5, 1.0231793d5, 5.6470682d4, 2.4395192d4, 1.9596941d3,
     + -1.3312068d4, -2.3637825d4, -3.0459258d4, -3.4515605d4,
     + -3.6578637d4, -3.7443365d4, -3.7149876d4, -3.6338257d4,
     + -3.4863835d4, -3.3111052d4, -3.1145104d4, -2.9043449d4,
     + -2.6885533d4, -2.4726344d4, -2.2596815d4, -2.0526965d4,
     + -1.8544092d4, -1.6660392d4, -1.4885294d4, -1.3231326d4,
     + -1.1701756d4, -1.0295944d4, -9.0171090d3, -7.8634916d3,
     + -6.8284414d3, -5.9078165d3, -5.0966428d3, -3.7636006d3,
     + -2.7628616d3, -2.0223506d3, -1.4808657d3, -1.0870202d3,
     + -8.0040990d2, -3.7971881d2, -2.1406920d2, -1.2354677d2,
     + -7.3047087d1, -4.4369826d1, -2.7787788d1, -1.7981324d1,
     + -1.2040085d1, -8.3205418d0, -5.9616288d0, -4.3832802d0,
     + -3.2969295d0, -2.5273773d0/
c 4Sigma-
      data v4sigm/
     + 2.4751421d5, 1.8139138d5, 1.3439404d5, 1.0112447d5, 7.7432732d4,
     + 6.0706176d4, 4.8720803d4, 4.0090637d4, 3.3877595d4, 2.9224464d4,
     + 2.5619846d4, 2.2641351d4, 2.0100238d4, 1.7862986d4, 1.5771186d4,
     + 1.3856426d4, 1.2115193d4, 1.0546096d4, 9.1510955d3, 7.9179511d3,
     + 6.8356141d3, 5.8941348d3, 5.0764433d3, 4.3674508d3, 3.7556713d3,
     + 3.2275277d3, 2.7710105d3, 2.3772374d3, 2.0374597d3, 1.7436381d3,
     + 1.4896922d3, 1.2705000d3, 9.1702975d2, 6.5387086d2, 4.5846367d2,
     + 3.1459108d2, 2.0974484d2, 1.3414048d2, 2.9410855d1,
     + -4.5874317d0, -1.7343722d1, -1.9999141d1, -1.8505908d1,
     + -1.5506654d1, -1.2379542d1, -9.6524466d0, -7.4591150d0,
     + -5.7558797d0, -4.4458407d0, -3.4479561d0, -2.6865380d0/
c 2Sigma-
       data v2sigm/
     + 2.5024247d5, 1.8418478d5, 1.3726391d5, 1.0408313d5, 8.0497109d4,
     + 6.3913854d4, 5.2121706d4, 4.3736318d4, 3.7822410d4, 3.3515841d4,
     + 3.0293856d4, 2.7721333d4, 2.5509078d4, 2.3502760d4, 2.1576630d4,
     + 1.9717792d4, 1.7923702d4, 1.6209201d4, 1.4596628d4, 1.3092983d4,
     + 1.1704412d4, 1.0435491d4, 9.2808998d3, 8.2346356d3, 7.2924554d3,
     + 6.4459730d3, 5.6865615d3, 5.0080279d3, 4.4031018d3, 3.8640594d3,
     + 3.3849887d3, 2.9604874d3, 2.2517096d3, 1.7003849d3, 1.2739011d3,
     + 9.4669448d2, 6.9755538d2, 5.0901962d2, 2.2060717d2, 1.0541064d2,
     + 4.5257609d1, 1.5216588d1, 1.1901581d0, -4.6124432d0,
     + -6.4052253d0, -6.3973047d0, -5.6806870d0, -4.8116157d0,
     + -3.9449311d0, -3.1858233d0, -2.5548735d0/ 
c 4Pi
      data v4pi/
     + 2.8720210d5, 2.2053430d5, 1.7294151d5, 1.3903277d5, 1.1463528d5,
     + 9.7023294d4, 8.3977159d4, 7.4127647d4, 6.6343836d4, 5.9844107d4,
     + 5.4148381d4, 4.8866766d4, 4.3917074d4, 3.9246646d4, 3.4879257d4,
     + 3.0851757d4, 2.7168453d4, 2.3837493d4, 2.0855787d4, 1.8197517d4,
     + 1.5839402d4, 1.3762078d4, 1.1934938d4, 1.0330530d4, 8.9290952d3,
     + 7.7062281d3, 6.6393077d3, 5.7121290d3, 4.9077430d3, 4.2095157d3,
     + 3.6049103d3, 3.0829839d3, 2.2428907d3, 1.6209023d3, 1.1620751d3,
     + 8.2607008d2, 5.8164257d2, 4.0471662d2, 1.5279154d2, 6.2270423d1,
     + 1.9670639d1, 1.0002616d0, -6.1864029d0, -8.1392314d0,
     + -7.8728743d0, -6.8386488d0, -5.6589243d0, -4.5761185d0,
     + -3.6663787d0, -2.9304163d0, -2.3450201d0/
* short and long range extrapolations for potential and derivatives
      data asr /
     + 13041115.93d0, 3757266.60d0, 4697775.11d0, 4813843.52d0, 
     + 7160354.38d0/
      data bsr/
     + 4.8396d0, 2.8570d0, 3.2582d0, 3.2975d0, 3.9382d0/
      data clr/
* 2Pi, 4Pi, 2Sigm, 4Sigm, 2Sig+
     + 2.3081d6, 1.8397d6, 2.0126d6, 2.3515d6, 2.1319d6/
      data dv_lr/
     +  1.7790d0, 1.4172d0, 1.5504d0, 1.8115d0, 1.6423d0/
      data dv_sr/ 
     + -8.1002d5, -8.2052d5, -8.1535d5, -8.1619d5, -8.1459d5/
      data dv_lr/
     +  1.7d0, 1.4219d0, 1.5555d0, 1.8115d0, 1.6423d0/
      data one, half /1d0,0.5d0/
*  set vv0 term to zero
      vv0 = 0.d0
*
*  spline fits
      if (ifirst .eq. 0) then
*  spline fit of the coefficients for the PE curves
*  using calculated derivatives
         do i=1,4
            der1=dv_sr(i)
            der2=dv_lr(i)
            if (i.eq.1)
     :         call dspline(rr,v2pi,npts,der1,der2,cs2p)
            if (i.eq.2)
     :         call dspline(rr,v4pi,npts,der1,der2,cs4p)
            if (i.eq.3)
     :         call dspline(rr,v2sigm,npts,der1,der2,cs2s)
            if (i.eq.4)
     :         call dspline(rr,v4sigm,npts,der1,der2,cs4s)
         enddo
*
        ifirst = 1
 
      endif
*
*  determine splined or extrapolated coefficients for PE curves at r=R
*  note: ordering is 1-4:  2Pi, 4Pi, 2Sig-, 4sig-
      if (r.lt.rr(1)) then
         do i=1,4 
            v(i)=asr(i)*exp(-bsr(i)*r)
         enddo
      else
*  spline term
         call dsplint(rr,v2pi,cs2p,npts,r,vvx)
         v(1)=vvx
         call dsplint(rr,v4pi,cs4p,npts,r,vvx)
         v(2)=vvx
         call dsplint(rr,v2sigm,cs2s,npts,r,vvx)
         v(3)=vvx
         call dsplint(rr,v4sigm,cs4s,npts,r,vvx)
         v(4) = vvx
*  long-range contribution
         do i=1,4 
		    vlr(i) = -clr(i)/r**6
         enddo
*  merge function
         alpha=1d0
         xmerge= -half*(tanh(alpha*(r-rr(npts)))-one)
         do i=1,4
           v(i)=xmerge*v(i)+(one-xmerge)*vlr(i)
         enddo
      endif
*
*  IN THIS VERSION, SET SPIN-ORBIT MATRIX ELEMENTS TO THEIR
*  ASYMPTOTIC VALUES AT ALL R
*  CBS value 2(E(J=0)-E(J=1))
      a = 151.2
      apar = a / 2.0d0
      aperp = a / 2.0d0
      v(5) = apar
      v(6) = aperp
*
*      if (chlist) then
*         write (9,120) r, xmerge, (v(i),i=1,4)
*120      format(7e15.5)
*      endif
      call dcopy(6,v,1,vvl,1)
* convert to hartree
      econv=1.d0/219474.6d0
      call dscal(6,econv,vvl,1)
      return
      end
c  -----------------------------------------------------------------
      subroutine ground(wf, r, nch, nphoto, mxphot)
c  -----------------------------------------------------------------
*  driven equation source vector for oh

*  current revision date: 30-dec-2018
c  -----------------------------------------------------------------
*     variables in call list:

*     wf        array of dimension nch*nphoto, containing, on return,
*               ground state wavefunction in each of nch components
*               nphoto is number of different ground state wavefunctions
*     r         value of separation coordinate
*     nch       total number of channels (row dimension of q)
*     nphoto    number of different wavefunctions calculated
*               column index of q vector
*     mxphot    maximum size of q vector (mxphot .ge. nch*nphoto)
*  -------------------------------------------------------
*  variable in common block /cosysi/
*    nvib:     initial vibrational level in A state
*  variable in common block /coj12/
*    j12:      array containing vector sum of ja + j (similar
*              here j12 is lower-case j in alexander, dagdigian, klos notes
*              situation with molecule-molecule scattering,
*              see hibastpln basis routine)
*  variable in commmon block /coja/
*     jja:      atomic 3P atom electronic angular momentum for each channel for 3P+2S scattering
*  variable in common block /coel/
*     ll:       end-over-end HS orbital angular momentum
*  variable in common block /corv_A/
*     rpsi:     knots for spline integration of A state vibrationaal wavefunctions
*               values given in blockdata psiv
*  variable in common block /copsiv_A/
*     psiv_0 -- psiv_7  vibrational wavefunctions for all 3 vibrational levels of the A state
*               values at r=rpsi given in blockdata psiv
*               wavefunctions normalized int(psiv_n*psiv_n)=1
      implicit double precision (a-h,o-z)
      logical grndtst, grndtst1, ljunk
      common /cosysi/ junk(3), nvib
      common /colpar/ ljunk(12), noprin
*  grndtst=.true.  if noprin=.false.
      common /cojtot/ jjtot, jjlpar
      common /coj12/  j12(9)
      common /coja/  jja(9)
      common /coel/  ll(9)
      common /corv_A/ rpsi(157)
      common /copsiv_A/ 
     + psiv_0(157), psiv_1(157), psiv_2(157), psiv_3(157),
     + psiv_4(157), psiv_5(157), psiv_6(157), psiv_7(157),
     + psiv_8(157)
      dimension psi_v(157)
      dimension c2sig(5), c4sig(5), cpi(5)
*  9 states of a given parity
      dimension wf(9), sc(9), const_pi(9), const_2sig(9), 
     +   const_4sig(9), term(9),
     +   term1(9), const(9)
      dimension rrh(39), hpi(39), h2sigm(39), h4sigm(39), 
     +   csh2sigm(39), csh4sigm(39), cshpi(39)
      dimension cspsi(157)
      common /coered/ ered, rmu
      econv=1.d0/219474.6d0
      data nl, npsi / 39, 157/
*  knots for spline expansion of Hsigm and Hpi spin-orbit coupling terms
      data rrh/
     + 7.5589d-1, 9.4486d-1, 1.1338d0, 1.3228d0, 1.5118d0, 1.7008d0,
     + 1.8897d0, 2.0787d0, 2.2677d0, 2.4566d0, 2.5322d0, 2.6456d0,
     + 2.7401d0, 2.8346d0, 2.9291d0, 3.0236d0, 3.1180d0, 3.2125d0,
     + 3.4015d0, 3.5905d0, 3.7795d0, 3.9684d0, 4.1574d0, 4.5353d0,
     + 4.9133d0, 5.2912d0, 5.6692d0, 6.0471d0, 6.4251d0, 6.8030d0,
     + 7.5589d0, 8.3148d0, 9.0707d0, 1.0393d1, 1.1338d1,
     + 1.2283d1, 1.3228d1, 1.4173d1, 1.5118d1/
c H2sigm -i <^2\Sigma^-(xy),1/2|Hso|A^2\Sigma^+,1/2> from av5z calculations
      data h2sigm /
     + 2.3157d0, 3.2162d0, 3.9084d0, 5.4456d0, 8.0002d0, 1.1876d+1,
     + 1.7218d+1, 2.1466d+1, 2.4120d+1, 2.5727d+1, 2.6134d+1,
     + 2.6510d+1, 2.6619d+1, 2.6551d+1, 2.6320d+1, 2.5945d+1,
     + 2.5444d+1, 2.4868d+1, 2.3588d+1, 2.2294d+1, 2.1131d+1,
     + 2.0153d+1, 1.9362d+1, 1.8264d+1, 1.7620d+1, 1.7263d+1,
     + 1.7067d+1, 1.6965d+1, 1.6912d+1, 1.6886d+1, 1.6864d+1,
     + 1.6855d+1, 1.6855d+1, 1.6857d+1, 1.6857d+1, 1.6857d+1,
     + 1.6856d+1, 1.6855d+1, 1.6855d+1/
c H4sigm -i <^4\Sigma^-(xy),1/2|Hso|A^2\Sigma^+,1/2> from av5z calculations
      data h4sigm /
     + 4.0113d0, 5.5704d0, 6.7683d0, 9.4287d0, 1.3849d+1, 2.0552d+1,
     + 2.9835d+1, 3.7180d+1, 4.1776d+1, 4.4560d+1, 4.5264d+1,
     + 4.5916d+1, 4.6104d+1, 4.5985d+1, 4.5583d+1, 4.4932d+1,
     + 4.4065d+1, 4.3066d+1, 4.0847d+1, 3.8605d+1, 3.6591d+1,
     + 3.4899d+1, 3.3530d+1, 3.1628d+1, 3.0516d+1, 2.9898d+1,
     + 2.9561d+1, 2.9383d+1, 2.9293d+1, 2.9247d+1, 2.9209d+1,
     + 2.9194d+1, 2.9193d+1, 2.9197d+1, 2.9198d+1, 2.9197d+1,
     + 2.9195d+1, 2.9195d+1, 2.9194d+1/
c Hpi -\sqrt 2 < ^4\Pi_x,1/2|Hso| A,-1/2> from av5z calculation
      data hpi /
     + 3.2748d0, 4.5484d0, 5.5272d0, 7.7012d0, 1.1314d1, 1.6795d1,
     + 2.4350d1, 3.0358d1, 3.4111d1, 3.6384d1, 3.6959d1, 3.7491d1,
     + 3.7646d1, 3.7549d1, 3.7222d1, 3.6691d1, 3.5984d1, 3.5169d1,
     + 3.3358d1, 3.1528d1, 2.9883d1, 2.8501d1, 2.7383d1, 2.5829d1,
     + 2.4919d1, 2.4413d1, 2.4137d1, 2.3992d1, 2.3918d1, 2.3880d1,
     + 2.3849d1, 2.3837d1, 2.3836d1, 2.3839d1, 2.3840d1, 2.3839d1,
     + 2.3838d1, 2.3837d1, 2.3837d1/
* CPi,C2sig, and C4sig(Eqs. 248-252)
      data cpi /
     +  0d0, -7.4536d-1, -1.6667d0, 4.7140d-1, 6.66667d-1/
* change sign from previous versions 11/24/2018
      data c2sig/
     +  0d0,  -7.4536d-1, 3.33333d-1, 4.7140d-1, -3.33333d-1/
      data c4sig/
     +  0.7746d0, -0.10541d0, -0.23570d0, -0.333333d0, -0.47105d0/
      data zero, half, one, two /0d0, 0.5d0, 1.d0,2.d0/
*
* conversion to hartree
      econv=1.d0/219474.6d0
      sq2=sqrt(2d0)

      grndtst1=.false.
      grndtst=.not.noprin
      
      
      interp=1
      if (interp .eq. 1) then
*  here for first pass through source subroutine
         if (grndtst) then
             write(6,14)
14           format(/,'** GROUND TEST, SEE OUTPT',/)
         endif
*  interpolate the coefficients for spline expansion of spin-orbit coupling matrix elts
*  and vibrational wavefunction in A state
*  assume zero derivatives at initial and final endpoints
        
         call dspline(rrh,h2sigm,nl,0d0,0d0,csh2sigm)
         call dspline(rrh,h4sigm,nl,0d0,0d0,csh4sigm)
         call dspline(rrh,hpi,nl,0d0,0d0,cshpi)
* now vibrational wavefunction 
         if (nvib.eq.0) call dcopy(npsi,psiv_0,1,psi_v,1)
         if (nvib.eq.1) call dcopy(npsi,psiv_1,1,psi_v,1)
         if (nvib.eq.2) call dcopy(npsi,psiv_2,1,psi_v,1)
         if (nvib.eq.3) call dcopy(npsi,psiv_3,1,psi_v,1)
         if (nvib.eq.4) call dcopy(npsi,psiv_4,1,psi_v,1)
         if (nvib.eq.5) call dcopy(npsi,psiv_5,1,psi_v,1)
         if (nvib.eq.6) call dcopy(npsi,psiv_6,1,psi_v,1)
         if (nvib.eq.7) call dcopy(npsi,psiv_7,1,psi_v,1)
         if (nvib.eq.8) call dcopy(npsi,psiv_8,1,psi_v,1)
         if ((nvib.lt.0).or.(nvib.gt.8)) then
            write (6,15) nvib
            write (9,15) nvib
15            format(/'nvib =',i3,' OUT OF RANGE, ABORT ***')
            call exit
         endif
***  write out wf array (for testing purposes)
         if (grndtst1) then
            write(9,*) 'GROUND VIBRATIONAL WAVEFUNCTION'
            do irr = 1, 157
               write(9,336) rpsi(irr),psi_v(irr)
336            format(2(1pe16.6))
            enddo
            write(9,*) 'WAVEFUNCTION END'
         endif

         call dspline(rpsi,psi_v,npsi,0d0,0d0,cspsi)
          
* sum over all channels to get the pi and sigma coefficients of Eq. 229
         if (grndtst) then
             write(9,*) 'JJTOT    JJA    J12    L'
             do ich=1,nch
                write (9,30) jjtot, jja(ich), j12(ich), ll(ich) 
30              format(4i4)
             enddo
         endif
         do i=1,nch
             xll=ll(i)
             xj12=j12(i)+0.5
             xj=jjtot+0.5
             if (grndtst) write (9,34) i, ll(i), jjtot
34           format(3i3)
                const1=sq2*sqrt(2*ll(i)+1d0)*((-1)**(jjtot+1))
                const3j= xf3j(xll,xj12,xj,zero,-half,half)
                const(i)=const1*const3j
                if (grndtst) then
                   write(9,35) xll, xj12, xj, const1, const3j, const(i)
35                 format('3J PARAMETERS:',/,3f6.2,3g16.8)
                endif
                if (jja(i).eq.2) then
                   if (j12(i).eq.2) then
                      const_pi(i)=cpi(1)
                      const_2sig(i)=c2sig(1)
                      const_4sig(i)=c4sig(1)
                   elseif (j12(i).eq.1) then
                      const_pi(i)=cpi(2)
                      const_2sig(i)=c2sig(2)
                      const_4sig(i)=c4sig(2)
                   else
                      write (6,20) i, jja(i), j12(i), ll(i)
                       write (9,20) i, jja(i), j12(i), ll(i)
20                    format(/'WRONG VALUES OF ANG MOMENTA IN GROUND',
     +                       ', ICH, JJA, J!2, LL: ',i3,2i2,i4) 
                      call exit
                   endif
                elseif (jja(i).eq.1) then
                   if (j12(i).eq.1) then
                      const_pi(i)=cpi(3)
                      const_2sig(i)=c2sig(3)
                      const_4sig(i)=c4sig(3)
                   elseif (j12(i).eq.0) then
                      const_pi(i)=cpi(4)
                      const_2sig(i)=c2sig(4)
                      const_4sig(i)=c4sig(4)
                   else
                      write (6,20) i, jja(i), j12(i), ll(i)
                      write (9,20) i, jja(i), j12(i), ll(i)
                      call exit
                   endif
                elseif (jja(i).eq.0) then
                   if (j12(i).eq.0) then
                      const_pi(i)=cpi(5)
                      const_2sig(i)=c2sig(5)
                      const_4sig(i)=c4sig(5)
                   else
                      write (9,20) i, jja(i), j12(i), ll(i)
                      call exit
                   endif
                else
                      write (6,20) i, jja(i), j12(i), ll(i)
                      write (9,20) i, jja(i), j12(i), ll(i)
                      call exit
                endif

         enddo
         if (grndtst) then
            write (9,*)
     :  ' I    CONST      CONST_PI    CONST_2SIG     CONST_4SIG'
            do i=1,nch
                write(9,21) i, const(i), const_pi(i), const_2sig(i), 
     :            const_4sig(i)
21              format(i3,4g14.4)
            enddo
         endif
         interp=1
      endif
* rmlmda is 2mu/hbar^2 in atomic units
      if (grndtst1) then
         xmuconv=1d0/5.4857990d-4
         rmlmda = two*xmuconv*(32*1)/33
      else
         rmlmda = two*rmu
      endif
* determine A state vibrational wavefunction, Hpi, and Hsigma-
* here for long range extrapolation
      if (r .gt. rpsi(npsi)) then
         hhpi=0.d0
         hh2sigm=0.d0
         hh4sigm=0.d0
         psivib=0d0
      else
* here for spline interpolation
         call dsplint (rrh,hpi(1),cshpi(1),nl,r,hhpi)
         call dsplint (rrh,h2sigm(1),csh2sigm(1),nl,r,hh2sigm)
         call dsplint (rrh,h4sigm(1),csh4sigm(1),nl,r,hh4sigm)
* now vibrational wavefunction
         call dsplint(rpsi,psi_v,cspsi,npsi,r,psivib)
      end if
      if (grndtst) then
        write (9,25) nvib, hhsigm, hhpi, psivib
25      format(' HSIGM         HPI           PSI (V=', i3,'):' 3g14.4)
      endif
      hh2sigm=hh2sigm*econv
      hh4sigm=hh4sigm*econv
      hhpi=hhpi*econv
* create source vector (eq. 229)
      if (grndtst) then 
         write (9,40) rmlmda
40       format('RMLMDA = ',d12.5)
         write (9,41)nch
41       format('SOURCE VECTOR, NCH = ',i3)
         write (9,*) 'RMLMDA, R: =', rmlmda, r
         write (9,411) rmlmda*psivib
411      format('rmlmda='g13.4,'psivib='g13.4)
      endif
      do i=1,nch
*   next line is eq. 242
*   for reference, the old hcl command is:    wf(6)=rmlmda*p*dipfunc
*   where p is the vibrational wavefunction
         term(i)=(const_pi(i)*hhpi+const_2sig(i)*hh2sigm+
     :       const_4sig(i)*hh4sigm)
         wf(i)=rmlmda*const(i)*term(i)*psivib
      enddo
      if (grndtst) then
         write (9,*) 'TERM'
         write (9, 42) term(i)
         write (9,*) 'CONST'
         write (9, 42) 
     :  ( term(i),i=1,nch)
42       format(11g14.5)
         write (9,*) 'WF VECTOR'
         write (9, 42) (wf(i),i=1,nch)
      endif
      entry spline_fit(ifirst)
* dummy initialization subroutine
      return
      entry wfintern(wf,yymin,nnvib,nny)
*  dummy ground subroutine
      return
      end

* --------------------------------------------------------------------------
      block data psiv
* --------------------------------------------------------------------------

*     revision date: 30-dec-2018
* --------------------------------------------------------------------------

*      data for OH A2sig+ state vibrational wavefunctions dawes calculation
      implicit double precision (a-h,o-z)
      common /corv_A/ rpsi(157)
      common /copsiv_A/ 
     + psiv_0(157), psiv_1(157), psiv_2(157), psiv_3(157),
     + psiv_4(157), psiv_5(157), psiv_6(157), psiv_7(157),
     + psiv_8(157)
* knots for spline interpolation of vibrational wavefunctions psi(v=0:7) in A state
      data rpsi /
     + 1.045d0, 1.075d0, 1.105d0, 1.135d0, 1.165d0, 1.195d0, 1.225d0,
     + 1.255d0, 1.285d0, 1.315d0, 1.345d0, 1.375d0, 1.405d0, 1.435d0,
     + 1.465d0, 1.495d0, 1.525d0, 1.555d0, 1.585d0, 1.615d0, 1.645d0,
     + 1.675d0, 1.705d0, 1.735d0, 1.765d0, 1.795d0, 1.825d0, 1.855d0,
     + 1.885d0, 1.915d0, 1.945d0, 1.975d0, 2.005d0, 2.035d0, 2.065d0,
     + 2.095d0, 2.125d0, 2.155d0, 2.185d0, 2.215d0, 2.245d0, 2.275d0,
     + 2.305d0, 2.335d0, 2.365d0, 2.395d0, 2.425d0, 2.455d0, 2.485d0,
     + 2.515d0, 2.545d0, 2.575d0, 2.605d0, 2.635d0, 2.665d0, 2.695d0,
     + 2.725d0, 2.755d0, 2.785d0, 2.815d0, 2.845d0, 2.875d0, 2.905d0,
     + 2.935d0, 2.965d0, 2.995d0, 3.025d0, 3.055d0, 3.085d0, 3.115d0,
     + 3.145d0, 3.175d0, 3.205d0, 3.235d0, 3.265d0, 3.295d0, 3.325d0,
     + 3.355d0, 3.385d0, 3.415d0, 3.445d0, 3.475d0, 3.505d0, 3.535d0,
     + 3.565d0, 3.595d0, 3.625d0, 3.655d0, 3.685d0, 3.715d0, 3.745d0,
     + 3.775d0, 3.805d0, 3.835d0, 3.865d0, 3.895d0, 3.925d0, 3.955d0,
     + 3.985d0, 4.015d0, 4.045d0, 4.075d0, 4.105d0, 4.135d0, 4.165d0,
     + 4.195d0, 4.225d0, 4.255d0, 4.285d0, 4.315d0, 4.345d0, 4.375d0,
     + 4.405d0, 4.435d0, 4.465d0, 4.495d0, 4.525d0, 4.555d0, 4.585d0,
     + 4.615d0, 4.645d0, 4.675d0, 4.705d0, 4.735d0, 4.765d0, 4.795d0,
     + 4.825d0, 4.855d0, 4.885d0, 4.915d0, 4.945d0, 4.975d0, 5.005d0,
     + 5.035d0, 5.065d0, 5.095d0, 5.125d0, 5.155d0, 5.185d0, 5.215d0,
     + 5.245d0, 5.275d0, 5.305d0, 5.335d0, 5.365d0, 5.395d0, 5.425d0,
     + 5.455d0, 5.485d0, 5.515d0, 5.545d0, 5.575d0, 5.605d0, 5.635d0,
     + 5.665d0, 5.695d0, 5.725d0/
* data for spline interpolation of Dawes vibrational wavefunctions psi(v=0:8) in A state
* from numerov calculations an a grid of r=[1:0.005:20], splined to the rpsi grid
      data psiv_0/
     + 3.8059884d-7, 1.3511044d-6, 4.4395603d-6, 1.3546473d-5,
     + 3.8530188d-5, 1.0259563d-4, 2.5653849d-4, 6.0374613d-4,
     + 1.3407989d-3, 2.8188900d-3, 5.6272499d-3, 1.0685961d-2,
     + 1.9337251d-2, 3.3421085d-2, 5.5315263d-2, 8.7830654d-2,
     + 1.3395681d-1, 1.9654199d-1, 2.7798862d-1, 3.7974638d-1,
     + 5.0153031d-1, 6.4099308d-1, 7.9397443d-1, 9.5490211d-1,
     + 1.1161466d0, 1.2685342d0, 1.4031556d0, 1.5135232d0, 1.5952879d0,
     + 1.6449777d0, 1.6602124d0, 1.6400219d0, 1.5857088d0, 1.5017422d0,
     + 1.3946334d0, 1.2713782d0, 1.1383697d0, 1.0016101d0,
     + 8.6661373d-1, 7.3789423d-1, 6.1868824d-1, 5.1105199d-1,
     + 4.1610788d-1, 3.3415281d-1, 2.6479989d-1, 2.0716594d-1,
     + 1.6008112d-1, 1.2223309d-1, 9.2271626d-2, 6.8890454d-2,
     + 5.0889605d-2, 3.7209430d-2, 2.6940626d-2, 1.9322303d-2,
     + 1.3732701d-2, 9.6749765d-3, 6.7592009d-3, 4.6842783d-3,
     + 3.2212755d-3, 2.1987840d-3, 1.4901930d-3, 1.0030978d-3,
     + 6.7082678d-4, 4.4582030d-4, 2.9451709d-4, 1.9345650d-4,
     + 1.2638463d-4, 8.2138616d-5, 5.3118497d-5, 3.4189740d-5,
     + 2.1907986d-5, 1.3978507d-5, 8.8829971d-6, 5.6232848d-6,
     + 3.5468651d-6, 2.2295104d-6, 1.3968866d-6, 8.7252156d-7,
     + 5.4341627d-7, 3.3752505d-7, 2.0910421d-7, 1.2923068d-7,
     + 7.9685332d-8, 4.9030234d-8, 3.0107832d-8, 1.8453370d-8,
     + 1.1290239d-8, 6.8962364d-9, 4.2058132d-9, 2.5612842d-9,
     + 1.5576671d-9, 9.4610561d-10, 5.7397127d-10, 3.4782394d-10,
     + 2.1056080d-10, 1.2734265d-10, 7.6944447d-11, 4.6453029d-11,
     + 2.8022518d-11, 1.6891909d-11, 1.0175374d-11, 6.1255108d-12,
     + 3.6852963d-12, 2.2159354d-12, 1.3317183d-12, 7.9993291d-13,
     + 4.8027994d-13, 2.8823540d-13, 1.7291156d-13, 1.0368979d-13,
     + 6.2157539d-14, 3.7248330d-14, 2.2314324d-14, 1.3363895d-14,
     + 8.0013395d-15, 4.7893856d-15, 2.8660970d-15, 1.7147555d-15,
     + 1.0256986d-15, 6.1340747d-16, 3.6677104d-16, 2.1926149d-16,
     + 1.3105559d-16, 7.8321102d-17, 4.6798925d-17, 2.7959562d-17,
     + 1.6701890d-17, 9.9757336d-18, 5.9575953d-18, 3.5575154d-18,
     + 2.1240991d-18, 1.2681114d-18, 7.5700199d-19, 4.5185145d-19,
     + 2.6968412d-19, 1.6094519d-19, 9.6042923d-20, 5.7308510d-20,
     + 3.4193283d-20, 2.0400081d-20, 1.2170090d-20, 7.2598532d-21,
     + 4.3304721d-21, 2.5829570d-21, 1.5405463d-21, 9.1877472d-22,
     + 5.4792485d-22, 3.2674697d-22, 1.9484159d-22, 1.1618019d-22,
     + 6.9272946d-23, 4.1302569d-23, 2.4624811d-23, 1.4680859d-23,
     + 8.7521027d-24, 5.2173901d-24, 3.1100428d-24/
      data psiv_1/
     + 1.3650547d-6,  4.7674506d-6,  1.5395725d-5,  4.6116777d-5, 
     + 1.2861036d-4,  3.3533368d-4,  8.1990154d-4,  1.8838809d-3, 
     + 4.0777007d-3,  8.3402116d-3,  1.6164655d-2,  2.9736445d-2, 
     + 5.1999150d-2,  8.6604621d-2,  1.3770227d-1,  2.0932166d-1, 
     + 3.0443327d-1,  4.2400276d-1,  5.6630035d-1,  7.2607162d-1, 
     + 8.9358601d-1,  1.0551026d0,  1.1946767d0,  1.2962456d0, 
     + 1.3441303d0,  1.3253013d0,  1.2326579d0,  1.0672744d0, 
     + 8.3719185d-1,  5.5535762d-1,  2.3873033d-1, -9.2955067d-2,
     + -4.1889766d-1, -7.1966213d-1, -9.7945315d-1, -1.1871387d0,
     + -1.3362253d0, -1.4251557d0, -1.4570707d0, -1.4386087d0,
     + -1.3784593d0, -1.2863067d0, -1.1721209d0, -1.0452774d0,
     + -9.1388368d-1, -7.8444394d-1, -6.6187065d-1, -5.4954417d-1,
     + -4.4944477d-1, -3.6237546d-1, -2.8825252d-1, -2.2637429d-1,
     + -1.7563501d-1, -1.3470633d-1, -1.0218607d-1, -7.6709287d-2,
     + -5.7013463d-2, -4.1974771d-2, -3.0624411d-2, -2.2151012d-2,
     + -1.5890719d-2, -1.1310730d-2, -7.9908674d-3, -5.6053134d-3,
     + -3.9053128d-3, -2.7033891d-3, -1.8599329d-3, -1.2721764d-3,
     + -8.6533009d-4, -5.8549841d-4, -3.9418587d-4, -2.6413029d-4,
     + -1.7618966d-4, -1.1702883d-4, -7.7421302d-5, -5.1024763d-5,
     + -3.3507470d-5, -2.1929526d-5, -1.4306438d-5, -9.3053502d-6,
     + -6.0354171d-6, -3.9041415d-6, -2.5191765d-6, -1.6217208d-6,
     + -1.0416964d-6, -6.6774431d-7, -4.2720785d-7, -2.7282468d-7,
     + -1.7393835d-7, -1.1071835d-7, -7.0371951d-8, -4.4666108d-8,
     + -2.8313668d-8, -1.7926273d-8, -1.1336856d-8, -7.1620343d-9,
     + -4.5201478d-9, -2.8501651d-9, -1.7956171d-9, -1.1303369d-9,
     + -7.1101061d-10, -4.4693082d-10, -2.8075015d-10, -1.7625161d-10,
     + -1.1058521d-10, -6.9347130d-11, -4.3465271d-11, -2.7230283d-11,
     + -1.7051817d-11, -1.0673586d-11, -6.6785576d-12, -4.1773205d-12,
     + -2.6119524d-12, -1.6326542d-12, -1.0202187d-12, -6.3733859d-13,
     + -3.9804492d-13, -2.4853395d-13, -1.5514494d-13, -9.6826380d-14,
     + -6.0417030d-14, -3.7691178d-14, -2.3509295d-14, -1.4661000d-14,
     + -9.1414651d-15, -5.6990186d-15, -3.5523854d-15, -2.2140085d-15,
     + -1.3796880d-15, -8.5966222d-16, -5.3557824d-16, -3.3363285d-16,
     + -2.0781073d-16, -1.2942635d-16, -8.0600019d-17, -5.0188864d-17,
     + -3.1249373d-17, -1.9455338d-17, -1.2111597d-17, -7.5392974d-18,
     + -4.6927634d-18, -2.9207622d-18, -1.8177528d-18, -1.1312166d-18,
     + -7.0393137d-19, -4.3801558d-19, -2.7253644d-19, -1.6956506d-19,
     + -1.0549355d-19, -6.5628743d-20, -4.0826461d-20, -2.5396246d-20,
     + -1.5797119d-20, -9.8257724d-21, -6.1113103d-21, -3.8008038d-21,
     + -2.3636063d-21/
      data psiv_2/
     + -3.4465028d-6, -1.1851888d-5, -3.7647335d-5, -1.1080161d-4,
     + -3.0324541d-4, -7.7493192d-4, -1.8544077d-3, -4.1637023d-3,
     + -8.7917645d-3, -1.7508367d-2, -3.2970992d-2, -5.8793893d-2,
     + -9.9393665d-2, -1.5955117d-1, -2.4366134d-1, -3.5432827d-1,
     + -4.9065709d-1, -6.4696974d-1, -8.1243784d-1, -9.7104300d-1,
     + -1.1019702d0, -1.1824598d0, -1.1922644d0, -1.1177531d0,
     + -9.5357089d-1, -7.0513214d-1, -3.9014607d-1, -3.7044430d-2, 
     + 3.1966436d-1,  6.4409228d-1,  9.0302977d-1,  1.0698366d0, 
     + 1.1280901d0,  1.0742152d0,  9.1700660d-1,  6.7499245d-1, 
     + 3.7307690d-1,  3.9520221d-2, -2.9746318d-1, -6.1283516d-1,
     + -8.8675663d-1, -1.1056416d0, -1.2625132d0, -1.3564550d0,
     + -1.3914266d0, -1.3748497d0, -1.3163531d0, -1.2265027d0,
     + -1.1156713d0, -9.9321640d-1, -8.6705217d-1, -7.4346088d-1,
     + -6.2703966d-1, -5.2080571d-1, -4.2642603d-1, -3.4450758d-1,
     + -2.7485711d-1, -2.1671975d-1, -1.6899121d-1, -1.3039825d-1,
     + -9.9627090d-2, -7.5408819d-2, -5.6574935d-2, -4.2090335d-2,
     + -3.1066450d-2, -2.2758353d-2, -1.6554133d-2, -1.1960444d-2,
     + -8.5865223d-3, -6.1273060d-3, -4.3476175d-3, -3.0682847d-3,
     + -2.1544023d-3, -1.5054614d-3, -1.0472418d-3, -7.2539184d-4,
     + -5.0043627d-4, -3.4393487d-4, -2.3553514d-4, -1.6076198d-4,
     + -1.0938160d-4, -7.4202578d-5, -5.0198202d-5, -3.3871172d-5,
     + -2.2799023d-5, -1.5311225d-5, -1.0260638d-5, -6.8623280d-6,
     + -4.5809863d-6, -3.0527245d-6, -2.0309798d-6, -1.3491512d-6,
     + -8.9494838d-7, -5.9286801d-7, -3.9226300d-7, -2.5923395d-7,
     + -1.7113397d-7, -1.1286050d-7, -7.4359595d-8, -4.8949480d-8,
     + -3.2195939d-8, -2.1160226d-8, -1.3897193d-8, -9.1209754d-9,
     + -5.9824915d-9, -3.9216404d-9, -2.5692996d-9, -1.6824332d-9,
     + -1.1011619d-9, -7.2039127d-10, -4.7108813d-10, -3.0793817d-10,
     + -2.0121623d-10, -1.3143489d-10, -8.5825421d-11, -5.6025703d-11,
     + -3.6562222d-11, -2.3853891d-11, -1.5558738d-11, -1.0145758d-11,
     + -6.6144775d-12, -4.3113480d-12, -2.8095868d-12, -1.8305797d-12,
     + -1.1924939d-12, -7.7669268d-13, -5.0579187d-13, -3.2932739d-13,
     + -2.1439801d-13, -1.3955767d-13, -9.0830111d-14, -5.9108788d-14,
     + -3.8461221d-14, -2.5023346d-14, -1.6278761d-14, -1.0588958d-14,
     + -6.8872070d-15, -4.4791238d-15, -2.9127609d-15, -1.8940010d-15,
     + -1.2314616d-15, -8.0062360d-16, -5.2048018d-16, -3.3833715d-16,
     + -2.1992075d-16, -1.4294040d-16, -9.2900316d-17, -6.0374553d-17,
     + -3.9234328d-17, -2.5495001d-17, -1.6566133d-17, -1.0763786d-17,
     + -6.9933735d-18, -4.5434406d-18, -2.9515840d-18, -1.9172835d-18,
     + -1.2452329d-18/
      data psiv_3/
     + -7.2525103d-6, -2.4578086d-5, -7.6863502d-5, -2.2248018d-4,
     + -5.9812111d-4, -1.4995280d-3, -3.5155040d-3, -7.7211199d-3,
     + -1.5919926d-2, -3.0898149d-2, -5.6584769d-2, -9.7883813d-2,
     + -1.6007032d-1, -2.4773074d-1, -3.6332716d-1, -5.0504066d-1,
     + -6.6470198d-1, -8.2704084d-1, -9.7085576d-1, -1.0711290d0,
     + -1.1020303d0, -1.0426618d0, -8.8338097d-1, -6.2989405d-1,
     + -3.0340451d-1,  6.0518761d-2,  4.1623344d-1,  7.1601600d-1, 
     + 9.1857767d-1,  9.9537163d-1,  9.3508504d-1,  7.4566075d-1, 
     + 4.5359480d-1,  9.9912631d-2, -2.6720816d-1, -6.0022802d-1,
     + -8.5863571d-1, -1.0138954d0, -1.0522978d0, -9.7498712d-1,
     + -7.9577560d-1, -5.3782336d-1, -2.2975248d-1,  9.8503247d-2, 
     + 4.1923457d-1,  7.0947110d-1,  9.5248228d-1,  1.1382805d0, 
     + 1.2632170d0,  1.3289513d0,  1.3411928d0,  1.3083545d0, 
     + 1.2402005d0,  1.1466916d0,  1.0371522d0,  9.1976492d-1, 
     + 8.0125029d-1,  6.8674705d-1,  5.7986950d-1,  4.8290360d-1, 
     + 3.9703320d-1,  3.2256839d-1,  2.5917312d-1,  2.0608101d-1, 
     + 1.6227549d-1,  1.2662033d-1,  9.7956796d-2,  7.5173543d-2, 
     + 5.7253506d-2,  4.3295850d-2,  3.2522694d-2,  2.4276875d-2, 
     + 1.8014527d-2,  1.3293279d-2,  9.7582256d-3,  7.1281575d-3, 
     + 5.1829344d-3,  3.7522147d-3,  2.7054086d-3,  1.9432302d-3, 
     + 1.3907926d-3,  9.9206936d-4,  7.0543543d-4,  5.0014813d-4, 
     + 3.5362918d-4,  2.4939026d-4,  1.7545485d-4,  1.2316192d-4, 
     + 8.6273697d-5,  6.0315521d-5,  4.2090453d-5,  2.9322159d-5, 
     + 2.0394696d-5,  1.4164231d-5,  9.8234655d-6,  6.8041353d-6, 
     + 4.7071332d-6,  3.2527506d-6,  2.2453662d-6,  1.5484450d-6, 
     + 1.0668555d-6,  7.3441861d-7,  5.0516581d-7,  3.4721514d-7, 
     + 2.3848401d-7,  1.6369525d-7,  1.1229165d-7,  7.6985676d-8, 
     + 5.2752046d-8,  3.6128569d-8,  2.4731944d-8,  1.6922878d-8, 
     + 1.1574717d-8,  7.9136693d-9,  5.4086298d-9,  3.6952846d-9, 
     + 2.5238787d-9,  1.7232860d-9,  1.1763117d-9,  8.0273255d-10, 
     + 5.4765783d-10,  3.7354586d-10,  2.5473033d-10,  1.7367009d-10, 
     + 1.1838101d-10,  8.0678239d-11,  5.4973398d-11,  3.7451981d-11, 
     + 2.5510966d-11,  1.7374508d-11,  1.1831369d-11,  8.0555950d-12, 
     + 5.4840751d-12,  3.7329751d-12,  2.5407114d-12,  1.7290463d-12, 
     + 1.1765525d-12,  8.0051891d-13,  5.4461489d-13,  3.7048192d-13, 
     + 2.5200315d-13,  1.7139892d-13,  1.1656684d-13,  7.9269890d-14, 
     + 5.3902549d-14,  3.6650469d-14,  2.4918404d-14,  1.6940747d-14, 
     + 1.1516423d-14,  7.8284584d-15,  5.3211935d-15,  3.6167287d-15, 
     + 2.4580771d-15,  1.6704889d-15,  1.1351418d-15,  7.7124448d-16, 
     + 5.2386412d-16/
      data psiv_4/
     + -1.3289212d-5, -4.4430761d-5, -1.3695501d-4, -3.9032639d-4,
     + -1.0320893d-3, -2.5418027d-3, -5.8458555d-3, -1.2576412d-2,
     + -2.5356557d-2, -4.8030168d-2, -8.5656803d-2, -1.4393117d-1,
     + -2.2795022d-1, -3.4044227d-1, -4.7975553d-1, -6.3736255d-1,
     + -7.9623909d-1, -9.3175283d-1, -1.0154410d0, -1.0200039d0,
     + -9.2492864d-1, -7.2378827d-1, -4.2971138d-1, -7.6006498d-2,
     + 2.8851075d-1, 6.0709986d-1, 8.2608205d-1, 9.0708396d-1,
     + 8.3519331d-1, 6.2156596d-1, 3.0195431d-1, -6.9268689d-2,
     + -4.2876867d-1, -7.1623839d-1, -8.8557918d-1, -9.1207196d-1,
     + -7.9483050d-1, -5.5562569d-1, -2.3376335d-1, 1.2183464d-1,
     + 4.6102240d-1, 7.3978745d-1, 9.2578321d-1, 1.0012980d0,
     + 9.6354684d-1, 8.2284358d-1, 5.9946407d-1, 3.1968838d-1,
     + 1.1786316d-2, -2.9740903d-1, -5.8502823d-1, -8.3373264d-1,
     + -1.0323144d0, -1.1754650d0, -1.2629942d0, -1.2987730d0,
     + -1.2894939d0, -1.2434465d0, -1.1694660d0, -1.0761661d0,
     + -9.7137482d-1, -8.6175795d-1, -7.5264405d-1, -6.4803448d-1,
     + -5.5072104d-1, -4.6243280d-1, -3.8401625d-1, -3.1563356d-1,
     + -2.5696080d-1, -2.0734574d-1, -1.6593461d-1, -1.3177449d-1,
     + -1.0389548d-1, -8.1365779d-2, -6.3323354d-2, -4.8993937d-2,
     + -3.7699713d-2, -2.8860631d-2, -2.1988477d-2, -1.6678058d-2,
     + -1.2597402d-2, -9.4780336d-3, -7.1051162d-3, -5.3082121d-3,
     + -3.9531874d-3, -2.9353318d-3, -2.1735209d-3, -1.6052787d-3,
     + -1.1827525d-3, -8.6948512d-4, -6.3785198d-4, -4.6701398d-4,
     + -3.4131170d-4, -2.4902127d-4, -1.8139859d-4, -1.3194420d-4,
     + -9.5840744d-5, -6.9527096d-5, -5.0377682d-5, -3.6461704d-5,
     + -2.6362316d-5, -1.9041834d-5, -1.3741684d-5, -9.9083767d-6,
     + -7.1387293d-6, -5.1394685d-6, -3.6975597d-6, -2.6584636d-6,
     + -1.9102174d-6, -1.3717930d-6, -9.8460914d-7, -7.0635460d-7,
     + -5.0649899d-7, -3.6303120d-7, -2.6009422d-7, -1.8627283d-7,
     + -1.3335529d-7, -9.5438215d-8, -6.8280170d-8, -4.8835459d-8,
     + -3.4918194d-8, -2.4960376d-8, -1.7837748d-8, -1.2744555d-8,
     + -9.1035494d-9, -6.5013503d-9, -4.6420323d-9, -3.3138228d-9,
     + -2.3652187d-9, -1.6878672d-9, -1.2042982d-9, -8.5913667d-10,
     + -6.1281066d-10, -4.3704820d-10, -3.1165514d-10, -2.2221021d-10,
     + -1.5841677d-10, -1.1292450d-10, -8.0487301d-11, -5.7361560d-11,
     + -4.0876248d-11, -2.9125911d-11, -2.0751441d-11, -1.4783559d-11,
     + -1.0531090d-11, -7.5012355d-12, -5.3426751d-12, -3.8049802d-12,
     + -2.7096591d-12, -1.9295055d-12, -1.3738716d-12, -9.7816731d-13,
     + -6.9637327d-13, -4.9570380d-13, -3.5280221d-13, -2.5102915d-13,
     + -1.7852982d-13/
      data psiv_5/
     + -2.1486110d-5, -7.0953574d-5, -2.1583631d-4, -6.0647816d-4,
     + -1.5793783d-3, -3.8263778d-3, -8.6459413d-3, -1.8247665d-2,
     + -3.6033454d-2, -6.6722115d-2, -1.1606743d-1, -1.8975116d-1,
     + -2.9148293d-1, -4.2064857d-1, -5.7009798d-1, -7.2398410d-1,
     + -8.5747919d-1, -9.4009903d-1, -9.4237137d-1, -8.4324001d-1,
     + -6.3687474d-1, -3.3888138d-1,  1.2270799d-2,  3.6073652d-1, 
     + 6.4405790d-1,  8.0644953d-1,  8.1288786d-1,  6.5959019d-1, 
     + 3.7566973d-1,  1.6933528d-2, -3.4515671d-1, -6.3803881d-1,
     + -8.0338058d-1, -8.0976090d-1, -6.5842931d-1, -3.8121207d-1,
     + -3.2428501d-2,  3.2257062d-1,  6.2045425d-1,  8.1145628d-1, 
     + 8.6694864d-1,  7.8242625d-1,  5.7607344d-1,  2.8344013d-1,
     + -5.0062773d-2, -3.7709042d-1, -6.5538487d-1, -8.5297934d-1,
     + -9.5101243d-1, -9.4418273d-1, -8.3935584d-1, -6.5285215d-1,
     + -4.0700389d-1, -1.2666849d-1,  1.6380669d-1,  4.4291814d-1, 
     + 6.9354259d-1,  9.0364678d-1,  1.0662960d0,  1.1791934d0, 
     + 1.2438815d0,  1.2647501d0,  1.2480256d0,  1.2008940d0, 
     + 1.1307923d0,  1.0448211d0,  9.4933413d-1,  8.4971079d-1, 
     + 7.5029872d-1,  6.5441836d-1,  5.6442621d-1,  4.8183043d-1, 
     + 4.0744902d-1,  3.4156196d-1,  2.8403915d-1,  2.3445630d-1, 
     + 1.9220007d-1,  1.5655635d-1,  1.2677033d-1,  1.0208957d-1, 
     + 8.1795231d-2,  6.5224608d-2,  5.1782048d-2,  4.0942025d-2, 
     + 3.2248233d-2,  2.5310535d-2,  1.9799915d-2,  1.5441679d-2, 
     + 1.2008507d-2,  9.3138786d-3,  7.2060655d-3,  5.5624720d-3, 
     + 4.2846011d-3,  3.2937393d-3,  2.5273326d-3,  1.9359125d-3, 
     + 1.4805124d-3,  1.1305490d-3,  8.6210607d-4,  6.5654861d-4, 
     + 4.9939814d-4,  3.7943457d-4,  2.8798444d-4,  2.1836006d-4, 
     + 1.6541622d-4,  1.2520179d-4,  9.4687982d-5,  7.1557160d-5, 
     + 5.4038873d-5,  4.0782548d-5,  3.0759242d-5,  2.3186039d-5, 
     + 1.7467978d-5,  1.3153415d-5,  9.8998333d-6,  7.4477163d-6, 
     + 5.6006148d-6,  4.2099460d-6,  3.1634129d-6,  2.3762011d-6, 
     + 1.7842968d-6,  1.3394167d-6,  1.0051636d-6,  7.5411480d-7, 
     + 5.6561946d-7,  4.2413444d-7,  3.1796615d-7,  2.3832075d-7, 
     + 1.7858769d-7,  1.3379951d-7,  1.0022483d-7,  7.5061658d-8, 
     + 5.6206546d-8,  4.2080916d-8,  3.1500427d-8,  2.3576749d-8, 
     + 1.7643742d-8,  1.3201995d-8,  9.8771864d-9,  7.3888058d-9, 
     + 5.5266883d-9,  4.1334020d-9,  3.0910379d-9,  2.3113034d-9, 
     + 1.7280938d-9,  1.2919240d-9,  9.6575535d-10,  7.2186935d-10, 
     + 5.3952491d-10,  4.0320387d-10,  3.0129709d-10,  2.2512059d-10, 
     + 1.6817908d-10,  1.2561449d-10,  9.3793311d-11,  6.9997683d-11, 
     + 5.2194152d-11/
      data psiv_6/
     + 3.1079947d-5, 1.0150516d-4, 3.0513103d-4, 8.4653356d-4,
     + 2.1745003d-3, 5.1908447d-3, 1.1542988d-2, 2.3942838d-2,
     + 4.6393057d-2, 8.4140507d-2, 1.4305854d-1, 2.2801245d-1,
     + 3.4041695d-1, 4.7560509d-1, 6.2089266d-1, 7.5437810d-1,
     + 8.4654020d-1, 8.6608214d-1, 7.8881923d-1, 6.0607645d-1,
     + 3.3062174d-1, -6.4703778d-4, -3.3109631d-1, -5.9613261d-1,
     + -7.3838935d-1, -7.2292558d-1, -5.4890428d-1, -2.5274751d-1,
     + 1.0000633d-1, 4.3015474d-1, 6.6296105d-1, 7.4553880d-1,
     + 6.5956397d-1, 4.2574859d-1, 9.7735806d-2, -2.5172223d-1,
     + -5.4769827d-1, -7.2935027d-1, -7.6217907d-1, -6.4339444d-1,
     + -4.0000546d-1, -8.1224739d-2, 2.5267002d-1, 5.4252497d-1,
     + 7.4077037d-1, 8.1840996d-1, 7.6793251d-1, 6.0214178d-1,
     + 3.4976468d-1, 4.9216519d-2, -2.5808666d-1, -5.3379675d-1,
     + -7.4725402d-1, -8.7821104d-1, -9.1769180d-1, -8.6729654d-1,
     + -7.3730248d-1, -5.4409016d-1, -3.0740626d-1, -4.7854154d-2,
     + 2.1519583d-1, 4.6512173d-1, 6.8893689d-1, 8.7763217d-1,
     + 1.0260924d0, 1.1326664d0, 1.1985235d0, 1.2269390d0, 1.2226324d0,
     + 1.1911305d0, 1.1382099d0, 1.0694660d0, 9.9004697d-1,
     + 9.0448476d-1, 8.1658387d-1, 7.2938512d-1, 6.4520376d-1,
     + 5.6571600d-1, 4.9203796d-1, 4.2481041d-1, 3.6429052d-1,
     + 3.1044976d-1, 2.6305245d-1, 2.2171694d-1, 1.8596740d-1,
     + 1.5527951d-1, 1.2911570d-1, 1.0694728d-1, 8.8269699d-2,
     + 7.2613337d-2, 5.9550878d-2, 4.8699606d-2, 3.9720916d-2,
     + 3.2318478d-2, 2.6235830d-2, 2.1253072d-2, 1.7182953d-2,
     + 1.3867021d-2, 1.1172022d-2, 8.9865837d-3, 7.2180493d-3,
     + 5.7896438d-3, 4.6379963d-3, 3.7110124d-3, 2.9660123d-3,
     + 2.3681256d-3, 1.8889366d-3, 1.5053513d-3, 1.1986484d-3,
     + 9.5368132d-4, 7.5821861d-4, 6.0240042d-4, 4.7829351d-4,
     + 3.7952449d-4, 3.0098002d-4, 2.3856296d-4, 1.8899473d-4,
     + 1.4965486d-4, 1.1845100d-4, 9.3714057d-5, 7.4113823d-5,
     + 5.8591140d-5, 4.6303324d-5, 3.6580402d-5, 2.8890080d-5,
     + 2.2809739d-5, 1.8004048d-5, 1.4207081d-5, 1.1208059d-5,
     + 8.8400037d-6, 6.9706980d-6, 5.4954961d-6, 4.3316069d-6,
     + 3.4135556d-6, 2.6895822d-6, 2.1187829d-6, 1.6688431d-6,
     + 1.3142420d-6, 1.0348303d-6, 8.1470394d-7, 6.4131321d-7,
     + 5.0475745d-7, 3.9722788d-7, 3.1256696d-7, 2.4592026d-7,
     + 1.9346144d-7, 1.5217504d-7, 1.1968511d-7, 9.4119842d-8,
     + 7.4004840d-8, 5.8178852d-8, 4.5727357d-8, 3.5930172d-8,
     + 2.8220136d-8, 2.2150534d-8, 1.7369410d-8, 1.3599304d-8/
      data psiv_7/
     + 3.1079947d-5, 1.0150516d-4, 3.0513103d-4, 8.4653356d-4,
     + 2.1745003d-3, 5.1908447d-3, 1.1542988d-2, 2.3942838d-2,
     + 4.6393057d-2, 8.4140507d-2, 1.4305854d-1, 2.2801245d-1,
     + 3.4041695d-1, 4.7560509d-1, 6.2089266d-1, 7.5437810d-1,
     + 8.4654020d-1, 8.6608214d-1, 7.8881923d-1, 6.0607645d-1,
     + 3.3062174d-1, -6.4703778d-4, -3.3109631d-1, -5.9613261d-1,
     + -7.3838935d-1, -7.2292558d-1, -5.4890428d-1, -2.5274751d-1,
     + 1.0000633d-1, 4.3015474d-1, 6.6296105d-1, 7.4553880d-1,
     + 6.5956397d-1, 4.2574859d-1, 9.7735806d-2, -2.5172223d-1,
     + -5.4769827d-1, -7.2935027d-1, -7.6217907d-1, -6.4339444d-1,
     + -4.0000546d-1, -8.1224739d-2, 2.5267002d-1, 5.4252497d-1,
     + 7.4077037d-1, 8.1840996d-1, 7.6793251d-1, 6.0214178d-1,
     + 3.4976468d-1, 4.9216519d-2, -2.5808666d-1, -5.3379675d-1,
     + -7.4725402d-1, -8.7821104d-1, -9.1769180d-1, -8.6729654d-1,
     + -7.3730248d-1, -5.4409016d-1, -3.0740626d-1, -4.7854154d-2,
     + 2.1519583d-1, 4.6512173d-1, 6.8893689d-1, 8.7763217d-1,
     + 1.0260924d0, 1.1326664d0, 1.1985235d0, 1.2269390d0, 1.2226324d0,
     + 1.1911305d0, 1.1382099d0, 1.0694660d0, 9.9004697d-1,
     + 9.0448476d-1, 8.1658387d-1, 7.2938512d-1, 6.4520376d-1,
     + 5.6571600d-1, 4.9203796d-1, 4.2481041d-1, 3.6429052d-1,
     + 3.1044976d-1, 2.6305245d-1, 2.2171694d-1, 1.8596740d-1,
     + 1.5527951d-1, 1.2911570d-1, 1.0694728d-1, 8.8269699d-2,
     + 7.2613337d-2, 5.9550878d-2, 4.8699606d-2, 3.9720916d-2,
     + 3.2318478d-2, 2.6235830d-2, 2.1253072d-2, 1.7182953d-2,
     + 1.3867021d-2, 1.1172022d-2, 8.9865837d-3, 7.2180493d-3,
     + 5.7896438d-3, 4.6379963d-3, 3.7110124d-3, 2.9660123d-3,
     + 2.3681256d-3, 1.8889366d-3, 1.5053513d-3, 1.1986484d-3,
     + 9.5368132d-4, 7.5821861d-4, 6.0240042d-4, 4.7829351d-4,
     + 3.7952449d-4, 3.0098002d-4, 2.3856296d-4, 1.8899473d-4,
     + 1.4965486d-4, 1.1845100d-4, 9.3714057d-5, 7.4113823d-5,
     + 5.8591140d-5, 4.6303324d-5, 3.6580402d-5, 2.8890080d-5,
     + 2.2809739d-5, 1.8004048d-5, 1.4207081d-5, 1.1208059d-5,
     + 8.8400037d-6, 6.9706980d-6, 5.4954961d-6, 4.3316069d-6,
     + 3.4135556d-6, 2.6895822d-6, 2.1187829d-6, 1.6688431d-6,
     + 1.3142420d-6, 1.0348303d-6, 8.1470394d-7, 6.4131321d-7,
     + 5.0475745d-7, 3.9722788d-7, 3.1256696d-7, 2.4592026d-7,
     + 1.9346144d-7, 1.5217504d-7, 1.1968511d-7, 9.4119842d-8,
     + 7.4004840d-8, 5.8178852d-8, 4.5727357d-8, 3.5930172d-8,
     + 2.8220136d-8, 2.2150534d-8, 1.7369410d-8, 1.3599304d-8/
      data psiv_8/
     + 4.5158683d-5, 1.4499591d-4, 4.2797679d-4, 1.1642381d-3,
     + 2.9278450d-3, 6.8307556d-3, 1.4816630d-2, 2.9911607d-2,
     + 5.6262959d-2, 9.8754576d-2, 1.6191269d-1, 2.4775952d-1,
     + 3.5316535d-1, 4.6769778d-1, 5.7309405d-1, 6.4440928d-1,
     + 6.5454036d-1, 5.8241488d-1, 4.2213979d-1, 1.8894892d-1,
     + -8.0118549d-2, -3.3156281d-1, -5.0797249d-1, -5.6445544d-1,
     + -4.8264130d-1, -2.7887291d-1, -3.0852742d-3, 2.7341432d-1,
     + 4.7775464d-1, 5.5552387d-1, 4.8569745d-1, 2.8664590d-1,
     + 1.1421069d-2, -2.6694697d-1, -4.7596984d-1, -5.6293054d-1,
     + -5.0753737d-1, -3.2583450d-1, -6.4490751d-2, 2.1249222d-1,
     + 4.4018724d-1, 5.6791724d-1, 5.6986592d-1, 4.4896530d-1,
     + 2.3394981d-1, -2.8805142d-2, -2.8619362d-1, -4.8946842d-1,
     + -6.0313206d-1, -6.1015615d-1, -5.1315647d-1, -3.3190912d-1,
     + -9.8199563d-2, 1.5059718d-1, 3.7787961d-1, 5.5312231d-1,
     + 6.5546138d-1, 6.7530176d-1, 6.1409570d-1, 4.8271604d-1,
     + 2.9887987d-1, 8.4153756d-2, -1.3895732d-1, -3.4948907d-1,
     + -5.2988075d-1, -6.6712756d-1, -7.5323346d-1, -7.8508663d-1,
     + -7.6393004d-1, -6.9452713d-1, -5.8417725d-1, -4.4172798d-1,
     + -2.7669501d-1, -9.8511338d-2, 8.4068617d-2, 2.6337742d-1,
     + 4.3307135d-1, 5.8823552d-1, 7.2537253d-1, 8.4229624d-1,
     + 9.3796822d-1, 1.0123287d0, 1.0661079d0, 1.1006242d0,
     + 1.1175938d0, 1.1189847d0, 1.1069013d0, 1.0834722d0, 1.0507578d0,
     + 1.0106895d0, 9.6504522d-1, 9.1542401d-1, 8.6322837d-1,
     + 8.0965856d-1, 7.5572594d-1, 7.0226841d-1, 6.4996132d-1,
     + 5.9933102d-1, 5.5077344d-1, 5.0457542d-1, 4.6093000d-1,
     + 4.1995007d-1, 3.8168189d-1, 3.4611967d-1, 3.1321692d-1,
     + 2.8289502d-1, 2.5505094d-1, 2.2956459d-1, 2.0630515d-1,
     + 1.8513537d-1, 1.6591502d-1, 1.4850386d-1, 1.3276435d-1,
     + 1.1856330d-1, 1.0577304d-1, 9.4272193d-2, 8.3946395d-2,
     + 7.4688640d-2, 6.6399319d-2, 5.8986135d-2, 5.2363960d-2,
     + 4.6454632d-2, 4.1186620d-2, 3.6494652d-2, 3.2319331d-2,
     + 2.8606760d-2, 2.5308155d-2, 2.2379440d-2, 1.9780869d-2,
     + 1.7476668d-2, 1.5434694d-2, 1.3626110d-2, 1.2025074d-2,
     + 1.0608462d-2, 9.3556078d-3, 8.2480602d-3, 7.2693630d-3,
     + 6.4048522d-3, 5.6414729d-3, 4.9676125d-3, 4.3729486d-3,
     + 3.8483119d-3, 3.3855621d-3, 2.9774772d-3, 2.6176529d-3,
     + 2.3004123d-3, 2.0207252d-3, 1.7741354d-3, 1.5566966d-3,
     + 1.3649137d-3, 1.1956911d-3, 1.0462866d-3, 9.1427000d-4,
     + 7.9748597d-4, 6.9402114d-4, 6.0217450d-4, 5.2043106d-4/
      end
