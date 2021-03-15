*  System:  S(3P) + H(2S)
*
*   calculation of potential energy curves and spin-orbit matrix
*   elements, taken from m. alexander's notes (2018)
*
*   MRCI+Q calculation of richard dawes, oct 2018
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
      potnam='S(3P)-H(2S) DAWES 2018'
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
      potnam='S(3P)-H(2S) DAWES 2018'
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
* author:  paul dagdigian
* latest revision date:  18-sep-2018
* ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical ljunk, chlist
      include "common/parbas"
      common /covvl/ vvl(6)
      common /colpar/ ljunk(4), chlist
      parameter (npsi = 170, npsi_mat=13*npsi)
      dimension rr(29),v2pi(29),v4pi(29),v2sigm(29),v4sigm(29),
     :          cs2p(29), cs4p(29), cs2s(29), cs4s(29)
      dimension asr(5), bsr(5), clr(5), dV_sr(5), dv_lr(5)
      dimension vl(29), vec(29), v(6)

      dimension d1(4), d2(4), vlr(4)
      data npts / 29/
*
c  number of r values= 29
      data rr /
     + 1.5d0, 1.508508509d0, 1.704204204d0, 1.891391391d0,
     + 2.078578579d0, 2.265765766d0, 2.452952953d0, 2.529529530d0,
     + 2.648648649d0, 2.742242242d0, 2.835835836d0, 3.023023023d0,
     + 3.210210210d0, 3.397397397d0, 3.593093093d0, 3.780280280d0,
     + 3.967467467d0, 4.154654655d0, 4.537537538d0, 4.911911912d0,
     + 5.294794795d0, 5.669169169d0, 6.043543544d0, 6.426426426d0,
     + 6.800800801d0, 7.183683684d0, 7.940940941d0, 8.689689690d0, 
     + 10.d0 /
c  X2Pi
      data v2pi /
     + 66662.39d0, 64165.04d0, 18426.22d0, -6418.24d0, -20049.89d0,
     + -26812.76d0, -29453.55d0, -29824.84d0, -29455.35d0, -28827.72d0,
     + -27819.05d0, -25298.13d0, -22388.30d0, -19364.19d0, -16278.63d0,
     + -13520.20d0, -11025.49d0, -8835.33d0, -5348.96d0, -3122.51d0,
     + -1761.47d0, -1008.75d0, -587.36d0, -345.53d0, -213.03d0,
     + -135.17d0, -63.83d0, -37.81d0, -22.97d0/
c  4Sigma-
      data v4sigm /
     + 122056.47d0, 119521.60d0, 73093.33d0, 47786.80d0, 33571.22d0,
     + 24682.16d0, 17582.59d0, 15361.81d0, 12509.22d0, 10675.15d0,
     + 9122.26d0, 6693.10d0, 4933.54d0, 3653.09d0, 2681.68d0,
     + 2002.60d0, 1497.44d0, 1116.99d0, 595.62d0, 295.98d0, 119.49d0,
     + 26.29d0, -19.24d0, -37.61d0, -40.99d0, -37.91d0, -26.38d0,
     + -16.56d0, -7.02d0/
c  2Sigma-
      data v2sigm /
     + 125108.02d0, 122558.33d0, 75801.36d0, 50291.88d0, 35951.85d0,
     + 27721.00d0, 22747.10d0, 21200.02d0, 18604.38d0, 16723.82d0,
     + 15055.70d0, 12205.74d0, 9880.46d0, 7971.63d0, 6338.55d0,
     + 5062.18d0, 4017.68d0, 3166.36d0, 1896.19d0, 1108.88d0, 605.48d0,
     + 311.67d0, 141.88d0, 46.42d0, -1.99d0, -25.72d0, -37.86d0,
     + -34.83d0, -26.23d0/
c  4Pi
      data v4pi /
     + 157991.68d0, 155477.65d0, 109275.53d0, 83774.19d0, 66691.60d0,
     + 52272.82d0, 42313.86d0, 38856.07d0, 33994.47d0, 30566.40d0,
     + 27429.32d0, 21951.98d0, 17413.09d0, 13694.65d0, 10562.89d0,
     + 8175.89d0, 6282.62d0, 4793.65d0, 2695.51d0, 1487.99d0, 780.49d0,
     + 396.19d0, 187.84d0, 77.35d0, 24.16d0, -0.75d0, -13.31d0,
     + -11.49d0, -5.85d0/
* short and long range extrapolations for potential and derivatives
      data asr /
     + 55880312.70d0, 2671423.43d0, 4718057.40d0, 4936013.32d0,
     + 8382787.01d0/
      data bsr/
     + 4.4876d0, 1.8852d0, 2.4200d0, 2.4666d0, 2.9359d0/
      data clr/
     + 22971092.67d0, 5850140.10d0, 26225332.18d0, 7024329.37d0,
     + 6591685.00d0/
      data dv_sr/ 
     +  -293511.7d0, -295472.4d0, -299663.6d0, -297921.6d0, -297281.5d0/
      data dv_lr/
     +  11.325d0,4.3074d0,6.5674d0,7.2766d0, 7.1211d0/
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
     :         call dspline(rr,v4pi,npts,der1,der2,cs2p)
            if (i.eq.3)
     :         call dspline(rr,v2sigm,npts,der1,der2,cs2p)
            if (i.eq.4)
     :         call dspline(rr,v4sigm,npts,der1,der2,cs2p)
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
      a = 394.14
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
*  driven equation source vector for hs

*  current revision date: 24-nov-2018
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
      implicit double precision (a-h,o-z)
      logical grndtst, grndtst1, ljunk
      common /cosysi/ junk(3), nvib
      common /colpar/ ljunk(12), noprin
*  grndtst=.true.  if noprin=.false.
      common /cojtot/ jjtot, jjlpar
      common /coj12/  j12(9)
      common /coja/  jja(9)
      common /coel/  ll(9)
      dimension c2sig(5), c4sig(5), cpi(5)
*  9 states of a given parity
      dimension wf(9), sc(9), const_pi(9), const_2sig(9), 
     +   const_4sig(9), term(9),
     +   term1(9), const(9)
      dimension rrh(33), hpi(33), h2sigm(33), h4sigm(33), 
     +   csh2sigm(33), csh4sigm(33), cshpi(33), csA(33)
      common /coered/ ered, rmu
      econv=1.d0/219474.6d0
      data nl, npsi, npsi_mat / 33, 170, 2210/
*  knots for spline expansion of Hsigm and Hpi spin-orbit coupling terms and A2Sig PEC
      data rrh /
     +  0.9449d0, 1.1338d0, 1.3228d0, 1.5118d0, 1.7008d0,
     +  1.8897d0, 2.0787d0, 2.2677d0, 2.4566d0, 2.5322d0,
     +  2.6456d0, 2.7401d0, 2.8346d0, 3.0236d0, 3.2125d0,
     +  3.4015d0, 3.5905d0, 3.7795d0, 3.9684d0, 4.1574d0,
     +  4.5353d0, 4.9133d0, 5.2912d0, 5.6692d0, 6.0471d0,
     +  6.4251d0, 6.8030d0, 7.5589d0, 8.3148d0, 9.0707d0,
     + 10.3935d0, 11.3384d0, 12.2832d0/ 
c H2sigm -i <^2\Sigma^-(xy),1/2|Hso|A^2\Sigma^+,1/2> from av6z calculations
      data h2sigm /
     + 5.2459d0, 8.0486d0, 7.6619d0, 5.2722d0, 2.0133d0, -1.1732d0,
     + -4.6307d0, -8.0046d0, -7.7865d0, -4.7367d0, 0.4415d0, 5.1114d0,
     + 10.1196d0,  21.1390d0,  33.0981d0,  45.5772d0,  58.0521d0,
     + 69.9698d0,  80.8735d0,  90.4798d0, 105.5272d0, 115.5981d0,
     + 121.9763d0, 125.8701d0, 128.1871d0, 129.5435d0, 130.3252d0,
     + 131.0295d0, 131.2498d0, 131.3188d0, 131.3384d0, 131.3408d0,
     + 131.3415d0/
c H4sigm -i <^4\Sigma^-(xy),1/2|Hso|A^2\Sigma^+,1/2> from av6z calculations
      data h4sigm /
     + 10.6002d0, 14.8989d0, 14.5825d0, 11.2491d0, 6.5017d0, 1.7531d0,
     + -3.4712d0, -8.8849d0, -10.4581d0, -6.5316d0, 0.3785d0, 6.7404d0,
     + 13.6415d0, 29.0300d0, 45.8875d0, 63.5848d0, 81.3368d0,
     + 98.3267d0, 113.8835d0, 127.5904d0, 149.0497d0, 163.3914d0,
     + 172.4616d0, 177.9926d0, 181.2796d0, 183.2037d0, 184.3114d0,
     + 185.3010d0, 185.6134d0, 185.7116d0, 185.7477d0, 185.7424d0,
     + 185.7432d0/
c Hpi -\sqrt 2 < ^4\Pi_x,1/2|Hso| A,-1/2>
      data hpi /
     +  6.3271d0, 1.7398d0, -0.0438d0, 0.4088d0, 2.8213d0,
     +  7.5917d0,  16.7281d0,  33.8591d0,  64.8611d0,  69.4033d0, 
     +  74.6481d0,  78.1138d0,  80.8622d0,  84.4940d0,  85.9448d0, 
     +  85.4480d0,  83.3610d0,  80.1383d0,  76.2572d0,  72.1573d0, 
     +  64.4852d0,  58.4735d0,  54.2613d0,  51.5008d0,  49.7688d0, 
     +  48.7156d0,  48.0900d0,  47.5172d0,  47.3308d0,  47.2732d0, 
     +  47.2537d0,  47.2527d0,  47.2532d0/
c A2Sig+ Potential.  zero of energy is S(3P,J=2)+H
      data VA/
     + 266654.14d0, 221410.00d0, 160230.15d0, 99226.39d0, 54322.35d0,
     + 28349.92d0, 13824.40d0, 6129.27d0, 2384.20d0, 1676.20d0,
     + 1325.85d0, 1284.46d0, 1573.00d0, 2517.13d0, 3706.86d0,
     + 4874.85d0, 5893.17d0, 6716.23d0, 7348.11d0, 7818.08d0,
     + 8410.62d0, 8726.20d0, 8902.00d0, 9007.16d0, 9074.35d0,
     + 9119.41d0, 9150.39d0, 9187.20d0, 9205.58d0, 9215.03d0,
     + 9221.59d0, 9223.49d0, 9224.55d0/
      
* derivatives of VA for spline extrapolation
      dinf=(VA(nl)-VA(nl-1))/(rrh(nl)-rrh(nl-1))
      dzero=(VA(2)-VA(1))/(rrh(2)-rrh(1))
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
*  interpolate the coefficients for spline expansion of spin-orbit coupling matrix elts
*  and vibrational wavefunction in A state
*  assume zero derivatives at initial and final endpoints, except for A state PEC
        
         call dspline(rrh,h2sigm,nl,0d0,0d0,csh2sigm)
         call dspline(rrh,h4sigm,nl,0d0,0d0,csh4sigm)
         call dspline(rrh,hpi,nl,0d0,0d0,cshpi)
         call dspline(rrh,vA,nl,dzero,dinf,csA)
          
* sum over all channels to get the pi and sigma coefficients of Eq. 229
* THIS WILL HAVE TO BE MODIFIED SLIGHTLY SINCE YOU HAVE THE A STATE ALSO 
         do i=1,nch
             xll=ll(i)
             xj12=j12(i)+0.5
             xj=jjtot+0.5
                const1=sq2*sqrt(2*ll(i)+1d0)*((-1)**(jjtot+1))
                const3j= xf3j(xll,xj12,xj,zero,-half,half)
                const(i)=const1*const3j
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
         interp=1
      endif
* rmlmda is 2mu/hbar^2 in atomic units
      if (grndtst1) then
         xmuconv=1d0/5.4857990d-4
         rmlmda = two*xmuconv*(32*1)/33
      else
         rmlmda = two*rmu
      endif
* determine A state potential,  Hpi, and Hsigma-
* here for long range extrapolation
      if (r .gt. rr(33)) then
         hhpi=0.d0
         hh2sigm=0.d0
         hh4sigm=0.d0
* r-6 exptrapolation of A state potentiol at long range
         vvA=9226.26-5.88663d6*r**(-6)
      else
* here for spline interpolation
         call dsplint (rrh,hpi(1),cshpi(1),nl,r,hhpi)
         call dsplint (rrh,h2sigm(1),csh2sigm(1),nl,r,hh2sigm)
         call dsplint (rrh,h4sigm(1),csh4sigm(1),nl,r,hh4sigm)
         call dspline(rrh,vA(1),nl,dzero,dinf,csA,nl,r,vvA)
      end if
      hh2sigm=hh2sigm*econv
      hh4sigm=hh4sigm*econv
      hhpi=hhpi*econv
      vvA=vvA*econv
      end
