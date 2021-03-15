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
      common /corv_A/ rpsi(170)
      common /copsiv_A/ 
     + psiv_0(170), psiv_1(170), psiv_2(170), psiv_3(170),
     + psiv_4(170), psiv_5(170), psiv_6(170), psiv_7(170)
      dimension psi_v(170)
      dimension c2sig(5), c4sig(5), cpi(5)
*  9 states of a given parity
      dimension wf(9), sc(9), const_pi(9), const_2sig(9), 
     +   const_4sig(9), term(9),
     +   term1(9), const(9)
      dimension rrh(33), hpi(33), h2sigm(33), h4sigm(33), 
     +   csh2sigm(33), csh4sigm(33), cshpi(33)
      dimension cspsi(170)
      common /coered/ ered, rmu
      econv=1.d0/219474.6d0
      data nl, npsi, npsi_mat / 33, 170, 2210/
*  knots for spline expansion of Hsigm and Hpi spin-orbit coupling terms
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
         if ((nvib.lt.0).or.(nvib.gt.7)) then
            write (6,15) nvib
            write (9,15) nvib
15            format(/'nvib =',i3,' OUT OF RANGE, ABORT ***')
            call exit
         endif
***  write out wf array (for testing purposes)
         if (grndtst1) then
            write(9,*) 'GROUND VIBRATIONAL WAVEFUNCTION'
            do irr = 1, 170
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

*     revision date: 19-oct-2018
* --------------------------------------------------------------------------

*      data for HS A2sig+ state vibrational wavefunctions dawes calculation
      implicit double precision (a-h,o-z)
      common /corv_A/ rpsi(170)
      common /copsiv_A/ 
     + psiv_0(170), psiv_1(170), psiv_2(170), psiv_3(170),
     + psiv_4(170), psiv_5(170), psiv_6(170), psiv_7(170)
* knots for spline interpolation of vibrational wavefunctions psi(v=0:7) in A state
      data rpsi /
     + 1.8141d0, 1.8897d0, 1.9653d0, 2.0409d0, 2.1165d0, 2.1921d0,
     + 2.2677d0, 2.3433d0, 2.4188d0, 2.4944d0, 2.5700d0, 2.6456d0,
     + 2.7212d0, 2.7968d0, 2.8724d0, 2.9480d0, 3.0236d0, 3.0992d0,
     + 3.1747d0, 3.2503d0, 3.3259d0, 3.4015d0, 3.4771d0, 3.5527d0,
     + 3.6283d0, 3.7039d0, 3.7795d0, 3.8550d0, 3.9306d0, 4.0062d0,
     + 4.0818d0, 4.1574d0, 4.2330d0, 4.3086d0, 4.3842d0, 4.4598d0,
     + 4.5353d0, 4.6109d0, 4.6865d0, 4.7621d0, 4.8377d0, 4.9133d0,
     + 4.9889d0, 5.0645d0, 5.1401d0, 5.2156d0, 5.2912d0, 5.3668d0,
     + 5.4424d0, 5.5180d0, 5.5936d0, 5.6692d0, 5.7448d0, 5.8204d0,
     + 5.8959d0, 5.9715d0, 6.0471d0, 6.1227d0, 6.1983d0, 6.2739d0,
     + 6.3495d0, 6.4251d0, 6.5007d0, 6.5762d0, 6.6518d0, 6.7274d0,
     + 6.8030d0, 6.8786d0, 6.9542d0, 7.0298d0, 7.1054d0, 7.1810d0,
     + 7.2565d0, 7.3321d0, 7.4077d0, 7.4833d0, 7.5589d0, 7.6345d0,
     + 7.7101d0, 7.7857d0, 7.8613d0, 7.9369d0, 8.0124d0, 8.0880d0,
     + 8.1636d0, 8.2392d0, 8.3148d0, 8.3904d0, 8.4660d0, 8.5416d0,
     + 8.6172d0, 8.6927d0, 8.7683d0, 8.8439d0, 8.9195d0, 8.9951d0,
     + 9.0707d0, 9.1463d0, 9.2219d0, 9.2975d0, 9.3730d0, 9.4486d0,
     + 9.5242d0, 9.5998d0, 9.6754d0, 9.7510d0, 9.8266d0, 9.9022d0,
     + 9.9778d0, 1.0053d+1, 1.0129d+1, 1.0205d+1, 1.0280d+1, 1.0356d+1,
     + 1.0431d+1, 1.0507d+1, 1.0582d+1, 1.0658d+1, 1.0734d+1,
     + 1.0809d+1, 1.0885d+1, 1.0960d+1, 1.1036d+1, 1.1112d+1,
     + 1.1187d+1, 1.1263d+1, 1.1338d+1, 1.1414d+1, 1.1490d+1,
     + 1.1565d+1, 1.1641d+1, 1.1716d+1, 1.1792d+1, 1.1867d+1,
     + 1.1943d+1, 1.2019d+1, 1.2094d+1, 1.2170d+1, 1.2245d+1,
     + 1.2321d+1, 1.2397d+1, 1.2472d+1, 1.2548d+1, 1.2623d+1,
     + 1.2699d+1, 1.2775d+1, 1.2850d+1, 1.2926d+1, 1.3001d+1,
     + 1.3077d+1, 1.3152d+1, 1.3228d+1, 1.3304d+1, 1.3379d+1,
     + 1.3455d+1, 1.3530d+1, 1.3606d+1, 1.3682d+1, 1.3757d+1,
     + 1.3833d+1, 1.3908d+1, 1.3984d+1, 1.4060d+1, 1.4135d+1,
     + 1.4211d+1, 1.4286d+1, 1.4362d+1, 1.4438d+1, 1.4513d+1, 1.4589d+1/
* data for spline interpolation of CBS vibrational wavefunctions psi(v=0:7) in A state
* from numerov calculations an a grid of r=[1:0.005:20], splined to the rpsi grid
      data psiv_0/
     + 1.2460d-4,7.1119d-4,3.2401d-3,1.1981d-2,3.6722d-2,9.4640d-2,
     + 2.0815d-1,3.9486d-1,6.5239d-1,9.4817d-1,1.2199d0,1.4123d0,
     + 1.4869d0,1.4306d0,1.2684d0,1.0450d0,8.0426d-1,5.8248d-1,
     + 3.9934d-1,2.6063d-1,1.6278d-1,9.7735d-2,5.6668d-2,3.1852d-2,
     + 1.7421d-2,9.3009d-3,4.8617d-3,2.4945d-3,1.2594d-3,6.2688d-4,
     + 3.0821d-4,1.4992d-4,7.2240d-5,3.4528d-5,1.6386d-5,7.7279d-6,
     + 3.6248d-6,1.6922d-6,7.8666d-7,3.6436d-7,1.6821d-7,7.7435d-8,
     + 3.5556d-8,1.6289d-8,7.4474d-9,3.3988d-9,1.5486d-9,7.0453d-10,
     + 3.2010d-10,1.4526d-10,6.5849d-11,2.9820d-11,1.3491d-11,
     + 6.0989d-12,2.7549d-12,1.2435d-12,5.6095d-13,2.5289d-13,
     + 1.1394d-13,5.1313d-14,2.3097d-14,1.0392d-14,4.6736d-15,
     + 2.1011d-15,9.4425d-16,4.2421d-16,1.9052d-16,8.5544d-17,
     + 3.8398d-17,1.7232d-17,7.7311d-18,3.4679d-18,1.5553d-18,
     + 6.9737d-19,3.1264d-19,1.4014d-19,6.2806d-20,2.8144d-20,
     + 1.2610d-20,5.6492d-21,2.5305d-21,1.1334d-21,5.0761d-22,
     + 2.2731d-22,1.0178d-22,4.5572d-23,2.0403d-23,9.1337d-24,
     + 4.0886d-24,1.8301d-24,8.1913d-25,3.6661d-25,1.6407d-25,
     + 7.3425d-26,3.2857d-26,1.4703d-26,6.5789d-27,2.9437d-27,
     + 1.3171d-27,5.8928d-28,2.6364d-28,1.1795d-28,5.2767d-29,
     + 2.3606d-29,1.0561d-29,4.7243d-30,2.1134d-30,9.4540d-31,
     + 4.2290d-31,1.8917d-31,8.4617d-32,3.7849d-32,1.6929d-32,
     + 7.5721d-33,3.3867d-33,1.5148d-33,6.7747d-34,3.0299d-34,
     + 1.3551d-34,6.0601d-35,2.7101d-35,1.2120d-35,5.4199d-36,
     + 2.4237d-36,1.0838d-36,4.8463d-37,2.1671d-37,9.6898d-38,
     + 4.3327d-38,1.9372d-38,8.6617d-39,3.8727d-39,1.7315d-39,
     + 7.7412d-40,3.4609d-40,1.5473d-40,6.9172d-41,3.0923d-41,
     + 1.3824d-41,6.1797d-42,2.7625d-42,1.2349d-42,5.5198d-43,
     + 2.4673d-43,1.1029d-43,4.9295d-44,2.2033d-44,9.8477d-45,
     + 4.4014d-45,1.9672d-45,8.7917d-46,3.9292d-46,1.7560d-46,
     + 7.8475d-47,3.5070d-47,1.5672d-47,7.0034d-48,3.1295d-48,
     + 1.3984d-48,6.2489d-49,2.7922d-49,1.2476d-49,5.5746d-50,
     + 2.4908d-50,1.1129d-50,4.9721d-51,2.2214d-51,9.9248d-52,
     + 4.4340d-52,1.9809d-52/
      data psiv_1/
     + 3.7305d-4,2.0312d-3,8.7681d-3,3.0465d-2,8.6851d-2,2.0556d-1,
     + 4.0846d-1,6.8490d-1,9.6952d-1,1.1518d0,1.1177d0,8.2945d-1,
     + 3.3764d-1,-2.3938d-1,-7.6174d-1,-1.1249d0,-1.2869d0,-1.2684d0,
     + -1.1241d0,-9.1740d-1,-6.9996d-1,-5.0481d-1,-3.4708d-1,-2.2906d-1,
     + -1.4594d-1,-9.0194d-2,-5.4298d-2,-3.1952d-2,-1.8437d-2,
     + -1.0458d-2,-5.8458d-3,-3.2261d-3,-1.7607d-3,-9.5174d-4,
     + -5.1017d-4,-2.7149d-4,-1.4355d-4,-7.5483d-5,-3.9498d-5,
     + -2.0579d-5,-1.0682d-5,-5.5262d-6,-2.8504d-6,-1.4664d-6,
     + -7.5259d-7,-3.8544d-7,-1.9703d-7,-1.0055d-7,-5.1232d-8,
     + -2.6068d-8,-1.3247d-8,-6.7237d-9,-3.4092d-9,-1.7269d-9,
     + -8.7400d-10,-4.4197d-10,-2.2333d-10,-1.1277d-10,-5.6905d-11,
     + -2.8699d-11,-1.4465d-11,-7.2874d-12,-3.6695d-12,-1.8470d-12,
     + -9.2924d-13,-4.6734d-13,-2.3495d-13,-1.1808d-13,-5.9329d-14,
     + -2.9800d-14,-1.4964d-14,-7.5123d-15,-3.7705d-15,-1.8920d-15,
     + -9.4924d-16,-4.7615d-16,-2.3880d-16,-1.1974d-16,-6.0035d-17,
     + -3.0095d-17,-1.5085d-17,-7.5599d-18,-3.7883d-18,-1.8982d-18,
     + -9.5098d-19,-4.7640d-19,-2.3864d-19,-1.1953d-19,-5.9862d-20,
     + -2.9979d-20,-1.5012d-20,-7.5171d-21,-3.7638d-21,-1.8844d-21,
     + -9.4343d-22,-4.7230d-22,-2.3643d-22,-1.1835d-22,-5.9242d-23,
     + -2.9653d-23,-1.4842d-23,-7.4284d-24,-3.7178d-24,-1.8607d-24,
     + -9.3125d-25,-4.6605d-25,-2.3324d-25,-1.1672d-25,-5.8410d-26,
     + -2.9229d-26,-1.4626d-26,-7.3188d-27,-3.6622d-27,-1.8324d-27,
     + -9.1685d-28,-4.5873d-28,-2.2952d-28,-1.1483d-28,-5.7450d-29,
     + -2.8741d-29,-1.4379d-29,-7.1931d-30,-3.5984d-30,-1.8001d-30,
     + -9.0043d-31,-4.5041d-31,-2.2530d-31,-1.1269d-31,-5.6366d-32,
     + -2.8192d-32,-1.4100d-32,-7.0522d-33,-3.5270d-33,-1.7639d-33,
     + -8.8214d-34,-4.4115d-34,-2.2061d-34,-1.1032d-34,-5.5166d-35,
     + -2.7585d-35,-1.3793d-35,-6.8970d-36,-3.4485d-36,-1.7243d-36,
     + -8.6209d-37,-4.3102d-37,-2.1549d-37,-1.0773d-37,-5.3860d-38,
     + -2.6926d-38,-1.3460d-38,-6.7288d-39,-3.3636d-39,-1.6814d-39,
     + -8.4046d-40,-4.2010d-40,-2.0998d-40,-1.0495d-40,-5.2457d-41,
     + -2.6218d-41,-1.3103d-41,-6.5488d-42,-3.2729d-42,-1.6356d-42,
     + -8.1738d-43,-4.0847d-43,-2.0412d-43,-1.0200d-43,-5.0968d-44,
     + -2.5468d-44/
      data psiv_2/
     + -7.6644d-4,-4.0004d-3,-1.6447d-2,-5.3985d-2,-1.4392d-1,
     + -3.1431d-1,-5.6576d-1,-8.3638d-1,-9.9793d-1,-9.1510d-1,
     + -5.3475d-1,4.0228d-2,5.9840d-1,9.2332d-1,9.0092d-1,5.6051d-1,
     + 3.7286d-2,-5.0031d-1,-9.2061d-1,-1.1607d0,-1.2219d0,-1.1460d0,
     + -9.8820d-1,-7.9837d-1,-6.1211d-1,-4.4949d-1,-3.1842d-1,
     + -2.1883d-1,-1.4658d-1,-9.6057d-2,-6.1782d-2,-3.9105d-2,
     + -2.4412d-2,-1.5059d-2,-9.1946d-3,-5.5640d-3,-3.3409d-3,
     + -1.9925d-3,-1.1814d-3,-6.9679d-4,-4.0911d-4,-2.3924d-4,
     + -1.3940d-4,-8.0973d-5,-4.6901d-5,-2.7097d-5,-1.5620d-5,
     + -8.9856d-6,-5.1595d-6,-2.9577d-6,-1.6929d-6,-9.6761d-7,
     + -5.5236d-7,-3.1495d-7,-1.7939d-7,-1.0208d-7,-5.8032d-8,
     + -3.2965d-8,-1.8711d-8,-1.0613d-8,-6.0154d-9,-3.4076d-9,
     + -1.9292d-9,-1.0917d-9,-6.1742d-10,-3.4904d-10,-1.9724d-10,
     + -1.1141d-10,-6.2910d-11,-3.5510d-11,-2.0038d-11,-1.1304d-11,
     + -6.3748d-12,-3.5942d-12,-2.0260d-12,-1.1418d-12,-6.4333d-13,
     + -3.6241d-13,-2.0412d-13,-1.1495d-13,-6.4721d-14,-3.6436d-14,
     + -2.0509d-14,-1.1543d-14,-6.4958d-15,-3.6551d-15,-2.0564d-15,
     + -1.1569d-15,-6.5078d-16,-3.6604d-16,-2.0587d-16,-1.1578d-16,
     + -6.5107d-17,-3.6610d-17,-2.0585d-17,-1.1574d-17,-6.5067d-18,
     + -3.6579d-18,-2.0563d-18,-1.1559d-18,-6.4973d-19,-3.6520d-19,
     + -2.0527d-19,-1.1537d-19,-6.4843d-20,-3.6444d-20,-2.0482d-20,
     + -1.1511d-20,-6.4686d-21,-3.6351d-21,-2.0427d-21,-1.1479d-21,
     + -6.4499d-22,-3.6241d-22,-2.0363d-22,-1.1441d-22,-6.4281d-23,
     + -3.6115d-23,-2.0289d-23,-1.1398d-23,-6.4033d-24,-3.5971d-24,
     + -2.0206d-24,-1.1350d-24,-6.3756d-25,-3.5811d-25,-2.0114d-25,
     + -1.1297d-25,-6.3451d-26,-3.5635d-26,-2.0013d-26,-1.1239d-26,
     + -6.3116d-27,-3.5444d-27,-1.9903d-27,-1.1176d-27,-6.2754d-28,
     + -3.5236d-28,-1.9784d-28,-1.1108d-28,-6.2365d-29,-3.5013d-29,
     + -1.9657d-29,-1.1035d-29,-6.1949d-30,-3.4776d-30,-1.9521d-30,
     + -1.0958d-30,-6.1507d-31,-3.4524d-31,-1.9377d-31,-1.0876d-31,
     + -6.1040d-32,-3.4258d-32,-1.9226d-32,-1.0789d-32,-6.0548d-33,
     + -3.3977d-33,-1.9066d-33,-1.0699d-33,-6.0032d-34,-3.3684d-34,
     + -1.8900d-34,-1.0604d-34,-5.9493d-35,-3.3378d-35,-1.8726d-35,
     + -1.0505d-35,-5.8932d-36,-3.3059d-36/
      data psiv_3/
     + -1.2356d-3,-6.2202d-3,-2.4521d-2,-7.6607d-2,-1.9251d-1,
     + -3.9106d-1,-6.4201d-1,-8.3816d-1,-8.2789d-1,-5.2156d-1,1.1187d-2,
     + 5.3624d-1,7.9241d-1,6.4905d-1,1.8583d-1,-3.6662d-1,-7.6207d-1,
     + -8.5543d-1,-6.4162d-1,-2.2241d-1,2.6109d-1,6.8739d-1,9.8477d-1,
     + 1.1333d0,1.1507d0,1.0729d0,9.3917d-1,7.8275d-1,6.2723d-1,
     + 4.8671d-1,3.6774d-1,2.7172d-1,1.9703d-1,1.4060d-1,9.8967d-2,
     + 6.8848d-2,4.7413d-2,3.2365d-2,2.1926d-2,1.4755d-2,9.8717d-3,
     + 6.5707d-3,4.3538d-3,2.8733d-3,1.8894d-3,1.2385d-3,8.0949d-4,
     + 5.2773d-4,3.4324d-4,2.2278d-4,1.4432d-4,9.3327d-5,6.0256d-5,
     + 3.8848d-5,2.5012d-5,1.6084d-5,1.0331d-5,6.6292d-6,4.2496d-6,
     + 2.7217d-6,1.7418d-6,1.1138d-6,7.1170d-7,4.5448d-7,2.9005d-7,
     + 1.8500d-7,1.1794d-7,7.5147d-8,4.7860d-8,3.0469d-8,1.9389d-8,
     + 1.2334d-8,7.8435d-9,4.9862d-9,3.1689d-9,2.0133d-9,1.2789d-9,
     + 8.1211d-10,5.1560d-10,3.2728d-10,2.0771d-10,1.3179d-10,
     + 8.3611d-11,5.3035d-11,3.3635d-11,2.1329d-11,1.3524d-11,
     + 8.5734d-12,5.4346d-12,3.4446d-12,2.1830d-12,1.3834d-12,
     + 8.7658d-13,5.5539d-13,3.5187d-13,2.2291d-13,1.4120d-13,
     + 8.9438d-14,5.6648d-14,3.5877d-14,2.2721d-14,1.4389d-14,
     + 9.1118d-15,5.7700d-15,3.6537d-15,2.3135d-15,1.4649d-15,
     + 9.2748d-16,5.8721d-16,3.7177d-16,2.3536d-16,1.4900d-16,
     + 9.4319d-17,5.9705d-17,3.7792d-17,2.3921d-17,1.5140d-17,
     + 9.5826d-18,6.0647d-18,3.8382d-18,2.4289d-18,1.5371d-18,
     + 9.7266d-19,6.1547d-19,3.8944d-19,2.4641d-19,1.5590d-19,
     + 9.8635d-20,6.2402d-20,3.9477d-20,2.4974d-20,1.5798d-20,
     + 9.9931d-21,6.3210d-21,3.9981d-21,2.5288d-21,1.5994d-21,
     + 1.0115d-21,6.3969d-22,4.0454d-22,2.5582d-22,1.6177d-22,
     + 1.0229d-22,6.4677d-23,4.0894d-23,2.5856d-23,1.6347d-23,
     + 1.0335d-23,6.5333d-24,4.1301d-24,2.6108d-24,1.6503d-24,
     + 1.0432d-24,6.5936d-25,4.1674d-25,2.6339d-25,1.6646d-25,
     + 1.0520d-25,6.6482d-26,4.2012d-26,2.6548d-26,1.6775d-26,
     + 1.0600d-26,6.6972d-27,4.2314d-27,2.6734d-27,1.6890d-27,
     + 1.0670d-27,6.7405d-28,4.2580d-28/
      data psiv_4/
     + -1.5991d-3,-7.8238d-3,-2.9834d-2,-8.9591d-2,-2.1461d-1,
     + -4.1062d-1,-6.2306d-1,-7.2604d-1,-5.8674d-1,-1.8837d-1,3.0493d-1,
     + 6.2293d-1,5.7535d-1,1.8769d-1,-3.1017d-1,-6.3831d-1,-6.3290d-1,
     + -3.1810d-1,1.4196d-1,5.4775d-1,7.5642d-1,7.2155d-1,4.8169d-1,
     + 1.2193d-1,-2.6555d-1,-6.0834d-1,-8.6417d-1,-1.0192d0,-1.0801d0,
     + -1.0652d0,-9.9669d-1,-8.9562d-1,-7.7937d-1,-6.6078d-1,-5.4831d-1,
     + -4.4688d-1,-3.5873d-1,-2.8427d-1,-2.2279d-1,-1.7294d-1,
     + -1.3314d-1,-1.0176d-1,-7.7289d-2,-5.8372d-2,-4.3868d-2,
     + -3.2822d-2,-2.4461d-2,-1.8166d-2,-1.3448d-2,-9.9266d-3,
     + -7.3084d-3,-5.3680d-3,-3.9343d-3,-2.8778d-3,-2.1012d-3,
     + -1.5317d-3,-1.1148d-3,-8.1026d-4,-5.8815d-4,-4.2641d-4,
     + -3.0881d-4,-2.2341d-4,-1.6147d-4,-1.1661d-4,-8.4137d-5,
     + -6.0663d-5,-4.3708d-5,-3.1471d-5,-2.2646d-5,-1.6286d-5,
     + -1.1707d-5,-8.4105d-6,-6.0398d-6,-4.3354d-6,-3.1108d-6,
     + -2.2313d-6,-1.5999d-6,-1.1468d-6,-8.2177d-7,-5.8870d-7,
     + -4.2162d-7,-3.0189d-7,-2.1611d-7,-1.5467d-7,-1.1068d-7,
     + -7.9183d-8,-5.6640d-8,-4.0509d-8,-2.8967d-8,-2.0711d-8,
     + -1.4806d-8,-1.0583d-8,-7.5640d-9,-5.4055d-9,-3.8625d-9,
     + -2.7598d-9,-1.9716d-9,-1.4085d-9,-1.0061d-9,-7.1860d-10,
     + -5.1323d-10,-3.6653d-10,-2.6176d-10,-1.8692d-10,-1.3348d-10,
     + -9.5310d-11,-6.8053d-11,-4.8588d-11,-3.4689d-11,-2.4764d-11,
     + -1.7678d-11,-1.2619d-11,-9.0075d-12,-6.4292d-12,-4.5886d-12,
     + -3.2748d-12,-2.3371d-12,-1.6678d-12,-1.1901d-12,-8.4916d-13,
     + -6.0587d-13,-4.3227d-13,-3.0839d-13,-2.2001d-13,-1.5694d-13,
     + -1.1195d-13,-7.9853d-14,-5.6955d-14,-4.0621d-14,-2.8970d-14,
     + -2.0660d-14,-1.4732d-14,-1.0505d-14,-7.4906d-15,-5.3408d-15,
     + -3.8078d-15,-2.7147d-15,-1.9353d-15,-1.3796d-15,-9.8337d-16,
     + -7.0093d-16,-4.9959d-16,-3.5606d-16,-2.5376d-16,-1.8084d-16,
     + -1.2887d-16,-9.1826d-17,-6.5429d-17,-4.6618d-17,-3.3214d-17,
     + -2.3662d-17,-1.6857d-17,-1.2008d-17,-8.5537d-18,-6.0927d-18,
     + -4.3395d-18,-3.0907d-18,-2.2011d-18,-1.5675d-18,-1.1163d-18,
     + -7.9486d-19,-5.6597d-19,-4.0298d-19,-2.8691d-19,-2.0426d-19,
     + -1.4541d-19,-1.0351d-19,-7.3685d-20,-5.2449d-20,-3.7331d-20/
      data psiv_5/
     + -1.6210d-3,-7.7779d-3,-2.8988d-2,-8.4701d-2,-1.9621d-1,
     + -3.5976d-1,-5.1522d-1,-5.4928d-1,-3.6831d-1,-4.8465d-3,3.5735d-1,
     + 4.9441d-1,3.1311d-1,-7.1449d-2,-4.1206d-1,-5.0164d-1,-2.9861d-1,
     + 6.9498d-2,4.0505d-1,5.5392d-1,4.6803d-1,2.0139d-1,-1.3659d-1,
     + -4.3484d-1,-6.1679d-1,-6.5248d-1,-5.5235d-1,-3.5180d-1,
     + -9.5516d-2,1.7448d-1,4.2545d-1,6.3618d-1,7.9613d-1,9.0326d-1,
     + 9.6125d-1,9.7717d-1,9.5951d-1,9.1685d-1,8.5705d-1,7.8677d-1,
     + 7.1140d-1,6.3505d-1,5.6070d-1,4.9037d-1,4.2533d-1,3.6625d-1,
     + 3.1337d-1,2.6660d-1,2.2566d-1,1.9015d-1,1.5957d-1,1.3341d-1,
     + 1.1117d-1,9.2357d-2,7.6514d-2,6.3229d-2,5.2129d-2,4.2886d-2,
     + 3.5213d-2,2.8860d-2,2.3615d-2,1.9293d-2,1.5739d-2,1.2824d-2,
     + 1.0435d-2,8.4823d-3,6.8876d-3,5.5872d-3,4.5283d-3,3.6670d-3,
     + 2.9672d-3,2.3993d-3,1.9387d-3,1.5656d-3,1.2635d-3,1.0192d-3,
     + 8.2173d-4,6.6220d-4,5.3340d-4,4.2948d-4,3.4567d-4,2.7812d-4,
     + 2.2370d-4,1.7987d-4,1.4458d-4,1.1619d-4,9.3349d-5,7.4980d-5,
     + 6.0212d-5,4.8343d-5,3.8806d-5,3.1144d-5,2.4991d-5,2.0051d-5,
     + 1.6084d-5,1.2901d-5,1.0346d-5,8.2963d-6,6.6518d-6,5.3327d-6,
     + 4.2748d-6,3.4265d-6,2.7464d-6,2.2011d-6,1.7640d-6,1.4136d-6,
     + 1.1327d-6,9.0752d-7,7.2706d-7,5.8245d-7,4.6656d-7,3.7370d-7,
     + 2.9930d-7,2.3969d-7,1.9194d-7,1.5369d-7,1.2306d-7,9.8518d-8,
     + 7.8868d-8,6.3132d-8,5.0531d-8,4.0443d-8,3.2366d-8,2.5900d-8,
     + 2.0725d-8,1.6582d-8,1.3266d-8,1.0613d-8,8.4895d-9,6.7905d-9,
     + 5.4310d-9,4.3434d-9,3.4734d-9,2.7774d-9,2.2207d-9,1.7754d-9,
     + 1.4193d-9,1.1346d-9,9.0690d-10,7.2485d-10,5.7930d-10,4.6294d-10,
     + 3.6992d-10,2.9557d-10,2.3615d-10,1.8866d-10,1.5071d-10,
     + 1.2038d-10,9.6151d-11,7.6792d-11,6.1326d-11,4.8971d-11,
     + 3.9102d-11,3.1220d-11,2.4925d-11,1.9897d-11,1.5883d-11,
     + 1.2677d-11,1.0118d-11,8.0749d-12,6.4438d-12,5.1418d-12,
     + 4.1026d-12,3.2731d-12,2.6112d-12,2.0830d-12,1.6615d-12,
     + 1.3252d-12,1.0569d-12,8.4283d-13/
      data psiv_6/
     + -1.2445d-3,-5.9089d-3,-2.1752d-2,-6.2625d-2,-1.4245d-1,
     + -2.5512d-1,-3.5364d-1,-3.5767d-1,-2.1085d-1,4.6173d-2,2.7183d-1,
     + 3.1817d-1,1.5005d-1,-1.2009d-1,-3.1118d-1,-3.0281d-1,-1.0803d-1,
     + 1.5073d-1,3.3008d-1,3.4611d-1,2.0384d-1,-2.4482d-2,-2.4517d-1,
     + -3.8400d-1,-4.0681d-1,-3.1921d-1,-1.5406d-1,4.4835d-2,2.3605d-1,
     + 3.8846d-1,4.8403d-1,5.1696d-1,4.9082d-1,4.1521d-1,3.0269d-1,
     + 1.6641d-1,1.8537d-2,-1.3065d-1,-2.7311d-1,-4.0307d-1,-5.1679d-1,
     + -6.1225d-1,-6.8879d-1,-7.4678d-1,-7.8733d-1,-8.1198d-1,
     + -8.2256d-1,-8.2103d-1,-8.0933d-1,-7.8930d-1,-7.6268d-1,
     + -7.3102d-1,-6.9571d-1,-6.5794d-1,-6.1872d-1,-5.7891d-1,
     + -5.3921d-1,-5.0017d-1,-4.6224d-1,-4.2574d-1,-3.9091d-1,
     + -3.5792d-1,-3.2687d-1,-2.9780d-1,-2.7073d-1,-2.4563d-1,
     + -2.2244d-1,-2.0110d-1,-1.8152d-1,-1.6361d-1,-1.4727d-1,
     + -1.3240d-1,-1.1889d-1,-1.0665d-1,-9.5569d-2,-8.5564d-2,
     + -7.6541d-2,-6.8415d-2,-6.1106d-2,-5.4541d-2,-4.8650d-2,
     + -4.3369d-2,-3.8640d-2,-3.4409d-2,-3.0627d-2,-2.7247d-2,
     + -2.4231d-2,-2.1540d-2,-1.9140d-2,-1.7002d-2,-1.5098d-2,
     + -1.3403d-2,-1.1895d-2,-1.0554d-2,-9.3622d-3,-8.3028d-3,
     + -7.3617d-3,-6.5260d-3,-5.7840d-3,-5.1255d-3,-4.5413d-3,
     + -4.0232d-3,-3.5637d-3,-3.1563d-3,-2.7952d-3,-2.4750d-3,
     + -2.1911d-3,-1.9396d-3,-1.7167d-3,-1.5192d-3,-1.3442d-3,
     + -1.1892d-3,-1.0520d-3,-9.3041d-4,-8.2279d-4,-7.2752d-4,
     + -6.4319d-4,-5.6855d-4,-5.0251d-4,-4.4408d-4,-3.9239d-4,
     + -3.4667d-4,-3.0623d-4,-2.7047d-4,-2.3886d-4,-2.1091d-4,
     + -1.8621d-4,-1.6437d-4,-1.4508d-4,-1.2804d-4,-1.1298d-4,
     + -9.9679d-5,-8.7932d-5,-7.7560d-5,-6.8401d-5,-6.0316d-5,
     + -5.3180d-5,-4.6881d-5,-4.1323d-5,-3.6419d-5,-3.2093d-5,
     + -2.8277d-5,-2.4911d-5,-2.1943d-5,-1.9327d-5,-1.7019d-5,
     + -1.4986d-5,-1.3193d-5,-1.1614d-5,-1.0222d-5,-8.9960d-6,
     + -7.9158d-6,-6.9645d-6,-6.1267d-6,-5.3890d-6,-4.7395d-6,
     + -4.1677d-6,-3.6644d-6,-3.2215d-6,-2.8318d-6,-2.4889d-6,
     + -2.1872d-6,-1.9218d-6,-1.6885d-6,-1.4832d-6,-1.3028d-6,
     + -1.1441d-6,-1.0047d-6,-8.8212d-7,-7.7441d-7/
      data psiv_7/
     + -7.1063d-4,-3.3605d-3,-1.2311d-2,-3.5240d-2,-7.9585d-2,
     + -1.4122d-1,-1.9322d-1,-1.9124d-1,-1.0638d-1,3.5161d-2,
     + 1.5301d-1,1.6816d-1,6.8017d-2,-7.9170d-2,-1.7318d-1,
     + -1.5407d-1,-3.8175d-2,1.0114d-1,1.8568d-1,1.7595d-1,
     + 8.3443d-2,-4.5582d-2,-1.5774d-1,-2.1506d-1,-2.0455d-1,
     + -1.3580d-1,-3.2160d-2,7.9252d-2,1.7508d-1,2.3969d-1,
     + 2.6602d-1,2.5439d-1,2.1029d-1,1.4218d-1,5.9477d-2,
     + -2.8778d-2,-1.1482d-1,-1.9255d-1,-2.5768d-1,-3.0754d-1,
     + -3.4090d-1,-3.5768d-1,-3.5864d-1,-3.4516d-1,-3.1903d-1,
     + -2.8220d-1,-2.3673d-1,-1.8460d-1,-1.2772d-1,-6.7808d-2,
     + -6.4128d-3,5.5127d-2,1.1568d-1,1.7429d-1,2.3022d-1,
     + 2.8289d-1,3.3186d-1,3.7683d-1,4.1763d-1,4.5418d-1,
     + 4.8648d-1,5.1460d-1,5.3867d-1,5.5885d-1,5.7534d-1,
     + 5.8836d-1,5.9814d-1,6.0492d-1,6.0895d-1,6.1047d-1,
     + 6.0971d-1,6.0692d-1,6.0230d-1,5.9607d-1,5.8841d-1,
     + 5.7952d-1,5.6957d-1,5.5870d-1,5.4707d-1,5.3480d-1,
     + 5.2201d-1,5.0882d-1,4.9533d-1,4.8161d-1,4.6776d-1,
     + 4.5383d-1,4.3989d-1,4.2600d-1,4.1221d-1,3.9855d-1,
     + 3.8506d-1,3.7178d-1,3.5872d-1,3.4591d-1,3.3338d-1,
     + 3.2112d-1,3.0917d-1,2.9751d-1,2.8617d-1,2.7514d-1,
     + 2.6443d-1,2.5404d-1,2.4397d-1,2.3422d-1,2.2477d-1,
     + 2.1562d-1,2.0677d-1,1.9820d-1,1.8992d-1,1.8191d-1,
     + 1.7418d-1,1.6672d-1,1.5952d-1,1.5258d-1,1.4588d-1,
     + 1.3943d-1,1.3322d-1,1.2724d-1,1.2149d-1,1.1596d-1,
     + 1.1064d-1,1.0553d-1,1.0062d-1,9.5909d-2,9.1387d-2,
     + 8.7049d-2,8.2890d-2,7.8904d-2,7.5086d-2,7.1429d-2,
     + 6.7928d-2,6.4578d-2,6.1374d-2,5.8311d-2,5.5383d-2,
     + 5.2586d-2,4.9915d-2,4.7364d-2,4.4931d-2,4.2609d-2,
     + 4.0395d-2,3.8285d-2,3.6274d-2,3.4358d-2,3.2534d-2,
     + 3.0798d-2,2.9146d-2,2.7574d-2,2.6080d-2,2.4659d-2,
     + 2.3309d-2,2.2027d-2,2.0810d-2,1.9654d-2,1.8557d-2,
     + 1.7516d-2,1.6529d-2,1.5594d-2,1.4707d-2,1.3867d-2,
     + 1.3072d-2,1.2318d-2,1.1605d-2,1.0930d-2,1.0292d-2,
     + 9.6886d-3,9.1180d-3,8.5788d-3,8.0693d-3,7.5881d-3/
      end
