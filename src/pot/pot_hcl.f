*system: H(2S)+Cl(2P)
*reference:  M. H. Alexander, B. Pouilly, and T. Duhoo, J. Chem. Phys. 99, 175
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(14)
      include "common/parpot"
      potnam='ALEXANDER-POUILLY-DUHOO HCl(X,A,a,triplet)'
      print *, potnam
      print *, 'r, v1sig, v1pi, v3sig, v3pi'
      do ir=1,100
        r=1.5+(ir-1)*0.25d0
        if (r .gt. 20) goto 99
      call pot(vv0,r)
      write (6, 100) r,(vvl(i), i=1,4)
100   format(f8.3,5(1pe16.8))
      enddo
99    end

      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      common /coselb/ ibasty
      include "common/parpot"
      potnam='ALEXANDER-POUILLY-DUHOO HCl(X,A,a,triplet)'
      ibasty=10
      lammin(1)=1
      lammax(1)=14
      mproj(1)=0
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* --------------------------------------------------------------------------
      subroutine pot (vv0, rz)
*  subroutine to calculate the r-dependent coefficients in the
*  hcl potentials and matrix elements of l^2 and ly of alexander and pouilly.
*  in atomic units
* ----------------------------------------------------------------------
*  on entry:
*    rz:      interparticle distance
*  on return:
*  vv0        (not used)
*  variable in common block /conlam/ used here
*    nlam:    the number of angular coupling terms actually used
*  variable in common block /covvl/
*    vvl:     vector of length 14 to store r-dependence of each term
*             in potential expansion
*     vvl(1):  energy of 1 sigma state
*     vvl(2):  energy of 1 pi state
*     vvl(3):  energy of 3 sigma state
*     vvl(4):  energy of 3 pi state
*     vvl(5):  l^2 matrix element in 1 sigma state
*     vvl(6):  l^2 matrix element in 1 pi state
*     vvl(7):  l^2 matrix element in 3 sigma state
*     vvl(8):  l^2 matrix element in 3 pi state
*     vvl(9):  ly matrix element between 1 sigma state and 1 pi state
*     vvl(10):  ly matrix element between 3 sigma state and 3 pi state
*     vvl(11):  aso (for constant spin-orbit term)
*     vvl(12):  aso * ly matrix element between 1 sigma state and 1 pi state
*     vvl(13):  aso * ly matrix element between 3 sigma state and 3 pi state
*     vvl(14):  B off diagonal J.S matrix element
*  subroutines called:
*   bint4,bvalu  (spline routines from cmlib bspline package)

*   splint  (from numerical recipes)
*      to return a cubic-spline interpolated value for a given value rz
* author:  millard alexander
* latest revision date:  27-nov-1995 by mha
* ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      double precision l21si, l21pi, l23si, l23pi, ly1, ly3
      parameter (nl = 25, nn=27, n3sig=24)
      parameter (nnp6=nn+6, n3p6=n3sig+6, nw=5*nnp6)
      common /coered/ ered, rmu
      common /cosysr/ junk1, junk2, aso
      common /covvl/ vvl(14)
      common /cotest/ rc
      common /coconv/ econv, xmconv

      dimension rr(nn),e1sig(nn),e1pi(nn),e3pi(nn),e3sig(n3sig),
     :          r3sig(n3sig)
      dimension work(nw)
      dimension t1sig(nnp6),t1pi(nnp6),t3sig(n3p6),t3pi(nnp6),
     :   b1sig(nnp6),b1pi(nnp6),b3sig(n3p6),b3pi(nnp6)
      dimension rl(nl), l21si(nl), l21pi(nl),
     :  l23si(nl), l23pi(nl), ly1(nl), ly3(nl)
      dimension dl21si(nl), dl21pi(nl),
     :  dl23si(nl), dl23pi(nl), dly1(nl), dly3(nl)
      dimension clr(4,2)
      data rr /
     : 1.1d0, 1.3d0, 1.7d0, 1.9d0, 2.1d0, 2.3d0, 2.4d0, 2.5d0, 2.6d0,
     : 2.8d0, 3.d0, 3.5d0, 4.d0, 4.3d0, 4.8d0, 5.d0, 5.5d0, 6.d0, 6.5d0,
     : 7.d0, 8.d0, 8.5d0, 9.d0, 10.d0, 11.d0, 12.d0, 15.d0/
      data e1sig /
     : 1.3072525d0, 6.131322d-1, 1.865020d-2, -9.194390d-2,
     : -1.452875d-1, -1.655991d-1, -1.680491d-1, -1.670238d-1,
     : -1.634184d-1, -1.510889d-1, -1.349651d-1, -9.149480d-2,
     : -5.487870d-2, -3.831790d-2, -1.955770d-2, -1.465090d-2,
     : -6.9195d-3, -3.2336d-3, -1.5512d-3, -7.859d-4, -2.715d-4,
     : -1.913d-4, -1.499d-4,-1.1489d-4, -1.0109d-4, -9.4692d-5,
     : -8.8392d-5/
      data e1pi /
     : 1.6571828d0, 9.588772d-1, 3.580657d-1, 2.416018d-1,
     : 1.771567d-1, 1.371320d-1, 1.220484d-1, 1.086497d-1,
     : 9.642810d-2, 7.547560d-2, 5.877380d-2, 3.078770d-2,
     : 1.547170d-2, 9.9711d-3, 4.5333d-3, 3.2281d-3, 1.2566d-3,
     : 3.848d-4, 3.26d-5, -8.86d-5, -1.083d-4, -9.23d-5, -7.64d-5,
     : -5.3554d-5, -4.0954d-5, -3.4354d-5, -2.7554d-5/
      data e3pi /
     : 1.6512248d0, 9.537590d-1, 3.539502d-1, 2.361912d-1,
     : 1.679616d-1, 1.231962d-1, 1.062502d-1, 9.160190d-2,
     : 7.905790d-2, 5.904410d-2, 4.421870d-2, 2.161490d-2,
     : 1.046290d-2, 6.6549d-3, 3.0025d-3, 2.15d-3, 5.07d-4, 1.89d-5,
     : -1.424d-4, -1.708d-4, -1.255d-4, -9.99d-5, -7.97d-5,
     : -5.4055d-5, -4.1055d-5, -3.4355d-5, -2.7555d-5/
      data r3sig /
     : 1.75d0, 2.d0, 2.25d0, 2.4d0, 2.5d0, 2.6d0, 2.8d0, 3.d0, 3.5d0,
     : 4.d0, 4.3d0, 4.8d0, 5.d0, 5.5d0, 6.d0, 6.5d0, 7.d0, 8.d0, 8.5d0,
     : 9.d0, 10.d0, 11.d0, 12.d0, 15.d0/
      data e3sig/
     : 4.9292d-1, 3.5971d-1, 2.6668d-1, 2.1738d-1, 1.9267d-1, 1.7026d-1,
     : 1.3205d-1, 1.0161d-1, 5.1451d-2, 2.2331d-2, 1.3762d-2, 5.8248d-3,
     : 4.0326d-3, 1.464d-3, 3.972d-4, -9.1d-6, -1.418d-4, -1.631d-4,
     : -1.478d-4, -1.329d-4,-1.1189d-4, -1.0059d-4, -9.4594d-5,
     : -8.8394d-5/
      data clr /-2.2943d3, -6.9067d2, -2.2949d3, -6.9067d2,
     :      2.8966d5, 8.4665d4, 2.8978d5, 8.4665d4/
      data rl/ 1.8, 2.0, 2.2, 2.3, 2.4, 2.5, 2.6, 2.8, 3.0, 3.25,
     :  3.5, 3.75,
     : 4.0, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 8.25, 9.25,
     :  10.25,
     : 11.25, 30/
      data l21si / 2.1140321, 2.2347, 2.4371105,
     : 2.6236865, 2.8773171, 3.1601399, 3.455894,
     : 4.0947695, 4.8195224, 5.8714889, 7.1018307,
     : 8.5146526, 1.0098107d+1, 1.1825017d+1, 1.5565214d+1,
     : 1.9500822d+1, 2.3560678d+1, 2.7786238d+1, 3.2239602d+1,
     : 3.6966848d+1, 4.7342488d+1, 5.901643d+1, 7.2014782d+1,
     : 8.6343416d+1, 6.0176944d+2/
      data l21pi / 4.2338236, 4.5001031, 5.2762324, 6.0420134,
     : 6.8951097, 7.5936377, 8.1263882, 8.9423566,
     : 9.6547034, 1.0547238d+1, 1.1493884d+1, 1.2516124d+1,
     : 1.362624d+1, 1.4832768d+1, 1.7553607d+1, 2.0689396d+1,
     : 2.4225675d+1, 2.8143244d+1, 3.2426054d+1, 3.7062276d+1,
     : 4.7366337d+1, 5.9022062d+1, 7.2016064d+1, 8.6343709d+1,
     : 6.0176944d+2/
      data ly1 / 2.095848d-2, 5.388128d-2, 1.6594713d-1, 2.3644219d-1,
     : 3.0225764d-1, 3.5541769d-1, 3.9853086d-1, 4.6831081d-1,
     : 5.2699103d-1, 5.9275509d-1, 6.5385714d-1, 7.1168691d-1,
     : 7.6601837d-1, 8.1566616d-1, 8.954997d-1, 9.464112d-1,
     : 9.7420198d-1, 9.8797953d-1, 9.9446727d-1, 9.9744003d-1,
     : 9.9937759d-1, 9.99758d-1, 9.998302d-1, 9.998438d-1,
     : 9.9984671d-1/
      data l23si /3.4625324, 4.0096686, 5.736907,
     : 6.5460681, 7.1138628, 7.5426417, 7.9072724,
     : 8.583406, 9.2772285, 1.0213651d+1, 1.12315d+1,
     : 1.2328155d+1, 1.3504083d+1, 1.4762723d+1, 1.7545065d+1,
     : 2.0700701d+1, 2.423753d+1, 2.815058d+1, 3.24294d+1,
     : 3.7063273d+1, 4.7365921d+1, 5.9021729d+1, 7.2015914d+1,
     : 8.6343653d+1, 6.0176944d+2/
      data l23pi / 4.1900632, 4.7714766, 6.5123339,
     : 7.3124023, 7.8613708, 8.2627621, 8.5927349,
     : 9.1836094, 9.7813095, 1.0599331d+1, 1.1514184d+1,
     : 1.2529012d+1, 1.3643736d+1, 1.4858379d+1, 1.7588707d+1,
     : 2.0720197d+1, 2.4246155d+1, 2.8154381d+1, 3.2431069d+1,
     : 3.7064012d+1, 4.7366083d+1, 5.9021765d+1, 7.2015929d+1,
     : 8.6343664d+1, 6.0176944d+2/
      data ly3 /9.3631183d-1, 9.3155984d-1, 9.2871575d-1, 9.339911d-1,
     : 9.4121738d-1, 9.4840456d-1, 9.5487882d-1, 9.6543725d-1,
     : 9.7347887d-1, 9.8116926d-1, 9.8683204d-1, 9.909724d-1,
     : 9.938994d-1, 9.9590691d-1, 9.9813269d-1, 9.9908365d-1,
     : 9.9949166d-1, 9.9967344d-1, 9.9975929d-1, 9.9980177d-1,
     : 9.9983d-1, 9.9984037d-1, 9.9984363d-1, 9.9984495d-1,
     : 9.998467d-1/
      data interp / 1/
      data zero,half,one,two /0.d0,0.5d0,1.d0,2.d0/
      data inb1sig, inb3sig, inb1pi, inb3pi /1,1,1,1/
*
      if (interp .eq. 1) then
*
*  here for first pass through pot subroutine
*  define spline initial conditions
      ibcl = 1
      ibcr = 2
      fbcl = zero
      fbcr = zero
      kntopt = 2
      inbv=1
      dr=rr(2)-rr(1)
*  compute initial derivatives
      d1sig=(e1sig(2)-e1sig(1))/dr
      d1pi=(e1pi(2)-e1pi(1))/dr
      d3pi=(e3pi(2)-e3pi(1))/dr
      dr=r3sig(2)-r3sig(1)
      d3sig=(e3sig(2)-e3sig(1))/dr


*  interpolate the  coefficients for the energies
          call bint4 (rr, e1sig, nn, ibcl,ibcr,d1sig,fbcr,kntopt,
     :                t1sig, b1sig,
     :                nc1sig,kord, work)
          call bint4 (rr, e1pi, nn, ibcl,ibcr,d1pi,fbcr,kntopt,
     :                t1pi, b1pi,
     :                nc1pi,kord, work)
          call bint4 (r3sig, e3sig, n3sig, ibcl,ibcr,d3sig,fbcr,kntopt,
     :                t3sig, b3sig,
     :                nc3sig,kord, work)
          call bint4 (rr, e3pi, nn, ibcl,ibcr,d3pi,fbcr,kntopt,
     :                t3pi, b3pi,
     :                nc3pi,kord, work)

*  interpolate the  coefficients for the matrix elements

      nlm1 = nl -1
        dl2f1si =( l21si(nl) - l21si(nlm1)) / (rl(nl) - rl(nlm1))
        dl2i1si =( l21si(2) - l21si(1)) / (rl(2) - rl(1))
        call dspline( rl, l21si(1), nl, dl2i1si, dl2f1si, dl21si)
        dl2f1pi =( l21pi(nl) - l21pi(nlm1)) / (rl(nl) - rl(nlm1))
        dl2i1pi =( l21pi(2) - l21pi(1)) / (rl(2) - rl(1))
        call dspline( rl, l21pi(1), nl, dl2i1pi, dl2f1pi, dl21pi)
        dl2f3sio=( l23si(nl) - l23si(nlm1)) / (rl(nl) - rl(nlm1))
        dl2f3pi =( l23pi(nl) - l23pi(nlm1)) / (rl(nl) - rl(nlm1))
        dl2i3pi =( l23pi(2) - l23pi(1)) / (rl(2) - rl(1))
        call dspline( rl, l23pi(1), nl, dl2i3pi, dl2f3pi, dl23pi)
        dly1f =( ly1(nl) - ly1(nlm1)) / (rl(nl) - rl(nlm1))
        dly1i =( ly1(2) - ly1(1)) / (rl(2) - rl(1))
        call dspline( rl, ly1(1), nl, dly1i, dly1f, dly1)
        dly3f =( ly3(nl) - ly3(nlm1)) / (rl(nl) - rl(nlm1))
        dly3i =( ly3(2) - ly3(1)) / (rl(2) - rl(1))
        call dspline( rl, ly3(1), nl, dly3i, dly3f, dly3)

        interp = 0
      endif
* here is  the entry point for subsequent entries into pot subroutine
      call dscal(14,zero,vvl,1)
      vv0 = 0.d0
* here to determine potential
      kord=4
      izero=0
      if (rz .gt. rr(nn)) then
*  here for inverse power extrapolation at large r
        r8=rz**8
        r2=rz*rz
        do 50 ii=1,4
          vvl(ii)=(clr(ii,1)*r2+clr(ii,2))/r8
50      continue
      else
* here for spline range
        v1sig= bvalu (t1sig, b1sig, nc1sig, kord,
     :              izero, rz, inb1sig, work)
* exponential extrapolation of 3sigma at short r
        if (rz.gt.r3sig(1)) then
          v3sig= bvalu (t3sig, b3sig, nc3sig, kord,
     :              izero, rz, inb3sig, work)
        else
          v3sig=exp(-1.2601949512382d0*rz+ 1.49793277475464d0)
        endif
        v1pi= bvalu (t1pi, b1pi, nc1pi, kord,
     :             izero, rz, inb1pi, work)
        v3pi= bvalu (t3pi, b3pi, nc3pi, kord,
     :             izero, rz, inb3pi, work)
        vvl(1)=v1sig
        vvl(2)=v1pi
        vvl(3)=v3sig
        vvl(4)=v3pi
      endif
      if (rz .gt. 11.25d0) then
* here for long range extrapolation of l^2 matrix elements
        vvl(5) = l21si(nl-1)*(rz/11.25d0)**1.98
        vvl(7) = l23si(nl-1)*(rz/11.25d0)**1.98
        vvl(6) = l21pi(nl-1)*(rz/11.25d0)**1.98
        vvl(8) = l23pi(nl-1)*(rz/11.25d0)**1.98
      else
* here for spline interpolation of l^2 matrix elements
        call dsplint (rl, l21si(1), dl21si(1), nl, rz, vvl(5))
        call dsplint (rl, l23si(1), dl23si(1), nl, rz, vvl(7))
        call dsplint (rl, l21pi(1), dl21pi(1), nl, rz, vvl(6))
        call dsplint (rl, l23pi(1), dl23pi(1), nl, rz, vvl(8))
      endif
      if (rz .gt. 5.d0) then
* here for long-range extrapolation of l coupling matrix elements
        vl=-37436.d0*exp(-2.9869*rz)+227.26*exp(-1.5714*rz)
        vvl(9)=one-vl
        vl=2.1103*exp(-1.4561*rz)
        vvl(10)=one-vl
      else
* here for spline interpolation of l coupling matrix elements
        call dsplint (rl, ly1(1), dly1(1), nl, rz, vvl(9))
        call dsplint (rl, ly3(1), dly3(1), nl, rz, vvl(10))
      end if
* subtract unity from l^2 matrix elements of pi states to get l-perp^2
      vvl(6)=vvl(6)-one
      vvl(8)=vvl(8)-one
* multiply ly matrix elements by -1/sqrt(2) to conform to brigitte's
* convention of l+
      fact=-one/sqrt(two)
      vvl(9)=fact*vvl(9)
      vvl(10)=fact*vvl(10)
* constant in vvl(11)
      vvl(11)=one
* duplicate vvl(9) in vvl(12) and vvl(10) in vvl(13)
      vvl(12)=vvl(9)
      vvl(13)=vvl(10)
* and multiply by spin-orbit constant (in hartree)
      fact=aso/econv
      call dscal(3,fact,vvl(11),1)
* multiply l coupling matrix elements by b
      b=one/(two*rmu*rz*rz)
      call dscal(6,b,vvl(5),1)
* vvl(14) just equal to B
      vvl(14)=b
      return
      end

c  -----------------------------------------------------------------
      subroutine ground(wf, r, nch, nphoto, mxphot)
c  -----------------------------------------------------------------
*  driven equation source vector for hcl
*  author:  millard alexander
*  current revision date:  29-jan-1993
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
      parameter (nl =  25)
      dimension wf(6), sc(6)
      common /cotrans/ t(6,6)
      common /cojtot/ j,jlpar
      common /coered/ ered, rmu
      dimension rl(nl), dip1(nl), ddip1(nl), dip3(nl)
* these are the morse parameters for hcl (avqz-mr-acpf with CASSCF
* orbitals optimized for each state
* note that g0 is more commonly denoted beta
      data de, r0, g0 /0.16811d0, 2.4154d0, 0.99471d0/
*  knots for spline expansion of 1Sig-1Pi and 3Sig-3Pi
*  transition moments
      data rl/ 1.8, 2.0, 2.2, 2.3, 2.4, 2.5, 2.6, 2.8, 3.0, 3.25,
     :  3.5, 3.75,
     : 4.0, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 8.25, 9.25,
     :  10.25,
     : 11.25, 30/
      data dip1 /6.2614624e-1, 5.7032945e-1, 4.5637519e-1, 3.6713607e-1,
     : 2.7828866e-1, 2.1192959e-1, 1.6771475e-1, 1.1915488e-1,
     : 9.602213e-2, 8.099689e-2, 7.163568e-2, 6.391378e-2,
     : 5.640394e-2, 4.87717e-2, 3.394984e-2, 2.177649e-2,
     : 1.336152e-2, 8.14642e-3, 5.08626e-3, 3.31558e-3,
     : 1.63871e-3, 9.6572e-4, 6.3707e-4, 4.3484e-4,
     : 8.38e-6/
      data dip3 /1.725617d-2, 1.4428889d-1, 1.6576706d-1, 1.555682d-1,
     : 1.4045505d-1, 1.244462d-1, 1.0909787d-1, 8.239537d-2,
     : 6.160453d-2, 4.292095d-2, 3.046575d-2, 2.229494d-2,
     : 1.695148d-2, 1.338722d-2, 9.15858d-3, 6.7703d-3,
     : 5.19709d-3, 4.07935d-3, 3.25825d-3, 2.64638d-3,
     : 1.39839d-3, 9.1192d-4, 6.1812d-4, 4.2624d-4,
     : 8.25d-6/
      data one, two /1.d0,2.d0/
      data interp / 1/
      if (interp .eq. 1) then
*  here for first pass through source subroutine
*  interpolate the coefficients for spline expansion of 1Sig-1Pi
*  transition moment
        call tcasea22(j,jlpar)
        nlm1 = nl -1
        ddip1f =( dip1(nl) - dip1(nlm1)) / (rl(nl) - rl(nlm1))
        ddip1i =( dip1(2) - dip1(1)) / (rl(2) - rl(1))
        call dspline( rl, dip1(1), nl, ddip1f, ddip1i, ddip1)
        interp=0
      endif
* determine morse function for ground vibrational state
      rmlmda = two*rmu
      d = dsqrt(rmlmda*de)/g0
      twod = two*d
      twodm1 = twod-one
      dln2d = log(twod)
      gln2d = gammln(twodm1)
      z = g0*(r-r0)
      arg = twodm1*(dln2d-z)-gln2d-twod*exp(-z)
      p = sqrt(g0*dexp(arg))
* determine transition moment
* here for long range extrapolation
      if (r .gt. 11.25d0) then
        dipfunc=0.d0
      else
* here for spline interpolation
        call dsplint (rl,dip1(1),ddip1(1),nl,r,dipfunc)
      end if
* initialize source vector to zero
      call dscal(6,zero,wf,1)
* since 1Pi state is channel 6, projection of ground state is only onto
* this channel
      wf(6)=rmlmda*p*dipfunc
      nch=6
      call mxma(t,1,nch,wf,1,nch,sc,1,nch,nch,nch,1)
      call dcopy(6,sc,1,wf,1)
      entry wfintern(wf,yymin,nnvib,nny)
*  dummy ground subroutine
      return
      end
* -----------------------------------------
      function gammln (xx)
      implicit double precision (a-h,o-z)
c
c     -----------------------------------------------------------------
c     natural log of the gamma function (numerical recipes)
c     -----------------------------------------------------------------
c
      dimension cof(6)
      data cof,stp /76.18009173d0,-86.50532033d0,24.01409822d0,
     + -1.231739516d0,0.120858003d-2,-0.536382d-5,2.50662827465d0/
c
      x = xx-1.0d0
      tmp = x+5.5d0
      tmp = (x+0.5d0)*dlog(tmp)-tmp
      ser = 1.0d0
      do 11 j = 1,6
         x = x+1.0d0
         ser = ser+cof(j)/x
  11  continue
      gammln = tmp+dlog(stp*ser)
      return
      end
* --------------------------------------------------------------------
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
      return
      end
* --------------------------------------------------------------------
      subroutine bausr (j, l, is, jhold, ehold, ishold, nlevel, nlevop,
     :                  sc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu,
     :                  csflag, clist, bastst, ihomo, nu, numin, jlpar,
     :                  n, nmax, ntop)
      return
      end
* --------------------------------------------------
