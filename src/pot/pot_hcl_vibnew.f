*system: H(2S)+Cl(2P)
*reference:  M. H. Alexander, B. Pouilly, and T. Duhoo, 
* J. Chem. Phys. 99, 1752 (1993)


*  to run this pot routine, the eispack routine RS must be replaced
*  with the lapack routine DSYEV, to be compatible with the ifort
*  compiler and MKL library


      subroutine driver
      implicit double precision (a-h,o-z)
      common /comxm/ncache, mxmblk
      common /cojtot/ jtot,jlpar
      common /covvl/ vvl(16)
      common /cosysi/ nscode,isicod,iscod(2),nvib
      dimension wf(6)
      include "common/parpot"
      if (ncache .eq. 0) then
        ncache=4096
        mxmblk=64
      endif
      jlpar=1
      jtot=2
      print *, 'input vibrational quantum number'
      read (5, *, end=99) nvib
      potnam='ALEXANDER-POUILLY-DUHOO HCl(X,A,a,triplet)'
      print *, potnam
      print *, 'r, v1sig, v1pi, v3sig, v3pi'
      do ir=1,100
        r=1.5+(ir-1)*0.25d0
        if (r .gt. 20) goto 99
      call pot(vv0,r)
      write (6, 100) r,(vvl(i), i=1,4)
      nch=1
      nphoto=1
      mxphot=1
      call ground(wf, r, nch, nphoto, mxphot)
100   format(f8.3,6(1pe16.8))
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
*     vvl(15):  1sigma dipole moment
*     vvl(16):  1pi dipole moment
*  subroutines called:
*   bint4,bvalu  (spline routines from cmlib bspline package)

*   splint  (from numerical recipes)
*      to return a cubic-spline interpolated value for a given value rz
* author:  millard alexander
* latest revision date:  3-mar-1997 by mha
* ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      double precision l21si, l21pi, l23si, l23pi, ly1, ly3
      parameter (nl = 25, nn=27, n3sig=24)
      parameter (nnp6=nn+6, n3p6=n3sig+6)
      parameter (ndip=32, ndipp6=ndip+6, nwdip=5*ndipp6)
      parameter (nw=5*max(nnp6,ndipp6))
      common /coered/ ered, rmu
      common /cosysr/ junk1, junk2, aso
      common /covvl/ vvl(16)
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
      dimension rdip(ndip), dippi(ndip), dipsig(ndip),ddipsig(ndip),
     :          ddippi(ndip)
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
* data for permanent dipole moments
      data rdip /
     : 2.0d0, 2.2d0, 2.4d0, 2.6d0, 2.8d0, 3.0d0, 3.2d0, 3.4d0, 3.6d0,
     : 3.8d0, 4.0d0, 4.25d0, 4.5d0, 4.75d0, 5.0d0, 5.25d0, 5.5d0,
     : 5.75d0, 6.0d0, 6.25d0, 6.5d0, 6.75d0, 7.0d0, 7.25d0, 7.5d0,
     : 7.75d0, 8.0d0, 8.5d0, 9.0d0, 9.5d0, 10d0, 10.5d0/
       data dipsig /
     : 3.9485816d-1, 4.0386723d-1, 4.1847469d-1, 4.4081023d-1,
     : 4.6286588d-1, 4.8070598d-1, 4.9110068d-1, 4.9119111d-1,
     : 4.7893768d-1, 4.5368216d-1, 4.1651174d-1, 3.5767255d-1,
     : 2.9247929d-1, 2.2891638d-1, 1.728664d-1, 1.2708756d-1,
     : 9.177669d-2, 6.562553d-2, 4.678357d-2, 3.343993d-2, 2.412478d-2,
     : 1.757989d-2, 1.300383d-2, 9.7897d-3, 7.51356d-3, 5.88761d-3,
     : 4.69332d-3, 3.16786d-3, 2.27427d-3, 1.71468d-3, 1.34235d-3,
     : 1.0804d-3/
       data dippi /
     : -1.3540803, -1.1640610, -9.1776846d-1, -7.5380637d-1,
     : -6.2626371d-1, -5.1740307d-1, -4.2248918d-1, -3.4049466d-1,
     : -2.7103326d-1, -2.1342039d-1, -1.6661273d-1, -1.2099645d-1,
     : -8.716796d-2, -6.246361d-2, -4.466574d-2, -3.195369d-2,
     : -2.293758d-2, -1.657437d-2, -1.209015d-2, -8.93597d-3,
     : -6.71633d-3, -5.14756d-3, -4.02968d-3, -3.22258d-3, -2.62615d-3,
     : -2.17699d-3, -1.85158d-3, -1.38057d-3, -1.07202d-3, -8.6224d-4,
     : -7.0728d-4, -5.8544d-4/

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
*  interpolate the  coefficients for the permanent dipole moments
        ndipm1=ndip-1
        dsigf=(dipsig(ndip)-dipsig(ndipm1))/(rdip(ndip)-rdip(ndipm1))
        dsigi=(dipsig(2)-dipsig(1))/(rdip(2)-rdip(1))
        call dspline( rdip, dipsig(1), ndip, dsigi, dsigf, ddipsig)
        dpif=(dippi(ndip)-dippi(ndipm1))/(rdip(ndip)-rdip(ndipm1))
        dpii=(dippi(2)-dippi(1))/(rdip(2)-rdip(1))
        call dspline( rdip, dipsig(1), ndip, dsigi, dsigf, ddippi)
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
* multiply ly matrix elements by -1/sqrt(2) to conform to brigittes
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
* now for dipole moment function
* constant extrapolation at short range
      if (rz .lt. 2d0) then
        vvl(15)=dipsig(1)
        vvl(16)=dippi(1)
* exponential extrapolation at long-range
      elseif (rz .gt. 10.5) then
        vvl(15)=1.0315d-01*exp(-0.43418*rz)
        vvl(16)=-3.1029e-02*exp(-0.37813*rz)
      else
        call dsplint (rdip, dipsig(1), ddipsig(1), ndip, rz, vvl(15))
        call dsplint (rdip, dippi(1), ddippi(1), ndip, rz, vvl(16))
      endif
      return
      end

c  -----------------------------------------------------------------
      subroutine ground(wf, r, nch, nphoto, mxphot)
c  -----------------------------------------------------------------
*  driven equation source vector for hcl, now includes v=5
*  author:  millard alexander
*  current revision date:  17-mar-1997
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
      character*1 forma, formb
      parameter (nl =  25, npsi = 41)
      dimension wf(6), sc(6)
      common /cotrans/ t(6,6)
      common /cojtot/ j,jlpar
      common /cosysi/ nscode,isicod,iscod(2),nvib,ibran
      common /codim/ mairy,mmax
      common /coeint/ eint(6)
      common /cocent/ cent(6)
      common /coered/ ered, rmu
      common /covvl/ vvl(16)
      dimension w(49),eignow(7),vecnow(49),scr(49),hp(49)
      dimension nsc(6)
      dimension rl(nl), dip1(nl), ddip1(nl), dip3(nl)
      dimension r5(npsi), rpsi(npsi), rx(npsi), psi(npsi),
     :  psi0(npsi), psi1(npsi), psi2(npsi), psi3(npsi),
     :  psi4(npsi), psi5(npsi), psi6(npsi), dpsi(npsi)
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
      data r5 /
     : 1.59d0, 1.65d0, 1.71d0, 1.77d0, 1.83d0, 1.89d0, 1.95d0, 2.01d0,
     : 2.07d0, 2.13d0, 2.19d0, 2.25d0, 2.31d0, 2.37d0, 2.43d0, 2.49d0,
     : 2.55d0, 2.61d0, 2.67d0, 2.73d0, 2.79d0, 2.85d0, 2.91d0, 2.97d0,
     : 3.03d0, 3.09d0, 3.15d0, 3.21d0, 3.27d0, 3.33d0, 3.39d0, 3.45d0,
     : 3.51d0, 3.57d0, 3.63d0, 3.69d0, 3.75d0, 3.81d0, 3.87d0, 3.93d0,
     : 3.99d0/
      data psi5 /
     : 1.3613206d-3, 6.9015243d-3, 2.7837145d-2, 9.0163368d-2,
     : 2.3531518d-1, 4.9316917d-1, 8.1964719d-1, 1.0460307d0,
     : 9.3621255d-1, 3.8629551d-1, -3.8153211d-1, -8.8047518d-1,
     : -7.2651757d-1, -1.0673546d-2, 7.1203495d-1, 8.6546769d-1,
     : 3.3282918d-1, -4.6102774d-1, -9.1172720d-1, -7.1045865d-1,
     : -2.8307611d-2, 6.7834718d-1, 1.0043884d0, 8.1468443d-1,
     : 2.4857236d-1, -4.2319931d-1, -9.6018497d-1, -1.2419150d0,
     : -1.2702391d0, -1.1204855d0, -8.8594939d-1, -6.4167271d-1,
     :  -4.3172392d-1, -2.7251899d-1, -1.6260803d-1, -9.2262472d-2,
     : -5.0023546d-2, -2.6025082d-2, -1.3038555d-2, -6.3102129d-3,
     : -2.9583018d-3/
      data rpsi /
     : 1.65d0, 1.71d0, 1.77d0, 1.83d0, 1.89d0, 1.95d0, 2.01d0,
     : 2.07d0, 2.13d0, 2.19d0, 2.25d0, 2.31d0, 2.37d0, 2.43d0,
     : 2.49d0, 2.55d0, 2.61d0, 2.67d0, 2.73d0, 2.79d0, 2.85d0,
     : 2.91d0, 2.97d0, 3.03d0, 3.09d0, 3.15d0, 3.21d0, 3.27d0,
     : 3.33d0, 3.39d0, 3.45d0, 3.51d0, 3.57d0, 3.63d0, 3.69d0,
     : 3.75d0, 3.81d0, 3.87d0, 3.93d0, 3.99d0, 4.05d0/
      data psi0 /
     : 8.9050635d-5, 4.5942723d-4, 1.9711971d-3, 7.1206810d-3,
     : 2.1880341d-2, 5.7763556d-2, 1.3229343d-1, 2.6525127d-1,
     : 4.6955283d-1, 7.3991452d-1, 1.046066d0, 1.3361546d0,
     : 1.551097d0, 1.6458369d0, 1.6062232d0, 1.4492814d0,
     : 1.2108703d0, 9.3739353d-1, 6.7483860d-1, 4.5438537d-1,
     : 2.8772006d-1, 1.7204959d-1, 9.7470354d-2, 5.2451413d-2,
     : 2.6873413d-2, 1.3139500d-2, 6.1453752d-3, 2.7559085d-3,
     : 1.1877863d-3, 4.9307445d-4, 1.9754518d-4, 7.6529342d-5,
     : 2.8719857d-5, 1.0458662d-5, 3.7018944d-6, 1.2755980d-6,
     : 4.2855709d-7, 1.4058597d-7, 4.5093854d-8, 1.4161605d-8,
     : 4.3599509d-9/
      data psi1 /
     : -3.3859285d-4, -1.6619428d-3, -6.7421868d-3, -2.2860243d-2,
     : -6.5350050d-2, -1.5875237d-1, -3.2998218d-1, -5.8986793d-1,
     : -9.0897257d-1, -1.2056693d0, -1.3640752d0, -1.2811489d0,
     : -9.2024165d-1, -3.4008707d-1, 3.2326411d-1, 9.1206791d-1,
     : 1.3034964d0, 1.4473317d0, 1.3733763d0, 1.1602507d0,
     : 8.9244187d-1, 6.3350222d-1, 4.1873410d-1, 2.5937406d-1,
     : 1.5131810d-1, 8.3502813d-2, 4.3755892d-2, 2.1849483d-2,
     : 1.0431020d-2, 4.7748069d-3, 2.1011902d-3, 8.9103650d-4,
     : 3.6492889d-4, 1.4464379d-4, 5.5591976d-5, 2.0755974d-5,
     : 7.5413904d-6, 2.6708752d-6, 9.2347715d-7, 3.1218436d-7,
     : 1.0332848d-7/
      data psi2 /
     : -9.1072072d-4, -4.2511793d-3, -1.6292527d-2, -5.1767178d-2,
     : -1.3728310d-1, -3.0537393d-1, -5.7118138d-1, -8.9645733d-1,
     : -1.1681515d0, -1.2278504d0, -9.5733545d-1, -3.7211497d-1,
     : 3.4829855d-1, 9.3012540d-1, 1.142078d0, 9.0913101d-1,
     : 3.3522771d-1, -3.5877134d-1, -9.4763507d-1, -1.295778d0,
     : -1.3805679d0, -1.2602794d0, -1.0255367d0, -7.5984824d-1,
     : -5.1951100d-1, -3.3088494d-1, -1.9776667d-1, -1.1158776d-1,
     : -5.9736325d-2, -3.0468183d-2, -1.4859721d-2, -6.9520161d-3,
     : -3.1288261d-3, -1.3581288d-3, -5.6992198d-4, -2.3171672d-4,
     : -9.1465337d-5, -3.5118936d-5, -1.3139675d-5, -4.7985735d-6,
     : -1.7132024d-6/
      data psi3 /
     : -2.0100744d-3, -8.9350507d-3, -3.2383784d-2, -9.6471202d-2,
     : -2.3719565d-1, -4.8181652d-1, -8.0518028d-1, -1.090710d0,
     : -1.1509948d0, -8.3776654d-1, -1.8303416d-1, 5.5250577d-1,
     : 9.9986223d-1, 9.1093278d-1, 3.2993708d-1, -4.2806727d-1,
     : -9.6759768d-1, -1.0345498d0, -6.3258315d-1, 3.0084876d-2,
     : 6.8860809d-1, 1.1482752d0, 1.3363828d0, 1.285712d0,
     : 1.0841362d0, 8.2483040d-1, 5.7623256d-1, 3.7409163d-1,
     : 2.2768195d-1, 1.3079482d-1, 7.1305232d-2, 3.7057873d-2,
     : 1.8430711d-2, 8.8017101d-3, 4.0481286d-3, 1.7979508d-3,
     : 7.7305351d-4, 3.2249989d-4, 1.3080920d-4, 5.1685691d-5,
     : 1.9929845d-5/
      data psi4/
     : 3.9053479d-3, 1.6531280d-2, 5.6634330d-2, 1.5796786d-1,
     : 3.5903296d-1, 6.6186260d-1, 9.7501552d-1, 1.1035578d0,
     : 8.5238097d-1, 2.1214809d-1, -5.3587095d-1, -9.5017981d-1,
     : -7.4564143d-1, -4.2201160d-2, 6.9011841d-1, 9.6395250d-1,
     : 6.1002223d-1, -1.2911149d-1, -7.9807330d-1, -1.0386623d0,
     : -7.6900449d-1, -1.5840569d-1, 5.1937555d-1, 1.035976d0,
     : 1.2860232d0, 1.283332d0, 1.1102909d0, 8.6259377d-1,
     : 6.1413571d-1, 4.0605120d-1, 2.5169335d-1, 1.4732164d-1,
     : 8.1895333d-2, 4.3441861d-2, 2.2078239d-2, 1.0788247d-2,
     : 5.0841280d-3, 2.3172024d-3, 1.0239437d-3, 4.3968213d-4,
     : 1.8384753d-4/
      data psi6/
     : 1.1225189d-2, 4.3172601d-2, 1.3223463d-1, 3.2268823d-1,
     : 6.2181257d-1, 9.2391050d-1, 9.9504605d-1, 6.2531388d-1,
     : -1.0772864d-1, -7.4985666d-1, -7.8611689d-1, -1.5053160d-1,
     : 6.1569308d-1, 8.1848975d-1, 2.7062039d-1, -5.2479437d-1,
     : -8.5277355d-1, -4.3774239d-1, 3.4803178d-1, 8.6472451d-1,
     : 7.4127390d-1, 1.0123974d-1, -6.1063001d-1, -9.7290053d-1,
     : -8.2683945d-1, -2.9407832d-1, 3.6333940d-1, 9.0541955d-1,
     : 1.2058113d0, 1.2575266d0, 1.1281568d0, 9.0684262d-1,
     : 6.6807868d-1, 4.5766866d-1, 2.9453534d-1, 1.7944346d-1,
     : 1.0412677d-1, 5.7835927d-2, 3.0877507d-2, 1.5901912d-2,
     : 7.9246039d-3/

      data zero,one,two,three,xmin1 /0.d0,1.d0,2.d0,3.d0,-1.d0/
      data interp / 1/
*  reset for first pass if vibrational quantum number has changed
      if (nvib.ne.nvibhold) interp=1
      call tcasea22(j,jlpar)
      if (interp .eq. 1) then
*  here for first pass through source subroutine
*  interpolate the coefficients for spline expansion of 1Sig-1Pi
*  transition moment

        nlm1 = nl -1
        ddip1f =( dip1(nl) - dip1(nlm1)) / (rl(nl) - rl(nlm1))
        ddip1i =( dip1(2) - dip1(1)) / (rl(2) - rl(1))
        call dspline( rl, dip1(1), nl, ddip1f, ddip1i, ddip1)
*  interpolate the coefficients for spline expansion of psi
        if (nvib .eq. 5) then
           call dcopy(npsi,r5,1,rx,1)
           call dcopy(npsi,psi5,1,psi,1)
        else
           call dcopy(npsi,rpsi,1,rx,1)
           if (nvib .eq. 0) call dcopy(npsi,psi0,1,psi,1)
           if (nvib .eq. 1) call dcopy(npsi,psi1,1,psi,1)
           if (nvib .eq. 2) call dcopy(npsi,psi2,1,psi,1)
           if (nvib .eq. 3) call dcopy(npsi,psi3,1,psi,1)
           if (nvib .eq. 4) call dcopy(npsi,psi4,1,psi,1)
           if (nvib .eq. 6) call dcopy(npsi,psi6,1,psi,1)
        endif
        npsim1 = npsi -1
        dpsif =( psi(npsi) - psi(npsim1)) / (rx(npsi) - rx(npsim1))
        dpsii =( psi(2) - psi(1)) / (rx(2) - rx(1))
        call dspline( rx, psi(1), npsi, dpsif, dpsii, dpsi)
        interp=0
        nvibhold=nvib
      endif
      rmlmda = two*rmu
      if (r.lt.rx(1)) then
        p=0.d0
      elseif (r. gt. rx(npsi)) then
        p=0.d0
      else
        call dsplint(rx,psi,dpsi,npsi,r,p)
      endif
* for diagnostic purposes uncomment the next line
*      if (r .lt. 4) write (6,100) p
*100   format (f8.4)
* determine transition moment
* here for long range extrapolation
      if (r .gt. 11.25d0) then
        dipfunc=0.d0
      else
* here for spline interpolation
        call dsplint (rl,dip1(1),ddip1(1),nl,r,dipfunc)
      end if
      nch=6
      nmax=mmax
      nmax=nch+1
      nmaxm1=nmax-1
      if (abs(r-2.4) .lt. 0.01) then
        xxx=1
      endif
      if (abs(r-3.37) .lt. 0.01) then
        xxx=1
      endif
      call potmath(w,r,nch,nmax)
*  since potmat returns negative of lower triangle of w(r) matrix (eq.(3) of
*  m.h. alexander, "hybrid quantum scattering algorithms ..."),
*  next loop changes its sign
      ipt = 1
      do 100 icol = 1, nch
*  nrow is the number of (diagonal plus subdiagonal) elements in column icol
*  ipt points to the diagonal element in column icol for a matrix stored in
*  packed column form
        nrow = nch - icol + 1
        call dscal (nrow, xmin1, w(ipt), 1)
        ipt = ipt + nmaxp1
 100  continue
*  next loop fills in upper triangles of w
      if (nch .gt. 1) then
        ipt = 2
        do 110 icol = 1, nch -1
*  ipt points to the first subdiagonal element in column icol
*  nrow is the number of subdiagonal elements in column icol
          nrow = nch - icol
          call dcopy (nrow, w(ipt), 1, w(ipt + nmaxm1), nmax)
          ipt = ipt + nmaxp1
110     continue
      end if
*  diagonalize w at rnow, the vecnow matrix is identical to the tn matrix
*  of eq.(6) of m.h. alexander, "hybrid quantum scattering algorithms ..."
*  now call eispack eigenvalue and eigenvector routine (hp is used as
*  a scratch vector here)
      call rs (nmax, nch, w, eignow, ione, vecnow, scr, hp, ierr)
      if (ierr .ne. 0) then
        write (6, 115) ierr
        write (9, 115) ierr
115     format (' *** IERR =',i3,' IN AIRPRP/POTENT/RS;  ABORT ***')
        write (9, 120) (eignow (i), i=1, nch)
120     format (' EIGENVALUES ARE:',/,8(1pe16.8) )
        call exit
      end if
* sort eigenvectors so that they correspond to nominal case(a)
* ordering (see hiba22p.f)
      do 125 i=1,6
        call dcopy(nmax,vecnow(i),nmax,sc,1)
        nnn=idamax(6,sc,1)
        nsc(i)=nnn
125   continue
      call dscal(nmax*nmax,zero,w,1)
      do 130 i=1,6
        ind=(nsc(i)-1)*nmax+1
        call dcopy(6,vecnow(ind),1,w((i-1)*nmax+1),1)
130   continue
* w now contains sorted eigenvectors

* construct 1Pi contribution to source vector
* NB hoenl-london factors assume different initial j's, for a given final
* j.  see herzberg, diatomics, p. 208 using the expressions with j"
      if (ibran .eq. -1) then
        hlfactor=j-1
        hlfacsig=j
      end if
      if (ibran .eq. 0) then
        hlfactor=2*j+1
        hlfacsig=zero
      end if
      if (ibran .eq. 1) then
        hlfactor=j+2
        hlfacsig=j+1
      end if
      hlfactor=sqrt(hlfactor/(2*j+1))
      hlfacsig=sqrt(hlfacsig/(2*j+1))
      call dset(6,zero,wf,1)
* sixth row is the 1Pi component; this is the source vector
      wf(6)=hlfactor*dipfunc

* add 1sig-1sig component if e-levels
      if (jlpar .eq. -1) then
         call dcopy(6,w(5),nmax,hp,1)
         hp(5)=zero
         call mxva(w,1,7,hp,1,vecnow,1,6,6)
* if e-labelled states, then add the component of the 1sig state
* multiplied by the 1sig moment and hlfactor for sig-sig
* transitions
         call daxpy(6,hlfacsig*vvl(15),vecnow,1,wf,1)
      endif
      call dscal(6,rmlmda*p,wf,1)
      nch=6
      call mxma(t,1,nch,wf,1,nch,sc,1,nch,nch,nch,1)
      call dcopy(6,sc,1,wf,1)
      entry wfintern(wf,yymin,nnvib,nny)
*  dummy ground subroutine
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
*
*  subroutine tcasea22 called by this pot routine.  this subr is part
*  of the basisi routine ba22p in hiba22p.f
*
*  this subr was deleted from this file
*  p. dagdigian, 23-jan-2011
*
* -----------------------------------------
      subroutine potmath (w, r, nch, nmax)
*  determine negative of lower triangle of w(r) matrix.  see eq. (3) of
*  m.h. alexander, "hybrid quantum scattering algorithms ..."
*  just for Hcl (2p-2s), returns w(r) in Hund's case (a) basis
*  author:  millard alexander
*  latest revision date:  10-mar-1997 by mha
*  -------------------------------------------------------------------
*  variables in call list:
*    w:        matrix of maximum row dimension nmax
*              on output contains negative of lower triangle of w(r)
*    r:        value of interparticle distance at which -w(r) is to
*              be evaluated
*    nch:      actual size of matrix w
*    nmax:     maximum row and column dimension of matrix w
*  variable in common block /cov2/
*    nv2max:    maximum core memory allocated for the v2 matrix
*    ndummy:    dummy variable for alignment
*    v2:        lower triangle of nonzero elements of coupling matrix stored
*               in packed row form that is (1,1), (2,1), (2,2), (3,1), (3,2),
*               (3,3)  ...
*  variable in common block /coiv2/
*   lamnum:     number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx
*   iv2:        row+column index of v2 matrix for each non-zero element
*  variable in common block /cocent/
*    cent:      array containing centrifugal barrier of each channel
*  variable in common block /coeint/
*    eint:      array containing channel energies (in hartree)
*  variables in common block /coered/
*    ered:      collision energy in atomic units (hartrees)
*    rmu:       collision reduced mass in atomic units (mass of electron = 1)
*  variable in common block /conlam/
*    nlammx:    the maximum number of angular coupling terms allowed
*    nlam:      the total number of angular coupling terms used
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx
*  variables in common block /copmat/
*    rtmn,rtmx: minimum and maximum turning points (not used here)
*    iflag:     variable used in determination of turning points (not used her
*           iflag = 0 if all channels are in classically forbidden region
*           iflag = 1 if some channels are open
*           iflag = 2 if all asymptotically open channels are open at r
*  variable in common block /covvl/
*    vvl        array to store r-dependence of each angular term in the
*               potential
*  subroutines called:
*    pot:      returns r-dependence of each angular term in the potential
*    daxpy:    blas routine
*    vsmul:    multiplies vector by scalar and stores result in another
*              vector
*  -------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c      real*4 cent, eint, v2, vvl, w
c      real*4 ered, factor, r, r2, rmu, rtmn, rtmx, vv0, wmax, wmin,
c     :     one, two, zero
      integer i, iflag, ilam, lamnum, nch, nlam, nlammx, nmax, nv2max
      integer iv2
      common /cputim/ cpuld,cpuai,cpupot,cpusmt,cpupht
      common /cov2/ nv2max, ndummy, v2(1)
      common /coiv2/ iv2(100)
      common /conlam/ nlam, nlammx, lamnum(1)
      common /copmat/ rtmn, rtmx, iflag
      common /cocent/ cent(1)
      common /coeint/ eint(1)
      common /coered/ ered, rmu
      common /covvl/ vvl(1)
      common /coselb/ ibasty
      dimension w(49)
      zero = 0.d0
      one = 1.d0
      two = 2.d0
cABER only for testing potential-matrix
c     if (r.le.3.41) icount=0
cABER
*  calculate coefficients of each angular term
      cpupot=cpupot-second()
      call pot( vv0, r)
*  vv0 is the isotropic term in the potential
*  the coefficients for each angular term in the coupling potential
*  [ vvl(i) for i = 1, nlam ] are returned in common block vvlcontains
*  multiply all vvl terms by twice the reduced mass
      twormu = two * rmu
      call dscal(nlam, twormu, vvl, 1)
*  now loop over angular coupling matrix to calculate the potential matrix
*    w(ij) = 2 * rmu * vv0 + sum [ 2 * rmu * vvl(ilam) * v2(ij,ilam) ]
*  first zero out lower triangle of potential matrix
      iwpt = 1
      do 20 icol = 1, nch
        ncol = nch - icol + 1
        call dset(ncol, zero, w(iwpt), 1)
        iwpt = iwpt + nmax + 1
20    continue
      ioff = 0
cABER
c       call druckq(w,nmax,nch,'potential matrix',icount)
cABER
      do  40  ilam = 1, nlam
           do  30 i = 1, lamnum(ilam)
             w(iv2(ioff+i)) = w(iv2(ioff+i)) + v2(ioff+i) * vvl(ilam)
30         continue
         ioff = ioff + lamnum(ilam)
40    continue
* here for 2s-2p scattering
* fill in upper triangle of w matrix
*  first fill in upper half of original matrix
      if (ibasty .eq. 10) then
        nmaxp1 = nmax + 1
        icolpt = 2
        irowpt = nmaxp1
        do  150  icol = 1, nch - 1
*  icolpt points to first sub-diagonal element in column icol
*  irowpt points to first super-diagonal element in row icol
*  ncol is number of subdiagonal elements in column icol
          ncol = nch - icol
          call dcopy (ncol, w(icolpt), 1, w(irowpt), nmax)
          icolpt = icolpt + nmaxp1
          irowpt = irowpt + nmaxp1
150     continue
      endif
      return
      end
