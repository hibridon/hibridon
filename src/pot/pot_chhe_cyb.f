* System:  CH(X 2Pi)+He
* New Cybulski, Chalasinski, Szczesniak PES
* Reference: S. M. Cybulski, G. Chalasinski, M. M. Szczesniak, J. Chem.
*    Phys. 105, 9525 (1996)
*
      include "common/syusr"
* ------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(19)
      econv=219474.6d0
1      print *, ' r (bohr)'
      read (5, *, end=99) r
      call pot(vv0,r)
      if (r.le.0) goto 99
      write (6, 100) r,vv0,vvl
100   format(' vsum',/,5(1pe16.8),/,5e16.8,/,e16.8,/,
     :    '  vdif',/,5e16.8,/,4e16.8)
      goto 1
99    rr=3.0d0
      dr=0.1d0
      open (unit=12,file='chhe_cyb_vsum.dat')
      write(12,109)
      open (unit=13,file='chhe_cyb_vdif.dat')
      write(13,109)
109   format(' %R/bohr V00  ...')
      do i=1,250
        call pot(vv0,rr)
        write (12,110) rr,vv0*econv,(econv*vvl(j),j=1,9)
110     format(f7.2,10(1pe16.8))
        write (13,112) rr,(econv*vvl(j),j=10,18)
112     format(f7.2,9(1pe16.8))
        rr = rr + dr
      enddo
      close(12)
      close(13)
      end
      include "common/bausr"
      include "common/ground"
      subroutine loapot(iunit,filnam)
* ------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      common /conlam/ nlam, nlammx, lamnum(1)
      potnam='Cybulski He-CH(X)'
      lammin(1)=1
      lammax(1)=6
      lammin(2)=2
      lammax(2)=6
      nlam=6+6-2+1
      mproj(1)=0
      mproj(2)=2
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* ----------------------------------------------------------------------
      subroutine pot ( vv0, rr)
*
*  subroutine to calculate the r-dependent coefficients in the
*  He-CH(A) potentials. deepened and shifted
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:        interparticle distance (in atomic units)
*  on return:
*    vv0:      contains isotropic term (d00 term in vsum)
*  variable in common block /covvl/
*    vvl:      vector of length 11 to store r-dependence of each term
*              in potential expansion
*    vvl(1-6)  expansion coefficients in dl0 (l=1:10) of vsum
*    vvl(7-11) expansion coefficients in dl2 (l=11:19) of vdif

* uses linear least squares routines from cmlib

* author:  millard alexander
* note:    modified for he-ch by mark cybulski 7/9/96
*          slightly revised by mha 1/29/97
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension xlam1(22),xlam2(22),r0(22),c1(22),c2(22),c3(22),
     :          clr(22),vsum(11),xsum(11),vdif(11),xdif(11),
     :          vap(11),va2p(11),d0(121),d2(81),aa(121)
      dimension kpvt(12),qraux(12),work(77),rsd(11), re(22)

      common /covvl/ vvl(19)
      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
      data onemin /-1.d0/

* mark cybulski's parameters:
*     npoint is the number of points for which data are available
*     nparm is the number of parameters in the fitted function
      data nparm, npoint / 7, 11/

*
*     for distances beyond rmax difference potential is damped
*
      data rmax /20d0/
*
*     re's for ab initio potentials
*
      data re /
*     for the A' state
     :    7.7359, 7.7391, 7.8032, 7.8586, 7.8483, 7.7926, 7.7016,
     :    7.5445, 7.7116, 8.0488, 8.1561,
*     for the A" state
     :    7.7359, 7.7079, 7.5580, 7.1826, 6.2665, 5.6028, 5.2618,
     :    6.0808, 7.3292, 7.9613, 8.1561 /

*      coefficicients for d0 rotation matrices
*      stored (by column) for each of 11 angles and for L=0:6
*      angles are 0 20 40 60 80 90 100 120 140 160 180

      data d0/
*      L = 0
     *  1.0000000d+00,  1.0000000d+00,  1.0000000d+00,  1.0000000d+00,
     *  1.0000000d+00,  1.0000000d+00,  1.0000000d+00,  1.0000000d+00,
     *  1.0000000d+00,  1.0000000d+00,  1.0000000d+00,
*      L = 1
     *  1.0000000d+00,  9.3969262d-01,  7.6604444d-01,  5.0000000d-01,
     *  1.7364818d-01,  5.0532155d-16, -1.7364818d-01, -5.0000000d-01,
     * -7.6604444d-01, -9.3969262d-01, -1.0000000d+00,
*      L = 2
     *  1.0000000d+00,  8.2453333d-01,  3.8023613d-01, -1.2500000d-01,
     * -4.5476947d-01, -5.0000000d-01, -4.5476947d-01, -1.2500000d-01,
     *  3.8023613d-01,  8.2453333d-01,  1.0000000d+00,
*      L = 3
     *  1.0000000d+00,  6.6488473d-01, -2.5233334d-02, -4.3750000d-01,
     * -2.4738193d-01, -7.5798232d-16,  2.4738193d-01,  4.3750000d-01,
     *  2.5233334d-02, -6.6488473d-01, -1.0000000d+00,
*      L = 4
     *  1.0000000d+00,  4.7497774d-01, -3.1900435d-01, -2.8906250d-01,
     *  2.6590161d-01,  3.7500000d-01,  2.6590161d-01, -2.8906250d-01,
     * -3.1900435d-01,  4.7497774d-01,  1.0000000d+00,
*      L = 5
     *  1.0000000d+00,  2.7149175d-01, -4.1968205d-01,  8.9843750d-02,
     *  2.8101754d-01,  9.4747791d-16, -2.8101754d-01, -8.9843750d-02,
     *  4.1968205d-01, -2.7149175d-01, -1.0000000d+00,
*      L = 6
     *  1.0000000d+00,  7.1903002d-02, -3.2357073d-01,  3.2324219d-01,
     * -1.3212134d-01, -3.1250000d-01, -1.3212134d-01,  3.2324219d-01,
     * -3.2357073d-01,  7.1903002d-02,  1.0000000d+00,
*      L = 7
     *  1.0000000d+00, -1.0722616d-01, -1.0060171d-01,  2.2314453d-01,
     * -2.8347992d-01, -1.1053909d-15,  2.8347992d-01, -2.2314453d-01,
     *  1.0060171d-01,  1.0722616d-01, -1.0000000d+00,
*      L = 8
     *  1.0000000d+00, -2.5183943d-01,  1.3862680d-01, -7.3638916d-02,
     *  2.3307850d-02,  2.7343750d-01,  2.3307850d-02, -7.3638916d-02,
     *  1.3862680d-01, -2.5183943d-01,  1.0000000d+00,
*      L = 9
     *  1.0000000d+00, -3.5169654d-01,  2.9001295d-01, -2.6789856d-01,
     *  2.5962717d-01,  1.2435648d-15, -2.5962717d-01,  2.6789856d-01,
     * -2.9001295d-01,  3.5169654d-01, -1.0000000d+00,
*      L =10
     *  1.0000000d+00, -4.0126914d-01,  2.9734522d-01, -1.8822861d-01,
     *  6.4682128d-02, -2.4609375d-01,  6.4682128d-02, -1.8822861d-01,
     *  2.9734522d-01, -4.0126914d-01,  1.0000000d+00 /

*     coefficicients for d2 rotation matrices
*     stored (by column) for each of 9 angles and for L=2:6
*     angles are 20 40 60 80 90 100 120 140 160

      data d2/
*      L = 2
     *  7.1633967d-02,  2.5301754d-01,  4.5927933d-01,  5.9390715d-01,
     *  6.1237244d-01,  5.9390715d-01,  4.5927933d-01,  2.5301754d-01,
     *  7.1633967d-02,
*      L = 3
     *  1.5051848d-01,  4.3340069d-01,  5.1348990d-01,  2.3060769d-01,
     *  6.9194003d-16, -2.3060769d-01, -5.1348990d-01, -4.3340069d-01,
     * -1.5051848d-01,
*      L = 4
     *  2.3957418d-01,  5.0756736d-01,  2.2234765d-01, -3.0244624d-01,
     * -3.9528471d-01, -3.0244624d-01,  2.2234765d-01,  5.0756736d-01,
     *  2.3957418d-01,
*      L = 5
     *  3.2835759d-01,  4.3600553d-01, -1.6982082d-01, -2.7746877d-01,
     * -9.1535062d-16,  2.7746877d-01,  1.6982082d-01, -4.3600553d-01,
     * -3.2835759d-01,
*      L = 6
     *  4.0592179d-01,  2.3830028d-01, -3.4523418d-01,  1.5131756d-01,
     *  3.2021721d-01,  1.5131756d-01, -3.4523418d-01,  2.3830028d-01,
     *  4.0592179d-01,
*      L = 7
     *  4.6231022d-01, -1.3906540d-02, -1.9131358d-01,  2.8490318d-01,
     *  1.0854723d-15, -2.8490318d-01,  1.9131358d-01,  1.3906540d-02,
     * -4.6231022d-01,
*      L = 8
     *  4.8973053d-01, -2.2700359d-01,  1.1374299d-01, -3.5240961d-02,
     * -2.7731624d-01, -3.5240961d-02,  1.1374299d-01, -2.2700359d-01,
     *  4.8973053d-01,
*      L = 9
     *  4.8345441d-01, -3.2461587d-01,  2.7905801d-01, -2.6334950d-01,
     * -1.2296697d-15,  2.6334950d-01, -2.7905801d-01,  3.2461587d-01,
     * -4.8345441d-01,
*      L =10
     *  4.4236809d-01, -2.7891371d-01,  1.6870456d-01, -5.7117496d-02,
     *  2.4836194d-01, -5.7117496d-02,  1.6870456d-01, -2.7891371d-01,
     *  4.4236809d-01 /

*     coefficients for expansion of V(A')(first 11 entries) and
*     for V(A") (entries 12:22)

      data xlam1/
     *  8.3415180d-01,  8.3655230d-01,  7.9039590d-01,  7.4855320d-01,
     *  7.2383030d-01,  7.0726850d-01,  6.9536850d-01,  6.9769030d-01,
     *  7.3286850d-01,  2.1551890d+00,  1.6566600d+00,
     *  8.3415180d-01,  1.3942260d+00,  7.6448150d-01,  7.1896530d-01,
     *  7.1925910d-01,  1.6412760d+00,  1.6733940d+00,  9.6108760d-01,
     *  7.7301580d-01,  1.7703380d+00,  1.6566600d+00 /
      data xlam2/
     *  2.0216150d+00,  2.0195960d+00,  2.0136340d+00,  1.9935080d+00,
     *  1.9661290d+00,  1.9607250d+00,  1.9688170d+00,  2.0482980d+00,
     *  2.1215750d+00,  2.1623680d+00,  2.0507600d+00,
     *  2.0216150d+00,  1.9949090d+00,  2.0706990d+00,  2.1256560d+00,
     *  2.2576140d+00,  2.3850670d+00,  2.4708780d+00,  2.4482290d+00,
     *  2.2170650d+00,  2.0571940d+00,  2.0507600d+00 /
      data r0 /
     *  3.8134270d+00,  3.8824980d+00,  3.9309310d+00,  3.8733090d+00,
     *  3.8090830d+00,  3.7973730d+00,  3.8033260d+00,  3.9391050d+00,
     *  2.8216570d+00,  3.1157650d+00,  3.1894430d+00,
     *  3.8134270d+00,  3.9445850d+00,  3.9429750d+00,  3.7645640d+00,
     *  3.4670520d+00,  3.4330330d+00,  3.3721860d+00,  3.5290920d+00,
     *  3.0846890d+00,  3.6840700d+00,  3.1894430d+00 /
      data c1 /
     * -1.5337840d+04, -1.4781490d+04, -1.1291030d+04, -8.2402870d+03,
     * -5.8078090d+03, -5.1548940d+03, -4.5375600d+03, -5.6384730d+03,
     * -9.9391440d+03, -5.4486500d+06, -1.3765500d+07,
     * -1.5337840d+04, -5.2244120d+05, -8.4470820d+03, -5.6583110d+03,
     * -3.1215130d+03, -2.5114390d+05, -8.9747930d+05, -1.0270760d+04,
     * -9.7752400d+03, -4.3667750d+07, -1.3765500d+07 /
      data c2/
     * -2.0204580d+07, -2.3924940d+07, -2.6047680d+07, -2.0237690d+07,
     * -1.1297860d+07, -7.2559300d+06, -4.9401140d+06, -1.1761180d+07,
     * -5.3292430d+07, -1.5990940d+08, -7.1468870d+07,
     * -2.0204580d+07, -2.6972790d+07, -4.4253780d+07, -4.6427270d+07,
     * -4.3504530d+07, -6.9556960d+07, -1.0156260d+08, -1.5983540d+08,
     * -1.0833420d+08, -7.3120690d+07, -7.1468870d+07 /
      data c3/
     *  2.3471920d+07,  2.3657030d+07,  2.2258200d+07,  1.7765330d+07,
     *  1.2376670d+07,  1.0316850d+07,  9.4428070d+06,  1.5325570d+07,
     *  4.2364170d+07,  1.0769020d+08,  8.8954890d+07,
     *  2.3471920d+07,  2.3226170d+07,  2.6351970d+07,  2.3146810d+07,
     *  2.0659610d+07,  3.0879900d+07,  4.6075670d+07,  6.9639620d+07,
     *  6.1803210d+07,  9.3372390d+07,  8.8954890d+07 /
       data clr/
     * -3.9696770d+06, -3.9877570d+06, -3.3100680d+06, -2.6855830d+06,
     * -2.5623320d+06, -2.3367800d+06, -2.3216730d+06, -2.0898290d+06,
     * -2.2441540d+06, -7.0960400d+06, -7.4588120d+06,
     * -3.9696770d+06, -5.8114580d+06, -3.1735720d+06, -2.5450750d+06,
     * -3.1654150d+06, -4.3745160d+06, -4.3923340d+06, -4.3510600d+06,
     * -3.3539990d+06, -6.8638170d+06, -7.4588120d+06 /
*
      npnt2 = npoint - 2
      ntots = npoint * npoint
      ntotd = npnt2 * npnt2

* rshift is factor for shift of atan midpoint in potential correction
* rrshift is inward shift of ab initio potential
* xfact is factor for additional potential correction

*     rshift = 0.5d0
*     rrshift = 0.4d0
*     xfact = 0.6667d0

* unmodified potential

      rshift  = zero
      rrshift = zero
      xfact   = zero

*     determine A' and A" potentials at angles
*     first shift radius out

      r = rr + rrshift
      rm3 = one / ( r**3 )
      rm6 = rm3 * rm3
      do 100 i = 1, npoint
         vap( i ) = c1( i ) * exp( -xlam1( i ) * r ) +
     :            ( c2( i ) + c3( i ) * r) * exp( - xlam2( i ) * r )
     :              + ( tanh( alph * ( r - r0( i ) ) ) + one )
     :              * clr( i ) * rm6

*        determine switching function (turned off by smc)

*        fact=xfact*half*(tanh(alph*(r-re(i)+rshift))+1)
*        vap(i)=vap(i)+fact*vap(i)

         j = i + npoint
         va2p( i ) = c1( j ) * exp( -xlam1( j ) * r ) +
     :            ( c2( j ) + c3( j ) * r ) * exp( -xlam2( j ) * r )
     :            + ( tanh( alph * ( r - r0( j ) ) ) + 1 )
     :            * clr( j ) * rm6

*        determine switching function (turned off by smc )

*        fact=xfact*half*(tanh(alph*(r-re(j)+rshift))+1)
*        va2p(i)=va2p(i)+fact*va2p(i)

         vsum( i ) = half * ( va2p( i ) + vap( i ) )

*        don't compute vdif for colinear geometries

         if( i .ne. 1 .and. i .ne. npoint ) then
            vdif( i - 1 ) = half * ( va2p( i ) - vap( i ) )

*        for long range damp out difference potential

            if( r .gt. rmax ) then
               damp = - ( tanh( one * ( r - rmax ) ) - one )
               vdif( i - 1 ) = vdif( i - 1 ) * damp
            end if

         end if

100   continue

* solve simultaneous equations for solutions
* first for vsigma

      tol = 1.0d-10
      call dcopy( ntots, d0, 1, aa, 1 )

      call dqrank(aa,npoint,npoint,npoint,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,npoint,npoint,npoint,kr,vsum,xsum,rsd,kpvt,qraux)

*     convert to hartree (potential is microhartrees)

      conv = 1.0d-6
      call dscal( npoint, conv, xsum, 1 )
      vv0 = xsum( 1 )
      call dcopy( npoint - 1, xsum( 2 ), 1, vvl, 1 )

*     solve simultaneous equations for solutions
*     for vdif

      call dcopy( ntotd, d2, 1, aa, 1 )
      call dqrank( aa, npnt2, npnt2, npnt2, tol, kr, kpvt, qraux, work )
      call dqrlss( aa, npnt2, npnt2, npnt2, kr, vdif, xdif, rsd, kpvt,
     *             qraux )
* convert to hartree and multiply by minus one so that
* Vdiff is defined as VA"-VA'
      call dscal( npnt2, onemin*conv, xdif, 1 )
      call dcopy( npnt2, xdif, 1, vvl( npoint ), 1 )
      return
      end
