* System:  CH(X 2Pi)+He, original ab initio CEPA PES's
* Reference: 
* Albert F. Wagner, Thom H. Dunning, and Randall A. Kok,
* J. Chem. Phys. 100, 1326 (1994) ; 
* M. H. Alexander, W. Kearney, and A. F. Wagner,
* J. Chem. Phys. 100, 1338 (1994).
*
*  Note:  this pot routine requires a data file to be in hibxx/bin/progs/potdata:
*         chhe_abin.dat
*
      include "common/syusr"
      subroutine driver
      implicit double precision (a-h,o-z)
      character*60 filnam, filnm1
      character*80 delete
      logical lpar, lpar1, batch, readpt
      common /colpar/ lpar(3), batch,lpar1(10),readpt
      common /covvl/ vvl(11)
      econv=219474.6d0
      readpt=.true.
      batch=.true.
      filnam='chhe_abin.dat'
      open(unit=9,file='fort.9',status='unknown',
     :      access='sequential')
      delete='rm fort.9'
      filnm1 = 'potdata/'//filnam
      call loapot(10,filnam)
1      print *, ' r (bohr)'
      read (5, *, end=99) r
      call pot(vv0,r)
      if (r.le.0) goto 99
      write (6, 100) vv0,vvl
100   format(' vsum',/,8(1pe16.8),/,
     :    '  vdif',/,4e16.8)
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
        write (12,110) rr,vv0*econv,(econv*vvl(j),j=1,7)
110     format(f7.2,19(1pe16.8))
        write (13,110) rr,vv0*econv,(econv*vvl(j),j=8,11)
        rr = rr + dr
      enddo
      close(12)
      close(13)
      call system(delete)
      end
      include "common/bausr"
      include "common/ground"
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
*  subroutine to calculate the r-dependent coefficients in the
*  CH-He 2Pi potential

*  authors:  al wagner, millard alexander, and bill kearney
*  latest revision date: 28-may-1997 by mha
* ----------------------------------------------------------------------

*  on return:
*  vv0 is v-average for lambda=0
*  the coefficients which multiply each reduced rotation matrix element
*  for each angular term in the coupling potential
*  [ vvl(i) for i = 1, nlam=nsum+ndif (+1 for variable spin-orbit)]
*         are returned in common block /covvl/
*  vvl(i=1,nsum) contains the potential expansion coefficients for
*  the average pi potential for lam = 1, nsum
*  vvl(i=nsum, nsum+ndif-1) contains the potential expansion coefficients for
*  the difference pi potential for lam = 2, 1 + ndif
*  vvl(nsum+ndif+1) contains the spin-orbit contribution

*  variable in common block /conlam/
*    nlammx:    the maximum number of angular coupling terms
*    nlam:      the number of angular coupling terms actually used
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx
*    vvl        array to store r-dependence of each angular term in the

*  variable in common block /covvl/
*               potential

*  variables in common block /cosysi/

*    nscode:   total number of system dependent variables
*    isicod:   total number of integer system dependent variables
*    nterm:    number of different diabatic potentials, this is 2 here
*    jmax:     the maximum rotational angular momenta for each 2pi channel
*              in each spin-orbit manifold with convention
*              omega .le. j .le. jmax+0.5
*    igupi:    permutation inversion symmetry of 2pi electronic state
*              igu=1 for gerade states, igu=-1 for ungerade states
*              for heteronuclear molecules igu should be +1
*    isa:      s/a label for molecular states. If ihomo=.true. then only
*              s states will be included if isa=1 and only a states if
*              isa=-1
*    nparpi:   number of 2pi symmetry doublets included (npar=2 will ensure
*              both lambda doublets)


      implicit double precision (a-h,o-z)
      logical existf, readpt, lpar, lpar1, batch
      character*(*) filnam
      character*75 title, filnm1
      include "common/parbas"
      include "common/parpot"

*  nvmax is maximum number of variables
*  npmax is maximum number of anisotropic terms
*  ntmax is maximum number of terms (nterm) in wagner expansion
*   this is the sum of the number of terms for sum potential and dif. potentia
      parameter (nvmax=3,  npmax=15,  ntmax=90)
      parameter (ksum = 7, kdif = 4)
      dimension ynow(npmax), nlegp(npmax), ncor(ntmax),
     :          mineold(nvmax), maxeold(nvmax)
      dimension nexp(nvmax, ntmax), ispf(nvmax), shif(nvmax),
     :          xprod(nvmax,-30:20), iatm(5,nvmax), ipm(nvmax),
     :          minexp(nvmax), maxexp(nvmax), fk(ntmax), fke(ntmax),
     :          f(nvmax), x(nvmax), scale(nvmax)
      dimension iconv(2:1 + kdif)
      common /cochhe1/ynow
      common /cochhe2/nlegp
      common /cochhe3/ncor
      common /cochhe4/mineold
      common /cochhe5/maxeold
      common /cochhe6/nexp
      common /cochhe7/ispf
      common /cochhe8/shif
      common /cochhe9/xprod
      common /cochhe10/iatm
      common /cochhe11/ipm
      common /cochhe12/minexp
      common /cochhe13/maxexp
      common /cochhe14/fk
      common /cochhe15/fke
      common /cochhe16/f
      common /cochhe17/x
      common /cochhe19/scale

      common /colpar/ lpar(3), batch,lpar1(10),readpt
      common /conlam/ nlam, nlammx, lamnum(1)
      common /covvl/vvl(1)
      common /cosysr/ isrcod, junkr, brot, aso, p, q, rsm,asod,bso
      common /cosysi/ nscode, isicod, nterm, jmax, igupi, isa, nparpi


*  nsum is maximum number of legendre polynomials in current expansion
*  of ansisotropy of sum potential (0 .le. lambda .le. nsum)
*  ndif is maximum number of legendre polynomials in current expansion of
*  difference potential ( 2 .le. lambda .le. ndif+1)
      zero = 0.d0
      nsum = ksum
      ndif = kdif
      potnam='WAGNER et al CI CHHE PES'

      lammin(1)=1
      lammax(1)=7
      lammin(2)=2
      lammax(2)=5
      nlam=7+5-2+1
      mproj(1)=0
      mproj(2)=2
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
* return at this point if readpt is false
      if (.not.readpt) return

*  iconv is sqrt[(l+2)!/(l-2)!] to convert from expansion in legendre
*  polynomials of order 2 to reduced rotation matrix elements

      do  10 i = 2, 1 + ndif
        iconv(i) = (i+2)*(i+1)*i*(i-1)
10    continue

c input is handled here
c first read in the .5(a''+a') surface parameters
      filnm1 = 'potdata/'//filnam
      inquire(file=filnm1,exist=existf)
      if (.not. existf) then
        write (6, 11) filnam
11      format
     :  (/'   *** INPUT FILE ',(a),' NOT IN LOCAL DIRECTORY; ABORT ***')
        if (batch) call exit
        return
      endif
      open(unit=iunit,file=filnm1,status='old',
     :      access='sequential')
      read(iunit,180)title
180   format(a)
      read (iunit, *)npts,nterms,nvar,natoms,nspf
      if(nspf.ne.0) read (iunit, *) (ispf(i), i=1,nvar)
      read (iunit, *) (minexp(i),i=1,nvar)
      read (iunit, *) (maxexp(i),i=1,nvar)
      do 190 i=1,nterms
        read (iunit, *) fk(i),fke(i),(nexp(j,i),j=1,nvar)
190   continue
      read (iunit, *)      ((iatm(i,j),i=1,5),j=1,nvar)
      read (iunit, *)      (ipm(j),           j=1,nvar)
      read (iunit, *)      (shif(j),          j=1,nvar)
      read (iunit, *)      (scale(j),         j=1,nvar)
      read (iunit, *)      (f(j),             j=1,nvar)
      write(9,195) filnam,title
      write(6,195) filnam,title
195   format(/,' SUM POTENTIAL SUCCESSFULLY LOADED FROM FILE:',
     :      2x,(a),/,'  TITLE: ',(a))

      do 200 i = 1,nvar
        mineold(i) = minexp(i)
        maxeold(i) = maxexp(i)
200   continue
c now read in the .5(a''-a') surface

      nterm0 = nterms+1
      read(iunit,180)title
      read (iunit, *)npts,nterms,nvar,natoms,nspf
      nterms = nterm0-1+nterms
      if(nspf.ne.0) read (iunit, *) (ispf(i), i=1,nvar)
      read (iunit, *) (minexp(i),i=1,nvar)
      read (iunit, *) (maxexp(i),i=1,nvar)
      do 240 i=nterm0,nterms
        read (iunit, *) fk(i),fke(i),(nexp(j,i),j=1,nvar)
240   continue
      read (iunit, *)      ((iatm(i,j),i=1,5),j=1,nvar)
      read (iunit, *)      (ipm(j),           j=1,nvar)
      read (iunit, *)      (shif(j),          j=1,nvar)
      read (iunit, *)      (scale(j),         j=1,nvar)
      read (iunit, *)      (f(j),             j=1,nvar)
      close (iunit)
      write(9,250) filnam,title
      write(6,250) filnam,title
250   format(' DIFFERENCE POTENTIAL SUCCESSFULLY LOADED FROM FILE:',
     :      2x,(a),/,'  TITLE: ',(a))
      do 260 i = 1,nvar
        minexp(i) = min0(minexp(i),mineold(i))
        maxexp(i) = max0(maxexp(i),maxeold(i))
260   continue

*  nlegp is a unique index for each legendre term, for congruence with
*  wagner's data files

*  for the sum potential (regular legendre polynomials)
*     nlegp = l * (l+1) / 2
*  for the difference potential (associated legendre polynomials of order2)
*     nlegp = l * (l+1) / 2 + 2

      do  280  i = 1, nsum + 1
        ll = i - 1
        nlegp(i) = ( ll * (ll + 1) ) / 2
280   continue

      do 300  i = 1, ndif
        ll = i + 1
        nlegp(i + nsum + 1) = ( ll * (ll + 1) ) / 2 + 2
300   continue

      nleg = nsum + ndif + 1

c order this in the most compact way for multiplication of the
c legendre polynomials

      do 350 k = 1,nterms
        nnow = nexp(3,k)
        do 340 j = 1,nleg
          if(nlegp(j).eq.nnow) ncor(k) = j
340     continue
350   continue
c
c assume values for vibrational coordinate x(1) and center of mass angle
c   x(3).  these must be provided by input.
c assume the distance of the atom from the center of mass of the
c   vibrator is x(2).  the value provided here will by overriden by
c   the r in the input argument.
c
*  r(ch) and dummy values of Rch-he and theta
*  r(ch) is set here to the outer turning point of the v=0 state

      f(1) = 2.110755d0
      f(2) = 3.d0
      f(3) = 90.d0
      do 380 i=1,nvar
        x(i) = (f(i) / scale(i) ) - shif(i)
380   continue

      return
      entry pot(vv0, r)


c assume values for vibrational coordinate x(1) and center of mass angle
c   x(3).  these must be provided by input via f(i).
c assume the distance of the atom from the center of mass of the
c   vibrator is x(2).  this is provided by the argument r which must
c   be provided in the final scaled and shifted form.

      x(2)=r
      do 450 i=1,2
      go to (410, 430, 450, 450, 450, 450), ispf(i)+1
      stop 'error in eval - loop 450'
410   xtemp=x(i)
      xprod(i,minexp(i))=xtemp**minexp(i)
      do 420 j=minexp(i)+1, maxexp(i)
        xprod(i,j)=xprod(i, j-1)*xtemp
420   continue
      go to 450

430   xtemp= x(i)/(x(i)+shif(i))
      xprod(i,minexp(i))=xtemp**minexp(i)
      do 440 j=minexp(i)+1, maxexp(i)
        xprod(i,j)=xprod(i,j-1)*xtemp
440   continue
450   continue
c
c compute the coefficients for the legendre expansion
c
      do 470 j = 1, nleg
        ynow(j) = zero
470   continue
      do 500 k=1, nterms
        jj = ncor(k)
        ykl = 1.0
        do 480 j=1, 2
          njk=nexp(j,k)
          if(njk.ne.0) ykl = ykl*xprod(j,njk)
480     continue
        ynow(jj) = ynow(jj) + fk(k)*ykl
500   continue
      vv0 = ynow(1)
      do 550 i = 2, nleg
        vvl(i-1) = ynow(i)
550   continue


      do 650 i = nsum + 1, nsum + ndif
         lam = i - nsum +1
           vvl(i) = vvl(i) * sqrt( float(iconv(lam)) )
650   continue
      return
      end
