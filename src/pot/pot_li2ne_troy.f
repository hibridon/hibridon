****************************************************************************
* pot_li2ne_low-v.f
*
* This is a special version of Millard H. Alexander's Li2(A)-Ne potential
* module (pot_li2ne.f).  I have made a few small modifications to the
* original source code to enable it to handle v=0-4.  (MHA's original
* version was set up for calculations in the range v=7-11.)
*
* All of my modifications are flagged with comments containing [TS].  If you
* modify this code further, please document your changes.
*
* Troy Stephens
* January 25, 1997
****************************************************************************
*
*system: Li2(A)-Ne using alexander-werner potential
*reference: M. H. Alexander and H.-J. Werner, J. Chem. Phys. 95, 6524 (1991).
*
*  Note:  this pot routine requires a data file to be in hibxx/bin/progs/potdata:
*         LI2NE.MHA
*
      include "common/syusr"
      subroutine driver
      implicit double precision (a-h,o-z)
      character*9 filnam
      character*80 delete
      logical lpar,readpt,lpar1,batch
      common /cosysi/ junk(3), nvibmn, nvibmx
      common /covvl/ vvl(30)
      common /colpar/ lpar(3), batch,lpar1(10),readpt
      include "common/parpot"
      readpt=.true.
      batch=.true.
      potnam='ALEXANDER-WERNER LI2(A)-NE'
      filnam='LI2NE.MHA'
      open(unit=9,file='fort.9',status='unknown',
     :      access='sequential')
      delete='rm fort.9'
      print *, potnam
      call loapot(10,filnam)
1      print *, ' enter r in bohr or control-D to quit'
      nvibmn=8
      nvibmx=9
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) (vvl(i), i=1,10)
100   format(' <v=8 | V | v=8>',/,10(1pe16.8))
      write (6, 101) (vvl(i), i=11,20)
101   format(' <v=9 | V | v=8>',/,10(1pe16.8))
      write (6, 102) (vvl(i), i=11,20)
102   format(' <v=9 | V | v=9>',/,10(1pe16.8))
      goto 1
99    call system(delete)
      end
      include "common/bausr"
      include "common/ground"
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
*
*   subroutine to load potential for li2-ne
*   author: m. alexander
*   current revision date: 12-jun-1991
*
*   on entry: iunit -> unit assigned to file containing the potential
*             filnam -> filename
*
*  variables in common block /conlam/
*    nlammx:    the maximum number of angular coupling terms used
*               here this should be nvib*(nvib+1)*5
*    nlam:      the number of angular coupling terms used for each <v | v'> bl
*               this should be 10
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx (this is not used here)
*
*  variable in common block /covvl/
*    vvl        array to store r-dependence of each angular term in the
*               potential
* variables in common block /cosysi/
*    nscode:    total number of system variables (not used here)
*    isicod:    total number of real variables
*    nterm:     number of different types of legendre terms (not used here)
*    nvibmn:    minimum vibrational quantum number
*    nvibmn:    maximum vibrational quantum number
*               the vibrational levels used range from nvibmn to nvibmx
* variable in common block /copt1/
*    inbv       parameter which should be 1 for first entry into bvalu
*    nlams      total number of vvl terms used (should be 30)
*    npoint     number of values of r used in spline interpolation
*    inclog     increment used in accumulating expansion coefficients for
*               logarthmic spline interpolation
*    inclin     increment used in accumulating expansion coefficients for
*               linear spline interpolation
* variable in common block /copt2/
*    rr         array of legth npoint containing values of r used in spline in
* variable in common block /copt3/ and /copt4/
*    aexp, bexp expansion coefficients for exponential extrapolation at short
* variable in common block /copt5/
*    c6         c-6 coefficient for extrapolation at long range
* variable in common block /copt6/ and /copt7/
*    aexpr
*    bexpr      expansion coefficients for exponential extrapolation at long r
* variable in common block /copt8/
*    logcut     index of last point used in logarithmic interpolation
* variable in common block /copt9/
*    cutlog     r value of last point used in logarithmic interpolation
* variable in common block /copt10/
*    linfirst   index of first point used in linear interpolation
* variable in common block /copt11/
*    lincut     index of last point used in linear interpolation
* variable in common block /copt12/
*    cutlin     r value of last point used in linear interpolation
* variable in common block /copt13/
*    tmat       knot array for logarithmic interpolation
* variable in common block /copt14/
*    ttmat      knot array for linear interpolation
* variable in common block /copt15/
*    coef       b-spline coefficient array for logarithmic interpolation
* variable in common block /copt16/
*    ccoef      b-spline coefficient array for linear interpolation
* variable in common block /copt17/
*    work       scratch array used in spline interpolations
* variable in common block /copt18/
*    ncoef      number of coefficients in logarithmic interpolation
* variable in common block /copt19/
*    nncoef     number of coefficients in linear interpolation
* variable in common block /copt20/
*    nlin       number of points in linear interpolation
*  bint4 and bvalu are derived from the bspline package,
*  part of cmlib (adapted for the vax 11/785).
*
* --------------------------------------------------------------------------
      implicit double precision(a-h,o-z)
* mlams is the maximum number of legendre terms in input
* melts is the maximum number of matrix elements to be calculated
* mpoint is the maximum number of values of r for each legendre term
* maxlog is the maximum number of points for logarithmic spline fit
* maxlin is the maximum number of points for the linear spline fit
* all these are specific to this potential and should not be changed
* without consulting mha
      parameter (mlams=30, maxlog=18, maxlin=11, mpoint=21)
      parameter (mt=mlams*(maxlog+6), mc=mlams*(maxlog+2),
     :           mtt=mlams*(maxlin+6), mcc=mlams*(maxlin+2),
     :           mw=5*(maxlog+2) )
      character*(*) filnam
      character*80 potlab, text, filnm1
      logical lpar,readpt,lpar1,batch,existf
      dimension tmat(mt), ttmat(mtt),  cmat(mc), ccmat(mcc),
     :          logcut(mlams), cutlog(mlams), linfirst(mlams),
     :          lincut(mlams), cutlin(mlams), aexp(mlams), bexp(mlams),
     :          c6(mlams), aexpr(mlams), bexpr(mlams), rr(mpoint),
     :          nlin(mlams), ncoef(mlams), nncoef(mlams),
     :          work(mw), v(mpoint), vlog(maxlog)
      include "common/parpot"
      include "common/parbas"
* comdeck parbas
*      parameter (maxtrm=4,maxvib=10,maxvb2=maxvib**2)
*  variables in common block /cobspt/
*    lammin:   array containing minimum value of lambda for each term
*    lammax:   array containing maximum value of lambda for each term
*    mproj:    array containing the order of the reduced rotation matrix
*              elements for each term.  lammin can not be less than mproj.
*              for homonuclear molecules, the allowed values of lambda for
*              each term range from lammin to lammax in steps of 2

*      common /cobspt/ lammin(maxtrm), lammax(maxtrm), mproj(maxtrm)
*      common /cobsp2/ ntv(maxtrm),ivcol(maxvb2,maxtrm),
*     :                ivrow(maxvb2,maxtrm)

      dimension vvl(1)
      common /cotrms/ ntrm,lmmin(10),lmmax(10),mld(10),nv(10)
      common /colpar/ lpar(3), batch,lpar1(10),readpt
      common /covvl/ vvl
      common /cosysi/ nscode, isicod, nterm, nvibmn, nvibmx
      common /conlam/ nlam, nlammx, lamnum
      common /copt1/ inbv, nlams, npoint, inclog, inclin
      common /coskip/ nskip,iskip
      common /copt2/ rr
      common /copt3/ aexp
      common /copt4/ bexp
      common /copt5/ c6
      common /copt6/ aexpr
      common /copt7/ bexpr
      common /copt8/ logcut
      common /copt9/ cutlog
      common /copt10/ linfirst
      common /copt11/ lincut
      common /copt12/ cutlin
      common /copt13/ tmat
      common /copt14/ ttmat
      common /copt15/ cmat
      common /copt16/ ccmat
      common /copt17/ work
      common /copt18/ ncoef
      common /copt19/ nncoef
      common /copt20/ nlin
      data zero /0.d0/
      potnam='ALEXANDER-WERNER LI2(A)+NE'
      ntrm=1
      iskip=2
      lmmin(1)=0
      lmmax(1)=18
      lammin(1)=0
      lammax(1)=18
      mproj(1)=0
      mld(1)=0
      nvib=nvibmx-nvibmn+1
      nv(1)=nvib*(nvib+1)/2
      ntv(1)=nv(1)
      index=0
      do 150 jvrow=nvibmn, nvibmx
        do 140 jvcol=nvibmn,jvrow
          index=index+1
          ivrow(index,1)=jvrow
          ivcol(index,1)=jvcol
140     continue
150   continue
* return at this point if readpt is false
      if (.not.readpt) return
*
      nlams=mlams
      npoint=mpoint
      inclog=maxlog
      inclin=maxlin
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
*      open(unit=9,file='file9.out',access='sequential',status='unknown')
      read(iunit,'(a)',err=800,iostat=ierr) potlab
      read(iunit,'(a)',err=800,iostat=ierr) text
      read(iunit,*,err=800,iostat=ierr) (rr(i), i=1, npoint)
* define spline initial conditions
      ibcl = 2
      ibcr = 2
      fbcl = zero
      fbcr = zero
      kntopt = 1
      inbv=1
* initialize running indices for coefficient arrays
      indt=1
      indtt=1
      indc=1
      indcc=1
       do 200 it = 1, nlams
* read in potential
        read(iunit,'(a)',err=800,iostat=ierr) text
        read(iunit,*,err=800,iostat=ierr) (v(i), i=1,7)
        read(iunit,*,err=800,iostat=ierr) (v(i), i=8,14)
        read(iunit,*,err=800,iostat=ierr) (v(i), i=15,21)
        read(iunit,*,err=800,iostat=ierr) aexp(it), bexp(it),
     :       aexpr(it), bexpr(it), c6(it)
        read(iunit,*,err=800,iostat=ierr) logcut(it), cutlog(it),
     :               linfirst(it), lincut(it), cutlin(it)
        nlin(it)=lincut(it)-linfirst(it)+1
        if (logcut(it) .gt. 0) then
* here under normal conditions
* convert necessary points to logarithms for logarithmic spline
          do 100 i=1,logcut(it)
            vlog(i)=log(v(i))
100       continue
*         write (6, *) 'before first bint4, it=',it
* do spline fit for logarithmic spline
          call bint4 (rr, vlog, logcut(it), ibcl, ibcr, fbcl,
     :                fbcr,kntopt, tmat(indt), cmat(indc),
     :                ncoef(it),kord, work)
* do spline fit for linear spline (if needed)
          nncoef(it)=0
          if (nlin(it) .ge. 2) then
            call bint4 (rr(linfirst(it)), v(linfirst(it)),
     :                  nlin(it), ibcl, ibcr, fbcl, fbcr, kntopt,
     :                  ttmat(indtt), ccmat(indcc), nncoef(it), kord,
     :                  work)
         endif
        else
* here if this particular term in potential is nul
          ncoef(it)=0
          nncoef(it)=0
        endif
* increment indices
        indt=indt+inclog+6
        indtt=indtt+inclin+6
        indc=indc+inclog+2
        indcc=indcc+inclin+2
200   continue
      close(iunit)
      write(6,205) filnam,potlab(:72)
      write(9,205) filnam,potlab(:72)
205   format(/,' POTENTIAL SUCCESSFULLY LOADED FROM FILE:',
     :      7x,(a),/,1x,'LABEL: ',(a))
      return
800   write(6,850) ierr
850   format(' *** READ ERROR IN LOAPOT, IOSTAT: ',i4,', ABORT ***')
      write (6, 855) text
855   format('     last potential label read was:',
     :       10x,(a) )
      return
      end
* ----------------------------------------------------------------------
              subroutine pot (vv0, rx)
***** calculates the coeff. vvl of each rotation matrix element,
* in atomic units (hartree).
* ----------------------------------------------------------------------
*  subroutine to calculate the r-dependent coefficients in the
*  li2-ne potential
*  author:  millard alexander
*  latest revision date:  18/12/89
* ----------------------------------------------------------------------
*  on return:
*  vv0 is v-average for lambda=0 (this is set equal to zero here!)
*  the coefficients for each angular term in the coupling potential
*  [ vvl(i) for i = 1, nlammx ] are returned in common block /covvl/
*  the vibrational quantum numbers range from nvibmn to nvibmn+nvib-1
*  the 10 angular coupling coefficients for each <v' | v> term are stored
*  in the array vvl in row-packed, lower triangular form, e.g.
*  <nvibmn | nvibmn> stored in vvl(1:10)
*  <nvibmn+1 | nvibmn> stored in vvl(11:20)
*  <nvibmn+1 | nvibmn+1> stored in vvl(21:30)
*  <nvibmn+2 | nvibmn> stored in vvl(31:40)
*  <nvibmn+2 | nvibmn+1> stored in vvl(41:50)
*  <nvibmn+2 | nvibmn+2> stored in vvl(51:60)
*  ...
*  variables in common block /conlam/
*    nlammx:    the maximum number of angular coupling terms used
*               here this should be at least nvib*(nvib+1)*5
*    nlam:      the number of angular coupling terms used for each <v | v'> bl
*               this should be 10
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx (this is not used here)
*
*  variable in common block /covvl/
*    vvl        array to store r-dependence of each angular term in the
*               potential
* variables in common block /cosysi/
*    nscode:    total number of system variables (not used here)
*    isicod:    total number of real variables
*    nterm:     number of different types of legendre terms (not used here)
*    nvibmn:    minimum vibrational quantum number
*    nvibmn:    maximum vibrational quantum number
*               the vibrational levels used range from nvibmn to nvibmx
* variable in common block /copt1/
*    inbv       parameter which should be 1 for first entry into bvalu
*    nlams      total number of vvl terms used
*    npoint     number of values of r used in spline interpolation
*    inclog     increment used in accumulating expansion coefficients for
*               logarthmic spline interpolation
*    inclin     increment used in accumulating expansion coefficients for
*               linear spline interpolation
* variable in common block /copt2/
*    rr         array of legth npoint containing values of r used in spline in
* variable in common block /copt3/ and /copt4/
*    aexp, bexp expansion coefficients for exponential extrapolation at short
* variable in common block /copt5/
*    c6         c-6 coefficient for extrapolation at long range
* variable in common block /copt6/ and /copt7/
*    aexpr
*    bexpr      expansion coefficients for exponential extrapolation at long r
* variable in common block /copt8/
*    logcut     index of last point used in logarithmic interpolation
* variable in common block /copt9/
*    cutlog     r value of last point used in logarithmic interpolation
* variable in common block /copt10/
*    linfirst   index of first point used in linear interpolation
* variable in common block /copt11/
*    lincut     index of last point used in linear interpolation
* variable in common block /copt12/
*    cutlin     r value of last point used in linear interpolation
* variable in common block /copt13/
*    tmat       knot array for logarithmic interpolation
* variable in common block /copt14/
*    ttmat      knot array for linear interpolation
* variable in common block /copt15/
*    coef       b-spline coefficient array for logarithmic interpolation
* variable in common block /copt16/
*    ccoef      b-spline coefficient array for linear interpolation
* variable in common block /copt17/
*    work       scratch array used in spline interpolations
* variable in common block /copt18/
*    ncoef      number of coefficients in logarithmic interpolation
* variable in common block /copt19/
*    nncoef     number of coefficients in linear interpolation
* variable in common block /copt20/
*    nlin       number of points in linear interpolation
*  bint4 and bvalu are derived from the bspline package,
*  part of cmlib (adapted for the vax 11/785).
* -----------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      dimension tmat(1), ttmat(1),  cmat(1), ccmat(1),
     :          logcut(1), cutlog(1), linfirst(1),
     :          lincut(1), cutlin(1), aexp(1), bexp(1),
     :          c6(1), aexpr(1), bexpr(1), rr(1),
     :          nlin(1), ncoef(1), nncoef(1),
     :          work(1)
      dimension rm66(5)
      dimension rmre(5,5), rmre2(5,5)
      common /cosysi/ nscode, isicod, nterm, nvibmn, nvibmx
      common /conlam/ nlam, nlammx, lamnum
      common /covvl/ vvl(1)
      common /copt1/ inbv, nlams, npoint, inclog, inclin
      common /copt2/ rr
      common /copt3/ aexp
      common /copt4/ bexp
      common /copt5/ c6
      common /copt6/ aexpr
      common /copt7/ bexpr
      common /copt8/ logcut
      common /copt9/ cutlog
      common /copt10/ linfirst
      common /copt11/ lincut
      common /copt12/ cutlin
      common /copt13/ tmat
      common /copt14/ ttmat
      common /copt15/ cmat
      common /copt16/ ccmat
      common /copt17/ work
      common /copt18/ ncoef
      common /copt19/ nncoef
      common /copt20/ nlin
      data conv /219474.6d0/
      data one /1.d0/
****************************************************************************
* The three tables below (rm66, rmre, and rmre2) are the main site of the
* low-v modifications.  I have calculated new numbers for v=0-4 to replace
* the old ones. [TS]
****************************************************************************
* corrected for bohr/angstrom confusion by mha 3/1/97

*  here are the diagonal matrix elements of (r-6.6) for v=0-4
*      data rm66 /-3.4741,-3.4385,-3.4026,-3.3663,-3.3298/
      data rm66 /-6.9290d-01, -6.2563d-01,-5.5779d-01,-4.8919d-01,
     :          -4.2022d-01/
*  here are the matrix elements of (r-re) for v=0-4
      data rmre/
     :  3.3796d-02, 2.5984d-01,-1.5685d-02, 1.7008d-03,-1.8897d-04,
     :  2.5984d-01, 1.0107d-01, 3.6850d-01,-2.7212d-02, 3.4015d-03,
     : -1.5685d-02, 3.6850d-01, 1.6891d-01, 4.5240d-01,-3.8550d-02,
     :  1.7008d-03,-2.7212d-02, 4.5240d-01, 2.3751d-01, 5.2364d-01,
     : -1.8897d-04, 3.4015d-03,-3.8550d-02, 5.2364d-01, 3.0648d-01/
*  here are the matrix elements of (r-re)^2 for v=0-4
      data rmre2/
     : 6.8611d-02, 2.9222d-02, 9.3802d-02,-1.3582d-02, 2.8277d-03,
     : 2.9222d-02, 2.1427d-01, 8.3780d-02, 1.6038d-01,-2.6901d-02,
     : 9.3802d-02, 8.3780d-02, 3.7014d-01, 1.5362d-01, 2.2368d-01,
     : -1.3582d-02, 1.6038d-01, 1.5362d-01, 5.3926d-01, 2.3683d-01,
     : 2.8277d-03,-2.6901d-02, 2.2368d-01, 2.3683d-01, 7.1769d-01/
      izero=0
      kord=4
      if (rx .lt. rr(1)) then
* here for exponential integration
          do 100 it=1, nlams
             vvl(it)=0.0d0
             if (ncoef(it) .ne. 0.)
     :         vvl(it)=exp(aexp(it)+bexp(it)*rx)/conv
100       continue
          return
        endif
* initialize running indices for coefficient arrays
      indt=1
      indtt=1
      indc=1
      indcc=1
      rm6 = one/rx**6
      do 200 it = 1, nlams
        if (ncoef(it) .eq. 0) then
          vv=0.0d0
        else
          if (rx .le. cutlog(it)) then
            vv= bvalu (tmat(indt), cmat(indc), ncoef(it), kord,
     :                  izero, rx, inbv, work)
            vv=exp(vv)
          else if (rx .le. cutlin(it)) then
	    if (nlin(it) .ge. 2) then
              vv= bvalu (ttmat(indtt), ccmat(indcc), nncoef(it),
     :                    kord, izero, rx, inbv, work)
            else
              vv=exp(aexpr(it)+bexpr(it)*rx)
            endif
          else if (rx .gt. cutlin(it)) then
            if (nlin(it) .ge. 2) then
              vv=0.
              if (c6(it) .lt. 0.) then
                vv=c6(it)*rm6
              else
                if (aexpr(it) .gt. 1.e-6)
     :              vv=exp(aexpr(it)+bexpr(it)*rx)
              endif
            else
              if (aexpr(it) .gt. 1.e-6)
     :            vv=exp(aexpr(it)+bexpr(it)*rx)
            endif
          endif
        endif
        vvl(it)=vv/conv
* increment indices
        indt=indt+inclog+6
        indtt=indtt+inclin+6
        indc=indc+inclog+2
        indcc=indcc+inclin+2
200   continue
* now calculate terms needed for vibrational averaging
* save v66 in work(1:10)
* save dv/dr |r=6.6 in work(11:20)
* save dv/dr |r=5.8733 in work(21:30)
* save d2v/dr2 |r=5.8733 in work(31:40)
      do 250 i=1, 10
        v66=vvl(i)
        v58=vvl(i+10)
        v50=vvl(i+20)
        work(i) = v66
        work(i+10) = (v66-v58)/0.72670
        work(i+20) = .730506*v66 -0.159982*v58
     :               -0.570523*v50
        work(i+30) = .888369*v66 -1.67346*v58
     :               +.785088*v50
250   continue
* calculate diagonal terms
      do 300 nvib=nvibmn, nvibmx
****************************************************************************
*  ind1 is now the displacement of the vibrational quantum number relative
*  to v=0.  (ind1=1 if v=0, 2 if v=1, etc.)  [TS]
****************************************************************************
*OLD COMMENT*  ind1 is displacement of vibrational quantum number relative
*OLD COMMENT*  to v=7 (ind1=1 if v=7, 2 if v=8, etc)
****************************************************************************
*        ind1 = nvib-6
        ind1 = nvib+1
*  ind2 is displacement of vibrational quantum number relative
*  to v=nvibmn
        ind2 = nvib-nvibmn + 1
*  kvib is the first element in the block corresponding to the diagonal elemen
*  (ind2,ind2)
        kvib=((ind2*(ind2+1))/2-1)*10+1
        do 290  i = 1, 10
          if (ind1 .le. 0 .or. ind1 .gt. 5) then
            write (6, *) 'ind1 out of range ',ind1
            stop
          endif
          if (kvib+i-1 .le. 0 .or. kvib+i-1 .gt. 150) then
            write (6, *) 'kvib+i-1 out of range ',kvib+i-1
            stop
          endif
          vvl(kvib+i-1)=work(i)+rm66(ind1)*work(i+10)
290     continue
300   continue
* now off diagonal terms
* only do this if more than two vibrational levels present
      if (nvibmx .gt. nvibmn) then
        do 350 nvib1=nvibmn+1, nvibmx
****************************************************************************
*  ind1 is now the displacement of the vibrational quantum number relative
*  to v=0.  (ind1=1 if v=0, 2 if v=1, etc.)  [TS]
****************************************************************************
*OLD COMMENT*  ind1 is displacement of vibrational quantum number relative
*OLD COMMENT*  to v=7 (ind1=1 if v=7, 2 if v=8, etc)
*          ind1 = nvib1-6
          ind1 = nvib1+1
*  jnd1 is displacement of vibrational quantum number relative
*  to v=nvibmn
        jnd1 = nvib1-nvibmn + 1
          do 320 nvib2=nvibmn, nvib1-1
****************************************************************************
*  ind2 is now the displacement of the vibrational quantum number relative
*  to v=0.  (ind2=1 if v=0, 2 if v=1, etc.)  [TS]
****************************************************************************
*OLD COMMENT*  ind2 is displacement of vibrational quantum number relative
*OLD COMMENT*  to v=7 (ind1=1 if v=7, 2 if v=8, etc)
*            ind2 = nvib2-6
            ind2 = nvib2+1
*  jnd2 is displacement of vibrational quantum number relative
*  to v=nvibmn
            jnd2 = nvib2-nvibmn + 1
*  indst is the running index of the particular block of terms
            indst=((jnd1*(jnd1-1))/2+jnd2-1)*10
            do 310  i=1, 10
              vvl(i+indst)=rmre(ind1,ind2)*work(i+20)+
     :                     rmre2(ind1,ind2)*work(i+30)
310         continue
320       continue
350     continue
      endif
      return
      end
