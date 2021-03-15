************************************************************************
*                                                                       *
*                         hibridon 2  library                           *
*                                                                       *
*************************************************************************
*                          routines included:                           *
*                                                                       *
*   1. hiblk       block data, default settings                         *
*   2. hinput      input driver                                         *
*   3. logdb       log derivativ propagator                             *
*   4. mxoutd (mxoutr) mxoutc      matrix print utility
*   5. prsg/aver1/xscpr1  integral cross section print
*   6. prsgpi/aver2/xscpr2 integral cross section print for sigma/pi
*                                                                       *
*************************************************************************
* ------------------------------------------------------------------
      subroutine default
*  current revision date:  22-jan-2008 by mha
* ------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      character*40 jobnam,input,output,savfil
      logical logwr, swrit, t2writ, wrpart, partw, airyfl, airypr,
     :        ipos, noprin, chlist, wrxsec, xsecwr, writs, csflag,
     :        flaghf, rsflag, t2test, logdfl, flagsu, batch,
     :        readpt, ihomo, bastst, twomol, nucros, photof, wavefl,
     :        boundc
      parameter (maxpar=150)
      include "common/parpot"
      common /cofile/ input, output, jobnam, savfil
      common /cosysi/ nscode, isicod, ispar(maxpar)
      common /cosysr/ isrcod, junkr, rspar(maxpar)
      common /cosysl/ islcod, lspar(maxpar)
      common /cosys/ scod(maxpar*2+3)
* nb if the nextcommon is changed, it should be also changed in common/parsys
      common /coiout/ niout, indout(20)
      common /cosout/ nnout, jout(21)
      common /coselb/ ibasty
      common /coipar/ jtot1,jtot2,jtotd,jlpar,nerg,numax,numin,nud,
     :                lscreen, iprint
      common /corpar/ fstfac, rincr, rcut, rendai, rendld, rstart, spac,
     :                tolai, xmu
      common /colpar/ airyfl, airypr, bastst, batch, chlist,
     :                csflag, flaghf, flagsu, ihomo, ipos, logdfl,
     :                logwr, noprin, partw, readpt, rsflag, swrit,
     :                t2test, t2writ, twomol, writs, wrpart, wrxsec,
     :                xsecwr, nucros, photof, wavefl, boundc
*  this sets the maximum number of energies
      common /coener/ energ(1)
c     This common block is kept in this subroutine for backward
c     compatibility.
c     The constants module should be used in new subroutines.
      common /coconv/ econv, xmconv, ang2c
*
      jtot1=20
      jtot2=20
      jtotd=20
      jlpar=1
      nerg=1
      numax=0
      numin=0
      nud=1
* lscreen is the number of lines available on your terminal screen
      lscreen=48
* iprint controls degree of print output in some routines
*     iprint=-1 (no print); iprint=0 (min print); iprint=1 some print, etc
      iprint=0
      energ(1)=208.509d0
      do i=2,25
        energ(i)=0.d0
      enddo
      fstfac=15d0
      rincr=8d0
      rcut=30d0
      rendai=25d0
      rendld=8d0
      rstart=5.6d0
      spac=0.15d0
      tolai=1.15d0
      xmu=16.47d0
* econv is the conversion from hatree->cm-1
* xmconv is the conversion from au (electron mass)-> amu, i.e. the mass of the
* in amu
* ang2c is bohr^2 in angstroms^2
* default basis type is 1 (singlet sigma)
      ibasty=1
      econv=219474.6d0
      xmconv=5.485930d-4
      ang2c=0.280002846d0
      airyfl=.true.
      airypr=.false.
      batch=.false.
      chlist=.true.
      csflag=.false.
      flaghf=.false.
      flagsu=.false.
      ihomo=.true.
      ipos=.false.
      logdfl=.true.
      logwr=.false.
      noprin=.false.
      partw=.false.
      readpt=.false.
      rsflag=.false.
      swrit=.false.
      t2test=.false.
      t2writ=.true.
      writs=.true.
      wrpart=.false.
      wrxsec=.false.
      nucros=.false.
      xsecwr=.true.
      bastst=.false.
      photof=.false.
      wavefl=.false.
      boundc=.false.
      twomol=.false.
      nnout=-2
      do i=1,20
        jout(i)=0
      enddo
      jout(2)=2
* the following values are for a singlet sigma molecule in v=0
      niout=1
      indout(1)=0
* the following values are for a doublet sigma type molecule
*      data niout/2/,indout(1),indout(2)/1,-1/
      nscode=0
      isicode=0
      islcod=0
      isrcod=0
*  set up default names for input and output files and jobname
      input='Inpt'
      output='Outpt'
      jobnam='Job'
* default labels
      label=' LABEL NOT YET ASSIGNED'
      potnam=' POTNAM NOT YET ASSIGNED'

      end
* ----------------------------------------------------------------------
      subroutine enord(energ,nerg)
      implicit double precision (a-h,o-z)
      dimension energ(1)
      do 502 i=1,nerg
      emax=energ(i)
      jmax=i
      do 501 j=i+1,nerg
      if(energ(j).gt.emax) then
        emax=energ(j)
        jmax=j
      end if
501   continue
      if(jmax.ne.i) then
        energ(jmax)=energ(i)
        energ(i)=emax
      end if
502   continue
      return
      end
* -------------------------------------------------------------------
      subroutine logdb (z, w, amat, bmat, nmax, wref, z1, z2,
     :                  scr1, scr2, nch, rmin, rmax, nsteps,
     :                  eshift, iread, iwrite, tl, tp, twf)
*     routine to initialise the log derivative matrix, y, at r = rmin,
*     and propagate it from rmin to rmax using the method described in
*     d.e.manolopoulos, j.chem.phys., 85, 6425 (1986)
*     references to this paper appear below in comments
*     author:  david manolopoulos and millard alexander
*     current revision date (the propagator): 28-nov-2007
*     revised on 30-mar-2012 by q. ma for stream I/O of wfu files
*     current revision: 20-apr-2012 by q. ma
***********************************************************************
****   this integrator is not released for general public use      ****
****   all use must be by specific prior arrangement with:         ****
****     millard alexander, department of chemistry,               ****
****     university of maryland, college park, md, 20742           ****
****     tel: 1.301.405.1823; email: mha@umd.edu                   ****
****   no part of this program may be copied or used for other     ****
****   purposes without the author's permission.                   ****
***********************************************************************
*  ------------------------------------------------------------------
*     variables in call list:
*     z             array of dimension (nmax,nch)
*                   on return z contains the log derivative matrix
*                   at r = rmax
*     w             scratch array of dimension (nmax,nch)
*                   used as workspace for the potential matrix
*     amat,bmat     scratch matrices of dimension (nmax,nch)
*                   used as workspace for the dissociation overlap
*     nmax          leading dimension of arrays z and w
*     wref          scratch array of dimension nch
*                   used as workspace for the reference potential
*     z1,z2         scratch arrays of dimension nch
*                   used as workspace for the homogeneous
*                   half-sector propagators
*     scr1,scr2     scratch arrays of dimension nch
*                   used as workspace in calls to smxinv
*     nch           number of coupled equations
*                   = actual order of matrices z and w
*     rmin,rmax     integration range limits
*                   rmin is assumed to lie inside the classically
*                   forbidden region for all channels
*     nsteps        number of sectors partitioning integration range
*                   this version uses a constant step size throughout
*                   the integration range, h = (rmax-rmin)/(2*nsteps)
*           note:  if nsteps=0, then this routine only initializes the logd
*                  matrix; no propagation is done
*     eshift        energy shift for diagonal elements of the coupling
*                   matrix at subsequent energies:
*                   new energy = first energy + eshift
*     iread,iwrite  logical variables for i/o of energy-independent
*                   information from/to unit 11
*     tl            elapsed time used in logd integration
*                   exclusive of calls to potential and ground state
*                   wavefunction (in photodissociation calculations)
*     tp            elapsed time used in calls to potential
*     twf           elapsed time used in determination of ground
*                   state wavefunction (this should be zero for scattering)
*                   exclusive of calls to potential
*                   this timing information comes from repeated calls
*                   of the form call mtime(elapsed) where "elapsed"
*                   is the current elapsed time in seconds
*     variables in common block /cophot/
*     photof        true if photodissociation calculation
*                   false if scattering calculation
*     wavefn        true if g(a,b) transformation matrices are saved
*                   to be used later in computing the wavefunction
*     variables in common block /coqvec/
*     mxphot        maximum dimension of q matrix
*     nphoto        actual collumn dimension of q matrix
*     q             accumulated overlap matrix with ground state
*                   only calculated if photof = .true.
*                   this is stored with each wavefunction as a column vector
*     variable in common block /cowave/
*     irec          record number of last written g(a,b) matrix
*     ifil          local unit number for g(a,b) file
*     blas routines daxpy and dscal are used in o(n*n) loops
*     symmetry of the matrices z and w is not exploited in these loops,
*     and blas routines are not used for o(n) loops
*  ------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer ich, icode, icol, idiag, ierr, ij, irow, istep, kstep,
     :        nch, ncol, ndiag, nmax, nrow, nsteps
*      real arg, d3, d6, eight, eshift, fac, h, half, hi, one, r, rmax,
*     :     rmin, sqrt, tan, tanh, t0, t0w, t1, t1w, tf, tfw, th, tl,
*     :     tn, tp,  zero, wdiag
*      real w, z
*      real scr1, scr2, wref, z1, z2
      logical iread, iwrite, photof, wavefn, boundf, wrsmat
      logical flagsu
*      external mtime, potmat, daxpy, smxinv, dscal
*     matrices z, w, amat, and bmat are stored column by column as one-dimensi
      dimension z(1), w(1), amat(1), bmat(1)
      dimension wref(1), q(1), z1(1), z2(1), scr1(1), scr2(1)
      common /coqvec/ mxphot, nphoto, q
      common /cophot/ photof, wavefn, boundf, wrsmat
      common /cowave/ irec, ifil, nchwfu, ipos2, ipos3, nrlogd, iendwv,
     $     inflev
      common /cosurf/ flagsu
      data zero,  one,  two, three,six,  eight
     :    / 0.d0, 1.d0, 2.d0, 3.d0, 6.d0, 8.d0 /
      data izero, ione /0, 1/
c     The following variables are for size-determination of (machine
c     dependent) built-in types
      integer int_t
      double precision dble_t
      character char_t
c
*     make sure that rmin, rmax and nsteps are the same if second
*     energy calculation
      if (iwrite) write (11) rmin, rmax, nsteps
      if (iread ) read  (11) rmin, rmax, nsteps
      if (nsteps .ne. 0) then
        h = (rmax - rmin) / (2 * nsteps)
        hi = one / h
        simpwt= h/three
        d3 = h * h / three
        d6 = - h * h / six
        half = h / two
        if (photof) nqmax=nch*nphoto
*       initialize q vector
cstart unix cray
        call dset(nqmax, zero, q, 1)
cend
      end if
*     row, column and diagonal increments for matrices z and w
      nrow = 1
      ncol = nmax
      ndiag = nmax + 1
*     obtain coupling matrix, w, at r = rmin
*     diagonal elements must be shifted at subsequent energies
      tp = zero
      tpw = zero
      twf = zero
      twfw = zero
      call mtime(tf,tfw)
      if (iread) then
         icol = 1
         do 5  ich = 1, nch
            read (11) (w(ij), ij = icol, icol + nch - 1)
            icol = icol + ncol
   5     continue
         idiag = 1
         do 10  ich = 1, nch
            w(idiag) = w(idiag) - eshift
            idiag = idiag + ndiag
  10     continue
         if (photof) read (11) (q(i), i=1, nqmax)
      else
         istep = 0
         r = rmin
         call mtime(t0,t0w)
         call potmat(w, r, nch, nmax)
         call mtime(t1,t1w)
         tp = tp + t1 - t0
         tpw = tpw + t1w - t0w
*     determine ground state wavefunction times dipole moment at beginning of
*     first sector
         if (photof) then
           call mtime(t0,t0w)
           call ground(q, r, nch, nphoto, mxphot)
*     wt is weight for simpson's rule quadrature, isimp is alternator for
*     simpson's rule quadrature
           wt = simpwt
           isimp=1
*     multiply psi(0) mu(0) by simpson's rule wt at first point
cstart unix cray
           call dscal (nqmax, wt, q, 1)
cend
           call mtime(t1,t1w)
           twf=twf+t1-t0
           twfw=twfw+t1w-t0w
         endif
*
         if (iwrite) then
            icol = 1
            do 15  ich = 1, nch
               write (11) (w(ij), ij = icol,icol+nch-1)
               icol = icol + ncol
  15        continue
            if (photof) write (11) (q(i), i=1, nqmax)
         endif
      endif
*     use diagonal approximation to wkb initial value for log
*     derivative matrix  (eqn 16)
*     rmin is assumed to lie inside the classically forbidden
*     region in all channels  (w(ii) > 0)
*     first zero out z matrix
      icol = 1
      do 20  ich = 1, nch
cstart unix cray mac
         call dscal (nch, zero, z(icol), nrow)
cend
         icol = icol + ncol
  20  continue
      idiag = 1
      do  30 ich = 1, nch
         wdiag = w(idiag)
         if (wdiag .le. 0) then
           write (9, 25) ich, wdiag
           write (6, 25) ich, wdiag
25         format (' *** WAVEVECTOR > 0 IN CHANNEL', i5,
     :        '     IN LOGD INITIALIZATION;  ABORT ***',
     :       /'     V - E =',1pe15.7)
           call exit
         end if
         z(idiag) = sqrt(wdiag)
         idiag = idiag + ndiag
  30  continue

      if (nsteps .le. 0) then
        if (photof) then
          write (9, 35) nsteps
          write (6, 35) nsteps
35        format (' *** WARNING: IN LODG PHOTOF=TRUE AND NSTEPS=',i3,
     :            ' SHOULD BE .GT. 0 ***')
          return
        endif
      endif
*     with a constant step size it is convenient to propagate
*     the matrix z = h*y rather than the log derivative matrix, y
*     (eqns 10, 12, 13 and 14 are simply multiplied through by h)
      fac = h
      icol = 1
      do 40  ich = 1, nch

cstart unix cray
         call dscal (nch, fac, z(icol), nrow)
cend
         icol = icol + ncol
  40  continue
*     propagate z matrix from rmin to rmax
      do 250  kstep = 1, nsteps
*     apply quadrature contribution at beginning of sector,
*     r = a
*     after this loop z contains z(a)+(h^2/3)w(a)
         fac = d3
         icol = 1
         do  50 ich = 1, nch
            call daxpy (nch, fac, w(icol), nrow, z(icol), nrow)
            icol = icol + ncol
  50     continue
*     the reference potential for the sector is the diagonal
*     of the coupling matrix evaluated at the centre of the
*     sector, r = c  (eqn 15)
         if (iread) then
            read (11) (wref(ich), ich = 1, nch)
            do 60  ich = 1, nch
               wref(ich) = wref(ich) - eshift
  60        continue
         else
            istep = istep + 1
            r = rmin + istep * h
            call mtime(t0,t0w)
            call potmat(w, r, nch, nmax)
            call mtime(t1,t1w)
            tp = tp + t1 - t0
            tpw = tpw + t1w - t0w
            idiag = 1
            do 70  ich = 1, nch
               wref(ich) = w(idiag)
               idiag = idiag + ndiag
  70        continue
            if (iwrite) then
               write (11) (wref(ich), ich = 1, nch)
            endif
         endif
*     adjust quadrature contribution at r = a to account for
*     sector reference potential  (eqn 11)
*     after this loop z contains z(a)+(h^2/3)w1(a)
         fac = - d3
         idiag = 1
         do  80 ich = 1, nch
            z(idiag) = z(idiag) + fac * wref(ich)
            idiag = idiag + ndiag
  80     continue
*     evaluate homogeneous half sector propagators  (eqn 10)
*     z(i) = h * y(i), i=1,2,3,4
         do 85  ich = 1, nch
            arg = half * sqrt(abs(wref(ich)))
               if (wref(ich) .lt. zero) then
                  tn = tan(arg)
                  z1(ich) = arg / tn - arg * tn
                  z2(ich) = arg / tn + arg * tn
               else
                  th = tanh(arg)
                  z1(ich) = arg / th + arg * th
                  z2(ich) = arg / th - arg * th
               endif
*           z3(ich) = z2(ich)
*           z4(ich) = z1(ich)
  85     continue
*     propagate z matrix across the first half sector
*     from r = a to r = c  (eqn 14)
         idiag = 1
         do  90 ich = 1, nch
            z(idiag) = z(idiag) + z1(ich)
            idiag = idiag + ndiag
  90     continue
*     z now contains z(a)+z1(a,c)
cstart .not.unix-darwin .and. .not. unix-x86
c;         call smxinv(z, nmax, nch, scr1, scr2, ierr)
cend
cstart unix-darwin unix-x86
        call syminv(z,nmax,nch,ierr)
cend
         if (ierr .ne. 0) then
            write (6, 9000) kstep, ierr
            write (9, 9000) kstep, ierr
            call exit
         endif
*     z now contains [z(a)+z1(a,c)]^-1
         icol = 1
         do  95 ich = 1, nch
            fac = z2(ich)
            call dscal (nch, fac, z(icol), nrow)
            icol = icol + ncol
  95     continue
*     z now contains [z(a)+z1(a,c)]^-1 z2(a,c)
*     if photodissociation calculation or wavefunction desired:
*     save this matrix, which is hg(a,m), in amat
         if (photof .or. wavefn)
     :     call matmov(z, amat, nch, nch, nmax, nmax)
         irow = 1
         do  110 ich = 1, nch
            fac = - z2(ich)
            call dscal (nch, fac, z(irow), ncol)
            irow = irow + nrow
 110     continue
*     z now contains - z3(a,c) [z(a)+z1(a,c)]^-1 z2(a,c)
         idiag = 1
         do  120 ich = 1, nch
            z(idiag) = z(idiag) + z1(ich)
            idiag = idiag + ndiag
 120     continue
*     evaluate quadrature contribution at sector mid-point,
*     r = c  (eqn 12)     (first energy calculation only)
         if (iread) then
            icol = 1
            do  125 ich = 1, nch
               read (11) (w(ij), ij = icol,icol+nch-1)
               icol = icol + ncol
 125        continue
         else
            fac = d6
            icol = 1
            do  130 ich = 1, nch
               call dscal (nch, fac, w(icol), nrow)
               icol = icol + ncol
 130        continue
            idiag = 1
            do  140 ich = 1, nch
               w(idiag) = one
               idiag = idiag + ndiag
 140        continue
cstart .not.unix-darwin .and. .not. unix-x86
c;            call smxinv(w, nmax, nch, scr1, scr2, ierr)
cend
cstart unix-darwin unix-x86
            call syminv(w,nmax,nch,ierr)
cend
            if (ierr .eq. 2) then
               icode = 2
               write (9,9000) kstep, icode
               call exit
            endif
            idiag = 1
            do  150 ich = 1, nch
               w(idiag) = w(idiag) - one
               idiag = idiag + ndiag
 150        continue
            if (iwrite) then
               icol = 1
               do  155 ich = 1, nch
                  write (11) (w(ij), ij = icol,icol+nch-1)
                  icol = icol + ncol
 155           continue
            endif
         endif
*     apply quadrature contribution at sector mid-point
*     corrections to z4(a,c) and z1(c,b) are applied
*     simultaneously  (eqn 13)
         fac = eight
         icol = 1
         do  160 ich = 1, nch
            call daxpy (nch, fac, w(icol), nrow, z(icol), nrow)
            icol = icol + ncol
 160     continue
*     propagate z matrix across the second half sector
*     from r = c to r = b   (eqn 14)
         idiag = 1
         do  170 ich = 1, nch
            z(idiag) = z(idiag) + z1(ich)
            idiag = idiag + ndiag
 170     continue
*     at this point z contains z(c) + z1(c,b)
cstart .not.unix-darwin .and. .not. unix-x86
c;         call smxinv(z, nmax, nch, scr1, scr2, ierr)
cend
cstart unix-darwin unix-x86
         call syminv(z,nmax,nch,ierr)
cend
         if (ierr .eq. 2) then
            icode = 3
            write (9,9000) kstep, icode
            call exit
         endif
*     at this point z contains [z(c) + z1(c,b)]^-1
         icol = 1
         do 175 ich = 1, nch
            fac = z2(ich)
            call dscal (nch, fac, z(icol), nrow)
            icol = icol + ncol
 175     continue
*     z now contains [z(c)+z1(c,b)]^-1 z2(c,b)
*     if photodissociation calculation or wavefunction desired:
         if (photof .or. wavefn) then
*     first save this matrix, which is g(c,b), in bmat
           call matmov(z, bmat, nch, nch, nmax, nmax)
*     use bmat and w as temporary storage here
           call mxma (amat, 1, nmax, bmat, 1, nmax, w, 1, nmax,
     :                nch, nch, nch)
*     w now contains the matrix g(a,m)g(m,b)=g(a,b)
*     if wavefunction desired, save this matrix
          if (wavefn .and. wrsmat) then
             irec = irec + 1
c     nrlogd is the number of LOGD records - used to seek the wfu file
             nrlogd = nrlogd + 1
             write (ifil, err=950) r - h, r, (w(i), i=1, nch)
             icol = 1
             do ich = 1, nch
                write (ifil, err=950) (w(icol - 1 + i), i=1, nch)
                icol = icol + nmax
             end do
             write (ifil, err=950) 'ENDWFUR', char(mod(irec, 256))
             lrlogd = (nchwfu ** 2 + nchwfu + 2) * sizeof(dble_t)
     $            + 8 * sizeof(char_t)
             iendwv = iendwv + lrlogd
          endif
        endif
*     accumulate overlap matrix with ground state
*     if photodissociation calculation
        if (photof) then
*     premultiply g(a,b) by wt psi(a) mu(a)
*     use bmat as temporary storage here
           call mxma(q,nch,1,w,1,nmax,bmat,nch,1,nphoto,nch,nch)
*     bmat now contains [...+wt*psi(a)mu(a)] g(a,b)
*       stored as successive columns with each initial
*       state corresponding to a column
*     now determine psi(b)mu(b), save this in q
           rnew = rmin + (istep+1) * h
           if (iread) then
             read (11) (q(i), i=1, nqmax)
           else
             call mtime(t0,t0w)
             call ground(q, rnew, nch, nphoto, mxphot)
* recalculate simpson's rule wt for this point
             wt=(3.d0+isimp)*simpwt
             isimp=-isimp
*     multiply psi(b) mu(b) by simpson's rule wt at r=rnew
             call dscal (nqmax, wt, q, 1)
             call mtime(t1,t1w)
             twf=twf+t1-t0
             twfw=twfw+t1w-t0w
             if (iwrite) write (11) (q(i), i=1, nqmax)
           endif
*     add wt*psi(b)mu(b) to bmat and resave
           fac=one
           call daxpy(nqmax, fac, bmat, 1, q, 1)
         endif
*     now premultiply [z(c)+z1(c,b)]^-1 z2(c,b) by z2(c,b)
         irow=1
         do  190 ich = 1, nch
           fac = - z2(ich)
           call dscal (nch, fac, z(irow), ncol)
           irow = irow + nrow
 190     continue
*      z now contains - z2(c,b) [z(c)+z1(c,b)]^-1 z2(c,b)
         idiag = 1
         do  195 ich = 1, nch
            z(idiag) = z(idiag) + z1(ich)
            idiag = idiag + ndiag
 195     continue
*     apply reference potential adjustment to quadrature
*     contribution at r = b  (eqn 11)
         fac = - d3
         idiag = 1
         do  200 ich = 1, nch
            z(idiag) = z(idiag) + fac * wref(ich)
            idiag = idiag + ndiag
  200    continue
*     obtain coupling matrix, w, at end of sector
         if (iread) then
            icol = 1
            do  205 ich = 1, nch
               read (11) (w(ij), ij = icol, icol + nch - 1)
               icol = icol + ncol
 205        continue
            idiag = 1
            do  210 ich = 1, nch
               w(idiag) = w(idiag) - eshift
               idiag = idiag + ndiag
 210        continue
         else
            istep = istep + 1
            r = rmin + istep * h
            call mtime(t0,t0w)
            call potmat(w, r, nch, nmax)
            call mtime(t1,t1w)
            tp = tp + t1 - t0
            tpw = tpw + t1w - t0w
            if (iwrite) then
               icol = 1
               do  215 ich = 1, nch
                  write (11) (w(ij), ij = icol,icol+nch-1)
                  icol = icol + ncol
 215           continue
            endif
         endif
*     apply quadrature contribution at r = b  (eqn 12)
         fac = d3
         icol = 1
         do 220  ich = 1, nch
            call daxpy (nch, fac, w(icol), nrow, z(icol), nrow)
            icol = icol + ncol
 220     continue
*     propagation loop ends here
 250  continue
*     recover the log derivative matrix, y, at r = rmax
      fac = hi
      icol = 1
      do 260  ich = 1, nch
         call dscal (nch, fac, z(icol), nrow)
         icol = icol + ncol
 260  continue
      call mtime(tl,tlw)
      tl = tl - tf - tp - twf
      tlw = tlw - tfw - tpw - twfw
      return
c
 950  write (0, *) ' *** ERROR WRITING WFU FILE (LOGD). ABORT'
      call exit()
      return
c
9000  format(' *** MATRIX INVERSION ERROR IN LOGDB AT KSTEP =',
     : i4,  /'            ERROR OCCURRED AT ROW =', i2,
     : '; ABORT ***')
      end
* ----------------------------------------------------------------------
      subroutine mxoutd (iunit, xmat, nn, nmax, isym, ipos)
*  to print out nn*nn matrix xmat of maximum dimension nmax use mxoutd
*  to print out nn*m matrix xmat of maximum dimension nmax use mxoutr
*  author:  millard alexander
*  current revision date: 1-may-90 by mha
*  -------------------------------------------------------------------------
*  variables in call list:
*    iunit:     logical unit number
*    xmat:      matrix to be printed out
*    nn:         actual size of matrix (nn x nn)
*    nmax:      maximum row dimension of matrix
*    isym:      if isym = 1,  matrix is symmetrical, only lower triangle
*               is printed
*               if isym = 0 full matrix is printed
*               note that isym must be either 0 or 1 in calling sequence in ma
*               program
*    ipos:      if .true., printer is assumed to have 132 positions/line
*               if .false., printer is assumed to have 80 positions/line
*  -------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical ipos
      dimension xmat(nmax,nmax), ind(10)
      m=nn
      n=nn
      goto 1
      entry mxoutr (iunit, xmat, nrow, ncol, nmax, isym, ipos)
      m=nrow
      n=ncol
*  abort if isym not equal to zero or one
1     if (isym. ne. 0 .and. isym .ne. 1) then
        write (iunit, 5)
        if(iunit.ne.6) write (6, 5)
5       format (/' ***  ISYM .NE. 0 OR 1 IN MXOUTD; ABORT ***')
        call exit
      end if
*  if 132 line printer, then 10 columns of matrix are printed simultaneously
*  if  80 line printer, then  6 columns of matrix are printed simultaneously
      iskip = 10
      if (.not. ipos) iskip = 6
*  jmax is the total number of the groups of columns to be printed
      jmax = n / iskip + 1
      jlow = 1
      jhigh = min (iskip, n)
*  loop over the groups of columns from column jlow to jhigh
      do  50   j = 1, jmax
*  ncol is the number of columns to be printed
        nc = jhigh - jlow + 1
        do  10   jj = 1, nc
*  the array ind contains the index of each column which will be printed
          ind(jj) = jj - 1 + jlow
10      continue
*  write as a heading the column index of each column to be printed
        write (iunit, 15) ( ind(i), i = 1, nc )
15      format (/5x, 10 (3x, i3, 6x) )
*  now loop through the rows of the matrix, which will be printed out
*  in groups of 10 with a blank line in between
*  lowrow is the index of the first row to be printed
*  if the full matrix is desired (isym = 0), then lowrow = 1
*  if the lower triangle is desire (isym = 1), then lowrow = jlow
        lowrow = isym * jlow + (1 - isym)
*  loop over the rows of the matrix which will be printed
        do  35   jrow =  lowrow, m
*  jtop points to the maximum column which will be printed in this row
          jtop = isym * min (jhigh, jrow) + (1 - isym) * jhigh
*  now write out the row index followed by the desired matrix elements
          jj=jrow/10
          write (iunit,20) jrow,( xmat(jrow,jcol),jcol = jlow,jtop)
20        format (i4, 1x, 10 (1pe12.4) )
*  if this row is an integer multiple of 10, add a blank line
          jj = jrow / 10
          if (10 * jj .eq. jrow .and. jrow .ne. m) write (9, 25)
25        format(1h )
35      continue
*  all rows have been printed out for this group of columns
*  move to the next group of columns
*  jlow becomes old jhigh
*  jhigh is set equal to jhigh plus skip distance
        if (jhigh .eq. n) go to 60
        jlow = jhigh + 1
        jhigh = min ( (jlow + iskip - 1), n)
50    continue
*  entire matrix has been printed, return
60    return
      end
* --------------------------------------------
      subroutine mxoutc (ifile,zmat,nlevop,nmax,ipos,csflag,flaghf,
     :                   twomol,numax,jlev,inlev)
      use constants
      implicit double precision (a-h,o-z)
      logical csflag,flaghf,twomol,ipos, csff
      character*40 form
      common /coisc2/ nj,jlist(1)
      common /coipar/ ipar(9), iprint
      dimension zmat(nmax,nlevop),inlev(1),jlev(1)
* routine to print integral cross sections
*  latest revision date:  24-jun-1991
      spin=0
      csff=csflag
      if (csflag .and. iprint .ge. 2) then
* here for full print of cs cross sections, even if not converged
        write (ifile,1) numax
1       format(/'** COUPLED STATES CROSS SECTIONS NOT CONVERGED',
     :         ' FOR J .GT. NUMAX = ',i3)
        csff=.false.
      endif
      if(twomol) then
        if(ipos) write(ifile,5) (j,j=1,nlevop)
        if(.not.ipos) write(ifile,6) (j,j=1,nlevop)
5       format(/1x,'  INITIAL STATE',t27,'FINAL STATE'/
     :          1x,'  N   J1   J2  INDEX',(t21,i7,9i11))
6       format(/1x,'  INITIAL STATE',t27,'FINAL STATE'/
     :          1x,'  N   J1   J2  INDEX',(t21,i7,4i11))
* don't remove blank in this format!!
        form='(1x,i3,2i5,i5, (t21,10(1pd11.3)))'
      else
        if(ipos) write(ifile,10) (j,j=1,nlevop)
        if(.not.ipos) write(ifile,11) (j,j=1,nlevop)
10      format(/1x,'  INITIAL STATE',t24,'FINAL STATE'/
     :          1x,'  N   J   INDEX ',(t18,i7,9i11))
11      format(/1x,'  INITIAL STATE',t24,'FINAL STATE'/
     :          1x,'  N   J   INDEX ',(t18,i7,4i11))
        if(flaghf) then
          spin=0.5d0
          form='(1x,i3,f5.1,i6,(t18,10(1pd11.3)))'
        else
* don't remove blank in this format!!
cjk          form='(1x,i3,i4,i7,  (t18,10(1pd11.3)))'
         form='(1x,i3,i4,i7,  (t18,10(1pes22.9e3)))'
        end if
      end if
      if(.not.ipos) form(21:22)=' 5'
      if(csff) then
        nj=0
        do 100 i=1,nlevop
        if(jlev(i).le.numax) then
          nj=nj+1
          jlist(nj)=i
        end if
100     continue
        do 120 ii=1,nj
        i=jlist(ii)
        if(twomol) then
          j2 = mod( jlev(i), 10)
          j1 = jlev(i) / 10
          write(ifile,form) i,j1,j2,inlev(i),(zmat(i,j),j=1,nlevop)
        else if(flaghf) then
          write(ifile,form) i,jlev(i)+spin,inlev(i),
     :                     (zmat(i,j),j=1,nlevop)
        else
          write(ifile,form) i,jlev(i),inlev(i),(zmat(i,j),j=1,nlevop)
        end if
120     continue
        if(ipos) write(ifile,10) (jlist(j),j=1,nj)
        if(.not.ipos) write(ifile,11) (jlist(j),j=1,nj)
        do 140 i=1,nlevop
        if(jlev(i).le.numax) goto 140
        if(twomol) then
          j2 = mod( jlev(i), 10)
          j1 = jlev(i) / 10
          write(ifile,form) i,j1,j2,inlev(i),(zmat(i,jlist(j)),j=1,nj)
        else if(flaghf) then
          write(ifile,form) i,jlev(i)+spin,inlev(i),
     :                     (zmat(i,jlist(j)),j=1,nj)
        else
          write(ifile,form) i,jlev(i),inlev(i),(zmat(i,jlist(j)),j=1,nj)
        end if
140     continue
      else
        do 150 i=1,nlevop
        if(twomol) then
          j2 = mod( jlev(i), 10)
          j1 = jlev(i) / 10
          write(ifile,form) i,j1,j2,inlev(i),(zmat(i,j),j=1,nlevop)
        else if(flaghf) then
          write(ifile,form) i,jlev(i)+spin,inlev(i),
     :                      (zmat(i,j),j=1,nlevop)
        else
          write(ifile,form) i,jlev(i),inlev(i),(zmat(i,j),j=1,nlevop)
        end if
150     continue
      end if
      return
      end
*  ------------------------------------------------------------------
      subroutine prsg (fname, a)
*  subroutine to write out selected integral cross sections
*  from file {fname1}.ics
*  author:  millard alexander
*  last bug fix:17-may-1993 by mha
*  latest revision date:  23-feb-2013 by p. dagdigian
*  ------------------------------------------------------------------
*  variables in call list:
*    zmat:    on return:  contains the nlevop x nlevop matrix of integral
*                         cross sections
*    jlev:   rotational angular momenta of each energetically open level
*    elev:   energy (in hartree) of each energetically open level
*    inlev:  additonal quantum index for each energetically open level
*    jfirst:  initial value of total angular momentum
*    jfinal:  final value of total angular momentum
*    jtotd:   step size for total angular momentum
*    ipos:    if .true., 132 column printer
*             if .false., 80 column printer
*    csflag:  if .true. coupled-states calculation
*             if .false., close-coupled calculation
*    flaghf:  if .true., then system with half-integer spin
*             if .false., then system with integer spin
*    twomol:  if .true., then molecule-molecule collision
*             if .false., then atom-molecule or molecule-surface collision
*    nlevop:  number of energetically distinct levels in channel basis which
*             are open asymptotically
*    numin:   initial value of the coupled-states projection index
*    numax:   final value of the coupled-states projection index
*    nud:     step in coupled-states projection index
*    jlpar:   parity of channels
*               eps * (-1)**(j+l-jtot) = jlpar
*             if jlpar = 0, then integral cross sections include contributions
*             of both parities
*    note!!!   if flaghf = .true.( see below), then the true values
*    of the rotational quantum numbers, the total angular momentum,
*    and the coupled-states projection index are equal to the values
*    stored in jlev, jtot, jfirst, jfinal, nu, numin, and numax plus 1/2
*    flagsu:    if .true., then molecule-surface collisons
*    xmu:       collision reduced mass in (c12) atomic mass units
*    econv:     conversion factor from cm-1 to hartrees
*  ------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      character*(*) fname
      character*20 cdate
      character*3 stat
      character*12 accs
      character*40 xnam1, xnam2
      logical csflag, flaghf, iprint, ipos, flagsu, twomol, existf,
     :        openfl,lpar(17),
     :        airyfl,airypr,bastst,batch,chlist,ihomo, eprint
      common /colpar/ airyfl, airypr, bastst, batch, chlist,
     :                csflag, flaghf, flagsu, ihomo, ipos,lpar
*     real econv, ener, xmu
*     real a, elev, scmat, zmat
      integer i, ienerg, iout, isize, j, jbegin, jend, jfinal,
     :        jfirst, jj1, jj2, jlpar, jtemp, jtotd, lenx, n, nlevel,
     :        nlevop, nnout, nout, numax, numin, nud, iaver
      integer inlev, ipoint, jlev, jout
      include "common/parpot"
      common /coz/ zmat(1)
      common /cow/ scmat(1)
      common /cojhld/ jlev(5)
      common /coisc1/ inlev(5)
      common /coisc2/ jpoint(5)
      common /colq/ ipoint(5)
      common /cosc1/ elev(5)
      common /cosout/ nnout, jout(21)
      common /coiout/ niout, indout(5)
      common /coselb/ ibasty
      dimension  a(3)
      data econv / 219474.6d0/
*  input parameters
      iprint=.true.
      eprint=.false.
      if (a(1) .lt. 0.0) iprint =  .false.
      if (a(1) .gt. 0.0) eprint = .true.
      iaver=nint(a(2))
      ienerg = a(3) + 0.1
      if (ienerg .le. 0) ienerg=1
*  open integral cross section file
      call gennam (xnam1, fname, ienerg, 'ics', lenx)
      inquire (file = xnam1, exist = existf)
      if (.not. existf) then
        write (6, 10) xnam1(1:lenx)
10      format(/' integral cross section file ',(a),' not found',/)
        return
      end if
      open (1, file = xnam1,status = 'old')
*  open output file for integral cross sections
      call gennam (xnam2, fname, ienerg, 'xsc', lenx)
      if (.not.iprint) write(6,20) xnam2
20    format(' PRINTING SELECTED CROSS SECTIONS; OUTPUT IN FILE ',
     :        (a))
* inquire file specifications
      inquire(file=xnam2, exist=existf, opened=openfl)
      accs='sequential'

      if (.not. openfl) then
        if (existf) then
          stat='old'
* make sure sequential formatted files are appended not overwritten
cstart unix-hp unix-dec unix-iris unix-sun unix-ifort unix-pgi
          accs='append'
cend
        else
          stat='new'
        end if
        open (3, file = xnam2, status = stat, access = accs)
      endif
      call version(3)
      read (1, 40) cdate
40    format (1x, a)
      read (1, 40) label
      read (1, 40) potnam
*  print job information
      write (3, 50) xnam1, cdate, label, potnam
      if (iprint) write (6, 50) xnam1, cdate, label, potnam
50    format(/' INTEGRAL CROSS SECTIONS READ FROM FILE ',(a)/
     +        ' WRITTEN:   ',(a)/
     +        ' LABEL:     ',(a)/,
     +        ' POT NAME:  ',(a) )
      read (1, 60) ener, xmu
60    format (f10.3, f15.11)
      read (1, 70) csflag, flaghf, flagsu, twomol
70    format (4l3)
      read (1, 75) jfirst, jfinal, jtotd, numin, numax,
     :             nud, jlpar
75    format (7i5)
      if (ipos) then
80      format (24i5)
        read (1, 80)  nlevel, nlevop
        read (1, 80) (jlev(i), inlev(i), i = 1, nlevel)
        read (1, 90) (elev(i), i = 1, nlevel)
90      format (8(1pe15.8))
*90      format (8f16.9)
      else
        read (1, *)  nlevel, nlevop
        read (1, *) (jlev(i), inlev(i), i = 1, nlevel)
        read (1, *) (elev(i), i = 1, nlevel)
      endif
*  read in matrix of cross sections, column by column
      do 95  i = 1, nlevop
*  jbegin and jend point to first element and last element of column i of
*  matrix packed column by column
        jbegin = (i - 1) * nlevop + 1
        jend = i * nlevop
        if (ipos) then
          read (1, 90) (zmat(j), j = jbegin, jend)
        else
          read (1, *) (zmat(j), j = jbegin, jend)
        endif
95    continue
      if (.not. flagsu) then
        if (.not.flaghf .or. ibasty.eq.12) then
          if (iprint) write (6, 100) ienerg, xmu, ener, jlpar,
     :                                jfirst, jfinal, jtotd
          write (3, 100) ienerg, xmu, ener, jlpar, jfirst,
     :                   jfinal, jtotd
100       format (/' ** INTEGRAL CROSS SECTIONS; IEN=', i2,' **',
     :            /'    RMU=', f9.4, '  E=', f9.2,'  JLPAR=', i2,
     :             '  JTOT-1=', i3,
     :             '  JTOT-2=', i4,'  JTOT-D=', i3)
        else
          if (iprint) write (6, 105) ienerg, xmu, ener, jlpar,
     :                      (jfirst+0.5), (jfinal+0.5), jtotd
          write (3, 105) ienerg, xmu, ener, jlpar, (jfirst+0.5),
     :                   (jfinal+0.5), jtotd
105       format (/' ** INTEGRAL CROSS SECTIONS; IEN=', i2,' **',
     :            /'    RMU=', f9.4, '  E=', f9.2,'  JLPAR=', i2,
     :             '  JTOT-1=', f5.1,
     :             '  JTOT-2=', f6.1,'  JTOT-D=', i3)
        end if
      else
        if (.not. flaghf) then
          if (iprint) write (6, 110) ienerg, xmu, ener, numin,
     :                               numax, nud
          write (3, 110) ienerg, xmu, ener, numin, numax, nud
110       format (/' ** SUMMED DEGENERACY AVERAGED TRANSITION',
     :             ' PROBABILITIES;  IEN=', i2,' **',
     :            /'    RMU=', f9.4, '  E=', f8.2,'  M-MIN=', i3,
     :             '  M-MAX=', i4, '  M-STEP=', i2)
        else
          if (iprint) write (6, 115) ienerg, xmu, ener, (numin+0.5),
     :                               (numax+0.5), nud
          write (3, 115) ienerg, xmu, ener, (numin+0.5),
     :                   (numax+0.5), nud
115       format (/' ** SUMMED DEGENERACY AVERAGED TRANSITION',
     :             ' PROBABILITIES;  IEN=', i2,' **',
     :            /'    RMU=', f9.4, '  E=', f8.2,'  M-MIN=', f5.1,
     :             '  M-MAX=', f6.1, '  M-STEP=', i2)
        end if
      end if
      if (.not. csflag) then
        if (jlpar .eq. 0) then
          write (3, 120)
          if (iprint) write (6, 120)
120       format (' ** CC CALCULATION, BOTH PARITIES **')
        else
          write (3, 125) jlpar
          if (iprint) write (6, 125) jlpar
125       format (' ** CC CALCULATION, JLPAR=', i2, ' **')
        end if
      else
        if (.not. flaghf) then
          write (3, 130) numin, numax, nud
          if (iprint) write (6, 130) numin, numax, nud
130       format (' ** CS CALCULATION, NUMIN=', i2,', NUMAX=',
     :               i2,' NUD=', i2, ' **')
        else
          write (3, 140) numin + 0.5, numax + 0.5, nud
          if (iprint) write (6, 140) numin+0.5, numax+ 0.5, nud
140       format (' ** CS CALCULATION, NUMIN=', f4.1, ', NUMAX=',
     :            f4.1,' NUD=', i2, ' **')
        end if
      end if
      if (eprint) then
        if (.not. twomol) then
          if (.not. flagsu) then
            if (iprint) write (6, 145)
            write (3, 145)
145         format
     :       (' LEVEL LIST FOR INTEGRAL CROSS SECTIONS',
     :        /'   N   J  INDEX  EINT(CM-1)',/)
          else
            if (iprint) write (6, 150)
            write (3, 150)
150         format
     :       (' LEVEL LIST FOR DEGENERACY AVERAGED',
     :         ' TRANSITION PROBABILITIES',
     :        /'   N   J  INDEX  EINT(CM-1)',/)
          end if
          do 170  i = 1, nlevop
            if (.not. flaghf) then
              if (iprint)
     :        write (6, 160) i, jlev(i), inlev(i), elev(i) * econv
              write (3, 160) i, jlev(i), inlev(i), elev(i) * econv
160           format (i4, i5, i6, f11.3)
            else
              if (iprint) write (6, 165) i, (jlev(i)+0.5), inlev(i),
     :                       elev(i) * econv
              write (3, 165) i, (jlev(i)+0.5), inlev(i),
     :                       elev(i) * econv
165           format (i4, f5.1, i6, f11.3)
            end if
170       continue
        else
          if (iprint) write (6, 175)
          write (3, 175)
175         format (' LEVEL LIST FOR INTEGRAL CROSS SECTIONS',
     :              /'   N   J1   J2  INDEX  EINT(CM-1)'/)
          do 190  i = 1, nlevop
            jj2 = mod( jlev(i), 10)
            jj1 = jlev(i) / 10
            if (iprint)
     :      write (6, 180) i, jj1, jj2, inlev(i), elev(i) * econv
            write (3, 180) i, jj1, jj2, inlev(i), elev(i) * econv
180         format (i4, 2i5, i6, f11.3)
190       continue
        end if
      endif
*  now sum and average over positive and negative values of index
*     check that number of levels with index negative equals number of levels
*     with index positive, if not abort
*     this test works by (1) checking that no index equals zero
*     and, if so, (2) adding all the values of the index.  this sum
*     should equal zero if there are as many levels with index negative
*     as index positive
      if (iaver .gt. 0) then
        isum = 0
        do  200 i = 1, nlevop
          if (inlev(i) .eq. 0) then
            write (6, 195) i
            write (6, 195) i
195         format(' *** INLEV(',i3,')=0;',
     :               ' AVERAGING MAY NOT WORK ***')
          else if (inlev(i) .ne. 0) then
            isum = isum + inlev(i)
          end if
200     continue
        if (isum .ne. 0) then
          write (6, 230) isum
          write (3, 230) isum
230       format (' *** SUM OF INDICES =',i3,
     :               ' AVERAGING MAY NOT WORK ***')
        end if
        if  (iaver .eq. 2) then
          if (iprint) write (6, 245)
          write (3, 245)
245       format
     :      (' ** CROSS SECTIONS SUMMED AND AVERAGED OVER INDEX **')
          call aver1 (zmat, scmat, nlevop)
         else if (iaver .eq. 1) then
          if (iprint) write (6, 250)
          write (3, 250)
250       format
     :      (' ** CROSS SECTIONS SUMMED OVER FINAL STATE INDEX **')
         end if
      end if
*  find all rows of cross sections matrix for which initial rotational
*  quantum number is equal to one of the values of jout
      isize = 0
      nout = abs (nnout)
      do  280  iout = 1, nout
        jtemp = jout(iout)
        do 270 n = 1, nlevop
          if (jlev(n) .eq. jtemp) then
            isize = isize + 1
            jpoint(isize) = n
            ipoint(isize) = n
          end if
270     continue
280   continue
caber
      insize=0
caber
      if (niout .gt. 0) then
        nout=abs(niout)
        do  282  iout = 1, nout
          indtemp=indout(iout)
          do 281 n = 1, isize
            if (inlev(jpoint(n)) .eq. indtemp) then
              insize = insize + 1
              ipoint(insize) = jpoint(n)
            end if
281       continue
282     continue
        isize=insize
      endif
*  isize is the number of cross sections to be printed
      if (isize .eq. 0) then
*  here if no initial states found
        write (6, 283)
        write (3, 283)
283     format (' ** NO INITIAL STATES FOUND; ABORT')
      else
*  now print out desired columns of cross section matrix
        write (3, 260)
        if (iprint) write (6, 260)
260     format (/' ** COLUMN HEADINGS ARE INITIAL STATES, ROW',
     :        ' HEADINGS ARE FINAL STATES **')
        call xscpr1(zmat, nlevop, isize, iaver, ipos, iprint, flaghf)
      endif
      close (1)
      close (3)
      return
      end
      subroutine aver1 (zmat, scmat, n)
      implicit double precision (a-h,o-z)
*  subroutine to sum and average cross section matrix over positive
*  and negative values of index
      integer i, ind, index, j, jnd, n, nn, inlev, jlev
*     real scmat, zmat
      dimension zmat(1), scmat(n,n)
      common /cojhld/ jlev(1)
      common /coisc1/ inlev(1)
*  first copy cross section matrix into scmat
      call matmov (zmat, scmat, n, n, n, n)
      nn = n / 2
      index = 0
      do  30  i = 1, nn
        ind = 2 * i
        do  20  j = 1, nn
          jnd = 2 * j
          index = index + 1
          zmat(index) = ( scmat(ind - 1, jnd - 1) +
     :                   scmat(ind - 1, jnd) +
     :                   scmat(ind, jnd - 1) +
     :                   scmat(ind, jnd) ) * 0.5
20      continue
         jlev (i) = jlev (ind - 1)
        inlev(i) = iabs(inlev(ind - 1))
30    continue
*  because of indexing zmat need to be transposed
      call transp (zmat, nn, nn)
*  reduce the size of the cross section matrix
      n = nn
      return
      end
* -------------------------------------------------
      subroutine xscpr1 (zmat, nlevop, isize, iaver,
     :                   ipos, iprint, flaghf)
      implicit double precision (a-h,o-z)
*      current revision date: 16-dec-2007
*  subroutine to print out specified columns of cross section matrix
*  if iaver = 1, then nth and (n+1) st rows are added before printing
      integer i, isize, iskip, j, jcol, jhigh, jj, jlow, jmax,
     :        jrow, ncol, nlevop, iaver
      integer ind, inlev, ipoint, jlev
*     real elev, zmat
      logical ipos, iprint, flaghf
      dimension zmat(nlevop,nlevop), ind(50)
      common /cojhld/ jlev(4)
      common /coisc1/ inlev(4)
      common /colq/ ipoint(4)
      common /cosc1/ elev(4)
      common /coselb/ ibasty
*     first transpose the cross section matrix so that initial
*     states are columns and final states are rows
      call transp (zmat, nlevop, nlevop)
*  if 132 line printer, then 13 columns of matrix are printed simultaneously
*  if  80 line printer, then  7 columns of matrix are printed simultaneously
      iskip = 13
      if (.not. ipos) iskip = 7
      irowsk = 1
      if (iaver .eq. 1) irowsk = 2
*  jmax is the total number of the groups of columns to be printed
      if (mod(isize, iskip) .eq. 0) then
        jmax = isize / iskip
      else
        jmax = isize / iskip + 1
      end if
      jlow = 1
      jhigh = min (iskip, isize)
*  loop over columns by groups of 6 or 10
      do  150   j = 1, jmax
*  ncol is the number of columns to be printed
        ncol = jhigh - jlow + 1
        do  10   jj = 1, ncol
*  the array ind contains the index of each column which will be printed
          ind(jj) = ipoint(jj - 1 + jlow)
10      continue
*  write as a heading the column index of each column to be printed
        if (.not. flaghf .or. ibasty.eq.12) then
          if (iprint) write (6, 15) ( jlev(ind(i)), i = 1,ncol)
          write (3, 15) ( jlev(ind(i)), i = 1, ncol )
15        format (/12x,'J=', i5, 2x, 12 (2x, i6, 2x) )
        else
          if (iprint) write (6, 30) ( jlev(ind(i))+0.5,
     :                                 i = 1,ncol)
          write (3, 30) ( jlev(ind(i))+0.5, i = 1,ncol)
30        format (/11x,'J= ', f5.1, 2x, 12 (2x, f6.1, 2x) )
        end if
        if (iprint) write (6, 40) ( inlev(ind(i)), i = 1,ncol)
        write (3, 40) ( inlev(ind(i)), i = 1,ncol)
40      format ('   J    I | I=', i5, 2x, 12 (2x, i6, 2x))
        if (iprint) write (6, 50)
        write (3, 50)
50      format (1h )
*  now loop through the rows of the matrix, which will be printed out
*  in groups of 10 with a blank line in between
*  loop over the rows of the matrix which will be printed
        do  95   jrow =  1 , nlevop, irowsk
*  now write out the row index followed by the desired matrix elements
          jj=jrow/10
*  inrow holds additional quantum index for this row
*  if iaver = 1, then this is positive, since both indices are summed
          inrow = inlev(jrow)
          if (iaver .ne. 1) then
            if (.not. flaghf .or. ibasty.eq.12) then
              if (iprint)
     :          write (6, 60) jlev(jrow),inrow,
     :            ( zmat(jrow,ind(jcol)), jcol = 1,ncol)
              write (3, 60) jlev(jrow),inrow,
     :          ( zmat(jrow,ind(jcol)), jcol = 1,ncol)
cjk 60            format (i5, i5, 2x, 13 (1pe10.3,1x) )
 60            format (i5, i5, 2x, 13 (1pes12.4e3,1x) )
            else
              if (iprint)
     :          write (6, 70) jlev(jrow)+0.5, inrow,
     :            ( zmat(jrow,ind(jcol)), jcol = 1,ncol)
              write (3, 70) jlev(jrow)+0.5, inrow,
     :          ( zmat(jrow,ind(jcol)), jcol = 1,ncol)
cjk 70            format (f5.1, i5, 2x, 13 (1pe10.3,1x) )
 70            format (f5.1, i5, 2x, 13 (1pes12.4e3,1x) )
            end if
          else
            if (.not. flaghf .or. ibasty.eq.12) then
              if (iprint)
     :          write (6, 60) jlev(jrow),inrow,
     :            (zmat(jrow,ind(jcol)) + zmat(jrow + 1,ind(jcol)),
     :                jcol = 1,ncol)
              write (3, 60) jlev(jrow),inrow,
     :            (zmat(jrow,ind(jcol)) + zmat(jrow + 1,ind(jcol)),
     :                jcol = 1,ncol)
            else
              if (iprint)
     :          write (6, 70) jlev(jrow)+0.5, inrow,
     :            (zmat(jrow,ind(jcol)) + zmat(jrow + 1,ind(jcol)),
     :                jcol = 1,ncol)
              write (3, 70) jlev(jrow)+0.5, inrow,
     :            (zmat(jrow,ind(jcol)) + zmat(jrow + 1,ind(jcol)),
     :                jcol = 1,ncol)
            end if
          end if
*  if this row is an integer multiple of 10, add a blank line
*  to output on unit 6 only
          jj = jrow / 10
          if (10 * jj .eq. jrow .and. jrow .ne. nlevop) then
            if (iprint) write (6, 80)
          end if
80        format(1h )
95      continue
*  all rows have been printed out for this group of columns
*  move to the next group of columns
*  jlow becomes old jhigh
*  jhigh is set equal to jhigh plus skip distance
        if (jhigh .eq. nlevop) go to 160
        jlow = jhigh + 1
        jhigh = min ( (jlow + iskip - 1), isize)
150   continue
*  entire matrix has been printed, return
160   return
      end
*  ------------------------------------------------------------------
      subroutine prsgpi (fname, a)
*  subroutine to write out selected integral cross sections
*  from file {fname1}.ics for sigma - pi transitions
*  author:  millard alexander
*  latest revision date:  5-apr-2004 by mha
*  ------------------------------------------------------------------
*  variables in call list:
*    zmat:    on return:  contains the nlevop x nlevop matrix of integral
*                         cross sections
*    jlev:   rotational angular momenta of each energetically open level
*    elev:   energy (in hartree) of each energetically open level
*    inlev:  additonal quantum index for each energetically open level
*    jfirst:  initial value of total angular momentum
*    jfinal:  final value of total angular momentum
*    jtotd:   step size for total angular momentum
*    ipos:    if .true., 132 column printer
*             if .false., 80 column printer
*    csflag:  if .true. coupled-states calculation
*             if .false., close-coupled calculation
*    flaghf:  if .true., then system with half-integer spin
*             if .false., then system with integer spin
*    twomol:  if .true., then molecule-molecule collision
*             if .false., then atom-molecule or molecule-surface collision
*    nlevop:  number of energetically distinct levels in channel basis which
*             are open asymptotically
*    numin:   initial value of the coupled-states projection index
*    numax:   final value of the coupled-states projection index
*    nud:     step in coupled-states projection index
*    jlpar:   parity of channels
*               eps * (-1)**(j+l-jtot) = jlpar
*             if jlpar = 0, then integral cross sections include contributions
*             of both parities
*    note!!!   if flaghf = .true.( see below), then the true values
*    of the rotational quantum numbers, the total angular momentum,
*    and the coupled-states projection index are equal to the values
*    stored in jlev, jtot, jfirst, jfinal, nu, numin, and numax plus 1/2
*    flagsu:    if .true., then molecule-surface collisons
*    xmu:       collision reduced mass in (c12) atomic mass units
*    econv:     conversion factor from cm-1 to hartrees
*  ------------------------------------------------------------------
*      implicit none
      use constants
      implicit double precision (a-h,o-z)
      character*(*) fname
      character*20 cdate
      character*40 xnam1, xnam2
      character*80 line
cmha
      character*3 stat
      character*12 accs
cmha
      logical csflag, flaghf, iprint, ipos, flagsu, twomol, existf,
     :        openfl,lpar(17),
     :        airyfl, airypr, bastst, batch, chlist, ihomo, eprint
      include "common/parpot"
      common /colpar/ airyfl, airypr, bastst, batch, chlist,
     :                csflag, flaghf, flagsu, ihomo, ipos,lpar
      common /coz/ zmat(1)
      common /cow/ scmat(1)
      common /coamat/ zbuf(1)
      common /cojhld/ jlev(1)
      common /coisc1/ inlev(1)
cmha
      common /coisc2/ jpoint(5)
cmha
      common /colq/ ipoint(1)
      common /cosc1/ elev(1)
      common /cosout/ nnout, jout(21)
      common /coiout/ niout, indout(5)
      common /coener/ energ(1)
      common /coselb/ ibasty
      dimension  a(4)
*  input parameters
cmha
      iprint=.true.
      eprint=.false.
      iaver = 0
cmha
cABER
      if (a(1) .lt. 0.0d0) iprint = .false.
cmha
      if (a(1) .gt. 0.d0) eprint = .true.
      if (a(2) .lt. 0.5) iaver = 0
      if (a(2).gt. 0.d0) iaver=1
      if (a(2) .gt. 1.5) iaver = 2
cABER
      ienerg = a(3) + 0.1
      xthresh=a(4)
      iener=energ(1)
      if (ienerg .le. 0) ienerg=1
*  open integral cross section file
      call gennam (xnam1, fname, ienerg, 'ics', lenx)
      inquire (file = xnam1, exist = existf)
      if (.not. existf) then
        write (6, 10) xnam1(1:lenx)
10      format(/' integral cross section file ',(a),' not found',/)
        return
      end if
      open (1, file = xnam1,status = 'old')
*  open output file for integral cross sections
      call gennam (xnam2, fname, ienerg, 'xsc', lenx)
* inquire file specifications
      inquire(file=xnam2, exist=existf, opened=openfl)
cmha .xsc file is appended if it already exists
      if (.not. openfl) then
        accs='sequential'
        if (existf) then
          stat='old'
* make sure sequential formatted files are appended not overwritten
cstart unix-hp unix-dec unix-iris unix-sun unix-ifort unix-pgi
          accs='append'
cend

        else
          stat='new'
        end if
        open (3, file = xnam2, status = stat,access = accs)
      endif
cmha
cmha (print out message only if no other print out)
      if (.not.iprint) write(6,20) xnam2
20    format(
     :  ' PRINTING SELECTED SIGMA-PI CROSS SECTIONS; OUTPUT IN FILE ',
     :        (a))
cmha
*      call version(3)
      read (1, 40) cdate
40    format (1x, a)
      read (1, 40) label
      read (1, 40) potnam
*  print job information
      write (3, 50) xnam1, cdate, label, potnam
      if (iprint) write (6, 50) xnam1, cdate, label, potnam
50    format(
cmha
     :  /'% INTEGRAL SIGMA-PI CROSS SECTIONS READ FROM FILE ',(a)/
cmha
     +        '% WRITTEN:   ',(a)/
     +        '% LABEL:     ',(a)/,
     +        '% POT NAME:  ',(a) )
      read (1, 60) ener, xmu
60    format (f10.3, f15.11)
      read (1, 70) csflag, flaghf, flagsu, twomol, ihomo
70    format (5l3)
      read (1, 80) jfirst, jfinal, jtotd, numin, numax,
     :             nud, jlpar, isa
80    format (24i5)
      if (ipos) then
        read (1, 80)  nlevel, nlevop
        read (1, 80) (jlev(i), inlev(i), i = 1, nlevel)
        read (1, 90) (elev(i), i = 1, nlevel)
90      format (8f16.9)
      else
        read (1, *)  nlevel, nlevop
        read (1, *) (jlev(i), inlev(i), i = 1, nlevel)
        read (1, *) (elev(i), i = 1, nlevel)
      endif
caber_begin
*  zero out zmat
      call dset(nlevop*nlevop,0.d0, zmat,1)
caber_end
*  read in matrix of cross sections, column by column
99    do 95  i = 1, nlevop
*  jbegin and jend point to first element and last element of column i of
*  matrix packed column by column
        jbegin = (i - 1) * nlevop + 1
        jend = i * nlevop
c            write (nxfile, 90) (zmat(j,i), j = 1, nlevop)
*  chnage format of reac statement (p. dagdigian, 7-mar-2012)
*        read (1, 90) (zbuf(j), j = jbegin, jend)
        read (1, *) (zbuf(j), j = jbegin, jend)
        do 96 j=jbegin,jend
          zmat(j)=zmat(j)+zbuf(j)
96      continue
95    continue
      read(1,'(a)',end=999) line
      if (line(1:14).eq.' ** RESTART **') then
        read(1,40) cdate
        read(1,40) label
        read (1,40) potnam
        read (1, 60) ener, xmu
        read (1, 70) csflag, flaghf, flagsu, twomol, ihomo
        read (1, 80) idum, jfinal, jtotd, numin, numax,
     :               jlpar, isa
        read (1, 80)  nlevel, nlevop
        read (1, 80) (jlev(i), inlev(i), i = 1, nlevel)
        read (1, 90) (elev(i), i = 1, nlevel)
        goto 99
      end if
999   if (.not. flagsu) then
        if (.not. flaghf) then
          if (iprint) write (6, 100) ienerg, xmu, ener, jlpar,
     :                                jfirst, jfinal, jtotd
          write (3, 100) ienerg, xmu, ener, jlpar, jfirst,
     :                   jfinal, jtotd
100       format ('% ** INTEGRAL CROSS SECTIONS; IEN=', i2,' **',
     :            /'%    RMU=', f9.4, '  E=', f9.2,'  JLPAR=', i2,
     :             '  JTOT-1=', i3,
     :             '  JTOT-2=', i4,'  JTOT-D=', i3)
        else
          if (iprint) write (6, 105) ienerg, xmu, ener, jlpar,
     :                      (jfirst+0.5), (jfinal+0.5), jtotd
          write (3, 105) ienerg, xmu, ener, jlpar, (jfirst+0.5),
     :                   (jfinal+0.5), jtotd
105       format ('% ** INTEGRAL CROSS SECTIONS; IEN=', i2,' **',
     :            /'%    RMU=', f9.4, '  E=', f9.2,'  JLPAR=', i2,
     :             '  JTOT-1=', f5.1,
     :             '  JTOT-2=', f6.1,'  JTOT-D=', i3)
        end if
      else
        if (.not. flaghf) then
          if (iprint) write (6, 110) ienerg, xmu, ener, numin,
     :                               numax, nud
          write (3, 110) ienerg, xmu, ener, numin, numax, nud
110       format ('% ** SUMMED DEGENERACY AVERAGED TRANSITION',
     :             ' PROBABILITIES;  IEN=', i2,' **',
     :            /'    RMU=', f9.4, '  E=', f8.2,'  M-MIN=', i3,
     :             '  M-MAX=', i4, '  M-STEP=', i2)
        else
          if (iprint) write (6, 115) ienerg, xmu, ener, (numin+0.5),
     :                               (numax+0.5), nud
          write (3, 115) ienerg, xmu, ener, (numin+0.5),
     :                   (numax+0.5), nud
115       format ('% ** SUMMED DEGENERACY AVERAGED TRANSITION',
     :             ' PROBABILITIES;  IEN=', i2,' **',
     :            /'    RMU=', f9.4, '  E=', f8.2,'  M-MIN=', f5.1,
     :             '  M-MAX=', f6.1, '  M-STEP=', i2)
        end if
      end if
      if (.not. csflag) then
        if (jlpar .eq. 0) then
          write (3, 120)
          if (iprint) write (6, 120)
120       format ('$ ** CC CALCULATION, BOTH PARITIES **')
        else
          write (3, 125) jlpar
          if (iprint) write (6, 125) jlpar
125       format ('% ** CC CALCULATION, JLPAR=', i2, ' **')
        end if
      else
        if (.not. flaghf) then
          write (3, 130) numin, numax
          if (iprint) write (6, 130) numin, numax, nud
130       format ('% ** CS CALCULATION, NUMIN=', i2,', NUMAX=',
     :               i2,', NUD=',i2,' **')
        else
          write (3, 140) numin + 0.5, numax + 0.5
          if (iprint) write (6, 140) numin+0.5, numax+ 0.5,nud
140       format ('% ** CS CALCULATION, NUMIN=', f4.1, ', NUMAX=',
     :            f4.1,', NUD=',i2,' **')
        end if
      end if
cmha
c print out of level list only for iprint .ge. 1
      if (eprint) then
        if (.not. twomol) then
          if (.not. flagsu) then
            if (iprint) write (6, 145)
            write (3, 145)
145         format
     :       (/'% LEVEL LIST FOR INTEGRAL CROSS SECTIONS',
     :        /'   N   J  INDEX  EINT(CM-1)',/)
          else
            if (iprint) write (6, 150)
            write (3, 150)
150         format
     :       (/'% LEVEL LIST FOR DEGENERACY AVERAGED',
     :         ' TRANSITION PROBABILITIES',
     :        /'   N   J  INDEX  EINT(CM-1)',/)
          end if
          do 170  i = 1, nlevop
            if (.not. flaghf) then
              if (iprint)
     :        write (6, 160) i, jlev(i), inlev(i), elev(i) * econv
              write (3, 160) i, jlev(i), inlev(i), elev(i) * econv
160           format (i4, i5, i6, f11.3)
            else
              if (iprint) write (6, 165) i, (jlev(i)+0.5), inlev(i),
     :                       elev(i) * econv
              write (3, 165) i, (jlev(i)+0.5), inlev(i),
     :                       elev(i) * econv
165           format (1x,'%',i4, f5.1, i6, f11.3)
            end if
170       continue
        else
          if (iprint) write (6, 175)
          write (3, 175)
175         format (/'% LEVEL LIST FOR INTEGRAL CROSS SECTIONS',
     :              /'%   N   J1   J2  INDEX  EINT(CM-1)'/)
          do 190  i = 1, nlevop
            jj2 = mod( jlev(i), 10)
            jj1 = jlev(i) / 10
            if (iprint)
     :      write (6, 180) i, jj1, jj2, inlev(i), elev(i) * econv
            write (3, 180) i, jj1, jj2, inlev(i), elev(i) * econv
180         format ('%',i4, 2i5, i6, f11.3)
190       continue
        end if
      endif
cmha
*  now sum and average over positive and negative values of index
*     check that number of levels with index negative equals number of levels
*     with index positive, if not abort
*     this test works by (1) checking that no index equals zero
*     and, if so, (2) adding all the values of the index.  this sum
*     should equal zero if there are as many levels with index negative
*     as index positive
      if (iaver .gt. 0) then
cABER
cmha
        if (iaver.gt.1.d0.and.ihomo) then
cmha
          write(6,194)
194       format(' *** HOMONUCLEAR MOLECULE : AVERAGING',
     :            ' MAKES NO SENSE ! ***')
          iaver=1
        end if
cABER
        isum = 0
        do  200 i = 1, nlevop
          if (inlev(i) .eq. 0) then
            write (6, 195) i
            write (6, 195) i
195         format(/' *** INLEV(',i4,')=0;',
     :               ' AVERAGING MAY NOT WORK ***')
          else if (inlev(i) .ne. 0) then
            isum = isum + inlev(i)
          end if
200     continue
        if (isum .ne. 0 .and. .not.flaghf) then
          write (6, 230) isum
          write (3, 230) isum
230       format ('% *** SUM OF INDICES =',i4,
     :               ' AVERAGING MAY NOT WORK ***')
        end if
        if  (iaver .eq. 2) then
          if (iprint) write (6, 245)
          write (3, 245)
245       format
     :      ('% ** CROSS SECTIONS SUMMED AND AVERAGED OVER INDEX **')
          call aver2 (zmat, scmat, nlevop)
        else if (iaver .eq. 1) then
          if (iprint) write (6, 250)
          write (3, 250)
250       format
     :      ('% ** CROSS SECTIONS SUMMED FINAL STATE INDEX **')
        end if
      end if
cmha
c 4 lines moved from here 2/27/92
cmha
*  find all rows of cross sections matrix for which initial rotational
*  quantum number is equal to one of the values of jout
      isize = 0
      insize=0
      nout = abs (nnout)
      do  280  iout = 1, nout
        jtemp = jout(iout)
        do 270 n = 1, nlevop
          if (jlev(n) .eq. jtemp) then
            isize = isize + 1
            jpoint(isize) = n
cmha
            ipoint(isize) = n
cmha
          end if
270     continue
280   continue
cmha
      if (niout .gt. 0) then
        nout=abs(niout)
        do  282  iout = 1, nout
          indtemp=indout(iout)
          do 281 n = 1, isize
            if (inlev(jpoint(n)) .eq. indtemp) then
              insize = insize + 1
              ipoint(insize) = jpoint(n)
            endif
281       continue
282     continue
        isize=insize
      endif
cmha
*  isize is the number of cross sections to be printed
      if (isize .eq. 0) then
*  here if no initial states found
        write (6, 283)
        write (3, 283)
283     format ('% ** NO INITIAL STATES FOUND; ABORT')
      else
        write (3, 290) xthresh
        if (iprint) write (6, 290) xthresh
290     format (/'% ** COLUMN HEADINGS ARE INITIAL STATES, ROW',
     :        ' HEADINGS ARE FINAL STATES **',
     :      /'%      CROSS SECTION PRINT THRESHOLD=',1pd8.1)
        if (iener.lt.10) then
            write (3,295) iener
        elseif (iener.lt.100) then
            write (3,296) iener
        elseif (iener.lt.1000) then
            write (3,297) iener
        elseif (iener.lt.10000) then
            write (3,298) iener
        elseif (iener.lt.100000) then
            write (3,299) iener
        endif
295     format('x',i1,'=[')
296     format('x',i2,'=[')
297     format('x',i3,'=[')
298     format('x',i4,'=[')
299     format('x',i5,'=[')
*  now print out desired columns of cross section matrix
        call xscpr2(zmat, xthresh, nlevop, isize, iaver, iprint, isa)
        write (3,300)
300     format('];')
      endif
      close (1)
      close (3)
      return
      end
      subroutine aver2 (zmat, scmat, n)
      implicit double precision (a-h,o-z)
*  subroutine to sum and average cross section matrix over positive
*  and negative values of index
      integer i, ind, index, j, jnd, n, nn, inlev, jlev
*     real scmat, zmat
      dimension zmat(1), scmat(n,n)
      common /cojhld/ jlev(1)
      common /coisc1/ inlev(1)
*  first copy cross section matrix into scmat
      call matmov (zmat, scmat, n, n, n, n)
      nn = n / 2
      index = 0
      do  30  i = 1, nn
        ind = 2 * i
        do  20  j = 1, nn
          jnd = 2 * j
          index = index + 1
          zmat(index) = ( scmat(ind - 1, jnd - 1) +
     :                   scmat(ind - 1, jnd) +
     :                   scmat(ind, jnd - 1) +
     :                   scmat(ind, jnd) ) * 0.5d0
20      continue
         jlev (i) = jlev (ind - 1)
        inlev(i) = iabs(inlev(ind - 1))
30    continue
*  because of indexing zmat need to be transposed
      call transp (zmat, nn, nn)
*  reduce the size of the cross section matrix
      n = nn
      return
      end
* ------------------------------------------------------
      subroutine xscpr2 (zmat, xthresh, nlevop, isize, iaver,
     :                   iprint,isa)
      implicit double precision (a-h,o-z)

*  current revision date:  10-oct-2001 by ab
*  subroutine to print out specified columns of cross section matrix
*  if iaver = 1, then nth and (n+1) st rows are added before printing
      integer i, isize, iskip, j, jcol, jhigh, jj, jlow, jmax,
     :        jrow, ncol, nlevop, iaver
      integer ind, inlev, ipoint, jlev, isa
      logical airyfl, airypr, bastst, batch, chlist, csflag, flaghf,
     :                flagsu, ihomo, ipos, iprint,lpar(17)
      dimension zmat(nlevop,nlevop), ind(50)
      common /colpar/ airyfl, airypr, bastst, batch, chlist,
     :                csflag, flaghf, flagsu, ihomo, ipos,lpar
      common /cojhld/ jlev(1)
      common /coisc1/ inlev(1)
      common /colq/ ipoint(1)
      common /cosc1/ elev(1)
*  first transpose cross section matrix so that initial states are
*  columns and final states are rows
      call transp (zmat, nlevop, nlevop)
*  if 132 line printer, then 13 columns of matrix are printed simultaneously
*  if  80 line printer, then  7 columns of matrix are printed simultaneously
      iskip = 13
      if (.not. ipos) iskip = 7
*  jamx is the total number of the groups of columns to be printed
      if (mod(isize, iskip) .eq. 0) then
        jmax = isize / iskip
      else
        jmax = isize / iskip + 1
      end if
      jlow = 1
      jhigh = min (iskip, isize)
*  loop ovr columns by groups of 6 or 10
      do  150   j = 1, jmax
*  ncol is the number of columns to be printed
        ncol = jhigh - jlow + 1
        do  10   jj = 1, ncol
*  the array ind contains the index of each column which will be printed
          ind(jj) = ipoint(jj - 1 + jlow)
10      continue
*  write as a heading the column index of each column to be printed
        if (.not. flaghf) then
          if (iprint) write (6, 15) ( jlev(ind(i)),i = 1,ncol)
          write (3, 15) ( jlev(ind(i)), i = 1,ncol)
15        format (/11x,'J= ', i4, 2x, 12 (2x, i5, 2x) )
        else
          if (iprint) write (6, 30) ( jlev(ind(i))+0.5,
     :                                 i = 1,ncol)
          write (3, 30) ( jlev(ind(i))+0.5, i = 1,ncol)
30        format (/'%',11x,'J= ', f4.1, 2x, 12 (2x, f5.1, 2x) )
        end if
        if (iprint) write (6, 40) ( inlev(ind(i)), i = 1,ncol)
        write (3, 40) ( inlev(ind(i)), i = 1,ncol)
40      format ('%   J    I | I=', i4, 2x, 12 (2x, i5, 2x))
        if (iprint) write (6, 50)
        write (3, 50)
50      format (1h )
*  now loop through the rows of the matrix, which will be printed out
*  in groups of 10 with a blank line in between
*  lopp over the rows of the matrix which will be printed
        nlev = 0
*  nlev is the n quantum number for sigma states
        do 95 jrow = 1 , nlevop
*  inrow holds additional quantum index for this row
          inrow = inlev(jrow)
cABER  summing will be done for HETERONUCLEAR case (PI and SIGMA)
cABER  and HOMONUCLEAR case (only SIGMA)
          if (iaver.eq.1.and.
     :        (.not.ihomo.or.aint(iabs(inrow)/100.0).eq.3)) then
cABER  find beginning of a new vibrational branch for SIGMA ( ==> nlev=0 )
            if (aint(iabs(inrow)/100.0).eq.3) then
              if (abs(inlev(jrow))-abs(inlev(jrow-1)).ne.0) nlev = 0
cABER  initialize SWITCH parameter depending on symmetry (ISA) and "parity"
cABER  of row index (jrow)
              if (nlev.eq.0) then
                if (isa.eq.1.or.isa.eq.0) then
                  if (mod(jrow,2).eq.0) iswtch=0
                  if (mod(jrow,2).ne.0) iswtch=1
                else if (isa.eq.-1) then
                  if (mod(jrow,2).eq.0) iswtch=1
                  if (mod(jrow,2).ne.0) iswtch=0
                end if
                if (isa.ne.-1) go to 94
              end if
            end if
cABER  if sum over final states desired, then sum and skip rows (depending on
cABER  SWITCH parameter)
            if (iswtch .eq. 0) then
              if (mod (jrow,2) .eq. 0) go to 95
            else if (iswtch .eq. 1) then
              if (mod (jrow,2) .ne. 0) go to 95
            end if
cABER  now write out the row index followed by the desired matrix elements
cABER  here for PI states
94          if (aint(iabs(inrow)/100.0).ne.3) then
              if (.not. flaghf) then
                if (iprint)
     :          write (6, 60) jlev(jrow),inrow,
     :            (zmat(jrow,ind(jcol)) + zmat(jrow + 1,ind(jcol)),
     :                jcol = 1,ncol)
                write (3, 60) jlev(jrow),inrow,
     :            (zmat(jrow,ind(jcol)) + zmat(jrow + 1,ind(jcol)),
     :                jcol = 1,ncol)
 60              format (i5, i5, 2x, 13 (1pe9.2,1x) )
              else
                if (iprint)
     :            write (6, 70) jlev(jrow)+0.5, inrow,
     :            (zmat(jrow,ind(jcol)) + zmat(jrow + 1,ind(jcol)),
     :                jcol = 1,ncol)
                write (3, 70) jlev(jrow)+0.5, inrow,
     :            (zmat(jrow,ind(jcol)) + zmat(jrow + 1,ind(jcol)),
     :                jcol = 1,ncol)
 70              format (f5.1, i5, 2x, 13 (1pe9.2,1x) )
              end if
cABER  here for SIGMA states
            else if (aint(iabs(inrow)/100.0).eq.3) then
              if (nlev.eq.0.and.isa.ne.-1) then
                if (iprint)
     :            write (6, 60) jlev(jrow),inrow,
     :            ( zmat(jrow,ind(jcol)), jcol = 1,ncol)
                write (3, 60) jlev(jrow),inrow,
     :           ( zmat(jrow,ind(jcol)), jcol = 1,ncol)
                nlev = 1
              else
                if (iprint)
     :          write (6, 60) jlev(jrow)+1,inrow,
     :            (zmat(jrow,ind(jcol)) + zmat(jrow + 1,ind(jcol)),
     :                jcol = 1,ncol)
                write (3, 60) jlev(jrow)+1,inrow,
     :            (zmat(jrow,ind(jcol)) + zmat(jrow + 1,ind(jcol)),
     :                jcol = 1,ncol)
                nlev = nlev + 1
              end if
            end if
cABER  end of summing over final states encountered
          else
            zmax=0d0
            do 75 jcol=1,ncol
               zz=zmat(jrow,ind(jcol))
               if (zz.gt.zmax) zmax=zz
75          continue
            if (zmax.gt.xthresh) then
               if (.not. flaghf) then
                 if (iprint)
     :             write (6, 60) jlev(jrow),inrow,
     :               ( zmat(jrow,ind(jcol)), jcol = 1,ncol)
                 write (3, 60) jlev(jrow),inrow,
     :             ( zmat(jrow,ind(jcol)), jcol = 1,ncol)
               else
                 if (iprint)
     :             write (6, 70) jlev(jrow)+0.5, inrow,
     :               ( zmat(jrow,ind(jcol)), jcol = 1,ncol)
                 write (3, 70) jlev(jrow)+0.5, inrow,
     :             ( zmat(jrow,ind(jcol)), jcol = 1,ncol)
               end if
            endif
          end if
          jj = jrow / 10
          if (10 * jj .eq. jrow .and. jrow .ne. nlevop) then
            if (iprint) write (6, 80)
          end if
80        format(1h )
95      continue
        if (jhigh .eq. nlevop) go to 160
        jlow = jhigh + 1
        jhigh = min ( (jlow + iskip - 1), isize)
150   continue
160   return
      end
