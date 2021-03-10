*system: NH3-H2 using canonical Valiron expansion of potential
      include "common/syusr"
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(45)
      include "common/parpot"
      potnam='RIST/VALIRON NH3-H2'
       print *, potnam
1      print *, ' r (bohr)'
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' v',/,45(1pe16.8))
      goto 1
99    end

      include "common/syusr"
      include "common/bausr"
      include "common/ground"
      subroutine loapot(iunit, filnam)
*  variables in common block /cobspt/
 
*    lammin:   array containing minimum value of lambda for each term
*    lammax:   array containing maximum value of lambda for each term
*    mproj:    array containing the order of the reduced rotation matrix
*              elements for each term. here, lammin is greater to mproj.
*  variables in common block /cobsptln/
*    lam2:     array containing the ordre of the reduced rotation matrix
*              elements for each term. in case of homonuclear molecule
*              is even.
*    m2proj:   array containing the order of the reduced rotation
*              matrix elements for each term.
*              here, lammin and lam2 are greater than m2proj .
*  variables in common block /cosysi/
 
*    nscod:    total number of variable names which are passed to hinput
*              nscod must equal isrcod + isicod + 3
*    isicod:   total number of integer system dependent variables
*    nterm:    number of different associated legendre terms in
*              expansion of potential
*    numpot:   the number of the potential used, this variable is passed
*              to the pot subroutine
*    ipotsy:   cylindrical symmetry of potential.  only terms with mu equal to
*              an integral multiple of ipotsy can be included in the potential.
*              example:  for nh3, ipotsy = 3
*    iop:      ortho/para label for molecular states. if ihomo=.true. then only
*              para states will be included if iop=1 and only ortho states if
*              iop=-1
*    ninv:     number of inversion doublets included
*              if ninv = +1, only + inversion levels included
*              if ninv = -1, only - inversion levels included
*              if ninv = 2, both inversion levels included
*    kmax:     the maximum projection quantum number included
*    jmax0:    the maximum rotational angular momenta for the k=0 stack
*    jmax1:    the maximum rotational angular momenta for the k=1 stack
*    jmax2:    the maximum rotational angular momenta for the k=2 stack
*    jmax3:    the maximum rotational angular momenta for the k=3 stack
*    jmax4:    the maximum rotational angular momenta for the k=4 stack
*    jmax5:    the maximum rotational angular momenta for the k=5 stack
*    jmax6:    the maximum rotational angular momenta for the k=6 stack
*    ipotsy2:  symmetry of potential. if linear molecule is homonuclear
*              then ipotsy=2 and only terms with lambda2  even can be
*              included in the potential,else ipotsy=1.
*    j2max:    the maximum rotational angular momentum for linear
*              molecule
*    j2min:    the minimum rotational angular momentum for linear
*              molecule
 
*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
 
      logical readpt
      include "common/parbas"
      include "common/parpot"
      common /cosysi/ nscode, isicod, nterm, numpot, ipotsy, iop, ninv,
     :                kmax, jmax0, jmax1, jmax2, jmax3, jmax4, jmax5,
     :                jmax6, ipotsy2, j2max, j2min
 
*  variable in common block /conlam/
*    nlammx:      the total number of angular coupling terms
*    nlam:    the maximum number of angular coupling terms allowed
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx (not used here)
 
      common /conlam/ nlam, nlammx, lamnum(1)
 
c  variables in common block /cutoff/
c  ----------------------------------
 
c      special cutoff : to check if the anisotropic are really significant
c      for the collision at long range. if they had not a negligible effect,
c      additional scf and disp points should be computed beyond 10 au.
 
c    rcut:       cutoff distance (au)
c    drcut:      cutoff width (au)
c    llcut:      logical flag, enables the cutoff of all the vij's terms
c                in the scf expansion excepted the v00 term.
c                the cutoff function is an hyperbolic tangent in order to
c                reach 1 and 0 faster than any polynomials.
c                the cutoff is applied on the original data before the spline
c                interpolation. thus the cutoff width should be much larger
c                than the grid step around rcut to prevent spurious spline
c                oscillations.
 
      logical llcut
      common /cutoff/ rcut,drcut,llcut
 
c  variables in common block /fisurf/
c  ----------------------------------
 
c    conv and econv : set the unit conversion for distances and energies.
c                these conversions are designed to match the units in use
c                in the collision code. all input data on units ifile and
c                ifile 2 must be expressed in au for distances and in cm-1
c                for distances. output data and listed values are converted
c                according to conv and econv.
c    conv:       atomic unit of distance expressed in user units.
c                set conv = 1.0 for atomic units
c                set conv = 0.529177 for angstroms
c    econv:      user energy unit expressed in cm-1.
c                set econv = 1.0 for cm-1.
c                set econv = 2 * 109737.312 for hartrees.
c    isurf:      select the potential surface
c                  1   uses scf data only
c                  2   uses disp data only
c                  3   uses both scf and disp ( default )
c    lsurf:      previous value of isurf. setting lsurf to 0 forces reading
c                again the data on ifile.
c    damped:     logical flag. applies only for paul wormer's data.
c                enforces a uniform damping on all terms.
c                the corresponding v00 term is read on ifile after the scf
c                and disp information with iflag set to 3.
 
*   ===========================>>>> check what's this damping is....
 
      logical damped
      common /fisurf/ conv, econv, isurf, lsurf, damped
 
 
c  variables in common blocks /fiunit/ and /finame/
c  ------------------------------------------------
 
c    iwrite       logical unit for listing output (on file).
c                 logical unit for display is assumed equal to 6.
c    ifile        logical unit to read the potential expansion.
c                 distances are read in bohrs and energies in cm-1.
c    ifile2       logical unit to read the additional potential data
c                 using the namelist /fit/.
c    datafl:      name of data file (pot. expansion) read on unit ifile
c    datanm:      name of namelist data file ( set with pot= )
 
      integer iwrite, ifile, ifile2
      common /fiunit/ iwrite,ifile,ifile2
      character*60 datafl
      character*40 datanm
      common /finame/ datafl, datanm
 
 
c  variables in common block /fiextr/
c  ----------------------------------
 
c     in order to avoid numerical problems, a simple extrapolation at
c     short range and at long range should be provided.
c      .   short range part of v00 is fitted to exponential
c      .   short range part of other vij's is damped exponentially
c      .   long range part of all vij's is damped exponentially
c    rmin:        short range radius (au)
c    rmax:        long range radius (au)
c    umin:        short range exponant (au-1)
c    umax:        long range exponant (au-1)
 
      double precision rmin, rmax, umin, umax
      common /fiextr/ rmin,rmax,umin,umax
 
      namelist /fit/  rcut, drcut, llcut
     1               ,conv, econv, isurf, damped
     1               ,rmin, rmax, umin, umax
     1               ,datafl,nterm,lammax
 
      character*(*) filnam
      parameter (nvmx = 45)
      dimension ivij(nvmx), jvij(nvmx), i2vij(nvmx), j2vij(nvmx),
     :          lambda(nvmx), mu(nvmx),
     :          lambda2(nvmx), mu2(nvmx)
      common /coloapot/ s4pi, ivv(nvmx), nvv
      potnam='RIST/VALIRON NH3-H2'
 
      nterm = 3
      do j=1, 4
        do i=1, 3
          k=3*(j-1)+i
          mproj(k) = 3*(i-1)
          if(j.le.2) m2proj(k) = 0
          if(j.eq.3) m2proj(k) = 1
          if(j.eq.4) m2proj(k) = 2
          if(j.eq.1) lam2(k) = 0
          if(j.ge.2) lam2(k) = 2
          lammax(k) = 6
          lammin(k) = max(mproj(k), m2proj(k))
        enddo
      enddo
      lammin(1) = 1
      if (filnam .eq. ' ') return
 
      datanm = filnam
 
200   continue
      iwrite = 9
      ifile  = 60
      ifile2 = 61
 
c     conv = 0.529177d0
      conv = 1.d0
c     econv = 1.0d0
      econv = 109737.312d0 * 2
 
      isurf = 3
      damped = .false.
 
      rmin = 4.0d0
      umin = 20.0d0
      rmax = 50.d0
      umax = 1.0d0
 
      llcut = .false.
      rcut  = 12.d0
      drcut = 1.5d0
 
c  read the namelist /fit/ on unit ifile2
      open (unit=ifile2, file=datanm, status='old')
      read(ifile2,fit)
      close (unit=ifile2)

*  calculate total number of anisotropic terms
      nlam = 0
      do 135  i = 1, nterm
        if (mproj(i) .gt. lammin(i) ) then
          write (6, 130)  mproj(i), lammin(i), i
          write (9, 130)  mproj(i), lammin(i), i
130       format (' *** mproj=',i2,' > lammin=',i2,
     :            ' for term',i2,'; abort ***')
          call abort
        end if
        if (m2proj(i) .gt. lammin(i) ) then
          write (6, 131)  m2proj(i), lammin(i),i
          write (9, 131)  m2proj(i), lammin(i), i
131       format (' *** m2proj=',i2,' > lammin=',i2,
     :            ' for term',i2,'; abort ***')
          call abort
        end if
        if (m2proj(i) .gt. lam2(i) ) then
          write (6, 132)  m2proj(i), lam2(i),i
          write (9, 132)  m2proj(i), lam2(i),i
132       format (' *** m2proj=',i2,' > lam2=',i2,
     :            ' for term',i2,'; abort ***')
          call abort
        endif
        nlam = nlam + lammax(i) - lammin(i) + 1
135   continue
      if (nlammx .lt. nlam) then
        write (6, 140) nlam, nlammx
        write (9, 140) nlam, nlammx
140     format (' ** total number of anisotropic terms=', i2,
     :          ' .gt. nlam=', i2,'; abort')
        call abort
      end if
c  convert to units of distance according to conv
      rcut = rcut * conv
      drcut = drcut * conv
      rmin = rmin * conv
      rmax = rmax * conv
      umin = umin / conv
      umax = umax / conv
 
c  display and list parameters
      write(6,1000) datanm, datafl
      write(6,1001) iwrite,ifile,ifile2
      write(6,1002) conv,econv,isurf,damped
      write(6,1003) rmin,umin,rmax,umax
      write(6,1004) rcut,drcut,llcut
 
      write(iwrite,1000) datanm, datafl
      write(iwrite,1001) iwrite,ifile,ifile2
      write(iwrite,1002) conv,econv,isurf,damped
      write(iwrite,1003) rmin,umin,rmax,umax
      write(iwrite,1004) rcut,drcut,llcut
 
1000  format(/'0fi2par -- initialise potential surface parameters'
     +       /'0 namelist file  ',a/
     +       /' potential data  ',a)
1001  format('0    iwrite     ifile    ifile2'/,1x,4i10)
1002  format('0      conv     econv      isurf    damped'
     +      /,1x,f12.6,f12.3,i10,l10)
1003  format('0      rmin      umin      rmax      umax'/,1x,4f10.6)
1004  format('0      rcut     drcut     llcut'/,1x,2f10.6,l10)
 
*  force reading the potential expansion ( with lsurf = 0 )
      lsurf = 0
      call vinit (1, rinit, v)
* match the order of potential terms in file potfil and in basis
* define in ivij, jvij, ivij2, jvij2 the order of potfil
      kv = 0
      do 210 i2v = 0, 2, 2
      do 210 j2v = 0, i2v
      do 210 iv = j2v, 6
      do 210 jv = 0, iv, 3
        kv = kv+1
        if (kv .gt. nvmx) then
          write(6,*) 'laopot -- dimension nvmx insuffisante'
          stop
        endif
        ivij(kv) = iv
        jvij(kv) = jv
        i2vij(kv) = i2v
        j2vij(kv) = j2v
210   continue
      nvij = kv
c      write(6,*)
c      write(6,*) 'termes du potentiel fitte :', nvij
c      write(6,211) (kv, ivij(kv), jvij(kv), kv=1, nvij)
211   format(1x, i6, '   v',2i1)
* define  in ivv a pointer from basis definition to the potfil
* definition.
      nvv = 0
      do 220 iterm = 1, nterm
        imu = mproj(iterm)
        ilam2 = lam2(iterm)
        imu2 = m2proj(iterm)
        do 220 ilamb = lammin(iterm), lammax(iterm)
          nvv = nvv+1
          lambda(nvv) = ilamb
          mu(nvv) = imu
          lambda2(nvv) = ilam2
          mu2(nvv) = imu2
          if (nvv .gt. nvmx) then
            write(6,*) 'sysdat -- probleme initialisation vecteur ivv'
            write(6,*) ilamb, imu, ilam2, imu2,
     :                 '   ilamb, imu, ilam2, imu2'
            write(6,*) nvv, '   nvv .gt. nvmx'
            stop
          endif
          do 230 kv = 1, nvij
            if ( (ivij(kv) .eq. ilamb) .and. (jvij(kv) .eq. imu)
     +      .and.(i2vij(kv) .eq. ilam2) .and. (j2vij(kv) .eq. imu2))
     +                                                         then
              ivv(nvv) = kv
              goto 240
            endif
230       continue
 
          write(6,*) 'sysdat -- probleme initialisation vecteur ivv'
          write(6,*) ilamb, imu, ilam2, imu2,
     +               '   ilamb, imu, ilam2, imu2   '
          write(6,*)
     +    'valeur non trouvee dans ivij, jviji, i2vij et j2vij'
          stop
240       continue
220   continue
      write(6,*)
      write(6,*) 'termes du potentiel vvl dans pot :', nvv
      write(6,250) (kvv, lambda(kvv), mu(kvv), lambda2(kvv),
     +              mu2(kvv), ivv(kvv),
     +              ivij(ivv(kvv)), jvij(ivv(kvv)),
     +              i2vij(ivv(kvv)), j2vij(ivv(kvv)), kvv=1, nvv)
250   format(1x, i6, '   vvl',4i2, '   as  ', i3, '  ( v',4i2,' )' )
 
      s4pir8 = sqrt ( 4.d0 * acos(-1.d0) )
      s4pi = s4pir8
c      write(6,*) s4pi, '    sqrt(4*pi)'
*
*  potential has been defined : force irpot=1
      irpot = 1
      write(6,260) datafl
260   format('Potential parameters successfully read from file:  ',
     :        a)
      return
      end
*  -----------------------------------------------------------------------
      subroutine pot( vv0, r)
*  -----------------------------------------------------------------------
*  -----------------------------------------------------------------------
*                .........   nh3 - para h2    .............
*  -----------------------------------------------------------------------
*  subroutine to calculate the r-dependent coefficients in the
*  collision of a symmetric top with a structureless target
*  in units of hartree for energy and bohr for distance
 
*  on return:
*  vv0 contains the isotropic term d1[0,0,0]*d2[0,0,0] in the potential
*  the coefficients for each angular term in the coupling potential
*  [ vvl(i) for i = 1, nlam ] are returned in common block /covvl/
*  in the following order:[l1,m1,l2,m2] associated with the angular
*  function drot(l1,m1,m2)*drot(l2,m2,0)
*    vvl(1) contains the [1,0,0,0] term
*    vvl(2) contains the [2,0,0,0] term
*    ....
*    vvl(lammax(1)) contains the [lammax(1),0,0,0] term
*    vvl(lammax(1)+1) contains the [lammin(2),3,0,0]
*    ....
*    vvl(lammax(1)+lammax(2) - 2) contains the [lammax(2),3,0,0]
*    etc... [l2,m2] sorted as: [0,0],[2,0],[2,1],[2,2]
*  variable in common block /conlam/
*    nlammx:      the total number of angular coupling terms
*    nlam:    the maximum number of angular coupling terms allowed
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx (not used here)
 
*  variable in common block /covvl/
*    vvl        array to store r-dependence of each angular term in the
*               potential
 
*  variable in common block /coloapot/
*  s4pi : facteur de normalisation du potentiel isotrope
*  ivv(1) corresponds to the [1,0,0,0] term
*  ivv(2) corresponds to the [2,0,0,0] term
*  ivv(3) corresponds to the [3,0,0,0] term
*  ivv(4) corresponds to the [3,3,0,0] term
*  [4,0,0,0],[4,3,0,0]...[6,6,0,0],[0,0,2,0]...[6,6,2,0],[1,0,2,1],...
*  [6,6,2,1],[2,0,2,2]...[6,6,2,2].
*  ivv(kv) correspond to the order in which the potential terms are read
*  nvv is the total number of angular coupling terms
 
      implicit double precision (a-h,o-z)
      integer ivv, nvv, nlam, nlammx, kv, kvv, nvmx
      parameter (nvmx=45)
      double precision  r, v, vv, vv0, vvl, s4pi
      common /conlam/ nlam, nlammx, lamnum(1)
      common /covvl/ vvl(nvmx)
      common /coloapot/ s4pi, ivv(nvmx), nvv
 
      if (nlam .ne. nvv) then
        write(6,*) 'entry pot --', nvv, nlam, '   nvv, expected nlam'
c        stop
      endif
      call vstar(1,r,v)
c
c  vv0 should be scaled as it must be the actual isotropic potential.
c
c  MAKE SURE that the scaling is consistent with the normalisation of the basis used.
c  SHOULD HAVE TO BE CHANGED when the new 4-index data are available.
c
      vv0 = v/s4pi
      do 400 kvv = 1, nvv
        kv = ivv(kvv)
        call vstar(kv,r,v)
        vvl(kvv) = v
400   continue
 
c      write(1,'(x,f12.6,1p5e20.12)') r, vv0, (vvl(kvv),kvv=1,min(nvv,4))
 
      return
      end
 
 
      subroutine vinit(i,r,v)
c***********************************************************************
c  file  vmol2/h2nh3sr - basically same as vmol/h2nh3ex except that
c
c      .   short range part of v00 is fitted to exponential
c      .   short range part of other vij's is damped exponentially
c      .   long range part of all vij's is damped exponentially
c      .   *** additional parameters are read on unit ifile2 ***
c
c  this potential is now well behaved from o to infinity.  pv 29/10/86
c
c  adapted to hibridon version 2, pv 23-jun-88
c  latest revision claire rist 18-jan-92
c***********************************************************************
      implicit double precision (a-h,o-z)
      integer i, isurf, lsurf, iwrite, ifile, ifile2
      double precision r, v, rcut, drcut, conv, econv, rmin, rmax, umin,
     :                 umax
      double precision bvs, avs, bv1s, av1s, av2s, bv2s, fcut, u, x,
     :                 rstart
      double precision v0, v1, v2
      logical llcut , damped
      common /cutoff/ rcut,drcut,llcut
      common /fisurf/ conv,econv,isurf,lsurf,damped
      common /fiunit/ iwrite,ifile,ifile2
      common /fiextr/ rmin,rmax,umin,umax
      save bvs, avs, bv1s, av1s, av2s, bv2s
 
c simple formula. second derivative non continuous.
      fcut (x) = (1.d0+u*x) * exp(-u*x)
 
      if (i.ne.1) return
      rstart = 10.d0
      call fi4vij (rstart,v0,1,1)
c
c set extrapolation for v00 at short range
c
      call fi4vij (rmin,v0,1,1)
      call fi4vij (rmin,v1,1,2)
      call fi4vij (rmin,v2,1,3)
      bvs = -v1/v0
      avs = v0 * exp(bvs*rmin)
      bv1s = -v2/v1
      av1s = v1 * exp(bv1s*rmin)
      av2s = -av1s * bv1s
      bv2s =  bv1s
      return
 
 
 
c=========================
      entry vstar(i,r,v)
c=========================
 
      if (r .lt. rmin) go to 12
      if (r .gt. rmax) go to 11
 10   call fi4vij (r,v,i,1)
      return
 11   u = umax
      call fi4vij (r,v0,i,1)
      v = fcut(r-rmax) * v0
      return
 12   if (i.gt.1) goto 13
        v = avs * exp(-bvs*r)
        return
 13   u = -umin
      call fi4vij (r,v0,i,1)
      v = fcut(r-rmin) * v0
      return
      end
        function fi2xyz (xx,yy,zz,nterms)
c
c       returns the value of the nh3-h2(j=0) potential for
c       cartesian coordinates of h2 (green conventions)
c       converts cartesian coordinates to polar (in radians)
c       and calls fi2rtp.
c       *** warning *** fi2xyz is much slower than fi2vij
c       original version                       p. valiron 14-june-84
c       modified arguments                     p. valiron 19-sept-87
c compute r, theta, phi as a function of xx, yy, zz
c
        implicit double precision (a-h,o-z)
        integer nterms
        double precision  xx, yy, zz, r, zzr, theta, phi, fi4rtp, fi2xyz
c
c no problem for r
        r     = sqrt (xx**2 + yy**2 + zz**2)
c take care of roundoff errors in argument of acos
        zzr   = zz/r
        if (zzr.gt.1)  zzr=1
        if (zzr.lt.-1) zzr=-1
c acos is identical to arcos
        theta = acos (zzr)
c force phi to 0 on z axis
        phi = 0
        if ((xx**2 + yy**2) .gt. 1.d-12) phi = atan2 (yy,xx)
        fi2xyz = fi4rtp (r,theta,phi,0,0,nterms)
        return
        end
 
 
 
        function fi4rtp (r,theta,phi,theta2, phi2, nterms)
c
c       returns the value of the nh3-h2) potential for
c       polar coordinates of h2 (green conventions).
c       calls the routine drot... for rotation matrix
c       calls the interpolation routine fi4vij.
c       *** warning *** theta and phi are expected to be in radians.
c       *** warning *** fi4rtp is much slower than fi4vij as it
c                       makes nterms calls to drot.
c       beware that the decomposition on rotation matrix used
c       here is inefficient to interpolate on random
c       (r,theta,phi,theta2,phi2).
c       moreover, the routine drot is not optimised in this case.
c       original version                       p. valiron 14-june-84
c       modified arguments                     p. valiron 19-sept-87
c       adaptation to NH3-H2                   c. rist    17-jan-92
        implicit double precision (a-h,o-z)
        integer nterms, id, jd, i2d, j2d, in, l, m, mul
        double precision pi, r, theta, phi, v, drot, fit, fi4rtp
        data pi /3.1415926535897932 d0/
        dimension v(45),id(45),jd(45),i2d(45), j2d(45)
c       table of correspondance between v index and related ylm
c
c               1     2     3     4     5     6     7     8
c               0000  1000  2000  3000  3300  4000  4300  5000
c               9     10   11   12
c               y5300 6000 6300 6600
c               13    14   15   16   17   18   19   20
c               0020  1020 2020 3020 3320 4020 4320 5020
c               21    22   23   24
c               5320  6020 6320 6620
c               25    26   27   28   29   30   31   32
c               1021  2021 3021 3321 4021 4321 5021 5321
c               33    34   35
c               6021  6321 6621
c               36    37   38   39   40   41   42   43
c               2022  3022 3322 4022 4322 5022 5322 6022
c               44    45
c               6322  6622
 
        data id /0,1,2,3,3,4,4,5,5,6,6,6,0,1,2,3,3,4,4,5,5,6,6,6,
     :           1,2,3,3,4,4,5,5,6,6,6,2,3,3,4,4,5,5,6,6,6/
        data jd /0,0,0,0,3,0,3,0,3,0,3,6,0,0,0,0,3,0,3,0,3,0,3,6,
     :           0,0,0,3,0,3,0,3,0,3,6,0,0,3,0,3,0,3,0,3,6/
        data i2d/0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,2,2,
     :           2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
        data j2d/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     :           1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2/
 
c
        do 10 in=1,nterms
            call vstar (in,r,v(in))
10      continue
        fit=0
        do 100 in=1,nterms
            l=id(in)
            m=jd(in)
            l2=i2d(in)
            m2=i2d(in)
*            mul=1
*            if (m.gt.0) mul=2
*c angle in degrees for ylm2 and in radians for cos
*            ylm = ylm2(theta*180/pi,l,m)
*     1                  * cos(m * phi) * mul
* expansion NH3-H2 on rotation matrix:
            coeff = 1
            if(m.eq.0)coeff=coeff/2
            if(m2.eq.0)coeff=coeff/2
            ylm=2*coeff*drott(l2,m2,0,theta2)*
     :          (cos(m*phi+m2*phi2)*drott(l1,m,m2,theta) +
     :        (-1)**m*cos(m*phi-m2*phi2)*drott(l1,-m,m2,theta))
            fit = fit + ylm * v(in)
100     continue
        fi4rtp = fit
        return
        end
 
 
 
        subroutine fi4vij (r,v,in,ider)
c
c       decomposition of the nh3-h2(j=0) potential on spherical
c       harmonics for a given intermolecular distance according
c       to s. green conventions
c
c       argument list
c
c       r       internuclear distance
c       v       vij coefficient returned by the routine
c       in      index of requested vij coefficient
c       ider    1 for vij, 2 and 3 for first or second derivatives
c
c       calls to fi4vij are optimised for a fixed value of isurf (set in
c       common /fisurf/). this parameter should not be frequently changed.
c               1       scf surface
c               2       disp surface
c               3       scf+disp
c
c       relation between v index and associated vij
c
c               1     2     3     4     5     6     7     8
c               0000  1000  2000  3000  3300  4000  4300  5000
c               9     10   11   12
c               y5300 6000 6300 6600
c               13    14   15   16   17   18   19   20
c               0020  1020 2020 3020 3320 4020 4320 5020
c               21    22   23   24
c               5320  6020 6320 6620
c               25    26   27   28   29   30   31   32
c               1021  2021 3021 3321 4021 4321 5021 5321
c               33    34   35
c               6021  6321 6621
c               36    37   38   39   40   41   42   43
c               2022  3022 3322 4022 4322 5022 5322 6022
c               44    45
c               6322  6622
c       the potential data is read on the logical unit ifile by  ddyini.
c       this file is read every time the value of isurf is changed.
c       informational messages are output on logical unit iwrite.
c
c       original version                       p. valiron 14-june-84
c       modified arguments                     p. valiron 19-sept-87
c       adaptation to new expansion NH3-H2     c. rist 19-jan-92
c       interpolation scheme
c       ====================
c               all scf distances are taken into account during the
c               interpolation process.
c               first, scf and disp energies are interpolated and
c               extrapolated for each configuration on a fine distance
c               grid.
c               then periodic bicubic splines are used to perform an
c               interpolation on (theta,phi) and the vij coefficients
c               are obtained by a least squares fit on the portion of
c               the sphere theta1=0...180 and phi1=0...120 .
c               the interpolation of the resulting vij coefficients
c               on the fine distance grid is made by usual cubic splines
c               this step provides an accurate evaluation of vij
c
c       parameter nddmx         max number of distances
c       parameter nvmx          max number of vij terms
        implicit double precision (a-h,o-z)
        integer nddmx, nvmx
        parameter (nddmx=85, nvmx=45)
        integer in, ider, isurf, lsurf, iwrite, ifile, ifile2, nv, ndd
     1          , k
        double precision r, v, dd, y, y1, y2, y3, yref, h, conv, econv
        dimension dd(nddmx),y(nddmx,nvmx),y1(nddmx,nvmx)
     1       ,y2(nddmx,nvmx),y3(nddmx,nvmx) , yref(nddmx)
        logical damped
        common /fisurf/ conv,econv,isurf,lsurf,damped
        common /fiunit/ iwrite,ifile,ifile2
        save k, nv, ndd, dd, y, y1, y2, y3
c
c compute surface
c extrapolated within 3.5 to 100 au.
c accurate evaluation of v00 to v66. coarse evaluation of v70 to v99.
c algorithm is fast but does not allow frequent changes in isurf value.
        if (isurf.eq.lsurf) goto 3000
c initialisation of the decomposition for all r according to isurf
        call ddyini (dd,y,y1,y2,y3,yref,nddmx,nvmx,ndd,nv)
        lsurf=isurf
c get the decomposition at r
3000    call splget (ndd,dd,r,k)
        h = r - dd(k)
        if (in .lt. 1 .or. in .gt. nv) goto 3010
        if (ider .lt. 1 .or. ider .gt. 3) goto 3020
        go to (111,222,333), ider
c y2 and y3 have been scaled for optimisation
 111    v = y(k,in)+h*(y1(k,in)+h*(y2(k,in)+h*y3(k,in)))
        return
 222    v = y1(k,in) + 2.d0*h*y2(k,in) + 3.d0*h*h*y3(k,in)
        return
 333    v = 2.d0*y2(k,in) + 6.d0*h*y3(k,in)
        return
c
3010    continue
        write (6,888)
        write (iwrite,888)
 888    format (/,'0fi2vij -- term in out of range',/)
        stop
3020    continue
        write (6,444)
        write (iwrite,444)
 444    format (/,'0fi2vij -- ider out of range',/)
        stop
        end
 
 
 
        subroutine ddyini (dd,y,y1,y2,y3,yref,nddmx,nvmx,ndd,nv)
c
c       reads the vij decomposition on logical unit ifile
c       sets up the coefficients of the spline for scf, disp or scf+disp
c       according to isurf
c
        implicit double precision (a-h,o-z)
        integer nddmx, nvmx, ndd, nv, isurf, lsurf
     +        , iwrite, ifile, ifile2, idd, iv, iflag, l
        logical llcut , damped
        double precision dd, y, y1, y2, y3, yref, rcut, drcut, conv,
     +                   econv, damp
        dimension dd(nddmx),y(nddmx,nvmx),y1(nddmx,nvmx)
     1         ,y2(nddmx,nvmx),y3(nddmx,nvmx) , yref(nddmx)
        common /cutoff/ rcut,drcut,llcut
        common /fisurf/ conv,econv,isurf,lsurf,damped
        common /fiunit/ iwrite,ifile,ifile2
        character*60 datafl
        common /finame/ datafl
c
        if (isurf.lt.1 .or. isurf.gt.3) then
            write (6,111)
            write (iwrite,111)
111         format (/,'0ddyini-- isurf out of range',/)
            stop
        endif
c
c read distances
c
        open (unit=ifile, file=datafl,  status='old')
        read(ifile,201) ndd,nv
201     format(/2i10)
        if (ndd.gt.nddmx) then
            write (6,222)
            write (iwrite,222)
222         format (/,'0ddyini-- increase dimension nddmx',/)
            stop
        endif
        if (nv.gt.nvmx) then
            write (6,333)
            write (iwrite,333)
333         format (/,'0ddyini-- increase dimension nvmx',/)
            stop
        endif
        read(ifile,202) (dd(idd),idd=1,ndd)
202     format(/(10f8.3))
c ***************************************************************
c convert input units(bohr) to conv units (angstroms for molscat)
c ***************************************************************
      do 100 idd=1,ndd
        dd(idd) = dd(idd) * conv
100   continue
c
c read scf interpolated coefficients
c
        read(ifile,203) iflag
203     format(i10)
        if (iflag.ne.1) then
            write (6,444)
            write (iwrite,444)
444         format (/,'0ddyini-- iflag ne 1',/)
            stop
        endif
        do 40 iv=1,nv
            read(ifile,204) (y1(idd,iv),idd=1,ndd)
204         format(/(4f20.12))
40      continue
c*********************************************************************
c  set llcut=.true. to enable the cutoff of all vij terms excepted v00
c  in the scf expansion.
c  rcut  =  cutoff distance
c  drcut =  cutoff width
c  hint **** the cutoff width should be much larger than the grid step
c            around rcut to prevent spurious spline oscillations.
c  the cutoff function is an hyperbolic tangent in order to reach
c  1 and 0 faster than a polynomial.
c*********************************************************************
      if (.not.llcut) goto 1000
      do 101 idd=1,ndd
        damp = (1.d0-tanh((dd(idd)-rcut)/drcut)) / 2.d0
c let v00 unchanged
        do 102 iv=2,nv
          y1(idd,iv)=y1(idd,iv) * damp
102     continue
101   continue
      write(6,1001) rcut,drcut
      write(iwrite,1001) rcut,drcut
1001  format(/,1x,30('* * '),/,' cutoff function for vij scf coeffic'
     1,'ients excepted v00 unchanged',
     1/' hyperbolic tangent -- rcut =',f12.6,'    drcut =',f12.6
     1,'   angstroms'/,1x,30('* * '))
1000  continue
c
c read disp interpolated coefficients
c
        read(ifile,203) iflag
        if (iflag.ne.2) then
            write (6,555)
            write (iwrite,555)
555         format (/,'0ddyini-- iflag ne 2',/)
            stop
        endif
        do 50 iv=1,nv
        read(ifile,204) (y2(idd,iv),idd=1,ndd)
50      continue
c
c damp the dispersion to values read with iflag=3
c
        if (.not.damped) goto 60
        read(ifile,203) iflag
        if (iflag.ne.3) then
            write (6,567)
            write (iwrite,567)
567         format (/,'0ddyini-- iflag ne 3',/)
            stop
        endif
        read(ifile,204) (yref(idd),idd=1,ndd)
        do 301 idd=1,ndd
          do 302 iv=2,nv
            y2(idd,iv)=y2(idd,iv) * (yref(idd)/y2(idd,1))
302       continue
          y2(idd,1)=yref(idd)
301     continue
        write(6,1003)
        write(iwrite,1003)
1003    format(/,1x,30('* * '),/,' special uniform damping applied to '
     1  ,'dispersion.' /,1x,30('* * '))
60      continue
c
c accumulate coefficients in array y according to isurf
c
        goto (61,62,63),isurf
61      do 261 iv=1,nv
        do 261 idd=1,ndd
            y(idd,iv)=y1(idd,iv)
261     continue
        write (6,1) dd(1),dd(ndd)
        write (iwrite,1) dd(1),dd(ndd)
1       format(/' fi2vij -- initialisation of scf surface for '
     1  ,'nh3-h2'/' interp=2 / distances ',f6.2,' to ',f6.2,' a  '
     1  ,'/ v00-v66 and estimation for v70-v99'/)
        goto 70
62      do 262 iv=1,nv
        do 262 idd=1,ndd
            y(idd,iv)=y2(idd,iv)
262     continue
        write (6,2) dd(1),dd(ndd)
        write (iwrite,2) dd(1),dd(ndd)
2       format(/' fi2vij -- initialisation of disp surface for '
     1  ,'nh3-h2'/' interp=2 / distances ',f6.2,' to ',f6.2,' a  '
     1  ,'/ v00-v66 and estimation for v70-v99'/)
        goto 70
63      do 263 iv=1,nv
        do 263 idd=1,ndd
            y(idd,iv)=y1(idd,iv)+y2(idd,iv)
263     continue
        write (6,3) dd(1),dd(ndd)
        write (iwrite,3) dd(1),dd(ndd)
3       format(/' fi2vij -- initialisation of scf+disp surface for '
     1  ,'nh3-h2'/' interp=2 / distances ',f6.2,' to ',f6.2,' a  '
     1  ,'/ v0000-v6622 '/)
c
c initialise spline coefficients
c
70      continue
        do 80 iv=1,nv
            call cubspl (ndd,dd,y(1,iv),y1(1,iv),y2(1,iv),y3(1,iv),0,0)
c scale y y1 y2 y3 to units of energy
c scale y2 and y3 to obtain directly taylor coefficients
            do 80 idd=1,ndd
                y (idd,iv) = y (idd,iv)        / econv
                y1(idd,iv) = y1(idd,iv)        / econv
                y2(idd,iv) = y2(idd,iv) / 2.d0 / econv
                y3(idd,iv) = y3(idd,iv) / 6.d0 / econv
80      continue
c
c end of initialisation
c
        close (unit=ifile)
        return
        end
 
 
 
      subroutine splget (n,x,t,k)
      integer n, k
      double precision x, t
      dimension x(n)
c
c***********************************************************************
c     *** the subroutine splget modifies the index k so that the
c     argument t lies within the interval @x(k)...x(k+1)|.
c     in case of extrapolation, k is forced to the value 1 or n-1.
c
c     n      number of data points (n is assumed .ge. 2).
c     (x(i), i=1,...,n) abcissae of the points
c            (x is assumed to be strictly increasing).
c     t      argument for which the spline function is to be determined.
c     k      initial guess for k.
c
c                                       p. valiron  8-june-84
c***********************************************************************
c
      if(k.lt.1) k=1
      if(k.gt.n-1) k=n-1
      if(t.le.x(k+1)) go to 11
   10 if(k.eq.n-1) goto 20
      k=k+1
      if(t.gt.x(k+1)) go to 10
      go to 20
   11 if(k.eq.1) goto 20
      if(t.ge.x(k)) go to 20
      k=k-1
      go to 11
   20 return
      end
 
 
 
      subroutine cubspl (n,tau,c1,c2,c3,c4,ibcbeg,ibcend)
c******  piecewise cubic spline interpolants computation; adapted from
c  'a practical guide to splines' , carl de boor , applied mathematical
c  sciences, springer-verlag, vol.27, p57-59 (1978).
c     ************************* input **************************
c     n = number of data points. assumed to be .ge. 2.
c     (tau(i), c1(i), i=1,...,n) = abscissae and ordinates of the
c        data points. tau is assumed to be strictly monotonous.
c     ibcbeg, ibcend = boundary condition indicators, and
c     c2(1) , c2(n)  = boundary condition information. specifically,
c        ibcbeg = 0  means no boundary condition at tau(1) is given.
c           in this case, the not-a-knot condition is used, i.e. the
c           jump in the third derivative across tau(2) is forced to
c           zero, thus the first and the second cubic polynomial pieces
c           are made to coincide.
c        ibcbeg = 1  means that the slope at tau(1) is made to equal
c           c2(1), supplied by input.
c        ibcbeg = 2  means that the second derivative at tau(1) is
c           made to equal c2(1), supplied by input.
c        ibcend = 0, 1, or 2 has analogous meaning concerning the
c           boundary condition at tau(n), with the additional infor-
c           mation taken from c2(n).
c     ********************** output ****************************
c     n, tau, c1, c2, ibcbeg, ibcend  are not altered by cubspl.
c     cj(i), j=1,...,4; i=1,...,l (= n-1) = the polynomial coefficients
c        of the cubic interpolating spline with interior knots (or
c        joints) tau(2), ..., tau(n-1). precisely, in the interval
c        (tau(i), tau(i+1)), the spline f is given by
c           f(x) = c1(i)+h*(c2(i)+h*(c3(i)+h*c4(i)/3.)/2.)
c        where h = x - tau(i).
c     in other words, for i=1,...,n, c2(i) and c3(i) are respectively
c        equal to the values of the first and second derivatives of
c        the interpolating spline, and c4(i) is equal to the third deri-
c        vative of the interpolating spline in the interval (tau(i),
c        tau(i+1)). c4(n) is meaningless and is set to 0. for clarity.
c     **********************************************************
      implicit double precision (a-h,o-z)
      integer n, ibcbeg, ibcend, i, l, m, jj, j
      double precision tau, c1, c2, c3, c4, taum1, g, dtau, divdf1,
     :                 divdf3
      dimension tau(n),c1(n),c2(n),c3(n),c4(n)
c***** a tridiagonal linear system for the unknown slopes s(i) of
c  f at tau(i), i=1,...,n, is generated and then solved by gauss elim-
c  ination, with s(i) ending up in c2(i), all i.
c     c3(.) and c4(.) are used initially for temporary storage.
c
check -- n.ge.2
      if (n.lt.2) write (6,111)
 111  format (/,'0cubspl -- less than two pivots',/)
      if (n.lt.2) stop
check -- tau strictly monotonous
      taum1=tau(2)
      if (tau(2)-tau(1)) 101,102,103
  101 if (n.eq.2) goto 200
      do 1 i=3,n
      if ((tau(i)-taum1).ge.0.d0) goto 102
    1 taum1=tau(i)
      goto 200
  102 write (6,222)
  222 format (/,'0cubspl -- non monotonous abscissae',/)
      stop
  103 if (n.eq.2) goto 200
      do 3 i=3,n
      if ((tau(i)-taum1).le.0.d0) goto 102
    3 taum1=tau(i)
c
  200 l = n-1
compute first differences of tau sequence and store in c3(.). also,
compute first divided difference of data and store in c4(.).
      do 10 m=2,n
         c3(m) = tau(m) - tau(m-1)
   10    c4(m) = (c1(m) - c1(m-1))/c3(m)
construct first equation from the boundary condition, of the form
c             c4(1)*s(1) + c3(1)*s(2) = c2(1)
      if (ibcbeg-1)                     11,15,16
   11 if (n .gt. 2)                     goto 12
c     no condition at left end and n = 2.
      c4(1) = 1.d0
      c3(1) = 1.d0
      c2(1) = 2.d0*c4(2)
                                        goto 25
c     not-a-knot condition at left end and n .gt. 2.
   12 c4(1) = c3(3)
      c3(1) = c3(2) + c3(3)
      c2(1) = ((c3(2)+2.d0*c3(1))*c4(2)*c3(3)+c3(2)**2*c4(3))/c3(1)
                                        goto 19
c     slope prescribed at left end.
   15 c4(1) = 1.d0
      c3(1) = 0.d0
                                        goto 18
c     second derivative prescribed at left end.
   16 c4(1) = 2.d0
      c3(1) = 1.d0
      c2(1) = 3.d0*c4(2) - c3(2)/2.d0*c2(1)
   18 if(n .eq. 2)                      goto 25
c  if there are interior knots, generate the corresp. equations and car-
c  ry out the forward pass of gauss elimination, after which the m-th
c  equation reads    c4(m)*s(m) + c3(m)*s(m+1) = c2(m).
   19 do 20 m=2,l
         g = -c3(m+1)/c4(m-1)
         c2(m) = g*c2(m-1) + 3.d0*(c3(m)*c4(m+1)+c3(m+1)*c4(m))
   20    c4(m) = g*c3(m-1) + 2.d0*(c3(m) + c3(m+1))
construct last equation from the second boundary condition, of the form
c           (-g*c4(n-1))*s(n-1) + c4(n)*s(n) = c2(n)
c     if slope is prescribed at right end, one can go directly to back-
c     substitution, since c array happens to be set up just right for it
c     at this point.
      if (ibcend-1)                     21,30,24
   21 if (n .eq. 3 .and. ibcbeg .eq. 0) goto 22
c     not-a-knot and n .ge. 3, and either n.gt.3 or also not-a-knot at
c     left endpoint.
      g = c3(l) + c3(n)
      c2(n) = ((c3(n)+2.d0*g)*c4(n)*c3(l)
     *            + c3(n)**2*(c1(l)-c1(n-2))/c3(l))/g
      g = -g/c4(l)
      c4(n) = c3(l)
                                        goto 29
c     either (n=3 and not-a-knot also at left) or (n=2 and not not-a-
c     knot at left endpoint).
   22 c2(n) = 2.d0*c4(n)
      c4(n) = 1.d0
                                        goto 28
c     second derivative prescribed at right endpoint.
   24 c2(n) = 3.d0*c4(n) + c3(n)/2.d0*c2(n)
      c4(n) = 2.d0
                                        goto 28
   25 if (ibcend-1)                     26,30,24
   26 if (ibcbeg .gt. 0)                goto 22
c     not-a-knot at right endpoint and at left endpoint and n = 2.
      c2(n) = c4(n)
                                        goto 30
   28 g = -1.d0/c4(l)
complete forward pass of gauss elimination.
   29 c4(n) = g*c3(l) + c4(n)
      c2(n) = (g*c2(l) + c2(n))/c4(n)
carry out back substitution
   30 do 40 jj=1,l
      j = l + 1 - jj
   40    c2(j) = (c2(j) - c3(j)*c2(j+1))/c4(j)
c****** generate cubic coefficients in each interval, i.e., the deriv.s
c  at its left endpoint, from value and slope at its endpoints.
      do 50 i=2,n
         dtau = c3(i)
         divdf1 = (c1(i) - c1(i-1))/dtau
         divdf3 = c2(i-1) + c2(i) - 2.d0*divdf1
         c3(i-1) = 2.d0*(divdf1 - c2(i-1) - divdf3)/dtau
   50    c4(i-1) = (divdf3/dtau)*(6.d0/dtau)
c****** compute in addition c3(n). set c4(n) to 0.
         c3(n) = c3(l) + c4(l)*dtau
         c4(n) = 0.d0
                                        return
      end
 
 
 
      function ylm2 (th, ll, m)
c
c     ylm2(th,l,m) is an update of ylm  from the package abme
c       (braithwaite, w.j. , comp. phys. commun. 5 (1973) 390)
c
c       th      angle (degrees)
c       ll,m    usual quantum numbers (not increased by one)
c
c       it is no more necessary to call the routine within a special
c       loop on l and m. however computer time is saved when ylm2
c       is called consecutively with same values of th and m and
c       increasing values of l.
c                                       p. valiron 11-apr-83
c
c updated to real*8 precision           p. valiron 21-may-84
c modified to fortran 4 syntax          p. valiron 14-jun-84
c
        double precision ylm2, th, c0, th0, root, factr, fp, arg, cx, sx
     +         , yx, ylmx, cf, y0
        integer ll, m, im, iroo, i, i2, nextl, m0, l, l2
        logical ifresh
        parameter (im=51,iroo=2*im+1)
        dimension root(iroo),factr(im)
c       dimension root(103),factr(51)
c       data iroo,im /103,51/
        data ifresh/.true./, c0/.2820947917738781d0/, th0/1.d37/
c
c       if(ifresh) then
        if (.not.ifresh) goto 501
                ifresh=.false.
                i2=2*im+1
                do 300 i=1,i2
300             root(i)=sqrt(float(i))
                fp=1.d0
                do 310 i=1,im
                i2=2*i
                fp=-fp*root(i2-1)/root(i2)/10.d0
310             factr(i)=c0*root(i2+1)*fp
c       endif
501     continue
c
c       if(th.ne.th0)then
        if(th.eq.th0) goto 502
                arg=th*0.1745329251994330d-01
                cx=cos(arg)
                sx=10.d0*sin(arg)
                th0=th
                nextl=m
                m0=m
c       endif
502     continue
c
c       if(m.ne.m0)then
        if(m.eq.m0) goto 503
                nextl=m
                m0=m
c       endif
503     continue
c
        if(ll.lt.m.or.m.lt.0) goto 100
        if(ll.lt.nextl) nextl=m
        do 2000 l=nextl,ll
        if(l.gt.0) goto 500
      yx=c0
      ylmx=0.
      go to 1000
  500 if(l.gt.m.or.m.eq.0) go to 900
      if(l.ge.im) write (6,111)
 111  format (/,'0ylm2 -- augmenter dimensions internes',/)
      if(l.ge.im) stop
      yx=factr(l)*sx**l
      ylmx=0.
      go to 1000
  900 l2=2*l
      cf=root(l2+1)/(root(l+m)*root(l-m))
      if(l.gt.m+1) go to 910
      yx=cf* root(l2-1)*cx*ylmx
      go to 1000
  910 yx=cf*(root(l2-1)*cx*ylmx-root(l+m-1)*root(l-m-1)*y0/root(l2-3))
 1000 y0=ylmx
      ylmx=yx
 2000 continue
      ylm2=ylmx
      nextl=ll+1
      return
  100 ylm2=0.
      return
      end
 
      function drott(l,m1,m2,theta)
c     routine bidon provisoire
      implicit double precision (a-h,o-z)
      drott = 1.d0
      return
      end
