*system:  NH3-H2 - Pierre Valiron/Claire Rist
*  details are published in Maret at el., MNRAS 399, 425 (2009)
*  subr provided by Claire Rist - Nov-2011
*
*  revised by p.j.dagdigian
*  current revision date:  18-nov-2011
*
*  Note:  this pot routine requires a data file to be in hibxx/bin/progs/potdata:
*         fitvij_bf_62.h2
*
      include "common/syusr"
      include "common/ground"
      include "common/bausr"
* ------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      character*60 filnam
      common /covvl/ vvl(300)
      common /fisurf/ conv,econv,lsurf,nv
      common /coloapot/ s4pi, nvv, ivv(300)
      common /conlam/ nlam, nlammx, lamnum
      include "common/parpot"
      potnam='RIST/VALIRON NH3-H2 2009'
      nlammx = 55
      nlam = 55
      call loapot(iunit, filnam)
      econv=219474.6d0
      conv=1.d0
      s4pi = sqrt ( 4.d0 * acos(-1.d0) )
      print *, potnam
1     print *, ' R (bohr)'
      read (5, *, end=99) r
      call pot(vv0,r)
*  potential is returned in atomic units (hartree)
*  convert from atomic units for printout
      write (6, 100) vv0*econv*s4pi, (vvl(iv)*econv,iv=1,nv)
100   format(' v',/,10(10(1pe16.8)/))
      goto 1
99    end
* ------------------------------------------------------------------------
               subroutine loapot(iunit, filnam)
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
*... jmax12:
*    ipotsy2:  symmetry of potential. if linear molecule is homonuclear
*              then ipotsy=2 and only terms with lambda2  even can be
*              included in the potential,else ipotsy=1.
*    j2max:    the maximum rotational angular momentum for linear
*              molecule
*    j2min:    the minimum rotational angular momentum for linear
*              molecule

*  variables in common block /cobspt/
*              Order of reduced rotation matrix d(theta1) as defined in 
*              eq (21) of reference J. Chem Phys. 98 (6), 1993 with    
*              (lambda, mproj) = (l1, m1).   
*    lammin:   array containing minimum value of lambda for each term
*    lammax:   array containing maximum value of lambda for each term
*    mproj:    array containing the order of the reduced rotation matrix
*              elements for each term. here, lammin is greater to mproj.
*  variables in common block /cobsptln/
*              Order of reduced rotation matrix d(theta2) as defined in 
*              eq (21) of reference J. Chem Phys. 98 (6), 1993 with
*              (lam2, m2proj) = (l2, m2).      
*    lam2:     array containing the order of the reduced rotation matrix
*              elements for each term. in case of homonuclear molecule
*              is even.
*    m2proj:   array containing the order of the reduced rotation
*              matrix elements for each term.
*              here, lammin and lam2 are greater than m2proj .
 
*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical readpt
      include "common/parbas"

      integer         nscode, isicod, nterm, numpot, ipotsy, iop, ninv,
     :                kmax, jmax0, jmax1, jmax2, jmax3, jmax4, jmax5,
     :                jmax6, jmax7, jmax8, jmax9, jmax10, jmax11, 
     :                jmax12, ipotsy2, j2max, j2min

      common /cosysi/ nscode, isicod, nterm, numpot, ipotsy, iop, ninv,
     :                kmax, jmax0, jmax1, jmax2, jmax3, jmax4, jmax5,
     :                jmax6, jmax7, jmax8, jmax9, jmax10, jmax11, 
     :                jmax12, ipotsy2, j2max, j2min
 
*  variable in common block /conlam/
*    nlammx:    the total number of angular coupling terms
*    nlam:      the maximum number of angular coupling terms allowed
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx (not used here)

*      integer         nlam, nlammx, lamnum(1)
 
      common /conlam/ nlam, nlammx, lamnum

c  variables in common block /fisurf/
c  ----------------------------------
 
c    conv and econv : set the unit conversion for distances and energies.
c                these conversions are designed to match the units in use
c                in the collision code. 
c                For NH3-H2, input data are read on unit ifile and
c                are expressed in au for distances and in cm-1
c                for energies. output data and listed values are converted
c                into au and Hartrees, according to conv and econv.
c    conv:       atomic unit of distance expressed in user units.
c                set conv = 1.0 for atomic units
c                set conv = 0.529177 for angstroms
c    econv:      user energy unit expressed in cm-1.
c                set econv = 1.0 for cm-1.
c                set econv = 2 * 109737.312 for hartrees.
c    lsurf:      Usually set to 1. Setting lsurf to 0 forces reading
c                again the data on ifile.
 
      common /fisurf/ conv, econv, lsurf
 
c  variables in common blocks /fiunit/ and /finame/
c  ------------------------------------------------
 
c    iwrite       logical unit for listing output (on file).
c                 logical unit for display is assumed equal to 6.
c    ifile        logical unit to read the potential expansion.
c                 distances are read in bohrs and energies in cm-1.
c    datafl:      name of data file (pot. expansion) read on unit ifile
c    datanm:      name of namelist data file ( set with pot= )
 
      common /fiunit/ iwrite,ifile
      character*(*) filnam
      include "common/parpot"
      parameter (nvmx = 300)
      integer   ivij(nvmx), jvij(nvmx), i2vij(nvmx), j2vij(nvmx),
     :          lambda(nvmx), mu(nvmx),
     :          lambda2(nvmx), mu2(nvmx)
   
*      common /coloapot/ s4pi, ivv(nvmx), nvv
      common /coloapot/ s4pi, nvv, ivv(300)
      logical is_wop
      
      potnam='RIST/VALIRON NH3-H2 2009'

      nterm = 16
      l2m2 = 6
      m1max = 6
      l1max = 6
      l2max = 2
      
      if (l2max .gt. 4) then
         write(6,*) 'l2max too high) -- set l2max < 6. stop'
         stop    
      if (l2max .eq. 0) l2m2 = 1
      if (l2max .eq. 2) l2m2 = 4
      if (l2max .eq. 4) l2m2 = 15  
      endif
      imax = m1max/3 +1
      nterm = imax * l2m2 - 2
      n = 0
      k = 0
      do j=1, l2m2
        do i=1, imax
          if (.not. (i.eq.1 .and. (j.eq.3 .or. j.eq.5))) then
          k = k + 1
          mproj(k) = 3*(i-1)
          if(j.eq.1) then 
              m2proj(k) = 0
              lam2(k) = 0
	  elseif (j.eq.2) then
              m2proj(k) = 0
              lam2(k) = 2
	  elseif (j.eq.3) then
              m2proj(k) = -1
              lam2(k) = 2
	  elseif (j.eq.4) then 
              m2proj(k) = 1
              lam2(k) = 2
	  elseif (j.eq.5) then
              m2proj(k) = -2
              lam2(k) = 2
	  elseif (j.eq.6) then
              m2proj(k) = 2 
              lam2(k) = 2
	  elseif (j.eq.7) then
              m2proj(k) = 0
              lam2(k) = 4
	  elseif (j.eq.8) then 
              m2proj(k) = -1 
              lam2(k) = 4
	  elseif (j.eq.9) then 
              m2proj(k) = 1 
              lam2(k) = 4
	  elseif (j.eq.10) then 
              m2proj(k) = -2
              lam2(k) = 4
	  elseif (j.eq.11) then
              m2proj(k) = 2
              lam2(k) = 4
	  elseif (j.eq.12) then 
              m2proj(k) = -3
              lam2(k) = 4
	  elseif (j.eq.13) then 
              m2proj(k) = 3
              lam2(k) = 4
	  elseif (j.eq.12) then 
              m2proj(k) = -4
              lam2(k) = 4
	  elseif (j.eq.13) then 
              m2proj(k) = 4
              lam2(k) = 4
          endif
          endif
          lammax(k) = l1max
          lammin(k) = max(mproj(k), iabs(m2proj(k)))
        enddo
      enddo
      lammin(1) = 1
200   continue

      iwrite = 9
      ifile  = 60
 
c     conv = 0.529177d0
      conv = 1.d0
c     econv = 1.0d0
      econv = 109737.312d0 * 2

*  calculate total number of anisotropic terms
      nlam = 0
      do 135  i = 1, nterm
        if (mproj(i) .gt. lammin(i) ) then
          write (6, 130)  mproj(i), lammin(i), i
          write (9, 130)  mproj(i), lammin(i), i
130       format (' *** mproj=',i2,' > lammin=',i2,
     :            ' for term',i2,'. stop ***')
          stop
        end if
        if (m2proj(i) .gt. lammin(i) ) then
          write (6, 131)  m2proj(i), lammin(i),i
          write (9, 131)  m2proj(i), lammin(i), i
131       format (' *** m2proj=',i2,' > lammin=',i2,
     :            ' for term',i2,'; abort ***')
          stop
        end if
        if (m2proj(i) .gt. lam2(i) ) then
          write (6, 132)  m2proj(i), lam2(i),i
          write (9, 132)  m2proj(i), lam2(i),i
132       format (' *** m2proj=',i2,' > lam2=',i2,
     :            ' for term',i2,'. stop ***')
          stop
        endif
        if(lammax(i) .lt.0)goto 135 
        nlam = nlam + lammax(i) - lammin(i) + 1
135   continue
      if (nlammx .lt. nlam) then
        write (6, 140) nlam, nlammx
        write (9, 140) nlam, nlammx
140     format (' ** total number of anisotropic terms=', i2,
     :          ' .gt. nlam=', i2,'. stop')
        call abort
      end if

c  display and list parameters
      write(6,1000)
*      write(6,1001) iwrite,ifile
      write(6,1002) conv,econv

      inquire (iwrite, opened=is_wop)
      if (is_wop) then
         write(iwrite,1000)
*     write(iwrite,1001) iwrite,ifile
         write(iwrite,1002) conv,econv
      end if
 
1000  format(' loapot -- initialize potential surface parameters'
     +       /' potential data:  ',
     +        'potdata/fitvij_bf_62.h2')
1001  format('     iwrite     ifile'/,1x,2i10)
1002  format('     conv =',f9.6,'   econv = ',f10.3)
 
*  force reading the potential expansion ( with lsurf = 0 )
      lsurf = 0
      call vinit (1, rinit, v)

* match the order of potential terms in file potfil and in basis
* define in ivij, jvij, ivij2, jvij2 the order of potfil
      kv = 0
      do 210 i2v = 0, l2max, 2
      do 210 j2v = 0, i2v
      do 210 iv = j2v, l1max
      jvmax = min(iv, m1max)
      do 210 jv = 0, jvmax, 3
        kv = kv+1
        if (kv .gt. nvmx) then
          write(6,*) 'laopot -- dimension nvmx too small'
          stop
        endif
        ivij(kv) = iv
        jvij(kv) = jv
        i2vij(kv) = i2v
        j2vij(kv) = -j2v
        if (jv.ne.0 .and. j2v.ne.0) then
          kv = kv + 1
          ivij(kv) = iv
          jvij(kv) = jv
          i2vij(kv) = i2v
          j2vij(kv) = j2v
        else
          j2vij(kv) = j2v
        end if
        if (kv .gt. nvmx) then
          write(6,*) 'laopot -- dimension nvmx too small'
          stop
        endif
210   continue
      nvij = kv

* define in ivv a pointer from basis definition to the potfil
* definition.
      nvv = 0
      do 220 iterm = 1, nterm
        imu = mproj(iterm)
        ilam2 = lam2(iterm)
        imu2 = m2proj(iterm)
        if (lammax(iterm).lt.0) goto 220
        do 240 ilamb = lammin(iterm), lammax(iterm)
          nvv = nvv+1
          lambda(nvv) = ilamb
          mu(nvv) = imu
          lambda2(nvv) = ilam2
          mu2(nvv) = imu2
          if (nvv .gt. nvmx) then
            write(6,*) 'sysdat -- ivv initialization problem  '
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
 
          write(6,*) 'sysdat -- ivv initialization problem '
          write(6,*) ilamb, imu, ilam2, imu2,
     +               '   ilamb, imu, ilam2, imu2   '
          write(6,*)
     +    ' Unfound index in ivij, jvij, i2vij, j2vij'
          stop
240     continue
220   continue
      write(6,*)
      write(6,*)'number of potential terms found in subroutine loapot :'
     +           , nvv
      write(6,250) (kvv, lambda(kvv), mu(kvv), lambda2(kvv),
     +              mu2(kvv), ivv(kvv),
     +              ivij(ivv(kvv)), jvij(ivv(kvv)),
     +              i2vij(ivv(kvv)), j2vij(ivv(kvv)), kvv=1, nvv)
250   format(1x, i6, '   vvl',4i2, '   as  ', i3, '  ( v',4i2,' )' )
 
      s4pir8 = sqrt ( 4.d0 * acos(-1.d0) )
      s4pi = s4pir8
*
*  potential has been defined : force irpot=1
      irpot = 1
*      write(6,260) datafl
*260   format('Potential parameters successfully read from file:  ',
*     :        a)
      write(6,260) 
260   format('Potential parameters successfully read from file:  ',
     :        'potdata/fitvij_bf_62.h2')
      return
      end
*  -----------------------------------------------------------------------
      subroutine pot(vv0, r)
*  -----------------------------------------------------------------------
*  -----------------------------------------------------------------------
*                .........   nh3 - h2    .............
*  -----------------------------------------------------------------------
*  subroutine to calculate the r-dependent coefficients in the
*  collision of a symmetric top with a linear target
*  in units of hartree for energy and bohr for distance
 
*  on return:
*  vv0 contains the isotropic term d1[0,0,0]*d2[0,0,0] in the potential
*  the coefficients for each angular term in the coupling potential
*  [ vvl(i) for i = 1, nlam ] are returned in common block /covvl/
*  in the following order:[l1,m1,l2,m2] associated with the angular
*  function d(l1,m1,m2)*d(l2,m2,0)
*    vvl(1) contains the [1,0,0,0] term
*    vvl(2) contains the [2,0,0,0] term
*    ....
*    vvl(lammax(1)) contains the [lammax(1),0,0,0] term
*    vvl(lammax(1)+1) contains the [lammin(2),3,0,0]
*    ....
*    vvl(lammax(1)+lammax(2) - 2) contains the [lammax(2),3,0,0]
*    etc... [l2,m2] sorted as: [0,0],[2,0],[2,-+1],[2,-+2]
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
*  [6,6,2,-1],[6,6,2,1],[2,0,2,2]...[6,6,2,-2],[6,6,2,2].
*  ivv(kv) correspond to the order in which the potential terms are read
*  nvv is the total number of angular coupling terms
 
*      implicit none
      implicit double precision (a-h,o-z)
*      integer ivv, nvv, nlam, nlammx, kv, kvv, nvmx, lamnum
      include "common/parpot"

*      double precision  r, v, vv, vv0, vvl, s4pi
      double precision  r, v, vv, vv0, s4pi
      common /conlam/ nlam, nlammx, lamnum(1)
*      common /covvl/ vvl(nvmx)
*      common /coloapot/ s4pi, ivv(nvmx), nvv
      common /covvl/ vvl(300)
      common /coloapot/ s4pi, nvv, ivv(300)
 
      if (nlam .ne. nvv) then
        write(6,*) ' entry pot --', nvv, nlam,'   nvv, expected nlam'
c        stop
      endif

      call vstar(1,r,v)

*****?
c  vv0 should be scaled as it must be the actual isotropic potential.
c  vv0 is to be multiplied by the isotropic matrix element
c  MAKE SURE that the scaling is consistent with the normalisation of 
c  the basis used.
c
      vv0 = v
c      vv0 = v/s4pi

      do 400 kvv = 1, nvv
        kv = ivv(kvv)
        call vstar(kv,r,v)
        vvl(kvv) = v
400   continue
 
c      write(1,'(x,f12.6,1p5e20.12)') r, vv0, (vvl(kvv),kvv=1,min(nvv,4))
 
      return
      end
*  -----------------------------------------------------------------------
      subroutine vinit(i,r,v)

      double precision  r, v, conv, econv
      double precision  rstart
      double precision v0, v1, v2
      common /fisurf/ conv,econv,lsurf
      common /fiunit/ iwrite,ifile

      if (i.ne.1) return
      rstart = 10.d0
      call fit4 (rstart,v0,1,1)
      return 
c=========================
      entry vstar(i,r,v)
c=========================
 
      call fit4 (r,v,i,1)

      return
      end
*  -----------------------------------------------------------------------
        subroutine fit4(r,v,in,ider)
c
c       decomposition of the nh3-h2 potential on rotation matrices
c       for a given intermolecular distance according
c       to equation (21) of reference J. Chem Phys. 98 (6), 1993.
c
c       argument list
c
c       r       internuclear distance
c       v       vij coefficient returned by the routine
c       in      index of requested vij coefficient
c       ider    1 for vij, 2 and 3 for first or second derivatives
c
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
c               25    26   27   28    29   30   31    32   33   34    35
c               1021  2021 3021 332-1 3321 4021 432-1 4321 5021 532-1 5321
c               36    37    38   39    40
c               6021  632-1 6321 662-1 6621
c               41    42   43    44   45   46    47   48   49    50   51
c               2022  3022 332-2 3322 4022 432-2 4322 5022 532-2 5322 6022
c               52    53    54    55    
c               632-2 6322  662-2 6622
c       the potential data is read on the logical unit ifile by  ddyini.
c       this file is read every time the value of lsurf is set to 0.
c       informational messages are output on logical unit iwrite.
c
c       interpolation scheme
c       ====================
c               Here the interpolation of the resulting vij coefficients
c               on the fine distance grid is made by usual cubic splines
c               this step provides an accurate evaluation of vij
c
c       parameter nddmx         max number of distances
c       parameter nvmx          max number of vij terms

*        implicit none
        implicit double precision (a-h,o-z)
*        integer nddmx, nvmx
        include "common/parpot"
*        parameter (nvmx = 45)
        parameter (nddmx=250)
        double precision r, v, dd, y, y1, y2, y3, yref, h, conv, econv
*        dimension dd(nddmx),y(nddmx,nvmx),y1(nddmx,nvmx)
*     1       ,y2(nddmx,nvmx),y3(nddmx,nvmx) , yref(nddmx)
        dimension dd(250),y(250,300),y1(250,300)
     1       ,y2(250,300),y3(250,300) , yref(250)
        common /fisurf/ conv,econv,lsurf
        common /fiunit/ iwrite,ifile
        parameter (nvmx = 300)
        save k, nv, ndd, dd, y, y1, y2, y3
c
c compute surface
c algorithm is fast but does not allow frequent changes in lsurf value.
        if (lsurf.eq.0) then 
c initialisation of the decomposition for all r when lsurf = 0
          call ddyini (dd,y,y1,y2,y3,yref,nddmx,nvmx,ndd,nv)
          lsurf=1
        endif

c get the decomposition at r
        call splget (ndd,dd,r,k)
        h = r - dd(k)

        if (in .lt. 1 .or. in .gt. nv) then
          write (6,888)
          write (iwrite,888)
 888      format (/,' fit4 -- term in out of range',/)
          stop
        endif

        if (ider .lt. 1 .or. ider .gt. 3) then
          write (6,444)
          write (iwrite,444)
 444      format (/,' fi2vij -- ider out of range',/)
          stop
        endif

c y2 and y3 have been scaled for optimisation

        if (ider .eq. 1) then
          v = y(k,in)+h*(y1(k,in)+h*(y2(k,in)+h*y3(k,in)))
        elseif (ider .eq. 2) then
          v = y1(k,in) + 2.d0*h*y2(k,in) + 3.d0*h*h*y3(k,in)
        elseif (ider .eq. 3) then
          v = 2.d0*y2(k,in) + 6.d0*h*y3(k,in)
        endif

        return
        end
*  -----------------------------------------------------------------------
        subroutine ddyini (dd,y,y1,y2,y3,yref,nddmx,nvmx,ndd,nv)
c
c       reads the vijkl decomposition on logical unit ifile
c       sets up the coefficients of the spline interpolation
c
       implicit double precision (a-h,o-z)
        double precision dd, y, y1, y2, y3, yref, conv,
     +                   econv, damp
*        dimension dd(nddmx),y(nddmx,nvmx),y1(nddmx,nvmx)
*     1         ,y2(nddmx,nvmx),y3(nddmx,nvmx) , yref(nddmx)
        dimension dd(250),y(250,300),y1(250,300)
     1       ,y2(250,300),y3(250,300) , yref(250)
        common /fisurf/ conv,econv,lsurf,nv1
        common /fiunit/ iwrite,ifile
*        character*(*) datafl
*        common /finame/ datafl
c
c
c read distances
c
        open (unit=ifile, file=
     :    'potdata/fitvij_bf_62.h2',
     :    status='old')
        read(ifile,*) ndd
*  comment out write statement
*        write(6,221) hdd
*221     format (' ddyini - reading distances.  ndd: ', i4)
201     format(/i10)
        if (ndd.gt.nddmx) then
            write (6,222)
            write (iwrite,222)
222         format (/,' ddyini-- increase dimension nddmx',/)
            stop
        endif
        read(ifile,*) (dd(idd),idd=1,ndd)
*  comment out write statements
*        write(6,*) 'ddyini - reading distances'
*        write(6, 202) (dd(idd),idd=1,ndd)
202     format(/(8f10.4))
c ***************************************************************
c convert input units(bohr) to conv units (angstroms for molscat)
c ***************************************************************
      do idd=1,ndd
        dd(idd) = dd(idd) * conv
      enddo

        read (ifile,*) nv
        nv1=nv
*  comment out write statement
*        write(6,332) nv
*332     format (' ddyini - reading potential.  nv: ',i4)
        if (nv.gt.nvmx) then
            write (6,333)
            write (iwrite,333)
333         format (/,' ddyini-- increase dimension nvmx',/)
            stop
        endif

        do iv=1,nv
            read(ifile,*) l1, m1, l2, m2
*  comment out write statement
*            write(6, '(a,i3,a,4i2)')'reading term', iv,' = ',
*     :                l1,m1,l2,m2
            read(ifile,*) (y1(idd,iv),idd=1,ndd)
204         format(/(4f20.12))
40      enddo

        do iv=1,nv
          do idd=1,ndd
            y(idd,iv)=y1(idd,iv)
          enddo
        enddo

        write (6,1) dd(1),dd(ndd)
1       format(' fi2vij -- initialization of surface for '
     1  ,'nh3-h2  distances ',f6.2,' to ',f6.2,' bohr '
     1  ,'/ v0000-v6622 ')

c initialise spline coefficients
c
        do iv=1,nv
            call cubspl (ndd,dd,y(1,iv),y1(1,iv),y2(1,iv),y3(1,iv),0,0)
c scale y y1 y2 y3 to units of energy
c scale y2 and y3 to obtain directly taylor coefficients
            do idd=1,ndd
                y (idd,iv) = y (idd,iv)        / econv
                y1(idd,iv) = y1(idd,iv)        / econv
                y2(idd,iv) = y2(idd,iv) / 2.d0 / econv
                y3(idd,iv) = y3(idd,iv) / 6.d0 / econv
            enddo
        enddo
c
c end of initialisation
c
        close (unit=ifile)
        return
        end
*  -----------------------------------------------------------------------
      subroutine splget (n,x,t,k)
      integer n, k
      double precision x, t
      dimension x(n)
c
c***********************************************************************
c     *** the subroutine splget modifies the index k so that the
c     argument t lies within the interval |x(k)...x(k+1)|.
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
*  -----------------------------------------------------------------------
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
*  -----------------------------------------------------------------------
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
*  -----------------------------------------------------------------------
 
