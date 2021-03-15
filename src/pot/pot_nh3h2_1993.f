*  system: NH3-H2 using canonical Valiron expansion of potential
*  this pot routine provided by claire wrist - aug-2011
*
*  reference:  c. rist, m.h. alexander, and p. valiron, j. chem. phys.
*  998, 4662 (1993)
*
*  V(l1,m1,l2,m2) coeffs over a grid of R values provided by
*  claire rist  -- aug-2011
*  pot routine revised to use these coeffs by p. dagdigian
*
*  Note:  this pot routine requires a data file to be in hibxx/bin/progs/potdata:
*         nh3h2_vijfit_vtab.dat
*
      include "common/syusr"
      include "common/ground"
      include "common/bausr"
* ------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      character*60 filnam
      common /covvl/ vvl(44)
      include "common/parpot"
      potnam='RIST/VALIRON NH3-H2'
      print *, potnam
1     print *, ' R (bohr)'
      read (5, *, end=99) r
      call pot(vv0,r)
*  potential is returned in atomic units (hartree)
*  convert from atomic units for printout
      econv=219474.6d0
      write (6, 100) vv0*econv, vvl*econv
100   format(' v',/,45(1pe16.8))
      goto 1
99    end
* ------------------------------------------------------------------------
      subroutine loapot(iunit, filnam)
      implicit double precision (a-h,o-z)
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
*
*  variables in common block /conlam/
*    nlammx:      the total number of angular coupling terms
*    nlam:    the maximum number of angular coupling terms allowed
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx (not used here)
      common /conlam/ nlam, nlammx, lamnum(1)
      common /coloapot/ s4pi, ivv(45), nvv
      common /co2mol/ twomol
*
      parameter (nvmx = 45)
      dimension ivij(nvmx), jvij(nvmx), i2vij(nvmx), j2vij(nvmx),
     :          lambda(nvmx), mu(nvmx),
     :          lambda2(nvmx), mu2(nvmx)
      potnam='RIST/VALIRON NH3-H2 1996'
*
      nterm = 12
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
* 
*  calculate total number of anisotropic terms
      nlam = 0
      do 135  i = 1, nterm
        if (mproj(i) .gt. lammin(i) ) then
          write (6, 130)  mproj(i), lammin(i), i
          write (9, 130)  mproj(i), lammin(i), i
130       format (' *** mproj=',i2,' > lammin=',i2,
     :            ' for term',i2,'; abort ***')
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
     :            ' for term',i2,'; abort ***')
          stop
        endif
        nlam = nlam + lammax(i) - lammin(i) + 1
135   continue
      if (nlammx .lt. nlam) then
        write (6, 140) nlam, nlammx
        write (9, 140) nlam, nlammx
140     format (' ** total number of anisotropic terms=', i2,
     :          ' .gt. nlam=', i2,'; abort')
        stop
      end if
*
* match the order of potential terms in basis subr
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
          stop
240       continue
220   continue
      write(6,221) nvv
221   format (/'TERMS IN THE ANGULAR EXPANSION OF THE POTENTIAL'/
     :  '   NUMBER OF ANISOTROPIC TERMS ',i4)
      write(6,250) (kvv, lambda(kvv), mu(kvv), lambda2(kvv),
     +              mu2(kvv), ivv(kvv),
     +              ivij(ivv(kvv)), jvij(ivv(kvv)),
     +              i2vij(ivv(kvv)), j2vij(ivv(kvv)), kvv=1, nvv)
250   format(1x, i6, '   vvl',4i2, '   as  ', i3, '  ( v',4i2,' )' )
*  
      return
      end
*  -----------------------------------------------------------------------
      subroutine pot( vv0, r)
*  subroutine to calculate the r-dependent coefficients
*  in atomic units (distance and energy)
*  -----------------------------------------------------------------------
*  on entry:
*    r:          interparticle distance
*  on return:
*    vv0         contains isotropic term -- d1[0,0,0]*d2[0,0,0]
*  variable in common block /covvl/
*    vvl:        vector of length 45 to store r-dependence of each term
*                in potential expansion
*    the coefficients for each angular term in the coupling potential
*    [ vvl(i) for i = 1, nlam ] are returned in common block /covvl/
*    in the following order:[l1,m1,l2,m2] associated with the angular
*    function drot(l1,m1,m2)*drot(l2,m2,0)
*      vvl(1) contains the [1,0,0,0] term
*      vvl(2) contains the [2,0,0,0] term
*      ....
*      vvl(lammax(1)) contains the [lammax(1),0,0,0] term
*      vvl(lammax(1)+1) contains the [lammin(2),3,0,0]
*      ....
*      vvl(lammax(1)+lammax(2) - 2) contains the [lammax(2),3,0,0]
*      etc... [l2,m2] sorted as: [0,0],[2,0],[2,1],[2,2]
*  variable in common block /coloapot/
*    s4pi : facteur de normalisation du potentiel isotrope
*    ivv(1) corresponds to the [1,0,0,0] term
*    ivv(2) corresponds to the [2,0,0,0] term
*    ivv(3) corresponds to the [3,0,0,0] term
*    ivv(4) corresponds to the [3,3,0,0] term
*    [4,0,0,0],[4,3,0,0]...[6,6,0,0],[0,0,2,0]...[6,6,2,0],[1,0,2,1],...
*    [6,6,2,1],[2,0,2,2]...[6,6,2,2].
*    ivv(kv) correspond to the order in which the potential terms are read
*    nvv is the total number of angular coupling terms
*
*  author:  paul dagdigian
*  major revision of routine written by p. valiron (mid 1990's)
*  current revision date:  11-aug-2011
*
      implicit double precision (a-h,o-z)
      dimension iwork(1000)
      dimension swork(12), work(1812)
      parameter (nvmx=45)
      common /conlam/ nlam, nlammx, lamnum(1)
      common /covvl/ vvl(44)
      common /coloapot/ s4pi, ivv(44), nvv
*
      common /vspline/ vsp(45)
      common /vvalue/ v(209,45), vvec(209)
      econv = 219474.6d0
*
*  call subroutine spline_nh3h2 to get values of vlam coefficients at distance r
      call spline_nh3h2(vsp,r)
*
*  copy to vvl array and convert to atomic units
      call dcopy(45,vsp(2),1,vvl,1)
      tohat = 1.d0/econv
      vv0 = vsp(1)*tohat
      call dscal(45,tohat,vvl,1)
*
      return
      end
*  -----------------------------------------------------------------------
      subroutine spline_nh3h2(vsp, r)
c  
c  based on spline_ch2he_3.f, using Jacek's code to do spline_fit of pot
c  example:
c    K=57
c    call splinej(K,Rp,XX,b0,c0,d0)
c    VNO=sevalj(K, r, Rp, XX, b0,c0,d0)
c  where,
c  K  the number of points to be fitted
c  Rp, XX   original data, Rp & XX should have the same dimensions
c  r  new distance, where pot need to be calculated by spline-Fitchbur
c  VNO  pot at r, output
      implicit double precision (a-h,o-z)   
      dimension b0(209,45),c0(209,45),a0(209,45),rr(209)
      dimension vsp(45), v(209,45),vvec(209,45)
      data ifirst /0/
      if (ifirst.eq.0) then
         open (unit=10,file=
     :     'potdata/nh3h2_vijfit_vtab.dat')
         do i=1,209
           read (10,*) rr(i),(v(i,j),j=1,45)
         enddo
         close(10)
*
*  determine spline coefficients for all 45 vlam coeffs
         do j=1,45
           call dcopy(209,v(1,j),1,vvec(1,j),1)
*  vvec now contains potentials at all 209 values of R (rows) for each 
*  45 vlam valuesi (columns)
           call spline(209,rr,vvec(1,j),a0(1,j),b0(1,j),c0(1,j))
         end do
         ifirst = 1
      endif
      do j=1,45
* using previously determined spline coefficients to determine 
* all 45 lam coefficients at distance r
          vsp(j)=seval(209,r,rr,vvec(1,j),a0(1,j),b0(1,j),c0(1,j))
      end do
      return
      end   
*  -----------------------------------------------------------------------
