*  -------------------------------------------------------------
      subroutine difcrs(fname1,a,ihomo,flaghf)
*  -------------------------------------------------------------
*  calculates differential cross sections
*  author h.-j. werner
*  addition for negative nnout 2-8-95 by moonbong yang
*  revised to include 8pole alignment:  7-apr-2003
*  revised to include A22+ differential alignment moment:
*     4-sep-2010 by mha
*  revised to allow differential cross sections for molecule-
*     molecule collisions (specifically ibasty=9):
*     1-feb-2012 by p.dagdigian
*  revised to allocate memory dynamically, 22-jul-2013 by q. ma
*  fixed a bug that only limited number of angles are supported,
*     22-jul-2014 by q. ma
*  use multi-dimensional arrays for atom-molecule DCSs for better
*     readability, 22-jul-2014 by q. ma
*  revised to include calculation of A21+ and A11- differential
*  polarization moments 28-jul-2014 by t. sharples
*  fixed bug in calculation of K=2 moments for j2=1,1.5 31-july-2014
*     by t. sharples
*  included option for molecule-molecule half-integer bases
*
*  current revision:  26-jun-2017 some half-int. mol-mol debugging
*                                 by j.klos
*---------------------------------------------------------------*
*                                                               *
*  IMPORTANT NOTE:  ONE SHOULD READ THE INPUT FILE FOR THE      *
*  SYSTEM BEFORE RUNNING DIFCRS, SO THAT THE LOGICAL PARAMETERS *
*  FOR THE BASIS TYPE ARE IN MEMORY                             *
*                                                               *
*---------------------------------------------------------------*
*  IMPORTANT CAVEAT:  difcrs has been extended to twomol = T basis
*  types, but not to other basis types requiring the j12 array
*  -------------------------------------------------------------
      use constants
      implicit double precision (a-h,o-z)
      character*(*) fname1
      character*20 cdate1
      character*20  cdate
      character*10  elaps, cpu
      character*40 xnam1,xnam2,xnam3
      character*1 m1string,m2string
      character*8 amplstring

      complex(8) :: stampl, stamplm
      logical existf,csflg1,flghf1,flgsu1,ihomo,flaghf,twomol,
     1        nucros,iprint,mflag,stflag,is_twomol
      include "common/parpot"
      common /coselb/ ibasty
      common /coz/ sreal1(1)
      common /cow/ simag1(1)
      common /cojq/ jq(1)
      common /colq/ lq(1)
      common /coinq/ inq(1)
      common /coj12/ j12(1)
      common /cojhld/ jlev(1)
      common /coisc1/ inlev(1)
      common /coisc2/ jout1(1)
      common /cosc1/ elev(1)
      common /cosc2/ jpack1(1)
      common /cosc3/ lpack1(1)
      common /cosc4/ ipack1(1)
      common /coj12p/ j12pk(1)
      common /codim/ mairy,mmax
      real(8), dimension(:), allocatable :: s, sm, sm6
* to store a22p,a21p and a11p amplitudes
      complex(8), dimension(:, :), allocatable :: fm1m2
      real(8), dimension(:), allocatable :: a22, a22fak, a21, a21fak,
     $     a11, a11fak
* to store rho and rho2 density matrices
      real(8), dimension(:, :), allocatable :: rho, rho2
c     scattering amplitudes
      complex(8), dimension(:, :, :), allocatable :: q, qm
c     m-dependent integral cross sections
      real(8), dimension(:, :), allocatable :: ximdep
      real(8), dimension(:), allocatable :: ytmp
      real(8), dimension(:, :, :), allocatable :: y
      dimension a(15)
c
      maxq=mmax*mmax/2
      maxy=nv2max
c
c.....input parameters
c
      j1=a(1)
      in1=a(2)
      j2=a(3)
      in2=a(4)
      ang1=a(5)
      ang2=a(6)
      dang=a(7)
      ienerg=a(8)+0.1d0
      jtot1=0
      jtot2=a(9)
      iprint = .false.
      if (a(10) .ne. 0.d0) iprint = .true.
      mflag=.false.
* if mflag = .true., calculate mdependent cross sections in collision frame
      if (a(11).ne.0) mflag=.true.
      stflag=.false.
* if stflag = .true., calculate steric effect cross sections
      if (a(12).ne.0) then
         stflag=.true.
         alph1=a(13)
         alphm1=a(14)
         msteric=a(15)
         in1=iabs(in1)
      endif
      if (stflag) then
         if (.not.flaghf) then
             write (6,5)
5            format
     :   (' *** STFLAG = .TRUE. NOT ALLOWED FOR FLAGHF = .FALSE.')
             return
         else if (ibasty.ne.3) then
             write (6,7) ibasty
7            format
     :   (' *** STFLAG = .TRUE. NOT ALLOWED FOR IBASTY =',i3)
             return
         endif
      endif
      if (flaghf) then
         xj2=j2+0.5d0
      else
         xj2=j2
      endif
*
      ang0=ang1
      if(ang2.eq.0) ang2=ang1+60.0d0
      if(dang.eq.0) dang=1.0d0
      if(ienerg.le.0) ienerg=1
c
c.....open s-matrix file
c
      call gennam(xnam1,fname1,ienerg,'smt',lenx)
      inquire (file=xnam1(1:lenx), exist=existf)
      if (.not. existf) then
        write (6, 10) xnam1(1:lenx)
10      format(/' FILE ',(a),' NOT FOUND')
        return
      end if
* initialize timer
      call mtime(cpu0,ela0)

      call openf(1,xnam1(1:lenx),'tu',0)
c
c.....open output file for differential cross sections
c
c     ifil=1
20    call gennam(xnam2,fname1,ienerg,'dcs',lenx)
      call openf(2,xnam2(1:lenx),'sf',0)
c.....open output file for final-m dependence of (m,m) and (m,m+2) density matrix
      if (mflag) then
         call gennam(xnam3,fname1,ienerg,'rho',lenx)
         call openf(3,xnam3(1:lenx),'sf',0)
      endif
* write a header (mha 10/22/08, don't write this for easier import into matlab)
*      call version(2)
      call dater(cdate)
      write (2,23) cdate
      if (mflag) write (3,23) cdate
      write (6,22) cdate
22    format (/,'**  DIFFERENTIAL CROSS SECTION (ANG^2/SR)',
     :    /,'    TODAYS DATE:  ',(a))
23    format (/,'% **  DIFFERENTIAL CROSS SECTION (ANG^2/SR)',
     :    /,'%   TODAYS DATE:  ',(a))
c
c.....read header of s-matrix file
c
      call rdhead(1,cdate1,ered1,rmu1,csflg1,flghf1,
     1   flgsu1,
     1   twomol,nucros,jfirst,jfinal,jtotd,numin,numax,nud,
     1   nlevel,nlevop,nnout,jlev,inlev,elev,jout1)
      if(csflg1) then
        write(6,30)
30      format(/' DIFFERENTIAL CROSS SECTIONS NOT IMPLEMENTED FOR',
     1    ' COUPLED STATES APPROXIMATION'/)
        goto 500
      end if
c
c.....print job information
c
      write(2,41) xnam1,cdate1,label,potnam,econv*ered1,xmconv*rmu1
      if (mflag)
     :     write(3,41) xnam1,cdate1,label,potnam,econv*ered1,xmconv*rmu1
      write(6,40) xnam1,cdate1,label,potnam,econv*ered1,xmconv*rmu1
40    format('    S-MATRICES READ FROM FILE ',(a),/,
     : '      WRITTEN:   ',(a),/,
     : '      LABEL:     ',(a)/,
     : '      POT NAME:  ',(a)/,
     : '    E-TOT:   ',f10.3,';  REDUCED MASS:  ',1pg11.4)
41    format('%    S-MATRICES READ FROM FILE ',(a),/,
     : '%     WRITTEN:   ',(a),/,
     : '%     LABEL:     ',(a)/,
     : '%     POT NAME:  ',(a)/,
     : '%    E-TOT:   ',f10.3,';  REDUCED MASS:  ',1pg11.4)
      if(jtot2.eq.0) jtot2=jfinal
      jtot2=min0(jfinal,jtot2)
      if(jtotd.ne.1) then
        write(6,50) jtotd
50      format(' *** JTOTD =',i2,'.NE.1')
        goto 500
      end if
c
c.....check whether initial level exists
c
      do 60 j=1,nlevel
60    if(jlev(j).eq.j1.and.inlev(j).eq.in1) goto 90

      write(6,70) 'j1=',j1,'in1=',in1
70    format(/1x,(a),i2,2x,(a),i2,' NOT FOUND IN LEVEL LIST'//
     1           '   N  J  INDEX   EINT(CM-1)')
      write(6,80) (j,jlev(j),inlev(j),elev(j)*econv,j=1,nlevel)
80    format(1x,2i3,i5,f12.3)
c     check if other symmetry doublet exists for steric effect
82    if (stflag) then
         do 85 j=1,nlevel
85       if(jlev(j).eq.j1.and.inlev(j).eq.-in1) goto 90
         write(6,70) 'j1=',j1,'in1=',-in1
         write(6,80) (j,jlev(j),inlev(j),elev(j)*econv,j=1,nlevel)
      endif
      goto 500
c
c.....ca is wavevector for initial state, ecol is collision energy
c
* initialize for integral cross sections
90    xint=0d0
      xintm=0d0
      ecol= ered1-elev(j)
      ca=sqrt(2d0*rmu1*ecol)
c
      do 100 j=1,iabs(nnout)
100   if(j1.eq.jout1(j))  goto 120
      write(6,110) 'j1=',j1,(jout1(j),j=1,iabs(nnout))
110   format(/1x,(a),i2,' NOT FOUND IN JOUT LIST'/
     1        1x,'JOUT: ',20i3)
      goto 500
c
c.....ckeck whether second level exists
c
c added for negative nnout by moonbong yang 2-8-95

120   continue
*     write(6,*) 'nnout=',nnout
      if(nnout.le.0) goto 160

      do 130 j=1,nlevel
130   if(jlev(j).eq.j2.and.inlev(j).eq.in2) goto 140
      write(6,70) 'j2=',j2,'in2=',in2
      write(6,80) (j,jlev(j),inlev(j),elev(j)*econv,j=1,nlevel)
      goto 500
140   do 150 j=1,iabs(nnout)
150   if(j2.eq.jout1(j)) goto 160
      write(6,110) 'j2=',j2,(jout1(j),j=1,iabs(nnout))
      goto 500
c
c.....print header
c
160   write(2,171) j1,in1,j2,in2,jtot1,jtot2,ca,1.8897*ca,
     :             econv*ecol,econv*ecol/8.065465
      if (mflag) write(3,171) j1,in1,j2,in2,jtot1,jtot2,ca,1.8897*ca,
     :             econv*ecol,econv*ecol/8.065465
      write(6,170) j1,in1,j2,in2,jtot1,jtot2,ca,1.8897*ca,
     :             econv*ecol,econv*ecol/8.065465
170   format('    TRANSITION:  J1=',i4,' IN1=',i4,
     1        '  ->  J2=',i4,' IN2=',i4/
     1        '    SUMMING PARTIAL WAVES FROM JTOT = ',i1,
     :        ' TO',i4/
     1        '    WAVE VECTOR IN INITIAL CHANNEL:  ',f10.4,' Bohr^(-1)',
     :        '; ',f10.4,' Angstroms^(-1)',/
     :        '    COLLISION ENERGY: ',f10.4,' cm-1; ',f10.4,' meV')
171   format('%   TRANSITION:  J1=',i4,' IN1=',i4,
     1        '  ->  J2=',i4,' IN2=',i4/
     1        '%   SUMMING PARTIAL WAVES FROM JTOT = ',i1,
     :        ' TO',i4/
     1        '%   WAVE VECTOR IN INITIAL CHANNEL:  ',f10.4,' Bohr^(-1)',
     :        '; ',f10.4,' Angstroms^(-1)',/
     :        '%   COLLISION ENERGY: ',f10.4,' cm-1; ',f10.4,' meV')
      if (stflag) then
         write (6,173) msteric+0.5, alph1,alphm1
         write (2,172) msteric+0.5, alph1,alphm1
      endif
172   format('%    STERIC (ORIENTED) CROSS SECTIONS: M = ',f4.1,/,
     :        '        ALPH(1) =',f7.4,
     :        ', ALPH(-1)=',f7.4)
173   format('    STERIC (ORIENTED) CROSS SECTIONS: M = ',f4.1,/,
     :        '        ALPH(1) =',f7.4,
     :        ', ALPH(-1)=',f7.4)

      if (.not.iprint) write (6, 175) xnam2(1:lenx)
175   format('    DIFFERENTIAL CROSS SECTIONS SAVED IN FILE: ',(a))
      if (.not.iprint.and.mflag) write (6, 176) xnam3(1:lenx)
176   format('    ROTATIONAL DENSITY MATRIX ELEMENTS SAVED IN FILE: ',
     :       (a))
c
c.....determine number of angles per batch
c
         l2max=jtot2+j2+1
         mlmax=j1+j2+1
         ideg1=2*j1+1
         ideg2=2*j2+1
         j1p=j1
         j2p=j2
         if(flaghf) then
c.....here for half-integer spin
           l2max=l2max+1
           mlmax=mlmax+1
           ideg1=2*j1+2
           ideg2=2*j2+2
           j1p=j1+1
           j2p=j2+1
         end if
         ideg=ideg1*ideg2
 179     nangle=(ang2-ang1)/dang+1.4d0
         if(nangle.eq.0) then
            write(6,180) maxq,maxy
 180        format(' *** NOT ENOUGH CORE IN DIFCRS. MAXQ=',i6,
     :           ' MAXY=',i6)
            return
         end if
c
c     Allocate and zero out m-dependent integral cross sections
      if (allocated(ximdep)) deallocate(ximdep)
      allocate(ximdep(-j1p:j1, -j2p:j2))
      ximdep = 0
c     Allocate and zero out amplitudes
      if (allocated(q)) deallocate(q)

      if (is_twomol(ibasty)) then
         iq1min = 1
         iq1max = (2 * j1 / 10 + 1) * (2 * mod(j1, 10) + 1)
     $        * (2 * j2 / 10 + 1) * (2 * mod(j2, 10) + 1) * nangle
c jk Added for half-integer molecule-molecule
      if(flaghf) then
         iq1max = (2 * j1 / 10 + 2) * (2 * mod(j1, 10) + 1)
     $        * (2 * j2 / 10 + 2) * (2 * mod(j2, 10) + 1) * nangle
      endif
c jk
         iq2min = 1
         iq2max = 1
         iq3 = 1
      else
         iq1min = -j1p
         iq1max = j1
         iq2min = -j2p
         iq2max = j2
         iq3 = nangle
         if (allocated(qm)) deallocate(qm)
         allocate(qm(iq1min:iq1max, iq2min:iq2max, iq3))
         qm = 0
         if (allocated(a22)) deallocate(a22)
         allocate(a22(nangle))
         a22 = 0
         if (allocated(a22fak)) deallocate(a22fak)
         allocate(a22fak(-j2p:j2p))
         if (allocated(a21)) deallocate(a21)
         allocate(a21(nangle))
         a21 = 0
         if (allocated(a21fak)) deallocate(a21fak)
         allocate(a21fak(-j2p:j2p))
         if (allocated(a11)) deallocate(a11)
         allocate(a11(nangle))
         a11 = 0
         if (allocated(a11fak)) deallocate(a11fak)
         allocate(a11fak(-j2p:j2p))
*     rho is array to accumulate |f(j1m1->j2m2)|^2 summed over m1
         if (allocated(rho)) deallocate(rho)
         allocate(rho(-j2p:j2p, nangle))
         rho = 0
*     rho2 is array to accumulate Ref(j1m1->j2m2)|^2 summed over m1
         if (allocated(rho2)) deallocate(rho2)
         allocate(rho2(-j2p:j2p, nangle))
         rho2 = 0
*     fm1m2 is array to accumulate f(j1m1->j2m2)
         if (allocated(fm1m2)) deallocate(fm1m2)
         allocate(fm1m2(-j2p:j2p, nangle))
         fm1m2 = 0
         if (allocated(sm)) deallocate(sm)
         allocate(sm(nangle))
         sm = 0
         if (allocated(sm6)) deallocate(sm6)
         allocate(sm6(nangle))
         sm6 = 0
      end if
      allocate(q(iq1min:iq1max, iq2min:iq2max, iq3))
      q = 0
      if (allocated(s)) deallocate(s)
      allocate(s(nangle))
      s = 0
c
      iydim1 = l2max - 1
      iydim2 = mlmax - 1
      if (allocated(y)) deallocate(y)
      allocate(y(0:iydim1, 0:iydim2, nangle))

c
c.....precalculate all required spherical harmonics
      if (allocated(ytmp)) deallocate(ytmp)
      allocate(ytmp(l2max))
      do ml = 0, mlmax-1
         angle=ang1
         do i = 1, nangle
            call sphn(ml, l2max - 1, angle, ytmp, 1)
            y(0:l2max-1, ml, i) = ytmp(1:l2max)
            angle=angle+dang
         end do
      end do
      if (allocated(ytmp)) deallocate(ytmp)
c
      jtlast=-1
      jplast=0
c
c.....read next s-matrix
c
250   call sread (0,sreal1, simag1, jtot, jlpar, nu1,
     :                  jq, lq, inq, ipack1, jpack1, lpack1,
     :                  1, mmax, nopen1, lengt1, ierr)
      if(ierr.eq.-1) then
         write(6,260) xnam1,jtlast,jplast
260      format(' END OF FILE DETECTED READING FILE ',(a),
     :     ' LAST JTOT,JLPAR PROCESSED:',2i5)
         goto 310
      end if
      if(ierr.lt.-1) then
        write(6,270) xnam1,jtlast,jplast
270     format(' ERROR READING FILE ',(a),
     :     ' LAST JTOT,JLPAR PROCESSED:',2i5)
        goto 310
      end if
c
c.....this assumes that jlpar=1 is stored first
c
      if(jtot.gt.jtot2) goto 300
c
c.....copy row labels into column labels if s-matrices are stored
c.....triangular
c
      if(jlpar.eq.jplast.and.jtot.ne.jtlast+1) write(6,275) jtot,jtlast
275   format(' *** WARNING: JTOT.NE.JTLAST+1:',2i4)
      jtlast=jtot
      jplast=jlpar
      if(nnout.gt.0) then
         do 290 i=1,lengt1
         inq(i)=ipack1(i)
         jq(i)=jpack1(i)
290      lq(i)=lpack1(i)
         nopen1=lengt1
      end if
c
c.....calculate contributions to amplitudes for present jtot
c
      call ampli(j1,in1,j2,in2,jtot,sreal1,simag1,mmax,jpack1,lpack1,
     :     ipack1,lengt1,jq,lq,inq,nopen1,y,q,l2max,nangle,ihomo,flaghf,
     $     iydim1,iydim2,iq1min,iq1max,iq2min,iq2max,iq3)
c.... calculate contributions for negative initial index
      if (stflag)
     :     call ampli(j1,-in1,j2,in2,jtot,sreal1,simag1,mmax,
     $     jpack1,lpack1,ipack1,lengt1,jq,lq,inq,nopen1,y,qm,
     :     l2max,nangle,ihomo,flaghf,iydim1,iydim2,iq1min,iq1max,
     $     iq2min,iq2max,iq3)
c
c.....loop back to next jtot/jlpar
c
300   if(jtot.lt.jtot2.or.jlpar.eq.1) goto 250
c
c.....print differential cross sections for this batch of angles
c
310   angle=ang1
      fak=ang2c/(ideg1*ca**2)
      faksq=sqrt(fak)
*
*  atom-molecule case
      if (.not. is_twomol(ibasty)) then
      do 343 mj1=-j1p,j1
      do 340 mj2=-j2p,j2
         if (flaghf) then
            xm=mj1+0.5d0
            xmj2=mj2+0.5d0
         else
            xm=mj1
            xmj2=mj2
         endif

* from C. H. Greene and R. N. Zare, J. Chem. Phys. 78, 6741 (1983)

* quadrupole alignment
*  A_0^(2)= <3m^2-j*(j+1)>/(j*(j+1))
* octupole alignment
*  A_0^(4)= <3(j*(j+1))^2-6(j*(j+1))-30m^2*j(j+1)+25m^2+35m^4>/8(j(j+1))^2
* A22+ moment
*\begin{eqnarray}
*A^{(2)+}_{2}(J)&=&(-1)^J \left [ \frac{2(2J-1)(2J+1)(2J+3)}{ J(J+1)}\right ]^{1/2}\nonumber \\
*&& \times \sum_{M=-J}^{J-2} (-1)^{-M}
*\left( {\begin{array}{*{20}c}
*   J & 2 & J  \\
*   { - M-2} & {  -2} & {M+2}  \\
*\end{array}} \right) {\CMcal{R}}\left({\rho_{M+2,M}}\right )
*\end{eqnarray}

* and
*\begin{equation}
*\rho_{M'M}= {{\sum\limits_{M''}^{} {f_{J''M'' \to JM'}^* } f_{J''M'' \to JM}} \mathord{\left/
* {\vphantom {{\sum\limits_{M''}^{} {f_{J''M'' \to JM'}^* } f_{J''M'' \to JM}} {\sum\limits_{M'',M'}^{} {\left|
*{f_{J''M'' \to JM'}} \right|^2 } }}} \right.
* \kern-\nulldelimiterspace} {\sum\limits_{M'',M'}^{} {\left| {f_{J''M'' \to JM'}} \right|^2 } }}
*\end{equation}

         fjjp1=xj2*(xj2+1)
         fjjp12=fjjp1*fjjp1
         xm2=xmj2*xmj2
         algfak=(3d0*xm2/fjjp1-1d0)
         octupfak=
     :      (3*fjjp12-6*fjjp1-30*fjjp1*xm2+25*xm2+35*xm2*xm2)/
     :      (8d0*fjjp12)
         xj2=j2
         if (flaghf) then
            xj2=xj2+0.5d0
         endif
         a2term=
     :       sqrt(2d0*(2*xj2-1)*(2*xj2+1)*(2*xj2+3)/(xj2*(xj2+1)))
         threej=xf3j(xj2,2d0,xj2,-xmj2,-2d0,xmj2+2d0)
         a22fak(mj2)=a2term*threej*(-1)**(xj2-xmj2)
	 a1term=sqrt(2d0*(2*xj2+1))
	 threej=xf3j(xj2,2d0,xj2,-xmj2,-1d0,xmj2+1d0)
         a21fak(mj2)=a2term*threej*(-1)**(xj2-xmj2)
         threej=xf3j(xj2,1d0,xj2,-xmj2,-1d0,xmj2+1d0)
         a11fak(mj2)=a1term*threej*(-1)**(xj2-xmj2)
* new (older than current revision)
         threej=xf3j(xj2,2d0,xj2,-xmj2,0d0,xmj2)
         algfak=threej*(a2term/sqrt(2d0))*(-1)**(xj2-xmj2)
         if (xm.lt.0d0) then
            msign=-1
         else
            msign=1
         endif
         if (mflag) then
           if (mj1.eq.-j1p .and. mj2.eq. -j2p) then
              write (2,331)
              if (iprint) write (6,331)
331           format('%',4x,'F-REAL (first line)',7x,
     :         'F-IMAG (second line); NOT divided by (2j+1)k^2 factor')
           endif
           write(m1string,'(i1)') abs(mj1)
           if (mj1.ge.0)
     :        amplstring='fm2_'//m1string
           if (mj1.lt.0)
     :        amplstring='fm2_'//'m'//m1string
           if (mj2.eq.-j2p) then
               write (2,333) amplstring
               if(iprint) write (6,333) amplstring
           endif
333        format((a),'=[')
         endif
* Comment(JK): Loop over DCS angular grid
         aangle=ang1
         do 339 i=1,nangle
         ii=ii+1
* space frame (z-axis is initial velocity vector) scattering amplitude
* returned in q,  ordering of loops:
*      inner loop:  angle
*      middle loop:  mj2
*      outer loop:  mj1
* can save m-resolved amplitudes separately (no longer done 9/23/10)
*        if (mflag) then
*           if (iprint)
*    :      write (6,336) aangle,faksq*dreal(q(ii)),faksq*imag(q(ii))
*           write (2,337) aangle,faksq*dreal(q(ii)),faksq*imag(q(ii))
*336         format(1x,f7.2,2g15.4)
*337         format(1x,f7.2,2g15.4)
*        endif
         sn=sin(aangle*pi/180d0)
* here for steric effect
         if (stflag) then
            stampl=alph1*q(mj1, mj2, i)
     $           + msign * alphm1 * qm(mj1, mj2, i)
            stamplm=alph1*q(mj1, mj2, i)
     $           - msign * alphm1 * qm(mj1, mj2, i)
            if (msteric.ge.0 .and. abs(xm).ne.(msteric+0.5d0)) then
* REPLACED 10/22/99
*           if (msteric.ge.0 .and. mj1.ne.msteric) then
                stampl=0d0
                stamplm=0d0
            endif
            term=dreal(stampl*conjg(stampl))
            termm=dreal(stamplm*conjg(stamplm))
            s(i)=s(i)+dreal(stampl*conjg(stampl))
            sm(i)=sm(i)+dreal(stamplm*conjg(stamplm))
            xint=xint+sn*term
            xintm=xintm+sn*termm
*            xint=xint+sn*dreal(stampl*conjg(stampl))
*            xintm=xintm+sn*dreal(stamplm*conjg(stamplm))
         else
            ximdep(mj1, mj2) = ximdep(mj1, mj2)
     $           + sn * dreal(q(mj1, mj2, i) * conjg(q(mj1, mj2, i)))
            dsigterm = dreal(q(mj1, mj2, i) * conjg(q(mj1, mj2, i)))
            s(i)=s(i)+dsigterm
* accumulate Re[f(j1m1->j2m2)f(j1m1->j2,m2+2)]+Re[f(j1m1->j2m2)f(j1m1->j2,m2+2)]
            fm1m2(mj2,i)=q(mj1, mj2, i)
            rho(mj2,i)=rho(mj2,i)+dsigterm
*           print *, 'mj1,mj2,mj2+j2p+1, fm1m2, rho(mj2,i):  ',
*    :         mj1,mj2, mj2+j2p+1,fm1m2(mj2,i),
*    :         rho(mj2,i)
            sm(i)=sm(i)+algfak*dsigterm
            sm6(i)=sm6(i)+octupfak*dsigterm
         endif
         aangle=aangle+dang
339      continue
340   continue
c       print *, 'mj1 is:  ', mj1
c      if (j2.ge.1) then
         do 342 mj2=-j2p,j2-2
            do 341 i=1,nangle
c            print *, ' fm1m2(mj2,i),fm1m2(mj2+2,i):  ',
c     :                   fm1m2(mj2,i),fm1m2(mj2+2,i)
            term=
     :            (dreal(fm1m2(mj2,i))*dreal(fm1m2(mj2+2,i))+
     :             dimag(fm1m2(mj2,i))*dimag(fm1m2(mj2+2,i)))
            rho2(mj2,i)=rho2(mj2,i)+term
            a22(i)=a22(i)+a22fak(mj2)*term
c            print *, 'mj1,mj2,term,a22fak(mj2):  ',
c     :     mj1,mj2,term,a22fak(mj2)
341         continue
342      continue
         do 1342 mj2=-j2p,j2-1
            do 1341 i=1,nangle
            term_21=
     :            (dreal(fm1m2(mj2,i))*dreal(fm1m2(mj2+1,i))+
     :             dimag(fm1m2(mj2,i))*dimag(fm1m2(mj2+1,i)))
            a21(i)=a21(i)+a21fak(mj2)*term_21
1341         continue
1342      continue
c      endif
c      if (j2.ge.0) then
         do 2342 mj2=-j2p,j2-1
            do 2341 i=1,nangle
            term_11=
     :            (dreal(fm1m2(mj2,i))*dimag(fm1m2(mj2+1,i))-
     :             dimag(fm1m2(mj2,i))*dreal(fm1m2(mj2+1,i)))
            a11(i)=a11(i)+a11fak(mj2)*term_11
2341         continue
2342      continue
c      endif
      if (mflag) then
         aangle=ang1
         do i=1,nangle
* print out real and imaginary parts of scattering amplitudes
            write (2,344) aangle,(dreal(fm1m2(jj,i)),jj=-j2p,j2)
            write (2,344) aangle,(dimag(fm1m2(jj,i)),jj=-j2p,j2)
            aangle=aangle+dang
         enddo
         write (2,345)
         if (iprint) write (6,345)
      endif
343   continue
344   format(f8.2,26(1pg15.4))
345   format('];')
*
      write (2,346)
      if (mflag) write (3,346)
346   format('  ')
      if (.not.stflag) then
         write (2,348)
         if (iprint) write (6,347)
347      format
     :  ('    DEGENERACY AVERAGED DXSC (ANG^2/SR) AND ',
     :       'A20, A40, A22+, A21+ AND A11- MOMENTS')
348      format
     :  ('%   DEGENERACY AVERAGED DXSC (ANG^2/SR) AND ',
     :       'A20, A40, A22+, A21+ AND A11- MOMENTS',/,'dcs=[')
         if (mflag) write (3,349)
349      format
     :  ('%   DIAGONAL (M-FINAL) DENSITY MATRIX (FIRST ROW)',
     :        ' SECOND SUPRA-DIAGONAL (M,M+2) DENSITY MATRIX (2ND ROW)',
     :        /,'rho=[')
         do 351 i=1,nangle
         if (iprint)
     :        write(6,350)
     :        angle,s(i)*fak,sm(i)/s(i),sm6(i)/s(i),a22(i)/s(i),
     $           a21(i)/s(i),a11(i)/s(i)
         write(2,350) angle,s(i)*fak,sm(i)/s(i),sm6(i)/s(i),a22(i)/s(i),
     $        a21(i)/s(i),a11(i)/s(i)
350      format(1x,f7.2,6g15.4)
         if (mflag) then
            write(3,352) angle,(rho(ij,i)/s(i),ij=-j2p,j2)
            write(3,352) angle,(rho2(ij,i)/s(i),ij=-j2p,j2)
         endif
351      angle=angle+dang
352      format(f7.2,25g15.4)
      else
         write (2,360)
         if (iprint) write (6,360)
360      format
     :  ('    DEGENERACY AVERAGED STERIC DXSC''S: HEADS AND TAILS')
         do 370 i=1,nangle
         if (iprint) write(6,350) angle,0.5d0*s(i)*fak,0.5d0*sm(i)*fak
         write(2,361) angle,0.5d0*s(i)*fak,0.5d0*sm(i)*fak
361      format(1x,'%',f7.2,3g15.4)
370      angle=angle+dang
      endif
*      ang1=angle
      ang1=angle-dang
c
c.....loop back if not all required angles are finished
c
      if(ang1.lt.ang2-0.001d0) then
         call rdhead(1,cdate1,ered1,rmu1,csflg1,flghf1,
     1   flgsu1,
     1   twomol,nucros,jfirst,jfinal,jtotd,numin,numax,nud,
     1   nlevel,nlevop,nnout,jlev,inlev,elev,jout1)
         ang1=ang1+dang
         goto 179
      end if
* determine integral m-dependent cross sections

      if (mflag) then
         write(6,375) ang0,dang,ang2
         write(2,376) ang0,dang,ang2
375      format
     :   (/,'    M-DEPENDENT (COLLISION FRAME) INTEGRAL CROSS SECTIONS',
     :    /,9x,'ANGULAR RANGE:',f6.2,':',f5.2,':',f6.2,
     :   /,'       M    M''    XSC (ANG^2)')
376      format
     :     ('];',/,'%   M-DEPENDENT (COLLISION FRAME) ',
     :     'INTEGRAL CROSS SECTIONS',
     :     /,'%',9x,'ANGULAR RANGE:',f6.2,':',f5.2,':',f6.2,
     :     /,'%','      M    M''    XSC (ANG^2)',/,' sigm1m2=[')
         algn=0d0
         xsctot=0d0
         do 385 mj1=-j1p,j1
         do 384 mj2=-j2p,j2
            if (flaghf) then
               xmj1=mj1+0.5d0
               xmj2=mj2+0.5d0
            else
               xmj1=mj1
               xmj2=mj2
            endif
c     (3m^2-1)/(j*(j+1))
            algfak=(3d0*xmj2**2/(xj2*(xj2+1d0))-1d0)
            ximdep(mj1, mj2)=ximdep(mj1, mj2)
     $           *fak*dang*2d0*pi*pi/180d0
            xsctot=xsctot+ximdep(mj1, mj2)
            algn=algn+algfak*ximdep(mj1, mj2)
* correct for mj degeneracy factor of initial state
            ximdep(mj1, mj2)=ximdep(mj1, mj2)*ideg1
            if (flaghf) then
               write (6,378) xmj1,xmj2,ximdep(mj1, mj2)
               write (2,378) xmj1,xmj2,ximdep(mj1, mj2)
378            format(4x,2f6.1,g13.4)
            else
               write (6,379) mj1,mj2,ximdep(mj1, mj2)
               write (2,379) mj1,mj2,ximdep(mj1, mj2)
379            format(4x,2i5,g13.4)
            endif
384      continue
385      continue
         write (2,387) xsctot, algn/xsctot
         write (6,386) xsctot, algn/xsctot
386      format(/,'    DEG. AVER. XSC =',g12.3,'; ALIGNMENT =',f6.3)
387      format(/,'%   DEG. AVER. XSC =',g12.3,'; ALIGNMENT =',f6.3)

      endif
* determine integral oriented (steric) cross sections
      if (stflag) then
         xint=0.5d0*fak*xint*2d0*pi*dang*pi/180d0	
         xintm=0.5d0*fak*xintm*2d0*pi*dang*pi/180d0	
         write (6,390) ang0,dang,ang2,xint,xintm,
     :            100d0*(xint-xintm)/(xint+xintm)
         write (2,391) ang0,dang,ang2,xint,xintm,
     :            100d0*(xint-xintm)/(xint+xintm)
390   format (/,
     : '    INTEGRAL ORIENTED CROSS SECTIONS; ANGLES:',f6.2,':',f5.2,
     :   ':',f6.2,/,'    XSC-HEADS =',g13.4,'  XSC-TAILS =',g13.4,
     :   '    ASYMMETRY(%) =',f8.2)
391   format(/,'%   DEG. AVER. XSC =',f8.3,'; ALIGNMENT =',f6.3)
      endif
*
*  molecule-molecule collisions
*
*  compute degeneracy averaged cross sections only
      else
        ii = 0
        j1_i = j1/10
        j1_ip=j1_i
        j2_i = mod(j1,10)
        j1_f = j2/10
        j1_fp=j1_f
        j2_f = mod(j2,10)
        ideg1 = (2*j1_i + 1)*(2*j2_i + 1)
c jk Added projections and degeneracy factor for half-integer molecule-molecule
        if (flaghf) then
        j1_ip=j1_i+1
        j1_fp=j1_f+1
        ideg1 = (2*j1_i + 2)*(2*j2_i + 1)
        endif
c jk
        fak = ang2c/(ideg1*ca**2)
        do 1500 mj1_i = -j1_ip, j1_i
        do 1500 mj2_i = -j2_i, j2_i
        do 1500 mj1_f = -j1_fp, j1_f
        do 1500 mj2_f = -j2_f, j2_f
          do 1400 i = 1, nangle
            ii =ii + 1
            dsigterm = dreal(q(ii, 1, 1) * conjg(q(ii, 1, 1)))
            s(i) = s(i) + dsigterm
1400      continue
1500    continue
*
        write(2,1348)
        if (iprint) write(6,1347)
1347    format(/'   DEGENERACY AVERAGED DXSC (ANG^2/SR)')
1348    format(/'%  DEGENERACY AVERAGED DXSC (ANG^2/SR)'/'dcs = [')
        angle = ang1
        do 1351 i = 1, nangle
          write(2,1350) angle, s(i)*fak
          if (iprint) write(6,1350) angle, s(i)*fak
1350      format(1x,f7.2,1g15.4)
          angle = angle + dang
1351    continue
      end if
*
 500  if (allocated(y)) deallocate(y)
      if (allocated(s)) deallocate(s)
      if (allocated(sm)) deallocate(sm)
      if (allocated(sm6)) deallocate(sm6)
      if (allocated(fm1m2)) deallocate(fm1m2)
      if (allocated(a22)) deallocate(a22)
      if (allocated(a22fak)) deallocate(a22fak)
      if (allocated(a21)) deallocate(a21)
      if (allocated(a21fak)) deallocate(a21fak)
      if (allocated(a11)) deallocate(a11)
      if (allocated(a11fak)) deallocate(a11fak)
      if (allocated(rho)) deallocate(rho)
      if (allocated(rho2)) deallocate(rho2)
      if (allocated(q)) deallocate(q)
      if (allocated(qm)) deallocate(qm)
      if (allocated(ximdep)) deallocate(ximdep)
c
      call mtime(cpu1,ela1)
      cpu1 = cpu1 - cpu0
      ela1 = ela1 - ela0
      call gettim(cpu1,cpu)
      call gettim(ela1,elaps)
      write(6,720) elaps, cpu
      write(2,721) elaps, cpu
720   format(/,
     : ' ** DIFFERENTIAL CROSS SECTION CALCULATION FINISHED:',
     :       /,'    ELAPSED TIME:',(a),'  CPU TIME: ',(a))
721   format(/,'];',/
     : '% ** DIFFERENTIAL CROSS SECTION CALCULATION FINISHED:',
     :       /,'%   ELAPSED TIME:',(a),'  CPU TIME: ',(a))
      close(1)
      close(2)
      if (mflag) then
         write(3,721) elaps, cpu
         close(3)
      endif
      return
      end
* -----------------------------------------------------------------------
      subroutine ampli(j1,in1,j2,in2,jtot,sreal,simag,mmax,jpack,lpack,
     :     ipack,length,jq,lq,inq,nopen,y,q,maxl2,nangle,ihomo,flaghf,
     $     iydim1,iydim2,iq1min,iq1max,iq2min,iq2max,iq3)
csubr  calculates scattering amplitudes for given jtot and set
csubr  of angles
c
*  author h.-j. werner
*  modified to allow calculation for molecule-molecule collisions
*  (ibasty=9):  1-feb-2012 by p. dagdigian
*  modified to dynamic arrays, simplified some statements: 22-jul-2014
*     by q. ma
*  included option for molecule-molecule half-integer bases
*
*  current revision date:  23-jun-2017 by p. dagdigian
*
c.....jpack,lpack,ipack: labels for rows
c.....jq,lq,inq:         labels for columns
*
      implicit double precision (a-h,o-z)
      complex*16 ai,yy,tmat
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      logical ihomo,flaghf,elastc,is_twomol
      common /coj12/ j12(1)
      common /coj12p/ j12pk(1)
      common /coselb/ ibasty
      dimension jpack(1),lpack(1),ipack(1),jq(1),lq(1),inq(1)
      dimension sreal(mmax,1),simag(mmax,1)
      integer, dimension(:), allocatable :: ilab1
      real(8), dimension(:), allocatable :: fak1
      complex(8), dimension(:), allocatable :: fak2, fak3
      real(8), dimension(0:iydim1, 0:iydim2, nangle) :: y
      complex(8), dimension(iq1min:iq1max, iq2min:iq2max, iq3) :: q
      sqpi=1.772453850905516d0
c
      ai=cmplx(zero,one)
      elastc=j1.eq.j2.and.in1.eq.in2
c
      if(flaghf) then
c.....here for half-integer spin
        spin=0.5d0
        if (.not. is_twomol(ibasty)) then
          fakj=sqpi*(two*dble(jtot) + two)*(-1)**(j1+j2+1)
          j1p=j1 + 1
          j2p=j2 + 1
        else
          j1_i = j1/10
          j1_ip = j1_i + 1
          j2_i = mod(j1,10)
          j1_f = j2/10
          j1_fp = j1_f + 1
          j2_f = mod(j2,10)
          fakj = sqpi*(two*dble(jtot) + two)
     :        *(-1)**(j1_i + j1_f + j2_i + j2_f)
        end if
      else
c.....here for integer spin
        spin=0.0d0
        if (.not. is_twomol(ibasty)) then
          fakj=sqpi*(two*dble(jtot) + one)*(-1)**(j1+j2)
          j1p=j1
          j2p=j2
        else
          j1_i = j1/10
          j1_ip=j1_i
          j2_i = mod(j1,10)
          j1_f = j2/10
          j1_fp=j1_f
          j2_f = mod(j2,10)
          fakj = sqpi*(two*dble(jtot) + one)
     :        *(-1)**(j1_i + j1_f + j2_i + j2_f)
        end if
      end if
      xjtot=dble(jtot) + spin
      if (.not. is_twomol(ibasty)) then
        xj1=dble(j1) + spin
        xj2=dble(j2) + spin
      else
        xj1_i = j1_i + spin
        xj2_i = j2_i
        xj1_f = j1_f + spin
        xj2_f = j2_f
      end if
c
      llmax = 0
      do ilab=1,length
         if(jpack(ilab).ne.j1.or.ipack(ilab).ne.in1) cycle
         llmax = llmax + 1
      end do
c
      allocate(ilab1(llmax))
      allocate(fak1(llmax))
      allocate(fak2(llmax))
      allocate(fak3(llmax))
      ll = 0
      do ilab=1,length
         if(jpack(ilab).ne.j1.or.ipack(ilab).ne.in1) cycle
         l1 = lpack(ilab)
         ll = ll + 1
         fak1(ll)=fakj*sqrt(two * dble(l1) + one)
         ilab1(ll)=ilab
      end do
c
      do 500 jlab=1,nopen
      if(jq(jlab).ne.j2.or.inq(jlab).ne.in2) goto 500
      l2=lq(jlab)
      xl2=l2
      if (is_twomol(ibasty)) then
        j12_f = j12(jlab)
        xj12_f = j12_f + spin
      end if
      do 60 ll=1,llmax
      ilab=ilab1(ll)
      l1=lpack(ilab)
      if (is_twomol(ibasty)) j12_i = j12pk(ilab)
c.....convert to t-matrix
      tmat=-cmplx(sreal(jlab,ilab),simag(jlab,ilab))
      if (.not.is_twomol(ibasty)) then
        if (elastc .and. l1.eq.l2) tmat=tmat+1.0d0
      else
        if (elastc .and. l1.eq.l2 .and.
     :    j12_i.eq.j12_f) tmat=tmat+1.0d0
      end if
      fak2(ll)=cmplx(fak1(ll)*(-1)**(l1+l2),zero)*(ai**(l1-l2))*tmat
60    continue
*
*  summation over magnetic quantum numbers
      if (.not. is_twomol(ibasty)) then
*  atom-molecule case here
        do 400 mj1=-j1p,j1
        xmj1=dble(mj1)+spin
c
        do 70 ll=1,llmax
        ilab=ilab1(ll)
        xl1=lpack(ilab)
70      fak3(ll)=fak2(ll)*xf3j(xj1,xl1,xjtot,xmj1,zero,-xmj1)
c
        do 300 mj2=-j2p,j2
        xmj2=dble(mj2)+spin
        ml2=mj1-mj2
        iyof=(iabs(ml2)*maxl2+l2)*nangle
        xml2=ml2
        fak=xf3j(xj2,xl2,xjtot,xmj2,xml2,-xmj1)
        if(ml2.gt.0) fak=fak*(-1)**ml2
c
        q(mj1, mj2, :) = q(mj1, mj2, :) + sum(fak3)
     $       * fak * y(l2, iabs(ml2), :)
300     continue ! mj2 (jk)
400     continue ! mj1 (jk)
*
*  molecule-molecule systems
      else
        ii = 0
        do 1400 mj1_i = -j1_ip, j1_i
        do 1400 mj2_i = -j2_i, j2_i
          xmj1_i = dble(mj1_i) + spin
          xmj2_i = mj2_i
          mj12_i = mj1_i + mj2_i
          xmj12_i = xmj1_i + xmj2_i
c
          do 1070 ll = 1, llmax
            ilab = ilab1(ll)
            xl1 = lpack(ilab)
            j12_i = j12pk(ilab)
            xj12_i = j12_i + spin
            fak3(ll) = fak2(ll) * (-1)**j12_i
     :        * sqrt(two * xj12_i + one)
     :        * xf3j(xj1_i,xj2_i,xj12_i,xmj1_i,xmj2_i,-xmj12_i)
     :        * xf3j(xj12_i,xl1,xjtot,xmj12_i,zero,-xmj12_i)
1070      continue
c
          do 1300 mj1_f = -j1_fp, j1_f
          do 1300 mj2_f = -j2_f, j2_f
            xmj1_f = dble(mj1_f) + spin
            xmj2_f = mj2_f
            mj12_f = mj1_f + mj2_f
            xmj12_f = xmj1_f + xmj2_f
            ml2 = mj12_i - mj12_f
            xml2 = ml2
            iyof = (iabs(ml2) * maxl2 + l2) * nangle
            fak = (-1)**(j12_f + ml2)
     :        * sqrt(two * xj12_f + one)
     :        * xf3j(xj1_f,xj2_f,xj12_f,xmj1_f,xmj2_f,-xmj12_f)
     :        * xf3j(xj12_f,xl2,xjtot,xmj12_f,xml2,-xmj12_i)
            if (ml2.gt.0) fak = fak * (-1)**ml2
c
            do 1200 iang = 1, nangle
              yy = fak * y(l2, iabs(ml2), iang)
              ii = ii + 1
              q(ii, 1, 1) = q(ii, 1, 1) + sum(fak3) * yy
1200        continue
1300      continue
1400    continue
      end if
500   continue ! jlab (jk)
      deallocate(fak1)
      deallocate(fak2)
      deallocate(fak3)
      deallocate(ilab1)
      return
      end
* -----------------------------------------------------------------------
      subroutine sphn(m,lmax,theta,yr,incy)
c
c  subroutine to calculate the
c  spherical harmonics y  (theta,0) for  m .le. l .le. lmax
c
c  -------------------------------------------------------------------
c  variables in call list:
c    m:     magnetic projection quantum number
c    lmax:  maximum value of l (minimum value of l is m)
c    theta: polar angle
c    yr:    on return: array containing real part of spherical harmonic
c    results in yr(l*incy+1) for given l
c    for l values less than m zero is returned
c
c  -------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension yr(1)
      mm = iabs(m)
c  determine array of associated legendre polynomials
      call plm(mm,lmax,theta,yr,incy)
c
c                      (l-mm)!  mm
c  on return:  yr(j) = ------- p        [cos(theta)]
c                      (l+mm)!  l=mm+j-1
c     fpi = sqrt(1 / 4 * pi)
      fpi = 2.820947917738781d-01
      if (m .lt. 0) go to 5
c  here for positive m
      a = (( - 1.d0) ** m) * fpi
      go to 10
c  here for negative m
5     a = fpi
10    ll=1
      do 11  l = 0, mm-1
      yr(ll) = 0
11    ll=ll+incy
      do 15  l = mm, lmax
      fact = sqrt(2.d0 * l + 1.d0)
      yr(ll) = fact * yr(ll) * a
15    ll = ll + incy
      return
      end
* -----------------------------------------------------------------------
      subroutine plm(m,lmax,theta,p,incp)
c
c  generates alexander-legendre polynomials for m.le.l.le.lmax
c
      implicit double precision (a-h,o-z)
      dimension p(1)
      data zero, one,  two,        rad
     :     /0.d0, 1.d0, 2.d0, 57.29577951308232d0/
      x=cos(theta/rad)
      ll = m*incp + 1
      if (m.ge.0) go to 1
      write (6,100)
100   format('  NEGATIVE M IN LEGENDRE ROUTINE:  ABORT')
      call exit
1     if (m.gt.0) go to 5
c  here for regular legendre polynomials
      p(ll)=one
      ll=ll+incp
      pm1=one
      pm2=zero
      do 2 l=1,lmax
      pp=((two*l-one)*x*pm1-(l-one)*pm2)/dble(l)
      p(ll)=pp
      ll=ll+incp
      pm2=pm1
2     pm1=pp
      return
c
c  here for alexander-legendre polynomials
c
5     imax=2*m
      rat=one
      do 6 i=2,imax,2
      ai=i
6     rat=rat*((ai-one)/ai)
      y=sin(theta/rad)
      pm1=sqrt(rat)*(y**m)
      p(ll)=pm1
      ll=ll+incp
      pm2=zero
      low=m+1
      do 10 l=low,lmax
      al=(l+m)*(l-m)
      al=one/al
      al2=((l+m-1)*(l-m-1))*al
      al=sqrt(al)
      al2=sqrt(al2)
      pp=(two*l-one)*x*pm1*al-pm2*al2
      p(ll)=pp
      ll=ll+incp
      pm2=pm1
10    pm1=pp
      return
      end
