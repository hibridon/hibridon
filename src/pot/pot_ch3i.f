* shapiro CH3I PES's modified by Guo and Schatz
* References:  M. Shapiro, J. Phys. Chem. 90, 3644 (1986);
*  H. Guo and G. C. Schatz, J. Chem. Phys. 93, 393 (1990);

* programmed by Claire Rist

*  on return:
*  vv0 contains the isotropic term (n=0) in the potential (this is zero here)
*  vvl(1) contains the isotropic (1,1) potential
*  vvl(2) contains the isotropic (2,2) potential
*  vvl(3) contains the isotropic (1,2)=(2,1) coupling potential
*  vvl(4) contains the first anisotropic (1,1) potential
*  vvl(5) contains the second anisotropic (1,1) potential
*  vvl(6) contains the first anisotropic (2,2) potential
*  vvl(7) contains the second anisotropic (2,2) potential
      subroutine driver
      implicit double precision (a-h,o-z)
      logical ifull
      dimension wf(16)
      common /covvl/ vvl(7)
      common /coered/ ered, rmu
      common /coground/ ifull
      common /covib/ie(50), iv(50)
      common /coconv/ econv, xmconv
      common /cosysi/ nscode,isicod,iscod(3)
      xmconv=5.485930d-4
      iscod(3)=1
      nch=16
      rmu=13.42/xmconv
1      print *, ' r (bohr)'
      rshift=0.5
      xfact=0.8
      do i=1,8
        ie(i)=1
        iv(i)=i-1
        iv(i+8)=i-1
        ie(i+8)=2
      enddo
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vvl',/,7(1pe16.8))
      ifull=.false.
      call ground(wf,r,nch,1,1)
      print *, wf
      goto 1
99    end
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
*  dummy loapot subroutine
       character*(*) filnam
      include "common/parpot"
      include "common/parbas"
      common /coselb/ ibasty
      potnam='SHAPIRO-GUO-SCHATZ 2D CH3I'
      ibasty=99
      lammin(1)=1
      lammax(1)=7
      mproj(1)=0
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* --------------------------------------------------------------------------
            subroutine ground(wf, rr, nch, nphoto, mxphot)

*  author: c. rist and m. alexander
*  test routine for photodissociation calculations
*  using guo/schatz two dimensional ch3-i model
*  latest revision date: 27-nov-1991
*  ---------------------------------
*     variables in call list:

*     wf        array of dimension nch*nphoto, containing, on return,
*               ground state wavefunction in each of nch components
*               nphoto is number of difference ground state wavefunctions
*     r         value of separation coordinate
*     nch       total number of channels (row dimension of q)
*     nphoto    number of different wavefunctions calculated
*               column index of q vector
*     mxphot    maximum size of q vector (mxphot .ge. nch*nphoto)
*  variables in common block /coconv/
*    econv:       conversion factor from cm-1 to hartree
*    xmconv:      conversion factor from amu to atomic units
*  variables in common /cosysr/
*    isrcod:  number of real parameters
*    junkr:   dummy variable for alignment
*    rshift:  shift of scattering potentials from nominal value
*    eel:     asymptotic electronic energy in au with respect to the
*             ground state.
*    evib:    asymptotic vibrational energy in each electronic channel (in
*             au)
*  variable in common /cosysi/
*    nscode:  total number of system dependent parameters
*             nscode = isicod + isrcod +3
*    isicod:  number of integer system dependent parameters
*    nterm:   number of different angular terms in potential
*             nterm should be 1
*    nphoto:  number of different initial ground states
*    ndip:    number of excited electronic channels (also number of dipole
*             coupling terms)
*    nvmin    minimum vibrationnal state in each electronic state
*    nvmax    maximum vibrationnal state in each electronic state
*  variables in common block /covib/
*    iv       vibrational quantum number for each asymptotic channel
*    ie       electronic quantum number for each channel
*  variables in common block /coered/
*    ered:      collision energy in atomic units (hartrees)
*    rmu:       collision reduced mass in atomic units (mass of electron = 1)

*  -------------------------------------------------------

      implicit double precision (a-h, o-z)
      parameter (ngr=21, nymx=500)
      logical ifull
      dimension wf(nch),en(ngr),f(ngr,nymx),psi(nymx),gr(ngr),
     :          dmu(2),q(4)
      include "common/parsys"
      include "common/parbas"
      common /coconv/ econv, xmconv
      common /coered/ ered, rmu
      common /cosysi/ nscode,isicod,iscod(maxpar)
      common /cosysr/ isrcod, junkr, rcod(maxpar)
      common /covib/ie(50), iv(50)
      common /coground/ ifull
      data pi4, sq2 / 0.75112554446494d0, 1.41421356237310d0/
      data y0, reg, q / 0.619702d0, 4.16799,
     :                  7.830d0, -0.1762d0, 0.6183d0, 4.939d0/
      data yymu, cf / 2.4d0, 0.036225d0/
      ymin=-2.d0
      dy=0.08
      ny=51
*     separate the variables into CH3 umbrella motion and C-I stretch
*     for each value of CH3--I distance, the component of the ground
*     state on the asymptotic vibrationnal basis must be determined.

*      if (nch*nphoto .gt. mxphot) then
*        write (6, 10) nch, nphoto,mxphot
*        write (9,10) nch, nphoto, mxphot
*10      format (' *** NCH=',i3,' * NPHOTO=',i3,' .GE. MXPHOT=',i3,
*     :          ' IN G; ABORT ***')
*        stop
*      endif
* shift ground state wavefunction
      rshift=rcod(1)
      r=rr-rshift
* define dipole moment as a function of intermolecular distance.
      ndip=iscod(3)
      dmu(1)=1.0/(1+exp(2.0*(r-9.8)))
      dmu(1)=2.d0*rmu/(1+exp(2.0*(r-9.8)))
      if (ndip .eq. 1) then
        dmu(2)=0.d0
      else
        dmu(2)=0.48304*rmu/(1+exp(2.0*(r-9.8)))
      endif

* Ground state wave-function is defined in Guo and Schatz (1990),93,393
* as a gaussian of normal coordinates: q1 and q2
* for each CH3-I distance we define the components of the ground state
* wave function in  asymptotic vibrationnal basis:

      ifull=.true.
*      entry wfintern(wf,yymin,nnvib,nny) ! original Rist statement
      entry wfintern(wf,yymin,nch,nny)
      if (.not.ifull) then
        ymin=yymin
*       nvib=nnvib/2 ! original Rist statement
        nvib=nch/2
        ny=nny
      endif

      ymu = yymu/xmconv
      alpha = (ymu * cf)**0.25
* en0 is the normalization factor for v=0
* en1 is the normalization factor for v=1
      en(1)=sqrt(alpha)*pi4
      en(2)=en(1)/sq2
      en(3)=en(1)/(2.d0*sq2)
      en(4)=en(1)/sqrt(48.d0)
      en(5)=en(1)/sqrt(384.d0)
      en(6)=en(1)/sqrt(3840.d0)
      en(7)=en(1)/sqrt(46080.d0)
      en(8)=en(1)/sqrt(645120.d0)
      en(9)=en(1)/sqrt(1.03219d+07)
      en(10)=en(1)/sqrt(1.8579d+08)
      en(11)=en(1)/sqrt(3.7159d+09)
      en(12)=en(1)/sqrt(8.1750d+10)
      en(13)=en(1)/sqrt(1.9620d+12)
      en(14)=en(1)/sqrt(5.1012d+13)
      en(15)=en(1)/sqrt(1.4283d+15)
*      en(16)=en(1)/sqrt(4.2850d+16)
*      en(17)=en(1)/sqrt(1.3712d+18)
*      en(18)=en(1)/sqrt(4.6620d+19)
*      en(19)=en(1)/sqrt(1.6782d+21)
*      en(20)=en(1)/sqrt(6.3762d+22)
*      en(21)=en(1)/sqrt(2.5493d+24)

      do i=1,ny
        y= ymin +(i-1)*dy
        q1=q(1)*(r-reg) + q(2)*(y-y0)
        q2=q(3)*(r-reg) + q(4)*(y-y0)
        psi(i)=exp(-0.5*(q1**2+q2**2))
* normalise psi:
        psi(i)=psi(i)/sqrt(0.0810)
* f(1) is the normalized h.o. function for v=0; f(2) for v=1

        xp=exp(- 0.5d0 * (alpha * y )**2)
        x = alpha*(y)
        f(1,i) = en(1) * xp
        f(2,i) = en(2) * 2.d0 * x * xp
        f(3,i) = en(3) * (4*x*x -2) * xp
        f(4,i) = en(4) * (8*x*x -12)*x * xp
        f(5,i) = en(5) * ((16*x*x -48)*x*x +12) * xp
        f(6,i) = en(6) * ((32*x*x -160)*x*x +120)*x* xp
        f(7,i) = en(7) * (((64*x*x -480)*x*x +720)*x*x-120) * xp
        f(8,i) = en(8) * (((128*x*x-1344)*x*x+3360)*x*x-1680)*x*xp
        f(9,i) = en(9) * ((((256*x*x  -3584)*x*x +13440)*x*x
     :                -13440)*x*x+1680) * xp
        f(10,i) = en(10)*((((512*x*x -9216)*x*x +48384)*x*x
     :                -80640)*x*x +30240)*x*xp
        f(11,i) = en(11)*(((((1024*x*x -23040)*x*x + 161280)*x*x
     :                -403200)*x*x + 302400)*x*x - 30240)*xp
        f(12,i) = en(12)*(((((2048*x*x - 56320)*x*x + 506880)*x*x
     :                -1774080)*x*x +2217600)*x*x -665280)*x*xp
        f(13,i) = en(13)*((((((4096*x*x - 135168)*x*x+1520640)*x*x
     :                - 7096320)*x*x + 13305600)*x*x -7983360)*x*x
     :                + 665280)*xp
        f(14,i)=en(14)*x*(17297280+x*x*(-69189120 + x*x*(69189120+
     :          x*x*(-26357760 + x*x*(4392960  +x*x*(- 319488
     :          + 8192*x*x))))))*xp
        f(15,i)=en(15)*( -17297280 + x*x*(242161920 +x*x*(-484323840+
     :          x*x*(322882560+x*x*(-92252160  + x*x*(12300288
     :          +x*x*(- 745472 + x*x*16384)))))))*xp
        f(16,i)=en(1)*x*(-2.5068229d0 + x*x*(1.1698507d+01 +x*x*(-
     :          1.4038208d+01 + x*x*(
     :          6.6848611d0  +x*x*(-1.4855247d0 +
     :          x*x*(1.6205724d-01+x*x*(
     :          - 8.3106277d-03 + x*x*1.5829767d-04)))))))*xp
        f(17,i)=en(1)*( 4.4314787d-01 +x*x*(- 7.0903660d0 +
     :           x*x*(1.6544187d+01 +x*x*(- 1.3235350d+01
     :          +x*x*(4.7269106d0 +x*x*(- 8.4033967d-01 +
     :          x*x*(7.6394515d-02 +x*x*(- 3.3580007d-03
     :           + x*x*5.5966678d-05))))))))*xp
        f(18,i)=en(1)*x*(2.5839961d0 +x*x*(- 1.3781312d1 +
     :          x*x*(1.9293837d1
     :         +x*x*(- 1.1025050d1 + x*x*(3.0625139d0*x**9
     :         +x*x*(- 4.4545657d-01
     :         +x*x*( 3.4265890d-02 +x*x*(- 1.3053672d-03+
     :           x*x*1.9196577d-05))))))))*xp
        f(19,i)=en(1)*( -4.3068141d-01 + x*x*(7.7522654d0+ x*x*(-
     :           2.0672708d1 +x*x*(
     :           1.9294527d1 +x*x*(-8.2690831d0+x*x*(1.8375740d0
     :           +x*x*(- 2.2273625d-01 + x*x*(1.4685906d-02 +x*x*(-
     :           4.8953021d-04 + x*x*6.3990877d-06)))))))))*xp
        f(20,i)=en(1)*x*(-2.6550984d+00 + x*x*(1.5930590d+01
     :          +x*x*(- 2.5488944d1 +x*x*(1.6992630d+01
     :          +x*x*(- 5.6642099d0 + x*x*(1.0298563d+00+
     :          x*x*(- 1.0562629d-01 +x*x*(
     :          6.0357881d-03 + x*x*(- 1.7752318d-04 +
     :          x*x*2.0762945d-06)))))))))*xp
        f(21,i)=en(1)*(4.1990506d-01 +x*x*(- 8.3981013d0
     :          + x*x*(2.5194304d+01 +x*x*( -2.6873924d+01  +
     :          x*x*( 1.3436962d+01  +x*x*(- 3.5831899d+00
     :          + x*x*(5.4290756d-01 +x*x*(- 4.7728137d-02 +
     :          x*x*( 2.3864068d-03 +x*x*(- 6.2389721d-05 +
     :           6.5673391d-07*x*x))))))))))*xp
      enddo
      if (ny .eq. 1) then
* here just for determination of vibrational wavefunctions
        call dcopy(nvib,f,1,wf,1)
        call dcopy(nvib,f,1,wf(nvib+1),1)
        return
      endif

      do j=1,ngr
        gr(j)=0.d0
        u=1.d0
        do i=1, ny
          wt=(3.d0+u)
          if (i.eq.1 .or. i.eq.ny) wt=1.d0
          wt=wt/3.d0
          gr(j)= gr(j) + wt*psi(i)*f(j,i)*dy
          u=-u
        enddo
      enddo
* channels correspond to electronic-vibrational states of CH3I
      do 30 i=1,nch
        iel=ie(i)
        ivb=iv(i)+1
        wf(i)= dmu(iel)*gr(ivb)
 30   continue
      return
      end
* --------------------------------------------------------------------------
      subroutine pot (vv0, r)
*  -----------------------------------------------------------------------

*  subroutine to calculate the r-dependent coefficients in the
*  model ch3i calculation to check the photodissociation program

*  on return:
*  vv0 contains the isotropic term (n=0) in the potential (this is zero here)
*  vvl(1) contains the isotropic (1,1) potential
*  vvl(2) contains the isotropic (2,2) potential
*  vvl(3) contains the isotropic (1,2)=(2,1) coupling potential
*  vvl(4) contains the first anisotropic (1,1) potential
*  vvl(5) contains the second anisotropic (1,1) potential
*  vvl(6) contains the first anisotropic (2,2) potential
*  vvl(7) contains the second anisotropic (2,2) potential
*
*  variable in common block /conlam/
*    nlam:      the total number of angular coupling terms
*    nlammx:    the maximum number of anisotropic terms
*  variable in common block /covvl/
*    vvl        array to store r-dependence of each angular term in the
*               potential
      implicit double precision (a-h, o-z)
      common /conlam/ nlam, nlammx
      common /covib/ ie(50), ic(50)
      common /covvl/ vvl(7)
*  -----------------------------------------------------------------------
      data a11, b11, a12, b12, b13, y0, r0 / 51.53d0, 1.64d0, 25.15d0,
     :          1.3d0,
     :          1.4d0, 0.619702d0, 4.16799d0/
      data cf / 0.036225d0/
      data a21, b21, b22, r02, a22, b23, a23, b24 / 0.71398d0,
     :          0.38597d0,
     :          1.5d0, 4.5d0, 0.82978d0, 1.5d0, 0.048149d0, 0.5d0/
      data bb12, g12, r012 / 0.002226d0, 0.5d0, 4.2d0/
      vvl(1) = a11*exp(-b11*r)+
     :       0.5*(cf+a12*exp(-b12*r))*y0**2*exp(-2*b13*(r-r0))
      vvl(4) = -(cf+a12*exp(-b12*r))*y0*exp(-b13*(r-r0))
*       vvl(4)=0.d0
      vvl(5) = 0.5*a12*exp(-b12*r)
*       vvl(5)=0.d0
      vvl(2) = a21*exp(-b21*r)/(1+exp(b22*(r-r02)))
      vvl(6) = a22*exp(-b23*r)
*      vvl(6)=0
      vvl(7) = a23*exp(-b24*r)
*      vvl(7)=0
*      vvl(3) = 0.00152*exp(-5.*(r-4.56)**2)
      vvl(3) = bb12*exp(-g12*(r-r012))
      return
      end
*  -----------------------------------------------------------------------
      subroutine syusr (irpot, readp, iread)
*  subroutine to read in system dependent parameters for
*  model photodissociation calculation
*  if iread = 1 read data from input file
*  if iread = 0 return after defining variable names
*  current revision date: 1-may-90
*  -----------------------------------------------------------------------
*  variables in common block /cobspt/
*    lammin:   array containing minimum value of lambda for each term
*    lammax:   array containing maximum value of lambda for each term
*    mproj:    array containing the order of the reduced rotation matrix
*              elements for each term.  lammin can not be less than mproj.
*              for homonuclear molecules, the allowed values of lambda for
*              each term range from lammin to lammax in steps of 2
*              the length of each of these arrays is limited to nterm
*              terms where nterm is defined in the common block /cosysi/
*  variables in common /cosysr/
*    isrcod:  number of real parameters
*    junkr:   dummy variable for alignment
*    rshift:  shift of scattering potentials from nominal value
*    eel:     asymptotic electronic energy in au with respect to the
*             ground state.
*    evib:    asymptotic vibrational energy in each electronic channel
*             in au)
*  variable in common /cosysi/
*    nscode:  total number of system dependent parameters
*             nscode = isicod + isrcod +3
*    isicod:  number of integer system dependent parameters
*    nterm:   number of different angular terms in potential
*             nterm should be 1
*    nphoto:  number of different initial ground states
*    ndip:    number of excited electronic channels (also number of dipole
*             coupling terms)
*    nvmin    minimum vibrationnal state in each electronic state
*    nvmax    maximum vibrationnal state in each electronic state
*  variable in common /cosys/
*    scod:    character*8 array of dimension lcode, which contains names
*             of all system dependent parameters.
*             note that the ordering of the variable names in scod must
*             correspond to the ordering of the variable names in cosysi
*             followed by the ordering of variable names in cosysr
*             followed by lammin, lammax, and mproj.

*  line 13:
*    nphoto:  total number of anistropic terms in potential
*    ndip:    total number of electronic states.
*  line 14 to 13+nel:
*    nvmin:   maximum molecular vibrationnal state included in
*              both electronic channels :
*    nvmax:   minimum molecular vibrationnal state included in
*              both electronic channels :
*    eel:     asymptotic electronic energy
*    evib:    asymptotic vibrationnal energy for each electronic
*              channel.
* line 14+nel:
*    rshift:  shift of potential surface with respect to the ground
*             state .
*  subroutines called: loapot(iunit,filnam)
*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include "common/parsys"
      integer irpot
      logical readp
      logical airyfl, airypr, logwr, swrit, t2writ, writs, wrpart,
     :        partw, xsecwr, wrxsec, noprin, chlist, ipos, flaghf,
     :        csflag, flagsu, rsflag, t2test, existf, logdfl, batch,
     :        readpt, ihomo, bastst, twomol
      character*1 dot
      character*4 char
      character*8 scod
      character*(*) fname
      character*40 filnam, line, potfil
      common /cosys/ scod(lencod)
      common /cosysi/ nscode,isicod,iscod(maxpar)
      common /cosysr/ isrcod, junkr, rcod(maxpar)
      common /conlam/ nlam
      include "common/parbas"
      common /cogrnd/ reg, caypot
      common /coqvec/ mxphot, nphoto, q
      common /colpar/ airyfl, airypr, bastst, batch, chlist, csflag,
     :                flaghf, flagsu, ihomo, ipos, logdfl, logwr,
     :                noprin, partw, readpt, rsflag, swrit,
     :                t2test, t2writ, twomol, writs, wrpart, wrxsec,
     :                xsecwr
      save potfil
      include "common/comdot"
      equivalence(iscod(1),nterm)

*     default number and names of system dependent parameters
      isicod = 7
      isrcod = 5
      readpt=.false.
      nscode = isicod + isrcod
      scod(1) = 'NTERM'
      scod(2) = 'NPHOTO'
      scod(3) = 'NDIP'
*  set default values for model dissociation problem
* typical values from Guo and Schatz (1990) 93,393, including 10
* vibrational levels for CH3.
      iscod(1)=1
      nterm = iscod(1)
      lammin(1) = 1
      lammax(1) = 7
      mproj(1) = 0
      nlam =7
      niout=20
      do i=1,niout
        if(i.le.10) then
           indout(i)=i
         else
           indout(i)=-(i-10)
         endif
      enddo
      iscod(2)=1
      nphoto=iscod(2)
      iscod(3)=2
      ndip=iscod(3)
      rcod(1) = 0.d0
      rshift= rcod(1)
      nel=2
      do i=1, nel
        iscod(3+i)=0
        iscod(4+i)=9
        rcod(2*i) = (2-i)*0.034642d0
        rcod(2*i+1) = 0.0028776d0
        char=' '
        if(nel.gt.1) then
          if(i.le.9) write(char,'(''('',i1,'')'')') i
          if(i.gt.9) write(char,'(''('',i2,'')'')') i
        end if
        scod(2+2*i)='VMIN'//char
        scod(3+2*i)='VMAX'//char
        scod(7+2*i)='EEL'//char
        scod(7+2*i+1)='EVIB'//char
      enddo
      scod(8)='RSHIFT'

      potfil = ' '
      if(iread.eq.0) return

*  read in number of ground state wavefunctions, maximum and minimum
*  vibrationnal states for each electronic channel
      if(iread .ne. 0) then
        read (8, *, err=888) nphoto, ndip
      endif
      iscod(2)=nphoto
      iscod(3)=ndip
      scod(1)='NTERM'
      scod(2)='NPHOTO'
      scod(3)='NDIP'
      isicod = 3
      isrcod = 0
      nel=2
      iofr = 2*nel + 3
      scod(iofr+1)='RSHIFT'
      iofr = iofr+1
      isrcod = 1
      do i=1, nel
        if(iread.ne.0) then
          read (8, *, err=888)iel, (iscod(isicod+j),j=1,2)
          read (8, *, err=888) ( rcod(isrcod+j),j=1,2)
        endif
        char=' '
        if(nel.gt.1) then
          if(i.le.9) write(char,'(''('',i1,'')'')') i
          if(i.gt.9) write(char,'(''('',i2,'')'')') i
        end if
        scod(isicod+1)='VMIN'//char
        scod(isicod+2)='VMAX'//char
        scod(iofr+1)='EEL'//char
        scod(iofr+2)='EVIB'//char
        iofr=iofr+2
        isicod=isicod+2
        isrcod=isrcod+2
      enddo
      if (iread .ne. 0) then
        read(8, *, err=888) rcod(1)
      endif
      if(isicod+isrcod+3.gt.lencod) stop 'lencod'
      nscode=isicod+isrcod
      goto 286
* here if read error occurs
888   write(6,1000)
1000  format(/'   *** ERROR DURING READ FROM INPUT FILE ***')
      return
* --------------------------------------------------------------
      entry ptrusr (fname,readp)
      readp = .true.
186   if (readp) then
* now call loapot(iunit,filnam) routine to read potential parameters
        filnam = ' '
        call loapot(1,filnam)
      endif
286   close (8)
      rsm = 0
      irpot=1
      return
* --------------------------------------------------------------
      entry savusr (readp)
*  save input parameters for model dissociation problem
      write (8, 290) iscod(2), iscod(3)
290   format(2i4,24x,' nphoto, ndip')
      iofi = 3
      iofr = 1
      nel = 2
      do i= 1, nel
        write (8, 295)i,(iscod(iofi+j),j=1,2)
295     format (3i4, t50,'iel, vmin, vmax')
        write (8, 296) (rcod(iofr+j),j=1,2)
296     format(2f15.8,t50,'eel, evib')
        iofi=iofi+2
        iofr=iofr+2
      enddo
      write (8, 300) rcod(1)
300   format(f11.5, t50,'rshift')
      return
      end
* --------------------------------------------------------------------

      subroutine bausr (j, l, is, jhold, ehold, ishold, nlevel, nlevop,
     :                  sc1, sc2, sc3, sc4, rcut, jtot, flaghf, flagsu,
     :                  csflag, clist, bastst, ihomo, nu, numin, jlpar,
     :                  n, nmax, ntop)
* --------------------------------------------------------------------
*  subroutine to determine angular coupling potential
*  for model photodissociation problem
*  authors:  millard alexander
*  current revision date:  28-apr-90
* --------------------------------------------------------------------
*  variables in call list:
*    j:        on return contains rotational quantum numbers for each
*              channel (dummy)
*    l:        on return contains orbital angular momentum for each
*              channel (zero)
*    is:       on return contains symmetry index for each channel
*  note that we have adopted the following convention for the symmetry
*  index:  is=iv for the first (I* channel) electronic state
*          is=-iv for the second (I channel) electronic state
*              where iv is the vibrational quantum number.
*    jhold:    on return contains rotational quantum numbers for each
*              rotational level (dummy)
*    ehold:    on return contains energy in hartrees of each level
*
*    ishold:   on return contains symmetry index of each rotational level
*    nlevel:   on return contains number of energetically distinct
*              rotational levels used in channel basis
*    nlevop:   on return contains number of energetically distinct
*              rotational levels used in channel basis which are open
*              asymptotically
*    sc1,sc2:  scratch vectors of length at least nmax
*    sc3,sc4:  scratch vectors of length at least nmax
*              these scratch vectors are not used here
*    rcut:     cut-off point for keeping higher energy channels
*              if any open channel is still closed at r=rcut, then
*              all closed channels as well any open channels which are
*              still closed at r=rcut are dropped from basis
*              note that this is ignored in molecule-surface collisions!!!
*    jtot:     total angular momentum
*    flaghf:   if .true., then system with half-integer spin
*              if .false., then system with integer spin (this is the case
*              here)
*    flagsu:   if .true., then molecule-surface collisons
*    csflag:   if .true., then coupled-states calculation
*              if .false., the close-coupled calculation
*    clist:    if .true., then quantum numbers and energies listed for
*              each channel
*    bastst:   if .true., then execution terminates after the first call
*              to basis
*    ihomo:    if .true. , then homonuclear molecule
*              if .false., then heteronuclear molecule
*              if the molecule is homonuclear (ihomo = .true.), the
*              rotational levels included go from jmin to jmax in steps
*              of 2 and only even lambda terms in the anisotropy are
*              included. Not used here.
*    nu:       coupled-states projection index
*    numin:    minimum coupled states projection index
*              for cc calculations nu and numin are both set = 0 by calling
*              program
*    jlpar:    total parity of included channels in cc calculation
*              only those channels are included for which
*                  (-1)**(parity+l-jtot)=jlpar
*              where parity is defined for each electronic, rotationnal
*              vibrationnal state, Here both sigma 2 surfaces have the
*              same symmetry and no rotation is taken into account.
*              We just define jlpar=1 for all channels involved in
*              calculation.
*              in cs calculation jlpar is set equal to 1 in calling
*              program
*    n:        on return equals number of channels
*    nmax:     maximum dimension of arrays
*    ntop:     maximum row dimension of all matrices passed to subroutines
*              propag and soutpt.  ntop is set in basis only if nu = numin
*              otherwise it is unchanged from the value supplied by the
*              calling program
*    note!!!   if flaghf = .true., then the true values of the rotational
*    quantum numbers, the total angular momentum, and the coupled-states
*    projection index are equal to the values stored in j, jtot, and nu
*    plus 1/2
*  variables in common block /cosysi/
*    nscode:   total number of system dependent variables
*    isicod:   total number of integer system dependent variables
*    nterm:    number of different associated legendre terms in
*              expansion of potential.  This
*              must be consistent with the pot routine.
*    nphoto:   number of vibrationnal ground states
*    ndip:      number of electronic channels
* for each electronic states:
*    nvmin, nvmax: minimum and maximum vibrational channels
*
*  variables in common block /cosysr/
*    isrcod:  number of real parameters
*    junkr:   dummy variable for alignment
*    rshift:  shift of scattering potentials from nominal value
*    eel:     asymptotic electronic energy ( in au with respect to the ground
*             state )
*    evib:    asymptotic vibrational energy in each electronic channel
*             (in au)
*  variables in common block /covib/
*    iv:      on return, contains vibr.channel for each
*             electronic channel
*    ie:      on return contains electronic channel
*  variables in common block /cobspt/
*    lammin:   array containing minimum value of lambda for each term
*    lammax:   array containing maximum value of lambda for each term
*    mproj:    array containing the order of the reduced rotation matrix
*              elements for each term.  lammin can not be less than mproj.
*              for homonuclear molecules, the allowed values of lambda for
*              each term range from lammin to lammax in steps of 2
*  variable in common block /cocent/
*    cent:      array containing centrifugal barrier of each channel
*  variable in common block /coint/
*    eint:      array containing channel energies (in hartree)
*  variables in common block /coered/
*    ered:      collision energy in atomic units (hartrees)
*    rmu:       collision reduced mass in atomic units
*               (mass of electron = 1)
*  variable in common block /conlam/
*    nlam:      the number of angular coupling terms actually used
*    nlammx:    the maximum number of angular coupling terms
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx
*  variable in common block /cov2/
*    nv2max:    maximum core memory allocated for the v2 matrix
*    v2:        lower triangle of nonzero angular coupling matrix elements
*               only nonzero elements are stored
*  variable in common block /coiv2/
*   iv2:        matrix address of v2 matrix for each non-zero element
* --------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical ihomo, flaghf, csflag, clist, flagsu, bastst
      logical clfl
      include "common/parbas"
      include "common/parsys"
      common /coipar/ iiipar(9), iprint
      common /cosysi/ nscode,isicod,iscod(maxpar)
      common /cosysr/ isrcod, junkr, rcod(maxpar)
      common /covib/ ie(50), iv(50)
      common /coicl/ clfl
      common /cov2/ nv2max, junkv, v2(1)
      common /coiv2/ iv2(1)
      common /conlam/  nlam, nlammx,lamnum(1)
      common /cocent/ cent(1)
      common /coeint/ eint(1)
      common /coered/ ered, rmu
      common /coconv/ econv, xmconv
      common /coiscl/ iscl(40)
      dimension j(1), l(1), jhold(1), ehold(1),
     :          sc1(1), sc2(1), sc3(1),
     :          sc4(1), ishold(1), is(1)
*   econv is conversion factor from cm-1 to hartrees
*   xmconv is conversion factor from amu to atomic units
      zero = 0.d0
      one = 1.d0
      two = 2.d0
      clfl = .true.
      nel= 2
      n=0
      nlevel=0
      do i=1, nel
        nvmin=iscod(2+2*i)
        nvmax=iscod(3+2*i)
        nvib=nvmax-nvmin+1
        do k=1, nvib
          n=n+1
          nlevel=nlevel+1
          ie(n)=i
          iv(n)=nvmin+k-1
          cent(n)=0.d0
          j(n)=0
          jhold(n)=j(n)
          l(n)=0
          is(n)=(3-2*i)*(iv(n)+1)
          ishold(n)=is(n)
          eel=rcod(2*i)
          evib=rcod(2*i+1)
          eint(n)=eel + (iv(n)+0.5d0)*evib
          ehold(n)=eint(n)
        enddo
      enddo
      if (n .gt. nmax) then
        write (9, 140) n, nmax
        write (6, 140) n, nmax
140     format(/' *** NCHANNELS=', i3,' .GT. MAX DIMENSION OF',
     :         i3,'; ABORT')
        if (bastst) then
          return
        else
          call exit
        end if
      end if
*  return if no channels
      if (n .eq. 0) return
*  now sort this list in terms of increasing energy
      if (n .gt. 1) then
        do 144 i1 = 1, n - 1
          esave = ehold(i1)
          iel = ie(i1)
          do i2 = i1 + 1, n
            if (iel .ne. ie(i2)) go to 144
            if (ehold(i2) .lt. esave) then
*  state i2 has a lower energy than state i1, switch them
              esave = ehold(i2)
              ehold(i2) = ehold(i1)
              ehold(i1) = esave
              esave = eint(i2)
              eint(i2) = eint(i1)
              eint(i1) = esave
              jsave = j(i2)
              j(i2) = j(i1)
              j(i1) = jsave
              jsave = jhold(i2)
              jhold(i2) = jhold(i1)
              jhold(i1) = jsave
              lsave = l(i2)
              l(i2) = l(i1)
              l(i1) = lsave
              csave = cent(i2)
              cent(i2) = cent(i1)
              cent(i1) = csave
              issave = is(i2)
              is(i2) = is(i1)
              is(i1) = issave
              issave = ishold(i2)
              ishold(i2) = ishold(i1)
              ishold(i1) = issave
              iesave = ie(i2)
              ie(i2) = ie(i1)
              ie(i1) = iesave
              ivsave = iv(i2)
              iv(i2) = iv(i1)
              iv(i1) = ivsave
           end if
          enddo
144     continue
      end if
*
* form list of all energetically open rotational levels included in the
*  calculations and their energies (same list as previously initialised
* as no angular momentum is included).
*  now sort this list to put closed levels at end
*  also determine number of levels which are open
      nlevop = 0
      do 150  i = 1, nlevel - 1
        if (ehold(i) .le. ered) then
          nlevop = nlevop + 1
        else
          do 145 ii = i + 1, nlevel
            if (ehold(ii) .le. ered) then
              nlevop = nlevop + 1
              ikeep = jhold(i)
              jhold(i) = jhold(ii)
              jhold(ii) = ikeep
              ikeep = ishold(i)
              ishold(i) = ishold(ii)
              ishold(ii) = ikeep
              ekeep = ehold(i)
              ehold(i) = ehold(ii)
              ehold(ii) = ekeep
              go to 150
            end if
145        continue
        end if
150    continue
* Claires fooling
      do i=1,nlevop
        iscl(i)=ishold(i)
*        write(6,*)i,ishold(i),iscl(i)
      enddo

      if (nlevop .le. 0) then
        write (9, 155)
        write (6, 155)
155     format('*** NO OPEN LEVELS IN BAUSR; ABORT')
        if (bastst) return
        call exit
      endif
      if (ehold(nlevel) .le. ered) nlevop = nlevop + 1

      if (nu .eq. numin) then
        ntop = max(n, nlevel)
*  ntop is the maximum row dimension of all matrices passed in the
*  call list of subroutines propag and soutpt.
*  for fps make sure this is an odd number, for faster bank access.
*  this has no effect on vax or cray
        if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax)
     :       ntop = ntop + 1
      end if
*  now list channels if requested
      if (clist .or. bastst) then
        write (6, 255)
        write (9, 255)
255     format(/' All channels:'//
     1          '   N   EL  V  IND     EINT(CM-1)')
256     format(/' Open channels:'//
     1          '   N   EL  V  IND     EINT(CM-1)')
        do 265  i = 1, nlevel
          nv = abs(is(i))-1
          ieps = is(i)/abs(is(i))
          iel = (3-ieps)/2
          write (6, 260) i, iel, nv, is(i), eint(i) * econv
          write (9, 260) i, iel, nv, is(i), eint(i) * econv
260       format (4i4, f13.3)
265     continue
        write (6, 256)
        do 266  i = 1, nlevop
          nv = abs(ishold(i))-1
          ieps=ishold(i)/abs(ishold(i))
          iel = (3-ieps)/2
          write (6, 261) i, iel, nv, ishold(i), ehold(i) * econv
261       format (4i4, f13.3)
266     continue
      end if
*  now calculate coupling matrix elements
      if (bastst .and. (iprint .gt. 1)) then
        write (6, 280)
        write (9, 280)
280     format (/' ILAM  LAMBDA  ICOL  IROW   I    IV2    VEE')
      end if
* i counts v2 elements
* inum counts v2 elements for given lambda
* ilam counts numver of v2 matrices
* ij is address of given v2 element in present v2 matrix
      ii=0
      call ptmatrix(0,0,0,vee,ii)
      ii=1
      i = 0
      ilam=0
      do 320 il = lammin(1), lammax(1), 1
        inum = 0
        ilam=ilam+1
        lb=il
        if(ilam.gt.nlammx) then
          write(6,311) ilam
311       format(/' ILAM.GT.NLAMMX IN BAUSR')
          call exit
        end if
        do 310  icol= 1, n
          do 300  irow = icol, n
            ij = ntop * (icol - 1) +irow
            vee=zero
            ier=ie(irow)
            iec=ie(icol)
            ivr=iv(irow) +1
            ivc=iv(icol) +1
            if (il .eq. 1) then
              if (ier .eq. 1 .and. iec .eq. 1) then
                if(ivr.eq.ivc) vee=1.d0
              endif
            endif
            if (il .eq. 2) then
              if (ier .eq. 2 .and. iec .eq. 2) then
                if(ivr.eq.ivc) vee=1.d0
              endif
            endif
            if (il .eq. 3) then
              if (ier .ne. iec) then
                if(ivr .eq. ivc) vee=1.d0
              endif
            endif
            if (il.eq.4) then
              if (ier.eq.1 .and. iec.eq.1) then
                call ptmatrix(il,ivr,ivc,vee,ii)
              endif
            endif
            if (il.eq.5) then
              if (ier.eq.1 .and. iec.eq.1) then
                call ptmatrix(il,ivr,ivc,vee,ii)
              endif
            endif
            if (il.eq.6) then
              if (ier.eq.2 .and. iec.eq.2) then
                call ptmatrix(il,ivr,ivc,vee,ii)
              endif
            endif
            if (il.eq.7) then
              if (ier.eq.2 .and. iec.eq.2) then
                call ptmatrix(il,ivr,ivc,vee,ii)
              endif
            endif
            if (il.gt.7) then
              write(6,289) il
289           format('LAMBDA.GT.',i2,' UNDEFINED')
            endif
            if (vee .eq. zero) goto 300
              i = i + 1
              inum = inum + 1
              if (i .gt. nv2max) goto 325
                v2(i) = vee
                iv2(i) = ij
                if (bastst .and. (iprint .gt. 1)) then
                  write (6, 290) ilam, lb, icol, irow, i, iv2(i),
     :                           vee
                  write (9, 290) ilam, lb, icol, irow, i, iv2(i),
     :                           vee
290               format (i4, 2i7, 2i6, i6, g17.8)
                endif
300       continue
310     continue
        lamnum(ilam) = inum
        if (bastst) then
          write (6, 315) ilam, lamnum(ilam)
          write (9, 315) ilam, lamnum(ilam)
315       format ('ILAM=',i3,' LAMNUM(ILAM) = ',i3)
        end if
320   continue
325   if ( i.gt. nv2max) then
        write (6, 350) i, nv2max
        write (9, 350) i, nv2max
350     format (' *** NUMBER OF NONZERO V2 ELEMENTS = ',i6,
     :           ' .GT. NV2MAX=',i6,'; ABORT ***')
        if (bastst) then
          return
        else
          call exit
        end if
      end if
      if (clist .or. bastst .and. (iprint .ge. 0)) then
        write (6, 360) i
        write (9, 360) i
360     format (' ** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS',
     :           i6)
      end if
      return
      end


      subroutine hof(y,f,ny)
* returns the coupling matrix elements for lambda=il=4,5,6,7
* ie. matrix elements of y, y**2, exp(b23*0.20218*y) and
* exp(b24*0.2218*y).
* between  states of vibrational quantum number ir-1 and ic-1.
      implicit double precision (a-h,o-z)

      parameter (ngr=21, nymx=500)
      dimension en(ngr), f(ngr, nymx), y(ny)

      common /coconv/ econv, xmconv
      data pi4, sq2 / 0.75112554446494d0, 1.41421356237310d0/
      data yymu, cf / 2.4d0, 0.036225d0/

      ymu = yymu/xmconv
      alpha = sqrt(sqrt((ymu * cf)))
* en0 is the normalization factor for v=0
* en1 is the normalization factor for v=1
      en(1)=sqrt(alpha)*pi4
      en(2)=en(1)/sq2
      en(3)=en(1)/(2.d0*sq2)
      en(4)=en(1)/sqrt(48.d0)
      en(5)=en(1)/sqrt(384.d0)
      en(6)=en(1)/sqrt(3840.d0)
      en(7)=en(1)/sqrt(46080.d0)
      en(8)=en(1)/sqrt(645120.d0)
      en(9)=en(1)/sqrt(1.03219d+07)
      en(10)=en(1)/sqrt(1.8579d+08)
      en(10)=en(1)/sqrt(1.8579d+08)
      en(11)=en(1)/sqrt(3.7159d+09)
      en(12)=en(1)/sqrt(8.1750d+10)
      en(13)=en(1)/sqrt(1.9620d+12)
*      en(16)=en(1)/sqrt(4.2850d+16)
*      en(17)=en(1)/sqrt(1.3712d+18)
*      en(18)=en(1)/sqrt(4.6620d+19)
*      en(19)=en(1)/sqrt(1.6782d+21)
*      en(20)=en(1)/sqrt(6.3762d+22)
*      en(21)=en(1)/sqrt(2.5493d+24)

      do i=1,ny
* f(1) is the normalized h.o. function for v=0; f(2) for v=1
        xp=exp(- 0.5d0 * (alpha * y(i) )**2)
        x = alpha*(y(i))
        x2=x*x
        x3=x2*x
        x4=x3*x
        x5=x4*x
        x6=x5*x
        x7=x6*x
        x8=x7*x
        x9=x8*x
        x10=x9*x
        x11=x10*x
        x12=x11*x
        x13=x12*x
        x14=x13*x
        x15=x14*x
        x16=x15*x
        x17=x16*x
        x18=x17*x
        x19=x18*x
        x20=x19*x
        f(1,i) = en(1) * xp
        f(2,i) = en(2) * 2.d0* x * xp
        f(3,i) = en(3) * (4.d0*x2 -2.d0) * xp
        f(4,i) = en(4) * (8.d0*x3 -12.d0*x) * xp
        f(5,i) = en(5) * (16.d0*x4 -48.d0*x2 +12.d0) * xp
        f(6,i) = en(6) * (32.d0*x5 -160.d0*x3 +120.d0*x)* xp
        f(7,i) = en(7) * (64.d0*x6 -480.d0*x4 +720.d0*x2-120.d0)
     :                 * xp
        f(8,i) = en(8) * (128.d0*x7 -1344.d0*x5 +3360.d0*x3
     :                   -1680.d0*x) * xp
        f(9,i) = en(9) * (256.d0*x8  -3584.d0*x6 +13440.d0*x4
     :                -13440.d0*x2+1680.d0) * xp
        f(10,i) = en(10)*(512.d0*x9 -9216.d0*x7 +48384.d0*x5
     :                -80640.d0*x3 +30240.d0*x)*xp
        f(11,i) = en(11)*(1024.d0*x10 -23040.d0*x8 + 161280.d0*x6
     :                -403200.d0*x4 + 302400.d0*x2 - 30240.d0)*xp
        f(12,i) = en(12)*(2048.d0*x11 - 56320.d0*x9 + 506880.d0*x7
     :                -1774080.d0*x5 +2217600.d0*x3 -665280.d0*x)*xp
        f(13,i) = en(13)*(4096.d0*x12 - 135168.d0*x10+1520640.d0*x8
     :                - 7096320.d0*x6 + 13305600.d0*x4
     :                -7983360.d0*x2
     :                + 665280.d0)*xp
        f(14,i)=en(14)*(17297280d0*x - 69189120d0*x3 + 69189120d0*x5 -
     :          26357760d0*x7 + 4392960d0*x9  - 319488d0*x11
     :          + 8192d0*x13)*xp
        f(15,i)=en(15)*( -17297280 + 242161920d0*x2 - 484323840d0*x4 +
     :          322882560d0*x6 -92252160d0*x8  + 12300288d0*x10
     :          - 745472d0*x12 + 16384d0*x14)*xp
        f(16,i)=en(1)*(-2.5068229d0*x + 1.1698507d+1*x3 -
     :          1.4038208d+1*x5 +
     :          6.6848611d0*x7  -1.4855247d0*x9 +
     :          1.6205724d-1*x11
     :          - 8.3106277d-3*x13 + 1.5829767d-4*x15)*xp
        f(17,i)=en(1)*( 4.4314787d-1 - 7.0903660d0*x2 +
     :           1.6544187d+1*x4 - 1.3235350d+1*x6
     :          +4.7269106d0*x8 - 8.4033967d-1*x10 +
     :          7.6394515d-2*x12 - 3.3580007d-3*x14
     :           + 5.5966678d-5*x16)*xp
        f(18,i)=en(1)*(2.5839961d0*x - 1.3781312d1*x3 +
     :          1.9293837d1*x5
     :         - 1.1025050d1*x7 + 3.0625139d0*x9 -
     :          4.4545657d-1*x11
     :         + 3.4265890d-2*x13 - 1.3053672d-3*x15+
     :           1.9196577d-5*x17)*xp
        f(19,i)=en(1)*( -4.3068141d-1 + 7.7522654d0*x2 -
     :           2.0672708d1*x4 +
     :           1.9294527d1*x6 -8.2690831d0*x8+1.8375740d0*x10
     :           - 2.2273625d-1*x12 + 1.4685906d-2*x14 -
     :           4.8953021d-4*x16 + 6.3990877d-6*x18)*xp
        f(20,i)=en(1)*(-2.6550984d0*x + 1.5930590d+1*x3
     :          - 2.5488944d1*x5 +1.6992630d+1*x7
     :          -  5.6642099d0*x9 + 1.0298563d0*x11
     :          - 1.0562629d-1*x13 +
     :          6.0357881d-3*x15 - 1.7752318d-4*x17 +
     :          2.0762945d-6*x19)*xp
        f(21,i)=en(1)*(4.1990506d-1 - 8.3981013d0*x2
     :          + 2.5194304d+1*x4  -2.6873924d+1*x6  +
     :           1.3436962d+1*x8  - 3.5831899d0*x10
     :          + 5.4290756d-1*x12 - 4.7728137d-2*x14 +
     :           2.3864068d-3*x16 - 6.2389721d-5*x18 +
     :           6.5673391d-7*x20)*xp

      enddo
      return
      end

      subroutine  ptmatrix(il,ir,ic,v,ii)

      implicit double precision (a-h,o-z)
      parameter (ngr=21, nymx=500)
      dimension f(ngr, nymx),y(nymx)

      data b23, b24 / 1.5d0, 0.5d0/
      data ymin, dy / -2.d0, 0.01d0/
      data ny / 401/
      save y, f

      if (ii .eq. 0) then
        do i=1,ny
          y(i)= ymin +(i-1)*dy
        enddo
        call hof(y,f,ny)
        ii=1
        return
      endif
* check that vibrational wave function has been defined:
      if( ir.gt.21) then
        write(6,369)ir-1
        return
      elseif (ic.gt.21) then
        write(6,369)ic-1
        return
      endif
369   format(' VIBRATIONAL WAVE FUNCTION V=',i3,' NOT DEFINED')

      v=0.d0
      fint=0.d0
      do i=2,ny
        fint=0.d0
        yu=y(i)
        yl=yu-dy
        if(il.eq.4) then
          fu=yu
          fl=yl
        elseif (il.eq.5) then
          fu=yu**2
          fl=yl**2
        elseif (il.eq.6) then
          fu=exp(b23*0.20218d0*yu)
          fl=exp(b23*0.20218d0*yl)
        elseif (il.eq.7) then
          fu=exp(b24*0.20218d0*yu)
          fl=exp(b24*0.20218d0*yl)
        else
          write(6,370) il
370       format('LAMBDA.GT.',i2,' UNDEFINED')
        endif
        fint=0.5d0*(f(ir,i)*f(ic,i)*fu + f(ir,i-1)*f(ic,i-1)*fl)
        v = v + fint*dy
      enddo
      return
      end

