* System:  OH(X 2Pi)+Ar,  ab initio UCCSD(T) PES's
* with Complete Basis Set extrapolation 
* obtained from 3D PES by averaging over 
* OH(v=1) wavefunction
* Calculations by  J. Klos
* Department of Chemistry
* University of Maryland, College Park, MD
* Prepared on 2013.08.26
*Reference fpor v=0: L. Scharfenberg, J. Klos, P. J. Dagdigian,
*           M. H. Alexander,G. Meijer, S. Y. T. van der Meerakker
*           Phys. Chem. Chem. Phys. 12, 10660-10670 (2010)
*
*  damp long-range part of PES to an isotropic C6 form, with
*  C6 determined from Voo term at R=10
*  (26-aug-2013, p.dagdigian)


      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(19)
      include "common/parpot"
      potnam=' J.Klos Ar-OH(X) UCCSDT CBS <v=1> PES'
      print *, potnam
1      print *, ' r (bohr)'
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vsum',/,11(1pe16.8),/,
     :    '  vdif',/,5e16.8)
      goto 1
99    end

      include "common/syusr"
      include "common/bausr"
      include "common/ground"
      subroutine loapot(iunit,filnam)
* ------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      potnam='J.Klos Ar-OH(X) UCCSDT CBS <v=1> PES'
      lammin(1)=1
      lammax(1)=10
      lammin(2)=2
      lammax(2)=10
      mproj(1)=0
      mproj(2)=2
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* ----------------------------------------------------------------------
      subroutine pot (vv0, r)
*  subroutine to calculate the r-dependent coefficients
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*  vv0        contains isotropic term (d00 term in vsum)
*  variable in common block /covvl/
*    vvl:     vector of length 6 to store r-dependence of each term
*             in potential expansion
*    vvl(1-10) expansion coefficients in dl0 (l=1:10) of vsum
*    vvl(11-19) expansion coefficients in dl2 (l=2:10) of vdif

* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  8-oct-1993
* revised for Ar-OH  Theta=0 for Ar-HO
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      common /covvl/ vvl(19)
      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
      call dscal(19,zero,vvl,1)
* convert to hartree
      conv=1.d0/219474.6
      vv0=   VL0_v1_0(r)*conv
      vvl(1)=VL0_v1_1(r)*conv
      vvl(2)=VL0_v1_2(r)*conv
      vvl(3)=VL0_v1_3(r)*conv
      vvl(4)=VL0_v1_4(r)*conv
      vvl(5)=VL0_v1_5(r)*conv
      vvl(6)=VL0_v1_6(r)*conv
      vvl(7)=VL0_v1_7(r)*conv
      vvl(8)=VL0_v1_8(r)*conv
      vvl(9)=VL0_v1_9(r)*conv
      vvl(10)=VL0_v1_10(r)*conv
      vvl(11)=VL2_v1_2(r)*conv
      vvl(12)=VL2_v1_3(r)*conv
      vvl(13)=VL2_v1_4(r)*conv
      vvl(14)=VL2_v1_5(r)*conv
      vvl(15)=VL2_v1_6(r)*conv
      vvl(16)=VL2_v1_7(r)*conv
      vvl(17)=VL2_v1_8(r)*conv
      vvl(18)=VL2_v1_9(r)*conv
      vvl(19)=VL2_v1_10(r)*conv
*
* switch to isotropic C6 potential at long range
* r^-6 fit to isotropic part of potential at R=20
      c6sum = -1.1384556e+07*conv
* switching function for long-range
      rsw = 25.d0
      switch_lr=0.5*(tanh(0.5*(r - rsw)) + 1.d0)
* kill anisotropic terms at large R
      do ilam=1,19
        vvl(ilam) = (1.d0 - switch_lr)*vvl(ilam)
      end do
* isotropic term - merge with C6 form
      vv0 = (1.d0 - switch_lr)*vv0 + switch_lr*c6sum/(r**6)
      return
      end
*
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Ar-OH(X) PES averaged over v=1 of OH
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_v1_0(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=35)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .       3.5000D0,
     .       3.7500D0,
     .       4.0000D0,
     .       4.2500D0,
     .       4.5000D0,
     .       4.7500D0,
     .       5.0000D0,
     .       5.2500D0,
     .       5.5000D0,
     .       5.7500D0,
     .       6.0000D0,
     .       6.2500D0,
     .       6.5000D0,
     .       6.7500D0,
     .       7.0000D0,
     .       7.2500D0,
     .       7.5000D0,
     .       7.7500D0,
     .       8.0000D0,
     .       8.2500D0,
     .       8.5000D0,
     .       8.7500D0,
     .       9.0000D0,
     .       9.2500D0,
     .       9.5000D0,
     .       9.7500D0,
     .      10.0000D0,
     .      11.0000D0,
     .      12.0000D0,
     .      13.0000D0,
     .      14.0000D0,
     .      16.0000D0,
     .      18.0000D0,
     .      20.0000D0,
     .      25.0000D0/
      DATA VCOEF/
     *    0.114849914576047363D+11,
     *   -0.129962041662471561D+11,
     *    0.423168241019788170D+10,
     *   -0.205333620800140071D+10,
     *    0.554474907850036025D+08,
     *   -0.703888838812267900D+09,
     *   -0.433523916034938693D+09,
     *   -0.463186569711734533D+09,
     *   -0.354139156637377858D+09,
     *   -0.268905877104441524D+09,
     *   -0.169950344297000051D+09,
     *   -0.820418822709082365D+08,
     *   -0.369843038808184862D+07,
     *    0.602537113177222908D+08,
     *    0.107739640526475668D+09,
     *    0.138048658399644494D+09,
     *    0.154105818608712375D+09,
     *    0.159309078512946814D+09,
     *    0.154482023201620162D+09,
     *    0.148810651668177128D+09,
     *    0.130447699862471700D+09,
     *    0.123813198772628069D+09,
     *    0.101672428214501619D+09,
     *    0.935537088370110989D+08,
     *    0.888645776702790260D+08,
     *    0.203383502623372078D+08,
     *    0.160201514521385193D+09,
     *    0.670836338168311119D+08,
     *   -0.880303461056900024D+07,
     *   -0.802305152424097061D+07,
     *   -0.119020196152300835D+08,
     *   -0.801927502817916870D+07,
     *   -0.133221909095780849D+08,
     *    0.225242431460022926D+06,
     *   -0.104214635694338083D+08/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_v1_0=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Ar-OH(X) PES averaged over v=1 of OH
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_v1_10(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=35)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .       3.5000D0,
     .       3.7500D0,
     .       4.0000D0,
     .       4.2500D0,
     .       4.5000D0,
     .       4.7500D0,
     .       5.0000D0,
     .       5.2500D0,
     .       5.5000D0,
     .       5.7500D0,
     .       6.0000D0,
     .       6.2500D0,
     .       6.5000D0,
     .       6.7500D0,
     .       7.0000D0,
     .       7.2500D0,
     .       7.5000D0,
     .       7.7500D0,
     .       8.0000D0,
     .       8.2500D0,
     .       8.5000D0,
     .       8.7500D0,
     .       9.0000D0,
     .       9.2500D0,
     .       9.5000D0,
     .       9.7500D0,
     .      10.0000D0,
     .      11.0000D0,
     .      12.0000D0,
     .      13.0000D0,
     .      14.0000D0,
     .      16.0000D0,
     .      18.0000D0,
     .      20.0000D0,
     .      25.0000D0/
      DATA VCOEF/
     *    0.177215103432925272D+10,
     *   -0.311728321535485744D+10,
     *    0.141434956468678164D+10,
     *   -0.322399763187381268D+09,
     *    0.191250742770919919D+09,
     *   -0.130825309630335644D+08,
     *    0.408237575537267029D+08,
     *    0.867940509000888653D+07,
     *    0.936945527461042628D+07,
     *    0.422026521140757576D+07,
     *    0.401598765113616735D+07,
     *    0.442781462301395088D+07,
     *    0.100560333452760568D+07,
     *    0.202263754814216495D+07,
     *   -0.477573066690370440D+06,
     *   -0.266743072949822247D+07,
     *    0.513821555010868609D+07,
     *   -0.431003085777941346D+07,
     *    0.298936182966720313D+07,
     *    0.373796926142794080D+07,
     *   -0.698136474107010663D+07,
     *    0.620253179662559927D+07,
     *   -0.319888361507406831D+07,
     *   -0.364738510852944851D+07,
     *    0.608614687293463945D+07,
     *   -0.777367732125163078D+07,
     *    0.687134747187960148D+07,
     *   -0.247328317390811443D+07,
     *    0.268481699474167824D+07,
     *   -0.569512425483489037D+07,
     *    0.686127835524809361D+07,
     *   -0.616374904850697517D+07,
     *    0.606137910760885477D+07,
     *   -0.337481435724860430D+07,
     *    0.453720726148068905D+06/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_v1_10=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Ar-OH(X) PES averaged over v=1 of OH
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_v1_1(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=35)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .       3.5000D0,
     .       3.7500D0,
     .       4.0000D0,
     .       4.2500D0,
     .       4.5000D0,
     .       4.7500D0,
     .       5.0000D0,
     .       5.2500D0,
     .       5.5000D0,
     .       5.7500D0,
     .       6.0000D0,
     .       6.2500D0,
     .       6.5000D0,
     .       6.7500D0,
     .       7.0000D0,
     .       7.2500D0,
     .       7.5000D0,
     .       7.7500D0,
     .       8.0000D0,
     .       8.2500D0,
     .       8.5000D0,
     .       8.7500D0,
     .       9.0000D0,
     .       9.2500D0,
     .       9.5000D0,
     .       9.7500D0,
     .      10.0000D0,
     .      11.0000D0,
     .      12.0000D0,
     .      13.0000D0,
     .      14.0000D0,
     .      16.0000D0,
     .      18.0000D0,
     .      20.0000D0,
     .      25.0000D0/
      DATA VCOEF/
     *    0.702118738250757408D+10,
     *   -0.968424474874148178D+10,
     *    0.361385080217848682D+10,
     *   -0.111332569144691133D+10,
     *    0.466051414522558570D+09,
     *   -0.180366900975352407D+09,
     *   -0.560328087198958099D+08,
     *   -0.161601168279991627D+09,
     *   -0.155485610138041854D+09,
     *   -0.155170032850356340D+09,
     *   -0.129208653042756468D+09,
     *   -0.947006994848804474D+08,
     *   -0.569853353742408454D+08,
     *   -0.183107056848988011D+08,
     *    0.102719818144303858D+08,
     *    0.336344691481464207D+08,
     *    0.484664737674354613D+08,
     *    0.608616364916947037D+08,
     *    0.705528638627163768D+08,
     *    0.766075441682578623D+08,
     *    0.807712670301932693D+08,
     *    0.751776370758002996D+08,
     *    0.679275837123830318D+08,
     *    0.516974515465584993D+08,
     *    0.512749610676114559D+08,
     *   -0.193736601238250732D+05,
     *    0.918599617057712078D+08,
     *    0.108025262631044388D+08,
     *   -0.821495713985729218D+07,
     *    0.499180908605861664D+07,
     *   -0.123334528919196129D+08,
     *    0.113561624907701015D+08,
     *   -0.130264008706103563D+08,
     *   -0.586432050278687477D+07,
     *   -0.161819261073172092D+05/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_v1_1=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Ar-OH(X) PES averaged over v=1 of OH
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_v1_2(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=35)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .       3.5000D0,
     .       3.7500D0,
     .       4.0000D0,
     .       4.2500D0,
     .       4.5000D0,
     .       4.7500D0,
     .       5.0000D0,
     .       5.2500D0,
     .       5.5000D0,
     .       5.7500D0,
     .       6.0000D0,
     .       6.2500D0,
     .       6.5000D0,
     .       6.7500D0,
     .       7.0000D0,
     .       7.2500D0,
     .       7.5000D0,
     .       7.7500D0,
     .       8.0000D0,
     .       8.2500D0,
     .       8.5000D0,
     .       8.7500D0,
     .       9.0000D0,
     .       9.2500D0,
     .       9.5000D0,
     .       9.7500D0,
     .      10.0000D0,
     .      11.0000D0,
     .      12.0000D0,
     .      13.0000D0,
     .      14.0000D0,
     .      16.0000D0,
     .      18.0000D0,
     .      20.0000D0,
     .      25.0000D0/
      DATA VCOEF/
     *    0.855547342836268044D+10,
     *   -0.105361646993892746D+11,
     *    0.358089640504435730D+10,
     *   -0.149349331029282951D+10,
     *    0.297791482617483020D+09,
     *   -0.358045993467365861D+09,
     *   -0.184521709249399006D+09,
     *   -0.260484952691158712D+09,
     *   -0.214621874960084826D+09,
     *   -0.178243015115898222D+09,
     *   -0.120937739320218384D+09,
     *   -0.669274621238223910D+08,
     *   -0.168320588444537818D+08,
     *    0.272657213520122617D+08,
     *    0.614313294789564237D+08,
     *    0.837831524167415202D+08,
     *    0.937870364367666990D+08,
     *    0.971327061425242126D+08,
     *    0.965520490800392032D+08,
     *    0.977611919465008378D+08,
     *    0.925453506089751720D+08,
     *    0.817663764482641220D+08,
     *    0.753185992407371998D+08,
     *    0.501179488562288284D+08,
     *    0.598946164654572010D+08,
     *    0.739481794965362549D+07,
     *    0.102192847847032547D+09,
     *   -0.141269444881777763D+08,
     *   -0.985315728152322769D+07,
     *    0.949354326984167099D+07,
     *   -0.242455710320062637D+08,
     *    0.108839673804066181D+08,
     *   -0.149206813160777092D+08,
     *   -0.151579723783016205D+07,
     *   -0.353630466480416059D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_v1_2=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Ar-OH(X) PES averaged over v=1 of OH
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_v1_3(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=35)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .       3.5000D0,
     .       3.7500D0,
     .       4.0000D0,
     .       4.2500D0,
     .       4.5000D0,
     .       4.7500D0,
     .       5.0000D0,
     .       5.2500D0,
     .       5.5000D0,
     .       5.7500D0,
     .       6.0000D0,
     .       6.2500D0,
     .       6.5000D0,
     .       6.7500D0,
     .       7.0000D0,
     .       7.2500D0,
     .       7.5000D0,
     .       7.7500D0,
     .       8.0000D0,
     .       8.2500D0,
     .       8.5000D0,
     .       8.7500D0,
     .       9.0000D0,
     .       9.2500D0,
     .       9.5000D0,
     .       9.7500D0,
     .      10.0000D0,
     .      11.0000D0,
     .      12.0000D0,
     .      13.0000D0,
     .      14.0000D0,
     .      16.0000D0,
     .      18.0000D0,
     .      20.0000D0,
     .      25.0000D0/
      DATA VCOEF/
     *    0.859192554475184250D+10,
     *   -0.125609422863674507D+11,
     *    0.485285662895839500D+10,
     *   -0.142725788113523340D+10,
     *    0.675275668764965653D+09,
     *   -0.159896123082239896D+09,
     *    0.117128172011925727D+08,
     *   -0.122351120107458651D+09,
     *   -0.109787835854080349D+09,
     *   -0.111999096574869514D+09,
     *   -0.902431887137089968D+08,
     *   -0.671145275146352053D+08,
     *   -0.380717732748788893D+08,
     *   -0.667208293719951808D+07,
     *    0.215669835093981028D+08,
     *    0.372439284504075646D+08,
     *    0.456626226829891205D+08,
     *    0.522630424903914481D+08,
     *    0.503791634772214442D+08,
     *    0.583762346972227395D+08,
     *    0.556177036558197737D+08,
     *    0.485622411301975250D+08,
     *    0.482893527577382326D+08,
     *    0.211983449592940807D+08,
     *    0.377898507066907883D+08,
     *   -0.278066608520050049D+08,
     *    0.131346341599422932D+09,
     *    0.256655076251201630D+08,
     *   -0.419634147122349739D+08,
     *    0.368849311817879677D+08,
     *   -0.370363812506012917D+08,
     *    0.929989744292688370D+07,
     *   -0.230516313281011581D+07,
     *   -0.621479897002601624D+07,
     *   -0.149550995318967104D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_v1_3=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Ar-OH(X) PES averaged over v=1 of OH
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_v1_4(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=35)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .       3.5000D0,
     .       3.7500D0,
     .       4.0000D0,
     .       4.2500D0,
     .       4.5000D0,
     .       4.7500D0,
     .       5.0000D0,
     .       5.2500D0,
     .       5.5000D0,
     .       5.7500D0,
     .       6.0000D0,
     .       6.2500D0,
     .       6.5000D0,
     .       6.7500D0,
     .       7.0000D0,
     .       7.2500D0,
     .       7.5000D0,
     .       7.7500D0,
     .       8.0000D0,
     .       8.2500D0,
     .       8.5000D0,
     .       8.7500D0,
     .       9.0000D0,
     .       9.2500D0,
     .       9.5000D0,
     .       9.7500D0,
     .      10.0000D0,
     .      11.0000D0,
     .      12.0000D0,
     .      13.0000D0,
     .      14.0000D0,
     .      16.0000D0,
     .      18.0000D0,
     .      20.0000D0,
     .      25.0000D0/
      DATA VCOEF/
     *    0.827060704770217323D+10,
     *   -0.119769956121775742D+11,
     *    0.395470477972368670D+10,
     *   -0.134208402107469010D+10,
     *    0.694102913803149462D+09,
     *   -0.119050998100738972D+08,
     *    0.151272395408098072D+09,
     *    0.216229679495098516D+08,
     *    0.138876868693691939D+08,
     *   -0.827191774428579211D+07,
     *   -0.870397642370319366D+07,
     *   -0.999420625145703554D+07,
     *   -0.466337343880242109D+07,
     *    0.494562081041175127D+07,
     *    0.189954859161335230D+08,
     *    0.218555250190041065D+08,
     *    0.254796198460209370D+08,
     *    0.172770967333560586D+08,
     *    0.184335809171520472D+08,
     *    0.211576863012557030D+08,
     *    0.226754009881662130D+08,
     *    0.259947388990315199D+08,
     *    0.258114345444073677D+08,
     *    0.976693898658800125D+07,
     *    0.178316603936629295D+08,
     *   -0.292996272630167007D+08,
     *    0.654209042357521057D+08,
     *    0.209651759095721245D+08,
     *   -0.349817530532126427D+08,
     *    0.368634074738025665D+08,
     *   -0.347212968522057533D+08,
     *    0.742028971953582764D+07,
     *   -0.127989297456741333D+05,
     *   -0.316997991765940189D+07,
     *   -0.204940882508063316D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_v1_4=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Ar-OH(X) PES averaged over v=1 of OH
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_v1_5(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=35)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .       3.5000D0,
     .       3.7500D0,
     .       4.0000D0,
     .       4.2500D0,
     .       4.5000D0,
     .       4.7500D0,
     .       5.0000D0,
     .       5.2500D0,
     .       5.5000D0,
     .       5.7500D0,
     .       6.0000D0,
     .       6.2500D0,
     .       6.5000D0,
     .       6.7500D0,
     .       7.0000D0,
     .       7.2500D0,
     .       7.5000D0,
     .       7.7500D0,
     .       8.0000D0,
     .       8.2500D0,
     .       8.5000D0,
     .       8.7500D0,
     .       9.0000D0,
     .       9.2500D0,
     .       9.5000D0,
     .       9.7500D0,
     .      10.0000D0,
     .      11.0000D0,
     .      12.0000D0,
     .      13.0000D0,
     .      14.0000D0,
     .      16.0000D0,
     .      18.0000D0,
     .      20.0000D0,
     .      25.0000D0/
      DATA VCOEF/
     *    0.673037022419480705D+10,
     *   -0.104370710558217239D+11,
     *    0.395932903222782135D+10,
     *   -0.121055890276764584D+10,
     *    0.611286371384364843D+09,
     *   -0.286040903918310925D+08,
     *    0.142133570830656528D+09,
     *    0.355363070605719611D+08,
     *    0.357210549455145374D+08,
     *    0.169739563732413054D+08,
     *    0.147919810231767893D+08,
     *    0.819475811795598269D+07,
     *    0.553300978581219912D+07,
     *    0.544109598136220872D+07,
     *    0.130002798645498455D+08,
     *    0.147781671192125678D+08,
     *    0.154992881973854303D+08,
     *    0.718304814085435867D+07,
     *    0.881007782097814232D+07,
     *   -0.676924356623291969D+05,
     *    0.836652145501983166D+07,
     *    0.151895305461496115D+08,
     *    0.103820771577725410D+08,
     *    0.340151155089330673D+07,
     *    0.188280646259493828D+08,
     *   -0.434892747954919338D+08,
     *    0.383632137878322601D+08,
     *    0.259321910313081741D+08,
     *   -0.353506375950245857D+08,
     *    0.270293555188598633D+08,
     *   -0.193188648770427704D+08,
     *    0.108650854430484772D+08,
     *   -0.184290113350756168D+08,
     *    0.152589071692034006D+08,
     *   -0.688156147998526692D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_v1_5=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Ar-OH(X) PES averaged over v=1 of OH
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_v1_6(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=35)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .       3.5000D0,
     .       3.7500D0,
     .       4.0000D0,
     .       4.2500D0,
     .       4.5000D0,
     .       4.7500D0,
     .       5.0000D0,
     .       5.2500D0,
     .       5.5000D0,
     .       5.7500D0,
     .       6.0000D0,
     .       6.2500D0,
     .       6.5000D0,
     .       6.7500D0,
     .       7.0000D0,
     .       7.2500D0,
     .       7.5000D0,
     .       7.7500D0,
     .       8.0000D0,
     .       8.2500D0,
     .       8.5000D0,
     .       8.7500D0,
     .       9.0000D0,
     .       9.2500D0,
     .       9.5000D0,
     .       9.7500D0,
     .      10.0000D0,
     .      11.0000D0,
     .      12.0000D0,
     .      13.0000D0,
     .      14.0000D0,
     .      16.0000D0,
     .      18.0000D0,
     .      20.0000D0,
     .      25.0000D0/
      DATA VCOEF/
     *    0.692035016093603134D+10,
     *   -0.109727630511868362D+11,
     *    0.407612957493325949D+10,
     *   -0.110862044851891136D+10,
     *    0.684447455905812025D+09,
     *    0.745860028305446357D+07,
     *    0.167353449591341436D+09,
     *    0.528648712054972947D+08,
     *    0.478656035369950309D+08,
     *    0.261296042910535485D+08,
     *    0.213783671980904341D+08,
     *    0.153720175573219210D+08,
     *    0.683268162119179498D+07,
     *    0.648272817802405357D+07,
     *    0.482987243984264135D+07,
     *    0.860769370375424623D+07,
     *    0.712119476335406303D+07,
     *    0.677976677635732293D+07,
     *    0.476031739449297264D+07,
     *   -0.673426524572861195D+07,
     *    0.902105489477342367D+07,
     *    0.460276881071269512D+07,
     *    0.501938429079115391D+07,
     *    0.163454113066759109D+08,
     *   -0.426871518761610985D+07,
     *   -0.720205258449649811D+07,
     *   -0.912340902913689613D+07,
     *    0.221441716880474091D+08,
     *   -0.186448728385405540D+08,
     *    0.134573548104190826D+08,
     *   -0.995780165848541260D+07,
     *    0.453404942611432076D+07,
     *   -0.522227722381687164D+07,
     *    0.395108042498874664D+07,
     *   -0.128741123310571909D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_v1_6=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Ar-OH(X) PES averaged over v=1 of OH
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_v1_7(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=35)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .       3.5000D0,
     .       3.7500D0,
     .       4.0000D0,
     .       4.2500D0,
     .       4.5000D0,
     .       4.7500D0,
     .       5.0000D0,
     .       5.2500D0,
     .       5.5000D0,
     .       5.7500D0,
     .       6.0000D0,
     .       6.2500D0,
     .       6.5000D0,
     .       6.7500D0,
     .       7.0000D0,
     .       7.2500D0,
     .       7.5000D0,
     .       7.7500D0,
     .       8.0000D0,
     .       8.2500D0,
     .       8.5000D0,
     .       8.7500D0,
     .       9.0000D0,
     .       9.2500D0,
     .       9.5000D0,
     .       9.7500D0,
     .      10.0000D0,
     .      11.0000D0,
     .      12.0000D0,
     .      13.0000D0,
     .      14.0000D0,
     .      16.0000D0,
     .      18.0000D0,
     .      20.0000D0,
     .      25.0000D0/
      DATA VCOEF/
     *    0.499976737055067825D+10,
     *   -0.785685272226776218D+10,
     *    0.285803244322528934D+10,
     *   -0.823625497227981806D+09,
     *    0.482908492749451160D+09,
     *    0.195256577118816301D+08,
     *    0.137838921107580125D+09,
     *    0.511471263072611466D+08,
     *    0.435459973673850298D+08,
     *    0.217386660318475217D+08,
     *    0.175862239092457294D+08,
     *    0.128141996912072152D+08,
     *    0.716341318713834975D+07,
     *    0.460019955717620254D+07,
     *    0.289555340943744779D+07,
     *    0.334669101450830698D+07,
     *    0.696939573441737890D+07,
     *    0.426919419272516668D+07,
     *    0.297418479805972427D+07,
     *   -0.776662999972969294D+06,
     *   -0.279505194629400969D+07,
     *    0.488311379208642244D+07,
     *   -0.513567759295916557D+07,
     *    0.163399677751637697D+08,
     *   -0.121381308076456785D+08,
     *    0.159442330775008202D+08,
     *   -0.184272574036288261D+08,
     *    0.121868667101514339D+08,
     *   -0.946912604580903053D+07,
     *    0.367599686305856705D+07,
     *   -0.201664034753322601D+06,
     *   -0.124884855790472031D+07,
     *    0.190673325796842575D+06,
     *    0.743971553496778011D+06,
     *   -0.469869782510846853D+06/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_v1_7=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Ar-OH(X) PES averaged over v=1 of OH
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_v1_8(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=35)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .       3.5000D0,
     .       3.7500D0,
     .       4.0000D0,
     .       4.2500D0,
     .       4.5000D0,
     .       4.7500D0,
     .       5.0000D0,
     .       5.2500D0,
     .       5.5000D0,
     .       5.7500D0,
     .       6.0000D0,
     .       6.2500D0,
     .       6.5000D0,
     .       6.7500D0,
     .       7.0000D0,
     .       7.2500D0,
     .       7.5000D0,
     .       7.7500D0,
     .       8.0000D0,
     .       8.2500D0,
     .       8.5000D0,
     .       8.7500D0,
     .       9.0000D0,
     .       9.2500D0,
     .       9.5000D0,
     .       9.7500D0,
     .      10.0000D0,
     .      11.0000D0,
     .      12.0000D0,
     .      13.0000D0,
     .      14.0000D0,
     .      16.0000D0,
     .      18.0000D0,
     .      20.0000D0,
     .      25.0000D0/
      DATA VCOEF/
     *    0.488489210328006363D+10,
     *   -0.830621337003962708D+10,
     *    0.353475920771725273D+10,
     *   -0.841574881278055429D+09,
     *    0.513568890016118884D+09,
     *   -0.219607552564159594D+08,
     *    0.116157909031362623D+09,
     *    0.309948858047876321D+08,
     *    0.335875657986976281D+08,
     *    0.154071312678235844D+08,
     *    0.136104990342693031D+08,
     *    0.829106003845654428D+07,
     *    0.704638960663027689D+07,
     *   -0.432943269062012434D+06,
     *    0.140879310586333275D+07,
     *    0.262805596688133478D+07,
     *   -0.469720233107656240D+06,
     *    0.689708851751047373D+07,
     *    0.213022842749938369D+07,
     *   -0.240278921629370749D+07,
     *    0.706009447846928239D+07,
     *   -0.586440131135314703D+07,
     *    0.921337751780390739D+06,
     *   -0.790394030342519283D+07,
     *    0.156486458822079897D+08,
     *   -0.499030632926702499D+07,
     *   -0.566669525288105011D+07,
     *    0.243972150150561333D+07,
     *    0.474531014162063599D+07,
     *   -0.102721909754810333D+08,
     *    0.635537694496965408D+07,
     *    0.176271032631492615D+07,
     *   -0.587085899862825871D+07,
     *    0.380888036167943478D+07,
     *   -0.568295308615565300D+06/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_v1_8=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Ar-OH(X) PES averaged over v=1 of OH
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_v1_9(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=35)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .       3.5000D0,
     .       3.7500D0,
     .       4.0000D0,
     .       4.2500D0,
     .       4.5000D0,
     .       4.7500D0,
     .       5.0000D0,
     .       5.2500D0,
     .       5.5000D0,
     .       5.7500D0,
     .       6.0000D0,
     .       6.2500D0,
     .       6.5000D0,
     .       6.7500D0,
     .       7.0000D0,
     .       7.2500D0,
     .       7.5000D0,
     .       7.7500D0,
     .       8.0000D0,
     .       8.2500D0,
     .       8.5000D0,
     .       8.7500D0,
     .       9.0000D0,
     .       9.2500D0,
     .       9.5000D0,
     .       9.7500D0,
     .      10.0000D0,
     .      11.0000D0,
     .      12.0000D0,
     .      13.0000D0,
     .      14.0000D0,
     .      16.0000D0,
     .      18.0000D0,
     .      20.0000D0,
     .      25.0000D0/
      DATA VCOEF/
     *    0.246497077365119648D+10,
     *   -0.391084677382457018D+10,
     *    0.139753282659399629D+10,
     *   -0.362642941624485731D+09,
     *    0.244666383797357112D+09,
     *    0.164919094994843677D+08,
     *    0.690455348235098720D+08,
     *    0.249770769255332164D+08,
     *    0.207969856705390066D+08,
     *    0.104250566939826645D+08,
     *    0.812327438830447197D+07,
     *    0.639834807990576327D+07,
     *    0.406441150533746742D+07,
     *    0.168480549799098074D+07,
     *   -0.474282882602199912D+06,
     *    0.673207717631280422D+06,
     *    0.854392132356479764D+06,
     *    0.226026724373819679D+07,
     *    0.180913822251722217D+06,
     *    0.401005345576913655D+07,
     *   -0.128488733884060383D+07,
     *   -0.200506082244688272D+07,
     *    0.380500601768904924D+07,
     *   -0.150644024410961866D+08,
     *    0.174665940317852497D+08,
     *   -0.644627301660919189D+07,
     *    0.449468370433330536D+06,
     *   -0.740559105869412422D+06,
     *    0.122583123555755615D+07,
     *   -0.666987386851310730D+06,
     *    0.737421438435554504D+06,
     *   -0.417653487955093384D+07,
     *    0.101223066502622366D+08,
     *   -0.902588704053232074D+07,
     *    0.282681829763810337D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_v1_9=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Ar-OH(X) PES averaged over v=1 of OH
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL2_v1_10(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=35)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .       3.5000D0,
     .       3.7500D0,
     .       4.0000D0,
     .       4.2500D0,
     .       4.5000D0,
     .       4.7500D0,
     .       5.0000D0,
     .       5.2500D0,
     .       5.5000D0,
     .       5.7500D0,
     .       6.0000D0,
     .       6.2500D0,
     .       6.5000D0,
     .       6.7500D0,
     .       7.0000D0,
     .       7.2500D0,
     .       7.5000D0,
     .       7.7500D0,
     .       8.0000D0,
     .       8.2500D0,
     .       8.5000D0,
     .       8.7500D0,
     .       9.0000D0,
     .       9.2500D0,
     .       9.5000D0,
     .       9.7500D0,
     .      10.0000D0,
     .      11.0000D0,
     .      12.0000D0,
     .      13.0000D0,
     .      14.0000D0,
     .      16.0000D0,
     .      18.0000D0,
     .      20.0000D0,
     .      25.0000D0/
      DATA VCOEF/
     *    0.635378309750134468D+09,
     *   -0.126097034383236670D+10,
     *    0.676062502210063100D+09,
     *   -0.114397074728420645D+09,
     *    0.734572981431997865D+08,
     *   -0.161409743197609112D+08,
     *    0.727666744314333517D+07,
     *   -0.157940157589579560D+07,
     *    0.783483200075506466D+06,
     *    0.796020311101283878D+06,
     *    0.783356283861771226D+04,
     *   -0.212647929998196661D+06,
     *   -0.341908220775331487D+06,
     *   -0.169273988468282670D+07,
     *    0.267104070800969750D+07,
     *   -0.279388865339043736D+07,
     *    0.304412250920380652D+07,
     *   -0.186138569601456076D+07,
     *    0.167310666243232042D+07,
     *   -0.283668909895319119D+07,
     *    0.165787936427389598D+07,
     *    0.255586302646353841D+07,
     *   -0.847706443634274602D+07,
     *    0.168696053560767621D+08,
     *   -0.206798828112342507D+08,
     *    0.150405100702401698D+08,
     *   -0.650302241154339910D+07,
     *    0.327344537652409077D+07,
     *   -0.515177378635001183D+07,
     *    0.664671502199876308D+07,
     *   -0.573750262959885597D+07,
     *    0.458614993874284625D+07,
     *   -0.372581125872316957D+07,
     *    0.190855824749347568D+07,
     *   -0.131233984172215685D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_v1_10=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Ar-OH(X) PES averaged over v=1 of OH
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL2_v1_2(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=35)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .       3.5000D0,
     .       3.7500D0,
     .       4.0000D0,
     .       4.2500D0,
     .       4.5000D0,
     .       4.7500D0,
     .       5.0000D0,
     .       5.2500D0,
     .       5.5000D0,
     .       5.7500D0,
     .       6.0000D0,
     .       6.2500D0,
     .       6.5000D0,
     .       6.7500D0,
     .       7.0000D0,
     .       7.2500D0,
     .       7.5000D0,
     .       7.7500D0,
     .       8.0000D0,
     .       8.2500D0,
     .       8.5000D0,
     .       8.7500D0,
     .       9.0000D0,
     .       9.2500D0,
     .       9.5000D0,
     .       9.7500D0,
     .      10.0000D0,
     .      11.0000D0,
     .      12.0000D0,
     .      13.0000D0,
     .      14.0000D0,
     .      16.0000D0,
     .      18.0000D0,
     .      20.0000D0,
     .      25.0000D0/
      DATA VCOEF/
     *    0.551562021546957016D+10,
     *   -0.613824728459838009D+10,
     *    0.178447426498553967D+10,
     *   -0.891456156499824286D+09,
     *    0.105571055364457279D+09,
     *   -0.213809198376985490D+09,
     *   -0.999763196322844177D+08,
     *   -0.137540345748657167D+09,
     *   -0.120477609061157703D+09,
     *   -0.114736093227298677D+09,
     *   -0.972327972597948611D+08,
     *   -0.794784295772194266D+08,
     *   -0.570306527247428000D+08,
     *   -0.364744780502110645D+08,
     *   -0.136872004963729940D+08,
     *    0.546936048281349242D+07,
     *    0.223143892697739676D+08,
     *    0.347789682111487389D+08,
     *    0.414534727845090628D+08,
     *    0.482183225255981088D+08,
     *    0.461176439503142834D+08,
     *    0.464021553735514879D+08,
     *    0.515470522718217373D+08,
     *    0.372778955298789740D+08,
     *    0.654120442537258863D+08,
     *   -0.106442460678696632D+07,
     *    0.125516750110187054D+09,
     *    0.755387204071583748D+08,
     *    0.639895609346985817D+07,
     *    0.236507452274203300D+07,
     *   -0.588307587450456619D+07,
     *   -0.173211910118567944D+07,
     *   -0.756718095421719551D+07,
     *    0.345694234806984663D+07,
     *   -0.385454541798073053D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_v1_2=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Ar-OH(X) PES averaged over v=1 of OH
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL2_v1_3(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=35)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .       3.5000D0,
     .       3.7500D0,
     .       4.0000D0,
     .       4.2500D0,
     .       4.5000D0,
     .       4.7500D0,
     .       5.0000D0,
     .       5.2500D0,
     .       5.5000D0,
     .       5.7500D0,
     .       6.0000D0,
     .       6.2500D0,
     .       6.5000D0,
     .       6.7500D0,
     .       7.0000D0,
     .       7.2500D0,
     .       7.5000D0,
     .       7.7500D0,
     .       8.0000D0,
     .       8.2500D0,
     .       8.5000D0,
     .       8.7500D0,
     .       9.0000D0,
     .       9.2500D0,
     .       9.5000D0,
     .       9.7500D0,
     .      10.0000D0,
     .      11.0000D0,
     .      12.0000D0,
     .      13.0000D0,
     .      14.0000D0,
     .      16.0000D0,
     .      18.0000D0,
     .      20.0000D0,
     .      25.0000D0/
      DATA VCOEF/
     *   -0.105180609565459251D+10,
     *    0.134226699039712906D+10,
     *   -0.308118047680962682D+09,
     *    0.170845238497830629D+09,
     *   -0.447638823194369972D+08,
     *    0.163596255175952911D+08,
     *   -0.624161552646588720D+07,
     *    0.891793263677228242D+05,
     *   -0.528315470756258629D+07,
     *   -0.983411843942059204D+07,
     *   -0.142023756878590845D+08,
     *   -0.176417594504508600D+08,
     *   -0.171137189698466472D+08,
     *   -0.175248185381578058D+08,
     *   -0.139984349510112405D+08,
     *   -0.155138777271148562D+08,
     *   -0.801300700098683871D+07,
     *   -0.103103988147014230D+08,
     *   -0.536062918855826557D+07,
     *   -0.649101622110903263D+06,
     *   -0.295653563189901412D+07,
     *    0.473934445518890023D+07,
     *   -0.430266475316023827D+07,
     *    0.869329316362112761D+07,
     *   -0.124561440345205069D+08,
     *    0.103149831427003443D+08,
     *    0.989525925547480583D+05,
     *    0.102948965141957402D+08,
     *    0.146257383129262924D+07,
     *    0.598732068878573179D+07,
     *   -0.801081026164627075D+07,
     *    0.771003576858651638D+07,
     *   -0.895080356262741983D+07,
     *    0.637089204081906378D+07,
     *   -0.370312625229605287D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_v1_3=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Ar-OH(X) PES averaged over v=1 of OH
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL2_v1_4(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=35)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .       3.5000D0,
     .       3.7500D0,
     .       4.0000D0,
     .       4.2500D0,
     .       4.5000D0,
     .       4.7500D0,
     .       5.0000D0,
     .       5.2500D0,
     .       5.5000D0,
     .       5.7500D0,
     .       6.0000D0,
     .       6.2500D0,
     .       6.5000D0,
     .       6.7500D0,
     .       7.0000D0,
     .       7.2500D0,
     .       7.5000D0,
     .       7.7500D0,
     .       8.0000D0,
     .       8.2500D0,
     .       8.5000D0,
     .       8.7500D0,
     .       9.0000D0,
     .       9.2500D0,
     .       9.5000D0,
     .       9.7500D0,
     .      10.0000D0,
     .      11.0000D0,
     .      12.0000D0,
     .      13.0000D0,
     .      14.0000D0,
     .      16.0000D0,
     .      18.0000D0,
     .      20.0000D0,
     .      25.0000D0/
      DATA VCOEF/
     *    0.280296749180084419D+10,
     *   -0.423715166621422863D+10,
     *    0.130668342841033983D+10,
     *   -0.361179554832498550D+09,
     *    0.270027237421253681D+09,
     *    0.385043771239758953D+08,
     *    0.884510140557588488D+08,
     *    0.439556348915084749D+08,
     *    0.341101485052585676D+08,
     *    0.174640164731534608D+08,
     *    0.807750349776053429D+07,
     *   -0.443362451739631593D+06,
     *   -0.296555121641772892D+07,
     *   -0.712434248914509267D+07,
     *   -0.467320795593266189D+07,
     *   -0.592101949977491796D+07,
     *   -0.443084479860486090D+07,
     *   -0.156924236139996350D+07,
     *   -0.163402912105099112D+07,
     *   -0.192088107952223718D+07,
     *    0.517708758782666922D+07,
     *   -0.249093282965362072D+06,
     *    0.411346188615417480D+07,
     *   -0.113048494182419777D+07,
     *    0.164216655673587322D+07,
     *   -0.607329975595593452D+07,
     *    0.253497460660624504D+07,
     *    0.122386699235056639D+08,
     *   -0.159514692003297806D+07,
     *    0.439675625094854832D+07,
     *   -0.403740273882317543D+07,
     *    0.509764931016552448D+07,
     *   -0.879584676325440407D+07,
     *    0.691141643166953325D+07,
     *   -0.225518159516486526D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_v1_4=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Ar-OH(X) PES averaged over v=1 of OH
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL2_v1_5(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=35)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .       3.5000D0,
     .       3.7500D0,
     .       4.0000D0,
     .       4.2500D0,
     .       4.5000D0,
     .       4.7500D0,
     .       5.0000D0,
     .       5.2500D0,
     .       5.5000D0,
     .       5.7500D0,
     .       6.0000D0,
     .       6.2500D0,
     .       6.5000D0,
     .       6.7500D0,
     .       7.0000D0,
     .       7.2500D0,
     .       7.5000D0,
     .       7.7500D0,
     .       8.0000D0,
     .       8.2500D0,
     .       8.5000D0,
     .       8.7500D0,
     .       9.0000D0,
     .       9.2500D0,
     .       9.5000D0,
     .       9.7500D0,
     .      10.0000D0,
     .      11.0000D0,
     .      12.0000D0,
     .      13.0000D0,
     .      14.0000D0,
     .      16.0000D0,
     .      18.0000D0,
     .      20.0000D0,
     .      25.0000D0/
      DATA VCOEF/
     *   -0.685476243681700706D+09,
     *    0.128236026529426455D+10,
     *   -0.590380945912590504D+09,
     *    0.989527135224429965D+08,
     *   -0.102175077016646236D+09,
     *    0.142812179720055871D+07,
     *   -0.119441531280001346D+08,
     *    0.380895900518381270D+07,
     *    0.357622784239137452D+07,
     *    0.427432978621018305D+07,
     *    0.231254461487428658D+07,
     *    0.519002718261368573D+06,
     *   -0.238336767109919991D+07,
     *   -0.159075014068454504D+06,
     *   -0.356383690207342803D+07,
     *    0.746642304831810296D+06,
     *   -0.331249036186648905D+07,
     *   -0.848729352378085256D+06,
     *    0.106568276522093639D+07,
     *   -0.447544660556192603D+07,
     *    0.274459258314989135D+07,
     *    0.324465587169443816D+07,
     *   -0.989076741491079330D+06,
     *   -0.114765303144901991D+08,
     *    0.310365752907201648D+08,
     *   -0.318280714081393778D+08,
     *    0.798413908687976003D+07,
     *    0.361187246263575554D+07,
     *    0.569346881421506405D+06,
     *    0.411843227418839931D+07,
     *   -0.530919168008351326D+07,
     *   -0.113711887067079544D+07,
     *    0.108026286229470670D+08,
     *   -0.998896362503595650D+07,
     *    0.196627993770673499D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_v1_5=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Ar-OH(X) PES averaged over v=1 of OH
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL2_v1_6(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=35)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .       3.5000D0,
     .       3.7500D0,
     .       4.0000D0,
     .       4.2500D0,
     .       4.5000D0,
     .       4.7500D0,
     .       5.0000D0,
     .       5.2500D0,
     .       5.5000D0,
     .       5.7500D0,
     .       6.0000D0,
     .       6.2500D0,
     .       6.5000D0,
     .       6.7500D0,
     .       7.0000D0,
     .       7.2500D0,
     .       7.5000D0,
     .       7.7500D0,
     .       8.0000D0,
     .       8.2500D0,
     .       8.5000D0,
     .       8.7500D0,
     .       9.0000D0,
     .       9.2500D0,
     .       9.5000D0,
     .       9.7500D0,
     .      10.0000D0,
     .      11.0000D0,
     .      12.0000D0,
     .      13.0000D0,
     .      14.0000D0,
     .      16.0000D0,
     .      18.0000D0,
     .      20.0000D0,
     .      25.0000D0/
      DATA VCOEF/
     *    0.222738066091999435D+10,
     *   -0.401478595241571856D+10,
     *    0.175911275846395636D+10,
     *   -0.273265595525565147D+09,
     *    0.263999080491364300D+09,
     *   -0.134175029666525926D+08,
     *    0.385652583605211526D+08,
     *    0.297124825037958287D+07,
     *    0.729532912326268200D+07,
     *    0.265153576216162741D+07,
     *    0.233476017776091397D+07,
     *    0.104347246258773655D+07,
     *   -0.818465491453744471D+04,
     *   -0.230230210275948048D+05,
     *   -0.191199606419456005D+07,
     *    0.554934809124469757D+04,
     *   -0.615150888803780079D+06,
     *   -0.158668485400207341D+07,
     *    0.432995348635965586D+07,
     *   -0.901018568258102611D+07,
     *    0.617189692412494123D+07,
     *   -0.424780843507787585D+07,
     *   -0.480761733083820343D+07,
     *    0.614307969862449169D+07,
     *    0.785453659217715263D+06,
     *    0.240645136721813679D+07,
     *   -0.154900284221386909D+07,
     *   -0.449390252791643143D+06,
     *   -0.756720911410331726D+06,
     *    0.208633795571231842D+07,
     *   -0.667303069689273834D+06,
     *   -0.147945605526924133D+07,
     *    0.359337010912597179D+07,
     *   -0.366914324782004952D+07,
     *    0.206088958025124669D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_v1_6=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Ar-OH(X) PES averaged over v=1 of OH
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL2_v1_7(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=35)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .       3.5000D0,
     .       3.7500D0,
     .       4.0000D0,
     .       4.2500D0,
     .       4.5000D0,
     .       4.7500D0,
     .       5.0000D0,
     .       5.2500D0,
     .       5.5000D0,
     .       5.7500D0,
     .       6.0000D0,
     .       6.2500D0,
     .       6.5000D0,
     .       6.7500D0,
     .       7.0000D0,
     .       7.2500D0,
     .       7.5000D0,
     .       7.7500D0,
     .       8.0000D0,
     .       8.2500D0,
     .       8.5000D0,
     .       8.7500D0,
     .       9.0000D0,
     .       9.2500D0,
     .       9.5000D0,
     .       9.7500D0,
     .      10.0000D0,
     .      11.0000D0,
     .      12.0000D0,
     .      13.0000D0,
     .      14.0000D0,
     .      16.0000D0,
     .      18.0000D0,
     .      20.0000D0,
     .      25.0000D0/
      DATA VCOEF/
     *   -0.437126404743376553D+09,
     *    0.102956969565840364D+10,
     *   -0.745319121515727758D+09,
     *    0.188647713477123797D+09,
     *   -0.678404078863447756D+08,
     *    0.301320236318493113D+08,
     *   -0.459083904672794417D+07,
     *    0.455337607124263421D+07,
     *   -0.147519700925404206D+06,
     *    0.182448412770959176D+07,
     *    0.105136011600062624D+07,
     *    0.183053870134465396D+07,
     *    0.142972125082314014D+06,
     *    0.273929436146281660D+05,
     *   -0.634408579372234643D+06,
     *   -0.232383690879891813D+07,
     *   -0.823256846051402390D+06,
     *    0.428534111442386359D+07,
     *   -0.484436155225405842D+07,
     *    0.292445001185289025D+07,
     *    0.175215894486885518D+06,
     *   -0.225851722531075217D+07,
     *   -0.182601237367584184D+07,
     *    0.382292995975905657D+07,
     *   -0.922426884500966966D+07,
     *    0.933861848416770995D+07,
     *   -0.156813947954565287D+06,
     *   -0.384026118923676014D+07,
     *    0.850831200910717249D+07,
     *   -0.131869604516389966D+08,
     *    0.991461954176604748D+07,
     *   -0.511732093551319838D+07,
     *    0.614241444741892815D+07,
     *   -0.508887276001676917D+07,
     *    0.181778024060413986D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_v1_7=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Ar-OH(X) PES averaged over v=1 of OH
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL2_v1_8(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=35)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .       3.5000D0,
     .       3.7500D0,
     .       4.0000D0,
     .       4.2500D0,
     .       4.5000D0,
     .       4.7500D0,
     .       5.0000D0,
     .       5.2500D0,
     .       5.5000D0,
     .       5.7500D0,
     .       6.0000D0,
     .       6.2500D0,
     .       6.5000D0,
     .       6.7500D0,
     .       7.0000D0,
     .       7.2500D0,
     .       7.5000D0,
     .       7.7500D0,
     .       8.0000D0,
     .       8.2500D0,
     .       8.5000D0,
     .       8.7500D0,
     .       9.0000D0,
     .       9.2500D0,
     .       9.5000D0,
     .       9.7500D0,
     .      10.0000D0,
     .      11.0000D0,
     .      12.0000D0,
     .      13.0000D0,
     .      14.0000D0,
     .      16.0000D0,
     .      18.0000D0,
     .      20.0000D0,
     .      25.0000D0/
      DATA VCOEF/
     *    0.171496552512664461D+10,
     *   -0.332403650476360559D+10,
     *    0.169613864979622293D+10,
     *   -0.273277733691165328D+09,
     *    0.199853590832758695D+09,
     *   -0.355187432365255058D+08,
     *    0.232011704451312646D+08,
     *   -0.403865689515417255D+07,
     *    0.307425737129962258D+07,
     *   -0.113642258112376556D+07,
     *    0.745473259140715003D+06,
     *    0.124948015924632549D+06,
     *    0.687178587664143182D+06,
     *    0.155657431936047971D+07,
     *   -0.180175784386733174D+07,
     *    0.116915476627191901D+07,
     *   -0.214665737533244491D+07,
     *    0.151606867232201993D+07,
     *   -0.516908199331867695D+07,
     *    0.864544615260211751D+07,
     *   -0.865243910440486670D+07,
     *    0.653088007795675099D+07,
     *   -0.327527266523277760D+07,
     *    0.709029137343972921D+06,
     *    0.168134389358377457D+07,
     *   -0.540442137274205685D+07,
     *    0.453438314145314693D+07,
     *   -0.394835648525118828D+06,
     *   -0.207120828942644596D+07,
     *    0.366828105682921410D+07,
     *   -0.204337808928084373D+07,
     *    0.340709457450389862D+06,
     *   -0.115116775850397348D+07,
     *    0.123772740782240033D+07,
     *   -0.293368428844437003D+06/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_v1_8=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Ar-OH(X) PES averaged over v=1 of OH
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL2_v1_9(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPT=35)
      DIMENSION VAB(NPT),VCOEF(NPT)
      DIMENSION R(NPT),THETA(NPT)
      DATA R/
     .       3.5000D0,
     .       3.7500D0,
     .       4.0000D0,
     .       4.2500D0,
     .       4.5000D0,
     .       4.7500D0,
     .       5.0000D0,
     .       5.2500D0,
     .       5.5000D0,
     .       5.7500D0,
     .       6.0000D0,
     .       6.2500D0,
     .       6.5000D0,
     .       6.7500D0,
     .       7.0000D0,
     .       7.2500D0,
     .       7.5000D0,
     .       7.7500D0,
     .       8.0000D0,
     .       8.2500D0,
     .       8.5000D0,
     .       8.7500D0,
     .       9.0000D0,
     .       9.2500D0,
     .       9.5000D0,
     .       9.7500D0,
     .      10.0000D0,
     .      11.0000D0,
     .      12.0000D0,
     .      13.0000D0,
     .      14.0000D0,
     .      16.0000D0,
     .      18.0000D0,
     .      20.0000D0,
     .      25.0000D0/
      DATA VCOEF/
     *   -0.191707191425623775D+09,
     *    0.512698867559487104D+09,
     *   -0.444021409676480532D+09,
     *    0.137136571594292730D+09,
     *   -0.333951673903162107D+08,
     *    0.196620350327256992D+08,
     *   -0.337827765082944138D+07,
     *    0.306982345831889799D+07,
     *   -0.241579473882267252D+06,
     *    0.105203191038373392D+07,
     *   -0.907743877035256475D+06,
     *    0.508129328343223780D+06,
     *   -0.103831654301066324D+07,
     *    0.316073559016839601D+06,
     *    0.956809096646308899D+04,
     *    0.181941479224451259D+07,
     *   -0.122284751581811905D+07,
     *    0.921827676135409623D+06,
     *   -0.361870435070526227D+07,
     *    0.371338914228350297D+07,
     *   -0.180831497700294107D+07,
     *   -0.806711127732962370D+06,
     *    0.157046628262043186D+07,
     *    0.201108689703261480D+06,
     *    0.101436686084882915D+07,
     *   -0.505873305366268009D+07,
     *    0.454906326777891070D+07,
     *   -0.313907645809917152D+07,
     *    0.743117002184826136D+07,
     *   -0.137959044164623022D+08,
     *    0.126723586843979955D+08,
     *   -0.598019523046575487D+07,
     *    0.167134975469499826D+06,
     *    0.278295169733827561D+07,
     *   -0.169003616804941371D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_v1_9=SUMA
       RETURN
       END

C------RKHS ROUTINES

C------------------------RKHS---------------------------------
C Fortran  Deck with distance-like and angular-like 
C kernels for RKHS interpolation method
C****************************************************
C* By J. KLOS  12/05/2007   University of Maryland  *
C* REVISON: 1.0                                     * 
C* jasztol@gmail.com                                *
C****************************************************
C-------------------------------------------------------------
C FUNCTIONS:
C
C RKHS_DK(X,Y,N,M): Reproducing kernel for  distance-like (DK) 
C                   variables, see Ho, Rabitz Eq(17) of
C                   J. Chem. Phys. 104, 2584 (1996) 
C                   Range = [0, inf]
C
C  INPUT: X,Y-DOUBLE PRECISION
C          N:INTEGER:ORDER OF SMOOTHNESS OF THE FUNCTION
C          M>=0:INTEGER:RECIPROCAL POWER OF THE WEIGHT W(X)=X**(-M)
C  OUTPUT: DOUBLE PRECISION
C -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -            
C RKHS_AK(X,Y,N): Reproducing kernel for  angular-like (AK) 
C                   variables, see Ho, Rabitz Eq(23) of
C                   J. Chem. Phys. 104, 2584 (1996) 
C                    Range = [0, 1]
C    
C  INPUT: X,Y-DOUBLE PRECISION
C          N:INTEGER:ORDER OF SMOOTHNESS OF THE FUNCTION
C  OUPUT: DOUBLE PRECISION
C  HINT: SCALE ANGLES TO [0,1] INTEGRAL FOR ANGLES, FOR EXAMPLE: 
C          X=(1.d0-cosd(ANGLE))/2.d0
C -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
C RKHS_OP(X,Y,N,M):Reproducing kernel for angular like 
C                    variables with a use of orthogonal 
C                    polynomials Plm(L,M,X) where X is 
C                    K=SUM_L=M^N PLM(X)PLM(X')
C                    variable transformed from angles as
C                    x=cos(theta)
C-------------------------------------------------------------
       DOUBLE PRECISION FUNCTION RKHS_DKN2(X,Y,M)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IF (X .LT. Y) THEN
           XB=X
           XU=Y
       ELSE
           XB=Y
           XU=X
       ENDIF

C       M=ABS(M)! FORCE M>=0
       NSQR=4
       XUPOWM=XU**(-(M+1))

C      TERM WITH BETA FUNCTION
C       CALL BETA(DFLOAT(M+1),DFLOAT(N),BETATERM)
C      N=2 case direct expression
       BETATERM=1.D0/DFLOAT((M+1)*(M+2))

C      GET TERM WITH 2F1 HYPERGEOMETRIC  GAUSS FUNCTION
C       CALL HYGFX(DFLOAT(-N+1),DFLOAT(M+1),DFLOAT(N+M+1),XB/XU,HF)
C      DIRECT EXPRESSION FOR 2F1 WITH N=2 case
       HF=1.D0-(XB/XU)*(DFLOAT(M+1)/DFLOAT(M+3))
       RKHS_DKN2=NSQR*XUPOWM*BETATERM*HF
       RETURN
       END

       DOUBLE PRECISION FUNCTION RKHS_DK(X,Y,N,M)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IF (X .LT. Y) THEN
           XB=X
           XU=Y
       ELSE
           XB=Y
           XU=X
       ENDIF

C       M=ABS(M)! FORCE M>=0
       NSQR=N*N
       XUPOWM=XU**(-(M+1))

C      TERM WITH BETA FUNCTION
       CALL BETA(DFLOAT(M+1),DFLOAT(N),BETATERM)

C      GET TERM WITH 2F1 HYPERGEOMETRIC  GAUSS FUNCTION
       CALL HYGFX(DFLOAT(-N+1),DFLOAT(M+1),DFLOAT(N+M+1),XB/XU,HF)
       RKHS_DK=NSQR*XUPOWM*BETATERM*HF
       RETURN
       END

C-------------------------------------------------------------
       DOUBLE PRECISION FUNCTION RKHS_EXPR(X,Y,XE,N)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER (ZERO=0.D0,ONE=1.D0)
       
       SUMA=ZERO
       Z1=(X-XE)/XE
       Z2=(Y-XE)/XE
       DO I=0,N
       SUMA=SUMA+(Z1**I)*(Z2**I)
       ENDDO

       RKHS_EXPR=SUMA
       RETURN
       END


C-------------------------------------------------------------
       DOUBLE PRECISION FUNCTION RKHS_AK(X,Y,N)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER (ZERO=0.D0,ONE=1.D0)
       
       IF ((X.EQ.ZERO).OR.(Y.EQ.ZERO)) THEN
          RKHS_AK=ONE
          RETURN
       ENDIF
       IF (X .LT. Y) THEN
           XB=X
           XU=Y
       ELSE
           XB=Y
           XU=X
       ENDIF

       SUMA=ZERO
C      GET POLYNOMIAL TERM
       DO I=0,N-1
       SUMA=SUMA+(XU**I)*(XB**I)
       END DO

C      GET TERM WITH 2F1 HYPERGEOMETRIC  GAUSS FUNCTION
       CALL HYGFX(DFLOAT(1),DFLOAT(-N+1),DFLOAT(N+1),XB/XU,HF)
       
       RKHS_AK=SUMA+N*(XB**N)*(XU**(N-1))*HF
       RETURN
       END

       DOUBLE PRECISION FUNCTION RKHS_AKN2(X,Y)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER (ZERO=0.D0,ONE=1.D0)
       
       IF ((X.EQ.ZERO).OR.(Y.EQ.ZERO)) THEN
          RKHS_AKN2=ONE
          RETURN
       ENDIF 
       IF (X .LT. Y) THEN
           XB=X
           XU=Y
       ELSE
           XB=Y
           XU=X
       ENDIF

C      GET POLYNOMIAL TERM
       SUMA=1.D0+XU*XB

C      GET TERM WITH 2F1 HYPERGEOMETRIC  GAUSS FUNCTION
C       CALL HYGFX(DFLOAT(1),DFLOAT(-N+1),DFLOAT(N+1),XB/XU,HF)
       HF=1.D0-(XB/XU)/3.D0
       
       RKHS_AKN2=SUMA+2.D0*(XB**2)*XU*HF
       RETURN
       END


       DOUBLE PRECISION FUNCTION RKHS_OP(X,Y,N,M)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER (ZERO=0.D0,ONE=1.D0)

       SUMA=ZERO
C      GET ORTHOGONAL POLYNOMIAL TERM
       DO I=ABS(M),N
       SUMA=SUMA+PLGNDR(I,M,X)*PLGNDR(I,M,Y)
       END DO
       RKHS_OP=SUMA
       RETURN
       END

       DOUBLE PRECISION FUNCTION RKHS_OP2L(X,Y,N,M)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER (ZERO=0.D0,ONE=1.D0)
 
       SUMA=ZERO
C      GET ORTHOGONAL POLYNOMIAL TERM
       DO I=ABS(M),N,2
       SUMA=SUMA+PLGNDR(I,M,X)*PLGNDR(I,M,Y)
       END DO
       RKHS_OP2L=SUMA
       RETURN
       END


            


      FUNCTION PLGNDR(L,M,x)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.)pause
     *'bad arguments in plgndr'
      pmm=1.
      if(m.gt.0) then
        somx2=dsqrt((1.-x)*(1.+x))
        fact=1.
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.
11      continue
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
         do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndr=pll
        endif
      endif
      return
      END



        SUBROUTINE BETA(P,Q,BT)
C
C       ==========================================
C       Purpose: Compute the beta function B(p,q)
C       Input :  p  --- Parameter  ( p > 0 )
C                q  --- Parameter  ( q > 0 )
C       Output:  BT --- B(p,q)
C       Routine called: GAMMA for computing B(x)
C       ==========================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        CALL GAMMA(P,GP)
        CALL GAMMA(Q,GQ)
        PPQ=P+Q
        CALL GAMMA(PPQ,GPQ)
        BT=GP*GQ/GPQ
        RETURN
        END


        SUBROUTINE GAMMA(X,GA)
C
C       ==================================================
C       Purpose: Compute gamma function (x)
C       Input :  x  --- Argument of ( x )
C                       ( x is not equal to 0,-1,-2,...)
C       Output:  GA --- (x)
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,
     &          -0.6558780715202538D0, -0.420026350340952D-1,
     &          0.1665386113822915D0,-.421977345555443D-1,
     &          -.96219715278770D-2, .72189432466630D-2,
     &          -.11651675918591D-2, -.2152416741149D-3,
     &          .1280502823882D-3, -.201348547807D-4,
     &          -.12504934821D-5, .11330272320D-5,
     &          -.2056338417D-6, .61160950D-8,
     &          .50020075D-8, -.11812746D-8,
     &          .1043427D-9, .77823D-11,
     &          -.36968D-11, .51D-12,
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END



        SUBROUTINE HYGFX(A,B,C,X,HF)
C
C       ====================================================
C       Purpose: Compute hypergeometric function F(a,b,c,x)
C       Input :  a --- Parameter
C                b --- Parameter
C                c --- Parameter, c <> 0,-1,-2,...
C                x --- Argument   ( x < 1 )
C       Output:  HF --- F(a,b,c,x)
C       Routines called:
C            (1) GAMMA for computing gamma function
C            (2) PSI_AAA for computing psi function
C       ====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        LOGICAL L0,L1,L2,L3,L4,L5
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        L0=C.EQ.INT(C).AND.C.LT.0.0
        L1=1.0D0-X.LT.1.0D-15.AND.C-A-B.LE.0.0
        L2=A.EQ.INT(A).AND.A.LT.0.0
        L3=B.EQ.INT(B).AND.B.LT.0.0
        L4=C-A.EQ.INT(C-A).AND.C-A.LE.0.0
        L5=C-B.EQ.INT(C-B).AND.C-B.LE.0.0
        IF (L0.OR.L1) THEN
           WRITE(*,*)'The hypergeometric series is divergent'
           RETURN
        ENDIF
        EPS=1.0D-15
        IF (X.GT.0.95) EPS=1.0D-8
        IF (X.EQ.0.0.OR.A.EQ.0.0.OR.B.EQ.0.0) THEN
           HF=1.0D0
           RETURN
        ELSE IF (1.0D0-X.EQ.EPS.AND.C-A-B.GT.0.0) THEN
           CALL GAMMA(C,GC)
           CALL GAMMA(C-A-B,GCAB)
           CALL GAMMA(C-A,GCA)
           CALL GAMMA(C-B,GCB)
           HF=GC*GCAB/(GCA*GCB)
           RETURN
        ELSE IF (1.0D0+X.LE.EPS.AND.DABS(C-A+B-1.0).LE.EPS) THEN
           G0=DSQRT(PI)*2.0D0**(-A)
           CALL GAMMA(C,G1)
           CALL GAMMA(1.0D0+A/2.0-B,G2)
           CALL GAMMA(0.5D0+0.5*A,G3)
           HF=G0*G1/(G2*G3)
           RETURN
        ELSE IF (L2.OR.L3) THEN
           IF (L2) NM=INT(ABS(A))
           IF (L3) NM=INT(ABS(B))
           HF=1.0D0
           R=1.0D0
           DO 10 K=1,NM
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X
10            HF=HF+R
           RETURN
        ELSE IF (L4.OR.L5) THEN
           IF (L4) NM=INT(ABS(C-A))
           IF (L5) NM=INT(ABS(C-B))
           HF=1.0D0
           R=1.0D0
           DO 15 K=1,NM
              R=R*(C-A+K-1.0D0)*(C-B+K-1.0D0)/(K*(C+K-1.0D0))*X
15            HF=HF+R
           HF=(1.0D0-X)**(C-A-B)*HF
           RETURN
        ENDIF
        AA=A
        BB=B
        X1=X
        IF (X.LT.0.0D0) THEN
           X=X/(X-1.0D0)
           IF (C.GT.A.AND.B.LT.A.AND.B.GT.0.0) THEN
              A=BB
              B=AA
           ENDIF
           B=C-B
        ENDIF
        IF (X.GE.0.75D0) THEN
           GM=0.0D0
           IF (DABS(C-A-B-INT(C-A-B)).LT.1.0D-15) THEN
              M=INT(C-A-B)
              CALL GAMMA(A,GA)
              CALL GAMMA(B,GB)
              CALL GAMMA(C,GC)
              CALL GAMMA(A+M,GAM)
              CALL GAMMA(B+M,GBM)
              CALL PSI_AAA(A,PA)
              CALL PSI_AAA(B,PB)
              IF (M.NE.0) GM=1.0D0
              DO 30 J=1,ABS(M)-1
30               GM=GM*J
              RM=1.0D0
              DO 35 J=1,ABS(M)
35               RM=RM*J
              F0=1.0D0
              R0=1.0D0
              R1=1.0D0
              SP0=0.D0
              SP=0.0D0
              IF (M.GE.0) THEN
                 C0=GM*GC/(GAM*GBM)
                 C1=-GC*(X-1.0D0)**M/(GA*GB*RM)
                 DO 40 K=1,M-1
                    R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(K-M))*(1.0-X)
40                  F0=F0+R0
                 DO 45 K=1,M
45                  SP0=SP0+1.0D0/(A+K-1.0)+1.0/(B+K-1.0)-1.0/K
                 F1=PA+PB+SP0+2.0D0*EL+DLOG(1.0D0-X)
                 DO 55 K=1,250
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0))
                    SM=0.0D0
                    DO 50 J=1,M
50                     SM=SM+(1.0D0-A)/((J+K)*(A+J+K-1.0))+1.0/
     &                    (B+J+K-1.0)
                    RP=PA+PB+2.0D0*EL+SP+SM+DLOG(1.0D0-X)
                    R1=R1*(A+M+K-1.0D0)*(B+M+K-1.0)/(K*(M+K))*(1.0-X)
                    F1=F1+R1*RP
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 60
55                  HW=F1
60               HF=F0*C0+F1*C1
              ELSE IF (M.LT.0) THEN
                 M=-M
                 C0=GM*GC/(GA*GB*(1.0D0-X)**M)
                 C1=-(-1)**M*GC/(GAM*GBM*RM)
                 DO 65 K=1,M-1
                    R0=R0*(A-M+K-1.0D0)*(B-M+K-1.0)/(K*(K-M))*(1.0-X)
65                  F0=F0+R0
                 DO 70 K=1,M
70                  SP0=SP0+1.0D0/K
                 F1=PA+PB-SP0+2.0D0*EL+DLOG(1.0D0-X)
                 DO 80 K=1,250
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0))
                    SM=0.0D0
                    DO 75 J=1,M
75                     SM=SM+1.0D0/(J+K)
                    RP=PA+PB+2.0D0*EL+SP-SM+DLOG(1.0D0-X)
                    R1=R1*(A+K-1.0D0)*(B+K-1.0)/(K*(M+K))*(1.0-X)
                    F1=F1+R1*RP
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 85
80                  HW=F1
85               HF=F0*C0+F1*C1
              ENDIF
           ELSE
              CALL GAMMA(A,GA)
              CALL GAMMA(B,GB)
              CALL GAMMA(C,GC)
              CALL GAMMA(C-A,GCA)
              CALL GAMMA(C-B,GCB)
              CALL GAMMA(C-A-B,GCAB)
              CALL GAMMA(A+B-C,GABC)
              C0=GC*GCAB/(GCA*GCB)
              C1=GC*GABC/(GA*GB)*(1.0D0-X)**(C-A-B)
              HF=0.0D0
              R0=C0
              R1=C1
              DO 90 K=1,250
                 R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(A+B-C+K))*(1.0-X)
                 R1=R1*(C-A+K-1.0D0)*(C-B+K-1.0)/(K*(C-A-B+K))
     &              *(1.0-X)
                 HF=HF+R0+R1
                 IF (DABS(HF-HW).LT.DABS(HF)*EPS) GO TO 95
90               HW=HF
95            HF=HF+C0+C1
           ENDIF
        ELSE
           A0=1.0D0
           IF (C.GT.A.AND.C.LT.2.0D0*A.AND.
     &         C.GT.B.AND.C.LT.2.0D0*B) THEN
              A0=(1.0D0-X)**(C-A-B)
              A=C-A
              B=C-B
           ENDIF
           HF=1.0D0
           R=1.0D0
           DO 100 K=1,250
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X
              HF=HF+R
              IF (DABS(HF-HW).LE.DABS(HF)*EPS) GO TO 105
100           HW=HF
105        HF=A0*HF
        ENDIF
        IF (X1.LT.0.0D0) THEN
           X=X1
           C0=1.0D0/(1.0D0-X)**AA
           HF=C0*HF
        ENDIF
        A=AA
        B=BB
        IF (K.GT.120) WRITE(*,115)
115     FORMAT(1X,'Warning! You should check the accuracy')
        RETURN
        END



        SUBROUTINE PSI_AAA(X,PS)
C
C       ======================================
C       Purpose: Compute Psi function
C       Input :  x  --- Argument of psi(x)
C       Output:  PS --- psi(x)
C       ======================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XA=DABS(X)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        S=0.0D0
        IF (X.EQ.INT(X).AND.X.LE.0.0) THEN
           PS=1.0D+300
           RETURN
        ELSE IF (XA.EQ.INT(XA)) THEN
           N=XA
           DO 10 K=1 ,N-1
10            S=S+1.0D0/K
           PS=-EL+S
        ELSE IF (XA+.5.EQ.INT(XA+.5)) THEN
           N=XA-.5
           DO 20 K=1,N
20            S=S+1.0/(2.0D0*K-1.0D0)
           PS=-EL+2.0D0*S-1.386294361119891D0
        ELSE
           IF (XA.LT.10.0) THEN
              N=10-INT(XA)
              DO 30 K=0,N-1
30               S=S+1.0D0/(XA+K)
              XA=XA+N
           ENDIF
           X2=1.0D0/(XA*XA)
           A1=-.8333333333333D-01
           A2=.83333333333333333D-02
           A3=-.39682539682539683D-02
           A4=.41666666666666667D-02
           A5=-.75757575757575758D-02
           A6=.21092796092796093D-01
           A7=-.83333333333333333D-01
           A8=.4432598039215686D0
           PS=DLOG(XA)-.5D0/XA+X2*(((((((A8*X2+A7)*X2+
     &        A6)*X2+A5)*X2+A4)*X2+A3)*X2+A2)*X2+A1)
           PS=PS-S
        ENDIF
        IF (X.LT.0.0) PS=PS-PI*DCOS(PI*X)/DSIN(PI*X)-1.0D0/X
        RETURN
        END
