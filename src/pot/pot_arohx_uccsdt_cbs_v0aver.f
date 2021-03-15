* System:  OH(X 2Pi)+Ar,  ab initio UCCSD(T) PES's
* with Complete Basis Set extrapolation 
* obtained from 3D PES by averaging over 
* OH(v=0) wavefunction
* Calculations by  J. Klos
*Reference: L. Scharfenberg, J. Klos, P. J. Dagdigian,
*           M. H. Alexander,G. Meijer, S. Y. T. van der Meerakker
*           Phys. Chem. Chem. Phys. 12, 10660-10670 (2010)



      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(19)
      include "common/parpot"
      potnam=' J.Klos Ar-OH(X) UCCSDT CBS <v=0> PES'
      print *, potnam
1      print *, ' r (bohr)'
      rshift=0.5
      xfact=0.8
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vsum',/,11(1pe16.8),/,
     :    '  vdif',/,9e16.8)
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
      potnam='J.Klos Ar-OH(X) UCCSDT CBS <v=0> PES'
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
      vv0=VL0_v0_0(r)*conv
      vvl(1)=VL0_v0_1(r)*conv
      vvl(2)=VL0_v0_2(r)*conv
      vvl(3)=VL0_v0_3(r)*conv
      vvl(4)=VL0_v0_4(r)*conv
      vvl(5)=VL0_v0_5(r)*conv
      vvl(6)=VL0_v0_6(r)*conv
      vvl(7)=VL0_v0_7(r)*conv
      vvl(8)=VL0_v0_8(r)*conv
      vvl(9)=VL0_v0_9(r)*conv
      vvl(10)=VL0_v0_10(r)*conv
      vvl(11)=VL2_v0_2(r)*conv
      vvl(12)=VL2_v0_3(r)*conv
      vvl(13)=VL2_v0_4(r)*conv
      vvl(14)=VL2_v0_5(r)*conv
      vvl(15)=VL2_v0_6(r)*conv
      vvl(16)=VL2_v0_7(r)*conv
      vvl(17)=VL2_v0_8(r)*conv
      vvl(18)=VL2_v0_9(r)*conv
      vvl(19)=VL2_v0_10(r)*conv
      end
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_v0_0(RR)
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
     *    0.112769332852735691D+11,
     *   -0.126647346189749107D+11,
     *    0.409916073448995924D+10,
     *   -0.202589109441241336D+10,
     *    0.206695833714786023D+08,
     *   -0.714198222289084435D+09,
     *   -0.442740401665376663D+09,
     *   -0.462404513109031856D+09,
     *   -0.347865698193598330D+09,
     *   -0.258782730228064865D+09,
     *   -0.158096944755403101D+09,
     *   -0.695077742723271102D+08,
     *    0.784579730021914747D+07,
     *    0.704698607612643689D+08,
     *    0.114332483571669981D+09,
     *    0.142940609357136428D+09,
     *    0.155488227046947300D+09,
     *    0.159338057012939453D+09,
     *    0.153750874074849248D+09,
     *    0.145148756987732917D+09,
     *    0.129305275533037633D+09,
     *    0.116953005655155569D+09,
     *    0.997719271255163997D+08,
     *    0.857610084757263660D+08,
     *    0.878855699107456356D+08,
     *    0.184011059621722177D+08,
     *    0.151624036057432383D+09,
     *    0.595207045838386193D+08,
     *   -0.642825186206332222D+07,
     *   -0.905006757790911943D+07,
     *   -0.121827843649483752D+08,
     *   -0.641676340051302128D+07,
     *   -0.124400923536501732D+08,
     *   -0.192953345586136589D+07,
     *   -0.713830947204540763D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_v0_0=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_v0_10(RR)
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
     *    0.114548370692776799D+10,
     *   -0.209921103071221566D+10,
     *    0.102723225450452948D+10,
     *   -0.228616196421086162D+09,
     *    0.131501254372664884D+09,
     *   -0.138358041085489001D+08,
     *    0.242828400872150473D+08,
     *    0.189488805947855767D+07,
     *    0.302142289737880230D+07,
     *    0.131351742733157519D+07,
     *    0.253721124455423048D+07,
     *    0.277833634461998614D+07,
     *    0.767675459626374184D+06,
     *    0.108130653562400350D+07,
     *   -0.129206036096900841D+07,
     *   -0.173963849246999552D+07,
     *    0.367631534656860493D+07,
     *   -0.692358786058444181D+06,
     *   -0.217533868820894789D+07,
     *    0.548403739819580223D+07,
     *   -0.296480939109582175D+07,
     *   -0.157529079368427116D+07,
     *    0.308937992910572561D+07,
     *   -0.716575173168069869D+07,
     *    0.969230625995718688D+07,
     *   -0.119522578195944522D+08,
     *    0.945343314237251878D+07,
     *   -0.362711037103645690D+07,
     *    0.288884275303166406D+07,
     *   -0.255451553983186884D+07,
     *    0.170176151509777573D+07,
     *   -0.725384609973672777D+06,
     *    0.613138903888004948D+06,
     *   -0.639456234643412172D+06,
     *    0.264337067871226056D+06/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_v0_10=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_v0_1(RR)
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
     *    0.623097518731232738D+10,
     *   -0.844750330984467983D+10,
     *    0.313354038237594843D+10,
     *   -0.991944534705090284D+09,
     *    0.356889510071761131D+09,
     *   -0.204953411520108134D+09,
     *   -0.864284425933138132D+08,
     *   -0.163605574707998067D+09,
     *   -0.143704922190601379D+09,
     *   -0.132394434606944144D+09,
     *   -0.102441581627078399D+09,
     *   -0.663220074658539221D+08,
     *   -0.301889414207166098D+08,
     *    0.472534087697147112D+07,
     *    0.270156549244231656D+08,
     *    0.422729832397154197D+08,
     *    0.518828951115522608D+08,
     *    0.586633874445496127D+08,
     *    0.678756006167284697D+08,
     *    0.700690344058011323D+08,
     *    0.725030386729202420D+08,
     *    0.688983615867289454D+08,
     *    0.525900510112705752D+08,
     *    0.464993982339586467D+08,
     *    0.338715753022838682D+08,
     *    0.870894875123459660D+07,
     *    0.668469093999893069D+08,
     *   -0.616221947250695061D+07,
     *   -0.100729268325794779D+07,
     *   -0.101273860472375771D+07,
     *   -0.918382440661622584D+07,
     *    0.130734011114610303D+08,
     *   -0.138485548212141469D+08,
     *   -0.377375988961985335D+07,
     *   -0.595296107912251027D+06/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_v0_1=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_v0_2(RR)
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
     *    0.727345055571242332D+10,
     *   -0.852344175421495628D+10,
     *    0.282005501419316769D+10,
     *   -0.128158525912614131D+10,
     *    0.136394002208706826D+09,
     *   -0.390279974750951529D+09,
     *   -0.233852613605053902D+09,
     *   -0.271120033372511387D+09,
     *   -0.208619852365734279D+09,
     *   -0.157936346151298761D+09,
     *   -0.955858394845910966D+08,
     *   -0.385207659033918902D+08,
     *    0.105596887667294890D+08,
     *    0.510394863277869076D+08,
     *    0.804272392264462411D+08,
     *    0.927956723526225835D+08,
     *    0.972851680765383095D+08,
     *    0.960151346646435857D+08,
     *    0.929602865039862245D+08,
     *    0.932429042352715880D+08,
     *    0.840294665039381832D+08,
     *    0.709130367549279332D+08,
     *    0.640875903798206598D+08,
     *    0.378064033179711401D+08,
     *    0.492781805148246139D+08,
     *    0.121461007977146208D+08,
     *    0.808935577295990288D+08,
     *   -0.349515338742027432D+08,
     *   -0.130992984676782088D+07,
     *    0.484353606550122611D+07,
     *   -0.220063213055212870D+08,
     *    0.134306776349800527D+08,
     *   -0.154525401562347617D+08,
     *    0.302200467051492073D+06,
     *   -0.384482713039138727D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_v0_2=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_v0_3(RR)
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
     *    0.677517552303577042D+10,
     *   -0.980039939734265137D+10,
     *    0.388092604823027802D+10,
     *   -0.110902266297047520D+10,
     *    0.484620841709214568D+09,
     *   -0.188375559125150055D+09,
     *   -0.516986507877669260D+08,
     *   -0.144402328185978651D+09,
     *   -0.117633187439502150D+09,
     *   -0.105333082040884659D+09,
     *   -0.788939870263220221D+08,
     *   -0.520658024563920945D+08,
     *   -0.217799734963362888D+08,
     *    0.923372842667368054D+07,
     *    0.345881319939893261D+08,
     *    0.438914133654673249D+08,
     *    0.463901418529420123D+08,
     *    0.506890395224836394D+08,
     *    0.469617026386186108D+08,
     *    0.582019451724176630D+08,
     *    0.487657961643912122D+08,
     *    0.412047084775665924D+08,
     *    0.421191957380312383D+08,
     *    0.883516366176164895D+07,
     *    0.324348039676273838D+08,
     *   -0.269159825429413319D+08,
     *    0.117634312221577808D+09,
     *    0.913031930310125276D+07,
     *   -0.329246689470835291D+08,
     *    0.324900763548029736D+08,
     *   -0.352933558405040205D+08,
     *    0.115364924574743174D+08,
     *   -0.284161116143819923D+07,
     *   -0.588762103919353150D+07,
     *   -0.531222740620499593D+06/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_v0_3=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_v0_4(RR)
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
     *    0.617158326677711010D+10,
     *   -0.876063612691217422D+10,
     *    0.280884661798855066D+10,
     *   -0.961224387789221406D+09,
     *    0.488556669215859294D+09,
     *   -0.320467886226678118D+08,
     *    0.829779428168724775D+08,
     *   -0.649383917237082962D+07,
     *   -0.371355495360971428D+07,
     *   -0.120680893978045192D+08,
     *   -0.834119152991232276D+07,
     *   -0.725292673809417244D+07,
     *    0.531328921472493443D+06,
     *    0.117321832209778856D+08,
     *    0.254008871145194657D+08,
     *    0.254828451813950092D+08,
     *    0.240303044807130136D+08,
     *    0.157308422561277952D+08,
     *    0.164946967136654966D+08,
     *    0.202811841985239051D+08,
     *    0.206871551806152575D+08,
     *    0.233925983380205855D+08,
     *    0.176232441700203978D+08,
     *    0.108969911082269344D+08,
     *    0.280124340460098349D+07,
     *   -0.162211796709186472D+08,
     *    0.563171452923030853D+08,
     *    0.775761870216983557D+07,
     *   -0.262152393801022097D+08,
     *    0.300900810334736258D+08,
     *   -0.297764102398454547D+08,
     *    0.953168155650613271D+07,
     *   -0.496925966189104039D+07,
     *    0.113670522746274859D+06,
     *   -0.114013152925336431D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_v0_4=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_v0_5(RR)
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
     *    0.461217429490902138D+10,
     *   -0.720072992572966099D+10,
     *    0.282968898492834091D+10,
     *   -0.835152253153198004D+09,
     *    0.411522121840802550D+09,
     *   -0.450631510549229085D+08,
     *    0.757892976085763276D+08,
     *    0.638196247747072857D+07,
     *    0.146493739683052655D+08,
     *    0.716628571847380605D+07,
     *    0.921113259960973263D+07,
     *    0.430164409282412846D+07,
     *    0.470318622743794974D+07,
     *    0.715598296528425533D+07,
     *    0.160907240548391957D+08,
     *    0.161576600836003739D+08,
     *    0.130808148333667889D+08,
     *    0.737537025283837225D+07,
     *    0.486875468627173081D+07,
     *    0.175485647632014239D+07,
     *    0.844673615471202694D+07,
     *    0.123673073449564800D+08,
     *    0.714485929590797611D+07,
     *    0.289763669591167197D+07,
     *    0.976464155344773643D+07,
     *   -0.375649770910991803D+08,
     *    0.395165245553495288D+08,
     *    0.161645709536783136D+08,
     *   -0.298660178148569278D+08,
     *    0.255232092964598909D+08,
     *   -0.178020259094199613D+08,
     *    0.749960824746781401D+07,
     *   -0.109187351068161726D+08,
     *    0.950040437117647380D+07,
     *   -0.530256198599853087D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_v0_5=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_v0_6(RR)
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
     *    0.476982566173464108D+10,
     *   -0.764089782335866165D+10,
     *    0.288585213508764601D+10,
     *   -0.735602629176923752D+09,
     *    0.484828024624088585D+09,
     *   -0.380727063451127335D+07,
     *    0.102522561931984693D+09,
     *    0.238973851913663149D+08,
     *    0.244348312885826714D+08,
     *    0.160830920621736664D+08,
     *    0.136789378769209553D+08,
     *    0.102309636567730159D+08,
     *    0.288665419589183666D+07,
     *    0.350663934261692083D+07,
     *    0.609579811801462248D+07,
     *    0.859546791330591775D+07,
     *    0.645824100696412567D+07,
     *    0.387015109229033114D+07,
     *    0.628654841033376567D+07,
     *   -0.699833533521516528D+07,
     *    0.856996786028444767D+07,
     *    0.659069177904529031D+07,
     *   -0.740416932598389685D+06,
     *    0.186136804986725040D+08,
     *   -0.104356731931949183D+08,
     *   -0.729570059654864389D+07,
     *   -0.351729726463873778D+07,
     *    0.183722397784283459D+08,
     *   -0.185493915427700058D+08,
     *    0.131614388896130715D+08,
     *   -0.700321298561874405D+07,
     *    0.435567735666584922D+06,
     *    0.810239028492451878D+06,
     *   -0.105656035522515280D+07,
     *    0.915750263339141849D+06/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_v0_6=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_v0_7(RR)
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
     *    0.277901011020553398D+10,
     *   -0.432924045097288132D+10,
     *    0.154705044179952264D+10,
     *   -0.457148905305670440D+09,
     *    0.273730099607205451D+09,
     *    0.135100273636536561D+08,
     *    0.748703744977397025D+08,
     *    0.248921885559969768D+08,
     *    0.196700518752739616D+08,
     *    0.119146517786198929D+08,
     *    0.909855356097545289D+07,
     *    0.841431001305358484D+07,
     *    0.347648569439832540D+07,
     *    0.677077706042855862D+06,
     *    0.267254517226238688D+07,
     *    0.422661446138642728D+07,
     *    0.512134843562315311D+07,
     *    0.266634756300450489D+07,
     *    0.492508409705484100D+07,
     *   -0.200177518122392870D+07,
     *   -0.185203100258929608D+07,
     *    0.410105080715974886D+07,
     *   -0.576590639027249347D+07,
     *    0.116108658339938894D+08,
     *   -0.173715315358961560D+07,
     *   -0.431911717493474856D+07,
     *   -0.211823356899501150D+07,
     *    0.805353729835367110D+07,
     *   -0.784038081177155301D+07,
     *    0.222363244509718847D+07,
     *    0.129162495551474183D+07,
     *   -0.289212475147323497D+07,
     *    0.311935924776605237D+07,
     *   -0.138734374336399860D+07,
     *   -0.857513754458879121D+05/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_v0_7=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_v0_8(RR)
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
     *    0.320864431608361101D+10,
     *   -0.563169270738546562D+10,
     *    0.254222541811669350D+10,
     *   -0.575927825751942515D+09,
     *    0.355748202944453299D+09,
     *   -0.262399398007970415D+08,
     *    0.686783122963022590D+08,
     *    0.112945834423579723D+08,
     *    0.151503218255797550D+08,
     *    0.843505787626402453D+07,
     *    0.771434319779101852D+07,
     *    0.563035301971080527D+07,
     *    0.452182460926722549D+07,
     *   -0.391057136623464292D+07,
     *    0.194741026572634303D+07,
     *    0.166576025280055637D+07,
     *   -0.143120995610675076D+06,
     *    0.373929509850142663D+07,
     *    0.371712927579702670D+07,
     *    0.730860444473255426D+06,
     *   -0.290642214079938945D+06,
     *    0.747029262501635589D+06,
     *   -0.145825986197703797D+07,
     *   -0.578496589973646961D+07,
     *    0.872956111242619716D+07,
     *   -0.188703095389100281D+07,
     *   -0.438050550612784829D+07,
     *    0.393695749108811188D+07,
     *    0.182451897074817680D+06,
     *   -0.492785952148701996D+07,
     *    0.439573608627993800D+07,
     *   -0.172890472617617669D+07,
     *    0.624952017235872103D+06,
     *    0.335426825780328130D+06,
     *   -0.963694283633786370D+06/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_v0_8=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl0
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL0_v0_9(RR)
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
     *    0.116110190002285719D+10,
     *   -0.178595051431960750D+10,
     *    0.575421810832156658D+09,
     *   -0.157676181460162252D+09,
     *    0.119918032082672641D+09,
     *    0.163644292342365663D+08,
     *    0.340056479540262893D+08,
     *    0.112542224383108523D+08,
     *    0.705608653342142422D+07,
     *    0.509657711239284091D+07,
     *    0.384486827707689581D+07,
     *    0.420586972172959149D+07,
     *    0.273735939154163888D+07,
     *   -0.172356074323194392D+06,
     *   -0.128481029345794674D+07,
     *    0.139081558777726395D+07,
     *    0.104687822993886797D+07,
     *    0.118399617279725440D+06,
     *    0.187783372781600757D+07,
     *    0.195081589648232865D+07,
     *    0.187805537570319371D+06,
     *   -0.344112788324010512D+07,
     *    0.334646757539036870D+07,
     *   -0.686455725890562870D+07,
     *    0.182449413066698518D+07,
     *    0.664816098228275590D+07,
     *   -0.429329055613289960D+07,
     *    0.196225201969405425D+05,
     *    0.824888452987971716D+06,
     *   -0.573027329980814713D+06,
     *   -0.311558089166373771D+06,
     *    0.623507510665765149D+06,
     *    0.532100253096414381D+06,
     *   -0.184530488380961702D+07,
     *    0.158503166291752667D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL0_v0_9=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL2_v0_10(RR)
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
     *    0.636905674842750072D+09,
     *   -0.129020697860196137D+10,
     *    0.728086829708801627D+09,
     *   -0.139604087938575774D+09,
     *    0.772398842716064304D+08,
     *   -0.192882860177438706D+08,
     *    0.804270904064994864D+07,
     *   -0.237175641008630255D+07,
     *    0.117992043983302335D+07,
     *    0.740481951537411776D+06,
     *   -0.544944764983693720D+06,
     *    0.112412050083140726D+07,
     *   -0.210879613432518486D+07,
     *   -0.514175667547786943D+05,
     *    0.163960784702887689D+07,
     *   -0.235102160919488082D+07,
     *    0.185452746911219694D+07,
     *    0.170934568757311185D+07,
     *   -0.439664738431229163D+07,
     *    0.372275813094899617D+07,
     *   -0.186742223866412509D+07,
     *    0.894791361819264479D+06,
     *   -0.219858040185595769D+07,
     *    0.970227502550709248D+07,
     *   -0.207851562874825858D+08,
     *    0.241078484712500647D+08,
     *   -0.130684720964627042D+08,
     *    0.355675680862437654D+07,
     *   -0.284543057253148686D+07,
     *    0.176176475236889115D+07,
     *   -0.710326767893341603D+06,
     *   -0.774906505015448900D+06,
     *    0.319263535919516208D+07,
     *   -0.259430005250991089D+07,
     *   -0.311620898073950375D+06/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_v0_10=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL2_v0_2(RR)
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
     *    0.557699880277947044D+10,
     *   -0.624725253163322353D+10,
     *    0.183698054561760306D+10,
     *   -0.903119159534327984D+09,
     *    0.112976947530167639D+09,
     *   -0.215139309460971862D+09,
     *   -0.997884559704813659D+08,
     *   -0.138540958398813516D+09,
     *   -0.121511767950127214D+09,
     *   -0.115707253255967826D+09,
     *   -0.974583642795259953D+08,
     *   -0.799256073500085473D+08,
     *   -0.562524257009606212D+08,
     *   -0.356633330246923715D+08,
     *   -0.132075785806477722D+08,
     *    0.681858241375164595D+07,
     *    0.230661680167050846D+08,
     *    0.353992892715934813D+08,
     *    0.422474129659290835D+08,
     *    0.488929223395244107D+08,
     *    0.464942243732484207D+08,
     *    0.462958624816433340D+08,
     *    0.527252680441062450D+08,
     *    0.344666711472032517D+08,
     *    0.692305272251941711D+08,
     *   -0.594196776374670863D+07,
     *    0.127279190500899762D+09,
     *    0.739503163195448369D+08,
     *    0.649913047076319903D+07,
     *    0.288207891404262651D+07,
     *   -0.675446411433657259D+07,
     *   -0.152453449137396808D+07,
     *   -0.712887277563616727D+07,
     *    0.339568938362048427D+07,
     *   -0.420900327364548109D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_v0_2=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL2_v0_3(RR)
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
     *   -0.110056161941127276D+10,
     *    0.139297265185199976D+10,
     *   -0.293928137711912811D+09,
     *    0.172241332222457439D+09,
     *   -0.458231757407660484D+08,
     *    0.136924414012196884D+08,
     *   -0.941938014178624935D+07,
     *   -0.337478033980979118D+07,
     *   -0.853541897901558317D+07,
     *   -0.131853914140430447D+08,
     *   -0.162635126710114703D+08,
     *   -0.197326887114649229D+08,
     *   -0.181319008924542665D+08,
     *   -0.175039624190588370D+08,
     *   -0.150864036514796354D+08,
     *   -0.136192880245118830D+08,
     *   -0.984522038665108755D+07,
     *   -0.696467178588531725D+07,
     *   -0.764926466068208497D+07,
     *    0.282353067775653629D+07,
     *   -0.495253098717514798D+07,
     *    0.734155814480163809D+07,
     *   -0.659241206232764013D+07,
     *    0.140427448538027070D+08,
     *   -0.191381703774914481D+08,
     *    0.127502383502124958D+08,
     *    0.264004126865256950D+07,
     *    0.878064026898320019D+07,
     *    0.331871862052073237D+07,
     *    0.207387980032845214D+07,
     *   -0.335821747802917706D+07,
     *    0.400314675359732425D+07,
     *   -0.657550130831148755D+07,
     *    0.567614711112878006D+07,
     *   -0.339746800484404340D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_v0_3=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL2_v0_4(RR)
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
     *    0.290181304534570789D+10,
     *   -0.444428089981859779D+10,
     *    0.142894737065864015D+10,
     *   -0.388540497263005137D+09,
     *    0.286613146879595637D+09,
     *    0.362200195253765807D+08,
     *    0.905263302080151439D+08,
     *    0.431532067038534433D+08,
     *    0.335582934047516808D+08,
     *    0.165768073412668370D+08,
     *    0.656780694775183778D+07,
     *   -0.394856183818586229D+06,
     *   -0.400799377285283245D+07,
     *   -0.687636700767159648D+07,
     *   -0.445959255295562278D+07,
     *   -0.619554911549653392D+07,
     *   -0.313345963565925928D+07,
     *   -0.288559037207851745D+07,
     *    0.391594406602506817D+06,
     *   -0.230750336074081296D+07,
     *    0.481025127826596517D+07,
     *    0.227784737203444354D+06,
     *    0.673077815095109027D+07,
     *   -0.718220032057878282D+07,
     *    0.798034902998306043D+07,
     *   -0.105543310076045617D+08,
     *    0.461710103691845573D+07,
     *    0.111857739377570879D+08,
     *   -0.295576541132044338D+06,
     *    0.230852322416862333D+07,
     *   -0.148645913888945081D+07,
     *   -0.646770360172940127D+05,
     *    0.517614202881207515D+06,
     *   -0.640974469140543370D+06,
     *    0.910405308111979859D+06/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_v0_4=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL2_v0_5(RR)
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
     *   -0.711433542553594351D+09,
     *    0.128577243013030815D+10,
     *   -0.536586261684719026D+09,
     *    0.739839227737180144D+08,
     *   -0.100931000244962811D+09,
     *   -0.305342127776865056D+07,
     *   -0.125627495974449757D+08,
     *    0.265479481444472773D+07,
     *    0.317707401770288777D+07,
     *    0.371535876716061728D+07,
     *    0.154903738269069721D+07,
     *    0.685143213361830276D+05,
     *   -0.209746291253840318D+07,
     *   -0.204025722765863268D+07,
     *   -0.439405721434761363D+06,
     *   -0.212782936426520022D+07,
     *   -0.162586875185792847D+07,
     *   -0.413018745899144924D+05,
     *    0.445170770038122195D+06,
     *   -0.499587834191429801D+07,
     *    0.405951778953305818D+07,
     *    0.423985846568766888D+07,
     *   -0.100762269271380827D+08,
     *    0.102157847959735096D+08,
     *    0.287670840342416288D+07,
     *   -0.899278125448114052D+07,
     *   -0.126857746978211543D+07,
     *    0.486829432165133674D+07,
     *    0.369889548559936578D+06,
     *    0.537684892155467765D+06,
     *   -0.303689389792799775D+06,
     *   -0.164705868021505815D+07,
     *    0.328070595084952563D+07,
     *   -0.862846464547522832D+06,
     *   -0.191817725468137418D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_v0_5=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL2_v0_6(RR)
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
     *    0.227369013256841803D+10,
     *   -0.416008266739053249D+10,
     *    0.190505865984697390D+10,
     *   -0.323448373229954898D+09,
     *    0.275910751819743514D+09,
     *   -0.208677390370645188D+08,
     *    0.390221351274639592D+08,
     *    0.180083102439951920D+07,
     *    0.691583141768497135D+07,
     *    0.249449343829181930D+07,
     *    0.154798084386450564D+07,
     *    0.156170992998979893D+07,
     *   -0.580060451508121216D+06,
     *    0.331540630694518622D+06,
     *   -0.218963066917348467D+07,
     *    0.789416940417197882D+06,
     *   -0.129372480657914979D+07,
     *   -0.161318620340735232D+07,
     *    0.583742910900342744D+07,
     *   -0.110099797371312492D+08,
     *    0.763265242229388282D+07,
     *   -0.338973951322273724D+07,
     *   -0.923097080578689277D+07,
     *    0.119949601437342316D+08,
     *   -0.665329240788710304D+05,
     *   -0.323550397363493405D+07,
     *    0.356313247810042556D+07,
     *   -0.292743930485270265D+07,
     *    0.266384143365175789D+07,
     *   -0.193686716483137244D+07,
     *    0.174020768274276285D+07,
     *   -0.238267625864995550D+07,
     *    0.439108725246513728D+07,
     *   -0.421476838349953108D+07,
     *    0.240803499960590200D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_v0_6=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL2_v0_7(RR)
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
     *   -0.492245390906152189D+09,
     *    0.107619744980248785D+10,
     *   -0.685082151729213715D+09,
     *    0.141977040869015723D+09,
     *   -0.659701368867858052D+08,
     *    0.241932852405317798D+08,
     *   -0.467666994646243006D+07,
     *    0.376105969846729701D+07,
     *   -0.122299208940234603D+06,
     *    0.133537579078195011D+07,
     *    0.143712048738935334D+07,
     *    0.123587392362282518D+07,
     *    0.673746467049226339D+05,
     *    0.117257047571065207D+07,
     *   -0.329522308900472801D+07,
     *    0.162428484145748150D+07,
     *   -0.400078565159763023D+07,
     *    0.434959553847656958D+07,
     *   -0.146205799100256758D+07,
     *   -0.976560175721534295D+06,
     *    0.274228799515426392D+07,
     *   -0.285854211622587405D+07,
     *   -0.103017778669942846D+07,
     *   -0.826676925246696919D+04,
     *    0.278544083862312604D+06,
     *   -0.373783783494386775D+07,
     *    0.676582559455903154D+07,
     *   -0.329393174449633248D+07,
     *    0.408292578011030378D+07,
     *   -0.466629884943075292D+07,
     *    0.213622151849628054D+07,
     *    0.158964191425781697D+07,
     *   -0.331246426807687059D+07,
     *    0.188607903451994970D+07,
     *    0.152573242663507699D+06/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_v0_7=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL2_v0_8(RR)
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
     *    0.172732318260711479D+10,
     *   -0.341174099693365908D+10,
     *    0.182859450039141393D+10,
     *   -0.331994105728986561D+09,
     *    0.209120363079807341D+09,
     *   -0.430553380543953180D+08,
     *    0.238269865153793208D+08,
     *   -0.470140877083909977D+07,
     *    0.322452168589429324D+07,
     *   -0.172954167355673504D+07,
     *    0.172457190904876986D+07,
     *   -0.178055431929770904D+07,
     *    0.292807286047906941D+07,
     *   -0.782486265392334200D+06,
     *    0.141575887750610855D+06,
     *    0.140892047419081762D+06,
     *   -0.756130568219428300D+06,
     *   -0.203264433404238382D+07,
     *   -0.241880217389070516D+06,
     *    0.485807382966472208D+07,
     *   -0.625345826071981713D+07,
     *    0.441694769712568820D+07,
     *    0.304604568356459960D+06,
     *   -0.546153080329078622D+07,
     *    0.121530063273045309D+08,
     *   -0.180657830482467525D+08,
     *    0.114502930761773512D+08,
     *   -0.272119734962271946D+07,
     *    0.131801883099487913D+07,
     *    0.536718498029134236D+06,
     *   -0.146559889397471747D+07,
     *    0.253053755112899514D+07,
     *   -0.484181423460538872D+07,
     *    0.412360180674585840D+07,
     *   -0.148037445306411362D+07/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_v0_8=SUMA
       RETURN
       END
C PROGRAM TO GENERATE 1-D RKHS FITS OF Vl2
C DISTANCE IS IN BOHR
C EXPANSION COEFFFICIENTS IN CM-1
C REPRODUCING KERNEL HILBERT SPACE METHOD 
C AB INITIO: UCCSD(T)/CBS BSSE CORRECTED 
C Contact: Jacek Klos jklos@umd.edu
C---OPTIMIZED FOR EXECUTION TIME BY DIRECT
C---EXPRESIONS OF HYPERGEOMETRIC AND BETA FUNCTIONS
C
      FUNCTION VL2_v0_9(RR)
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
     *   -0.233075681908708483D+09,
     *    0.554175249230122328D+09,
     *   -0.409575793071995616D+09,
     *    0.105285998951568171D+09,
     *   -0.321772887701118663D+08,
     *    0.159813032627533283D+08,
     *   -0.326394303656015033D+07,
     *    0.269342684258913063D+07,
     *   -0.148166431411186029D+06,
     *    0.962699819266577950D+06,
     *   -0.834382461785805062D+06,
     *    0.292867644573944504D+06,
     *   -0.907584332995778532D+06,
     *    0.350270016857351220D+05,
     *    0.750923817363097682D+06,
     *    0.692077219892069814D+06,
     *   -0.314350670520343294D+06,
     *    0.202632029775614676D+06,
     *   -0.170114158321288181D+07,
     *    0.101285998101547547D+07,
     *   -0.316416069266698789D+06,
     *   -0.231470919912292575D+07,
     *    0.738932870421846304D+07,
     *   -0.100116393775518443D+08,
     *    0.112503250246319771D+08,
     *   -0.123687024990233537D+08,
     *    0.743677171301179007D+07,
     *   -0.235028332283700211D+07,
     *    0.345428463328313641D+07,
     *   -0.542322744519502018D+07,
     *    0.467867079135303013D+07,
     *   -0.221146410487418389D+07,
     *    0.539825523205072619D+06,
     *   -0.270459325132685990D+03,
     *    0.470409551634216099D+06/
C Kernel parameters
      M=5
C
C      EVALUATE PES    
       SUMA=ZERO
       DO I=1,NPT
        SUMA=SUMA+VCOEF(I)*RKHS_DKN2(R(I),RR,M)
       END DO
       VL2_v0_9=SUMA
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
