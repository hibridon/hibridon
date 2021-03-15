c   Define the sizes of grids
c       V2MAX: maximum value of v2 included in the pot file
c       V2TMAX: number of (v2, v2') combination, C(V2MAX+1, 2)
c       NVLM: number of v_lm terms for each (v2, v2') combination
c       NVVL: total number of v_lm terms, for all (v2, v2') blocks
c       NTHETA, NPHI: number of theta/phi's in the ab initio calculation
c       NANGLE: number of (theta, phi) tuples
c       NDIST: number of distances included in the ab initio calculation
      integer V2MAX, V2TMAX, NVLM, NTHETA, NPHI, NANGLE, NVVL, NDIST
      parameter (V2MAX=3, V2TMAX=(V2MAX+1)*(V2MAX+2)/2)
      parameter (NVLM=12, NVVL=NVLM*V2TMAX)
      parameter (NTHETA=19, NPHI=7, NANGLE=NTHETA*NPHI)
      parameter (NDIST=19)
c
c   Conversion factor
c       ECONV: hartree to wavenumber
c       XMCONV: amu to electrom mass
      double precision ECONV, XMCONV
      parameter (ECONV=219474.6d0, XMCONV=1822.88848477d0)
c
c   Max number of channels
      integer KMAX
      parameter (KMAX=10000)
c
c   Lengths of cod array, 
c       ICOD, IRCOD, LENCOD: lenghts of cod array
      integer ICOD, IRCOD, LENCOD
      parameter (ICOD=5, IRCOD=4, LENCOD=ICOD+IRCOD+3)
c
c   ch3he block: data used only by this pot/basis combination
c       brot, crot: rotational constants of CH3 for each vibrational level
c       evib: vibrational level energies
c       nlamsi: number of v_lm terms for each (v2, v2') combination
c       lamsym, musym: list of lambda/mu's used for the coupling potential symmetric to theta = 90 deg
c       lamasy, muasy: list of lambda/mu's used for the coupling potential anti-symmetric to theta = 90 deg
c
c     Source of rotational constants:
c     Yamada, C., et. al., JCP, 75, 5256
c     Amano, T., et. al., JCP, 77, 5284
c
      common /ch3he/ brot, crot, evib,
     +               lamsym, musym, lamasy, muasy
      double precision brot(V2MAX+1), crot(V2MAX+1), evib(V2MAX+1)
      integer lamsym(NVLM), musym(NVLM), lamasy(NVLM), muasy(NVLM)
      data brot /9.57789d0, 9.25814d0, 8.93320d0, 8.60974d0/
      data crot /4.74275d0, 4.811643d0, 4.871213d0, 4.92655d0/
      data evib /0d0, 6.064531d2, 1.28809d3, 2.0191657d3/
c   lambda/mu's used in this pot routine
c
c   When <v2'|V|v2> is symmetric about theta=90 deg (even v2 + v2')
c       lambda =  0  2  4  6  8  3  5  7  9  6  8  9
c       mu     =  0  0  0  0  0  3  3  3  3  6  6  9
c   When <v2'|V|v2> is anti-symmetric about theta=90 deg (odd v2 + v2')
c       lambda =  1  3  5  7  9  4  6  8 10  7  9 10
c       mu     =  0  0  0  0  0  3  3  3  3  6  6  9
      data lamsym /0, 2, 4, 6, 8, 3, 5, 7, 9, 6, 8, 9/
      data musym /0, 0, 0, 0, 0, 3, 3, 3, 3, 6, 6, 9/
      data lamasy /1, 3, 5, 7, 9, 4, 6, 8, 10, 7, 9, 10/
      data muasy /0, 0, 0, 0, 0, 3, 3, 3, 3, 6, 6, 9/
c
c
c   covvl block
c       vvl: r-dependence of each term in potential expansion
c   represent <v_2'|v_{\lambda\mu}(Q_2, R)|v_2> in the following order:
c       <0|v_00|0>, <0|v_20|0>, ..., <0|v_99|0>, <0|v_10|1>, ...
c   where v2 >= v2'
      common /covvl/ vvl
      double precision vvl(NVVL)
c
c   cosysi block
c       nscod: total number of variable names which are passed to HINPUT, nscod must equal isrcod + isicod + 3
c       isicod: total number of integer system dependent variables
c       nterm: number of different associated legendre terms in expansion of potential
c       numpot: the number of the potential used, this variable is passed to the pot subroutine
c       ipotsy: cylindrical symmetry of potential. Should be set to 3 for CH3.
c       iop: ortho/para label for molecular states. Only para states are included if iop=1 and only ortho states if iop=-1.
      common /cosysi/ nscode, isicod, nterm, ipotsy, iop, jmax, vmax
      integer nscode, isicod, nterm, ipotsy, iop, jmax, vmax
c
c   cosys block
c       scod: contains names of all system dependent parameters
      common /cosys/ scod
      character*8 scod(LENCOD)
c
c   cosysr block
c       isrcod: total number of real system dependent variables
c       junkr: junk variable (required by hibridon)
c       vmax: maximum value of v2 (starts from zero) included in the calculation
c       emax0, emax1, emax2, emax3: maximum total energy of a level to be included in the channel basis, for four vibrational levels
      common /cosysr/ isrcod, junkr, emax0, emax1, emax2, emax3
      integer isrcod, junkr
      double precision emax0, emax1, emax2, emax3
c
c   coiout block
c       niout: number of level indeces included in the output of hibridon
c       indout: level indeces included in the output of hibridon
      common /coiout/ niout, indout
      integer niout, indout(100)
c
c   conlam block
c       nlam: the number of angular coupling terms actually used
c       nlammx: the maximum number of angular coupling terms
c       lamnum: number of non-zero v2 matrix elements for each lambda
      common /conlam/ nlam, nlammx, lamnum
      integer nlam, nlammx, lamnum(NVVL)
c
c   coered block
c       ered: collision energy in atomic units (hartrees)
c       rmu: collision reduced mass in atomic units
      common /coered/ ered, rmu
      double precision ered, rmu
c
c   cov2 block
c       nv2max: maximum core memory allocated for the v2 matrix
c       junkv: junk variable, required by hibridon
c       v2: lower triangle of nonzero angular coupling matrix elements stored in packed column form that is (1,1), (2,1), (3,1) ... (n,1), (2,2), (3,2) ... (n,2), etc.
      common /cov2/ nv2max, junkv, v2
      integer nv2max, junkv
      double precision v2(1)
c
c   coiv2 block
c       iv2: row+column index of v2 matrix for each non-zero element
      common /coiv2/ iv2
      integer iv2(1)
c
c   cocent block
c       cent: array containing centrifugal barrier of each channel
      common /cocent/ cent
      double precision cent(KMAX)
c
c   coeint block
c       eint: array containing channel energies (in hartree). The zero of energy is assumed to be the v2=0, j=0, k=0 level
      common /coeint/ eint
      double precision eint(KMAX)
c
c   coipar block
c       junkip: not refered to
c       iprint: level of printing in calculation
      common /coipar/ junkip, iprint
      integer junkip(9), iprint
