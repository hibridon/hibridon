* ----------------------------------------------------------------------
      subroutine basphtp (j, l, is, jhold, ehold, ishold, nlevel,
     :       nlevop, etemp, fjtemp, fitemp, fistmp, rcut, jtot,
     :       flaghf, flagsu, csflag, clist, bastst, ihomo,
     :       nu, numin, jlpar, n, nmax, ntop)
* --------------------------------------------------------------------
*  subroutine to determine angular coupling potential for collision
*  of an spherical top molecule with a structureless atom or with an
*  uncorrugated surface
*
*  use procedure from heilmen et al. [j.c.p. 107, 902 (1997)] for
*  picking tetrahedral combinations of symmetric top wave functions
*  coefficients obtained from output of a program provided by
*  ad van der avoird.
*
*  this version includes levels for jmax .le. 16 for A levels
*                                   jmax .le. 10 for E and F levels
*
*  author:  paul dagdigian
*  revision:  17-feb-2016 (p.dagdigian)
* --------------------------------------------------------------------
*  variables in call list:
*    j:        on return contains rotational quantum number for each
*              channel
*    l:        on return contains orbital angular momentum for each
*              channel
*    is:       on return contains state index for each channel,
*              set to the order in energy of the near degenerate levels
*              of the given j
*    jhold:    on return contains rotational quantum numbers for each
*              rotational level
*    ehold:    on return contains energy in hartrees of each rotational
*              level
*    ishold:   on return contains symmetry index of each energetically
*              distinct level
*    nlevel:   on return contains number of energetically distinct
*              rotational levels used in channel basis
*    nlevop:   on return contains number of energetically distinct
*              rotational levels used in channel basis which are open
*              asymptotically
*    etemp:    scratch array used to create channel list
*    fjtemp:   scratch array used to create channel list
*    fitemp:   scratch array used to create channel list
*    fistmp:   scratch array used to create channel list
*    rcut:     cut-off point for keeping higher energy channels
*              if any open channel is still closed at r=rcut, then
*              all closed channels as well any open channels which are
*              still closed at r=rcut are dropped from basis
*              note that this is ignored in molecule-surface collisions!!!
*    jtot:     total angular momentum
*              in cc calculation jtot is the total angular momentum
*              in cs calculation jtot is the l-bar quantum number
*    flaghf:   if .true., then system has half-integer spin
*              if .false., then system has integer spin
*    flagsu:   if .true., then molecule-surface collisons
*    csflag:   if .true., then coupled-states calculation
*              if .false., the close-coupled calculation
*    clist:    if .true., then quantum numbers and energies listed for
*              each channel
*    bastst:   if .true., then execution terminates after the first call
*              to basis
*    nu:       coupled-states projection index
*    numin:    minimum coupled states projection index
*              for cc calculations nu and numin are both set = 0 by calling
*              program
*    jlpar:    total parity of included channels in cc calculation
*              only those channels are included for which
*                  (-1)**(is+j+kp+l-jtot)=jlpar
*              where parity designates the parity of the molecular state
*              in cs calculation jlpar is set equal to 1 in calling program
*    n:        on return equals number of channels
*    nmax:     maximum dimension of arrays
*    ntop:     maximum row dimension of all matrices passed to subroutines
*              propag and soutpt.  ntop is set in basis only if nu = numin
*              otherwise it is unchanged from the value supplied by the
*              calling program
*  variables in common bloc /cosysr/
*    isrcod:   total number of real system dependent variables
*    brot:     rotational constant for the spherical top
*  variables in common block /cosysi/
*    nscode:   total number of variable names which are passed to HINPUT
*              nscod must equal isrcod + isicod + 3
*    isicod:   total number of integer system dependent variables
*    nterm:    number of different tetrahedral terms in expansion of the potential
*    numpot:   the number of the potential used, this variable is passed
*              to the pot subroutine
*    iop:      nuclear spin modification for molecular states:
*                = 1 for A levels (A1 or A2)
*                = 2 for E levels
*                = 3 for F levels (F1 or F2)
*    jmax:     the maximum rotational angular momentum for the spherical top
*  variable in common block /cocent/
*    cent:      array containing centrifugal barrier of each channel
*  variable in common block /coeint/
*    eint:      array containing channel energies (in hartree)
*  variables in common block /coered/
*    ered:      collision energy in atomic units (hartrees)
*    rmu:       collision reduced mass in atomic units
*               (mass of electron = 1)
*  variables in common block /conlam/
*    nlam:      the number of angular coupling terms actually used
*    nlammx:    the maximum number of angular coupling terms
*    lamnum:    number of non-zero v2 matrix elements for each lambda
*               lamnum is an array of dimension nlammx
*  variables in common block /cov2/
*    nv2max:    maximum core memory allocated for the v2 matrix
*    v2:        lower triangle of nonzero angular coupling matrix elements
*               stored in packed column form that is (1,1), (2,1), (3,1) ...
*               (n,1), (2,2), (3,2) ... (n,2), etc.
*  variable in common block /coiv2/
*   iv2:        row+column index of v2 matrix for each non-zero element
*  variables in common block /coconv/
*   econv:      conversion factor from cm-1 to hartrees
*   xmconv:     converson factor from amu to atomic units
*  variables in common block /coatpi/
*   narray:     maximum size of symmetric top basis fn expansion
*               set to 200 in himain
*  variable in common block /coatpr/
*   c:          expansion coefficients for spherical top rotor wave fns.
*               in a signed-k symmmetric top basis |j,k,m>.
*               if c(+k) = c(-k), then eps = +1, or if c(+k) = -c(-k), eps = -1.
*               the expansion coefficients for a given eps stored in the array c are:
*               c(k=-j), c(k=-j+2), ... c(k=j-10, c(k=j).
*  variable in common block /coatp1/
*   ctemp:      temporary storage for rot. e.fn coeffs.
*  variable in common block /coatp2/
*   chold:      ditto
*  subroutines called:
*   vlmats:     returns angular coupling coefficient for particular
*               choice of channel index
*   prmats:     computes primitive cc and cs v-lambda matrix elements
*               between signed-k symmetric top basis fns.
* --------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical flaghf, csflag, clist, flagsu, ihomo, bastst
      include "common/parbas"
      include "common/parbasl"
      common /cosysi/ nscode, isicod, nterm, iop, jmax
      common /coipar/ iiipar(9), iprint
      common /cosysr/ isrcod, junkr, brot
      common /cov2/ nv2max, junkv, v2(1)
      common /coiv2/ iv2(1)
      common /conlam/ nlam, nlammx, lamnum(1)
      common /cocent/ cent(1)
      common /coeint/ eint(1)
      common /coered/ ered, rmu
      common /coconv/ econv, xmconv
      common /coatpi/ narray
      common /coatpr/ c(1)
      common /coatp1/ ctemp(1)
      common /coatp2/ chold(1)
      dimension j(1), l(1), is(1), jhold(1), ehold(1),
     :          ishold(1), etemp(1), fjtemp(1), fitemp(1),
     :          fistmp(1)
*  scratch arrays
      dimension sc1(narray), ieps(narray), ips(narray)
      dimension na(17), ne(11), nf(11), ipa(11), nepsa(25),
     :     nepse(20), nepsf(30), ua(25,17), ue(20,11),
     :     uf(30,11)
*
*  number of levels of A, E, F type for each j (=0 to 12 for A levels,
*  to 10 for E and F levels)
      data na / 1, 0, 0, 1, 1, 0, 2, 1, 1, 2, 2, 1, 3, 2,
     :     2, 3, 3 /,
     :     ne / 0, 0, 2, 0, 2, 2, 2, 2, 4, 2, 4 /,
     :     nf / 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5 /,
     :     jmax0 /10/, jmaxa /16/
      data nepsa / 1, 1, 1, -1, 1, 1, 1, 1, -1, 1, -1, -1,
     :     1, -1, 1, -1, 1, -1, 1, 1, -1, -1, 1, -1, 1 / ,
     :     nepse / 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1,
     :     1, -1, 1, -1, 1, -1, 1, -1 /,
     :     nepsf / 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, -1, 1,
     :     1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1,
     :     1, -1, 1, -1, 1 /
*  tetrahedral eigenfunctions
      data (ua(1,ii), ii=1,11) /1.d0, 10*0.d0/                            ! j=0 A
      data (ua(2,ii), ii=1,11) /.707107d0, 0.d0, .707107d0, 8*0.d0/       ! j=3 A
      data (ua(3,ii), ii=1,11) /-.456435d0, 0.d0, .763763d0, 0.d0,        ! j=4 A
     :   -.456435d0, 6*0.d0/
      data (ua(4,ii), ii=1,11) /-.395285d0, 0.d0, -.586302d0, 0.d0,       ! j=6 A
     :   .586302d0, 0.d0, .395285d0, 4*0.d0/
      data (ua(5,ii), ii=1,11) /0.d0, .661438d0, 0.d0, .353553d0, 0.d0,   ! j=6 A (2)
     :    .661438d0, 5*0.d0/
      data (ua(6,ii), ii=1,11) /-.478714d0, 0.d0, .520416d0, 0.d0,        ! j=7 A
     :    .520416d0, 0.d0, -.478714d0, 4*0.d0/
      data (ua(7,ii), ii=1,11) /.411425d0, 0.d0, -.270031d0, 0.d0,        ! j=8 A
     :    .718070d0, 0.d0, -.270031d0, 0.d0, .411425d0, 2*0.d0/
      data (ua(8,ii), ii=1,11) /0.d0, .637377d0, 0.d0, .306186d0, 0.d0,   ! j=9 A
     :    .306186d0, 0.d0, .637377d0, 0.d0, 2*0.d0/
      data (ua(9,ii), ii=1,11) /-.381881d0, 0.d0, -.595119d0, 0.d0,       ! j=9 A (2)
     :    0.d0, 0.d0, .595119d0, 0.d0, .381881d0, 2*0.d0/
      data (ua(10,ii), ii=1,11) /0.d0, -.493447d0, 0.d0, .414578d0,       ! j=10 A
     :    0.d0, .414578d0, 0.d0, .414578d0, 0.d0, -.493447d0, 0.d0/
      data (ua(11,ii), ii=1,11) /.407450d0, 0.d0, .111220d0, 0.d0,        ! j=10 A (2)
     :    -.567111d0, 0.d0, .567111d0, 0.d0, -.111220d0,
     :    0.d0, -.407450d0/
      data (ua(12,ii), ii=1,11) /-0.416146d0, 0.d0, 0.324760d0, 0.d0,     ! j=11 A
     :    -0.470483d0, 0.d0, -0.470483d0, 0.d0, 0.324760d0, 0.d0,
     :    -0.416146d0/
      data (ua(13,ii), ii=1,13) / 0.408103d0, 0.d0, -0.66480d0, 0.d0,     ! j=12 A (2)
     :    0.329221d0, 0.d0, -0.664297d0, 0.d0, 0.329221d0, 0.d0,
     :    -0.066480d0, 0.d0, 0.408103d0/
      data (ua(14,ii), ii=1,13) / 0.d0, -0.368874d0, 0.d0, -0.584634d0,   ! j=12 A
     :    0.d0, -0.148780d0, 0.d0, 0.148780d0, 0.d0, 0.584634d0, 0.d0,
     :    0.368874d0, 0.d0/
      data (ua(15,ii), ii=1,13) / 0.016765d0, 0.d0, 0.617513d0, 0.d0,     ! j=12 A (2)
     :    0.311737d0, 0.d0, 0.205993d0, 0.d0, 0.311737d0, 0.d0,
     :    0.617513d0, 0.d0, 0.0167765d0/
      data (ua(16,ii), ii=1,13) / 0.405849d0, 0.d0, 0.161374d0, 0.d0,     ! j=13 A
     :    -0.556098d0, 0.d0, 0.d0, 0.d0, 0.556098d0, 0.d0, -0.161374d0,
     :    0.d0, -0.405849d0/
      data (ua(17,ii), ii=1,13) / 0.d0, 0.504537d0, 0.d0, -0.348919d0,    ! j=13 A (2)
     :    0.d0, -0.351707d0, 0.d0, -0.351707d0, 0.d0, -0.348919d0,
     :    0.d0, 0.504537d0, 0.d0/
      data (ua(18,ii), ii=1,15) / -0.408198d0, 0.d0,  -0.014264d0,        ! j=14 A
     :    0.d0, 0.175738d0, 0.d0, -0.54980d0, 0.d0, 0.549806d0,
     :    0.d0, -0.175738d0, 0.d0, 0.014264d0, 0.d0, 0.408198d0/
      data (ua(19,ii), ii=1,15) / 0.d0, 0.421682d0, 0.d0, -0.347283d0,    ! j=14 A (2)
     :    0.d0, 0.323630d0, 0.d0, 0.440096d0, 0.d0, 0.323630d0, 0.d0,
     :    -0.347283d0, 0.d0, 0.421682d0, 0.d0/
      data (ua(20,ii), ii=1,15) / -0.407422d0, 0.d0, 0.120456d0, 0.d0,    ! j=15 A
     :    -0.384727d0, 0.d0, 0.414105d0, 0.d0, 0.414105d0, 0.d0,
     :    -0.384727d0, 0.d0, 0.120456d0, 0.d0, -0.407422d0/
      data (ua(21,ii), ii=1,15) / 0.d0, -0.356305d0, 0.d0, -0.572822d0,   ! j=15 A (2)
     :    0.d0, -0.211948d0, 0.d0, 0.d0, 0.d0, 0.211948d0, 0.d0,
     :    0.572822d0, 0.d0, 0.356305d0, 0.d0/
      data (ua(22,ii), ii=1,15) / 0.038397d0, 0.d0, 0.636059d0, 0.d0,     ! j=15 A (3)
     :    0.298110d0, 0.d0, 0.071305d0, 0.d0, -0.071305d0, 0.d0,
     :    -0.298110d0, 0.d0, -0.636059d0, 0.d0, -0.038397d0/
      data (ua(23,ii), ii=1,17) / -0.408150d0, 0.d0, 0.016655d0, 0.d0,    ! j=16 A
     :     -0.087764d0, 0.d0, 0.366930d0, 0.d0, -0.617731d0, 0.d0,
     :     0.366930d0, 0.d0, -0.087764d0, 0.d0, 0.016655d0, 0.d0,
     :     -0.408150d0/
      data (ua(24,ii), ii=1,17) / 0.d0, 0.403537d0, 0.d0, 0.198259d0,     ! j=16 A (2)
     :     0.d0, -0.509113d0, 0.d0, -0.196610d0, 0.d0, 0.196610d0,
     :     0.d0, 0.509113d0, 0.d0, -0.198259d0, 0.d0, -0.403537, 0.d0/
      data (ua(25,ii), ii=1,17) / -0.009488d0, 0.d0, -0.513140d0, 0.d0,   ! j=16 A (3)
     :     0.298222d0, 0.d0, 0.326082d0, 0.d0, 0.287511d0,0.d0,
     :     0.326082d0, 0.d0, 0.298222d0, 0.d0, -0.513140d0, 0.d0,
     :     -0.009488d0/
*
      data (ue(1,ii), ii=1,11) / 0.d0, 1.0d0, 0.d0,                       ! j=2 E
     :    8*0.d0/
      data (ue(2,ii), ii=1,11) /-.707107d0, 0.d0, .707107d0,              ! j=2 E (2)
     :    8*0.d0/
      data (ue(3,ii), ii=1,11) / .540061d0, 0.d0, .645497d0,              ! j=4 E
     :    0.d0, .540061d0, 6*0.d0/
      data (ue(4,ii), ii=1,11) /0.d0, .707107d0, 0.d0, -.707107d0,        ! j=4 E (2)
     :     0.d0, 6*0.d0/
      data (ue(5,ii), ii=1,11) / 0.d0, .707107d0, 0.d0,                   ! j=5 E
     :    .707107d0, 0.d0, 6*0.d0/
      data (ue(6,ii), ii=1,11) / -.707107d0, 0.d0, 0.d0, 0.d0,            ! j=5 E (2)
     :    .707107d0, 6*0.d0/
      data (ue(7,ii), ii=1,11) / 0.d0, -.250000d0, 0.d0, .935415d0,       ! j=6 E
     :    0.d0, -.250000d0, 0.d0, 4*0.d0/
      data (ue(8,ii), ii=1,11) / 0.58630196d0, 0.d0, -.395285d0, 0.d0,    ! j=6 E (2)
     :    .395285d0, 0.d0, -0.58630196d0, 4*0.d0/
      data (ue(9,ii), ii=1,11) / -.520416d0, 0.d0, -.478714d0, 0.d0,      ! j=7 E
     :    -.478714d0, 0.d0, -.520416d0, 4*0.d0/
      data (ue(10,ii), ii=1,11) / 0.d0, -.707107d0, 0.d0, 0.d0, 0.d0,     ! j=7 E (2)
     :    .707107d0, 0.d0, 4*0.d0/
      data (ue(11,ii), ii=1,11) /-.424490d0, 0.d0, .278607d0, 0.d0,       ! j=8 E
     :    .695970d0, 0.d0, .278606d0, 0.d0, -.424490d0, 2*0.d0/
      data (ue(12,ii), ii=1,11) / 0.d0, -.536942d0, 0.d0, .460102d0,      ! j=8 E (2)
     :    0.d0, -.460102d0, 0.d0, .536942d0, 0.d0, 2*0.d0/
      data (ue(13,ii), ii=1,11) / .387992d0, 0.d0, .591154d0, 0.d0,       ! j=8 E (3)
     :    0.d0, 0.d0, .591154d0, 0.d0, .387992d0, 2*0.d0/
      data (ue(14,ii), ii=1,11) / 0.d0, -.460102d0, 0.d0, -.536942d0,     ! j=8 E (4)
     :    0.d0, .536942d0, 0.d0, .460102d0, 0.d0, 2*0.d0/
      data (ue(15,ii), ii=1,11) / 0.d0, -.306186d0, 0.d0, .636378d0,      ! j=9 E
     :    0.d0, .637378d0, 0.d0, -.306186d0, 0.d0, 2*0.d0/
      data (ue(16,ii), ii=1,11) / .595119d0, 0.d0, -.381882d0, 0.d0,      ! j=9 E (2)
     :    0.d0, 0.d0, .381882d0, 0.d0, -.595119d0, 2*0.d0/
      data (ue(17,ii), ii=1,11) / 0.d0, .222741d0, 0.d0, -.187140d0,      ! j=10 E
     :    0.d0, .911444d0, 0.d0, -.187140d0, 0.d0, .222741d0, 0.d0/
      data (ue(18,ii), ii=1,11) /-.531788d0, 0.d0, .343800d0, 0.d0,       ! j=10 E (2)
     :    -.314648d0, 0.d0, .314648d0, 0.d0, -.3438000d0, 0.d0,
     :     .531788d0/
      data (ue(19,ii), ii=1,11) / 0.d0, .454959d0, 0.d0, .541391d0,       ! j=10 E (3)
     :    0.d0, 0.0, 0.d0, .541391d0, 0.d0, .454859d0, 0.d0/
      data (ue(20,ii), ii=1,11) /.226242d0, 0.d0, .607809d0, 0.d0,        ! j=10 E (4)
     :    .281748d0, 0.d0, -.281748d0, 0.d0, -.607809d0, 0.d0,
     :   -.226242d0/
*
      data (uf(1,ii), ii=1,11) /1.00000d0, 10*0.d0/                       ! j=1 F1
      data (uf(2,ii), ii=1,11) /.707107d0, 0.d0, .707107d0,               ! j=2 F2
     :    8*0.d0/
      data (uf(3,ii), ii=1,11) /0.d0, 1.d0, 0.d0, 8*0.d0/                 ! j=3 F1
      data (uf(4,ii), ii=1,11) /-.707107d0, 0.d0, .707107d0, 8*0.d0/      ! j=3 F2
      data (uf(5,ii), ii=1,11) /-.707107d0, 0.d0, 0.d0, 0.d0,             ! j=4 F1
     :    .707107d0, 6*0.d0/
      data (uf(6,ii), ii=1,11) /0.d0, .707107d0, 0.d0, .707107d0,         ! j=4 F2
     :    0.d0, 6*0.d0/
      data (uf(7,ii), ii=1,11) /0.d0, -.707107d0, 0.d0,.707107d0,         ! j=5 F1
     :    0.d0, 6*0.d0/
      data (uf(8,ii), ii=1,11) /0.d0, 0.d0, 1.000000d0, 0.d0, 0.d0,       ! j=5 F2
     :    6*0.d0/
      data (uf(9,ii), ii=1,11) /.707107d0, 0.d0, 0.d0, 0.d0,              ! j=5 F2
     :    .707107d0, 6*0.d0/
      data (uf(10,ii), ii=1,11) /0.d0, 0.d0, .707107d0, 0.d0,             ! j=6 F2
     :    .707107d0, 0.d0, 0.d0, 4*0.d0/
      data (uf(11,ii), ii=1,11) /0.d0, -.707107d0, 0.d0, 0.d0,            ! j=6 F1
     :     0.d0, .707107d0, 0.d0, 4*0.d0/
      data (uf(12,ii), ii=1,11) /.707107d0, 0.d0, 0.d0, 0.d0, 0.d0,       ! j=6 F2
     :     0.d0, .707107d0, 4*0.d0/
      data (uf(13,ii), ii=1,11) /0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0,      ! j=7 F1
     :     0.d0, 4*0.d0/
      data (uf(14,ii), ii=1,11) /0.d0, 0.d0, -.707107d0, 0.d0,            ! j=7 F1
     :     .707107d0, 0.d0, 0.d0, 4*0.d0/
      data (uf(15,ii), ii=1,11) /0.d0, .70710d0, 0.d0, 0.d0, 0.d0,        ! j=7 F2
     :     .707107d0, 0.d0, 4*0.d0/
      data (uf(16,ii), ii=1,11) /-.707107d0, 0.d0, 0.d0, 0.d0, 0.d0,      ! j=7 F1
     :     0.d0, .707107d0, 4*0.d0/
      data (uf(17,ii), ii=1,11) /0.d0, 0.d0, 0.d0, .707107d0, 0.d0,       ! j=8 F2
     :     .707107d0, 0.d0, 0.d0, 0.d0, 2*0.d0/
      data (uf(18,ii), ii=1,11) /0.d0, 0.d0, -.707107d0, 0.d0,            ! j=8 F1
     :     0.d0, 0.d0, .707107d0, 0.d0, 0.d0, 2*0.d0/
      data (uf(19,ii), ii=1,11) /0.d0, .707107d0, 0.d0, 0.d0, 0.d0,       ! j=8 F2
     :     0.d0, 0.d0, .707107d0, 0.d0, 2*0.d0/
      data (uf(20,ii), ii=1,11) /-.707107d0, 0.d0, 0.d0, 0.d0, 0.d0,      ! j=8 F1
     :     0.d0, 0.d0, 0.d0, .707107d0, 2*0.d0/
      data (uf(21,ii), ii=1,11) /0.d0, 0.d0, 0.d0, 0.d0, 1.000000d0,      ! j=9 F2
     :     0.d0, 0.d0, 0.d0, 0.d0, 2*0.d0/
      data (uf(22,ii), ii=1,11) /0.d0, 0.d0, 0.d0, -.707107d0, 0.d0,      ! j=9 F1
     :     .707107d0, 0.d0, 0.d0, 0.d0, 2*0.d0/
      data (uf(23,ii), ii=1,11) /0.d0, 0.d0, .707107d0, 0.d0, 0.d0,       ! j=9 F2
     :     0.d0, .707107d0, 0.d0, 0.d0, 2*0.d0/
      data (uf(24,ii), ii=1,11) /0.d0, -.707107d0, 0.d0, 0.d0, 0.d0,      ! j=9 F1
     :    0.d0, 0.d0, .707107d0, 0.d0, 2*0.d0/
      data (uf(25,ii), ii=1,11) /.707107d0, 7*0.d0, .707107d0,            ! j=9 F2
     :    2*0.d0/
      data (uf(26,ii), ii=1,11) /4*0.d0, .707107d0, 0.d0, .707107d0,      ! j=10 F2
     :    4*0.d0/
      data (uf(27,ii), ii=1,11) /3*0.d0, -.707107d0, 3*0.d0,              ! j=10 F1
     :    .707107d0, 3*0.d0/
      data (uf(28,ii), ii=1,11) /0.d0, 0.d0, .707107d0, 5*0.d0,           ! j=10 F2
     :    .7071207d0, 0.d0, 0.d0/
      data (uf(29,ii), ii=1,11) /0.d0, -.707107d0, 7*0.d0,                ! j=10 F1
     :     .707107d0, 0.d0/
      data (uf(30,ii), ii=1,11) /.707107d0, 9*0.d0, .707107d0/            ! j=10 F2
*
*  check for consistency in the values of flaghf and csflag
      if (flaghf) then
        write (6, 5)
        write (9, 5)
5      format (' *** FLAGHF = .TRUE. FOR SINGLET SYSTEM; ABORT ***')
        stop
      end if
      if (flagsu .and. .not. csflag) then
        write (6, 6)
        write (9, 6)
6       format
     :   ('  *** CSFLAG = .FALSE. FOR SURFACE CALCULATION; ABORT ***')
        stop
      end if
      if (flagsu) then
        write (6, 10)
        write (9, 10)
10      format ('  *** FLAGSU = .TRUE. FOR SPHERICAL TOP',
     :     ' COLLISIONS; NOT IMPLEMENTED.  ABORT ***')
        call exit
      end if
      jmx = jmax0
      if (iop .eq. 1) jmx = jmaxa
      if (jmax .gt. jmx) then
        write (6,51) jmax, jmx, iop
        write (9,51) jmax, jmx, iop
51      format(' *** JMAX =',i3,' .GT. THAN ALLOWED IN BASPHTP: ',i3,
     :    ' FOR IOP =',i3,'  ABORT *** ')
        stop
      endif
      xjtot = jtot
      nsum = 0
      do 35  i = 1, nterm
        if (mproj(i) .gt. lammin(i) ) then
          write (6, 30) mproj(i), lammin(i), i
          write (9, 30) mproj(i), lammin(i), i
30        format (' *** MPROJ=',i2,' > LAMMIN=',i2,
     :            ' FOR TERM',i2,'; ABORT ***')
          stop
        end if
        nsum = nsum + (lammax(i) - lammin(i)) + 1
35    continue
      if (nlammx .lt. nsum) then
        write (6, 40) nsum, nlammx
        write (9, 40) nsum, nlammx
40      format (' ** TOTAL NUMBER OF ANISOTROPIC TERMS=', i2,
     :          ' .GT. NLAMMX=', i2,'; ABORT')
        stop
      end if
      if (nsum .ne. nlam) then
        write (6, 45) nsum, nlam
        write (9, 45) nsum, nlam
45      format (' *** NLAM IN BASIS=', i3,' .NE. NLAM FROM SYSDAT=',
     :           i3, '; ABORT ***')
        stop
      end if
      if (bastst) write (6, 46) nsum
      write (9, 46) nsum
46    format (' ** TOTAL NUMBER OF ANISOTROPIC TERMS IN POTENTIAL =',
     :        i3)
      nlam = nsum
      if (iop.lt.1 .or. iop.gt.3) then
        write(6,49) iop
        write(9,49) iop
49      format (' ** IOP =',i3,' OUT OF RANGE')
        stop
      endif
      if (bastst) then
        if (csflag) then
          write (6,75) rmu * xmconv, brot, dj, dk,
     :      iop, ered * econv, jtot, nu
          write (9,75) rmu * xmconv, brot, crot,
     :      iop, ered * econv, jtot, nu
75        format(/,' **  CS SPHERICAL TOP **',
     :      /,'     RMU=', f9.4,'  BROT=',f7.3/
     :      '     DJ=',e13.6,'     DK=', e13.6/
     :      '     IOP=',i2,'  E=', f7.2,'  LBAR=', i5, 2x,' NU=', i3)
        else
          write (6,80) rmu * xmconv, brot, dj, dk,
     :      iop, ered * econv, jtot, jlpar
          write (9,80) rmu * xmconv, brot, dj, dk,
     :      iop, ered * econv, jtot, jlpar
80        format(/,' **  CC SPHERICAL TOP **',
     :      /,'     RMU=', f9.4,'  BROT=',f7.3/
     :      '     DJ=',e13.6,'     DK=', e13.6/
     :      '     IOP=',i2,
     :      '  E=', f7.2, '  JTOT=', i4, 2x,' JLPAR=', i2)
        end if
        if (.not. flagsu) write (9,90) rcut
90      format (/' OPEN CHANNELS ELIMINATED WHICH ARE CLOSED AT R=',
     :          f6.2)
      end if
*
*  first set up list of all rotational levels included in basis
      nlist = 0
      icnt = 0
      do 4000 ji = 0, jmax
*  set up wave functions for A, E, or F levels
        erot = brot * ji * (ji + 1.d0)
        isize = 2 * (ji/2) + 1
*
*  choose nuclear spin modification
        goto (1000,2000,3000), iop
*
*  A levels
1000    if (na(ji+1) .gt. 0) then
          do ia = 1, na(ji+1)
            nlist = nlist + 1
            etemp(nlist) = erot
            fjtemp(nlist) = ji
            fistmp(nlist) = ia
            icnt = icnt + 1
            ieps(nlist) = nepsa(icnt)
            do nn = 1, isize
              isub = (nlist - 1)*narray + nn
              ctemp(isub) = ua(icnt,nn)
            enddo
          enddo
        endif
        goto 4000
*
*  E levels
2000    if (ne(ji+1) .gt. 0) then
          do ie = 1, ne(ji+1)
            nlist = nlist + 1
            etemp(nlist) = erot
            fjtemp(nlist) = ji
            fistmp(nlist) = ie
            icnt = icnt + 1
            ieps(nlist) = nepse(icnt)
            do nn = 1, isize
              isub = (nlist - 1)*narray + nn
              ctemp(isub) = ue(icnt,nn)
            enddo
          enddo
        endif
        goto 4000
*
*  F levels
3000    if (nf(ji+1) .gt. 0) then
          do if = 1, nf(ji+1)
            nlist = nlist + 1
            etemp(nlist) = erot
            fjtemp(nlist) = ji
            fistmp(nlist) = if
            icnt = icnt + 1
            ieps(nlist) = nepsf(icnt)
            do nn = 1, isize
              isub = (nlist - 1)*narray + nn
              ctemp(isub) = uf(icnt,nn)
            enddo
          enddo
        endif
4000  continue
*
*  now sort this list in terms of increasing energy
      if (nlist .gt. 1) then
        do 120 i1 = 1, nlist - 1
          esave = etemp(i1)
          do 115 i2 = i1 + 1, nlist
            if (etemp(i2) .lt. esave) then
*  state i2 has a lower energy than state i1, switch them
              esave = etemp(i2)
              etemp(i2) = etemp(i1)
              etemp(i1) = esave
              fjsave = fjtemp(i2)
              fjtemp(i2) = fjtemp(i1)
              fjtemp(i1) = fjsave
              fissav = fistmp(i2)
              fistmp(i2) = fistmp(i1)
              fistmp(i1) = fissav
              issav = ieps(i2)
              ieps(i2) = ieps(i1)
              ieps(i1) = issav
*  also move e.fn coeffs (don't worry about size of e.fn)
              isize = 2 * (ji/2) + 1
              do mm = 1, isize
                isub1 = (i1 - 1)*narray + mm
                isub2 = (i2 - 1)*narray + mm
                sc1(mm) = ctemp(isub2)
                ctemp(isub2) = ctemp(isub1)
                ctemp(isub1) = sc1(mm)
              end do
            end if
115       continue
120     continue
      end if
*
*  now set up channel and level list
*  print this list if bastst = .true. or if clist = .true.
      if (bastst .or. clist) then
        write (6, 130)
130     format (/,2x,
     :   'LEVEL LIST SORTED BY ENERGY',/,'   N   J  ',
     :     'IS  EPS  EINT(CM-1)  COEFFS')
      end if
      n = 0
      nlevel = 0
      do 170  njk = 1, nlist
        ji = fjtemp(njk)
        isi = fistmp(njk)
*
*  here if this state is to be included
        nlevel = nlevel + 1
        ehold(nlevel) = etemp(njk) / econv
        jhold(nlevel) = ji
        ishold(nlevel) = isi
        fitemp(nlevel) = ieps(njk)
*  also move e.fn coeffs
        isize = 2 * (ji/2) + 1
        do mm = 1, isize
          isub = (nlevel - 1)*narray + mm
          isub1 = (njk - 1)*narray + mm
          chold(isub) = ctemp(isub1)
        end do
*
*  print this level if bastst = .true.
        if (bastst) then
          ecm = ehold(nlevel) * econv
          isize = 2 * jhold(nlevel) + 1
          isub = (nlevel - 1) * narray
          write (6, 1352) nlevel, jhold(nlevel), ishold(nlevel),
     :      int(fitemp(nlevel)), ecm,
     :          (chold((nlevel-1)*narray + mm),mm=1,2*(ji/2)+1)
1352      format (3i4, i5, f10.3,2x,6f8.4/10(29x,6f8.4/))
        endif
*
*  here for cs calculations; include state only if j at least equal to coupled
*  states projection index
        if (csflag) then
          if (ji .ge. nu) then
            n = n + 1
            if (n .gt. nmax) then
              write (6, 150) n, nmax
              write (9, 150) n, nmax
150           format(/' *** NCHANNELS=', i5,
     :             ' .GT. MAX DIMENSION OF',i5,' ABORT ***')
              stop
            end if
            is(n) = ishold(nlevel)
            j(n) = ji
            eint(n) = ehold(nlevel)
            l(n) = jtot
            cent(n) = jtot * (jtot + 1)
*  move e.fn also
            do mm = 1, narray
              isub = (n - 1)*narray + mm
              isub1 = (nlevel - 1)*narray + mm
              c(isub) = chold(isub1)
            end do
          end if
        else if (.not. csflag) then
*
*  here for cc calculations.  first calculate range of orbital angular
*  momentum quantum numbers allowed for this state
*  determine parity of molecular state [Eq. (A21) of S. Green, J. Chem. Phys.
*  73, 2740 (1980) and Townes and Schawlow, Microwave Spectroscopy, Eq. (3-27),
*  p. 64.]  See Eq. (A3) of Green.
*
*  modified determination of parity of rotational wave functions to deal with
*  spherical top.
*
*  first determine whether wave function is composed of eps = +1 or -1 basis fns.
*  note:  only even k basis functions included in spherical top wave functions
          ipar = (-1) ** (ji + 0) * ieps(nlevel)
          lmax = jtot + ji
          lmin = iabs (jtot - ji)
          do 155  li = lmin, lmax
            ix = ipar * (-1) ** (li - jtot)
            if (ix .eq. jlpar) then
              n = n + 1
              if (n .gt. nmax) then
                write (6, 150) n, nmax
                write (9, 150) n, nmax
                stop
              end if
              is(n) = ishold(nlevel)
              j(n) = ji
              eint(n) = ehold(nlevel)
              l(n) = li
              cent(n) = li * (li + 1)
              ips(n) = ieps(nlevel)
*  move e.fn also
             do mm = 1, narray
                isub = (n - 1)*narray + mm
                isub1 = (nlevel - 1)*narray + mm
                c(isub) = chold(isub1)
              end do
            end if
155       continue
        end if
170   continue
*
*  also determine number of levels which are open
      nlevop = 0
      do 250  i = 1, nlevel
        if (ehold(i) .le. ered) then
          nlevop = nlevop + 1
        end if
250   continue
*  now check to see if any of the open channels are closed at r=rcut
*  this is not done for molecule-surface collisions or for rcut < 0
      if (.not.flagsu .and. rcut .gt. 0.d0 .and. .not.boundc) then
        emin = 1. e+7
        do 290  i = 1, n
          if (eint(i) .le. ered) then
*  here if channel is open asymptotically
            if ( jtot * (jtot + 1) / (two * rmu * rcut * rcut)
     :          .gt. (ered - eint(i)) ) then
*  here if channel is open asymptotically but closed at r = rcut
              if (eint(i) .lt. emin) emin = eint(i)
*  emin now contains the lowest channel energy for which this condition is met
            end if
          end if
290     continue
*  now eliminate all channels with eint .ge. emin if any of the channels
*  are open asymptotically but closed at r = rcut
        if (emin.lt.ered) then
          nn = 0
          do 300 i = 1, n
            if (eint(i) .lt. emin) then
*  here if this channel is to be included
              nn = nn + 1
              eint(nn) = eint(i)
              j(nn) = j(i)
              is(nn) = is(i)
              cent(nn) = cent(i)
              l(nn) = l(i)
              ips(nn) = ips(n)
              do mm = 1, narray
                isub = (nn - 1)*narray + mm
                isub1 = (i - 1)*narray + mm
                c(isub) = c(isub1)
              end do
            end if
300       continue
*  reset number of channels
          n = nn
        end if
      end if
*  return if no channels
      if (n .eq. 0) return
      if (nu .eq. numin) then
        ntop = max(n, nlevop)
*  ntop is the maximum row dimension of all matrices passed in the
*  call list of subroutines propag and soutpt.
*  for fps make sure this is an odd number, for faster bank access.
*  this has no effect on vax or cray
        if (mod(ntop,2) .eq. 0 .and. ntop .lt. nmax)
     :       ntop = ntop + 1
      else
        if (n.gt.ntop) then
          write (6, 303) nu, n, ntop
          write (9, 303) nu, n, ntop
303       format (' *** NCH = ',i3, ' AT NU = ',i2, ' .GT. NTOP = ',i3,
     :            '; ABORT **',/,
     :    '     CHECK RCUT')
          call exit
        endif
      end if
*  now list channels if requested
      if (bastst) then
        if (.not.csflag) then
          write (6,305)
          write (9,305)
305       format
     :    (/,2x,'CHANNEL LIST SORTED BY ENERGY',/,
     :     '   N   J  IS  EPS  L  EINT(CM-1)')
        else
          write (6,310) nu
          write (9,310) nu
310       format
     :     (/,2x,'CHANNEL LIST SORTED BY ENERGY',/,
     :     '   N   J  IS   L   EINT(CM-1) ** NU = ',i2)
        end if
        do 330  i = 1, n
          ecm = eint(i) * econv
          if (bastst .or. clist) then
            if (.not.csflag) then
              write (6, 320) i, j(i), is(i), ips(i), l(i), ecm
              write (9, 320) i, j(i), is(i), ips(i), l(i), ecm
320           format (5i4, f10.3)
            else
              write (6, 1320) i, j(i), is(i), l(i), ecm
              write (9, 1320) i, j(i), is(i), l(i), ecm
1320          format (4i4, f10.3)
            endif
          end if
330     continue
      end if
*
*  now calculate coupling matrix elements
*  i counts v2 elements
*  inum counts v2 elements for given lambda
*  ilam counts number of v2 matrices
*  ij is address of given v2 element in present v2 matrix
      i = 0
      if (bastst.and. iprint.ge. 2) then
        write (6, 340)
        write (9, 340)
340     format (/' ILAM  LAMBDA  ICOL  IROW    I    IV2    VEE')
      end if
      lamsum = 0
      ilam = 0
      do 400 iterm = 1, nterm
        lb = lammin(iterm)
*  ilam is the index for the next term in the potential matrix
*  lb is the actual value of lambda
*  mu is not used for Td molecule
          mu = 0
          ilam = ilam + 1
          inum = 0
          ij = 0
          do 355  icol = 1, n
            do 350  irow = icol, n
              ij = ntop * (icol - 1) + irow
              lrow = l(irow)
              if (csflag) lrow = nu
              call vlmats (j(irow), lrow, j(icol), l(icol), jtot,
     :          lb, vee, csflag, irow, icol)
*  check for nonzero v2 matrix element
              if (abs(vee) .gt. 1.d-15) then
                i = i + 1
                if (i .le. nv2max) then
                  inum = inum + 1
                  v2(i) = vee
                  iv2(i) = ij
                  if (bastst.and. iprint.ge.2) then
                    write (6, 345) ilam, lb, icol, irow, i, iv2(i),
     :                             vee
                    write (9, 345) ilam, lb, icol, irow, i, iv2(i),
     :                             vee
345                 format (i4, 2i7, 2i6, i6, g17.8)
                  end if
                end if
              end if
350         continue
355       continue
          if (i .le. nv2max) lamnum(ilam) = inum
          if (bastst) then
            write (6, 370) ilam, lamnum(ilam)
            write (9, 370) ilam, lamnum(ilam)
370         format ('ILAM=',i3,' LAMNUM(ILAM) = ',i6)
          end if
          lamsum = lamsum + lamnum(ilam)
400   continue
      if (i.gt. nv2max) then
         write (6, 450) i, nv2max
         write (9, 450) i, nv2max
450      format (' *** NUMBER OF NONZERO V2 ELEMENTS = ',i6,
     :           ' .GT. NV2MAX=',i6,'; ABORT ***')
         stop
      end if
      if (bastst) then
        write (6, 460) lamsum
        write (9, 460) lamsum
460     format (' ** TOTAL NUMBER OF NONZERO V2 MATRIX ELEMENTS IS ',
     :           i5)
      end if
      return
      end
* --------------------------------------------------------------------
      subroutine vlmats (jp, lp, j, l, jtot, lambda,
     :                   v, csflag, irow, icol)
*  subroutine to calculate v-lambda matrices for close-coupled or coupled-states
*  treatment of collisions of a spherical top with an atom
*  the angular dependence of each lambda term in the tetrahedral expansion
*  potential is given in Table I of Heijmen et al., jcp 107, 902 (1997).
*  the terms are linear combinations of the functions S(lambda,mu), where
*       S(lambda,mu) = [C(lambda,mu) + C(lambda,-mu)]/(1 + delta(mu,0))
*  the C(lambda,mu) are Racah normalized spherical harmonics.
*  the primitive cc and cs matrix elements of normalized spherical
*  harmonics are given in eqs. (26) and (32), respectively,
*  of s. green, j. chem. phys. 64, 3463 (1976)
*  the expressions for the full matrix elements for definite-parity symmetric
*  top wavefunctions are given in eqs. (46-48) and (50) of the same article.
*  note, that in this article the bra indices (left side of matrix elements)
*  are primed, while in the conventions of the present subroutine the bra
*  indices are primed and the ket indices (right-hand side of matrix elements),
*  unprimed.
*
*  author of vlmats:  paul dagdigian
*  revised from vlmstp subr for symmetric top levels by millard alexander
*  appropriate for j/jp <= 10 and lambda <= 7
*
*  current revision date:  18-feb-2016 by pjd
*  -----------------------------------------------------------------------
*  variables in call list:
*    jp:       rotational quantum number of left side of matrix element (bra)
*    lp:       orbital angular momentum of left side of matrix element (bra)
*    j:        rotational quantum number of right side of matrix element (ket)
*    l:        orbital angular momentum of right side of matrix element (ket)
*    jtot:     total angular momentum
*    lambda:   order of legendre term in tetrahedral expansion of potential
*    v:        on return, contains desired matrix element
*    csflag:   if .true., then cs calculation; in this case:
*                jtot in call list is cs lbar
*                lp is nu (cs projection index)
*                l is not used
*              if .false., then cc calculation; in this case:
*                jtot is total angular momentum and l and lp correspond
*                to the orbital angular momenta
*    irow:     number of left-side level in channel list
*    icol:     number of right-side level in channel list
*  variable in common block /coatpi/
*    narray:   maximum size of asymmetric top basis fn expansion
*              set to 200
*  variable in common block /coatpr/
*    c:        expansion coefficients for asymmetric top rotor wave fns.
*  subroutines called:
*    xf3j, xf6j, prmats
*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical csflag
      common /coatpi/ narray
      common /coatpr/ c(1)
      dimension kmin(17)
      data kmin / 0, 0, -2, -2, -4, -4, -6, -6, -8, -8, -10,
     :  -10, -12, -12, -14, -14, 16 /
      sq514 = sqrt(5.d0 / 14.d0)
      sq72 = sqrt(7.d0 / 2.d0)
      sq1113 = sqrt(11.d0 / 13.d0)
*
      v = 0.d0
      isize = 2 * (j/2) + 1
      isizp = 2 * (jp/2) + 1
      do 100 mm = 1,isize
      do 100 nn = 1,isizp
        isub = (icol - 1) * narray + mm
        isubp = (irow - 1) * narray + nn
        if (abs(c(isubp) * c(isub)) .lt. 1.d-8) goto 100
        k = kmin(j+1) + 2 * (mm - 1)
        kp = kmin(jp+1) + 2 * (nn - 1)
        goto (100,100,300,400,100,600,700),lambda
300       call prmats (jp, lp, j, l, jtot, kp, k,
     :          3, 2, vprm32, csflag)
          call prmats (jp, lp, j, l, jtot, kp, k,
     :          3, -2, vpr3m2, csflag)
          vprm = vprm32 + vpr3m2
          goto 900
400       call prmats (jp, lp, j, l, jtot, kp, k,
     :          4, 0, vprm40, csflag)
          call prmats (jp, lp, j, l, jtot, kp, k,
     :          4, 4, vprm44, csflag)
          call prmats (jp, lp, j, l, jtot, kp, k,
     :          4, -4, vpr4m4, csflag)
          vprm = vprm40 - sq514 * (vprm44 + vpr4m4)
          goto 900
600       call prmats (jp, lp, j, l, jtot, kp, k,
     :          6, 0, vprm60, csflag)
          call prmats (jp, lp, j, l, jtot, kp, k,
     :          6, 4, vprm64, csflag)
          call prmats (jp, lp, j, l, jtot, kp, k,
     :          6, -4, vpr6m4, csflag)
          vprm = vprm60 + sq72 * (vprm64 + vpr6m4)
          goto 900
700       call prmats (jp, lp, j, l, jtot, kp, k,
     :          7, 2, vprm72, csflag)
          call prmats (jp, lp, j, l, jtot, kp, k,
     :          7, -2, vpr7m2, csflag)
          call prmats (jp, lp, j, l, jtot, kp, k,
     :          7, 6, vprm76, csflag)
          call prmats (jp, lp, j, l, jtot, kp, k,
     :          7, -6, vpr7m6, csflag)
          vprm = (vprm72 + vpr7m2)
     :      - sq1113 * (vprm76 + vpr7m6)
900     v = v + vprm * c(isubp) * c(isub)
100   continue
      return
      end
* ----------------------------------------------------------------------
      subroutine prmats (jp, lp, j, l, jtot, kp, k, lambda, mu,
     :                   vprm, csflag)
*  subroutine to calculate primitive v-lambda matrix elements for close-coupled
*  treatment of collisions of a symmetric top with an atom
*  the primitive cc and cs matrix elements are given in eqs. (26) and (32),
*  respectively, of s. green, j. chem. phys. 64, 3463 (1976)
*  note, that in this article the bra indices are unprimed and the ket indices
*  primed, while in the conventions of the present subroutine the bra indices
*  are primed and the ket indices, unprimed.
*
*  this subroutine is appropriate for calculation of matrix element of
*  unnormalized spherical harmonic terms C(lambda,mu), used in the angular
*  expansion of a tetrahedral molecule - atom interaction
*  [see Heijmen et al., JCP 107, 902 (1997)].
*
*  duplicate of prmstp
*  author:  millard alexander
*  current revision date:  19-aug-2015 by p.j.dagdigian
*  -----------------------------------------------------------------------
*  variables in call list:
*    jp:       rotational quantum number of left side of matrix element (bra)
*    lp:       orbital angular momentum of left side of matrix element (bra)
*    j:        rotational quantum number of right side of matrix element (ket)
*    l:        orbital angular momentum of right side of matrix element (ket)
*    jtot:     total angular momentum
*    kp:       k quantum number of bra
*    k:        k quantum number of ket
*    lambda:   order of legendre term in expansion of potential
*    mu:       index of legendre term in expansion of potential
*    vrpm:     on return, contains primitive matrix element
*    csflag:   if .true., then cs calculation; in this case:
*                jtot in call list is cs lbar
*                lp is nu (cs projection index)
*                l is not used
*              if .false., then cc calculation; in this case:
*                jtot is total angular momentum
*  subroutines called:
*     xf3j, xf6j
*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical csflag
      data one, two, zero / 1.d0, 2.d0, 0.d0 /
      vprm = zero
      xjp = jp
      xj = j
      xkp = kp
      xk = k
      xjtot = jtot
      if (csflag) then
        nu = lp
        xnu = nu
      end if
      xlp = float(lp)
      xl = float(l)
      xlamda = float(lambda)
      xmu = float(mu)
*  xnorm is appropriate for matrix element of Racah normalized spherical
*  harmonic C(lambda,mu)
      xnorm = (two * xjp + one) * (two * xj + one)
*  the desired matrix element is constructed in x
      if (.not. csflag) then
*  here for cc matrix elements

        x = xf3j (xlp, xl, xlamda, zero, zero, zero)
*  Hutson and Thornley
*        x = xf3j (xlp, xlamda, xl, zero, zero, zero)

        if (x .eq. zero) return

        x = x * xf3j (xjp, xj, xlamda, xkp, - xk, xmu)
*  Hutson and Thornley
*        x = x * xf3j (xjp, xlamda, xj, xkp, xmu, - xk)

        if (x .eq. zero) return

        x = x * xf6j (xj, xl, xjtot, xlp, xjp, xlamda)
*  Hutson and Thornley
*        x = x * xf6j (xl, xlamda, xlp, xj, xjtot, xjp)

        if (x .eq. zero) return

        iphase = jp + j + kp - jtot
*  Hutson and Thornley
*        iphase = lp + l + jtot - kp

        xnorm = xnorm * (two * lp + one) * (two * l + one)
      else if (csflag) then
*  here for cs matrix elements
        iphase = - k - nu
        x = xf3j (xjp, xlamda, xj, -xnu, zero, xnu)
        if (x .eq. zero) return
        x = x * xf3j (xjp, xlamda, xj, -xkp, xmu, xk)
      end if
      vprm = ( (-1.d0) ** iphase) * x * sqrt(xnorm)
      return
      end
* -------------------------------eof------------------------------------
