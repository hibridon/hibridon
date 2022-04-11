!*************************************************************************
!                                                                        *
!                 system dependent routines library                      *
!                                                                        *
!*************************************************************************
!                          routines included:                            *
!  0.  baschk      check that basistype is allowed                       *
!  1.  sysdat      dispatcher to select specific sysdat routine          *
!  2.  syssav      dispatcher to select specific syssav routine          *
!  3.  ptread      dispatcher to select specific ptread routine          *
!*************************************************************************
! -----------------------------------------------------------------------
subroutine baschk(ival)
! -------------------------------------------------------
! subroutine to check that allowed basis is being called for
! this should be updated as new bases are added
! allowed basis types currently are:
!  1:  singlet sigma + atom
!  2:  doublet sigma + atom
!  3:  doublet pi + atom
!  4:  sigma | pi + atom
!  5:  general pi + atom
!  6:  symmetric top + atom (with inversion doubling)
!  7:  1/3 P atom + atom
!  8:  1sigma + 1sigma
!  9:  symmetric top + 1sigma (with inversion doubling)
!  10: 2S atom + 2P atom
!  11: singlet delta + atom
!  12: homonuclear diatomic + 2P atom
!  13: homonuclear diatomic + 3P atom
!  14: doublet delta + atom
!  15: heteronuclear diatomic + 2P atom
!  16: asymmetric top + atom
!  17: CH2(X 3B1) (0,v2,0) bender level
!  18: symmetric top + 1sigma (with no inversion doubling)
!  19: 2sigma | 2pi + atom (no perturbations)
!  20: doublet pi + singlet sigma
!  21: symmetric top + singlet sigma
!  22: 1D/3P atom + atom
!  23: 3P atom + 2S atom
!  24: spherical top + atom
!  25: 1sigma + 1sigma (different molecules)
!  26: 2sigma + 1sigma
!  27: C2v asymmetric top + atom
!  28: 3sigma + 1sigma
!  29. chiral asymmetric top + atom
!  99 or higher: user defined basis
! author:  millard alexander
! revisions  by p. dagdigian
!
! added routine for collision of 1D/3P atom + atom
!   (p.dagdigian, nov-2013)
! added routine for 2pi - 1sigma system (q. ma, jun-2013)
! new 2sigma-2pi routine added by p.dagdigian (dec-2011)
! symmetric top routine for no inversion doubling added
!   by p. dagdigian (mar-2011)
! asymmetric top option added by p. dagdigian (aug-2009)
! CH2(X 3B1) (0,v2,0) optaion added by p. dagdigian (jun-2010)
! changed input variables for bastp basis routine to match
!     this for astp1 (p. dagdigian)
! modified systpln routine (p. dagdigian, aug-2011):
!     changed term from 3 to 12
! changed directory in which file filnam to be read, to:
!     /hibxx/bin/progs/potdata (p.j.dagdigian, 2-jan-2012)
! identify ibasty >99 as user-defined basis
!     (q. ma, 8-oct-2012)
! added 2pi--1sigma basis routine (q. ma, 10-jul-2013)
! added 3P atom + 2S atom basis routine (p.dagdigian, 18-sep-2014)
! added spherical top + atom basis routine (pjd, 21-jul-2105)
! added two unlike 1sigma molecules basis routine (pjd, 23-may-2017)
! added 2sigma + 1sigma basis routine (pjd, 8-jun-2017)
! added C2v asymmetric top basis routine (pjd, sep-2017)
! added 3sigma + 1sigma basis routine (pjd, 30-jun-2018)
! added chiral asymmetric top basis routine (pjd, 16-jan-2019)
!
!     COMMON BLOCK COMXBS DEFINED IN HIMAIN
! -------------------------------------------------------
use mod_comxbs, only: maxbas
common /coselb/ ibasty
logical icheck
icheck=.false.
do 100 i=1, maxbas
  if (ival .eq. i) icheck= .true.
100 continue
if (ival .ge. 99) icheck= .true.
if (.not. icheck) then
    write(6,50) ival
50     format(' *** IBASTY =',i3,' NOT YET INCLUDED,', &
           ' BASIS ROUTINES AVAILABLE:'/, &
         '   1 ->  SINGLET SIGMA + ATOM',/, &
         '   2 ->  DOUBLET SIGMA + ATOM',/, &
         '   3 ->  DOUBLET PI + ATOM',/, &
         '   4 ->  SIGMA | PI + ATOM',/, &
         '   5 ->  GENERAL PI + ATOM',/, &
         '   6 ->  SYMMETRIC TOP + ATOM   ',/, &
         '   7 ->  1/3 P ATOM + ATOM   ',/, &
         '   8 ->  TWO 1SIGMA MOLECULES   ',/, &
         '   9 ->  SYM TOP + 1SIGMA MOLECULE   ',/, &
         '   10 ->  2P ATOM + 2S ATOM   ',/, &
         '   11 ->  SINGLET DELTA + ATOM',/, &
         '   12 ->  HOMONUCLEAR + 2P ATOM',/, &
         '   13 ->  HOMONUCLEAR + 3P ATOM',/, &
         '   14 ->  DOUBLET DELTA	+ ATOM',/, &
         '   15 ->  HETERONUCLEAR + 2P ATOM',/, &
         '   16 ->  ASYMMETRIC TOP + ATOM   ',/, &
         '   17 ->  CH2(X 3B1) (0,V2,0) + ATOM ',/ &
         '   18 ->  SYMMETRIC TOP + ATOM (NO INV)  ',/, &
         '   19 ->  2SIGMA-2PI (NO PERTURBATIONS)  ',/, &
         '   20 ->  DOUBLET PI + SINGLET SIGMA', /, &
         '   21 ->  SYMMETRIC TOP + SINGLET SIGMA', /, &
         '   22 ->  1D/3P ATOM + ATOM  ',/, &
         '   23 ->  3P ATOM + 2S ATOM  ',/, &
         '   24 ->  SPHERICAL TOP + ATOM  '/, &
         '   25 ->  TWO DIFFERENT 1SIGMA MOLECULES   ',/ &
         '   26 ->  2SIGMA - 1SIGMA MOLECULES   ',/ &
         '   27 ->  C2V ASYMMETRIC TOP + ATPM   ',/ &
         '   28 ->  3SIGMA - 1SIGMA MOLECULES   ',/ &
         '   29 ->  CHIRal ASYM TOP + ATOM      ',/ &
         '   99+ -> USER DEFINED BASIS  ',/)
endif
return
end
! -----------------------------------------------------------------------
subroutine sysdat (irpot, readpt, iread)
!   dispatcher to select correct sysdat routine
!   the correct routine is selected according to value of ibasty
!   the following sysdat routines are currently available:
!
!  ibasty         sysdat routine          kind of problem
!    1              sy1sg             singlet sigma scattering
!    2              sy2sg             doublet sigma scattering
!    3              sy2pi             doublet pi scattering
!    4              sysgpi            doublet sigma/pi scattering
!    5              sypi              general pi scattering
!    6              systp             symmetric top scattering - w. inversion doubling
!    7              sy13p             1/3 P atom scattering
!    8              sy2mol            two 1sigma molecules
!    9              systpln           symetric top + 1 sigma molecule
!    10             sy22p             2/2 P atom scattering
!    11             sy1del            singlet delta atom scattering
!    12             syh2p             homonuclear + 2p atom
!    13             syh3p             homonuclear + 3p atom
!    14             syh2del           doublet delta + atom
!    15             sydiat2p          heteronuclear + 2P atom
!    16             syastp            asymmetric top scattering
!    17             sych2x            CH2(X 3B1) (0,v2,0) bender level + atom
!    18             systp1            symmetric top scattering - no inversion doubling
!    19             sysgpi1           2sigma | 2pi + atom (no perturbations)
!    20             sy2pi1sg          doublet pi + singlet sigma
!    21             systp1sg          symmetric top + singlet sigma
!    22             sy1d3p            1D/3P atom + closed-shell atom
!    23             sy3p2s            3P + 2S atom
!    24             sysphtp           sphericAl top + atom
!    25             sav1sg1sg         two different 1sigma molecules
!    26             sav2sg1sg         2sigma + 1sigma molecules
!    27             syastp1           C2v asymmetric top scattering
!    28             sav3sg1sg         3sigma + 1sigma molecules
!    29             syastp2           chiral asymmetric top scattering
!    30.            syastp3           C2v asym top - linear mo;ecule scattering
!    99             syusr             user supplied routine
!
!
!  if iread = 1 read data from input file
!  if iread = 0 return after defining variable names
!  author: b. follmeg
!  current revision date: 30-jul-2018 (p. dagdigian)
!  -----------------------------------------------------------------------
use mod_hiba1sg, only: sy1sg
use mod_hiba03_2pi, only: sy2pi
use mod_hiba12_h2p, only: syh2p
use mod_hiba26_2sg1sg, only: sys2sg1sg
integer ibasty, irpot, iread
logical readpt
common /coselb/ ibasty
#include "common/parbas.F90"
! set default for vibrational quantum numbers to zero for each term
do 10 it=1,maxtrm
ivrow(1,it)=0
ivcol(1,it)=0
10 ntv(it)=1
if (ibasty .ge. 99) then
!  user supplied routine
  call syusr(irpot, readpt, iread)
  return
endif
goto (100,200,300,400,500,600,700,800,900,1000,1100,1200, &
      1300,1400,1500,1600,1700,1800,1900,2000,2100,2200, &
      2300,2400,2500,2600,2700,2800,2900,3000) &
     ibasty
!  singlet sigma variables
100 call sy1sg(irpot, readpt, iread)
return
!  doublet sigma variables
200 call sy2sg(irpot, readpt, iread)
return
!  doublet pi variables
300 call sy2pi(irpot, readpt, iread)
return
!  sigma/pi variables
400 call sysgpi(irpot, readpt, iread)
return
!  general pi variables
500 call sypi(irpot, readpt, iread)
return
!  symmetric top variables - w. inversion doubling
600 call systp(irpot, readpt, iread)
return
!  1/3 P atom variables
700 call sy13p(irpot, readpt, iread)
return
! two 1 sigma molecules
800 call sy2mol(irpot, readpt, iread)
return
! symmetric top + 1 sigma molecule
900 call systpln(irpot, readpt, iread)
return
! 2S atom + 2P atom
1000 call sy22p(irpot, readpt, iread)
return
! singlet delta variables
1100 call sy1del(irpot, readpt, iread)
return
! homonuclear + 2P atom variables
1200 call syh2p(irpot, readpt, iread)
return
! homonuclear + 3P atom variables
1300 call syh3p(irpot, readpt, iread)
return
! doublet delta variable
1400 call sy2del(irpot, readpt, iread)
return
! heteronuclear + 2P atom variables
1500 call sydiat2p(irpot, readpt, iread)
return
! asymmetric top variables
1600 call syastp(irpot, readpt, iread)
return
! CH2(X 3B1) (0,v2,0) bender level variables
1700 call sych2x(irpot, readpt, iread)
return
! symmetric top variables - no inversion doubling
1800  call systp1(irpot, readpt, iread)
return
! 2sigma | 2pi + atom (no perturbations)
1900 call sysgpi1(irpot, readpt, iread)
return
! 2Pi + 1Sigma
2000 call sy2pi1sg(irpot, readpt, iread)
return
! Symmetric top + 1Sigma
2100 call systp1sg(irpot, readpt, iread)
return
! 1D/3P atom + closed-shell atom
2200 call sy1d3p(irpot, readpt, iread)
return
!  3P atom + 2S atom
2300 call sy3p2s(irpot, readpt, iread)
return
!  spherical top + atom
2400 call sysphtp(irpot, readpt, iread)
return
! two different 1sigma molecules
2500 call sy1sg1sg(irpot, readpt, iread)
return
! 2sigma + 1sigma molecules
2600 call sys2sg1sg(irpot, readpt, iread)
return
! C2v asymmetric top variables
2700 call syastp1(irpot, readpt, iread)
return
! 3sigma + 1sigma molecules
2800 call sys3sg1sg(irpot, readpt, iread)
return
! chiral asymmetric top variables
2900 call syastp2(irpot, readpt, iread)
return
! C2v asymmetric top +linear molecule variables
3000 call syastp3(irpot, readpt, iread)
return
end
! -----------------------------------------------------------------------
subroutine syssav (readpt)
!   dispatcher to select correct syssav routine
!   the correct routine is selected according to value of ibasty
!   the following savdat routines are currently available:
!
!  ibasty         syssav routine          kind of problem
!    1              sav1sg            singlet sigma scattering
!    2              sav2sg            doublet sigma scattering
!    3              sav2pi            doublet pi scattering
!    4              savsgpi           sigma/pi scattering
!    5              savpi             general pi scattering
!    6              savstp            symmetric top scattering - w. inversion doubling
!    7              sav13p            1/3 P atom scattering
!    8              sav2mol           1sigma+1sigma
!    9              savstpln          symetric top + 1 sigma molecule
!    10             sav22p            2/2 P atom scattering *
!    11             sav1del           singlet delta scattering
!    12             savh2p            homonuclear + 2p atom
!    13             savh3p            homonuclear + 3p atom
!    14             sav2del           doublet delta + atom
!    15             savdiat2p         heteronuclear +2P atom
!    16             savastp           asymmetric top scattering
!    17             savch2x           CH2(X 3B1) (0,v2,0) bender level + atom
!    18             savstp1           symmetric top scattering - no inversion doubling
!    19             savsgpi1          2sigma | 2pi + atom (no perturbations)
!    20             sav2pi1sg         doublet pi + singlet sigma
!    21             savstp1st         symmetric top + singlet sigma
!    22             sav1d3p           1D/3P atom + closed-shell atom
!    23             sav3p2s           3P atom + 2S atom
!    24             savsphtp          spherical top + atom
!    25             sav1sg1sg         two different 1sigma molecules
!    26             sav2sg1sg         2sigma + 1sigma
!    27             savastp1          C2v asymmetric top scattering
!    28             sav3sg1sg         3sigma + 1sigma
!    29             savastp2          chiral asymmetric top scattering
!    30             savastp3          C2v asymmetric top + linear molecule scattering
!    99             savusr            user supplied routine
!
!  author: b. follmeg
!  current revision date: 20-jun-2019 (p.dagdigian)
!  -----------------------------------------------------------------------
use mod_hiba1sg, only: sav1sg
use mod_hiba03_2pi, only: sav2pi
use mod_hiba12_h2p, only: savh2p
use mod_hiba26_2sg1sg, only: sav2sg1sg

integer ibasty
logical readpt
common /coselb/ ibasty
if (ibasty .ge. 99) then
!  user supplied routine
   call savusr(readpt)
   return
endif
goto (100,200,300,400,500,600,700,800,900,1000,1100,1200, &
      1300,1400,1500,1600,1700,1800,1900,2000,2100,2200, &
      2300,2400,2500,2600,2700,2800,2900,3000) &
     ibasty
!  singlet sigma variables
100 call sav1sg(readpt)
return
!  doublet sigma variables
200 call sav2sg(readpt)
return
!  doublet pi variables
300 call sav2pi(readpt)
return
!  sigma/pi variables
400 call savsgpi(readpt)
return
!  general pi variables
500 call savpi(readpt)
return
!  symmetric top variables - w. inversion doubling
600 call savstp(readpt)
return
!  1/3 P atom variables
700 call sav13p(readpt)
return
!  1sigma+1sigma variables
800 call sav2mol(readpt)
return
!  symmetric top + 1 sigma molecule
!900   call savstpln(irpot, readpt, iread) -- change call (pjd)
900 call savstpln(readpt)
return
!  2/2 P atom variables
1000 call sav22p(readpt)
return
!  singlet delta variables
1100 call sav1del(readpt)
return
!  homonuclear + 2P atom variables
!1200  call savh2p(irpot, readpt, iread) -- change call (pjd)
1200 call savh2p(readpt)
return
!  homonuclear + 3P atom variables
!1300  call savh3p(irpot, readpt, iread) -- change call (pjd)
1300 call savh3p(readpt)
return
!  doublet-delta + atom variables
!1400  call sav2del(irpot, readpt, iread) -- change call (pjd)
1400 call sav2del(readpt)
return
!  heteronuclear + 2P atom variables
!1500  call savdiat2p(irpot, readpt, iread) -- change call (pjd)
1500 call savdiat2p(readpt)
return
!  asymmetric top variables
1600 call savastp(readpt)
return
!  CH2(X 3B1) (0,v2,0) bender level variables
1700 call savch2x(readpt)
return
!  symmetric top variables - w/o. inversion doubling
1800 call savstp1(readpt)
return
!  2sigma | 2pi + atom (no perturbations)
1900 call savsgpi1(readpt)
return
!  2Pi + 1Sigma
2000 call sav2pi1sg(readpt)
return
!  Symmetric top + 1Sigma
2100 call savstp1sg(readpt)
return
!  1D/3P atom + closed-shell atom
2200 call sav1d3p(readpt)
return
!  3P atom + 2S atom
2300 call sav3p2s(readpt)
return
!  spherical top + atom
2400 call savsphtp(readpt)
return
!  two different 1sigma molecules
2500 call sav1sg1sg(readpt)
return
!  2sigma + 1sigma molecules
2600 call sav2sg1sg(readpt)
return
!  C2v asymmetric top variables
2700 call savastp1(readpt)
return
!  3sigma + 1sigma molecules
2800 call sav3sg1sg(readpt)
return
!  chiral asymmetric top variables
2900 call savastp2(readpt)
return
!  C2v asymmetric top + linear molecule variables
3000 call savastp3(readpt)
return
end
! -----------------------------------------------------------------------
subroutine ptread (filnam, readpt)
!   dispatcher to select correct ptread routine
!   the correct routine is selected according to value of ibasty
!   the following ptread routines are currently available:
!
!  ibasty         ptread routine          kind of problem
!    1              ptr1sg            singlet sigma scattering
!    2              ptr2sg            doublet sigma scattering
!    3              ptr2pi            doublet pi scattering
!    4              ptrsgpi           sigma/pi scattering
!    5              ptrpi             general pi scattering
!    6              ptrstp            symmetric top scattering - w. inversion doubling
!    7              ptr13p            1/3 P atom scattering
!    8              ptr2mol           1sigma+1sigma
!    9              ptrstpln          symmetric top + 1 sigma molecule
!    10             ptr22p            2/2 P atom scattering
!    11             ptr1del           singlet delta scattering
!    12             pth2p             homonuclear + 2p atom
!    13             pth3p             homonuclear + 3p atom
!    14             pt2del            doublet delta + atom
!    15             ptrdiat2p         heteronuclear + 2P atom
!    16             ptrastp           asymmetric top + atom
!    17             ptrch2x           CH2(X 3B1) (0,v2,0) bender level + atom
!    18             ptrstp1           symmetric top scattering - no inversion doubling
!    19             ptrsgpi1          2sigma | 2pi + atom (no perturbations)
!    20             ptr2pi1sg         doublet pi + singlet sigma
!    21             ptrstp1sg         symmetric top + singlet sigma
!    22             ptr1d3p           1D/3P atom + closed-shell atom
!    23             ptr3p2s           3P atom + 2S atom
!    24             ptrsphtp          spherical top + atom
!    25             ptr1sg1sg         two different 1sigma molecules
!    26             ptr2sg1sg         2sigma + 1sigma molecules
!    27             ptrastp1          C2v asymmetric top + atom
!    28             ptr3sg1sg         3sigma + 1sigma molecules
!    29             ptrastp2          chiral asymmetric top + atom
!    30             ptrastp3          C2v asymmetric top + linear molecule
!    99             ptrusr            user supplied routine
!
!  author: b. follmeg
!  current revision date: 20-jun-2019 (p.dagdigian)
!  -----------------------------------------------------------------------
use mod_hiba1sg, only: ptr1sg
use mod_hiba03_2pi, only: ptr2pi
use mod_hiba12_h2p, only: ptrh2p
use mod_hiba26_2sg1sg, only: ptr2sg1sg
integer ibasty
logical readpt
character*(*) filnam
common /coselb/ ibasty
if (ibasty .ge. 99) then
!  user supplied routine
   call ptrusr(filnam,readpt)
   return
endif
goto (100,200,300,400,500,600,700,800,900,1000,1100,1200, &
      1300,1400,1500,1600,1700,1800,1900,2000,2100,2200, &
      2300,2400,2500,2600,2700,2800,2900,3000) &
     ibasty
!  singlet sigma potential
100 call ptr1sg(filnam,readpt)
return
!  doublet sigma potential
200 call ptr2sg(filnam,readpt)
return
!  doublet pi potential
300 call ptr2pi(filnam,readpt)
return
!  sigma/pi potential
400 call ptrsgpi(filnam,readpt)
return
!  general pi variables
500 call ptrpi(filnam,readpt)
return
!  symmetric top variables
600 call ptrstp(filnam,readpt)
return
!  1/3 P atom variables
700 call ptr13p(filnam,readpt)
return
!  1sigma+1sigma variables
800 call ptr2mol(filnam,readpt)
return
! symmetric top + 1 sigma molecule
!900   call ptrstpln(irpot, readpt, iread) -- change call (pjd)
900 call ptrstpln(filnam,readpt)
return
!  2/2 P atom variables
1000 call ptr22p(filnam,readpt)
return
!  singlet delta variables
1100 call ptr1del(filnam,readpt)
return
! homonuclear + 2P atom variables
!1200  call ptrh2p(irpot, readpt, iread) -- change call (pjd)
1200 call ptrh2p(filnam,readpt)
return
! homonuclear + 3P atom variables
!1300  call ptrh3p(irpot, readpt, iread) -- change call (pjd)
1300 call ptrh3p(filnam,readpt)
return
! doublet delta + atom variables
!1400  call ptr2del(irpot, readpt, iread) -- change call (pjd)
1400 call ptr2del(filnam,readpt)
return
! heteronuclear + 2P atom variables
!1500  call ptrdiat2p(irpot, readpt, iread) -- change call (pjd)
1500 call ptrdiat2p(filnam,readpt)
return
! asymmetric top variables
1600 call ptrastp(filnam, readpt)
return
! CH2(X 3B1) (0,v2,0) bender level variables
1700 call ptrch2x(filnam, readpt)
return
!  symmetric top variables
1800 call ptrstp1(filnam,readpt)
return
!  2sigma | 2pi + atom (no perturbations) variables
1900 call ptrsgpi1(filnam,readpt)
return
!  2Pi + 1Sigma
2000 call ptr2pi1sg(filnam, readpt)
return
!  Symmetric top + 1Sigma
2100 call ptrstp1sg(filnam, readpt)
return
!  1D/3P atom + closed-shell atom
2200 call ptr1d3p(filnam, readpt)
return
!  3P atom + 2S atom
2300 call ptr3p2s(filnam, readpt)
return
!  spherical top + atom
2400 call ptrsphtp(filnam, readpt)
return
!  two different 1sigma molecules
2500 call ptr1sg1sg(filnam, readpt)
return
!  2sigma + 1sigma molecules
2600 call ptr2sg1sg(filnam, readpt)
return
! C2v asymmetric top variables
2700 call ptrastp1(filnam, readpt)
return
!  3sigma + 1sigma molecules
2800 call ptr3sg1sg(filnam, readpt)
return
! C2v asymmetric top variables
2900 call ptrastp2(filnam, readpt)
return
! C2v asymmetric top + linear molecule variables
3000 call ptrastp3(filnam, readpt)
return
end
! ---------------------------eof--------------------------------
