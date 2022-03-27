!***********************************************************************
!                                                                       *
!                         hibridon 4  library                           *
!                                                                       *
!************************************************************************
!                          routines included:                           *
!                                                                       *
!   2. sprint         prints s-matrices on the screen                   *
!   3. spropn         this subroutine calculates the diagonal matrices  *
!                     to propagate the log-derivative matrix through    *
!                     the current interval                              *
!   4. steppr         determines matrix to transform log-deriv matrix   *
!                     into new interval                                 *
!   5. transp         subroutine to carry out in place transposition    *
!                     of n x n matrix a                                 *
!   6. turn           function, determines classical turning point      *
!   7. wavevc         sets up wavevector matrix and diagonalizes it     *
!   8. xwrite         subroutine to write out integral cross sections   *
!   9. waverd         writes and reads header file for wavefunction     *
!  10. psiasym/psi    to determine wavefunction
!  11. flux           to determine fluxes
!  12. transmt        print out transformation matrix at rout
!************************************************************************
! ----------------------------------------------------------------------
subroutine sprint (fname, ia, basis)
!   reads and prints s-matrices
!   author: h.-j. werner
!   now print out j12 for molecule-molecule systems
!   and for 2P atom + molecule system
!   current revision date: 19-jun-2015 by p.dagdigian
! ----------------------------------------------------------------------
use mod_cosout, only : nnout, jout
use mod_coj12, only: j12
use mod_coj12p, only: j12pk
use constants
use mod_codim, only: nmax => mmax
use mod_cojq, only: jq ! jq(1)
use mod_colq, only: lq ! lq(1)
use mod_coinq, only: inq ! inq(1)
use mod_coinhl, only: jlev => inhold ! jlev(1)
use mod_coisc1, only: inlev => isc1 ! inlev(1)
use mod_coisc2, only: jpack => isc2 ! jpack(1)
use mod_coisc3, only: lpack => isc3 ! lpack(1)
use mod_coisc4, only: ipack => isc4 ! ipack(1)
use mod_cosc1, only: elev => sc1 ! elev(1)
use mod_coz, only: sreal => z_as_vec ! sreal(1)
use mod_cozmat, only: simag => zmat_as_vec ! simag(1)
use mod_hibrid5, only: sread
use mod_basis, only: ab_basis
implicit double precision (a-h,o-z)
class(ab_basis), intent(in) :: basis
character*(*) fname
character*20 cdate
character*40 xname
logical  existf, csflag, flaghf, flagsu, twomol, nucros
dimension ia(4)
#include "common/parpot.F90"
common /coselb/ ibasty

!
!.....jtota: first jtot to be printed
!.....jtotb: last jtot to be printed
!.....jlp:   jlpar to be printed. if zero, both parities are searched
!.....ienerg: number of energy
jtota = ia(1)
jtotb = ia(2)
jlp = ia(3)
ienerg = ia(4)
if(ienerg.eq.0) ienerg = 1
if(jtotb.eq.0) jtotb = jtota
call gennam(xname, fname, ienerg, 'smt', lenx)
inquire (file = xname,  exist = existf)
if (.not. existf) then
  write (6,  20) xname(1:lenx)
20   format(/' FILE ', (a), ' NOT FOUND')
  return
end if
call openf(1, xname, 'tu', 0)
call rdhead(1,cdate, &
   ered,rmu,csflag,flaghf,flagsu,twomol, &
   nucros,jfirst,jfinal,jtotd,numin,numax,nud,nlevel,nlevop,nnout, &
   jlev,inlev,elev,jout)
  if(csflag) jlp=0
  write (6, 70) xname, cdate, label, potnam, ered*econv
70   format(/' S-MATRICES READ FROM FILE ', (a)/ &
          ' WRITTEN: ', (a), /, ' LABEL:   ', (a),/, &
          ' POT NAME:  ', (a) &
         /,' ENERGY: ', f10.3)
  if (.not. twomol) then
    write (6, 80)
80     format ('   N   J   INDEX   EINT(CM-1)')
    do  95  i  =  1,  nlevop
      if (ibasty.ne.12) then
        if (.not. flaghf) then
          write (6,  85) i,  jlev(i),  inlev(i),  elev(i) * econv
85           format (i4,  i5,  i6,  f11.3)
        else
          write (6,  90) i,  (jlev(i)+0.5d0),  inlev(i), &
                           elev(i) * econv
90           format (i4,  f5.1,  i6,  f11.3)
        end if
      else
        write (6,91) i, jlev(i),inlev(i)+0.5,elev(i)*econv
91         format (i4,i5,f6.1,f11.3)
      endif
95     continue

  else
    write (6, 100)
100     format('   N   J1   J2  INDEX  EINT(CM-1)')
      do 135  i  =  1,  nlevop
        jj1 = jlev(i) / 10
        jj2 = mod(jlev(i), 10)
        write (6,  130) i,  jj1,  jj2,  inlev(i), &
                        elev(i) * econv
130         format (i4,  2i5,  i6,  f11.3)
135       continue
  end if
jtotb = min0(jtotb, jfinal)
!
iadr=0
30 nopen = 0
call sread (iadr, sreal, simag, jtot, jlpar, nu, &
                  jq, lq, inq, ipack, jpack, lpack, &
                   1, nmax, nopen, length, ierr, basis)
if(csflag) jlpar=0
if(ierr.eq.-1) goto 400
if(ierr.lt.-1) then
  write(6,35) xname
35   format(' ERROR READING FILE ',(a))
  goto 400
end if
!.....assume that jlpar=1 is stored first
220 if (jlpar .eq.1.and.jlp.eq.-1) goto 30
if(jtot.lt.jtota) goto 30
if(jtot.gt.jtotb) then
   if(jlp.eq.jlpar.or.jlp.eq.-1) goto 400
   jlp=-1
   goto 220
end if
if (.not.flaghf) then
  write (6, 230) jtot, jlpar, nu, nnout
230   format(/' JTOT=', i4, '  JLPAR=', i2, '  NU=', i3, &
          '  NNOUT=', i3)
else
  write (6, 240) jtot+0.5d0, jlpar, nu+0.5d0, nnout
240   format(/' JTOT=', f6.1, '  JLPAR=', i2, '  NU=', f5.1, &
          '  NNOUT=', i3)
end if
if (nnout.le.0) then
  write (6, 250)
250   format(/' COLUMN INDICES:')
  write (6, 290) 'N    ', (j, j=1, nopen)
  if (.not. twomol) then
    if (flaghf) then
      write (6, 260) 'J    ', (jq(j)+0.5d0, j=1, nopen)
260       format (1x, (a), (t10, 20(f6.1)) )
    else
      write (6, 290) 'J    ', (jq(j), j=1, nopen)
    end if
  else
    write (6, 270) 'J1/J2', (jq(j), j=1, nopen)
270     format(1x, (a), (t10, 10i6) )
  end if
  write (6, 290) 'L    ', (lq(j), j=1, nopen)
  write (6, 290) 'INDEX', (inq(j), j=1, nopen)
end if
  write (6, 280)
280   format(/' ROW INDICES:')
  write (6, 290) 'N    ', (j, j=1, length)
  if (.not. basis%uses_j12()) then
    if (flaghf) then
      write (6, 260) 'J    ', (jpack(j)+0.5d0, j=1, length)
    else
      write (6, 290) 'J    ', (jpack(j), j=1, length)
    end if
    write (6, 290) 'IS   ', (ipack(j), j=1, length)
  else
    if (ibasty.eq.12 .or. ibasty.eq.15) then
      write (6, 290) 'J    ', (jpack(j), j=1, length)
      write (6, 260) 'JA   ', (ipack(j)+0.5d0, j=1, length)
      write (6, 260) 'J12  ', (j12pk(j)+0.5d0, j=1, length)
    else
      write (6, 270) 'J1/J2', (jpack(j), j=1, length)
      write (6, 270) 'J12', (j12pk(j), j=1, length)
      write (6, 270) 'IS', (ipack(j), j=1, length)
    endif
  end if
  write (6, 290) 'L    ', (lpack(j), j=1, length)
290   format(1x, (a), (t10, 10i6))
ncol = nopen
if(nnout.gt.0) ncol = length
write (6, 300) 'REAL PART OF THE S-MATRIX'
300 format(/1x, (a))
do 330 ja = 1, length, 10
  je=min0(ja+9,length)
  write (6, 310) (j, j = ja, je)
310   format(10i12)
  ij=1-nmax
  do 320 i=1,ncol
    if(nnout.gt.0) je=min0(ja+9,i)
    if (je.ge.ja) write (6, 315) i, (sreal(ij+j*nmax), j=ja,je)
315     format(1x, i3, 10(1pe12.4))
320   ij=ij+1
330 continue
write (6, 300) 'IMAGINARY PART OF THE S-MATRIX'
do 350 ja = 1, length, 10
  je=min0(ja+9,length)
  write (6, 310) (j, j = ja, je)
  ij=1-nmax
  do 340 i = 1, ncol
    if(nnout.gt.0) je=min0(ja+9,i)
    if (je.ge.ja) write (6, 315) i, (simag(ij+j*nmax), j=ja,je)
340   ij=ij+1
350 continue
goto 30
400 call closf(1)
close(1)
return
end
! -----------------------------------------------------------------------
subroutine spropn (rnow, width, eignow, hp, y1, y4, y2, &
                   gam1, gam2, nch)
!-----------------------------------------------------------------------------
!  this subroutine calculates the diagonal matrices to propagate the
!  log-derivative matrix through the current interval
!  also calculated are the ihomogeneous propagators (explained below)
!-----------------------------------------------------------------------------
!  variables in call list:
!    rnow:       midpoint of the current interval
!    width:      width of the current interval
!    eignow:     array containing the wavevectors
!                these are defined by eq. (6) of m.alexander,
!                j. chem. phys. 81,4510 (1984)
!    hp:         array containing the negative of diagonal elements of the
!                derivative of the wavevector matrix at the center of the
!                current interval [see eq. (9) of m.alexander,
!                j. chem. phys. 81,4510 (1984)
!                this array thus contains the derivative of the diagonal
!                elements of the transformed hamiltonian matrix
!    y1, y2, y4: on entry, contain the desired diagonal elements of the
!                homogeneous propagator
!    gam1, gam2: on return, if photof .true, contain the desired diagonal
!                elements of the ihomogeneous propagators
!                otherwise gam1 and gam2 are returned as zero
!    nch:        the number of channels, this equals the dimensions of the
!                eignow, hp, y1, y2, y2, gam1, and gam2 arrays
!-----------------------------------------------------------------------------
!  the key equations, reproduced below, are taken from
!  m. alexander and d. manolopoulos, "a stable linear reference potential
!  algorithm for solution ..."
!  each uncoupled equation can be written as:
!         2    2
!     [ d / dr + eignow - hp * r ] f(r) = 0
!     where r is the distance from the midpoint of the current interval
!  the linearly indepedent solutions are the airy functions ai(x) and bi(x)
!  where  x = alpha (r + beta)
!                   1/3
!  with   alpha = hp   , and beta = (-eignow) / hp
!  the three diagonal elements of the cauchy propagator necessary to propagate
!  the log-derivative matrix are:
!    b = pi [ ai(x ) bi(x ) - ai(x )bi(x ) ] / alpha
!                 1      2        2     1
!    a = pi [ - ai'(x ) bi(x ) + ai(x ) bi'(x ) ]
!                    1      2        2       1
!    d = pi [ ai(x ) bi'(x ) - ai'(x ) bi(x ) ]
!                 1       2         2      1
!    where x  = alpha ( beta + width / 2) and
!           2
!          x  = alpha ( beta - width / 2)
!           1
!  here "width" denotes the width of the interval
!  the diagonal elements of the "imbedding type" propagator are given in terms
!  of the diagonal elements of the cauchy propagator by:
!     y = a/b     y = y = 1/b    and   y = d/b
!      1           2   3                4
!  for the calculation sof the homogeneous propagators
!  the airy functions are defined in terms of their moduli and phases
!  for negative x these definitions are:
!      ai(-x) = m(x) cos[theta(x)]
!      bi(-x) = m(x) sin[theta(x)]
!      ai'(-x) = n(x) cos[phi(x)]
!      bi'(-x) = n(x) sin[phi(x)]
!  in other words
!          2              2        2
!      m(x)  = sqrt[ ai(x)  + bi(x)  ]
!          2               2         2
!      n(x)  = sqrt[ ai'(x)  + bi'(x)  ]
!      theta(x) = atan [ bi(x) / ai(x) ]
!      phi(x)   = atan [ bi'(x) / ai'(x) ]
!  for positive x the moduli and phases are defined by:
!      ai(x) = m(x) sinh[theta(x)]
!      bi(x) = m(x) cosh[theta(x)]
!      ai'(x) = n(x) sinh[phi(x)]
!      bi'(x) = n(x) cosh[phi(x)]
!  in other words
!          2              2        2
!      m(x)  = sqrt[ bi(x)  - ai(x)  ]
!          2               2         2
!      n(x)  = sqrt[ bi'(x)  - ai'(x)  ]
!      theta(x) = atanh [ ai(x) / bi(x) ]
!      phi(x)   = atanh [ ai'(x) / bi'(x) ]
!  here the the exponentially scaled airy functions
!  ai(x), ai'(x), bi(x), bi'(x) are:
!      ai(x)  = ai(x)  * exp[zeta]
!      ai'(x) = ai'(x) * exp[zeta]
!      bi(x)  = bi(x)  * exp[-zeta]
!      bi'(x) = bi'(x) * exp[-zeta]
!                          3/2
!      where zeta = (2/3) x
!  note that for positive x the phases are labeled chi and eta in
!  m. alexander and d. manolopoulos, "a stable linear reference potential
!  algorithm for solution ..."
!-----------------------------------------------------------------------------
!  for both x  and x  negative
!            1      2
!  (this corresponds to a channel which is classically open at both ends of th
!  interval)
!  we find:
!  y     = 1 / { m  m  sin[theta -theta ] }
!   2             1  2          2      1
!          n  sin[phi -theta ]
!           1        1      2
!  y    = ----------------------
!   1      m  sin[theta - theta ]
!           1          2       1
!          n  sin[phi -theta ]
!           2        2      1
!  y    = ----------------------
!   4      m  sin[theta - theta ]
!           2          2       1
!  here the subscripts 1 and 2 imply the moduli and phases evaluated at x = x
!                                                                            1
!  and x = x  , respectively
!           2
!-----------------------------------------------------------------------------
!  for both x  and x  positive
!            1      2
!  (this corresponds to a channel which is classically closed at both ends of
!  the interval)
!  we find:
!  1 / y  =  m  m  cosh[z -z ] { sinh[theta -theta ]
!       2     1  2       2  1              1      2
!                         + tanh[z -z ] sinh[theta +theta ] }
!                                 2  1            1      2
!                     3/2
!   where z  = (2/3) x     and similarly for z
!          1          1                       2
!          n  { sinh [theta -phi ] - tanh[z -z ] sinh[theta +phi ] }
!           1              2    1          2  1            2    1
!  y    = --------------------------------------------------------
!   1      m  { sinh [theta -theta ] + tanh[z -z ] sinh[theta +theta ] }
!           1              2      1          2  1            2      1
!          n  { sinh [theta -phi ] - tanh[z -z ] sinh[theta +phi ] }
!           2              1    2          2  1            1    2
!  y     = --------------------------------------------------------
!   4      m  { sinh [theta -theta ] + tanh[z -z ] sinh[theta +theta ] }
!           2              2      1          2  1            2      1
!-----------------------------------------------------------------------------
!  for x  positive and x  negative we find:
!       1               2
!  1 / y  = m  m  cosh[z ] cosh[theta ] { - cos[theta ] (1 + tanh[z ])
!       2    1  2       1            1               2             1
!                                   + tanh[theta ] sin[theta ] (1 - tanh[z ])
!                                               1           2             1
!      n { cos[theta ](1 + tanh[z ]) - tanh[phi ] sin[theta ] (1 - tanh[z ]) }
!       1           2            1             1           2             1
! y = ------------------------------------------------------------------------
!  1   m {-cos[theta ](1 + tanh[z ]) + tanh[theta ] sin[theta ] (1 - tanh[z ])
!       1           2            1               1           2             1
!      n {-cos[phi ](1 + tanh[z ]) + tanh[theta ] sin[phi ] (1 - tanh[z ]) }
!       2         2            1               1         2             1
! y  = -----------------------------------------------------------------------
!  4   m {-cos[theta ](1 + tanh[z ]) + tanh[theta ] sin[theta ] (1 - tanh[z ])
!       2           2            1               1           2             1
!-----------------------------------------------------------------------------
!  for x  negative and x  positive we find:
!       1               2
! 1 / y  = m  m  cosh[z ] cosh[theta ] { cos[theta ] (1 + tanh[z ])
!      2    1  2       2            2             1             2
!                                  - tanh[theta ] sin[theta ] (1 - tanh[z ]) }
!                                              2           1             2
!      n {-cos[phi ](1 + tanh[z ]) + tanh[theta ] sin[phi ] (1 - tanh[z ]) }
!       1         1            2               2         1             2
! y  = -----------------------------------------------------------------------
!  1   m {cos[theta ](1 + tanh[z ]) - tanh[theta ] sin[theta ] (1 - tanh[z ])
!       1          1            2               2           1             2
!      n { cos[theta ](1 + tanh[z ]) - tanh[phi ] sin[theta ] (1 - tanh[z ]) }
!       2           1            2             2           1             2
! y  = -----------------------------------------------------------------------
!  4   m {cos[theta ](1 + tanh[z ]) - tanh[theta ] sin[theta ] (1 - tanh[z ])
!       2          1            2               2           1             2
!-----------------------------------------------------------------------------
!  for the special case of a constant reference potential (hp=0)
!  then the propagators are:
!  for eignow .gt. 0 (the classically allowed region)
!    y1 = y4 = k cot (k width)
!    y2 = k / sin (k width)
!    where k = sqrt (eignow)
!  for eignow .lt. 0 (the classically forbidden region)
!    y1 = y4 = kap coth (kap width)
!    y2 = kap / sinh (kap width)
!
!    where kap = sqrt (-eignow)
!-----------------------------------------------------------------------------
!  this subroutine also calculates the diagona linhomogeneous log-derivative
!  propagators.  key equations are 9.13, 9.14, and 9.25 of the ph. d.
!  thesis of d. manolopoulos
!  the diagonal elements of the two inhomogeneous propagators are defined
!  in terms of the linearly intependent solutions psi+ and psi-, which
!  are
!  psi+ = [ - bi(x2)ai(x) + ai(x2)bi(x) ] y2 / w
!  and
!  psi- = [ - bi(x1)ai(x) + ai(x1)bi(x) ] y2 / w
!  where w is the wronskian (1/pi)
!  in the determination of these inhomogeneous propagators
!  the airy functions are defined as follows:
!  for negative x :
!      ai(-x) = ai(x) cos[th] + bi(x) sin[th]
!      bi(-x) = bi(x) cos[th] - ai(x) cos[th]
!                            3/2
!      where zeta = (2/3) |x|
!  for positive x :
!      ai(x) = exp(-zeta) ai(x)
!      bi(x) = exp(zeta) bi(x)

!  the integrals of the airy functions are defined as:

!  for a > 0
!      int[ai(x),{0,a}] = 1/3 - exp(-zeta)iai(a)
!      int[bi(x),{0,a}] = exp(zeta)ibi(a)
!  and for a < 0
!      int[ai(x),{a,0}] = int[ai(-x),{0,-a}]
!                       = 2/3 - q(a) cos(th) + p(a) sin(th)
!      int[bi(x),{a,0}] = int[bi(-x),{0,-b}] =
!                       = p(a) cos(th) + q(a) sin(th)
!  we further assume that the ground state wavefunction times the dipole
!  moment function can be expanded as phi(x) = phi0 + phi1 x
!  where phi1 = d[phi,rmid]/alpha and
!        phi0 = phi(rmid) - beta * d[phi,rmid]
!-----------------------------------------------------------------------------
!  for both x  and x  negative
!            1      2
!  (this corresponds to a channel which is classically open at both ends of th
!  interval)
!  we find:
!  y     = 1 / { m  m  sin[theta -theta ] }
!   2             1  2          2      1
!          n  sin[phi -theta ]
!           1        1      2
!  y    = ----------------------
!   1      m  sin[theta - theta ]
!           1          2       1
!          n  sin[phi -theta ]
!           2        2      1
!  y    = ----------------------
!   4      m  sin[theta - theta ]
!           2          2       1
!  here the subscripts 1 and 2 imply the moduli and phases evaluated at x = x
!                                                                            1
!  and x = x  , respectively
!           2
!-----------------------------------------------------------------------------
!  for both x  and x  positive
!            1      2
!  (this corresponds to a channel which is classically closed at both ends of
!  the interval)
!  we find:
!  1 / y  =  m  m  cosh[z -z ] { sinh[theta -theta ]
!       2     1  2       2  1              1      2
!                         + tanh[z -z ] sinh[theta +theta ] }
!                                 2  1            1      2
!                     3/2
!   where z  = (2/3) x     and similarly for z
!          1          1                       2
!          n  { sinh [theta -phi ] - tanh[z -z ] sinh[theta +phi ] }
!           1              2    1          2  1            2    1
!  y    = --------------------------------------------------------
!   1      m  { sinh [theta -theta ] + tanh[z -z ] sinh[theta +theta ] }
!           1              2      1          2  1            2      1
!          n  { sinh [theta -phi ] - tanh[z -z ] sinh[theta +phi ] }
!           2              1    2          2  1            1    2
!  y     = --------------------------------------------------------
!   4      m  { sinh [theta -theta ] + tanh[z -z ] sinh[theta +theta ] }
!           2              2      1          2  1            2      1
!-----------------------------------------------------------------------------
!  for x  positive and x  negative we find:
!       1               2
!  1 / y  = m  m  cosh[z ] cosh[theta ] { - cos[theta ] (1 + tanh[z ])
!       2    1  2       1            1               2             1
!                                   + tanh[theta ] sin[theta ] (1 - tanh[z ])
!                                               1           2             1
!      n { cos[theta ](1 + tanh[z ]) - tanh[phi ] sin[theta ] (1 - tanh[z ]) }
!       1           2            1             1           2             1
! y = ------------------------------------------------------------------------
!  1   m {-cos[theta ](1 + tanh[z ]) + tanh[theta ] sin[theta ] (1 - tanh[z ])
!       1           2            1               1           2             1
!      n {-cos[phi ](1 + tanh[z ]) + tanh[theta ] sin[phi ] (1 - tanh[z ]) }
!       2         2            1               1         2             1
! y  = -----------------------------------------------------------------------
!  4   m {-cos[theta ](1 + tanh[z ]) + tanh[theta ] sin[theta ] (1 - tanh[z ])
!       2           2            1               1           2             1
!-----------------------------------------------------------------------------
!  for x  negative and x  positive we find:
!       1               2
! 1 / y  = m  m  cosh[z ] cosh[theta ] { cos[theta ] (1 + tanh[z ])
!      2    1  2       2            2             1             2
!                                  - tanh[theta ] sin[theta ] (1 - tanh[z ]) }
!                                              2           1             2
!      n {-cos[phi ](1 + tanh[z ]) + tanh[theta ] sin[phi ] (1 - tanh[z ]) }
!       1         1            2               2         1             2
! y  = -----------------------------------------------------------------------
!  1   m {cos[theta ](1 + tanh[z ]) - tanh[theta ] sin[theta ] (1 - tanh[z ])
!       1          1            2               2           1             2
!      n { cos[theta ](1 + tanh[z ]) - tanh[phi ] sin[theta ] (1 - tanh[z ]) }
!       2           1            2             2           1             2
! y  = -----------------------------------------------------------------------
!  4   m {cos[theta ](1 + tanh[z ]) - tanh[theta ] sin[theta ] (1 - tanh[z ])
!       2          1            2               2           1             2
!-----------------------------------------------------------------------------
!  for the special case of a constant reference potential (hp=0)
!  then the propagators are:
!  for eignow .gt. 0 (the classically allowed region)
!    y1 = y4 = k cot (k width)
!    y2 = k / sin (k width)
!    where k = sqrt (eignow)
!  for eignow .lt. 0 (the classically forbidden region)
!    y1 = y4 = kap coth (kap width)
!    y2 = kap / sinh (kap width)
!
!    where kap = sqrt (-eignow)
!-----------------------------------------------------------------------------
!  written by:  millard alexander
!  current revision date (algorithm):  30-dec-1994
!-----------------------------------------------------------------------------
use mod_coqvec2, only: q => q2

implicit double precision (a-h,o-z)
!      implicit none
double precision a, b, bfact, cs, cs1, cs2, csh, dalph2, dalpha, &
    darg, dbeta, dcay, delzet, denom, dhalf, dkap, dlzeta, &
    dmmod1, dmmod2, dnmod1, dnmod2, doneth, dphi1, dphi2, &
    dpi, droot, dslope, dthet1, dthet2, dtnhfm, &
    dtnhfp, dx1, dx2, dzeta1, dzeta2, emz1, emz2, &
    ez1, ez2, fact, oflow, one, rnow, scai1, scai2, scbi1, scbi2, &
    sn, sn1, sn2, snh, tnhfac, width, x1, x2, xairy1, xairy2, &
    xbiry1, xbiry2, zero
double precision eignow, gam1, gam2, hp, y1, y2, y4
double precision xinpt, fprop
integer i, nch, mxphot, nphoto
logical photof, wavefn, boundf, wrsmat
common /cophot/ photof, wavefn, boundf, wrsmat
dimension eignow(1), hp(1), y1(1), y2(1), y4(1), gam1(1), gam2(1)
data     doneth,        dhalf &
  / 0.333333333333333d0, 0.5d0 /
data zero, one /0.d0, 1.d0/
data  dpi / 3.1415926535897932d0 /
!  the parameter oflow is the largest value of x for which exp(x)
!  does not cause a single precision overflow
!                                     n
!  a reasonable value is x = [ ln(2) 2 ] - 5, where n is the number of bits in
!  the characteristic of a floating point number
data oflow / 83.d0 /
if (.not. photof) then
 call dset(nch,zero,gam1,1)
 call dset(nch,zero,gam2,1)
endif
!     now determine propagators for all nch channels
do 10  i = 1, nch
  dslope = hp(i)
! activate next statement for constant reference potential
! force slope to equal zero, to force constant potential
!        dslope=0.d0
  darg = 1.e+10
  if (dslope .ne. 0.d0) &
    darg = log (abs(eignow(i))) - log (abs(dslope))
  if (darg .gt. 20.d0 .or. width .lt. 1.d-5) then
!  here if the relative slope in the wavevector matrix is less than 1.**(-20)
!  in magnitude, or sector width less than 1.e-5 bohr,
!  in which case the potential is assumed to be constant
    if (eignow(i) .gt. 0) then
!  here for classically allowed region (sines and cosines as reference
!  solutions)
      dcay = sqrt (eignow(i))
      darg = dcay * width
      sn=sin(darg)
      y1(i) = dcay / tan (darg)
      y4(i) = y1(i)
      y2(i) = dcay / sn
!  here for inhomogeneous propagators
      if (photof) then
        cs=cos(darg)
        b=rnow+width*dhalf
        a=rnow-width*dhalf
        denom=dcay*sn
        fact=(one-cs)*(q(i)-q(nch+i)*rnow)
        gam1(i)=(fact+(b-a*cs-sn/dcay)*q(nch+i))/denom
        gam2(i)=(fact+(a-b*cs+sn/dcay)*q(nch+i))/denom
      endif
    else
!  here for classically forbidden region (hyperbolic sines and cosines as
!  reference solutions)
      dkap = sqrt ( - eignow(i))
      darg = dkap * width
      snh=sinh(darg)
      y1(i) = dkap / tanh (darg)
      y4(i) = y1(i)
      y2(i) = dkap / snh
!  here for inhomogeneous propagators
      if (photof) then
        csh=cosh(darg)
        b=rnow+width*dhalf
        a=rnow-width*dhalf
        denom=dkap*snh
        fact=(-one+csh)*(q(i)-q(nch+i)*rnow)
        gam1(i)=(fact+(-b+a*csh+snh/dkap)*q(nch+i))/denom
        gam2(i)=(fact+(-a+b*csh-snh/dkap)*q(nch+i))/denom
      endif
    end if
  else
!  here if the relative slope in the wavevector matrix is greater than
!  1.**(-20) in magnitude, in which case a linear reference potential is used,
!  with airy functions as reference solutions
    droot = ( abs (dslope) ) ** doneth
    dalpha   = sign (droot, dslope)
    dbeta = - eignow(i) / dslope
    dx1 = dalpha * ( dbeta - width * dhalf)
    dx2 = dalpha * ( dbeta + width * dhalf)
    call airymp (dx1, dthet1, dphi1, dmmod1, dnmod1,scai1, scbi1, &
             dzeta1)
    call airymp (dx2, dthet2, dphi2, dmmod2, dnmod2,scai2, scbi2, &
             dzeta2)
    if (photof) then
! determine required airy integrals
      call intairy(dx1, xairy1, xbiry1)
      call intairy(dx2, xairy2, xbiry2)
! convert ground state wavefunction and its derivative from r as
! independent variable to x
      q(i)=q(i)-dbeta*q(nch+i)
      q(nch+i)=q(nch+i)
    endif

    x1 = dx1
    x2 = dx2
!-----------------------------------------------------------------------------
    if (x1 .gt. zero .and. x2 .gt. zero) then
!  here for both x  and x  positive
!                 1      2
      tnhfac = tanh(dzeta2 - dzeta1)
      bfact = sinh(dthet1 - dthet2) + &
              tnhfac * sinh(dthet1 + dthet2)
      dlzeta = dzeta2 - dzeta1
      y2(i) = zero
      if (abs(dlzeta) .le. oflow) then
        b = dmmod1 * dmmod2 * cosh(dzeta2 - dzeta1) * bfact
        y2(i) = 1. / b
      end if
      y1(i) = dnmod1 * (sinh(dthet2 - dphi1) &
            - tnhfac * sinh(dthet2 + dphi1) ) / (dmmod1 * bfact)
      y4(i) = dnmod2 * (sinh(dthet1 - dphi2) &
            + tnhfac * sinh(dthet1 + dphi2) ) / (dmmod2 * bfact)
      if (photof) then
        gam1(i)=-scbi2*xairy2-scai2*xbiry2 &
           +exp(dlzeta)*scbi2*xairy1+exp(-dlzeta)*scai2*xbiry1
        gam2(i)=-scbi1*xairy1-scai1*xbiry1 &
           +exp(dlzeta)*scai1*xbiry2+exp(-dlzeta)*scbi1*xairy2
      endif
!-----------------------------------------------------------------------------
    else if (x1 .le. zero .and. x2 .le. zero) then
!  here for both x  and x  negative
!                 1      2
      b =  dmmod1 * dmmod2 * sin(dthet2 - dthet1)
      y2(i) = 1. / b
      y1(i) = dnmod1 * sin(dphi1 - dthet2) &
            / (dmmod1 * sin(dthet2 - dthet1) )
      y4(i) = dnmod2 * sin(dphi2 - dthet1) &
            / (dmmod2 * sin(dthet2 - dthet1) )
      if (photof) then
        delzet=dzeta2-dzeta1
        cs=cos(delzet)
        sn=sin(delzet)
        gam1(i)=-scai2*xairy2+scbi2*xbiry2 &
               +cs*(scai2*xairy1-scbi2*xbiry1) &
               +sn*(scai2*xbiry1+scbi2*xairy1)
        gam2(i)=-scai1*xairy1+scbi1*xbiry1 &
               +cs*(scai1*xairy2-scbi1*xbiry2) &
               -sn*(scai1*xbiry2+scbi1*xairy2)
      endif
!-----------------------------------------------------------------------------
    else if (x1 .gt. zero .and. x2 .le. zero) then
!  here for x  positive and x  negative
!            1               2
      dtnhfp = 1 + tanh(dzeta1)
      dtnhfm = 1 - tanh(dzeta1)
      bfact = cosh(dthet1) * ( - cos(dthet2) * dtnhfp &
            + tanh(dthet1) * sin(dthet2) * dtnhfm)
      y2(i) = zero
      if (abs(dzeta1) .le. oflow) then
        y2(i) = cosh(dzeta1) * (dmmod1 * dmmod2 * bfact)
        y2(i) = one / y2(i)
      end if
      y1(i) = (dnmod1 * cosh(dphi1) * ( cos(dthet2) * dtnhfp &
            - tanh(dphi1) * sin(dthet2) * dtnhfm) ) &
            / (dmmod1 * bfact)
      y4(i) = (dnmod2 * cosh(dthet1) * ( - cos(dphi2) * dtnhfp &
            + tanh(dthet1) * sin(dphi2) * dtnhfm) ) &
            / (dmmod2 * bfact)
      if (photof) then
        cs2=cos(dzeta2)
        sn2=sin(dzeta2)
        ez1=exp(dzeta1)
        emz1=one/ez1
        gam1(i)=scbi2*(-cs2+xbiry2) &
               +scai2*(sn2-xairy2) &
                +emz1*xairy1*(scbi2*cs2-scai2*sn2) &
                +ez1*xbiry1*(scai2*cs2+scbi2*sn2)
        gam2(i)=-scai1*xbiry1-scbi1*xairy1 &
                +emz1*scai1*(xairy2*cs2-xbiry2*sn2) &
                +ez1*scbi1*(one-xbiry2*cs2-xairy2*sn2)
      endif
!-----------------------------------------------------------------------------
    else if (x2 .gt. zero .and. x1 .le. zero) then
!  here for x  positive and x  negative
!            2               1
      dtnhfp = 1 + tanh(dzeta2)
      dtnhfm = 1 - tanh(dzeta2)
      bfact = cosh(dthet2) * ( cos(dthet1) * dtnhfp &
            - tanh(dthet2) * sin(dthet1) * dtnhfm)
      y2(i) = zero
      if (abs(dzeta2) .le. oflow) then
        y2(i) =  cosh(dzeta2) * (dmmod1 * dmmod2 * bfact)
        y2(i) = one / y2(i)
      end if
      y4(i) = (dnmod2 * cosh(dphi2) * ( cos(dthet1) * dtnhfp &
            - tanh(dphi2) * sin(dthet1) * dtnhfm) ) &
            / (dmmod2 * bfact)
      y1(i) = (dnmod1 * cosh(dthet2) * ( - cos(dphi1) * dtnhfp &
            + tanh(dthet2) * sin(dphi1) * dtnhfm) ) &
            / (dmmod1 * bfact)
      if (photof) then
        ez2=exp(dzeta2)
        emz2=one/ez2
        cs1=cos(dzeta1)
        sn1=sin(dzeta1)
        gam1(i)=-scai2*xbiry2-scbi2*xairy2 &
                +emz2*scai2*(xairy1*cs1-xbiry1*sn1) &
                +ez2*scbi2*(one-xbiry1*cs1-xairy1*sn1)
! bug corrected here 4/14/94
        gam2(i)=scbi1*(-cs1+xbiry1) &
               +scai1*(sn1-xairy1) &
                +emz2*xairy2*(scbi1*cs1-scai1*sn1) &
                +ez2*xbiry2*(scai1*cs1+scbi1*sn1)
! here is the old, incorrect code`
!             gam1(i)=scbi1*(-cs1+xbiry1)
!    :               +scai1*(sn1-xairy1)
!    :                +emz2*xairy2*(scbi1*cs1-scai1*sn1)
!    :                +ez2*xbiry2*(scai1*cs1+scbi1*sn1)
      endif
    end if
!-----------------------------------------------------------------------------
    y1(i) = dalpha * y1(i)
    y4(i) = dalpha * y4(i)
    y2(i) = dalpha * y2(i) / dpi
    if (photof) then
        dalph2=dalpha*dalpha
        gam1(i)=q(i)*gam1(i)*y2(i)*dpi &
                +q(nch+i)*(y1(i)-y2(i))/dalpha
        gam2(i)=q(i)*gam2(i)*y2(i)*dpi &
                +q(nch+i)*(y4(i)-y2(i))/dalpha
        gam1(i)=gam1(i)/dalph2
        gam2(i)=gam2(i)/dalph2
    endif
!  at this point the y1, y2, and y4 propagators correspond identically to
!  eqs. (38)-(44) of m. alexander and d. manolopoulos, "a stable linear
!  reference potential algorithm for solution ..."
  end if
10 continue
return
end
! -----------------------------------------------------------------------
subroutine steppr (vecnow, vecnew, tmat, nmax, n)
!  determine matrix to transform log-deriv matrix into new interval
!  see eq. (22) of m.h. alexander, "hybrid quantum scattering algorithms ..."
! --------------------------------------------------------------------------
!  variables in call list:
!    vecnow:     on entry: matrix of eigenvectors of wavevector matrix in
!                current interval
!                on return: matrix of eigenvectors of wavevector matrix in
!                new interval - this is the matrix tn in eq. (22) of
!                m.h. alexander, "hybrid quantum scattering algorithms ..."
!    vecnew:     on entry:  contains matrix of eigenvectors in next interval
!    tmat:       on return: contains transformation matrix pn in eq. (22)
!    n:          number of channels
!    nmax:       maximum row dimension of matrices
!  subroutines called:
!     rgmmul:    generalized matrix multiply, called here to evaluate
!                a.b-transpose
! --------------------------------------------------------------------------
implicit double precision (a-h,o-z)
!      real vecnow, vecnew, tmat
integer n, nmax, isw
!  matrices of maximum row dimension nmax, stored in packed column form
dimension vecnow(1), vecnew(1), tmat(1)
data isw / 0/
#if defined(HIB_NONE)
call mxma (vecnew, 1, nmax, vecnow, nmax, 1, tmat, 1, nmax, &
            n, n, n)
#endif
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
call dgemm('n','t',n,n,n,1.d0,vecnew,nmax,vecnow,nmax, &
           0d0,tmat,nmax)
#endif
!  restore eigenvectors
call matmov (vecnew, vecnow, n, n, nmax, nmax)
return
end
! -----------------------------------------------------------------------
subroutine transp (a, n, nmax)
!  subroutine to carry out in place transposition of n x n matrix a
!  of maximum row dimension nmax stored in packed column form
!  uses blas routine dswap
!  written by:  millard alexander
!  current revision date: 23-sept-87
implicit double precision (a-h,o-z)
integer icol, icolpt, irowpt, n, nmax, nrow
dimension a(1)
icolpt = 2
irowpt = nmax + 1
do 100 icol = 1, n - 1
!  icolpt points to first sub-diagonal element in column icol
!  irowpt points to first super-diagonal element in row icol
!  nrow is number of subdiagonal elements in this column
  nrow = n - icol
  call dswap (nrow, a(icolpt), 1, a(irowpt), nmax)
  icolpt = icolpt + nmax + 1
  irowpt = irowpt + nmax + 1
100 continue
return
end
! -----------------------------------------------------------------------
function turn(e)
! current revision date: 23-sept-87
use constants
implicit double precision (a-h,o-z)
ee = e/econv
r = 3.0d0
dr = 0.5d0
10 r = r+dr
call pot(vv0,r)
if(vv0-ee) 20,50,30
20 if(dr.lt.0) goto 10
goto 40
30 if(dr.gt.0) goto 10
40 dr = -dr*0.5d0
if(abs(dr).gt.0.01d0) goto 10
50 turn = r
return
end
! -----------------------------------------------------------------------
subroutine wavevc (w, eignow, scr1, scr2, rnow, nch, nmax)
!  this subroutine first sets up the wavevector matrix at rnow
!  then diagonalizes this matrix
!  written by:  millard alexander
!  current revision date: 14-dec-2007
! ----------------------------------------------------------------
!  variables in call list:
!  w:           matrix of maximum row dimension nmax used to store
!               wavevector matrix
!  eignow:      on return:  array containing eigenvalues of wavevector matrix
!  scr1, scr2:  scratch vectors of dimension at least nch
!  rnow:        value of interparticle separation at which wavevector matrix
!               is to be evaluated
!  nch:         number of channels
!  nmax:        maximum number of channels
!  subroutines called:
!     potmat:         determines wavevector matrix
!     tred1,tqlrat:   eispack routines to obtain eigenvalues of real,
!                     matrix
!     dsyevr:         latest lapack eigenvalue routine
!     dscal, dcopy:   linpack blas routines
! ----------------------------------------------------------------
use mod_hibrid3, only: potmat
implicit double precision (a-h,o-z)
!      real rnow, xmin1
!      real eignow, scr1, scr2, w
integer icol, ierr, ipt, nch, nmax, nmaxm1, nmaxp1, nrow
external dscal, dcopy
!     external dscal, dcopy, potmat, tred1, tqlrat
!  square matrix (of row dimension nmax)
dimension w(1)
!  vectors dimensioned at least nch
dimension eignow(1), scr1(1), scr2(1)
!  local arrays (for lapack dsyevr)
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
dimension isuppz(2*nch),iwork(10*nch),work(57*nch)
#endif

! ------------------------------------------------------------------
data xmin1 / -1.d0/
nmaxp1 = nmax + 1
nmaxm1 = nmax - 1
call potmat (w, rnow, nch, nmax)
!  since potmat returns negative of lower triangle of w(r) matrix (eq.(3) of
!  m.h. alexander, "hybrid quantum scattering algorithms ..."),
!  next loop changes its sign
ipt = 1
do 100 icol = 1, nch
!  nrow is the number of (diagonal plus subdiagonal) elements in column icol
!  ipt points to the diagonal element in column icol for a matrix stored in
!  packed column form
  nrow = nch - icol + 1
  call dscal (nrow, xmin1, w(ipt), 1)
  ipt = ipt + nmaxp1
100 continue
!  next loop fills in upper triangle of w
if (nch .gt. 1) then
  ipt = 2
  do 110 icol = 1, nch -1
!  ipt points to the first subdiagonal element in column icol
!  nrow is the number of subdiagonal elements in column icol
    nrow = nch - icol
    call dcopy (nrow, w(ipt), 1, w(ipt + nmaxm1), nmax)
    ipt = ipt + nmaxp1
110   continue
end if
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
lwork=57*nch
liwork=10*nch
abstol=1.e-16
lsup=2*nch
vl = 0.0
vu = 0.0
call dsyevr_wrapper('N','A','L',nch,w,nmax,vl,vu,il,iu,abstol,m, &
   eignow,vecnow,nmax,isuppz,work,lwork,iwork,liwork,ierr)

if (ierr .ne. 0) then
  write (6, 115) ierr
  write (9, 115) ierr
115   format (' *** IERR =',i3,' IN WAVEVC/DSYEVR;  ABORT ***')
  write (9, 120) (eignow (i), i=1, nch)
120   format (' EIGENVALUES ARE:',/,8(1pe16.8) )
  call exit
end if
#endif
#if defined(HIB_UNIX) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
!  transform w to tridiagonal form
!  eignow, scr1 and scr2 are used as scratch vectors here
call tred1 (nmax, nch, w, eignow, scr1, scr2)
!  get eigenvalues of tridiagonal matrix
call tqlrat (nch, eignow, scr2, ierr)
if (ierr .ne. 0) then
  write (9, 130) ierr
  write (6, 130) ierr
130   format &
    (' *** TQLRAT IERR =', i3, ' .N.E. 0 IN WAVEVC; ABORT ***')
  call exit
end if
#endif
return
end
!
!     ------------------------------------------------------------------
function iwavsk(irecr)
!     ------------------------------------------------------------------
!     Function to return offset of wfu file for recrod #irec (stream IO)
!
!     Written by: Qianli Ma
!     Latest revision: 20-apr-2012
!
!     This function needs nchwfu, ipos2 and ipos3 from the cowave common
!     blocks.  These variables are set by waverd.
!
!     The stream IO counterpart for `dbrr(,,,irec)` is `read
!     (,pos=wavesk(irec))...
!     ------------------------------------------------------------------
implicit double precision (a-h, o-z)
integer iwavsk
common /cowave/ irec, ifil, nchwfu, ipos2, ipos3, nrlogd, iendwv, &
     inflev
!     The following variables are for size-determination of (machine
!     dependent) built-in types
integer int_t
double precision dble_t
character char_t
!
if (irecr .le. 0) then
   iwavsk = -1
   goto 100
end if
if (irecr .eq. 1) then
   iwavsk = 1
   goto 100
end if
if (irecr .eq. 2) then
   iwavsk = ipos2
   goto 100
end if
if (irecr .eq. 3) then
   iwavsk = ipos3
   goto 100
end if
!     Length for record 1 (at the beginning of the file)
lr1 = 136 * sizeof(char_t) + (10 + 3 * nchwfu) * sizeof(int_t) &
     + (4 + nchwfu) * sizeof(dble_t)
!     Length for each block written by the Airy and LOGDpropagator
if (inflev .eq. 0) then
   lrairy = (2 * nchwfu ** 2 + 6 * nchwfu + 2) &
        * sizeof(dble_t) + 8 * sizeof(char_t)
   lrlogd = (nchwfu ** 2 + nchwfu + 2) * sizeof(dble_t) &
        + 8 * sizeof(char_t)
else if (inflev .eq. 1) then
   lrairy = (nchwfu + 2) * sizeof(dble_t) + 8 * sizeof(char_t)
   lrlogd = 0
end if
!
if ((irecr - 3) .le. nrlogd) then
!     within the logd range of the file
   iwavsk = lr1 + (irecr - 4) * lrlogd + 1
   goto 100
else
!     airy range of the file
   iwavsk = lr1 + nrlogd * lrlogd &
        + (irecr - 4 - nrlogd) * lrairy + 1
   goto 100
end if
!     should never reach here if function called properly
write (0, *) '*** OOPS! ERROR SEEKING WFU FILE. ABORT.'
call exit()
100 continue
end
!
! -------------------------------------------------------------------------
subroutine wavewr(jtot,jlpar,nu,nch,nchtop,rstart,rendld)
! -------------------------------------------------------
!  subroutine to write initial header information on wavefunction file
!  (file jobname.WFU, logical unit 22), unit is opened in subroutine openfi
!     written by:  millard alexander
!     latest revision:  11-dec-2011
!     common blocks amat and bmat are used to store real and
!     imaginary parts of asymptotic wavefunction (only used in
!     read of wavefunction from saved file)
!     increased nchtop to 1000
!
!     Major revision: 16-mar-2012 by Q. Ma
!     Use stream I/O for smaller file size and better compatibility
!
!     current revision: 20-apr-2012 by q. ma
!
! -------------------------------------------------------
#define AMAT_AS_VEC_METHOD_DISTINCT 1
#define AMAT_AS_VEC_METHOD_POINTER 2
#define AMAT_AS_VEC_METHOD_NOVEC 3
#define AMAT_AS_VEC_METHOD AMAT_AS_VEC_METHOD_DISTINCT
use mod_coeint, only: eint
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_DISTINCT)
use mod_coamat, only: amat => psir ! amat(25) psir(nopen, nopen)
#endif
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_POINTER)
use, intrinsic :: ISO_C_BINDING
use mod_coamat, only: amat ! amat(25)
#endif
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_NOVEC)
use mod_coamat, only: amat ! amat(25)
#endif
use mod_cobmat, only: bmat => psii ! bmat(25), here bmat is used as a vector 
use mod_cotq1, only: dpsir ! dpsir(25)
use mod_cotq2, only: dpsii ! dpsii(25)
use mod_cojq, only: jq ! jq(1)
use mod_colq, only: lq ! lq(1)
use mod_coinq, only: inq ! inq(1)
use mod_coisc1, only: isc1 ! isc1(25)
use mod_cosc1, only: sc1 ! sc1(10)
use mod_cosc2, only: sc2 ! sc2(10)
use mod_cosc3, only: sc3 ! sc3(10)
use mod_cosc4, only: sc4 ! sc4(10)
use mod_cosc5, only: sc5 ! sc5(10)
use mod_cow, only: w => w_as_vec ! w(25)
use mod_cozmat, only: zmat => zmat_as_vec ! zmat(25)

implicit double precision (a-h,o-z)
character*48 oldlab, oldpot
character*20 cdate, olddat
logical         airyfl, airypr, bastst, batch, chlist, &
                csflag, flaghf, flagsu, ihomo, ipos, logdfl, &
                logwr, noprin, partw, readpt, rsflag, swrit, &
                t2test, t2writ, twomol, writs, wrpart, wrxsec, &
                xsecwr, nucros, photof, wavefl
#include "common/parpot.F90"
common /colpar/ airyfl, airypr, bastst, batch, chlist, &
                csflag, flaghf, flagsu, ihomo, ipos, logdfl, &
                logwr, noprin, partw, readpt, rsflag, swrit, &
                t2test, t2writ, twomol, writs, wrpart, wrxsec, &
                xsecwr, nucros, photof, wavefl
common /coered/ ered, rmu
common /cowave/ irec, ifil, nchwfu, ipos2, ipos3, nrlogd, iendwv, &
     inflev
dimension iword(32), word(4)
character csize8(8), csize4(4)
!     The three variables below are used to determine the (machine
!     dependent) size of the built-in types
character char_t
integer int_t
double precision dble_t
!
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_POINTER)
real, pointer :: amat_as_vec(:)
#endif


ifil=22
zero=0.d0
izero=0
!     Size limit of wfu file is not necessary when stream I/O is used.
!     However, if 32-bit integer is used, files exceeding 2GB cannot be
!     seeked properly.  With ifort, '-i8' option will force using 64-bit
!     integers.
!
!$$$      if (nchtop .gt. 1000) then
!$$$         write (6, 60) nchtop
!$$$ 60      format(/' *** NCHTOP=',i3,
!$$$     :        ' EXCEEDS MAX FOR STORAGE OF WAVEFUNCTION; ABORT ***')
!$$$         call exit
!$$$      endif
!
nchwfu = nch
ipos2 = -1
ipos3 = -1
nrlogd = 0
!     Mark the position of the EOF of the WFU file in order to by pass
!     (likely) a bug in the intel compiler that INQUIRE does not return
!     the proper offset
iendwv = 1
!     Write magic number
write (ifil, err=950) char(128), 'WFU'
if (writs) then
   write (ifil, err=950) char(0), char(2), char(0), char(0)
else
   write (ifil, err=950) char(1), char(2), char(0), char(0)
end if
!
write (ifil, err=950) ipos2, ipos3, nrlogd
call dater(cdate)
write (ifil, err=950) cdate, label, potnam
!     Four zero-bytes for alignment / C struct compatibility
write (ifil, err=950) char(0), char(0), char(0), char(0)
!
write (ifil, err=950) jtot, jlpar, nu, nch, csflag, flaghf, photof
write (ifil, err=950) ered, rmu, rstart, rendld
!
write (ifil, err=950) (jq(i), i=1, nch), (lq(i), i=1, nch), &
     (inq(i), i=1, nch)
write (ifil, err=950) (eint(i), i=1, nch)
!
write (ifil, err=950) 'ENDWFUR', char(1)
!
lr1 = 136 * sizeof(char_t) + (10 + 3 * nchwfu) * sizeof(int_t) &
     + (4 + nchwfu) * sizeof(dble_t)
iendwv = iendwv + lr1
irec=3
return
!
!     ------------------------------------------------------------------
!     reads header file for wavefunction (wfu file)
entry waverd(jtot,jlpar,nu,nch,npts,nopen,nphoto,jflux, &
     rstart,rendld,rinf)
!
ifil = 22 ! the wfu file is expected to be open using this unit
!     Read the magic number (from the start of the file)
read (ifil, pos=1, end=900, err=950) csize8
inflev = ichar(csize8(5))
!
read (ifil, end=900, err=950) ipos2, ipos3, nrlogd
!
read (ifil, end=900, err=950) olddat, oldlab, oldpot
label = oldlab
potnam = oldpot
!     Read four zero bytes
read (ifil, end=900, err=950) csize4
!
read (ifil, end=900, err=950) jtot, jlpar, nu, nch, csflag, &
     flaghf, photof
!     nchwfu is used in locating the position for records
nchwfu = nch
read (ifil, end=900, err=950) ered, rmu, rstart, rendld
write (6, 245) olddat
if (jflux .ne. 0) write (3, 245) olddat
if (jflux .eq. 0) write (2, 245) olddat
245 format('    FROM CALCULATION ON: ',(a))
if (jflux .ne. 0) write (3, 250) oldlab
if (jflux .eq. 0) write (2, 250) oldlab
write (6, 250) oldlab
250 format('    INITIAL JOB LABEL: ', (a))
if (jflux .ne. 0) write (3, 251) oldpot
if (jflux .eq. 0) write (2, 251) oldpot
write (6, 251) oldpot
251 format('    INITIAL POT NAME: ', (a))
!
!     Read in channel labels
read (ifil, end=900, err=950) (jq(i), i=1, nch), &
     (lq(i), i=1, nch), (inq(i), i=1, nch), &
     (eint(i), i=1, nch)
!
! start reading in information from record 2 here
read (ifil, end=900, err=950, pos=iwavsk(2)) nrecs, nopen, &
     nphoto, rinf
npts = nrecs - 3
! read in wavevectors, bessel functions j, j', n, n'
! first initialize to zero for all channels
call dset(nch,zero,sc1,1)
call dset(nch,zero,sc2,1)
call dset(nch,zero,sc3,1)
call dset(nch,zero,sc4,1)
call dset(nch,zero,sc5,1)
read (ifil, end=900, err=950) (sc1(i), i=1, nopen), &
     (sc2(i), i=1, nopen), (sc3(i), i=1, nopen), &
     (sc4(i), i=1, nopen), (sc5(i), i=1, nopen)
nopsq = nopen ** 2
! read in sreal and simag, store in w and zmat
read (ifil, end=900, err=950) (w(i), i=1, nopsq), &
     (zmat(i), i=1, nopsq)
if (photof) then
! read in number of initial photodissociation states
!        call dbri(mphoto,1,ifil,izero)
!        nphoto=mphoto
! read in real part of photodissociation amplitude
! overlay sreal which is not needed for photodissociation problem
   read (ifil, end=900, err=950) (w(i), i=1, nphoto * nopen)
! read in imaginary part of photodissociation amplitude
! overlay simag which is not needed for photodissociation problem
   read (ifil, end=900, err=950) (zmat(i), i=1, nphoto * nopen)
endif
! read in channel packing list and real and imaginary parts
! of scattering wavefunction and derivative
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_DISTINCT)
read (ifil, end=900, err=950, pos=iwavsk(3)) &
     (isc1(i), i=1, nopen), (amat(i), i=1, nopsq), &
     (bmat(i), i=1, nopsq), (dpsir(i), i=1, nopsq), &
     (dpsii(i), i=1, nopsq)
#endif
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_POINTER)
! amat_as_vec is a view of the matrix amat(nopen, nopen) as a vector(nopen*nopen)
call C_F_POINTER (C_LOC(amat), amat_as_vec, [nopsq])
read (ifil, end=900, err=950, pos=iwavsk(3)) &
     (isc1(i), i=1, nopen), (amat_as_vec(i), i=1, nopsq), &
     (bmat(i), i=1, nopsq), (dpsir(i), i=1, nopsq), &
     (dpsii(i), i=1, nopsq)
#endif
#if (AMAT_AS_VEC_METHOD == AMAT_AS_VEC_METHOD_NOVEC)
read (ifil, end=900, err=950, pos=iwavsk(3)) &
     (isc1(i), i=1, nopen), amat, &
     (bmat(i), i=1, nopsq), (dpsir(i), i=1, nopsq), &
     (dpsii(i), i=1, nopsq)
#endif
irec = 3
return
!
900 continue
950 write (0, *) '*** ERROR READING/WRITING WFU FILE. ABORT.'
call exit
return
!
end
! ----------------------------------------------------------------------
subroutine psiasy(fj,fn,unit,sr,si,psir,psii,nopen,nmax)
! subroutine to determine real and imaginary part of asymptotic wavefunction o
! derivative of these
!  asmptotically, in the case of inelastic scattering, the wavefunction is
!  exp[-i(kr-l pi/2)] - S exp[i(kr-l pi/2)]
!  whereas in the case of photodissociation,
!  exp[-i(kr-l pi/2)] S - exp[i(kr-l pi/2)]
!  this is equivalent to, in the case of inelastic scattering
!  - yl (1-Sr) + jl Si + i [-jl(1+Sr)+yl Si]
!  and, for photodissociation,
!   yl (1-Sr) + jl Si + i [-jl(1+Sr)-yl Si]
!  written by:  millard alexander
!  current revision date:  16-jun-1990
! ---------------------------------------------------------------------
!  variables in call list:
!    fj             contains (for wavefunction calculation) the normalized
!                   ricatti bessel function jl
!                   contains (for derivative calculation) the derivative with
!                   respect to r of the normalized ricatti bessel function jl
!    fn             contains (for wavefunction calculation) the normalized
!                   ricatti bessel function yl
!                   contains (for derivative calculation) the derivative with
!                   respect to r of the normalized ricatti bessel function yl
!    unit           scratch vector
!    sr, si         matrices of order nmax x nmax which contain
!                   on input: real and imaginary parts of s-matrix
!                   on return: real and imaginary parts of asymptotic
!                   wavefunction
!    psir           on return contains nopen x nopen real part of
!                   asymptotic wavefunction (or derivative)
!    psii           on return contains nopen x nopen imag part of asymptotic
!                   wavefunction (or derivative)
!
!    nopen          number of open channels
!    nmax           row dimension of matrices
!  variables in common block /cophot/
!     photof        true if photodissociation calculation
!                   false if scattering calculation
!     wavefn        true if G(a,b) transformation matrices are saved
!                   to be used later in computing the wavefunction
! ----------------------------------------------------------------------------
implicit double precision (a-h,o-z)
logical photof, wavefn, boundf, wrsmat
common /cophot/ photof, wavefn, boundf, wrsmat
dimension fj(1), fn(1), unit(1), sr(nmax,nmax), si(nmax,nmax), &
          psii(nmax,nmax), psir(nmax,nmax)
one=1.d0
onemin=-1.d0
!   put unit vector into array unit
do 80  icol = 1, nopen
  unit(icol) = one
80 continue
! first we want to calculate real part of wavefunction at infinity
! that is   yl(kr) (Sr-1) + jl(kr) Si for scattering or
!         - yl(kr) (Sr-1) + jl(kr) Si for photodissociation
! first move Sreal into psii
  call matmov (sr, psii, nopen, nopen, nmax, nmax)
! now subtract unit matrix
  call daxpy_wrapper (nopen, onemin, unit, 1, psii(1, 1), nmax + 1)
! now premultiply by diagonal matrix -yl(kr) for photodissociation or
! +yl(kr) for scattering
  do 130 irow = 1, nopen
    fac=one*fn(irow)
    if (photof) fac=-fac
    call dscal(nopen, fac, psii(irow,1), nmax)
130   continue
! now store simag in psir
  call matmov(si, psir, nopen, nopen, nmax, nmax)
! premultiply by diagonal matrix jl(kr)
  do 140 irow = 1, nopen
    call dscal(nopen, fj(irow), psir(irow,1), nmax)
140   continue
! now evaluate +/- yl(kr) (Sr-1) + jl(kR) Si, save in psir
  do 150 icol = 1, nopen
    call daxpy_wrapper(nopen, one, psii(1, icol), 1, psir(1,icol), 1)
150   continue
! psir now contains real part of asymptotic scattering wavefunction
! now compute imaginary part of asymptotic wavefunction
! that is - jl(kr) (1+Sr) + yl(kr) Si for scattering or
!         - jl(kr) (1+Sr) - yl(kr) Si for photodissociation
! now move Sreal into psii
  call matmov (sr, psii, nopen, nopen, nmax, nmax)
! now add unit matrix
  call daxpy_wrapper (nopen, one, unit, 1, psii(1, 1), nmax + 1)
! now premultiply by diagonal matrix -jl(kr)
  do 157 irow = 1, nopen
    fac=-fj(irow)
    call dscal(nopen, fac, psii(irow,1), nmax)
157   continue
! replace real part of s matrix by real part of asymptotic wavefunction
  call matmov(psir,sr,nopen, nopen, nmax, nmax)
! premultiply Simag by diagonal matrix yl(kr) for scattering or by
! -yl(kr) for photodissociation
  do 159 irow = 1, nopen
    fac=fn(irow)
    if (photof) fac=-fac
    call dscal(nopen, fac, si(irow,1), nmax)
159   continue
! now evaluate - jl(kr) (1+Sr) +/- yl(kR) Si, save in psii
  do 161 icol = 1, nopen
    call daxpy_wrapper(nopen, one, si(1, icol), 1, psii(1,icol), 1)
161   continue
! replace imaginary part of s matrix by imaginary part of
! asymptotic wavefunction
  call matmov(psii,si,nopen, nopen, nmax, nmax)
return
end
! ------------------------------------------------------------------
subroutine psi(filnam,a)
!
! driver subroutine to calculate scattering wavefunction and fluxes
! from information stored in direct access file
!
! author: millard alexander
! current revision date (algorithm): 15-apr-1997 by mha
! revised on 30-mar-2012 by q. ma for stream I/O of wfu files
! current revision: 20-apr-2012 by q. ma
!
! special version for 13p collisions
!
! ------------------------------------------------------------------
use mod_cosout, only: nnout, jout
use mod_coiout, only: niout, indout
use constants
use mod_coqvec, only: nphoto
use mod_cocent, only: sc2 => cent
use mod_coeint, only: eint
use mod_coamat, only: psir ! psir(100) psir(nopen,nopen)
use mod_cobmat, only: psii ! psii(100) 
use mod_cotq1, only: dpsir ! dpsir(100)
use mod_cotq2, only: dpsii ! dpsii(100)
use mod_cojq, only: jq ! jq(60)
use mod_colq, only: lq ! lq(10)
use mod_coinq, only: inq ! inq(60)
use mod_coisc1, only: ipack => isc1 ! ipack(10)
use mod_coisc2, only: nlist => isc2 ! nlist(50)
use mod_coisc3, only: nalist => isc3 ! nalist(60)
use mod_coisc5, only: nblist  => isc5   ! nblist(60)
use mod_cosc1, only: pk  => sc1   ! pk(100)
use mod_cosc2, only: fj  => sc2   ! fj(10)
use mod_cosc3, only: fjp => sc3   ! fjp(10)
use mod_cosc4, only: fn  => sc4   ! fn(10)
use mod_cosc5, only: fnp => sc5   ! fnp(10)
use mod_cosc6, only: sc  => sc6   ! sc(100)
use mod_cosc7, only: sc1  => sc7   ! sc1(100)
use mod_cosysi, only: ispar
use mod_coz, only: scmat => z_as_vec ! scmat(100)
use mod_cow, only: sr => w_as_vec ! sr(100)
use mod_cozmat, only: si => zmat_as_vec ! si(100)
use mod_version, only : version
use mod_hibrid3, only: expand

implicit double precision (a-h,o-z)
character*(*) filnam
character*40  psifil, wavfil, amplfil, flxfil
character*20  cdate
character*10  elaps, cpu
character*5   s13p(12)
logical exstfl, batch, lpar(3), photof, wavefn, adiab, &
                ldum, csflag,kill,llpar(19),propf, sumf, &
                coordf
common /cowave/ irec, ifil, nchwfu, ipos2, ipos3, nrlogd, iendwv, &
     inflev
common /colpar/ lpar, batch,ldum,csflag,llpar,photof
common /cotrans/ ttrans(36)
! common for y1, y2, y4
common /coered/ ered, rmu
common /coselb/ ibasty
dimension a(7), sx(3)
data s13p /'3SG0f','3SG1f','3PI0f','3PI1f','3PI2f','1PI1f', &
           '3SG1e','3PI0e','3PI1e','3PI2e','1SG0e','1PI1e'/
!

integer, pointer :: ipol
ipol=>ispar(3)

one=1.d0
onemin=-1.d0
zero=0.d0
! initialize timer
call mtime(cpu0,ela0)
! input
iflux=a(1)
if (iflux .gt. 4 .or. iflux .lt. -3) then
  write (6, 2) iflux
2   format (' *** IFLUX =',i3,' MUST BE -3 ... 4  ***')
  return
endif
iprint=a(2)
inchj = a(5)
inchl = a(6)
inchi = a(7)
thresh=a(3)
factr=a(4)
coordf=.false.
sumf=.false.
adiab=.false.
jflux=iabs(iflux)
if (iflux .eq. -1 .or. iflux .eq. 2) adiab = .true.
if (iflux .eq. -2) then
  sumf=.true.
  jflux=1
  adiab = .false.
endif
if (iflux .eq. 3) then
  jflux=1
  adiab=.false.
  coordf=.true.
  ny=a(2)
  ymin=a(3)
  dy=a(4)
  iprint=a(5)
  thresh=a(6)
  factr=a(7)
endif
if (jflux .eq. 4) then
  rout=a(2)
  adiab=.true.
  sumf=.false.
  coordf=.true.
endif
if (photof) then
   if (thresh .eq. 0.d0) thresh=-1.d9
   if (iprint .eq. 0) iprint= 1
else
   if (factr .eq. 0.d0) factr=1.d0
endif
if (iprint .ne. 0) then
  kill=.false.
else
  kill=.true.
endif
!
! generate filename and check if it is present
ien=0
wavfil=filnam//'.wfu'
call gennam(wavfil,filnam,ien,'wfu',lenfs)
inquire(file = wavfil, exist = exstfl)
if (.not. exstfl) then
    write(6,10) wavfil(1:lenfs)
10     format(' ** WAVEFUNCTION INFORMATION FILE ',(a), &
           ' NOT FOUND **')
    return
end if
! open file which holds transformation data and asymptotic wavefunction
call openf(22, wavfil, 'TU', 0)
call dater(cdate)
! open file to save generated wavefunction
if (jflux .eq. 0) then
  call gennam(psifil,filnam,ien,'psi',lenft)
  call openf(2,psifil,'sf',0)
! write a header
  call version(2)
  write(2,11)
  write(6,11)
11   format(/' ** WAVEFUNCTION DETERMINATION ***',/)
  write (2, 12) wavfil
12   format('    INFORMATION FROM FILE: ',(a))
  write (2,13) cdate
13 format('    THIS CALCULATION ON: ',(a))
endif
if (jflux .ne. 0) then
! open file to save generated flux
  call gennam(flxfil,filnam,ien,'flx',lenft)
  call openf(3,flxfil,'sf',0)
! write a header
  call version(3)
  if (jflux .eq. 1) then
    if (photof) then
      write(3,14)
      write(6,14)
14       format(/, &
   ' ** DETERMINATION OF OUTGOING FLUX ***')
    else
      write(3,15)
      write(6,15)
15       format(/, &
   ' ** DETERMINATION OF INCOMING AND OUTGOING FLUX ***')
    endif
  endif
  if (jflux .eq. 2) then
    write(3,16)
    write(6,16)
16     format(/' ** ADIABATIC ENERGIES ***',/)
  endif
  if (jflux .eq. 4) then
    write(3,17)
    write(6,17)
17     format(/' ** TRANSFORMATION MATRIX **',/)
  endif
  write (3,18) cdate
18   format('    THIS CALCULATION ON: ',(a))
endif
! read header information, s matrix, and asymptotic wavefunction and
! derivative
call waverd(jtot,jlpar,nu,nch,npts,nopen,nphoto, &
            jflux,rstart,rendld,rinf)
if (inflev .ne. 0) then
   write (6, *) '** CALCULATION WITH WRSMAT=.T. REQUIRED.'
   goto 700
end if
if (photof) then
  write (6, 19)
  if (jflux .eq. 0) write (2, 19)
  if (jflux .ne. 0) write (3, 19)
19   format('    PHOTODISSOCIATION BOUNDARY CONDITIONS')
else
  write (6, 20)
  if (jflux .eq. 0) write (2, 20)
  if (jflux .ne. 0) write (3, 20)
20   format('    SCATTERING BOUNDARY CONDITIONS')
  photof=.false.
endif
if (adiab) then
  if (jflux .eq. 0)  write (2,21)
  if (jflux .ne. 0)  write (3,21)
  write (6,21)
21   format ('    ADIABATIC BASIS')
endif
if (.not.adiab .and. .not. sumf) then
  if (.not. coordf) then
    if (ibasty .ne. 7) then
      if (jflux .eq. 0) write (2,22)
      if (jflux .ne. 0)  write (3,22)
      write (6,22)
22       format ('    DIABATIC (ASYMPTOTIC) BASIS')
    else
      if (jflux .eq. 0) write (2,23)
      if (jflux .ne. 0)  write (3,23)
!  print flux even inside of closed region in molecular basis
      kill = .false.
      write (6,23)
23       format ('    MOLECULAR (CASE A) BASIS')
    endif
  else
    if (ny .gt. 0) then
      if (jflux .eq. 0) write (2,24)
      if (jflux .ne. 0) write (3,24)
      write (6,24)
24       format &
       ('    COORDINATE SPACE FLUX; POSITIVE INDEX CHOSEN')
    else
      if (jflux .eq. 0) write (2,25)
      if (jflux .ne. 0) write (3,25)
      write (6,25)
25       format &
       ('    COORDINATE SPACE FLUX; NEGATIVE INDEX CHOSEN')
    endif
  endif
endif
if (sumf) then
      if (jflux .eq. 0) write (2,26)
      if (jflux .ne. 0) write (3,26)
      write (6,26)
26       format &
   ('    DIABATIC (ASYMPTOTIC) BASIS;', &
    ' INELASTIC FLUX SUMMED OVER ROTATIONAL LEVELS')
endif
if (jflux .ne. 0) write (3, 12) wavfil
write (6, 12) wavfil
if (csflag) then
  if (jflux .ne. 0) &
    write(3,27) ered*econv, rmu*xmconv, jtot, nu
  if (jflux .eq.0) &
    write(2,27) ered*econv, rmu*xmconv, jtot, nu
  write(6,27) ered*econv, rmu*xmconv, jtot, nu
27   format('    ENERGY = ',f10.3,' cm(-1);  MASS = ',f9.4, &
     ' amu',/,'    CS CALCULATION:  JTOT = ',i3,'; NU =',i3)
else
  if (jflux .ne. 0) &
     write(3,29) ered*econv, rmu*xmconv, jtot, jlpar
  if (jflux .eq. 0) &
      write(2,29) ered*econv, rmu*xmconv, jtot, jlpar
  write(6,29) ered*econv, rmu*xmconv, jtot, jlpar
29   format('    ENERGY = ',f10.3,' cm(-1);  MASS = ',f9.4, &
     /,'    CC CALCULATION:  JTOT = ',i3,'; JLPAR =',i3)
endif
if (iabs(jflux).eq.1) then
  if (rendld .ge. rinf) then
    write (6, 30) rendld, rinf
30     format (' *** FLUX DESIRED; BUT RENDLD=',f7.3, &
            ' .GE. RINF=',f7.3)
    return
  endif
  write(3, 31) rendld, rinf
  write (6, 31) rendld, rinf
31   format ('    FLUXES DETERMINED FROM R = ', &
        f7.3,' TO R = ',f7.3)
  if (coordf) then
    write (6,34) ymin, dy, ymin+(iabs(ny)-1)*dy
    write (3,34) ymin, dy, ymin+(iabs(ny)-1)*dy
34     format ('                           R-INT = ',f5.2,':', &
          f5.2,':',f5.2)
  endif
  write(3,32) thresh
  write (6,32) thresh
32   format ('    CLOSED CHANNEL THRESHOLD = ',1pg10.3)
  write(3,33) factr
  write (6,33) factr
33   format ('    FACTOR FOR CLOSED CHANNEL DAMP = ',f5.2)
else if(jflux.eq.2) then
  write(3, 35) rendld, rinf
  write (6, 35) rendld, rinf
35   format ('    ADIABATIC ENERGIES DETERMINED FROM R = ', &
        f7.3,' TO R = ',f7.3)
else if(jflux.eq.4) then
  write (6,36) rout
  write (3,36) rout
36   format ( &
 '    DIABATIC->ADIABATIC TRANSFORMATION REQUESTED AT R = ', &
     f7.3)
else if(jflux.eq.0) then
  write(2, 37) rstart, rinf, npts
  write (6, 37) rstart, rinf, npts
37   format ( &
  /,'    WF DEFINED FROM R = ',f7.3,' TO R = ',f7.3,' AT ', &
   i4,' POINTS',/)
endif
inch=0
! check if initial channel is in list of channels
! not for photodissociation
if (.not. photof) then
  if (ibasty .ne. 7) then
    do 40 nn=1, nch
      j1 = jq(nn)
      l1 = lq(nn)
      i1 = inq(nn)



 write(6,443) nn,j1,l1,il,inchj,inchl,inchi
443  format(' ch#',i3,3i5,'  req:',3i5)



      if(j1.eq.inchj.and.l1.eq.inchl.and.i1.eq.inchi) then
        inch=nn
        goto 41
      endif
40     continue
  else
    inch=inchj
  endif
41   if ((jflux .ne. 2 .and. jflux .ne. 4).and. inch .eq. 0) then
    if (jflux.ne.0) write(3, 43) inchj, inchl, inchi
    write (6, 43) inchj, inchl, inchi
43     format( /,' ** CHANNEL (J, l, IN = ',2i3,i4,') NOT IN LIST')
    return
  endif
endif
do 45 i = 1, nch
45   nalist(i)=i
! reorder channels in increasing energy since this is eispack ordering
if (nch .gt. 1) then
  call dcopy(nch,eint,1,sc1,1)
  do 50 i = 1,nch-1
    do 48 j = i+1,nch
      if (sc1(j).lt.sc1(i)) then
!  switch
        ehold=sc1(i)
        sc1(i)=sc1(j)
        sc1(j)=ehold
        inhold=nalist(i)
        nalist(i)=nalist(j)
        nalist(j)=inhold
      endif
48     continue
50   continue
endif
do 51 i = 1, nch
  if (inch .eq. nalist(i)) then
    incha=i
    goto 52
  endif
51 continue
52 continue
! nalist now contains ordering of channels in energy
!   nalist(i) is the number in original list of the channel of ith
!   lowest energy
if (.not.coordf)  then
  if (jflux.ne.0) write(3, 55)
  write (6, 55)
55   format('    CHANNEL PARAMETERS:', &
       '   N   J   L   IN    ENERGY(CM-1)    SQRT(K)')
  if (jflux .ne. 2) then
    inchc=inch
    if (adiab) inchc=incha
    if (.not.photof) then
      if (ibasty .ne. 7 .or. &
          (ibasty .eq. 7 .and. (ipol .eq. 0 .or. adiab &
           .and. jlpar .eq. 1))) then
        if (jflux.eq.0) &
          write(2, 57) inchc, jq(inch), lq(inch), inq(inch), &
            econv*eint(inch)
        if (jflux.ne.0) &
          write(3, 57) inchc, jq(inch), lq(inch), inq(inch), &
            econv*eint(inch)
        write (6, 57) inchc, jq(inch), lq(inch), inq(inch), &
            econv*eint(inch)
57         format(/,15x,'INITIAL:',3i4,i5,f13.3)
      else if (ibasty .eq.7 .and. jlpar .eq. -1 &
              .and. ipol.ne.0) then
        if (jflux.eq.0) &
          write(2, 58) inchc, s13p(inch+6), &
            econv*eint(inch)
        if (jflux.ne.0) &
          write(3, 58) inchc, s13p(inch+6), &
            econv*eint(inch)
        write (6, 58) inchc, s13p(inch+6), &
            econv*eint(inch)
58         format(/,15x,'INITIAL:  ',i2,3x,a5,f18.3)
      endif
    else
! initial channel always 1 for photodissociation, since only one
! column in wavefunction
      inch=1
      incha=1
      inchc=1
    endif
  endif
endif
if (jflux .eq. 4) then
  write (3,55)
  do 60  i=1, nch
    write (3, 59) i, jq(i), lq(i), inq(i), econv*eint(i)
59     format(10x,4i4,f13.3)
60   continue
endif
if (photof) then
  inch=1
  incha=1
  inchc=1
endif
nchsq=nch*nch
nopsq=nopen*nopen
! make a list of pointers
! nlist is pointer to desired probed channels in full channel list
nj = 0
if (ibasty .ne. 7) then
  if (nnout .lt. 0) then
    write (6, 65) nnout
65     format ('  ** WARNING: NNOUT = ',i3, &
          ' .LE. O IN SUBROUTINE PSI')
  else if (nnout .eq. 0) then
    if(jflux.ne. 4) then
      write (6, 68)
68       format ('  ** NO PROBE STATES REQUESTED?  ABORT ***')
      go to 700
    endif
  endif
endif
do 120 i=1, iabs(nnout)
   jo = jout(i)
   if (jflux .eq. 2) then
     if (ibasty .ne. 4) then
       llo = iabs(indout(i))/100
       io = sign(iabs(indout(i))-100*llo,indout(i))
     else
       llo = iabs(indout(i))/1000
       io = sign(iabs(indout(i))-1000*llo,indout(i))
     endif
   endif
   do 100 nn=1, nch
     j1 = jq(nn)
     if(j1.ne.jo) goto 100
       if (jflux .ne. 2) then
         do 75 in = 1, iabs(niout)
           io = indout(in)
           if(inq(nn).ne.io) goto 75
           nj = nj + 1
           nlist(nj) =nn
75          continue
       else
          innq=inq(nn)
          if(innq.ne.io .or. lq(nn) .ne. llo) goto 100
! check to see if state has already been found
           ifound=0
           do 80 if=1,nj
             if (nn .eq. nlist(if)) ifound=1
80            continue
           if (ifound .eq. 0) then
             nj = nj + 1
             nlist(nj) = nn
           endif
       endif
100    continue
120 continue
! check if there had been any match
! if coordinate space flux or 13p scattering, include all states
if (coordf .or. ibasty .eq. 7 .or. sumf) then
  nj=nch
  do 121 i=1, nch
    nlist(i)=i
121   continue
endif
if(nj.eq.0) then
  if(jflux.ne.4) then
    if (jflux .ne. 0) write(3,130)
    if (jflux .eq. 0) write(2,130)
    write(6,130)
    if(.not. batch) write(6,130)
130       format(' *** NO PROBE STATES FOUND, ABORT ***')
    goto 700
  endif
end if
! reorder list in terms of energy
if (adiab) then
  if (nj .gt. 1) then
    do 135 i=1, nj-1
      do 134 j=i+1,nj
        if (eint(nlist(j)) .lt. eint(nlist(i))) then
          nhold=nlist(i)
          nlist(i)=nlist(j)
          nlist(j)=nhold
        endif
134       continue
135     continue
  endif
endif
! establish equivalence with adiabatic level list
do 138 i = 1, nj
  do 136 j =1,nch
    if(nalist(j) .eq. nlist(i)) then
      nblist(i)=j
      goto 138
    endif
136   continue
138 continue
if (.not. coordf) then
  if (.not.sumf) then
    do 142 i=1, nj
      nn=nlist(i)
      nnn=nn
      if (adiab) nnn=nblist(i)
      if (eint(nn) .le. ered) then
        sq=sqrt(2.d0*rmu*(ered-eint(nn)))
      else
        sq=-sqrt(2.d0*rmu*(-ered+eint(nn)))
      endif
      if (ibasty .ne. 7 .or. adiab) then
        if (jflux.eq.0) &
        write(2, 140) nnn, jq(nn), lq(nn), inq(nn), &
                 eint(nn)*econv, sq
        if (jflux.ne.0) write(3,140) nnn,jq(nn),lq(nn),inq(nn), &
                 eint(nn)*econv, sq
        write (6, 140) nnn, jq(nn), lq(nn), inq(nn), &
                 eint(nn)*econv, sq
140         format(16x,'PROBED:',3i4,i5,f13.3,f13.4)
      else
        if (jlpar .eq. -1) nn=nn+6
        if (jflux.eq.0) &
        write(2, 141) nnn, s13p(nn), &
                 eint(nnn)*econv, sq
        if (jflux.ne.0) write(3,141) nnn,s13p(nn), &
                 eint(nnn)*econv, sq
        write (6, 141) nnn, s13p(nn), &
                 eint(nnn)*econv, sq
141         format(16x,'PROBED:  ',i2,3x,a5,f18.3,f13.4)
      endif
142     continue
  else
    do 145 ni = 1, niout
      write (3, 143) indout(ni)
      write (6, 143) indout(ni)
143       format(16x,'PROBED:  INDEX =',i4)
145     continue
  endif
endif
if (jflux.eq.1) then
  if(.not.coordf) then
     write(3, 146)
     write (6, 146)
146      format(17x,'TOTAL:  LAST COLUMN')
  else
     write(3, 147)
     write (6, 147)
147      format('    TOTAL FLUX IN LAST COLUMN')
  endif
endif
! reverse order of adiabatic states, since eispack routines return
! highest energy first
do 148 i=1,nj
  nalist(i)=nch-nblist(i)+1
148 continue
npoint=(inch-1)*nch + 1
! npoint points to top of column inch of full wavefunction matrix
if (ibasty .eq.7 .and. jlpar.eq.-1 .and. ipol.eq.1) then
!         if (adiab) then
!           write (6, 149)
!           if (jflux.ne.0) write(3, 149)
!149        format('    STATE 5 IS PARALLEL POLARIZATION, STATE 6'm
!     :            ' IS PERPENDICULAR')
!         endif
   call waverot(jtot,nch)
 endif
! if 13p  or 2s-2p scattering, determine the matrix for case (a) -> case (e)
#ifdef BASIS_10_IS_INCLUDED
if (ibasty .eq.7) call tcasea(jtot,jlpar)
#endif
if (jflux .eq. 0) then
! here if wavefunction calculation
  if (.not.photof) then
!        if (photof) then
! here for scattering
! now expand psir, psii, dpsir, dpsii into nch x nch matrix, putting zeros
! as closed channel components
    if (nch .gt. nopen) then
      call expand(nopen,nopen,nch,nch,ipack, &
                psir,psii,scmat)
      call expand(nopen,nopen,nch,nch,ipack, &
                dpsir,dpsii,scmat)
    endif
! store desired wavefunction column (real and imaginary part) in 1st and 2nd
! columns of psir
    if (inch .ne. 1) call dcopy(nch,psir(npoint),1,psir,1)
    call dcopy(nch,psii(npoint),1,psir(nch+1),1)
    write(2, 150)
150     format(/' R (BOHR) AND REAL PART OF WAVEFUNCTION', &
          ' (R < 0 INDICATES AIRY PROPAGATION)',/)
    call psicalc(npts,nch,nchsq,nj)
! copy imaginary part of asymptotic wavefunction into first column of psir
! so we can use same loop as above
    call dcopy(nch,psir(nch+1),1,psir,1)
    write(2, 185)
185     format(/' R (BOHR) AND IMAGINARY PART OF WAVEFUNCTION', &
          '(R < 0 INDICATES AIRY PROPAGATION)',/)
    call psicalc(npts,nch,nchsq,nj)
  else
! here for photdissociation, in which case outgoing wavefunction is a
! given column of chi
! npoint points to which of the nphot ground state wavefunctions are to
! be selected
    if (nch .gt. nopen) then
      call expand(nopen,nopen,nch,nch,ipack, &
                  psir,psii,scmat)
      call expand(nopen,nopen,nch,nch,ipack, &
                  dpsir,dpsii,scmat)
    endif
! store desired wavefunction column (real and imaginary part) in 1st and
! 2nd columns of psir
    call dcopy(nch,psir(npoint),1,psir,1)
    call dcopy(nch,psii(npoint),1,psir(nch+1),1)
! store derivatives (real and imaginary)  in 3rd and 4th columns of psir
    call dcopy(nch,dpsir(npoint),1,psir(2*nch+1),1)
    call dcopy(nch,dpsii(npoint),1,psir(3*nch+1),1)
    ipoint=2*nch+1
    irec=npts+4
    write(2, 200)
200     format(/' R (BOHR) AND REAL PART OF CHI')
    iwf = 1
    propf=.true.
    call flux(npts,nch,nchsq,ipoint,nj,adiab,thresh,factr,kill, &
            photof,propf,sumf,inch,iwf,coordf,ny,ymin,dy)
    write(2, 210)
210     format(/' R (BOHR) AND IMAGINARY PART OF CHI')
! reread asymptotic information
    call waverd(jtot,jlpar,nu,nch,npts,nopen,nphoto, &
            jflux,rstart,rendld,rinf)
    if (nch .gt. nopen) then
      call expand(nopen,nopen,nch,nch,ipack, &
                  psir,psii,scmat)
      call expand(nopen,nopen,nch,nch,ipack, &
                  dpsir,dpsii,scmat)
    endif
! store desired wavefunction column (real and imaginary part) in 1st and
! 2nd columns of psir
    call dcopy(nch,psir(npoint),1,psir,1)
    call dcopy(nch,psii(npoint),1,psir(nch+1),1)
! store derivatives (real and imaginary)  in 3rd and 4th columns of psir
    call dcopy(nch,dpsir(npoint),1,psir(2*nch+1),1)
    call dcopy(nch,dpsii(npoint),1,psir(3*nch+1),1)
    iwf = -1
    irec=npts+4
    call flux(npts,nch,nchsq,ipoint,nj,adiab,thresh,factr,kill, &
            photof,propf,sumf,inch,iwf,coordf,ny,ymin,dy)
  endif
else if (jflux .eq. 2) then
  write(3, 300)
300   format(/' R (BOHR) AND ADIABATIC ENERGIES (CM-1)',/)
  irec=npts+4
  call eadiab(npts,nch,nchsq,nj)
else if (jflux .eq. 4) then
  call transmt(npts,nch,nchsq,rout)
else if (jflux .eq. 1) then
! here for flux calculation
! first for total flux (only for scattering)
! now expand psir, psii, dpsir, dpsii into nch x nch matrix, putting zeroz
! as closed channel components
  if (.not. photof) then
    if (nch .gt. nopen) then
      call expand(nopen,nopen,nch,nch,ipack, &
           psir,psii,scmat)
      call expand(nopen,nopen,nch,nch,ipack, &
           dpsir,dpsii,scmat)
    endif
! store desired wavefunction column (real and imaginary part) in 1st and 2nd
! columns of psir
    call dcopy(nch,psir(npoint),1,psir,1)
    call dcopy(nch,psii(npoint),1,psir(nch+1),1)
! store derivatives (real and imaginary)  in 3rd and 4th columns of psir
    call dcopy(nch,dpsir(npoint),1,psir(2*nch+1),1)
    call dcopy(nch,dpsii(npoint),1,psir(3*nch+1),1)
! now compute desired total fluxes
    write(3, 305)
305     format(/' R (BOHR) AND TOTAL FLUXES (UNITS OF HBAR/MU; ', &
            'NO DAMPING)'/)
! first initial flux (not if adiabatic or not if 13p scattering or
! not summed fluxes)
    scsum=0.d0
    if (.not. adiab .and. ibasty .ne. 7 .and. .not.sumf) then
      do 310 i=1, nj
        nni=nlist(i)
        sc(i)=psir(nni)*psir(3*nch+nni)-psir(nch+nni) &
              *psir(2*nch+nni)
      scsum=scsum+sc(i)
310       continue
      write(3, 320) rinf, (sc(i), i=1,nj), scsum
320       format(f10.4,30(1pe11.3))
    endif

    if ((.not. adiab .and. .not. coordf) .and. ibasty .ne.7) then
    scsum=0.d0
    if (.not.sumf) then
      do 321 i=1, nj
        nni=nlist(i)
        sc(i)=psir(nni)*psir(3*nch+nni)-psir(nch+nni) &
            *psir(2*nch+nni)
        scsum=scsum+sc(i)
321       continue
      nout=nj
    else
! here for sum of fluxes over desired indices
      do 322 i=1,niout
        sc(i)=0.d0
322       continue
      scsum=0.d0
      do 330 i=1,nch
        nni=nalist(i)
        scc=psir(nni)*psir(3*nch+nni)-psir(nch+nni) &
            *psir(2*nch+nni)
        do 327 ni=1, niout
          if (inq(nni) .eq. indout(ni)) sc(ni)=sc(ni)+scc
327         continue
330       continue
      scsum=dsum(niout,sc,1)
      nout=niout
    endif
    write(3, 320) rinf, (sc(i), i=1,nout), scsum
    endif
    irec=npts+4
    ipoint=2*nch+1
    iwf = 0
    propf=.true.
! plot out all fluxes for total flux which is numerically well behaved
    tthresh=-1.e9
    call flux(npts,nch,nchsq,ipoint,nj,adiab,thresh,factr,.false., &
            photof,propf,sumf,inch,iwf,coordf,ny,ymin,dy)
  endif
  if (.not. photof) then
! now for incoming flux (only for scattering)
! first construct real and imaginary parts of incoming waves
! here for scattering, in which case incoming flux is a diagonal matrix
    call dset(nchsq,zero,psir,1)
    call dset(nchsq,zero,psii,1)
    call dset(nchsq,zero,dpsir,1)
    call dset(nchsq,zero,dpsii,1)
    ipoint=1
    do 335 i=1, nch
      psir(ipoint)=-fn(i)
      psii(ipoint)=-fj(i)
      dpsir(ipoint)=-fnp(i)
      dpsii(ipoint)=-fjp(i)
      ipoint=ipoint+nch+1
335     continue
    write(3, 340)
340     format(/' R (BOHR) AND INCOMING FLUXES (UNITS OF HBAR/MU)',/)
    if (ibasty .eq.7 .and. jlpar.eq.-1 .and. ipol.eq.1) &
        call waverot(jtot,nch)
! store desired wavefunction column (real and imaginary part) in 1st and 2nd
! columns of psir
    call dcopy(nch,psir(npoint),1,psir,1)
    call dcopy(nch,psii(npoint),1,psir(nch+1),1)
! store derivatives (real and imaginary)  in 3rd and 4th columns of psir
    call dcopy(nch,dpsir(npoint),1,psir(2*nch+1),1)
    call dcopy(nch,dpsii(npoint),1,psir(3*nch+1),1)
! now compute desired incoming fluxes
! first initial flux (not if adiabatic or not if 13p scattering or
! not summed fluxes)
    if (.not. adiab .and. ibasty .ne. 7 ) then
      if (.not.sumf) then
        scsum=0.d0
        do 360 i=1, nj
          nni=nlist(i)
          sc(i)=psir(nni)*psir(3*nch+nni)-psir(nch+nni) &
              *psir(2*nch+nni)
        scsum=scsum+sc(i)
360         continue
        write(3, 320) rinf, (sc(i), i=1,nj), scsum
      else
! here for sum of fluxes over desired indices
        do 362 i=1,niout
          sc(i)=0.d0
362         continue
        scsum=0.d0
        do 365 i=1,nch
          nni=nalist(i)
          scc=psir(nni)*psir(3*nch+nni)-psir(nch+nni) &
            *psir(2*nch+nni)
          do 364 ni=1, niout
            if (inq(nni) .eq. indout(ni)) sc(ni)=sc(ni)+scc
364          continue
365         continue
        scsum=dsum(niout,sc,1)
        nout=niout
      endif
      write(3, 320) rinf, (sc(i), i=1,nout), scsum
    endif
    irec=npts+4
    ipoint=2*nch+1
    iwf = 0
    propf=.false.
    call flux(npts,nch,nchsq,ipoint,nj,adiab,thresh,factr,kill, &
            photof,propf,sumf,inch,iwf,coordf,ny,ymin,dy)
  endif
! now for outgoing flux
  if (.not.photof) then
! outgoing wave: move Sreal into psir and Simag into psii
    call matmov (sr, psir, nopen, nopen, nopen, nopen)
    call matmov (si, psii, nopen, nopen, nopen, nopen)
! premultiply Sr by diagonal matrix yl(kr) and Si by jl(kr)
    do 430 irow = 1, nopen
      fac1=fn(irow)
      fac2=fj(irow)
      call dscal(nopen, fac1, psir(irow), nopen)
      call dscal(nopen, fac2, psii(irow), nopen)
430     continue
! add together, resave in psir, this is real part of outgoing wave
    call daxpy_wrapper(nopsq, one, psii, 1, psir, 1)
! repeat for derivative of part of outgoing wave
    call matmov (sr, dpsir, nopen, nopen, nopen, nopen)
    call matmov (si, dpsii, nopen, nopen, nopen, nopen)
! premultiply Sr by diagonal matrix -ylp(kr) and Si by jlp(kr)
    do 440 irow = 1, nopen
      fac1=fnp(irow)
      fac2=fjp(irow)
      call dscal(nopen, fac1, dpsir(irow), nopen)
      call dscal(nopen, fac2, dpsii(irow), nopen)
440     continue
! add together, resave in dpsir, this is derivative of
! real part of outgoing wave
    call daxpy_wrapper(nopsq, one, dpsii, 1, dpsir, 1)
! repeat for imaginary part of outgoing wave
    call matmov (si, psii, nopen, nopen, nopen, nopen)
    call matmov (sr, scmat, nopen, nopen, nopen, nopen)
! premultiply Sr by diagonal matrix -jl(kr) and Si by yl(kr)
    do 450 irow = 1, nopen
      fac1=fn(irow)
      fac2=-fj(irow)
      call dscal(nopen, fac1, psii(irow), nopen)
      call dscal(nopen, fac2, scmat(irow), nopen)
450     continue
! add together, resave in psii, this is imaginary part of outgoing wave
    call daxpy_wrapper(nopsq, one,scmat, 1, psii, 1)
! repeat for derivative of imaginary part of outgoing wave
    call matmov (si, dpsii, nopen, nopen, nopen, nopen)
    call matmov (sr, scmat, nopen, nopen, nopen, nopen)
! premultiply Sr by diagonal matrix jl(kr) and Si by yl(kr)
    do 460 irow = 1, nopen
      fac1=fnp(irow)
      fac2=-fjp(irow)
      call dscal(nopen, fac1, dpsii(irow), nopen)
      call dscal(nopen, fac2, scmat(irow), nopen)
460     continue
! add together, resave in dpsii, this is imaginary part of derivative
! of outgoing wave
    call daxpy_wrapper(nopsq, one,scmat, 1, dpsii, 1)
! now expand psir, psii, dpsir, dpsii into nch x nch matrix, putting zeroz
! as closed channel components
    if (nch .gt. nopen) then
      call expand(nopen,nopen,nch,nch,ipack, &
                  psir,psii,scmat)
      call expand(nopen,nopen,nch,nch,ipack, &
                  dpsir,dpsii,scmat)
    endif
    if (ibasty .eq.7 .and. jlpar.eq.-1 .and. ipol.eq.1 &
      .and. .not.photof) &
       call waverot(jtot,nch)
  else
! here for photdissociation, in which case outgoing wavefunction is a
! given column of chi
! npoint points to which of the nphot ground state wavefunctions are to
! be selected
    if (nch .gt. nopen) then
      call expand(nopen,nopen,nch,nch,ipack, &
                  psir,psii,scmat)
      call expand(nopen,nopen,nch,nch,ipack, &
                  dpsir,dpsii,scmat)
    endif
  endif
! store desired wavefunction column (real and imaginary part) in 1st and 2nd
! columns of psir
  call dcopy(nch,psir(npoint),1,psir,1)
  call dcopy(nch,psii(npoint),1,psir(nch+1),1)
! store derivatives (real and imaginary)  in 3rd and 4th columns of psir
  call dcopy(nch,dpsir(npoint),1,psir(2*nch+1),1)
  call dcopy(nch,dpsii(npoint),1,psir(3*nch+1),1)


! now compute desired outgoing fluxes
  if (.not. photof) write(3, 550)
550   format(/' R (BOHR) AND OUTGOING FLUXES (UNITS OF HBAR/MU)',/)
  if (photof) write(3, 551)
551   format(/' R (BOHR) AND OUTGOING FLUXES (ATOMIC UNITS)',/)
! first initial flux
  if (coordf) then
    write (3, 554) (ymin+dy*(iy-1), iy=1,iabs(ny))
554     format('    R-INT:',f9.3,30(f11.3))
  endif
  if ((.not. adiab .and. .not. coordf) .and. ibasty .ne.7) then
    scsum=0.d0
    if (.not.sumf) then
      do 560 i=1, nj
       nni=nlist(i)
        sc(i)=psir(nni)*psir(3*nch+nni)-psir(nch+nni) &
            *psir(2*nch+nni)
        scsum=scsum+sc(i)
560       continue
      nout=nj
    else
! here for sum of fluxes over desired indices
        do 575 i=1,niout
          sc(i)=0.d0
575         continue
        scsum=0.d0
        do 580 i=1,nch
          nni=nalist(i)
          scc=psir(nni)*psir(3*nch+nni)-psir(nch+nni) &
            *psir(2*nch+nni)
          do 579 ni=1, niout
            if (inq(nni) .eq. indout(ni)) sc(ni)=sc(ni)+scc
579          continue
580         continue
        scsum=dsum(niout,sc,1)
        nout=niout
    endif
    if (photof) then
      call dscal(nj,1.d0/rmu,sc,1)
      scsum=scsum/rmu
    endif
    write(3, 320) rinf, (sc(i), i=1,nout), scsum
  endif
  ipoint=2*nch+1
  irec=npts+4
  iwf = 0
  if (photof) propf=.true.
  if (.not. photof) propf=.false.
  call flux(npts,nch,nchsq,ipoint,nj,adiab,thresh,factr,kill, &
            photof,propf,sumf,inch,iwf,coordf,ny,ymin,dy)
endif
700 if (photof .or. jflux .eq. 0) close (2)
if (jflux .ne. 0) close (3)
close (22)
call mtime(cpu1,ela1)
cpu1 = cpu1 - cpu0
ela1 = ela1 - ela0
call gettim(cpu1,cpu)
call gettim(ela1,elaps)
if(.not. batch) then
  if (iflux .eq. 0) write(6,720) elaps, cpu
720   format(/,' ** WAVEFUNCTION CALCULATION FINISHED:', &
       /,'    ELAPSED TIME:',(a),'  CPU TIME: ',(a))
  if (iflux .ne. 0) write(6,730) elaps, cpu
730   format(/,' ** FLUX CALCULATION FINISHED:', &
       /,'    ELAPSED TIME:',(a),'  CPU TIME: ',(a))
endif
return
end
! ------------------------------------------------------------------
subroutine flux(npts,nch,nchsq,ipoint,nj,adiab,thresh,factr,kill, &
                photof, propf, sumf,inch,iwf,coordf,nny,ymin,dy)
!
! subroutine to calculate fluxes
!
! author: millard alexander
! current revision date (algorithm): 19-apr-1996 by mha
! revised on 30-mar-2012 by q. ma for stream I/O of wfu files
! current revision: 20-apr-2012 by q. ma
!
! ------------------------------------------------------------------
use mod_coiout, only: niout, indout
use mod_cocent, only: sc2 => cent
use mod_coamat, only: psir ! psir(100) (4,nch)
use mod_cobmat, only: psii ! psii(100) Here psii is used as a vector
use mod_cotq2, only: scmat2 => dpsii ! scmat2(100)
use mod_cotq3, only: scmat3 => scmat ! scmat3(100)
use mod_coinq, only: inq ! inq(60)
use mod_cosc1, only: pk => sc1 ! pk(6)
use mod_coisc2, only: nlist => isc2 ! nlist(60)
use mod_coisc3, only: nalist => isc3 ! nalist(60)
use mod_cosc6, only: sc => sc6 ! sc(60)
use mod_cosc7, only: sc1 => sc7 ! sc1(6)
use mod_cosc8, only: sc8
use mod_cosc9, only: sc9
use mod_coz, only: scmat => z_as_vec ! scmat(100)
use mod_cozmat, only: tcoord => zmat_as_vec ! tcoord(100)
! steve, you may need more space, but i doubt it since tcoord is dimensioned n
implicit double precision (a-h,o-z)
logical adiab, kill, photof, propf, sumf, coordf, ifull
common /cowave/ irec, ifil, nchwfu, ipos2, ipos3, nrlogd, iendwv, &
     inflev

common /coered/ ered, rmu
common /coground/ ifull
common /cotrans/ ttrans(36)
common /coselb/ ibasty
dimension scc(100)
data zero, one, onemin /0.d0, 1.d0, -1.d0/
data ione, mone /1,-1/
! if propf = true then true back-subsititution for flux
! if propf = false then inward propagation
! noffset is start of 5th column of psir
  noffset=4*nch+1
!        open (unit=23, file='propagators.txt',status='unknown')
! determine coordinate matrix if coordf true
  if (coordf) then
    ind=1
    ny=iabs(nny)
    do 50 i= 1, nch
      sc(i)=zero
      if(nny .gt. 0 .and. inq(i).gt. 0) sc(i)=one
      if(nny .lt. 0 .and. inq(i).lt. 0) sc(i)=one
50     continue
! sc is now mask for those states for which index is desired
    do 100 iy = 1, ny
      y=ymin+(iy-1)*dy
      ifull=.false.
      call wfintern(scmat,y,nch,1)
! steve, you'll need to modify wfintern so that scmat returns both function an
! scmat is a vector of length nch containing the nch internal states
! evaluated at internal coordinate y
      call vmul(sc,1,scmat,1,tcoord(ind),1,nch)
! this masks the internal states depending on whether the index is
! positive or negative
      ind=ind+nch
100     continue
  endif
! tcoord now contains as rows internal states and as columns
! values of internal coordinate

! here beings loop over sectors (R), starting from outermost sector
  do 420 kstep=1, npts
    irec=irec-1
    read (ifil, end=900, err=950, pos=iwavsk(irec)) r
! branch out if we've gone beyond airy propagation region
    if (r .gt. 0) goto 450
    read (ifil, end=900, err=950) drnow
! read in local wavevectors
    read (ifil, end=900, err=950) (sc8(i), i=1, nch)
! read in first G(n,n+1) matrix,  then Tn transformation
! matrix into local interval
    read (ifil, end=900, err=950) (scmat3(i), i=1, nchsq), &
         (scmat(i), i=1, nchsq)
! read in propagators (y1=pk, y2=sc1, y4=sc2, gam1=sc9,
!     muab= 5th column of psi)
    read (ifil, end=900, err=950) (pk(i), i=1, nch), &
         (sc1(i), i=1, nch), (sc2(i), i=1, nch), &
         (sc9(i), i=1, nch), (psir(noffset - 1 + i), i=1, nch)
!          write (23, 299) -r, drnow, (scmat(ii), ii=1,4),
!     :      (pk(ii),ii=1,2),(sc1(ii),ii=1,2),
!     :      (sc2(ii),ii=1,2),(sc9(ii),ii=1,2),(psir(noffset+ii),ii=0,1)

299     format (2f16.12,14(1pe22.12e3))
! transform wave function into local basis
    if (propf) then
! scmat2(nch, 2) = scmat(nch, nch) * psir(nch, 2)
      call mxma(scmat,1,nch,psir,1,nch,scmat2,1,nch,nch,nch,2)
      call dcopy(2*nch,scmat2,1,psii(ipoint),1)
! here for back substitution for wavefunction
! 3rd and 4th columns of psii contain real and imaginary parts of function
! propagate functions
! evaluate G-tilde* psi-tilde(b)
       call mxma(scmat3,1,nch,psii(ipoint),1,nch,psii,1, &
               nch,nch,nch,2)
! 1st two column of psii now contain G-tilde(A-B)*psi-tilde(b)
! if photodissociation, subtract off mu-tilde(a,b) from real part
      if (photof) call vadd(mone,psii,1 &
                          ,psir(noffset),1,nch)
! at this point 1st two columns of psii contain real and imaginary
!   psi-tilde(a)
! 3rd and 4th columns still contain psi-tilde(b)
! propagate derivatives
      call vmul(pk,1,psii,1,scmat2,1,nch)
      call vmul(pk,1,psii(1+nch),1,scmat2(1+nch),1,nch)
! first and second columns of scmat2 now contain y1 Fa
      call vmul(sc1,1,psii(ipoint),1,scmat2(ipoint),1,nch)
      call vmul(sc1,1,psii(ipoint+nch),1, &
                scmat2(ipoint+nch),1,nch)
! 3rd and 4th columns of scmat2 now contain y2 Fb
      call daxpy_wrapper(2*nch,onemin,scmat2,1,scmat2(ipoint),1)
! 3rd and 4th columns of scmat2 now contains (-y1 Fa + y2 Fb)
! this is psip-tilde(a)  move to 2nd and 3rd columns of psir
      call dcopy(2*nch,scmat2(ipoint),1,psir(ipoint),1)
! if photodissociation subtract off gamma1 from real part of derivative
      if (photof) call vadd(mone,psir(ipoint),1,sc9,1,nch)
! for compatibility with previous version, copy derivative into 3rd and
! 4th columns of psii
      call dcopy(2*nch,psir(ipoint),1,psii(ipoint),1)
    else
      do 300 ii=1, nch
        sc1(ii)=onemin/sc1(ii)
300       continue
! scmat2(nch, 4) = scmat(nch, nch) * psir(nch, 4)
      call mxma(scmat,1,nch,psir,1,nch,scmat2,1,nch,nch,nch,4)
      call dcopy(2*nch,scmat2,1,psii(ipoint),1)
      call dcopy(2*nch,scmat2(ipoint),1,psii,1)
! 1st two columns of psii now contain real and imaginary part of derivative
! 3rd and 4th columns contain real and imaginary parts of function
! propagate functions
      call vmul(sc2,1,psii(ipoint),1,scmat2,1,nch)
      call vmul(sc2,1,psii(ipoint+nch),1,scmat2(1+nch),1,nch)
! first and second columns of scmat2 now contain y4 Fb
      call daxpy_wrapper(2*nch,onemin,scmat2,1,scmat2(ipoint),1)
! 3rd and 4th columns of scmat2 now contains (Fb' - y4 Fb)
! if photodissociation, subtract off gamma2 from real part
      if (photof) call vadd(mone,scmat2(2*nch+1),1 &
                          ,psir(noffset),1,nch)
! if photodissociation, 3rd column of scmat2 now contains
!     Re(Fb' - y4 Fb - gamma 2)
! multiply by - y2^-1 to get Fa
      call vmul(sc1,1,scmat2(ipoint),1,psii,1,nch)
      call vmul(sc1,1,scmat2(ipoint+nch),1,psii(1+nch),1,nch)
! 1st and 2nd columns of psii now contain Fa
      do 320 ii=1, nch
        sc1(ii)=onemin/sc1(ii)
320       continue
      call vmul(sc1,1,psii(ipoint),1,scmat2,1,nch)
      call vmul(sc1,1,psii(ipoint+nch),1,scmat2(1+nch),1,nch)
! 1st and 2nd columns of scmat2 now contain y2 Fb
! if photodissociation subtract off gamma1 from real part
      if (photof) call vadd(mone,scmat2,1,sc9,1,nch)
      call dscal(nch,onemin,pk,1)
      call vmul(pk,1,psii,1,psii(ipoint),1,nch)
      call vmul(pk,1,psii(1+nch),1,psii(ipoint+nch),1,nch)
! 3rd and 4th columns of psii now contain -y1 Fa
! add on y2 Fb and store in 3rd and 4th columns of psii
      call daxpy_wrapper(2*nch,one,scmat2,1,psii(ipoint),1)
! if wavevector matrix positive (lambda negative) channel is closed
! kill the corresponding components of psi and psi'
! thresh is the threshold for killing closed channel components
    do 330 ii=1, nch
      if (sc8(ii) .lt. thresh) then
        fact=exp(-sqrt(abs(sc8(ii)))*drnow*factr)
        psii(ii)=fact*psii(ii)
        psii(ii+nch)=fact*psii(ii+nch)
        psii(ii+2*nch)=fact*psii(ii+2*nch)
        psii(ii+3*nch)=fact*psii(ii+3*nch)
      endif
330     continue
    endif
! 1st two columns of psii now contain Fa
! 3rd and 4th columns of psii now contain Fa'
    if (adiab) then
! here for flux calculation in locally adiabatic basis
      scsum=0.d0
      do 360 i=1, nj
        nni=nalist(i)
        sc(i)=psii(nni)*psii(3*nch+nni)- &
              psii(nch+nni)*psii(2*nch+nni)
!  store real or imaginary parts of wf if desired
        if (iwf .eq. 1) pk(i)=psii(nni)
        if (iwf .eq. -1) pk(i)=psii(nch+nni)
! if wavevector matrix positive (lambda negative) channel is closed
! kill the corresponding components of psi and psi'
!              if (sc8(nni) .lt. thresh.and.kill) then
        if (sc8(nni) .lt. 0.d0.and.kill) then
          sc(i)=zero
          if (iwf .ne. 0) pk(i)=0.d0
        endif
        scsum=scsum+sc(i)
360       continue
      nout=nj
    endif
! transform wavefunction and derivative into asymptotic basis
    call mxma(scmat,nch,1,psii,1,nch,psir,1,nch,nch,nch,4)
! psir now contains this information
!          write (23, 299) -r, drnow, (psir(ii), ii=1,8)
    if (.not.adiab) then
! here for flux calculation in asymptotic basis
! calculate flux
      if (coordf) then
! here for coordinate space calculation of fluxes
        scsum=zero
! fluxes will be stored in vector sc, initialize to zero
        call dset(nch,zero,sc,1)
        do i=1,ny
          ind=(i-1)*nch+1
          scc1=ddot(nch,psir,1,tcoord(ind),1)
! scc1 contains sum(psi-real*phi_internal)
          scc2=ddot(nch,psir(nch+1),1,tcoord(ind),1)
! scc2 contains sum(psi-imag*phi_internal)
          scc3=ddot(nch,psir(2*nch+1),1,tcoord(ind),1)
! scc3 contains sum(dpsi-real*phi_internal)
          scc4=ddot(nch,psir(3*nch+1),1,tcoord(ind),1)
! scc4 contains sum(dpsi-imag*phi_internal)
          sc(i)=scc1*scc4-scc2*scc3
!steve, you'll need to append to this to calculate r-component of current dens
! try using sc9 as scratch storage for these
        enddo
        call dscal(ny,dy,sc,1)
! multiply by step width to recover J.R at ri times dri
        do i=1, ny
          scsum=scsum+sc(i)
        enddo
        nout=ny
      endif
      if (.not. coordf) then
! here for channel fluxes in diabatic basis
! transform wavefunction and derivative into molecular basis (only for
! 13p or 2s-2p
        if (sumf) call dset(niout,zero,scc,1)
        if (ibasty .eq. 7) then
! psii(nch, 4) = ttrans(nch, nch) * psir(nch, 4)
          call mxma(ttrans,nch,1,psir,1,nch, &
                    psii,1,nch,nch,nch,4)
! N. B.  as written, mxma multiplies psir by the transpose of ttrans
! psii now contains this transformed wavefunction and derivative
! copy these back to psir
          call dcopy(4*nch,psii,1,psir,1)
        endif
        do 370 i=1, nj
          mmi=nalist(i)
          nni=nlist(i)
          sc(i)=psir(nni)*psir(3*nch+nni)- &
            psir(nch+nni)*psir(2*nch+nni)
!  store real or imaginary parts of wf if desired
          if (iwf .eq. 1) pk(i)=psii(nni)
          if (iwf .eq. -1) pk(i)=psii(nch+nni)
! if wavevector matrix positive (lambda negative) channel is closed
! kill the corresponding components of psi and psi'
          if (sc8(mmi) .lt. 0.d0.and.kill) then
!              if (sc8(mmi) .lt. thresh.and.kill) then
            sc(i)=zero
            if (iwf .ne. 0) pk(i)=0.d0
          endif
370         continue
        if (sumf) then
          do 390 i=1, nj
            mmi=nalist(i)
            do 375 ni=1, niout
              if (inq(mmi) .eq. indout(ni)) then
                scc(ni)=scc(ni)+sc(mmi)
              endif
375             continue
390           continue
        endif
        if (.not.sumf) then
          nout=nj
          scsum=dsum(nj,sc,1)
        else
          nout=niout
          scsum=dsum(niout,scc,1)
        endif
! transform wavefunction and derivative back into asymptotic basis (only for
! 13p or 2s-2p
        if (ibasty .eq. 7) then
          call mxma(ttrans,1,nch,psir,1,nch, &
                    psii,1,nch,nch,nch,4)
! psii now contains this transformed wavefunction and derivative
! copy these back to psir
          call dcopy(4*nch,psii,1,psir,1)
        endif
      endif
    endif
    if (iwf .ne. 0) &
      write (2, 400) -r, (pk(i), i=1,nj)
    if (iwf .eq. 0) then
      if (photof) then
! for photodissociation, so, steve, you'll need to scale sc9 also
        call dscal(nout,1.d0/rmu,sc,1)
        scsum=scsum/rmu
        if (scsum .lt. 1.d-13) scsum=zero
      endif
      if (sumf) then
        write (3, 400) -r, (scc(i), i=1,nout), scsum
      else
        write (3, 400) -r, (sc(i), i=1,nout), scsum
! steve you'll also have to write out sc9
      endif
400       format(f10.4,30(1pe11.3))
    endif
420   continue

450   continue
!        close (23)
  return
!
900   continue
950   write (0, *) '*** ERROR READING WFU FILE (FLUX). ABORT'
  call exit()
  end
! ------------------------------------------------------------------
subroutine eadiab(npts,nch,nchsq,nj)
!
! subroutine to readin and print out adiabatic energies
!
! author: millard alexander
! current revision date (algorithm): 12-may-1997 by mha
! revised on 30-mar-2012 by q. ma for stream I/O of wfu files
!
! ------------------------------------------------------------------
use mod_cocent, only: sc2 => cent
use mod_coisc3, only: nalist => isc3 ! nalist(10)
use mod_cosc6, only: sc => sc6 ! sc(6)
use mod_cosc7, only: sc7  ! sc7(6)
use mod_cosc8, only: sc8
use mod_coz, only: scmat => z_as_vec ! scmat(100)
implicit double precision (a-h,o-z)
common /coered/ ered, rmu
common /cowave/ irec, ifil, nchwfu, ipos2, ipos3, nrlogd, iendwv, &
     inflev
! common for y1, y2, y4
data zero, one, onemin, two,   conv &
    /0.d0, 1.d0, -1.d0, 2.d0, 219474.6d0/
  do 420 kstep=1, npts
    irec=irec-1
    read (ifil, end=900, err=950, pos=iwavsk(irec)) r
! branch out if we've gone beyond airy propagation region
    if (r .gt. 0) return
    read (ifil, end=900, err=950) drnow
! read in local wavevectors
    read (ifil, end=900, err=950) (sc8(i), i=1, nch)
    do 370 i=1, nj
      nni=nalist(i)
      sc(i)=-conv*(sc8(nni)/(two*rmu)-ered)
370     continue
!          write (3, 170) -r+0.5*drnow, (nalist(i), i=1,nj),
    write (3, 170) -r+0.5*drnow, &
        (sc(i), i=1,nj)
!170       format(f10.4,2i4,20(1pe12.4))
170     format(f10.4,20(1pe12.4))
420   continue
  return
!
900   continue
950   write (0, *) '*** ERROR READING WFU FILE (EADIAB). ABORT'
  call exit()
  end
! ------------------------------------------------------------------
subroutine waverot(jtot,nch)
! special for singlet-triplet mixing;
! rearrange wavefunction to correspond to state assignment:
! psi-parallel = (Jtot+1)^1/2 |el=J+1> - Jtot^1/2 |el=J-1>
! psi-perpendicular = Jtot^1/2 |el=J+1> + (Jtot+1)^1/2 |el=J-1>
! with assignment that original column five is singlet with el=J-1 and
! column 6 is singlet with el=J+1
! use orthogonal plane rotation
use mod_coamat, only: psir ! psir(100) (4,nch)
use mod_cobmat, only: psii ! psii(100) Here psii is used as a vector
use mod_cotq1, only: dpsir ! dpsir(100) 
use mod_cotq2, only: dpsii ! dpsii(100)
implicit none
integer, intent(in) :: jtot
integer, intent(in) :: nch
integer :: jpoint
real(8) :: cs, sn ! cosine and sine of rotation angle
real(8) :: jtot_as_real
real(8) :: norm
jtot_as_real=jtot
cs=sqrt(jtot_as_real+1.d0)
sn=sqrt(jtot_as_real)
norm=sqrt(2.d0*jtot_as_real+1.d0)
cs=cs/norm
sn=sn/norm
jpoint=4*nch+1
call drot(nch,psir(jpoint),1,psir(jpoint+nch),1,cs,sn)
call drot(nch,psii(jpoint),1,psii(jpoint+nch),1,cs,sn)
call drot(nch,dpsir(jpoint),1,dpsir(jpoint+nch),1,cs,sn)
call drot(nch,dpsii(jpoint),1,dpsii(jpoint+nch),1,cs,sn)
return
end
! ------------------------------------------------------------------
subroutine psicalc(npts,nch,nchsq,nj)
!
! subroutine to propagate wavefunctions inward
!
! author: millard alexander
! current revision date (algorithm): 4-oct-1991 by mha
! revised on 30-mar-2012 by q. ma for stream I/O of wfu files
!
! ------------------------------------------------------------------
use mod_coamat, only: psir ! psir(100)  psir(nch, 1)
use mod_coisc2, only: nlist => isc2 ! nlist(6)
use mod_cosc6, only: sc => sc6 ! sc(6)
use mod_cosc7, only: sc1 => sc7 ! sc1(6)
use mod_cow, only: sr => w_as_vec ! sr(100)
use mod_cozmat, only: si => zmat_as_vec ! si(100)

implicit double precision (a-h,o-z)
common /cowave/ irec, ifil, nchwfu, ipos2, ipos3, nrlogd, iendwv, &
     inflev
  irec=npts+4
  do 180 kstep=1, npts
    irec=irec-1
    read (ifil, end=900, err=950, pos=iwavsk(irec)) r, drnow
! read in adiabtic energies (which we don't need)
    read (ifil, end=900, err=950) (sc1(i), i=1, nch)
! read in transformation matrices
    read (ifil, end=900, err=950) (sr(i), i=1, nchsq)
! si(nch, 1) = sr(nch, nch) * psir(nch, 1)
    call mxma(sr,1,nch,psir,1,nch,si,1,nch,nch,nch,1)
    do 160 i=1, nj
      sc(i)=si(nlist(i))
160     continue
    call dcopy(nch,si,1,psir,1)
    write (2, 170) r, (sc(i), i=1,nj)
170     format(f10.4,12(1pe10.2))
180   continue
return
!
900 continue
950 write (0, *) '*** ERROR READING WFU FILE (PSICALC). ABORT.'
call exit()
end
! ------------------------------------------------------------------
subroutine transmt(npts,nch,nchsq,rout)
!
! subroutine to print out transformation matrix at rout
!
! author: millard alexander
! current revision date (algorithm): 18-may-2008 by mha
! revised on 30-mar-2012 by q. ma for stream I/O of wfu files
! current revision: 20-apr-2012 by q. ma
!
! ------------------------------------------------------------------
use mod_coisc2, only: nlist => isc2 ! nlist(6)
use mod_cosc6, only: sc => sc6 ! sc(6)
use mod_cosc7, only: sc1 => sc7 ! sc1(6)
use mod_cow, only: sr => w_as_vec ! sr(100)
use mod_cozmat, only: si => zmat_as_vec ! si(100)
implicit double precision (a-h,o-z)
logical renormf
common /cowave/ irec, ifil, nchwfu, ipos2, ipos3, nrlogd, iendwv, &
     inflev
dimension scrvec(64)
irec=npts+4
delold=1.d+18
do 200 kstep=1, npts
    irec=irec-1
    read (ifil, end=900, err=950, pos=iwavsk(irec)) r, drnow
    del=abs(abs(r)-rout)
! read in adiabtic energies (which we don't need)
    read (ifil, end=900, err=950) (sc1(i), i=1, nch)
! read in first G(n,n+1) matrix,  then Tn transformation
! matrix into local interval
    read (ifil, end=900, err=950) (sr(i), i=1, nchsq), &
         (si(i), i=1, nchsq)
    if (rout .gt. 0.d0) then
      if (del .gt.delold) then
! transpose matrix is stored column by column
! after transposition, columns correspond to eigenvectors
        call transp (si, nch, nch)
        write (3,160) -r+0.5*drnow
160    format('    TRANSFORMATION MATRIX AT R = ',f10.6,' IS:',/)
        write (6,161) abs(r)
161    format('    TRANSFORMATION MATRIX DETERMINED AT R = ',f10.6)
! print out transpose matrix (reverse order since eispack
! determines highest eigenvalue first
        ind=nchsq-nch+1
! if more than 64 channels, then eliminate all elements with i>64, and renormalize vector
        nnch=nch
        if (nch.gt.64) then
           renormf=.true.
           write (3,162) nch
           write (6,162) nch
162    format('    CHANNEL EXPANSION TRUNCATED FROM NCH =',i4, &
             ' TO NCH = 64;')
           write (3,163)
           write (6,163)
163    format('       EIGENVECTOR RENORMALIZED')
           nnch=64
        else
           renormf=.false.
        endif
        call openf(4,'tmatrix.dat','sf',0)
        nstate=min(12,nch)
        do 190 ii=1,nstate
!             do 190 ii=1,nch
          if (renormf) then
             call dcopy(64,si(ind),1,scrvec,1)
             cnorm=ddot(64,scrvec,1,scrvec,1)
             cnorm=1d0/sqrt(cnorm)
             call dscal(64,cnorm,scrvec,1)
             call dcopy(64,scrvec,1,si(ind),1)
          endif
          write (3,175) -r+0.5*drnow, ii, &
            (si(ij), ij=ind,ind+nnch-1)
          write (4,175) -r+0.5*drnow, ii, &
            (si(ij), ij=ind,ind+nnch-1)
175         format(f7.4,i3,64f10.6)
          ind=ind-nch
190         continue
        close(4)
        return
      endif
    else
        call transp (si, nch, nch)
! print out transpose matrix (reverse order since eispack
! determines highest eigenvalue first
        ind=nchsq-nch+1

        call openf(4,'tmatrix.dat','sf',0)
        do 195 ii=1,nstate
          if (renormf) then
             call dcopy(64,si(ind),1,scrvec,1)
             cnorm=ddot(64,scrvec,1,scrvec,1)
             cnorm=1d0/sqrt(cnorm)
             call dscal(64,cnorm,scrvec,1)
             call dcopy(64,scrvec,1,si(ind),1)
          endif
!              do 195 ii=1,nch
          write (4,175) -r+0.5*drnow, ii, &
            (si(ij), ij=ind,ind+nnch-1)
          ind=ind-nch
195         continue
        close(4)
    endif
    call dcopy(nchsq,sr,1,si,1)
    delold=del
200 continue
! here if rout not reached at end of airy propagation region, print out
! last transformation matrix
call transp (si, nch, nch)
write (3,160)-r+0.5*drnow
write (6,160)-r+0.5*drnow
! print out transpose matrix
ind=nchsq-nch+1
do 250 ii=1,12
!      do 250 ii=1,nch
  write (3,175) -r+0.5*drnow,  ii, (si(ij), ij=ind,ind+nch-1)
  ind=ind-nch
250 continue
return
!
900 continue
950 write (0, *) '*** ERROR READING WFU FILE (TRANSMT). ABORT'
call exit()
end
!
!     ------------------------------------------------------------------
subroutine eadiab1(filnam, nchmin, nchmax)
!
!     Subroutine to readin and print out adiabatic energies from a wfu
!     file.  wfu files with only adiabatic energy information
!     (wrsmat=.f. when performaing wavefunction calculations) can be
!     supported.
!
!     Adiabatic energies as a function of r will be saved to
!     filnam.eadiab, with columns the adiabatic bender curves.  r is
!     stored in the first column.
!
!     author: qianli ma
!     current revision: 20-apr-2012
!
!     ------------------------------------------------------------------
use constants
use mod_cosc8, only: sc8
implicit double precision (a-h, o-z)
character*(*) filnam
logical exstfl
character*40 wavfil, eadfil
!
common /coered/ ered, rmu
common /cowave/ irec, ifil, nchwfu, ipos2, ipos3, nrlogd, iendwv, &
     inflev
!
double precision dble_t
!
if (nchmax .ne. 0 .and. nchmax .lt. nchmin) goto 990
ien = 0
wavfil = filnam // '.wfu'
call gennam(wavfil, filnam, ien, 'wfu', lenfs)
inquire (file=wavfil, exist=exstfl)
if (.not. exstfl) then
   write (6, 10) wavfil(1:lenfs)
10    format (' ** WAVEFUNCTION INFORMATION FILE ', (a), &
        ' NOT FOUND **')
   return
end if
call openf(22, wavfil, 'TU', 0)

eadfil = filnam // '.eadiab'
call gennam(eadfil, filnam, ien, 'eadiab', lenft)
call openf(2, eadfil, 'sf', 0)
write (6, 15) eadfil(1:lenft)
15 format (' ** WRITING ADIABATIC ENERGIES TO ', (a))
!
call waverd(jtot, jlpar, nu, nch, npts, nopen, nphoto, &
     0, rstart, rendld, rinf)
if (nchmin .gt. nch) goto 990
if (nchmax .eq. 0 .or. nchmax .gt. nch) nchmax = nch
nchpr = nchmax - nchmin + 1
noffst = (nch - nchmax + 2) * sizeof(dble_t)
!
write (2, 17) nchmin, nchmax
17 format (' ** ADIABATIC ENERGIES FROM NO.', i5, ' TO NO.', &
     i5, ' REQUESTED')
do i = 4 + nrlogd, npts + 3
   read (ifil, end=900, err=950, pos=iwavsk(i)) r, drnow
   read (ifil, end=900, err=950, pos=iwavsk(i)+noffst) &
        (sc8(j), j=1, nchpr)
   write (2, 20) -r + 0.5 * drnow
20    format (f10.5, 1x, $)
   do j = nchpr, 1, -1
      write (2, 30) -econv * (sc8(j) / (2d0 * rmu) - ered)
   end do
30    format (f13.5, 1x, $)
   write (2, 20)
end do
write (2, 20)
goto 990
!
900 continue
950 write (0, *) '*** ERROR READING WFU FILE'
990 close(ifil)
close(3)
return
end
