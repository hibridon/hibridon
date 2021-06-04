      subroutine airymp (x, ftheta, fphi, xmmod, xnmod, scai, scbi,
     :                    zeta)
*  subroutine to return the moduli and phases of the airy functions and their
*  derivatives
*  author:  millard alexander
*  current revision date: 4-nov-1991
* ---------------------------------------------------------------------------
*  variables in call list:
*    x     argument of airy functions
*    ftheta, xmmod         on return: contain the (double precision)
*                                 phase and modulus of ai(x) and bi(x) (see
*                                 below).
*    fphi, xnmod           on return: contain the (double precision)
*                                 phase and modulus of ai'(x) and Bi'(x) (see
*                                 below).
*    scai, scbi, zeta      on return:  for x .le. 0, returns constants such th
c     ai(x) = scai * cos(zeta) + scbi * sin(zeta)
c     bi(x) = scbi * cos(zeta) - scai * sin(zeta)
c     where zeta = (2/3) * (-x) ** (3/2) + pi/4
*                          and, for x .gt. 0
*     ai(x) = scai * exp(-zeta)
*     bi(x) = scbi * exp(+zeta)
*     where zeta = (2/3) * x ** (3/2)
*    N.B. on entry scai, scbi are identical to above values, except for
*    -5 .le. x .le. 0, where on entry scai = ai(x) and scbi = bi(x)

* ----------------------------------------------------------------------------
*  for negative x
* ----------------------------------------------------------------------------
*  the moduli and phases are defined by
*      ai(-x) = m(x) cos[theta(x)]
*      bi(-x) = m(x) sin[theta(x)]
*      ai'(-x) = n(x) cos[phi(x)]
*      bi'(-x) = n(x) sin[phi(x)]
*  in other words
*          2              2        2
*      m(x)  = sqrt[ ai(x)  + bi(x)  ]
*          2               2         2
*      n(x)  = sqrt[ ai'(x)  + Bi'(x)  ]
*      theta(x) = atan [ bi(x) / ai(x) ]
*      phi(x)   = atan [ bi'(x) / Ai'(x) ]
*  to determine these moduli and phases we use the subroutine
*  scairy, written by d. manolopoulos (sept. 1986)
*  this subroutine returns the following quantities:
*     scai, scbi, scaip, scpib, and zeta, where
c     for  x .lt. -5.0
c     ai(x) = scai * cos(zeta) + scbi * sin(zeta)
c     bi(x) = scbi * cos(zeta) - scai * sin(zeta)
c     ai'(x) = scaip * cos(zeta) + scbip * sin(zeta)
c     bi'(x) = scbip * cos(zeta) - scaip * sin(zeta)
c     where zeta = (2/3) * (-x) ** (3/2) + pi/4
c
c     for  -5.0 .le. x .le. 0.0
c
c     ai(x) = scai
c     bi(x) = scbi
c     ai'(x) = scaip
c     bi'(x) = scbip
c     and zeta = 0
* ----------------------------------------------------------------------------
*  for positive x ( x > 0)
* ----------------------------------------------------------------------------
*  the moduli and phases are defined by
*      ai(x) = m(x) sinh[theta(x)]
*      bi(x) = m(x) cosh[theta(x)]
*      ai'(x) = n(x) sinh[phi(x)]
*      bi'(x) = n(x) cosh[phi(x)]
*  in other words
*          2              2        2
*      m(x)  = sqrt[ bi(x)  - ai(x)  ]
*          2               2         2
*      n(x)  = sqrt[ bi'(x)  - ai'(x)  ]
*      theta(x) = atanh [ ai(x) / bi(x) ]
*      phi(x)   = atanh [ ai'(x) / bi'(x) ]
*  here the the exponentially scaled airy functions
*  ai(x), ai'(x), bi(x), bi'(x) are:
*      ai(x)  = ai(x)  * exp[zeta]
*      ai'(x) = ai'(x) * exp[zeta]
*      bi(x)  = bi(x)  * exp[-zeta]
*      bi'(x) = bi'(x) * exp[-zeta]
*  to determine these moduli and phases we use the subroutine
*  scairy, written by d. manolopoulos (sept. 1986)
*  this subroutine returns the following quantities:
*     scai, scbi, scaip, scpib, and zeta
*  in terms of which the exponentially scaled airy functions are defined by
*   ai(x) = scai * exp(-zeta)
*   bi(x) = scbi * exp(+zeta)
*   ai'(x) = scaip * exp(-zeta)
*   bi'(x) = scbip * exp(+zeta)
*   where zeta = (2/3) * x ** (3/2)
*
* ----------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      data twoth/0.666666666666666667d0/
      data pib4 / 7.85398163397448310d-01 /
      call scairy (x, scai, scbi, scaip, scbip, zeta)
      if ( x .le. 0.d0) then
        xmmod = sqrt( scai ** 2 + scbi ** 2)
        xnmod = sqrt( scaip ** 2 + scbip ** 2)
        ftheta = atan2 (scbi, scai)
        fphi = atan2 (scbip, scaip)
        if (x .lt. (-5.0d0) ) then
          ftheta = ftheta - zeta
          fphi = fphi - zeta
        else
          zeta=twoth*(-x)*sqrt(-x)+pib4
          sn=sin(zeta)
          cs=cos(zeta)
          ai=scai*cs-scbi*sn
          bi=scbi*cs+scai*sn
          scai=ai
          scbi=bi
        end if
      else
        xmmod = sqrt( - scai ** 2 + scbi ** 2)
        xnmod = sqrt( - scaip ** 2 + scbip ** 2)
        ratio = scai / scbi
        ftheta = 0.5 * log ( (1.d0 + ratio) / (1.d0 - ratio) )
        ratio = scaip / scbip
        fphi = 0.5 * log ( (1.d0 + ratio) / (1.d0 - ratio) )
      end if
      return
      end
