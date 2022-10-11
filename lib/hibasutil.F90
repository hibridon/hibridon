module mod_hibasutil
contains
!     ------------------------------------------------------------------
subroutine raise(mesg)
implicit none
character*(*) mesg
write (0, *) ' *** hibridon error: ', mesg
stop
end

! ----------------------------------------------------------------------
subroutine iswap(ia,ib)
implicit integer (a-z)
ii=ia
ia=ib
ib=ii
return
end
! ----------------------------------------------------------------------
subroutine rswap(a,b)
implicit double precision (a-h,o-z)
c=a
a=b
b=c
return
end
	
! ----------------------------------------------------------------------
subroutine prmatp (jp, lp, j, l, jtot, kp, k, lambda, mu, &
	vprm, csflag)
!  subroutine to calculate primitive v-lambda matrix elements for close-coupled
!  treatment of collisions of a symmetric top with an atom
!  the primitive cc and cs matrix elements are given in eqs. (26) and (32),
!  respectively, of s. green, j. chem. phys. 64, 3463 (1976)
!  note, that in this article the bra indices are unprimed and the ket indices
!  primed, while in the conventions of the present subroutine the bra indices
!  are primed and the ket indices, unprimed.
!
!  duplicate of prmstp
!  author:  millard alexander
!  current revision date:  27-mar-90
!  -----------------------------------------------------------------------
!  variables in call list:
!    jp:       rotational quantum number of left side of matrix element (bra)
!    lp:       orbital angular momentum of left side of matrix element (bra)
!    j:        rotational quantum number of right side of matrix element (ket)
!    l:        orbital angular momentum of right side of matrix element (ket)
!    jtot:     total angular momentum
!    kp:       k quantum number of bra
!    k:        k quantum number of ket
!    lambda:   order of legendre term in expansion of potential
!    mu:       index of legendre term in expansion of potential
!    vrpm:     on return, contains primitive matrix element
!    csflag:   if .true., then cs calculation; in this case:
!                jtot in call list is cs lbar
!                lp is nu (cs projection index)
!                l is not used
!              if .false., then cc calculation; in this case:
!                jtot is total angular momentum
!  subroutines called:
!     xf3j, xf6j
!  -----------------------------------------------------------------------
use mod_hiutil, only: xf3j, xf6j
implicit double precision (a-h,o-z)
logical csflag
data half, one, two, zero, four / 0.5d0, 1.d0, 2.d0, 0.d0, 4.d0/
data pi /3.1415926536d0 /
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
xnorm = (two * xjp + one) * (two * xj + one) * &
(two * lambda + one) / (four * pi)
!  special normalization for k and/or kp = 0
!  see Eq. (46) of S. Green, j. chem. phys. 64, 3463 (1976)
if (k .eq. 0) xnorm = xnorm * half
if (kp .eq. 0) xnorm = xnorm * half
!  the desired matrix element is constructed in x
if (.not. csflag) then
!  here for cc matrix elements
x = xf3j (xlp, xlamda, xl, zero, zero, zero)
if (x .eq. zero) return
x = x * xf3j (xjp, xj, xlamda, xkp, - xk, xmu)
if (x .eq. zero) return
x = x * xf6j (xj, xl, xjtot, xlp, xjp, xlamda)
if (x .eq. zero) return
iphase = jp + j + kp - jtot
xnorm = xnorm * (two * lp + one) * (two * l + one)
else if (csflag) then
!  here for cs matrix elements
iphase = - k - nu
x = xf3j (xjp, xlamda, xj, -xnu, zero, xnu)
if (x .eq. zero) return
x = x * xf3j (xjp, xlamda, xj, -xkp, xmu, xk)
end if
vprm = ( (-1) ** iphase) * x * sqrt(xnorm)
return
end
! ----------------------------------------------------------------------
subroutine prmstp (jp, lp, j, l, jtot, kp, k, lambda, mu, &
	vprm, csflag)
!  subroutine to calculate primitive v-lambda matrix elements for close-coupled
!  treatment of collisions of a symmetric top with an atom
!  the primitive cc and cs matrix elements are given in eqs. (26) and (32),
!  respectively, of s. green, j. chem. phys. 64, 3463 (1976)
!  note, that in this article the bra indices are unprimed and the ket indices
!  primed, while in the conventions of the present subroutine the bra indices
!  are primed and the ket indices, unprimed.
!  author:  millard alexander
!  current revision date:  27-mar-90
!  -----------------------------------------------------------------------
!  variables in call list:
!    jp:       rotational quantum number of left side of matrix element (bra)
!    lp:       orbital angular momentum of left side of matrix element (bra)
!    j:        rotational quantum number of right side of matrix element (ket)
!    l:        orbital angular momentum of right side of matrix element (ket)
!    jtot:     total angular momentum
!    kp:       k quantum number of bra
!    k:        k quantum number of ket
!    lambda:   order of legendre term in expansion of potential
!    mu:       index of legendre term in expansion of potential
!    vrpm:     on return, contains primitive matrix element
!    csflag:   if .true., then cs calculation; in this case:
!                jtot in call list is cs lbar
!                lp is nu (cs projection index)
!                l is not used
!              if .false., then cc calculation; in this case:
!                jtot is total angular momentum
!  subroutines called:
!     xf3j, xf6j
!  -----------------------------------------------------------------------
use mod_hiutil, only: xf3j, xf6j
implicit double precision (a-h,o-z)
logical csflag
data half, one, two, zero, four / 0.5d0, 1.d0, 2.d0, 0.d0, 4.d0/
data pi /3.1415926536d0 /
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
xnorm = (two * xjp + one) * (two * xj + one) * &
(two * lambda + one) / (four * pi)
!  special normalization for k and/or kp = 0
!  see Eq. (46) of S. Green, j. chem. phys. 64, 3463 (1976)
if (k .eq. 0) xnorm = xnorm * half
if (kp .eq. 0) xnorm = xnorm * half
!  the desired matrix element is constructed in x
if (.not. csflag) then
!  here for cc matrix elements
x = xf3j (xlp, xlamda, xl, zero, zero, zero)
if (x .eq. zero) return
x = x * xf3j (xjp, xj, xlamda, xkp, - xk, xmu)
if (x .eq. zero) return
x = x * xf6j (xj, xl, xjtot, xlp, xjp, xlamda)
if (x .eq. zero) return
iphase = jp + j + kp - jtot
xnorm = xnorm * (two * lp + one) * (two * l + one)
else if (csflag) then
!  here for cs matrix elements
iphase = - k - nu
x = xf3j (xjp, xlamda, xj, -xnu, zero, xnu)
if (x .eq. zero) return
x = x * xf3j (xjp, xlamda, xj, -xkp, xmu, xk)
end if
vprm = ( (-1) ** iphase) * x * sqrt(xnorm)
return
end
! ----------------------------------------------------------------------
subroutine prmstpln (jp, lp, j2p, j12p, j, l, j2, j12, &
	jtot, kp, k, lambda1, mu1,lambda2, mu2, &
		vprm, csflag)
!  subroutine to calculate primitive v-lambda matrix elements for
!  close-couple
!  treatment of collisions of a symmetric top with a linear molecule
!  author: claire rist
!  current revision date:  17-jan-1992
!
! -----------------------------------------------------------------------
!  variables in call list:
!    jp:       rotational quantum number of left side of matrix element
!              (bra) for symmetric top
!    lp:       orbital angular momentum of left side of matrix element
!              (bra)
!    j2p:      rotational quantum number of left side of matrix element
!              (bra) for linear molecule
!    j12p:     total molecular rotational quantum number of left side of
!              matrix element (bra)
!    j:        rotational quantum number of right side of matrix element
!              (ket) for symmetric top
!    l:        orbital angular momentum of right side of matrix element
!              (ket)
!    j2:       rotational quantum number of right side of matrix element
!              (ket) for linear molecule
!    j12:      total molecular rotational quantum number of right side
!              of matrix element
!    jtot:     total angular momentum
!    kp:       k quantum number of bra
!    k:        k quantum number of ket
!    lambda1:  order of rotation matrix term in expansion of potential
!    mu1:      index of rotation matrix term in expansion of potential
!    lambda2:  order of rotation matrix term in expansion of potential
!    mu2:      index of rotation matrix term in expansion of potential
!    vrpmstp:  on return, contains primitive matrix element
!    csflag:   if .true., then cs calculation; in this case:
!                jtot in call list is cs lbar
!                lp is nu (cs projection index)
!                l is not used
!              if .false., then cc calculation; in this case:
!                jtot is total angular momentum
!  subroutines called:
!     xf3j, xf6j
!     f9j (Claire's routines)
!
! -----------------------------------------------------------------------
use mod_hiutil, only: xf3j, xf6j, f9j
implicit double precision (a-h,o-z)
logical csflag
data half, one, two, zero, four / 0.5d0, 1.d0, 2.d0, 0.d0, 4.d0/
data pi /3.1415926536d0 /
vprm = zero
xjp = dble(jp)
xj = dble(j)
xkp = dble(kp)
xk = dble(k)
xj2p = dble(j2p)
xj2 = dble(j2)
xj12p = dble(j12p)
xj12 = dble(j12)
xjtot = dble(jtot)
xlp = dble(lp)
xl = dble(l)
xlambda1 = dble(lambda1)
xmu1 = dble(mu1)
xlambda2 = dble(lambda2)
xmu2 = dble(mu2)
! Last term in j12 added by Claire to reproduce the isotropic
! matrix elements
xnorm = (two * xjp + one) * (two * xj + one) * &
(two * xj2p + one) * (two * xj2 + one) * &
(two * xlp + one) * (two * xl + one)* &
(two * j12p + one) * (two * j12 + one)
!  special normalization for k and/or kp = 0
if (k .eq. 0) xnorm = xnorm * half
if (kp .eq. 0) xnorm = xnorm * half
!  special normalization for mu1 and/or mu2 = 0
!      if (mu1 .eq. 0) xnorm = xnorm * half * half
if (mu2 .eq. 0) xnorm = xnorm * half * half
!  the desired matrix element is constructed in x
!  here for cc matrix elements
x = xf3j (xjp, xlambda1, xj, xkp, xmu1, - xk)
x2 = xf3j (xj2p, xlambda2, xj2, zero, zero, zero)
x = x * xf3j (xj2p, xlambda2, xj2, zero, zero, zero)
if (x .eq. zero) return
lmin = max(iabs(lambda1-lambda2), iabs(l-lp))
lmax = min(lambda1+lambda2, l+lp)
sum = zero
do 500 lambda = lmin, lmax
xlambda = dble(lambda)
xs = xf3j (xlp, xlambda, xl, zero, zero, zero)
if (xs .eq. zero) go to 500
xs = xs * xf3j (xlambda1,xlambda2,xlambda,xmu2,-xmu2,zero)
if (xs .eq. zero) go to 500
xs = xs * xf6j (xj12, xl, xjtot, xlp, xj12p, xlambda)
if (xs .eq. zero) go to 500
!        xs4 = 1/sqrt((2*xjp+1)*(2*xj+1)*(2*xlambda+1))
!        xs = xs * xs4
xs = xs * f9j(lambda1, j, jp, lambda2, j2, j2p, &
lambda, j12, j12p)
if (xs .eq. zero) go to 500
xs = xs * (two * xlambda + one)
!        write(6,*)'lambda,  xs', lambda, xs
sum = sum + xs
500 continue
x = x * sum
!      iphase = jtot+k+mu2
iphase = jtot + k + mu2 + jp + j2p + j12p
vprm = ( (-1) ** iphase) * x * sqrt(xnorm)
! Claire factor ( 1+ eps*epsp*(-1)**(j+jp+lambda1+mu1))=2
vprm = two * vprm
!     provisoire: check with nh3h2(j=0) calculations
!     multiplicationfactor:
!      s4pi = sqrt ( 4.d0 * acos(-1.d0) )
!      vprm = vprm *sqrt(2*xlambda1+1) /s4pi
!      cnorm= two*sqrt(2*xlambda1+1)*sqrt(xnorm)/s4pi
return
end

! ----------------------------------------------------------------------
double precision function rotham(ji, ki, jf, kf)
!
!  subroutine to compute matrix elements of the asymmmetric top hamiltionian
!  in a prolate (case Ia) basis between unsymmetrized basis functions
!  (ji,ki) and (jf,kf)
!
!  author:  paul dagdigian
!  current revision date:  16-aug-2009
!  -----------------------------------------------------------------------
use mod_cosysr, only: isrcod, junkr, rspar
implicit double precision (a-h,o-z)

real(8), pointer :: arot, brot, crot, emax
arot=>rspar(1); brot=>rspar(2); crot=>rspar(3); emax=>rspar(4)

bpc = (brot + crot)*0.5d0  ! B Plus C
bmc = (brot - crot)*0.25d0 ! B Minus C
if (ji .ne. jf) goto 900
fjj1 = ji*(ji + 1.d0)
fk = ki
if (ki .ne. kf) goto 10
!
!     diagonal term
rotham = bpc*(fjj1 - fk**2) + arot*fk**2
goto 1000
!
!     off-diagonal terms
10 if (kf .ne. (ki + 2)) goto 20
rotham = bmc*sqrt((fjj1 - fk*(fk + 1.d0)) &
  *(fjj1 - (fk + 1.d0)*(fk + 2.d0)))
goto 1000
20 if (kf .ne. (ki - 2)) goto 900
rotham = bmc*sqrt((fjj1 - fk*(fk - 1.d0)) &
  *(fjj1 - (fk - 1.d0)*(fk - 2.d0)))
goto 1000
900 rotham = 0.d0
1000 continue
return
end

! ----------------------------------------------------------------------
double precision function rotham1(ji, ki, jf, kf)
!
!  subroutine to compute matrix elements of the asymmmetric top hamiltionian
!  in a prolate (case Ia) basis between unsymmetrized basis functions
!  (ji,ki) and (jf,kf)
!
!  here, the body-frame quantization axis is along the C2 axis
!  (b inertial axis) of the symmetrical molecule
!
!  modification of rotham subroutine
!
!  author:  paul dagdigian
!  current revision date:  18-sep-2017
!  -----------------------------------------------------------------------
use mod_cosysr, only: isrcod, junkr, rspar
implicit double precision (a-h,o-z)

real(8), pointer :: arot, brot, crot, emax
arot=>rspar(1); brot=>rspar(2); crot=>rspar(3); emax=>rspar(4)

aa = brot
bb = arot
cc = crot
bpc = (bb + cc)*0.5d0
bmc = (bb - cc)*0.25d0
if (ji .ne. jf) goto 900
fjj1 = ji*(ji + 1.d0)
fk = ki
if (ki .ne. kf) goto 10
!
!     diagonal term
rotham1 = bpc*(fjj1 - fk**2) + aa*fk**2
goto 1000
!
!     off-diagonal terms
10 if (kf .ne. (ki + 2)) goto 20
rotham1 = bmc*sqrt((fjj1 - fk*(fk + 1.d0)) &
  *(fjj1 - (fk + 1.d0)*(fk + 2.d0)))
goto 1000
20 if (kf .ne. (ki - 2)) goto 900
rotham1 = bmc*sqrt((fjj1 - fk*(fk - 1.d0)) &
  *(fjj1 - (fk - 1.d0)*(fk - 2.d0)))
goto 1000
900 rotham1 = 0.d0
1000 continue
return
end

!  -----------------------------------------------------------------------
subroutine vlm2sg (jp, lp, j, l, jtot, lambda, iepsp, ieps, &
	v, csflag)
!  subroutine to calculate v-lambda matrices for close-coupled and
!  coupled-states treatments of collisions of a molecule in a 2sigma
!  electronic state
!  latest revision date:  21-feb-89
!  the cc matrix elements are given in eq. (15) of m.h. alexander,
!  j. chem. phys. 76, 3637 (1982)
!  the cs matrix elements are given in eq. (24) of m.h. alexander,
!  j. chem. phys. 76, 3637 (1982), with corrections corresponding
!  to eq. 14 of t. orlikowski and m.h. alexander, j. chem. phys. 79, 6006
!  (1983)
!  note that for cc collisions of a 2sigma molecule with a flat surface, the
!  coupling matrix elements are assumed to be identical to the cs matrix
!  elements here
!  -----------------------------------------------------------------------
!  variables in call list:
!    jp:       rotational quantum number of left side of matrix element (bra)
!              minus 1/2
!    lp:       orbital angular momentum of left side of matrix element (bra)
!    j:        rotational quantum number of right side of matrix element (ket)
!              minus 1/2
!    l:        orbital angular momentum of right side of matrix element (ket)
!    jtot:     total angular momentum
!    lambda:   order of legendre term in expansion of potential
!    iepsp:    symmetry index of bra
!    ieps:     symmetry index of ket
!    v:        on return, contains matrix element
!    csflag:   if .true., then cs calculation; in this case:
!                jtot in call list is cs lbar
!                lp is nu (cs projection index) minus 1/2
!                l is not used
!              if .false., then cc calculation; in this case:
!                jtot is total angular momentum minus 1/2
!    for collisions of a 2sigma molecule with a surface, nu is equivalent
!    to m (the projection of j along the surface normal) minus 1/2
!  subroutines called:
!     xf3j, xf6j
!  -----------------------------------------------------------------------
use mod_hiutil, only: xf3j, xf6j
implicit double precision (a-h,o-z)
integer jp, j, jtot, lp, l, lambda, iepsp, ieps, nu, iphase
logical csflag
half = 0.5d0
zero = 0.d0
one = 1.d0
two = 2.d0
v = zero
xjp = jp + half
xj = j + half
xjtot = jtot + half
if (csflag) then
nu = lp
xnu = nu + half
end if
xlp = float(lp)
xl = float(l)
xlamda = float(lambda)
iphase = ieps * iepsp * ((-1) ** (jp + j + lambda + 1))
if (iphase .eq. 1) return
if (csflag) then
!  note true phase factor is nu - omega, where nu is a half/integer and
!  omega=1/2.  here, however nu is nu-true + 1/2, so nu-true - omega = nu
iphase = nu
xnorm = (two * xjp + one) * (two * xj + one)
xnu = nu + half
xnum = - xnu
!  indices in denominator of 3j correspond
!  to eq. 14 of t. orlikowski and m.h. alexander, j. chem. phys. 79, 6006
!  (1983) rather than eq. (24) of m.h. alexander,
!  j. chem. phys. 76, 3637 (1982)
x = xf3j (xjp, xlamda, xj, xnum, zero, xnu)
else
iphase = jp + j + 1 + jtot
xnorm = (two * xjp + one) * (two * xj + one) * (two * lp + one) &
* (two * l + one)
x = xf3j (xlp, xlamda, xl, zero, zero, zero)
if  (x .eq. zero) return
x = x * xf6j (xjp, xlp, xjtot, xl, xj, xlamda)
end if
if (x .eq. zero) return
x = x * xf3j (xjp, xlamda, xj, -half, zero, half)
iphase = (-1) ** iabs(iphase)
v = iphase * x * sqrt(xnorm)
return
end
! --------------------------------------------------------------------
subroutine vlmstp (jp, lp, j, l, jtot, kp, k, lambda, mu, &
	iepsp, ieps, v, csflag)
!  subroutine to calculate v-lambda matrices for close-coupled or coupled-stat
!  treatment of collisions of a symmetric top with an atom
!  the primitive cc and cs matrix elements are given in eqs. (26) and (32),
!  respectively, of s. green, j. chem. phys. 64, 3463 (1976)
!  the expressions for the full matrix elements are given in eqs. (46-48) and
!  (50) of the same article.
!  note, that in this article the bra indices (left side of matrix elements) u
!  are primed, while in the conventions of the present subroutine the bra
!  indices are primed and the ket indices (right side of matrix elements),
!  unprimed.
!
!  author:  millard alexander
!  current revision date:  27-mar-90
!  -----------------------------------------------------------------------
!  variables in call list:
!    jp:       rotational quantum number of left side of matrix element (bra)
!    lp:       orbital angular momentum of left side of matrix element (bra)
!    j:        rotational quantum number of right side of matrix element (ket)
!    l:        orbital angular momentum of right side of matrix element (ket)
!    jtot:     total angular momentum
!    kp:       k quantum number of bra
!    k:        k quantum number of ket
!    iepsp:    symmetry index of bra
!    ieps:     symmetry index of ket
!    lambda:   order of legendre term in expansion of potential
!    mu:       absolute value of index of legendre term in expansion of
!              potential
!    v:        on return, contains desired matrix element
!    csflag:   if .true., then cs calculation; in this case:
!                jtot in call list is cs lbar
!                lp is nu (cs projection index)
!                l is not used
!              if .false., then cc calculation; in this case:
!                jtot is total angular momentum and l and lp correspond
!                to the orbital angular momenta
!  subroutines called:
!     prmstp
!  -----------------------------------------------------------------------
implicit double precision (a-h,o-z)
logical csflag
data one, zero / 1.0d0, 0.0d0 /
v = zero
iphase = ieps * iepsp * (-1) ** (jp + j + lambda + mu)
if (iphase .eq. -1) return
kdif = k - kp
if (iabs(kdif) .eq. mu) then
omeg = one
if (csflag) then
!  the signed value of mu in the cs matrix elements is given by the
!  3j symbol (j' lambda j / -k' mu k) so that mu = k' - k = - kdif
musign = - kdif
!  the multiplicative factor is given by Eq. (52) of S. Green, j. chem. phys.
!  64, 3463 (1976)
if (kdif .gt. 0) omeg =  (-1) ** mu
else if (.not. csflag) then
!  the signed value of mu in the cc matrix elements is given by the
!  3j symbol (j' j lambda / k' -k mu) so that mu = k - k'= kdif
musign = kdif
!  the multiplicative factor is given by Eq. (48) of S. Green, j. chem. phys.
!  64, 3463 (1976)
if (kdif .lt. 0)  omeg =  (-1) ** mu
end if
!  contribution from (jp, kp, lp / Y(lambda, mu) / j, k, l), that is
!  the first term in Eq. (46) of S. Green, j. chem. phys. 64, 3463 (1976)
call prmstp (jp, lp, j, l, jtot, kp, k, lambda, musign, &
vprm, csflag)
v = v + omeg * vprm
end if
if (k + kp .eq. mu) then
!  n.b. for k = 0 and/or kp = 0, we recompute the same primitive matrix
!  element (here we follow MOLSCAT, although this might be somewhat inefficient
!  this is the second term in Eq. (46) of S. Green, j. chem. phys. 64, 3463
if (.not.csflag) then
!  cc contribution from (jp, -kp, lp / Y(lambda, mu) / j, k, l)
call prmstp (jp, lp, j, l, jtot, -kp, k, lambda, mu, &
  vprm, csflag)
v = v + vprm * iepsp
else if (csflag) then
!  cs contribution from (jp, kp, lp / Y(lambda, mu) / j, -k, l)
call prmstp (jp, lp, j, l, jtot, kp, -k, lambda, mu, &
  vprm, csflag)
v = v + vprm * ieps
end if
end if
return
end
! --------------------------------------------------------------------
subroutine vlmstpln (jp, lp, j, l, j2p, j2, &
	j12p, j12, jtot, kp, k, lambda, mu, &
	lambda2, mu2, &
	iepsp, ieps, v, csflag)

!  subroutine to calculate v-lambda matrices for close-coupled or
!  coupled-state treatment of collisions of a symmetric top with an atom
!  the primitive cc  elements are given in claire's thesis
!  note, that in this article the bra indices (left side of matrix
!  elements) u
!  are primed, while in the conventions of the present subroutine the
!  bra indices are primed and the ket indices (righ side of matrix
!  elements), unprimed.
!
!  author: claire rist
!  current revision date: 17-jan-1992
!
! -----------------------------------------------------------------------
!  variables in call list:
!    jp:       rotational quantum number of left side of matrix element
!              (bra) for symmetric top
!    lp:       orbital angular momentum of left side of matrix element
!              (bra)
!    j:        rotational quantum number of right side of matrix element
!              (ket) for symmetric top
!    l:        orbital angular momentum of right side of matrix element
!              (ket)
!    j2p:      rotational quantum number of left side of matrix element
!              (bra) for linear molecule
!    j2:       rotational quantum number of left side of matrix element
!              (ket) for linear molecule
!    j12p:     rotational quantum number of left side of matrix element
!              (bra) j12p = jp + j2p
!    j12:      rotational quantum number of left side of matrix element
!              (ket) j12 = j + j2
!    jtot:     total angular momentum
!    kp:       k quantum number of bra
!    k:        k quantum number of ket
!    iepsp:    symmetry index of bra
!    ieps:     symmetry index of ket
!    lambda:   order of rotation matrix term in expansion of potential
!    mu:       absolute value of index of rotation matrix term in expansion of
!              expansion of potential
!    lambda2:  order of rotation matrix term in expansion of potential
!    mu2:      absolute value of index of rotation matrix term in
!              expansion of potential
!    v:        on return, contains desired matrix element
!    csflag:   if .true., then cs calculation; in this case:
!                jtot in call list is cs lbar
!                lp is nu (cs projection index)
!                l is not used
!              if .false., then cc calculation; in this case:
!                jtot is total angular momentum and l and lp correspond
!                to the orbital angular momenta
!  subroutines called:
!     prmstp
!     f9j (Claire's routine)
!
! -----------------------------------------------------------------------
implicit double precision (a-h,o-z)
logical csflag
data one, zero, two  / 1.0d0, 0.0d0, 2.0d0 /
v = zero
iphase = ieps * iepsp * (-1) ** (jp + j + lambda + mu)
if (iphase .eq. -1) return
kdif = k - kp
if (iabs(kdif) .eq. mu) then
omeg = one
!  the signed value of mu in the cc matrix elements is given by the
!  3j symbol (j' lambda j/ k' mu -k) so that mu = k - k'= kdif
musign = kdif
if (kdif .lt. 0)  omeg =  (-1) ** mu
!  contribution from these 3j coeff (j' lambda j/k' mu -k)...
call prmstpln (jp, lp, j2p, j12p, j, l, j2, j12, &
jtot, kp, k, lambda, musign, lambda2, mu2, &
vprm, csflag)
v = v + omeg * vprm
end if
if ((k + kp) .eq. mu) then
!  n.b. for k = 0 and/or kp = 0, we recompute the same primitive matrix
!  element
!  cc contribution from  3j coeff (j' lambda j/ -k' mu, -k)
call prmstpln (jp, lp, j2p, j12p, j, l, j2, j12, &
jtot, -kp, k, lambda, mu, lambda2, mu2, &
vprm, csflag)
v = v + vprm * iepsp
end if
return
end
end module mod_hibasutil