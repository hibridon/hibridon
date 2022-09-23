!comdeck parbas
module mod_parbas
!  revised march 1992, c.r. 13-may-1997 by mha
!  revised by Q.Ma for explicit type declaration - 16-nov-2011
!  increased maxtrm from 12 to 18 (p.dagdigian, 17-nov-2011)
!
integer, parameter :: maxtrm=18
integer, parameter :: maxvib=10
integer, parameter :: maxvb2=maxvib**2

! cobsp2 parameters
!    ntv:      Number of vibrational blocks for each term. All of these
!              use same lammin, lammax, mproj. These numbers as well
!              as the corresponding ivcol, ivrow (see below) should be
!              set in loapot and must be consistent with the pot routine
!              Otherwise unrecognized chaos!!!
!              Each block corresponds to a pair of vibrational quantum
!              numbers (ivrow,ivcol), which must correspond to lower triangle
!              of potential matrix (i.e., ivrow.ge.ivcol for
!              iterm=1 (Vsig), iterm=2 (Vpi), iterm=4 (V2) and
!              for iterm=3 (V1) ivrow=ivs, ivcol=ivp)
!    ivrow(ivb,iterm): row vibrational state for vibrational block ivb
!                      in term iterm
!    ivcol(ivb,iterm): column vibrational state for vibrational block ivb
!                      in term iterm
integer :: ntv(maxtrm) = -1
integer :: ivcol(maxvb2,maxtrm) = -1
integer :: ivrow(maxvb2,maxtrm) = -1

! cobspt parameters
!  variables in obspt must be set in loapot!!
!              Order of reduced rotation matrix d(theta1) as defined in 
!              eq (21) of reference J. Chem Phys. 98 (6), 1993 with    
!              (lambda, mproj) = (l1, m1).   
!    lammin:   array containing minimum value of lambda for each term
!    lammax:   array containing maximum value of lambda for each term
!    mproj:    array containing the order of the reduced rotation matrix
!              elements for each term.  here, lammin is greater to mproj.
!              lammin can not be less than mproj.
!              for homonuclear molecules, the allowed values of lambda for
!              each term range from lammin to lammax in steps of 2
!              the length of each of these arrays is limited to nterm
!              terms where nterm is defined in the common block /cosysi/

integer :: lammin(maxtrm) = -1
integer :: lammax(maxtrm) = -1
integer :: mproj(maxtrm) = -1

! cobsptln parameters
!              Order of reduced rotation matrix d(theta2) as defined in 
!              eq (21) of reference J. Chem Phys. 98 (6), 1993 with
!              (lam2, m2proj) = (l2, m2).      
!    lam2:     array containing the order of the reduced rotation matrix
!              elements for each term. in case of homonuclear molecule
!              is even.
!    m2proj:   array containing the order of the reduced rotation
!              matrix elements for each term.
!              here, lammin and lam2 are greater than m2proj .
!
integer :: lam2(maxtrm) = -1
integer :: m2proj(maxtrm) = -1
end module
