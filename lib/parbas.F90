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
integer :: ntv(maxtrm)
integer :: ivcol(maxvb2,maxtrm)
integer :: ivrow(maxvb2,maxtrm)

! cobspt parameters
!    lammin:   array containing minimum value of lambda for each term
!    lammax:   array containing maximum value of lambda for each term
!    mproj:    array containing the order of the reduced rotation matrix
!              elements for each term.  lammin can not be less than mproj.
!              for homonuclear molecules, the allowed values of lambda for
!              each term range from lammin to lammax in steps of 2
!
integer :: lammin(maxtrm)
integer :: lammax(maxtrm)
integer :: mproj(maxtrm)

! cobsptln
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
integer :: lam2(maxtrm)
integer :: m2proj(maxtrm)
end module
