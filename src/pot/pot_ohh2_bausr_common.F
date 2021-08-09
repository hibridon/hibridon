C     Maximum number of R's, terms of B's and F's
C     The R grid does not need to be evenly spaced
      integer MAX_NR, MAX_NVB, MAX_NVF
      parameter (MAX_NR=400, MAX_NVB=100, MAX_NVF=70)
C     
C     Actual number of R's, terms of B's and F's, R, B, F
      common /pisg1/ nr, nvb, nvf
      common /pisg2/ rr, bcoef, fcoef
C     Arrays listing l1, l2, and l
      common /pisg3/ lam1b, lam2b, lamb, lam1f, lam2f, lamf
C     Arrays storing information for cubic spline evaluation
      common /pisgb/ splb_b, splb_c, splb_d
      common /pisgf/ splf_b, splf_c, splf_d
      integer nr, nvb, nvf
      integer lam1b(MAX_NVB), lam2b(MAX_NVB), lamb(MAX_NVB)
      integer lam1f(MAX_NVF), lam2f(MAX_NVF), lamf(MAX_NVF)
      double precision rr(MAX_NR), bcoef(MAX_NR, MAX_NVB),
     $     fcoef(MAX_NR, MAX_NVF)
      double precision splb_b(MAX_NR, MAX_NVB), splb_c(MAX_NR, MAX_NVB),
     $     splb_d(MAX_NR, MAX_NVB), splf_b(MAX_NR, MAX_NVF),
     $     splf_c(MAX_NR, MAX_NVF), splf_d(MAX_NR, MAX_NVF)
c
      common /covvl/ vvl
      double precision vvl(MAX_NVB+MAX_NVF)
c
      common /coconv/ econv, xmconv
      double precision econv, xmconv
      data econv /219474.6315343234d0/
      data xmconv /0.0005485799094979479d0/
C     Machine epsilon
      real(8), parameter :: machep=epsilon(0d0)
