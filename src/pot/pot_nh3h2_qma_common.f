      integer MAX_NR, MAX_NV
      parameter (MAX_NR=300, MAX_NV=300)
c
      common /stpln1/ nr, nv, rr, v_pot
      common /stpln2/ lb1, mu1, lb2, mu2
      common /stplnb/ spl_b
      common /stplnc/ spl_c
      common /stplnd/ spl_d
      integer nr, nv
      integer lb1(MAX_NV), mu1(MAX_NV), lb2(MAX_NV), mu2(MAX_NV)
      double precision rr(MAX_NR), v_pot(MAX_NR, MAX_NV),
     $     spl_b(MAX_NR, MAX_NV), spl_c(MAX_NR, MAX_NV),
     $     spl_d(MAX_NR, MAX_NV)
c
      common /covvl/ vvl
      double precision vvl(MAX_NV)
c
      common /coconv/ econv, xmconv
      double precision econv, xmconv
      data econv /219474.6315343234d0/
      data xmconv /0.0005485799094979479d0/
