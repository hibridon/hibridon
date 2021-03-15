c     Number of lambda terms used, and dimension of vvl array
c     To avoid naming confliction, use CS as prefix
      integer CSLMAX, CSNVVL
      parameter (CSLMAX = 11, CSNVVL = 2 * CSLMAX - 1)
c
c
c     Number of R's in the data file
      integer CSNR
      parameter (CSNR = 389)
c
c
c     covvl common block
c     potential expanded in D0 and D2 rotation matrices
c     vv0: v_00
c     vvl(1:CSLMAX): v_l0 (l from 1 to CSLMAX)
c     vvl(CSLMAX+1:): v_l2 (l from 2 to CSLMAX)
      common /covvl/ vvl
      double precision vvl(CSNVVL)
c
c
c     conlam common block
c     nlam: the number of angular coupling terms actually used
c     nlammx: the maximum number of angular coupling terms
c     lamnum: number of non-zero v2 matrix elements for each lambda
      common /conlam/ nlam, nlammx, lamnum
      integer nlam, nlammx, lamnum(CSNVVL)
