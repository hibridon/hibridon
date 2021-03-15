*  BAr(B2Sigma+) potential
*  E. Hwang, Y.-L. Huang, P.J.Dagdigian, and M.H.Alexander, 
*  J. Chem. Phys. 98, 8484 (1993).
*
      include "common/syusr"
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(1)
      include "common/parpot"
      potnam='BAR(B2Sigma+) HWANG ET AL.'
      print *, potnam
1      print *, ' r (bohr)'
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0*219474.6d0
100   format(' v +',1pe16.8)
      goto 1
99    end
      include "common/bausr"
      include "common/ground"
* --------------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parpot"
      include "common/parbas"
      potnam='BAR(B2Sigma+) HWANG ET AL.'
      lammin(1)=2
      lammax(1)=2
      mproj(1)=0
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end

      subroutine pot (vv0, r)
*  -----------------------------------------------------------------------

*  subroutine to calculate the r-dependent coefficients in the
*  collision of a homonuclear diatomic with a structureless target
*  in units of hartree for energy and bohr for distance

*  on return:
*  vv0 contains the isotropic term (n=0) in the potential
*  the coefficients for each angular term in the coupling potential
*  [ vvl(i) for i = 1, nlam ] are returned in common block /covvl/
*  vvl(1) contains the anisotropic (n=2) term in the potential

*  variable in common block /conlam/
*    nlammx:    the maximum number of anisotropic terms
*    nlam:      the total number of angular coupling terms
*  variable in common block /covvl/
*    vvl        array to store r-dependence of each angular term in the
*               potential

*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      common /conlam/ nlam,nlammx

      include "common/parbas"
      common /covvl/ vvl(1)

*  BAr(B2Sigma+) potential.  Taken from E.Hwang, Y.-L.Huang, P.J.Dagdigian,
*  and M.H.Alexander, JCP 98, 8484 (19993).
*
*  call subroutine spline_ch2he to get pot value at distance r      
      call spline_ch2he(vsp,r)
*
*  convert to hartree
      vv0 = vsp / 219474.6d0
      vvl(1) = 0.d0
      return
      end
* ----------------------------------------------------------------------
      subroutine spline_ch2he(vsp, r)
c  
c  use spline interpolation sub from hipotutil.f
c
      implicit double precision (a-h,o-z)   
      dimension b0(100),c0(100),a0(100),rr(23),vv(23)
      data rr /3.75,4.00,4.25,4.35,4.50,4.75,5.00,5.25,5.50,5.75,
     :  6.00,6.50,7.00,7.50,8.00,9.00,10.00,11.00,12.00,13.00,
     :  14.00,15.00,16.00/
      data vv /810.74,-36.03,-334.41,-413.23,-452.43,-429.08,
     :  -370.44,-276.86,-176.34,-80.72,2.74,120.68,174.23,
     :  179.44,156.22,84.64,28.16,-1.49,-12.53,-14.00,-11.86,
     :  -8.97,-6.44/
      data ifirst /0/
*
      if (ifirst.eq.0) then
*
* determine spline coefficients
         call spline(23,rr,vv,a0,b0,c0)
         ifirst=1
*
* determine C6 coeff. for large r
         c6 = -vv(23)*rr(23)**6
*
* determine C12 coeff. for small r
         c12 = vv(1)*rr(1)**12
      endif
*
* using previously determined spline coefficients to determine potential at r
* if r lies between 3.75 and 16.00
      if (r.ge.3.75d0 .and. r.le.16.00) then
        vsp = seval(23,r,rr,vv,a0,b0,c0)
      else
        if (r.gt.16.00) then
          vsp = -c6/r**6
        else
* here for r < 3.75
          vsp = c12/r**12
        end if
      end if
*
      return
      end   
