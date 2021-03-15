*  System:  O(3P) + He
*
*   calculation of potential energy curves.  potentials taken
*   from Krems et al., j. chem. phys. 116, 1457 (2002)
*
*   assume 3P spin-orbit matrix elements do not depend on R
*  
*   written by p. dagdigian (adapted from pot_oxe_3p1d.f)
*   current revision date:  24-sep-2014
*
      include "common/syusr"
      include "common/bausr"
      include "common/ground"
* --------------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      common /coselb/ ibasty
      include "common/parbas"
      include "common/parpot"
      potnam='O(3P)-He'
      npot=1
      ibasty=22
      lammin(1)=1
      lammax(1)=19
      mproj(1)=0
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* --------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysi/ junk(5), npot
      common /cosysr/ isrcod, junkr, en(4)
      common /covvl/ vvl(19)
      include "common/parpot"
      potnam='O(3P)-He'
      econv=219474.6d0
      print *, potnam
1     print *, 'R (bohr):'
      read (5, *, end=99) r
      if (r.le.0.d0) goto 99
      call pot(vv0,r)
      write (6, 100) econv*vvl(4),econv*vvl(5),
     +  econv*vvl(13),econv*vvl(14)
100   format('  ',2(1pe16.8))
      goto 1
99    rr=2.6d0
      dr=0.25d0
      open (unit=12,file='o3phe_vlms.txt')
      write(12,209)
209   format(' %R/bohr    3Sig-   3Pi   Axy(R=inf)   Azy(R=inf)')
      do i=1,201
        call pot(vv0,rr)
        write (12,110) rr,econv*vvl(4),econv*vvl(5),
     +  econv*vvl(13),econv*vvl(14)
110     format(f7.2,4(1pe16.8))
        rr = rr + dr
      enddo
      close(12)
      end
* --------------------------------------------------------------------------
      subroutine pot (vv0, r)
*  subroutine to calculate the r-dependent coefficients
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*    vv0        (not used)
*
*  the only nonzero elements here are vvl(4,5) and vvl(13,14)
*
*    vvl(1,3)   contains the 1Sigma+, 1Pi, 1Delta energies
*    vvl(4,5)   contains the 3Sigma- and 3Pi PE energies
*    vvl(6,7)   contains Axy and Azy, the spin-orbit matrix elements
*               within the 3P state
*    vvl(8,12)  contains Bss, Byx, Bxs, Bsy, and Bxd, the spin-orbit
*               matrix elements coupling the 1D and 3P states
*    vvl(13,19) values of the Axy, Azy, Bss, Byx, Bxs, Bsy, and Bxd
*               matrix elements at R=inf
*  variable in common block /conlam/ used here
*    nlam:      the number of angular coupling terms actually used
*  variable in common block /covvl/
*    vvl:       array to store r-dependence of each
*               angular term in the potential
*  variable in common block /coconv/
*    econv:     conversion factor from cm-1 to hartrees
*  variables in common /cosysr/
*    isrcod:    number of real parameters
*    en1d:      asymptotic energy of the 1D state (cm-1)
* 
* author:  paul dagdigian
* latest revision date:  24-sep-2014
* ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      include "common/parbas"
      common /covvl/ vvl(19)
      econv=219474.6d0
*
*  set vv0 term to zero
      vv0 = 0.d0
*  clear vvl array
      call dscal(19,0.d0,vvl,1)
      sigma = pot_heo_s(r) * 1.d-6
      pi = pot_heo_p(r) * 1.d-6
      vvl(4) = sigma
      vvl(5) = pi
*
*  enter R and R-inf values for Axy and Azy
* convert to hartree
      econv=219474.6d0
      vvl(6) = 76.425669d0 / econv
      vvl(7) = 53.870478d0 / econv
      vvl(13) = 76.425669d0 / econv
      vvl(14) = 53.870478d0 / econv
      return
      end
* --------------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION POT_HEO_S(RR)    
CSYSTEM:He-O SIGMA state 
CLEVEL:RCCSD(T) 
CRMS 0.075 microEh
C units R=a.u. E=microEh
C Citation: 
C R. V. Krems, A. A. Buchachenko, 
C M. M. Szczesniak, J. Klos, and G. Chalasinski
C J. Chem. Phys. 116, 1457 (2002); 

      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION XX(15)
      DATA XX/ 
     *  0.247052945146152503D+01,
     *  0.251844051379198852D+01,
     *  0.340445104396755691D+06,
     * -0.435976248635663243D+06,
     *  0.342494546743638231D+06,
     * -0.138458782770920021D+06,
     *  0.347369265242454421D+05,
     * -0.492620366671661395D+04,
     *  0.331966775853499371D+03,
     *  0.289184192423140888E-02,
     * -0.971006114259508779D+00,
     *  0.321811869155314483D+02,
     *  0.326267696066016470D+03,
     *  0.584138476146936323D+07,
     *  0.675054074614624213D+07/             
      TERM1 = DEXP(-XX(1)*(RR-XX(2)))
      TERM2 = XX(3) + XX(4)*RR+XX(5)*RR**2+
     *        XX(6)*RR**3+XX(7)*RR**4+XX(8)*RR**5
     *        +XX(9)*RR**6+XX(10)*RR**7+XX(11)*RR**8
      TERM3 = 0.5D0*
     *   (1.D0+ DTANH(XX(12)+XX(13)*RR))
      TERM4 = XX(14)/(RR**6)+XX(15)/(RR**8)
       POT_HEO_S = TERM1*TERM2-TERM3*TERM4
       RETURN
C LAST CARD IN POT_HEO_S FUNCTION 
       END
* --------------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION POT_HEO_P(RR)    
CSYSTEM:He-O PI state 
CLEVEL:RCCSD(T) 
CRMS 0.04 microEh
C units R=a.u. E=microEh
C Citation:
C R. V. Krems, A. A. Buchachenko, 
C M. M. Szczesniak, J. Klos, and G. Chalasinski
C J. Chem. Phys. 116, 1457 (2002);
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION XX(15)
      DATA XX/ 
     *  0.282900323803789755D+01,
     *  0.300548314220765889D+01,
     *  0.237226490676133340D+06,
     * -0.412752921752802096D+06,
     *  0.327277095970970870D+06,
     * -0.144766215678601089D+06,
     *  0.394032679819154509D+05,
     * -0.647739187808830047D+04,
     *  0.599292718657259798D+03,
     * -0.241741251231964469D+02,
     * -0.157188726404088319D-01,
     *  0.321811869155314483D+02,
     *  0.326267696066016470D+03,
     *  0.503668844517287519D+07,
     *  0.176724543567391448D+08/                      
      TERM1 = DEXP(-XX(1)*(RR-XX(2)))
      TERM2 = XX(3) + XX(4)*RR+XX(5)*RR**2+
     *        XX(6)*RR**3+XX(7)*RR**4+XX(8)*RR**5
     *        +XX(9)*RR**6+XX(10)*RR**7+XX(11)*RR**8
      TERM3 = 0.5D0*
     *   (1.D0+ DTANH(XX(12)+XX(13)*RR))
      TERM4 = XX(14)/(RR**6)+XX(15)/(RR**8)
       POT_HEO_P = TERM1*TERM2-TERM3*TERM4
       RETURN
C LAST CARD IN POT_HEO_P FUNCTION 
       END
*===================================eof====================================
