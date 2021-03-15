*  System:  n-C6H14-N2 spherically averaged splined potential
*
*   from A. Jasper (jun-2015)
*   computed vdW Rmin/eps and spherically averaged repulsivze potential
*   splined repulsive and vdW potentials
*  
*   written by p. dagdigian
*   current revision date:  23-jun-2015
      include "common/syusr"
      include "common/bausr"
      include "common/ground"
* ------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      common /covvl/ vvl(15)
      include "common/parpot"
      econv=219474.6d0
      potnam='nhexane-N2 splined pot (data from Jasper)'
      print *, potnam
1     print *, 'R (bohr):'
      read (5, *, end=99) r
      call pot(vv0,r)
      if (r.le.0.d0) goto 99
*     vlm coefficient is returned in atomic units (hartree)
*     convert from atomic units for printout
      write (6, 100) vv0*econv
100   format(1pe16.8)
      goto 1
99    rr=5.d0
      dr=0.2d0
      open (unit=12,file='nhex_n2_splined_pec.txt')
      write(12,109)
109   format(' %R/bohr V0')
      do i=1,250
        call pot(vv0,rr)
        write (12,110) rr,vv0*econv
110     format(f7.2,1pe16.8)
        rr = rr + dr
      enddo
      close(12)
      end
* --------------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      character*(*) filnam
      include "common/parbas"
      include "common/parpot" 
      common /conlam/ nlam, nlammx, lamnum(2)
      common /cosysi/ nscode, isicod, nterm
      potnam='nHexane-N2 splined pot (data from Jasper)'
      npot=1
      nterm=1
      lammin(1)=2
      lammax(1)=2
      mproj(1)=0
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* --------------------------------------------------------------------------
      subroutine pot (vv0, r)

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

* author:  paul dagdigian
* latest revision date:  23-jun-2015
* ----------------------------------------------------------------------
* CH4-N2 spherically averaged potential - exp6 fit
*
      implicit double precision (a-h,o-z)
      common /conlam/ nlam,nlammx
      include "common/parbas"
      common /covvl/ vvl(1)
      dimension rr(99),vl(99)
      dimension csplin(99),vec(99)
      data ifirst /0/
*
*  99 values of R
      data rr /
     + 4.8000000e+00, 4.9000000e+00, 5.0000000e+00, 5.1000000e+00, 
     + 5.2000000e+00, 5.3000000e+00, 5.4000000e+00, 5.5000000e+00, 
     + 5.6000000e+00, 5.7000000e+00, 5.8000000e+00, 5.9000000e+00, 
     + 6.0000000e+00, 6.1000000e+00, 6.2000000e+00, 6.3000000e+00, 
     + 6.4000000e+00, 6.5000000e+00, 6.6000000e+00, 6.7000000e+00, 
     + 6.8000000e+00, 6.9000000e+00, 7.0000000e+00, 7.1000000e+00, 
     + 7.2000000e+00, 7.2002360e+00, 7.3002360e+00, 7.4002360e+00, 
     + 7.5002360e+00, 7.6002360e+00, 7.7002360e+00, 7.8002360e+00, 
     + 7.9002360e+00, 8.0002360e+00, 8.1002360e+00, 8.2002360e+00, 
     + 8.3002360e+00, 8.4002360e+00, 8.5002360e+00, 8.6002360e+00, 
     + 8.7002360e+00, 8.8002360e+00, 8.9002360e+00, 9.0002360e+00, 
     + 9.1002360e+00, 9.2002360e+00, 9.3002360e+00, 9.4002360e+00, 
     + 9.5002360e+00, 9.6002360e+00, 9.7002360e+00, 9.8002360e+00, 
     + 9.9002360e+00, 1.0000236e+01, 1.0100236e+01, 1.0200236e+01, 
     + 1.0300236e+01, 1.0400236e+01, 1.0500236e+01, 1.0600236e+01, 
     + 1.0700236e+01, 1.0800236e+01, 1.0900236e+01, 1.1000236e+01, 
     + 1.1100236e+01, 1.1200236e+01, 1.1300236e+01, 1.1400236e+01, 
     + 1.1500236e+01, 1.1600236e+01, 1.1700236e+01, 1.1800236e+01, 
     + 1.1900236e+01, 1.2000236e+01, 1.2100236e+01, 1.2200236e+01, 
     + 1.2300236e+01, 1.2400236e+01, 1.2500236e+01, 1.2600236e+01, 
     + 1.2700236e+01, 1.2800236e+01, 1.2900236e+01, 1.3000236e+01, 
     + 1.3100236e+01, 1.3200236e+01, 1.3300236e+01, 1.3400236e+01, 
     + 1.3500236e+01, 1.3600236e+01, 1.3700236e+01, 1.3800236e+01, 
     + 1.3900236e+01, 1.5000000e+01, 2.0000000e+01, 2.5000000e+01, 
     + 3.0000000e+01, 3.5000000e+01, 4.0000000e+01 /
*  splined pot (data from Jasper)
      data vl /
     + 3.7380170e+05, 3.1762753e+05, 2.6989510e+05, 2.2933581e+05, 
     + 1.9487168e+05, 1.6558676e+05, 1.4070271e+05, 1.1955819e+05, 
     + 1.0159123e+05, 8.6324301e+04, 7.3351658e+04, 6.2328518e+04, 
     + 5.2961914e+04, 4.5002903e+04, 3.8239957e+04, 3.2493333e+04, 
     + 2.7610300e+04, 2.3461080e+04, 1.9935396e+04, 1.6939545e+04, 
     + 1.4393904e+04, 1.2230817e+04, 1.0392794e+04, 8.8309860e+03, 
     + 7.5038832e+03, 7.5010000e+03, 6.2794919e+03, 5.3608060e+03, 
     + 4.4791423e+03, 3.7518151e+03, 3.0957496e+03, 2.5389786e+03, 
     + 2.0848151e+03, 1.6898464e+03, 1.3485930e+03, 1.0716413e+03, 
     + 8.5671219e+02, 6.6457665e+02, 4.9158865e+02, 3.5398607e+02, 
     + 2.4812825e+02, 1.6584489e+02, 9.8965684e+01, 3.9320314e+01, 
     + -1.6079233e+01, -5.5885798e+01, -8.6313977e+01, -1.0922695e+02, 
     + -1.2611480e+02, -1.3818075e+02, -1.4639717e+02, -1.5155037e+02, 
     + -1.5427643e+02, -1.5509000e+02, -1.5440741e+02, -1.5256534e+02, 
     + -1.4983580e+02, -1.4643832e+02, -1.4254974e+02, -1.3831225e+02, 
     + -1.3383977e+02, -1.2922327e+02, -1.2453503e+02, -1.1983210e+02, 
     + -1.1515920e+02, -1.1055097e+02, -1.0603391e+02, -1.0162791e+02, 
     + -9.7347528e+01, -9.3202994e+01, -8.9201092e+01, -8.5345833e+01, 
     + -8.1639025e+01, -7.8080733e+01, -7.4669661e+01, -7.1403453e+01, 
     + -6.8278952e+01, -6.5292395e+01, -6.2439584e+01, -5.9716018e+01, 
     + -5.7117000e+01, -5.4637721e+01, -5.2273332e+01, -5.0018996e+01, 
     + -4.7869931e+01, -4.5821442e+01, -4.3868946e+01, -4.2007991e+01, 
     + -4.0234271e+01, -3.8543628e+01, -3.6932068e+01, -3.5395756e+01, 
     + -3.3931017e+01, -2.1649139e+01, -3.8768934e+00, -1.0166368e+00, 
     + -3.4048085e-01, -1.3502498e-01, -6.0598640e-02 /  
*
* spline fit
      if (ifirst .eq. 0) then
* spline fit of the vl coefficients
         call dcopy(99,vl,1,vec,1)
*    evaluate derivative at first point
         der1=(vec(2)-vec(1))/(rr(2)-rr(1))
         call dspline(rr,vec,99,der1,0d0,csplin)
         ifirst = 1
       end if
* r^-6 fit to at R = 30 bohr 
       c6sum = vl(97)*rr(97)**6
* switching function for long-range
       switch_lr=0.5*(tanh(0.5*(r - 25.d0)) + 1.d0)
* determine splined coefficients at R=r
       call dcopy(99,vl,1,vec,1)
       call dsplint(rr,vec,csplin,99,r,vvx)
* merge with asymptotic form
       if (r.gt.25.d0) then
         vvx = (1.d0 - switch_lr)*vvx 
     +     + switch_lr*c6sum/(r**6)
       endif
* convert to hartree
       econv=1.d0/219474.6d0
       call dscal(1,econv,vvl,1)
* isotropic term 
       vv0 = vvx*econv
       vvl(1) = 0.d0
*
       return
       end
*===========================eof===============================