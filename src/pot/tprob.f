*  ---------------------------------------------------------
      subroutine tprob(nmax, nch, nphoto, fluxn,flux_ja)
*  ---------------------------------------------------------

* subroutine to sort normalized photodissociation fluxes by atomic angular momentum is
* chalcogen hydrides (SH for example, basistyp=26)

* currewnt revision date:  13-oct-2018 by m. alexander

* variables in call list:
*   nmax:     maximum row dimension of matrices
*   nch:      number of states (should not exceed 9 here)
*   nphoto:   number of initial vibrational states (should be 1 here)
*   fluxn:    nphoto x nmax matrix of normalized photodissociation fluxes
*   flux_ja:  on return, total photodissociation flux in Ja=0, 1, and 2

*  ---------------------------------------------------------
*  variables in common cojtot
*     jjtot     total angular momentum jtot
*     jjlpar    parity parameter jlpar
*  variable in common block /coj12/
*    j12:      array containing vector sum of ja + j (similar
*              here j12 is lower-case j in alexander, dagdigian, klos notes
*  variable in commmon block /coja/
*     jja:      atomic electronic angular momentum for each channel for 3P+2S scattering
*  variable in common block /coel/
*     ll:       orbital angular momentum for each channel for 3P+2S scattering
* --------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical tprob_tst
      common /cojtot/ jjtot, jjlpar
      common /coj12/  j12(9)
      common /coja/  jja(9)
      common /coel/  ll(9)
      dimension fluxn(nmax,nmax)
      dimension flux_ja(3)


* initialize sorted fluxes
      flux0=0d0;
      flux1=0d0;
      flux2=0d0;
* sort normalized fluxes
      do i=1,nch
         if (ja(i).eq.0) flux0=flux0+fluxn(1,i)
         if (ja(i).eq.1) flux1=flux1+fluxn(1,i)
         if (ja(i).eq.2) flux2=flux2+fluxn(1,i)
      enddo
         flux_ja(1)=flux0
         flux_ja(2)=flux1
         flux_ja(3)=flux2

      write (9,30) flux0, flux, flux2
30    format(\,'** PHOTOFRAMENT FLUXES INTO V= 0, 1, 2: ',3f12.5)
      return
      end

