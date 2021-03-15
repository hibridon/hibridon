*system:  test routine for OH(A/X)-Kr system
*         to debug basgp1 basis routine
*
*   written by P.j.dagdigian (dec-2011)
*
      include "common/syusr"
      include "common/ground"
      include "common/bausr"
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      common /coselb/ ibasty
      common /covpot/ numvib,ivibpi(5)
      potnam='OH(A/X)-Kr TEST ROUTINE'
      ibasty=19
      numvib=2
*  consider only v=0 and 1 vim levels of 2pi state
      ivibpi(1)=0
      ivibpi(2)=1
      nterm=4 + 3*(numvib - 1)
*  for 2sigma state
      lammin(1)=0
      lammax(1)=4
      mproj(1)=0
*  for first 2pi level
      lammin(2)=0
      lammax(2)=4
      mproj(2)=0
      lammin(3)=2
      lammax(3)=4
      mproj(3)=2
      lammin(4)=1
      lammax(4)=4
      mproj(4)=1
*  for second 2pi level
      lammin(5)=0
      lammax(5)=4
      mproj(5)=0
      lammin(6)=2
      lammax(6)=4
      mproj(6)=2
      lammin(7)=1
      lammax(7)=4
      mproj(7)=1
      return
      end
* --------------------------------------------------------------------------
      subroutine driver
      implicit double precision (a-h,o-z)
      character *48 potnam
      logical csflag, ljunk, ihomo, lljunk
      include "common/parbas"
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(24)
      common /covpot/ numvib,ivibpi(5)
      potnam='OH(A/X)-Kr TEST ROUTINE'
      print *, potnam
*  consider only v=0 and 1 vib levels of 2pi state
      numvib = 2
      ivibpi(1)=0
      ivibpi(2)=1
      nterm=4 + 3*(numvib - 1)
      write (6,89) numvib, (ivibpi(i),i=1,numvpi)
89    format(' Number of 2pi vibrational levels =',i3/
     :  5x,'v =',10i3)       
      print *
1      print *, ' r (bohr) '
      read (5, *, end=93) r
      if (r .lt. 0.d0) go to 93
      call pot(vv0,r)
      write (6, 100) vvl
100   format(' vsig  ',1pe16.8/
     :  4(' vpisum',3(1pe16.8),/,
     :  ' vpi2  ',3(1pe16.8),/,
     :  ' vpi1  ',3(1pe16.8)/))
      goto 1
93    r=4

*      do i=1,100
*       call pot(vv0,r)
*       write(2,101) r,vvl
*101    format(f8.4,8(1pe16.8))
*       r=r+0.2
*      enddo

99    end
      subroutine pot (vv0, r)
* ----------------------------------------------------------------------
*  subroutine to calculate the r-dependent coefficients 
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*    vv0      dummy term here
*  variable in common block /covvl/
*    vvl:     vector to store r-dependence of each term
*             in potential expansion
*    vvl(1:x) expansion of vsigma in d(l,0) terms, for l=0,...
*    vvl(x+1:y) expansion of vsum for pi(v=0) in d(l,0) terms
*    vvl(y+1:z) expansion of vdif for pi(v=0) in d(l,2) terms
*    vvl(z+1:a) expansion of v1 coupling of sigma to pi(v=0) in d(l,1) terms
* É  and more terms if numvib > 0
*
* author:  paul dagdigian
* latest revision date:  24-jan-2012
* ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical ljunk, ihomo, csflag, lljunk
      include "common/parbas"
      common /colpar/ ljunk(5),csflag,lljunk(2),ihomo
      common /covvl/ vvl(1)
      common /covpot/ numvib,ivibpi(5)
*
*     vv0 term is a dummy
      vv0 = 0.
*
*  v for 2sigma state
      il = 0
      nlam = lammax(1) - lammin(1) + 1
      do ilam=1,nlam
        il = 0
        vvl(il) = ilam*10.
      end do
*  v's for 2pi levels
      do iv=1,numvib
        do it=1,3
          il = il + 1
          nlam = lammax(il) - lammin(il) + 1
          vvl(il) = 1000.*iv + ilam*10.
*         need to include terms for Vsum, Vdif, and V1 here
        end do
      end do
      return
      end



