!comdeck partens
! switched from jmax=14 to jmax=20 (pjd)
!.....for jmax=14 (49 channels, 29 buffers)
!      parameter (jmx=14,kmx=2*jmx+1,lmx=kmx,kkmx=3*kmx)
!      parameter (lbufs = 35525, lbuflb = 1421)
!      NOTE:  lbufs=5*5*7*7*29 and lbuflb=7*7*29
!      parameter (lbufs = 499851, lbuflb = 6171)
!      NOTE:  lbufs=9*9*11*11*51 and lbuflb=11*11*51
!
!  increase number of buffers from 41 to 51
! .....for jmax=20 (121 channels, 51 buffers)
!      parameter (lbufs = 302621, lbuflb=4961)
!      NOTE:  lbufs=11*11*41*61 and lbuflb=11*11*41
!      parameter (jmx=20,kmx=2*jmx+1,lmx=kmx,kkmx=3*kmx)
!
!...increase size of buffers - p. dagdigian, 22-jun-2010
!      parameter (jmx=26, kmx=2*jmx+1,lmx=kmx, kkmx=3*kmx)
!...increase size of buffers - p. dagdigian, 2-nov-2012
parameter (jmx=90, kmx=2*jmx+1,lmx=kmx, kkmx=3*kmx)
!
!   set lbuflb = 3*length and lbufs = 3*length(length+1)/2,
!   length = 1500
!      parameter (lbufs = 3377250, lbuflb = 4500)
!   length = 3000   (19-sep-2012   p.j.dagdigian)
parameter (lbufs = 13504500, lbuflb = 9000)
common /cosmrb/ srealp(lbufs)
common /cosmib/ simagp(lbufs)
common /coipp/  ipackp(lbuflb)
common /cojpp/  jpackp(lbuflb)
common /colpp/  lpackp(lbuflb)

