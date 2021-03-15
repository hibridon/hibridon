*Ar-OH(A^2Sigma+) PES of Fawzy and Heaven
* References: W. Fawzy and M. C. Heaven, J. Chem. Phys. 89, 7030 (1988)
* U. Schnupf, J. M. Bowman, and M. C. Heaven, Chem. Phys. Lett. 189,487 (1992)

      subroutine driver
      implicit double precision (a-h,o-z)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(10)
      include "common/parpot"
      potnam='Ar-OH(A^2Sigma+) PES of Fawzy and Heaven'
      print *, potnam
1      print *, ' r (bohr)'
      rshift=0.5
      xfact=0.8
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vsum',/,7(1pe16.8))
      goto 1
99    end

      include "common/syusr"
      include "common/bausr"
      include "common/ground"
      subroutine loapot(iunit,filnam)
* ------------------------------------------------------------------------
      character*(*) filnam
      include "common/parbas"
      include "common/parpot"
      potnam='Ar-OH(A^2Sigma+) PES of Fawzy and Heaven'
      lammin(1)=1
      lammax(1)=10
      mproj(1)=0
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* ----------------------------------------------------------------------
      subroutine pot (vv0, r)
*  subroutine to calculate the r-dependent coefficients
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*  vv0        contains isotropic term (d00 term in vsum)
*  variable in common block /covvl/
*    vvl:     vector of length 6 to store r-dependence of each term
*             in potential expansion
*    vvl(1-10) expansion coefficients in dl0 (l=1:11) of vsum

* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  8-oct-1993
* revised for He-NO(X) : 1-20-95 by Moonbong Yang
* revised for CCSD(T) PES: 2002.10.13 by Jacek Klos
*revised for Ar-OH(A^2Sigma+) 2005 J.Klos
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension xlam1(7),xlam2(7),r0(7),
     :          vsum(11),xsum(11),
     :          d0(121),aa(121),thta(11),cthta(11)
      dimension kpvt(11),qraux(11),work(121),rsd(11),re(14)

      common /covvl/ vvl(10)
      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
* for distances beyond rmax difference potential is damped
      data rmax /13d0/
* coefficicients for d0 rotation matrices
* stored (by column) for each of 7 angles and for l=0:6
* angles are 0 20 40 60 80 90 100 120 140 160 180
      data d0/
     :  1.d0,!l=0
     :  1.d0,
     :  1.d0,
     :  1.d0,
     :  1.d0,
     :  1.d0,
     :  1.d0,
     :  1.d0,
     :  1.d0,
     :  1.d0,
     :  1.d0,!l=0
     :                    1.d0, !l=1
     :    0.939692620785908d0,
     :    0.766044443118978d0,
     :                  0.5d0,
     :     0.17364817766693d0,
     :                  1d-16,
     :    -0.17364817766693d0,
     :                 -0.5d0,
     :   -0.766044443118978d0,
     :   -0.939692620785908d0,
     :                   -1d0, !l=1
     :                    1d0, !l=2
     :    0.824533332339233d0,
     :    0.380236133250197d0,
     :               -0.125d0,
     :   -0.454769465589431d0,
     :                 -0.5d0,
     :   -0.454769465589431d0,
     :               -0.125d0,
     :    0.380236133250197d0,
     :    0.824533332339233d0,
     :                    1d0, !l=2
     :                    1d0, !l=3
     :    0.664884732794716d0,
     :  -0.0252333338303835d0,
     :              -0.4375d0,
     :   -0.247381933374901d0,
     :                  1d-16,
     :    0.247381933374901d0,
     :               0.4375d0,
     :   0.0252333338303835d0,
     :   -0.664884732794716d0,
     :                   -1d0, !l=3
     : 1d0,                   !l=4
     :    0.474977735636283d0,
     :   -0.319004346471378d0,
     :           -0.2890625d0,
     :    0.265901610835095d0,
     :                0.375d0,
     :    0.265901610835095d0,
     :           -0.2890625d0,
     :   -0.319004346471378d0,
     :    0.474977735636283d0,
     :                    1d0,!l=4
     :  1.d0,                 !l=5
     :    0.271491745551255d0,
     :   -0.419682045437054d0,
     :           0.08984375d0,
     :    0.281017540988309d0,
     :                  1d-16,
     :   -0.281017540988309d0,
     :          -0.08984375d0,
     :    0.419682045437054d0,
     :   -0.271491745551255d0,
     :                   -1d0,!l=5
     :                   1.d0,!l=6
     :   0.0719030017842305d0,
     :   -0.323570725710931d0,
     :         0.3232421875d0,
     :   -0.132121338573299d0,
     :              -0.3125d0,
     :   -0.132121338573299d0,
     :         0.3232421875d0,
     :   -0.323570725710931d0,
     :   0.0719030017842305d0,
     :                    1d0, !l=6
     :                    1d0, !l=7
     :   -0.107226158692938d0,
     :   -0.100601708629502d0,
     :        0.22314453125d0,
     :   -0.283479918813435d0,
     :                  1d-16,
     :    0.283479918813435d0,
     :       -0.22314453125d0,
     :    0.100601708629502d0,
     :    0.107226158692938d0,
     :                   -1d0,!l=7
     :                    1d0,!l=8
     :   -0.251839432959275d0,
     :    0.138626797752243d0,
     :   -0.073638916015625d0,
     :   0.0233078500507821d0,
     :            0.2734375d0,
     :   0.0233078500507821d0,
     :   -0.073638916015625d0,
     :    0.138626797752243d0,
     :   -0.251839432959275d0,
     :                    1d0,!l=8
     :                    1d0, !l=9
     :   -0.351696543958561d0,
     :    0.290012951832139d0,
     :   -0.267898559570312d0,
     :    0.259627174131175d0,
     :                  1d-16,
     :   -0.259627174131175d0,
     :    0.267898559570312d0,
     :   -0.290012951832139d0,
     :    0.351696543958561d0,
     :                   -1d0, !l=9
     :                    1d0, !l=10
     :   -0.401269139852809d0,
     :    0.297345221371711d0,
     :   -0.188228607177734d0,
     :   0.0646821277096134d0,
     :          -0.24609375d0,
     :   0.0646821277096134d0,
     :   -0.188228607177734d0,
     :    0.297345221371711d0,
     :   -0.401269139852809d0,
     :                    1d0/!l=10
      thta(1)=0.d0
      thta(2)=20.d0
      thta(3)=40.d0
      thta(4)=60.d0
      thta(5)=80.d0
      thta(6)=90.d0
      thta(7)=100.d0
      thta(8)=120.d0
      thta(9)=140.d0
      thta(10)=160.d0
      thta(11)=180.d0
c      do i=1,11
c      cthta(i)=dcos(thta(i)*dacos(-1.d0)/180.d0)
c      enddo
      call prearoh
      do 100 i=1,11
        vsum(i)=pes_fh(r,thta(i))
100   continue
* solve simultaneous equations for solutions
      tol=1.e-10
      call dcopy(121,d0,1,aa,1)
      call dqrank(aa,11,11,11,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,11,11,11,kr,vsum,xsum,rsd,kpvt,qraux)
* convert to hartree
      conv=1.d0/219474.6
      call dscal(11,conv,xsum,1)
      vv0=xsum(1)
      call dcopy(10,xsum(2),1,vvl,1)
      end

c SYSTEM:Ar-OH(A)
c PES FROM PAPER BY  FAWZY AND
c MICHAEL HEAVEN J. Chem. Phys  94, 2226 (1990)
c Rbig=intermolecular distance in bohr
c THETA=Jacobi angle in degrees
c OUTPUT: cm-1

         function pes_fh(rbig,theta)
         implicit double precision (a-h, o-z)
         pi=dacos(-1.d0)
         gama=dabs(theta)*pi/180.d0
         call pes(rbig,1.9126D0,gama,pes_fh,vmorse)
         return
         end

          subroutine prearoh
           common/dnab/nab,ntrans
c           this must be called before calls to pes
           nab=1
           ntrans=0
           return
           end   
            subroutine pes (rt1,rt2,gamma,v,vmorse) 
***********************************************************
c oPPPtential surface for ArOH    
c eVVVrsion from 6/29/91, 13:25 
c orrrutinen by U.Schnupf, Emory University, Atlanta, GA 
***********************************************************
c    
c iLLLt : a) J.M.Bowman,P.Schafer,J.Mol.Struc., 224(1990)133-139.
c         b) J.M.Bowman,B.Gazdy,P.Schafer,M.C.Heaven,J.P.C. 94
c            (1990) 2226-2229.
c         c) U.Schnupf,J.M.Bowman,M.C.Heaven, unpublished. 
c     
c tnnnrans : integer, if ntrans=0 then ArOH, if ntrans=1 then ArOD
c            coordinates transformation 
c     
c accclls routine : pes1 
c     
        implicit real*8    (a-h,o-z)
        implicit integer (i-n) 
c     
       logical first          
c     
c*******common/dnab/ block from the calling program ! note you have to
c       give a value for ntrans !!!
c    
       common/dnab/   nab,ntrans   
c    
       common/daroh1/b1,ex2,alp1,alp2,r0,rpi,a,ex1,al01,umh,umd 
       common/daroh2/al02,sp,als,rs,dho,alho,r20,rone,drsw,d00,delpi 
       common/daroh3/d0p1,d0p2,dsp1,dsp2,als1,als2,rthr,rtwo,ds0,dpi0 
       common/daroh4/p1,p2,p3,p4 
       common/daroh5/t1,t2,t3,t4  
c     
       data rsw1,rsw2,gsw1,gsw2/3.0,10.0,0.0,3.141592653589793/   
       data fal1,fal2,fal3,fal4/8.055,14.175,9.185,2.065/ 
       data wntoau/219474.7/ 
       data first/.true./         
       data rmin,dal02/4.6,0.01 /   
c     
       if (first) then
          first=.false.
          if (ntrans.eq.1) write(6,10) 
          call pes1  
       end if
c     
c*******transformation of the coordinates      
c       convert from ArOD jacobi coordinates to internal
c    
        IF (NTRANS.EQ.1) THEN 
           r2 = rt2   
           f  = umd*r2
           fm = r2 - f
           gam= gamma  
           r1 = sqrt(rt1*rt1+f*f-2.0*rt1*f* cos(gam))
           r3 = sqrt(fm*fm+rt1*rt1+2.0*fm*rt1* cos(gam))
           che= (r1*r1+r2*r2-r3*r3)/(2.0*r1*r2)
           if (che.lt.-1.0) che=-1.0
           if (che.gt. 1.0) che= 1.0   
c     
c*******convert from internal coordinates to jacobi coordinates
c       for ArOH, watch out the masses !   
c    
          f     = r2*umh              
          cr2   = r2
          cr1   = sqrt(r1*r1+f*f-2.0*f*r1*che)
          cosga = (cr1*cr1+f*f-r1*r1)/(2.0*f*cr1)
          if (cosga .lt.-1.0) cosga = -1.0
          if (cosga .gt. 1.0) cosga =  1.0
          ga = acos (cosga)
          r1 = cr1 
       ELSE
          r1=rt1
          r2=rt2
          ga=gamma 
       END IF 
c     
c*******linear fit for d0 and ds parameter      
c     
       d0= d0p1  +  d0p2*(r2-r20)    
       ds= dsp1  +  dsp2*(r2-r20)    
c     
c*******r1 is the distance  of Ar to the cm of OH, r2 is roh, and ga is angle
c       between such that ga = 0 is ArHO.
c     
c*******dis is "d" in the JPC paper
c     
       dis = d0*d00                          
       dss = ds*ds0  
       dpi = (d0-delpi)*dpi0     
c     
                         xal = (r1 - rsw1)/drsw           
       if (r1 .gt. rsw2) xal = 1.0
       if (r1 .lt. rsw1) xal = 0.0
c     
       fal = fal1*xal**2-fal2*xal**3+fal3*xal**4-fal4*xal**5 
c     
       first_a = -2.0*exp(-als1*(r1 - rs))
       sec = exp(-als2*(r1 - rthr))
       third = first_a + sec
       vs = dss*third
c     
       if (ga .le. sp) then 
c     
       if ( r1.le.rmin) then
       xr = r1/rmin
       fr=1 - t1*xr*xr - t2*xr**3 - t3*xr**4 - t4*xr**5  
       al02=al02+dal02*(1-exp(-(r2-r20)))*fr
       endif  
c     
       v0 = dis*(-2.0*exp (-al01*(r1 - r0)) + exp (-al02*(r1 - rone)))
       x0 = (exp ((ga - gsw1)*ex1) - 1.0)/(exp ((sp - gsw1)*ex1) - 1.0)
       f0 = p1*(x0*x0) + p2*(x0**3) + p3*(x0**4) + p4*(x0**5) 
       v = (1.0 - f0)*v0 + f0*vs
c     
       else  
c    
       vp = dpi*(-2.0*exp (-alp1*(r1 - rpi)) + exp (-alp2*(r1 - rtwo)))
       xp = abs ((ga - sp)/(gsw2 - sp))**ex2
       fp = p1*(xp*xp) + p2*(xp**3) + p3*(xp**4) + p4*(xp**5) 
       v = (1.0 - fp)*vs + fp*vp
c     
       end if  
c     
c*******combine ArOH and morse potentential for OH  
c     
       vmorse=dho*(1.d0 - dexp(-alho*(r2-r20)))**2
c commented below by J. K. to leave only interaction part in cm-1
c       v = v/wntoau +vmorse 
        
c    
10     format(' I will do coordination transformation between ArOH/ArOD')
c    
       return
C           END OF ROUTINE PES  
       end
c~~~~~~~~~~~~~~~~~~~~~~
       subroutine pes1  
**********************************
c  IIIbelong to routine pes
c     
c accclls routine : pes2
c     
**********************************
       implicit real*8    (a-h,o-z)   
       implicit integer (i-n) 
c     
       logical first 
c     
       common/dnab/ nab, ntrans 
c    
       common/daroh1/b1,ex2,alp1,alp2,r0,rpi,a,ex1,al01,umh,umd 
       common/daroh2/al02,sp,als,rs,dho,alho,r20,rone,drsw,d00,delpi  
       common/daroh3/d0p1,d0p2,dsp1,dsp2,als1,als2,rthr,rtwo,ds0,dpi0 
       common/daroh4/p1,p2,p3,p4 
c    
       data first/.true./ 
       data ua/72820.69/,uh/1837.47/,uc/29164.37/,ud/3671.93/ 
       data rsw1,rsw2,gsw1,gsw2/3.0,10.0,0.0,3.141592653589793/   
c     
       data b1,ex2,alp1,alp2 /-1.0000, 1.00000, 2.50000, 0.86250/
       data r0,rpi,a,ex1     /5.29800, 4.10300, 9.69690, 1.39520/
       data al01,al02,sp,als1/1.09569, 2.98290, 1.37706, 0.90604/
       data rs,d0p1,d0p2,als2/7.17320, 1036.755,681.8721,1.78000/
       data dsp1,dsp2,delpi  /59.40297,175.474 ,123.0000        / 
       data dho,alho,r20     /0.123888,1.20967, 1.912600        /
       data ar,br / 10.0,-1./  
c     
c*******routine pes2 will be used for readin parameters from the calling
c       program, otherwise the data-parameter will be used ! 
c    
cu sss  call pes2 
c    
c*******write out parameters depending on nab=2
c     
       if(first) then
          if(nab.eq.2) then
             write(6,*) 
             write(6,*)' Potential parameters : (routine pes)'
             write(6,*)' b1,ex2,alp1,alp2 :',b1,ex2,alp1,alp2 
             write(6,*)' r0,rpi,a,ex1     :',r0,rpi,a,ex1
             write(6,*)' al01,al02,sp,als1:',al01,al02,sp,als1
             write(6,*)' rs,d0p1,d0p2,als2:',rs,d0p1,d0p2,als2 
             write(6,*)' dsp1,dsp2,delpi  :',dsp1,dsp2,delpi 
             write(6,*)' dho,alho,r20     :',dho,alho,r20 
             write(6,*)
           end if
           first=.false.
        end if  
c     
c*******Calculation of constant terms 
c     
       rone = r0 + (1.0/al02)*log (2.0*al01/al02)
       rtwo = rpi+ (1.0/alp2)*log (2.0*alp1/alp2)
       rthr = rs + (1.0/als2)*log (2.0*als1/als2)
c    
       d00 = 0.5/(1.0 - al01/al02) 
       ds0 = 0.5/(1.0 - als1/als2) 
       dpi0= 0.5/(1.0 - alp1/alp2) 
c    
       p1=0.5*a
            p2=(-1.5*a + 0.5*b1 + 10.0)
            p3=(1.5*a - b1 - 15.0)
            p4=(0.5*(b1 - a) + 6.0) 
c    
       T1 = 0.5*ar
           t2 = (-1.5*ar + 0.5*br + 10.0)
           t3 = (1.5*ar - br - 15.0)
           t4 = (0.5*(br - a ) + 6.0) 
       drsw = (rsw2 - rsw1)   
       umh  = uc/(uh + uc)   
       umd  = uc/(ud + uc)   
c     
       return
C           END OF ROUTINE PES1 
       end 
c~~~~~~~~~~~~~~~~~~~~~~~~~
	  subroutine pes2
********************************
c I belong to routine pes1 
********************************
c   
            implicit real*8    (a-h,o-z)
            implicit integer (i-n) 
c  
       parameter(npar=50)
c  
       common/dpar/ potpai(npar),potpaf(npar),potpar(npar),nloopi,nrand 
c   
       common/daroh1/b1,ex2,alp1,alp2,r0,rpi,a,ex1,al01,umh,umd 
       common/daroh2/al02,sp,als,rs,dho,alho,r20,rone,drsw,d00 
       common/daroh3/d0p1,d0p2,dsp1,dsp2,als1,als2,rthr,rtwo,ds0,dpi0 
       common/daroh4/p1,p2,p3,p4 
c 
c*****Set the parameter            
c
       if (b1.Eq.potpar(1)) then 
          write(6,*) 'OK-b1' 
        else
          write(6,*) 'b1',b1,potpar(1)  
       end if 
       b1    =potpar(1)
       if (ex2.Eq.potpar(3)) then 
          write(6,*) 'OK-ex2'   
        else
          write(6,*) 'ex2',ex2,potpar(3)
       end if 
       ex2   =potpar(3)
       if (alp1.Eq.potpar(5)) then 
          write(6,*) 'OK-alp1' 
        else
          write(6,*) 'alp1',alp1,potpar(5) 
       end if 
       alp1  =potpar(5)
       if (alp2.Eq.potpar(6)) then 
          write(6,*) 'OK-alp2' 
        else
          write(6,*) 'alp2',alp2,potpar(6) 
       end if 
       alp2  =potpar(6)
       if (r0.Eq.potpar(7)) then 
           write(6,*) 'OK-r0' 
        else
          write(6,*) 'r0',r0,potpar(7) 
       end if 
       r0    =potpar(7)
       if (rpI.Eq.potpar(9)) then 
           write(6,*) 'OK-rpi'
        else
          write(6,*) 'rpi',rpi,potpar(9)  
       end if 
       rpi   =potpar(9)
       if (a .Eq.potpar(10)) then 
           write(6,*) 'OK-a ' 
        else
          write(6,*) 'a',a,potpar(10) 
       end if 
       a     =potpar(10)
       if (ex1.Eq.potpar(11)) then 
           write(6,*) 'OK-ex1' 
        else
          write(6,*) 'ex1',ex1,potpar(11) 
       end if 
       ex1   =potpar(11)
       if (al01.Eq.potpar(13)) then 
           write(6,*) 'OK-al01' 
        else
          write(6,*) 'al01',al01,potpar(13) 
       end if 
       al01  =potpar(13)
       if (al02.Eq.potpar(14)) then 
           write(6,*) 'OK-al02' 
        else
          write(6,*) 'al02',al02,potpar(14) 
       end if 
       al02  =potpar(14)
       if (sp.Eq.potpar(15)) then 
           write(6,*) 'OK-sp' 
        else
          write(6,*) 'sp',sp,potpar(15) 
       end if 
       sp    =potpar(15)
       if (als1.Eq.potpar(17)) then 
           write(6,*) 'OK-als1'
        else
          write(6,*) 'als1',als1,potpar(17) 
       end if 
       als1  =potpar(17)
       if (rs.Eq.potpar(18)) then 
           write(6,*) 'OK-rs' 
        else
          write(6,*) 'rs',rs,potpar(18) 
       end if 
       rs    =potpar(18)
       if (d0p1.Eq.potpar(19)) then 
           write(6,*) 'OK-d0p1'
        else
          write(6,*) 'd0p1',d0p1,potpar(19) 
       end if 
       d0p1  =potpar(19)
       if (d0p2.Eq.potpar(20)) then 
           write(6,*) 'OK-d0p2'
        else
          write(6,*) 'd0p2',d0p2,potpar(20) 
       end if 
       d0p2  =potpar(20) 
       if (dhO.Eq.potpar(21)) then 
           write(6,*) 'OK-dho' 
        else
          write(6,*) 'dho',dho,potpar(21) 
       end if 
       dho   =potpar(21)
       if (alho.Eq.potpar(22)) then 
           write(6,*) 'OK-alho'
        else
          write(6,*) 'alho',alho,potpar(22) 
       end if 
       alho  =potpar(22)
       if (r20.Eq.potpar(23)) then 
           write(6,*) 'OK-r20' 
        else
          write(6,*) 'r20',r20,potpar(23) 
       end if 
       r20   =potpar(23)
       if (dsp1.Eq.potpar(24)) then 
           write(6,*) 'OK-dsp1' 
        else
          write(6,*) 'dsp1',dsp1,potpar(24) 
       end if 
       dsp1  =potpar(24) 
       if (dsp2.Eq.potpar(25)) then 
           write(6,*) 'OK-dsp2'
        else
          write(6,*) 'dsp2',dsp2,potpar(25) 
       end if 
       dsp2  =potpar(25) 
       if (delpi.Eq.potpar(26)) then 
           write(6,*) 'OK-delpi' 
        else
          write(6,*) 'delpi',delpi,potpar(26) 
       end if 
       als2  =potpar(30)
       if (als2.Eq.potpar(30)) then 
           write(6,*) 'OK-als2' 
        else
          write(6,*) 'als2',als2,potpar(30) 
       end if 
       als2  =potpar(30)
c
	  return
C         END OF ROUTINE PES2 
	  end 
c
c------------------------------end v--------------------------------
c   
      subroutine dpes(r,v)
      implicit double precision (a-h,o-z)
      dimension r(3)     
      return
      end
