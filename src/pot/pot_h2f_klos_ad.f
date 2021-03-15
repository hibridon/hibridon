*system: F-H2 using adiabatic klos potentials
*Reference
*J. Klos, G. Chalasinski and M. M. Szczesniak 
* Int. J. Quant. Chem., 90, 1038 (2002)
      include "common/syusr"
      include "common/bausr"
      include "common/ground"
* --------------------------------------------------------------------------
      subroutine driver
* --------------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      character *48 potnam
      common /capwern/ ipot
      common /covvl/ vvl(3)
      if (ipot.eq.1) then
         potnam='FIRST ADIABATIC KLOS F+H2'
      elseif (ipot.eq.2) then
         potnam='SECOND ADIABATIC KLOS F+H2'
      elseif (ipot.eq.3) then
         potnam='THIRD ADIABATIC KLOS F+H2'
      endif
      print *, potnam
1      print *, ' r (bohr)'
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vsum',/,7(1pe16.8))
      goto 1
99    end
* --------------------------------------------------------------------------
      block data cwad
      common /capwern/ ipot
      data ipot /3/
      end
* --------------------------------------------------------------------------
      subroutine loapot(iunit,filnam)
* --------------------------------------------------------------------------
      character*(*) filnam
      include "common/parpot"
      include "common/parbas"
      common /coselb/ ibasty
      common /capwern/ ipot
      if (ipot.eq.1) then
         potnam='FIRST ADIABATIC KLOS F+H2'
      elseif (ipot.eq.2) then
         potnam='SECOND ADIABATIC KLOS F+H2'
      elseif (ipot.eq.3) then
         potnam='THIRD ADIABATIC KLOS F+H2'
      endif
      ibasty=1
      lammin(1)=2
      lammax(1)=6
      mproj(1)=0
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      return
      end
* --------------------------------------------------------------------------
      subroutine pot (vv0, r)
*  -----------------------------------------------------------------------

*  subroutine to calculate the r-dependent coefficients in the
*  collision of a homonuclear diatomic with a structureless target
*  in units of hartree for energy and bohr for distance

*  on return:
*  vv0 contains the isotropic term (n=0) in the potential
*  the coefficients for each angular term in the coupling potential
*  [ vvl(i) for i = 1, nxl ] are returned in common block /covvl/
*  vvl(1:3) contains the anisotropic (n=2:2:6) terms in the potential

*  variable in common block /conxl/
*    nxlmx:    the maximum number of anisotropic terms
*    nxl:      the total number of angular coupling terms
*  variable in common block /covvl/
*    vvl        array to store r-dependence of each angular term in the
*               potential

*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension 
     : vau(6),
     : d0(16),
     :          vs(4),aa(16),xs(4),kpvt(4),qraux(4),work(50),rsd(4)
     
      common /conxl/ nxl,nxlmx
      common /capwern/ ipot

      include "common/parbas"
      common /covvl/ vvl(3)


*  adiabatic alexander-start-werner pes fitted to legendre expansion at r=1.4 bohr
* angles (not needed here, but included for completeness)
*      data beta / 0d0, 30d0, 60d0, 90d0/
      data zero, one, half /0.d0, 1.d0, 0.5d0/
* coefficients for d0 rotation matrices
* stored (by column) for each of 4 angles and for l=0,2,4,6
      data d0 /
     : 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 6.25d-1, -1.25d-1, -5.d-1,
     : 1.d0, 2.3437511d-2, -2.8906251d-1, 3.75d-1, 1.d0, -3.7402343d-1,
     : 3.2324218d-1, -3.125d-1/
      data iflag /1/
* determine potentials at angles
      do  i=1,4
        thjac=(i-1)*30d0
        call poth2f(r,rh2,thjac,vau,iflag)
           vs(i)=vau(ipot)
      enddo
* solve simultaneous equations for legendre expansion coefficients
      tol=1.e-10
      call dcopy(16,d0,1,aa,1)
      call dqrank(aa,4,4,4,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,4,4,4,kr,vs,xs,rsd,kpvt,qraux)
      vv0=xs(1)
      call dcopy(3,xs(2),1,vvl,1)
      vv0=vv0
      return
      end
c----------------------------------------------------------------------
      subroutine poth2f(r,rh2,th,vau,iflag)
      double precision au2ang,rh2,r,thjac,vau(3),rrh2
      complex hdmatr,work
      dimension hdmatr(3,3), work(5)
      real rwork(7),w(3), cmtohar
      character*1 jobz,uplo
      integer m,nh,lda,lwork,info
      real rhh,pi,DSO
      au2ang=0.529177249d0
      cmtohar=219474.6
      xmili2cm=1.0/4.556335
      jobz='N'
      uplo='U'
      nh=3
      lwork=5
      lda=nh
      pi=acos(-1.0)
      rrh2=1.448361259d0*au2ang
* NOTE, the rh2 distance is unchanged
      DSO=134.7
      hdmatr(1,1)=CMPLX(POT3DH11(r*au2ang,rrh2,th),0.0)
      hdmatr(1,2)=CMPLX(-POT3DH12(r*au2ang,rrh2,th)/sqrt(2.0),
     *    -sqrt(2.0)*DSO) 
      hdmatr(1,3)=CMPLX(POT3DH12(r*au2ang,rrh2,th)
     *       /sqrt(2.0),0.0)
      hdmatr(2,2)=CMPLX(POT3DVSUM(r*au2ang,rrh2,th)+DSO,0.0)
      hdmatr(2,3)=CMPLX(POT3DVDIF(r*au2ang,rrh2,th),0.0)
      hdmatr(3,3)=CMPLX(POT3DVSUM(r*au2ang,rrh2,th)-DSO,0.0)
      call cheev(jobz,uplo,nh,hdmatr,lda,w,work,lwork,rwork,info)
      if(info.ne.0) then
      write(6,*) 'ERROR in CHEEV:ABORT',info
      stop
      endif
      vau(1)=(w(1)+DSO)/cmtohar
      vau(2)=(w(2)+DSO)/cmtohar
      vau(3)=(w(3)-2.*DSO)/cmtohar
200   continue
      end 
      
C******************************************
C 3-D UCCSD(T) MODEL OF  F-H2 vdW SYSTEM
C   DIABATIC PESs 
C r, R in Angstroms
C T in degrees
C re(H-H)=0.74085A
C-----------------------------------------
C AUTHOR: Jacek Klos*) ,Grzegorz Chalasinski, 
C         Maria Szczesniak
C University of Warsaw
C Department of Chemistry
C ul. Pasteura 1 02-093
C Warsaw, POLAND
C jakl@tiger.chem.uw.edu.pl 
C *) Present address:
C Institute of Theoretical Chemistry 
C University of Nijmegen
C Toernooiveld 1 
C 6525 ED Nijmegen, The Netherlands
C------------------------------------------
C POT3DH11 - 3D H11 MODEL DIABAT
C POT3DH22 - 3D H22 MODEL DIABAT
C POT3DV3  - 3D H33 (A") MODEL DIABAT
C POT3DH12 - 3D H12 MRCI COUPLING ELEMENT 
C POT3DVSUM - 3D 0.5(H33+H22) MODEL DIABAT
C POT3DVDIF - 3D 0.5(H33-H22) MODEL DIABAT
C POT3DV1A - 1A' 3-D ADIABAT
C POT3DV2A - 2A' 3-D ADIABAT 
C*****************************************
       FUNCTION POT3DH11(RR,rsmall,T)
c     Vsigma diabat (Hxx or H11)
c     RR ,rsmall,in A, T in degrees
c     output in cm-1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PI=DACOS(-1.D0)
      XH2CM=1.d0/4.556335d0
      POT3DH11=(POT2DSG(RR,rsmall)*(DCOS(T*PI/180.D0))**2+
     * POT2DA1(RR,rsmall)*(DSIN(T*PI/180.D0))**2)*XH2CM
      RETURN
      END
C----------------------------------------
      FUNCTION POT3DH22(RR,rsmall,T)
c    Vpi diabat (Hyy or H22)
c     RR, rsmall in A, T in degrees
c     output in cm-1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PI=DACOS(-1.D0)
      XH2CM=1.d0/4.556335d0
      POT3DH22=(POT2DPI(RR,rsmall)*(DCOS(T*PI/180.D0))**2+
     * POT2DB2(RR,rsmall)*(DSIN(T*PI/180.D0))**2)*XH2CM
      RETURN
      END
C---------------------------------------
      FUNCTION POT3DVSUM(RR,rsmall,T)
c    Vsum=1/2(H33+H22)
c     RR, rsmall in A, T in degrees
c     output in cm-1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      POT3DVSUM=0.5d0*(POT3DV3(RR,rsmall,T)+POT3DH22(RR,rsmall,T))
      RETURN
      END
C--------------------------------------
      FUNCTION POT3DVDIF(RR,rsmall,T)
c    Vdif=1/2(H33-H22)
c     RR ,rsmall in A, T in degrees
c     output in cm-1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      POT3DVDIF=0.5d0*(POT3DV3(RR,rsmall,T)-POT3DH22(RR,rsmall,T))
      RETURN
      END
C-------------------------------------
      FUNCTION POT3DV3(RR,rsmall,T)
c    V3 diabat (H33, A")
c     RR,rsmall in A, T in degrees
c     output in cm-1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PI=DACOS(-1.D0)
      XH2CM=1.d0/4.556335d0
      POT3DV3=(POT2DPI(RR,rsmall)*(DCOS(T*PI/180.D0))**2+
     * POT2DB1(RR,rsmall)*(DSIN(T*PI/180.D0))**2)*XH2CM
      RETURN
      END
C******************************************************

C     3-D ADIABATIC F-H2 PESs

C******************************************************
      FUNCTION POT3DV1A(RR,rsmall,T)
c    V1A=1Aprim
c     RR in A, T in degrees
c     output in cm-1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      TERM1=(POT3DH11(RR,rsmall,T)+POT3DH22(RR,rsmall,T))*0.5D0
      TERM2=(POT3DH11(RR,rsmall,T)-POT3DH22(RR,rsmall,T))*0.5D0
      IF (RR.GE.2.0D0) THEN
      ROOT1=TERM1-DSQRT(TERM2**2+2.D0*(POT3DH12(RR,rsmall,T))**2)
      ROOT2=TERM1+DSQRT(TERM2**2+2.D0*(POT3DH12(RR,rsmall,T))**2)
      POT3DV1A=MIN(ROOT1,ROOT2)
                       ELSE
      ROOT1=TERM1-DSQRT(TERM2**2+2.D0*(POT3DH12(2.0D0,rsmall,T))**2)
      ROOT2=TERM1+DSQRT(TERM2**2+2.D0*(POT3DH12(2.0D0,rsmall,T))**2)
      POT3DV1A=MIN(ROOT1,ROOT2)
      ENDIF
      RETURN
      END

      FUNCTION POT3DV2A(RR,rsmall,T)
c    V2A=2Aprim
c     RR in A, T in degrees
c     output in cm-1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      TERM1=(POT3DH11(RR,rsmall,T)+POT3DH22(RR,rsmall,T))*0.5D0
      TERM2=(POT3DH11(RR,rsmall,T)-POT3DH22(RR,rsmall,T))*0.5D0
      IF (RR.GE.2.0D0) THEN
      ROOT1=TERM1-DSQRT(TERM2**2+2.D0*(POT3DH12(RR,rsmall,T))**2)
      ROOT2=TERM1+DSQRT(TERM2**2+2.D0*(POT3DH12(RR,rsmall,T))**2)
      POT3DV2A=MAX(ROOT1,ROOT2)
                       ELSE
      ROOT1=TERM1-DSQRT(TERM2**2+2.D0*(POT3DH12(2.0D0,rsmall,T))**2)
      ROOT2=TERM1+DSQRT(TERM2**2+2.D0*(POT3DH12(2.0D0,rsmall,T))**2)
      POT3DV2A=MAX(ROOT1,ROOT2)
      ENDIF
      RETURN
      END
C*******************************************************

      FUNCTION POT2DA1(R,rh2)
C FUNCTION FOR (R,r) DEPENDENCE OF
C F-H2 A1 STATE
C R,rh2 IN Angstrems
C OUTPUT IN MIKROHARTREE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension V0(9)
       CALL VlambdaA1(R,V0)
       re=0.74085D0
       z=(rh2-re)/re
       term1=V0(1)*z+VA1I(R,7)
       term2=V0(2)*z**2
       term3=V0(3)*z**3
       term4=V0(4)*z**4
       term5=V0(5)*z**5
       term6=V0(6)*z**6
       term7=V0(7)*z**7
       term8=V0(8)*z**8
       term9=V0(9)*z**9
       POT2DA1=term1+term2+term3+term4+
     * term5+term6+term7+term8+term9
       RETURN
      END
C
C
C
C
      subroutine VlambdaA1(R,V0)
      implicit double precision(a-h, o-z)
      dimension T(9,9),ipvt(9)
      dimension V0(9)
      pi=dacos(-1.d0)
      re=0.74085D0
      r1=0.42334D0
      r2=0.47626D0
      r3=0.52918D0
      r4=0.58210D0
      r5=0.63501D0
      r6=0.68793D0
      r7=re
      r8=0.79377D0
      r9=0.84668D0
      r10=0.89960D0
      t1=(r1-re)/re
      t2=(r2-re)/re
      t3=(r3-re)/re
      t4=(r4-re)/re
      t5=(r5-re)/re
      t6=(r6-re)/re
      t7=(r8-re)/re
      t8=(r9-re)/re
      t9=(r10-re)/re

      T(1,1)=t1
      T(1,2)=t1**2
      T(1,3)=t1**3
      T(1,4)=t1**4
      T(1,5)=t1**5
      T(1,6)=t1**6
      T(1,7)=t1**7
      T(1,8)=t1**8
      T(1,9)=t1**9
      T(2,1)=t2
      T(2,2)=t2**2
      T(2,3)=t2**3
      T(2,4)=t2**4
      T(2,5)=t2**5
      T(2,6)=t2**6
      T(2,7)=t2**7
      T(2,8)=t2**8
      T(2,9)=t2**9
      T(3,1)=t3
      T(3,2)=t3**2
      T(3,3)=t3**3
      T(3,4)=t3**4
      T(3,5)=t3**5
      T(3,6)=t3**6
      T(3,7)=t3**7
      T(3,8)=t3**8
      T(3,9)=t3**9
      T(4,1)=t4
      T(4,2)=t4**2
      T(4,3)=t4**3
      T(4,4)=t4**4
      T(4,5)=t4**5
      T(4,6)=t4**6
      T(4,7)=t4**7
      T(4,8)=t4**8
      T(4,9)=t4**9
      T(5,1)=t5
      T(5,2)=t5**2
      T(5,3)=t5**3
      T(5,4)=t5**4
      T(5,5)=t5**5
      T(5,6)=t5**6
      T(5,7)=t5**7
      T(5,8)=t5**8
      T(5,9)=t5**9
      T(6,1)=t6
      T(6,2)=t6**2
      T(6,3)=t6**3
      T(6,4)=t6**4
      T(6,5)=t6**5
      T(6,6)=t6**6
      T(6,7)=t6**7
      T(6,8)=t6**8
      T(6,9)=t6**9
      T(7,1)=t7
      T(7,2)=t7**2
      T(7,3)=t7**3
      T(7,4)=t7**4
      T(7,5)=t7**5
      T(7,6)=t7**6
      T(7,7)=t7**7
      T(7,8)=t7**8
      T(7,9)=t7**9
      T(8,1)=t8
      T(8,2)=t8**2
      T(8,3)=t8**3
      T(8,4)=t8**4
      T(8,5)=t8**5
      T(8,6)=t8**6
      T(8,7)=t8**7
      T(8,8)=t8**8
      T(8,9)=t8**9
      T(9,1)=t9
      T(9,2)=t9**2
      T(9,3)=t9**3
      T(9,4)=t9**4
      T(9,5)=t9**5
      T(9,6)=t9**6
      T(9,7)=t9**7
      T(9,8)=t9**8
      T(9,9)=t9**9
      V0(1)=VA1I(R,1)-VA1I(R,7)
      V0(2)=VA1I(R,2)-VA1I(R,7)
      V0(3)=VA1I(R,3)-VA1I(R,7)
      V0(4)=VA1I(R,4)-VA1I(R,7)
      V0(5)=VA1I(R,5)-VA1I(R,7)
      V0(6)=VA1I(R,6)-VA1I(R,7)
      V0(7)=VA1I(R,8)-VA1I(R,7)
      V0(8)=VA1I(R,9)-VA1I(R,7)
      V0(9)=VA1I(R,10)-VA1I(R,7)
      call dgesv(9,1,T,9,ipvt,V0,9,info)
C
      return
      end
C

      FUNCTION VA1I(R,N)
CSYSTEM: F-H2 A1
CLEVEL:RHF/UCCSD(T)/aug-cc-pvqz+332
C units R=A E=microEh
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(14),XX1(14),XX2(14),XX3(14),
     * XX4(14),XX5(14),XX6(14),XX7(14),XX8(14),
     * XX9(14),XX10(14)
C r=r1=0.42334 A
      DATA XX1/
     *  0.301418858448982929D+01,
     * -0.603705273584179025D+01,
     *  0.682915960700290307D+05,
     * -0.187666032260949782D+06,
     *  0.238301308334246220D+06,
     * -0.164244683880907105D+06,
     *  0.654163507393592809D+05,
     * -0.152247928229055105D+05,
     *  0.193569781326728071D+04,
     * -0.105288765445411940D+03,
     * -0.267482472856100245D-01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r2=0.47626 A
      DATA XX2/
     *  0.302063810440341030D+01,
     * -0.638833011929894568D+01,
     *  0.628740498834297978D+05,
     * -0.176096661305668618D+06,
     *  0.219680223613722919D+06,
     * -0.148219855928161065D+06,
     *  0.580248577947132362D+05,
     * -0.133253428905331384D+05,
     *  0.167675037673261022D+04,
     * -0.904316540418584083D+02,
     * -0.218200084608609485D-01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r3=0.52918 A
      DATA XX3/
     *  0.301685686644561724D+01,
     * -0.551775847247352313D+01,
     *  0.185856955473403039D+06,
     * -0.523721464309566189D+06,
     *  0.642625085935761104D+06,
     * -0.425782884896181815D+06,
     *  0.164142937796302431D+06,
     * -0.372112660932159197D+05,
     *  0.463005236787568083D+04,
     * -0.247007765851240379D+03,
     * -0.825900761957491514D-01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r4=0.58210 A
      DATA XX4/
     *  0.302848416668694886D+01,
     * -0.528421029366862793D+01,
     *  0.299444914273685019D+06,
     * -0.846202614471198991D+06,
     *  0.102471454895413003D+07,
     * -0.669979653675300884D+06,
     *  0.255722534873673751D+06,
     * -0.575729189852614654D+05,
     *  0.713228182899861167D+04,
     * -0.379708891873872460D+03,
     * -0.810470636018972429D-01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r5=0.63501 A
      DATA XX5/
     *  0.301903904174249416D+01,
     * -0.580166995255255191D+01,
     *  0.208814868503257632D+06,
     * -0.588130130126829026D+06,
     *  0.702723792203680612D+06,
     * -0.453320331730610400D+06,
     *  0.171024947143409197D+06,
     * -0.381150275527972772D+05,
     *  0.467918343905950132D+04,
     * -0.247067072552540992D+03,
     * -0.403090685500953963D-01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r6=0.68793 A
      DATA XX6/
     *  0.302492933202342806D+01,
     * -0.693082706820432204D+01,
     *  0.808806175727103546D+05,
     * -0.227090708047667606D+06,
     *  0.268582342694093648D+06,
     * -0.171656644661493512D+06,
     *  0.643261826310891774D+05,
     * -0.142737683146155960D+05,
     *  0.174882220895869500D+04,
     * -0.924922544120467336D+02,
     *  0.170244251061168403D-01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=re=0.74085 A
      DATA XX7/
     *  0.301980277505242789D+01,
     * -0.617615636080535335D+01,
     *  0.195263148647961119D+06,
     * -0.545196221413553925D+06,
     *  0.637964371158752940D+06,
     * -0.403625246676922659D+06,
     *  0.149959957821118820D+06,
     * -0.330330960155646899D+05,
     *  0.402097847767468920D+04,
     * -0.211193237785435201D+03,
     *  0.124209082604196828D-02,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r7=0.79377 A
      DATA XX8/
     *  0.301021338673093108D+01,
     * -0.624001964617968508D+01,
     *  0.202324783012565807D+06,
     * -0.561443567235070281D+06,
     *  0.650422295567868161D+06,
     * -0.407610947389787179D+06,
     *  0.150199530809161515D+06,
     * -0.328513252261865491D+05,
     *  0.397414152611111467D+04,
     * -0.207595239703559344D+03,
     *  0.282830888754052878D-02,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r8=0.84668 A
      DATA XX9/
     *  0.299847282996920583D+01,   
     * -0.655756498688726186D+01,
     *  0.159401107390797115D+06,
     * -0.439483762258756673D+06,
     *  0.504276070480560011D+06,
     * -0.313130611333813518D+06,
     *  0.114440752744809768D+06,
     * -0.248454294757396201D+05,
     *  0.298491272761102528D+04,
     * -0.154777916253962246D+03,
     * -0.203724640779095814D-01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r9=0.89960 A
      DATA XX10/
     *  0.380511452384088811D+01,
     * -0.908785189377794111D+01,
     * -0.122544933518701597D+06,
     *  0.458995191575567995D+06,
     * -0.732561192073940299D+06,
     *  0.652036493840829120D+06,
     * -0.355327266345501412D+06,
     *  0.121969156890437516D+06,
     * -0.258844174728603793D+05,
     *  0.311940518867429682D+04,
     * -0.164759837360351469D+03,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/

      IF(N.EQ.1) THEN
      DO I=1,14
      XX(I)=XX1(I)
      ENDDO
      ENDIF
      IF(N.EQ.2) THEN
      DO I=1,14
      XX(I)=XX2(I)
      ENDDO
      ENDIF
      IF(N.EQ.3) THEN
      DO I=1,14
      XX(I)=XX3(I)
      ENDDO
      ENDIF
      IF(N.EQ.4) THEN
      DO I=1,14
      XX(I)=XX4(I)
      ENDDO
      ENDIF
      IF(N.EQ.5) THEN
      DO I=1,14
      XX(I)=XX5(I)
      ENDDO
      ENDIF
      IF(N.EQ.6) THEN
      DO I=1,14
      XX(I)=XX6(I)
      ENDDO
      ENDIF
      IF(N.EQ.7) THEN
      DO I=1,14
      XX(I)=XX7(I)
      ENDDO
      ENDIF
      IF(N.EQ.8) THEN
      DO I=1,14
      XX(I)=XX8(I)
      ENDDO
      ENDIF
      IF(N.EQ.9) THEN
      DO I=1,14
      XX(I)=XX9(I)
      ENDDO
      ENDIF
      IF(N.EQ.10) THEN
      DO I=1,14
      XX(I)=XX10(I)
      ENDDO
      ENDIF


      TERM1 = DEXP(-XX(1)*R-XX(2))
      TERM2 = XX(3) + XX(4)*R+XX(5)*R**2+
     *        XX(6)*R**3+XX(7)*R**4+XX(8)*R**5
     *        +XX(9)*R**6+XX(10)*R**7+XX(11)*R**8
      TERM3=0.5D0*
     *   (1.D0+ DTANH(1.D0+XX(12)*R))
      TERM4 = XX(13)/(R**5)+XX(14)/(R**6)

       VA1I = (TERM1*TERM2-TERM3*TERM4)
       RETURN
       END
C
      FUNCTION POT2DB1(R,rh2)
C FUNCTION FOR (R,r) DEPENDENCE OF
C F-H2 B1 STATE
C R,rh2 IN Angstrems
C OUTPUT IN MIKROHARTREE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension V0(6)
       CALL VlambdaB1(R,V0)
       re=0.74085D0
       z=(rh2-re)/re
       term1=V0(1)*z+VB1I(R,4)
       term2=V0(2)*z**2
       term3=V0(3)*z**3
       term4=V0(4)*z**4
       term5=V0(5)*z**5
       term6=V0(6)*z**6
       POT2DB1=term1+term2+term3+term4
     * +term5+term6
       RETURN
      END
C
C
C
C
      subroutine VlambdaB1(R,V0)
      implicit double precision(a-h, o-z)
      dimension T(6,6),ipvt(6)
      dimension V0(6)
      pi=dacos(-1.d0)
      re=0.74085D0
      r1=0.42334D0
      r2=0.52918D0
      r3=0.63501D0
      r5=0.84668D0
      r6=0.95252D0
      r7=1.05835D0
      t1=(r1-re)/re
      t2=(r2-re)/re
      t3=(r3-re)/re
      t5=(r5-re)/re
      t6=(r6-re)/re
      t7=(r7-re)/re
      T(1,1)=t1
      T(1,2)=t1**2
      T(1,3)=t1**3
      T(1,4)=t1**4
      T(1,5)=t1**5
      T(1,6)=t1**6
      T(2,1)=t2
      T(2,2)=t2**2
      T(2,3)=t2**3
      T(2,4)=t2**4
      T(2,5)=t2**5
      T(2,6)=t2**6
      T(3,1)=t3
      T(3,2)=t3**2
      T(3,3)=t3**3
      T(3,4)=t3**4
      T(3,5)=t3**5
      T(3,6)=t3**6
      T(4,1)=t5
      T(4,2)=t5**2
      T(4,3)=t5**3
      T(4,4)=t5**4
      T(4,5)=t5**5
      T(4,6)=t5**6
      T(5,1)=t6
      T(5,2)=t6**2
      T(5,3)=t6**3
      T(5,4)=t6**4
      T(5,5)=t6**5
      T(5,6)=t6**6
      T(6,1)=t7
      T(6,2)=t7**2
      T(6,3)=t7**3
      T(6,4)=t7**4
      T(6,5)=t7**5
      T(6,6)=t7**6
      V0(1)=VB1I(R,1)-VB1I(R,4)
      V0(2)=VB1I(R,2)-VB1I(R,4)
      V0(3)=VB1I(R,3)-VB1I(R,4)
      V0(4)=VB1I(R,5)-VB1I(R,4)
      V0(5)=VB1I(R,6)-VB1I(R,4)
      V0(6)=VB1I(R,7)-VB1I(R,4)
      call dgesv(6,1,T,6,ipvt,V0,6,info)
C
      return
      end
C

      FUNCTION VB1I(R,N)
CSYSTEM: F-H2 B1 
CLEVEL:RHF/UCCSD(T)/aug-cc-pvqz+332
C units R=A E=microEh
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(14),XX1(14),XX2(14),XX3(14),
     * XX4(14),XX5(14),XX6(14),XX7(14)
C r=r1=0.42334 A
      DATA XX1/
     *  0.283357195276282159D+01,
     * -0.447996395278260007D+01,
     *  0.120474830125722481D+06,
     *  0.115015940999464419D+05,
     * -0.718971929144628084D+05,
     *  0.322034585517986015D+05,
     * -0.491576304350523515D+04,
     * -0.215086627651632853D+03,
     *  0.112713490885495432D+03,
     *  0.146634989562883871D-02,
     * -0.150112789126050239D+01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r2=0.52918 A
      DATA XX2/
     *  0.287565694372276637D+01,
     * -0.366724084316262822D+01,
     *  0.168661718634924589D+06,
     *  0.214034530736243818D+06,
     * -0.274929454983082716D+06,
     *  0.105687925716860424D+06,
     * -0.159634578179779619D+05,
     * -0.370851298572992278D+03,
     *  0.305974904741581042D+03,
     * -0.516782167836614750D-03,
     * -0.443307822670590834D+01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r3=0.63501 A
      DATA XX3/
     *  0.272357501448750794D+01,
     * -0.595235034416498454D+01,
     * -0.318407904728104040D+04,
     *  0.707250451042203640D+05,
     * -0.912179444956080115D+05,
     *  0.576819887385092807D+05,
     * -0.227460512564833298D+05,
     *  0.586382841600056418D+04,
     * -0.964449377738941280D+03,
     *  0.919700063136916839D+02,
     * -0.391895466702581707D+01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r4=0.74085 A
      DATA XX4/
     *  0.303641286457554482D+01,
     * -0.339987361360767126D+01,
     *  0.134927154486689280D+06,
     *  0.176307448709671007D+06,
     *  0.134958616753002862D+06,
     * -0.312351622227754153D+06,
     *  0.183894987760278251D+06,
     * -0.545020965378945184D+05,
     *  0.840349416828280118D+04,
     * -0.551336151607739453D+03,
     * -0.340617927152318768D-03,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r5=0.84668 A
      DATA XX5/
     *  0.314670355154628645D+01,
     * -0.421978869926994360D+01,
     *  0.380897949482160030D+02,
     *  0.203825493417480378D+06,
     * -0.799566701096517354D+05,
     *  0.396907951227411848D+00,
     * -0.691896869880620670D+04,
     *  0.997650545223350673D+04,
     * -0.425004450481778258D+04,
     *  0.803118691836126004D+03,
     * -0.606881967303535959D+02,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r6=0.95252 A
      DATA XX6/
     *  0.276653030787949072D+01,
     * -0.601197717883762461D+01,
     * -0.211874523650474330D+05,
     *  0.914183935555512144D+05,
     * -0.987666230306430953D+05,
     *  0.598198122497683580D+05,
     * -0.240408687476286286D+05,
     *  0.646133974122356358D+04,
     * -0.111220215202421628D+04,
     *  0.110682612636269241D+03,
     * -0.491308095624267693D+01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/
C r=r7=1.05835 A
      DATA XX7/
     *  0.317856973781594565D+01,
     * -0.475694080460601221D+01,
     * -0.741846712837660471D+04,
     *  0.824400855359893321D+05,
     * -0.910671521711559184D+00,
     * -0.594514654130900908D+04,
     * -0.163809303451785072D+05,
     *  0.139203513135412231D+05,
     * -0.485591071112742793D+04,
     *  0.823138186948582870D+03,
     * -0.576795108934895779D+02,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/

        IF(N.EQ.1) THEN
        DO I=1,14
        XX(I)=XX1(I)
        ENDDO
        ENDIF
        IF(N.EQ.2) THEN
        DO I=1,14
        XX(I)=XX2(I)
        ENDDO
        ENDIF
        IF(N.EQ.3) THEN
        DO I=1,14
        XX(I)=XX3(I)
        ENDDO
        ENDIF
        IF(N.EQ.4) THEN
        DO I=1,14
        XX(I)=XX4(I)
        ENDDO
        ENDIF
        IF(N.EQ.5) THEN
        DO I=1,14
        XX(I)=XX5(I)
        ENDDO
        ENDIF
        IF(N.EQ.6) THEN
        DO I=1,14
        XX(I)=XX6(I)
        ENDDO
        ENDIF
        IF(N.EQ.7) THEN
        DO I=1,14
        XX(I)=XX7(I)
        ENDDO
        ENDIF

      TERM1 = DEXP(-XX(1)*R-XX(2))
      TERM2 = XX(3) + XX(4)*R+XX(5)*R**2+
     *        XX(6)*R**3+XX(7)*R**4+XX(8)*R**5
     *        +XX(9)*R**6+XX(10)*R**7+XX(11)*R**8
      TERM3=0.5D0*
     *   (1.D0+ DTANH(1.D0+XX(12)*R))
      TERM4 = XX(13)/(R**5)+XX(14)/(R**6)

       VB1I = (TERM1*TERM2-TERM3*TERM4)
       RETURN
       END
C
      FUNCTION POT2DB2(R,rh2)
C FUNCTION FOR (R,r) DEPENDENCE OF
C F-H2 B2 STATE
C R,rh2 IN Angstrems
C OUTPUT IN MIKROHARTREE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension V0(6)
       CALL VlambdaB2(R,V0)
       re=0.74085D0
       z=(rh2-re)/re
       term1=V0(1)*z+VB2I(R,4)
       term2=V0(2)*z**2
       term3=V0(3)*z**3
       term4=V0(4)*z**4
       term5=V0(5)*z**5
       term6=V0(6)*z**6
       POT2DB2=term1+term2+term3+term4
     * +term5+term6
       RETURN
      END
C
C
C
C
      subroutine VlambdaB2(R,V0)
      implicit double precision(a-h, o-z)
      dimension T(6,6),ipvt(6)
      dimension V0(6)
      pi=dacos(-1.d0)
      re=0.74085D0
      r1=0.42334D0
      r2=0.52918D0
      r3=0.63501D0
      r5=0.84668D0
      r6=0.95252D0
      r7=1.05835D0
      t1=(r1-re)/re
      t2=(r2-re)/re
      t3=(r3-re)/re
      t5=(r5-re)/re
      t6=(r6-re)/re
      t7=(r7-re)/re
      T(1,1)=t1
      T(1,2)=t1**2
      T(1,3)=t1**3
      T(1,4)=t1**4
      T(1,5)=t1**5
      T(1,6)=t1**6
      T(2,1)=t2
      T(2,2)=t2**2
      T(2,3)=t2**3
      T(2,4)=t2**4
      T(2,5)=t2**5
      T(2,6)=t2**6
      T(3,1)=t3
      T(3,2)=t3**2
      T(3,3)=t3**3
      T(3,4)=t3**4
      T(3,5)=t3**5
      T(3,6)=t3**6
      T(4,1)=t5
      T(4,2)=t5**2
      T(4,3)=t5**3
      T(4,4)=t5**4
      T(4,5)=t5**5
      T(4,6)=t5**6
      T(5,1)=t6
      T(5,2)=t6**2
      T(5,3)=t6**3
      T(5,4)=t6**4
      T(5,5)=t6**5
      T(5,6)=t6**6
      T(6,1)=t7
      T(6,2)=t7**2
      T(6,3)=t7**3
      T(6,4)=t7**4
      T(6,5)=t7**5
      T(6,6)=t7**6
      V0(1)=VB2I(R,1)-VB2I(R,4)
      V0(2)=VB2I(R,2)-VB2I(R,4)
      V0(3)=VB2I(R,3)-VB2I(R,4)
      V0(4)=VB2I(R,5)-VB2I(R,4)
      V0(5)=VB2I(R,6)-VB2I(R,4)
      V0(6)=VB2I(R,7)-VB2I(R,4)
      call dgesv(6,1,T,6,ipvt,V0,6,info)
C
      return
      end
C

      FUNCTION VB2I(R,N)
CSYSTEM: F-H2 B2
CLEVEL:RHF/UCCSD(T)/aug-cc-pvqz+332
C units R=A E=microEh
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(14),XX1(14),XX2(14),XX3(14),
     * XX4(14),XX5(14),XX6(14),XX7(14)
C r=r1=0.42334 A
      DATA XX1/
     *  0.295390872587323061D+01,
     * -0.453947093032389226D+01,
     *  0.118156228584792130D+06,
     *  0.224167464154384543D+05,
     * -0.690554894091893366D+05,
     *  0.266362877211812229D+05,
     * -0.355625183723824512D+04,
     * -0.612217862939555546D-02,
     * -0.479362034306542313D+02,
     *  0.293949486242627032D+02,
     * -0.351786738641760000D+01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/

C r=r2=0.52918 A
      DATA XX2/
     *  0.275155005023126487D+01,
     * -0.615709853193104717D+01,
     *  0.945635569747667978D+04,
     *  0.341928142287404407D+05,
     * -0.520340841086358196D+05,
     *  0.333593620514415379D+05,
     * -0.128791165468607669D+05,
     *  0.321851674394995507D+04,
     * -0.512455089089091530D+03,
     *  0.473214849876440127D+02,
     * -0.195438603790498422D+01,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/

C r=r3=0.63501 A
      DATA XX3/
     *  0.333725803136194044D+01,
     * -0.267880823655760159D+01,
     *  0.596950710992069566D+06,
     *  0.163943950042255544D+02,
     *  0.867914629123724997D+06,
     * -0.981230853226367501D+06,
     *  0.450003677472476382D+06,
     * -0.103125263823762230D+06,
     *  0.814011750871661661D+04,
     *  0.934650894937899579D+03,
     * -0.174645165794070209D+03,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/

C r=r4=0.74085 A
      DATA XX4/
     *  0.287466199339493267D+01,
     * -0.740940376009246826D+01,
     * -0.219249015717386010D+04,
     *  0.197722690470349771D+05,
     * -0.226850242782013411D+05,
     *  0.136138783780608283D+05,
     * -0.526133850926898958D+04,
     *  0.133656314205604372D+04,
     * -0.214818339537281190D+03,
     *  0.197390619170379082D+02,
     * -0.801984425698802683D+00,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/

C r=r5=0.84668 A
      DATA XX5/
     *  0.344635486563916960D+01,
     * -0.226638385192963199D+01,
     *  0.124263179003213765D+06,
     *  0.107743554458265635D+07,
     *  0.867351491690326133D+06,
     * -0.133065270364417532D+07,
     *  0.697818993751274538D+06,
     * -0.195084434963865147D+06,
     *  0.243627862149680877D+05,
     *  0.111579050816088432D+01,
     * -0.233502575631271128D+03,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/

C r=r6=0.95252 A
      DATA XX6/
     *  0.334637810745046593D+01,
     * -0.687657540987928062D+01,
     * -0.817659920031785987D+04,
     *  0.352538520695417465D+05,
     * -0.240665573384614945D+05,
     *  0.959739155079603006D+04,
     * -0.220114290909493639D+04,
     *  0.866968678060150069D+02,
     *  0.575741998160614372D+02,
     * -0.889339266248669880D+01,
     * -0.655469724938011173D-04,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/

C r=r7=1.05835 A
      DATA XX7/
     *  0.349595516987215493D+01,
     * -0.371815832436441784D+01,
     *  0.301121766462376783D+02,
     * -0.115853613700391384D+04,
     *  0.778787608245793148D+06,
     * -0.872544104222255759D+06,
     *  0.506782723350841203D+06,
     * -0.176643484336668451D+06,
     *  0.355901030709365878D+05,
     * -0.375907062809635636D+04,
     *  0.143508214783077364D+03,
     * -0.289573768192804302D+02,
     *  0.876525683507813606D+06,
     *  0.387887406511697778D+07/

        IF(N.EQ.1) THEN
        DO I=1,14
        XX(I)=XX1(I)
        ENDDO
        ENDIF
        IF(N.EQ.2) THEN
        DO I=1,14
        XX(I)=XX2(I)
        ENDDO
        ENDIF
        IF(N.EQ.3) THEN
        DO I=1,14
        XX(I)=XX3(I)
        ENDDO
        ENDIF
        IF(N.EQ.4) THEN
        DO I=1,14
        XX(I)=XX4(I)
        ENDDO
        ENDIF
        IF(N.EQ.5) THEN
        DO I=1,14
        XX(I)=XX5(I)
        ENDDO
        ENDIF
        IF(N.EQ.6) THEN
        DO I=1,14
        XX(I)=XX6(I)
        ENDDO
        ENDIF
        IF(N.EQ.7) THEN
        DO I=1,14
        XX(I)=XX7(I)
        ENDDO
        ENDIF

      TERM1 = DEXP(-XX(1)*R-XX(2))
      TERM2 = XX(3) + XX(4)*R+XX(5)*R**2+
     *        XX(6)*R**3+XX(7)*R**4+XX(8)*R**5
     *        +XX(9)*R**6+XX(10)*R**7+XX(11)*R**8
      TERM3=0.5D0*
     *   (1.D0+ DTANH(1.D0+XX(12)*R))
      TERM4 = XX(13)/(R**5)+XX(14)/(R**6)

       VB2I = (TERM1*TERM2-TERM3*TERM4)
       RETURN
       END
C
      FUNCTION POT2DPI(R,rh2)
C FUNCTION FOR (R,r) DEPENDENCE OF
C F-H2 PI STATE
C R,rh2 IN Angstrems
C OUTPUT IN MIKROHARTREE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension V0(5)
       CALL VlambdaPI(R,V0)
       re=0.74085D0
       z=(rh2-re)/re
       term1=V0(1)*z+VPII(R,3)
       term2=V0(2)*z**2
       term3=V0(3)*z**3
       term4=V0(4)*z**4
       term5=V0(5)*z**5
       POT2DPI=term1+term2+term3
     *         +term4+term5 
       RETURN
      END
C
C
C
C
      subroutine VlambdaPI(R,V0)
      implicit double precision(a-h, o-z)
      dimension T(5,5),ipvt(5)
      dimension V0(5)
      pi=dacos(-1.d0)
      re=0.74085D0
      r1=0.52918D0
      r2=0.63501D0
      r3=0.84668D0
      r4=0.95252D0
      r5=1.05835D0
      t1=(r1-re)/re
      t2=(r2-re)/re
      t3=(r3-re)/re
      t4=(r4-re)/re
      t5=(r5-re)/re
      T(1,1)=t1
      T(1,2)=t1**2
      T(1,3)=t1**3
      T(1,4)=t1**4
      T(1,5)=t1**5
      T(2,1)=t2
      T(2,2)=t2**2
      T(2,3)=t2**3
      T(2,4)=t2**4
      T(2,5)=t2**5
      T(3,1)=t3
      T(3,2)=t3**2
      T(3,3)=t3**3
      T(3,4)=t3**4
      T(3,5)=t3**5
      T(4,1)=t4
      T(4,2)=t4**2
      T(4,3)=t4**3
      T(4,4)=t4**4
      T(4,5)=t4**5
      T(5,1)=t5
      T(5,2)=t5**2
      T(5,3)=t5**3
      T(5,4)=t5**4
      T(5,5)=t5**5



      V0(1)=VPII(R,1)-VPII(R,3)
      V0(2)=VPII(R,2)-VPII(R,3)
      V0(3)=VPII(R,4)-VPII(R,3)
      V0(4)=VPII(R,5)-VPII(R,3)
      V0(5)=VPII(R,6)-VPII(R,3)
      call dgesv(5,1,T,5,ipvt,V0,5,info)
C
      return
      end
C

      FUNCTION VPII(R,N)
CSYSTEM: F-H2 PI
CLEVEL:RHF/UCCSD(T)/aug-cc-pvqz+332
C units R=A E=microEh
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(14),XX1(14),XX2(14),XX3(14),
     * XX4(14),XX5(14),XX6(14)
C r=r1=0.52918 A
      DATA XX1/
     *   0.288101387672307796D+01,
     *  -0.354829457167905948D+01,
     *   0.291136354246576608D+06,
     *   0.369464696956159969D+05,
     *   0.237963337296140143D+01,
     *  -0.141753258082832588D+06,
     *   0.104161230265923543D+06,
     *  -0.333158711217429227D+05,
     *   0.519580324533040766D+04,
     *  -0.321386881706345662D+03,
     *  -0.241476449572259400D+01,
     *  -0.289573768192804302D+02,
     *   0.876525683507813606D+06,
     *   0.387887406511697778D+07/

C r=r2=0.63501 A
      DATA XX2/
     *   0.242584637231939126D+01,
     *  -0.445927465567822878D+01,
     *   0.441012349872669001D+05,
     *   0.167905792207280698D+06,
     *  -0.261928162331646425D+06,
     *   0.163923933653356362D+06,
     *  -0.593695738761498578D+05,
     *   0.136318252623030912D+05,
     *  -0.198993785386692002D+04,
     *   0.170439090478737825D+03,
     *  -0.665517431122651004D+01,
     *  -0.289573768192804302D+02,
     *   0.876525683507813606D+06,
     *   0.387887406511697778D+07/

C r=r3=0.74085 A
      DATA XX3/
     *   0.270604825982857111D+01,
     *  -0.379856573477150539D+01,
     *   0.219345198441248707D+06,
     *  -0.942320080548216938D+05,
     *   0.159497867944249418D+06,
     *  -0.178063998875649500D+06,
     *   0.837444404448538407D+05,
     *  -0.192982070065093685D+05,
     *   0.192860335836597301D+04,
     *   0.111114209799324945D+00,
     *  -0.107682709499718499D+02,
     *  -0.289573768192804302D+02,
     *   0.876525683507813606D+06,
     *   0.387887406511697778D+07/

C r=r4=0.84668 A
      DATA XX4/
     *   0.295616710734156696D+01,
     *  -0.582731330195263908D+01,
     *   0.119370392994385213D+06,
     *  -0.257451496362622100D+06,
     *   0.304612783285240177D+06,
     *  -0.186315242876922130D+06,
     *   0.599656510847001628D+05,
     *  -0.902007764541772121D+04,
     *  -0.399783258697925872D-01,
     *   0.175809304071151104D+03,
     *  -0.166037450690838355D+02,
     *  -0.289573768192804302D+02,
     *   0.876525683507813606D+06,
     *   0.387887406511697778D+07/

C r=r5=0.95252 A
      DATA XX5/
     *   0.249851276811401180D+01,
     *  -0.365805308591928435D+01,
     *   0.168273506377212878D+06,
     *  -0.525591948786162684D+05,
     *   0.144421037511982751D+06,
     *  -0.165331041044362937D+06,
     *   0.738001964887552749D+05,
     *  -0.157905495102099339D+05,
     *   0.144879681491836209D+04,
     *   0.589935257060922175D-01,
     *  -0.661101416330000635D+01,
     *  -0.289573768192804302D+02,
     *   0.876525683507813606D+06,
     *   0.387887406511697778D+07/

C r=r6=1.05835 A
      DATA XX6/
     *   0.255923875968635617D+01,
     *  -0.533068068602461409D+01,
     *   0.107771077313168193D+06,
     *  -0.219044839917762205D+06,
     *   0.265874971309535264D+06,
     *  -0.173347095176382485D+06,
     *   0.630010905801819827D+05,
     *  -0.131205994379011136D+05,
     *   0.148231350114034808D+04,
     *  -0.722857133485914716D+02,
     *   0.792093501325373955D-01,
     *  -0.289573768192804302D+02,
     *   0.876525683507813606D+06,
     *   0.387887406511697778D+07/



        IF(N.EQ.1) THEN
        DO I=1,14
        XX(I)=XX1(I)
        ENDDO
        ENDIF
        IF(N.EQ.2) THEN
        DO I=1,14
        XX(I)=XX2(I)
        ENDDO
        ENDIF
        IF(N.EQ.3) THEN
        DO I=1,14
        XX(I)=XX3(I)
        ENDDO
        ENDIF
        IF(N.EQ.4) THEN
        DO I=1,14
        XX(I)=XX4(I)
        ENDDO
        ENDIF
        IF(N.EQ.5) THEN
        DO I=1,14
        XX(I)=XX5(I)
        ENDDO
        ENDIF
        IF(N.EQ.6) THEN
        DO I=1,14
        XX(I)=XX6(I)
        ENDDO
        ENDIF



      TERM1 = DEXP(-XX(1)*R-XX(2))
      TERM2 = XX(3) + XX(4)*R+XX(5)*R**2+
     *        XX(6)*R**3+XX(7)*R**4+XX(8)*R**5
     *        +XX(9)*R**6+XX(10)*R**7+XX(11)*R**8
      TERM3=0.5D0*
     *   (1.D0+ DTANH(1.D0+XX(12)*R))
      TERM4 = XX(13)/(R**5)+XX(14)/(R**6)

       VPII = (TERM1*TERM2-TERM3*TERM4)
       RETURN
       END
C
      FUNCTION POT2DSG(R,rh2)
C FUNCTION FOR (R,r) DEPENDENCE OF
C F-H2 SIGMA STATE
C R,rh2 IN Angstrems
C OUTPUT IN MIKROHARTREE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension V0(3)
       CALL VlambdaSG(R,V0)
       re=0.74085D0
       z=(rh2-re)/re
       term1=V0(1)*z+VSGI(R,3)
       term2=V0(2)*z**2
       term3=V0(3)*z**3
       POT2DSG=term1+term2+term3
       RETURN
      END
C
C
C
C
      subroutine VlambdaSG(R,V0)
      implicit double precision(a-h, o-z)
      dimension T(3,3),ipvt(3)
      dimension V0(3)
      pi=dacos(-1.d0)
      re=0.74085D0
      r1=0.52918D0
      r2=0.63501D0
      r3=0.84668D0
      t1=(r1-re)/re
      t2=(r2-re)/re
      t3=(r3-re)/re
      T(1,1)=t1
      T(1,2)=t1**2
      T(1,3)=t1**3
      T(2,1)=t2
      T(2,2)=t2**2
      T(2,3)=t2**3
      T(3,1)=t3
      T(3,2)=t3**2
      T(3,3)=t3**3
      V0(1)=VSGI(R,1)-VSGI(R,3)
      V0(2)=VSGI(R,2)-VSGI(R,3)
      V0(3)=VSGI(R,4)-VSGI(R,3)
      call dgesv(3,1,T,3,ipvt,V0,3,info)
C
      return
      end
C

      FUNCTION VSGI(R,N)
CSYSTEM: F-H2 SIGMA
CLEVEL:RHF/UCCSD(T)/aug-cc-pvqz+332
C units R=A E=microEh
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(14),XX1(14),XX2(14),XX3(14),
     * XX4(14)
C r=r1=0.52918 A
      DATA XX1/
     *   0.347621098314397114D+01,
     *  -0.798710130301622012D+01,
     *  -0.467617238842874722D+05,
     *   0.491079091587328803D+05,
     *  -0.855869016233460206D+04,
     *  -0.101551258834595355D+05,
     *   0.477777240120719489D+04,
     *  -0.389575385098098886D-01,
     *  -0.487570912230031070D+03,
     *   0.125758103534125368D+03,
     *  -0.110594558700220453D+02,
     *   0.768592052646715018D+01,
     *   0.236725079764586641D+06,
     *  -0.173346530057967524D+07/

C r=r2=0.63501 A
      DATA XX2/
     *   0.333151909247473288D+01,
     *  -0.902874771605496029D+01,
     *  -0.418388904657136954D+05,
     *   0.538730352252237644D+05,
     *  -0.238232960317764409D+05,
     *  -0.612759984206418380D+03,
     *   0.404250979379386172D+04,
     *  -0.143924324941787609D+04,
     *   0.185531786841853489D+03,
     *  -0.235510589418366099D-05,
     *  -0.167298839352105655D+01,
     *   0.435400419555694447D+06,
     *   0.594406277535313508D+06,
     *  -0.432590569201788213D+07/

C r=r3=0.74085 A
      DATA XX3/
     *  0.246118595428434483D+01,
     * -0.106225896220179212D+02,
     * -0.144803130912030902D+05,
     *  0.297159808611018743D+05,
     * -0.280533690955581515D+05,
     *  0.154955056128223841D+05,
     * -0.545145928921563427D+04,
     *  0.124588772032276438D+04,
     * -0.180675232499912738D+03,
     *  0.152161667237378868D+02,
     * -0.575821782221406453D+00,
     *  0.435400419555694447D+06,
     * -0.486618092994658102D+06,
     * -0.718598600303590298D+07/

C r=r4=0.84668 A
      DATA XX4/
     *  0.340670754329443470D+01,
     * -0.144041266721032706D+02,
     *  0.747437940712431828D+04,
     * -0.190768193189354497D+05,
     *  0.223218687851440191D+05,
     * -0.151890424974961279D+05,
     *  0.654937928075852597D+04,
     * -0.182777406838075399D+04,
     *  0.323214681866157889D+03,
     * -0.332235709696224362D+02,
     *  0.154443880986455384D+01,
     *  0.435400419555694447D+06,
     * -0.542820067369836010D+07,
     *  0.393913160980229229D+08/

        IF(N.EQ.1) THEN
        DO I=1,14
        XX(I)=XX1(I)
        ENDDO
        ENDIF
        IF(N.EQ.2) THEN
        DO I=1,14
        XX(I)=XX2(I)
        ENDDO
        ENDIF
        IF(N.EQ.3) THEN
        DO I=1,14
        XX(I)=XX3(I)
        ENDDO
        ENDIF
        IF(N.EQ.4) THEN
        DO I=1,14
        XX(I)=XX4(I)
        ENDDO
        ENDIF

      TERM1 = DEXP(-XX(1)*R-XX(2))
      TERM2 = XX(3) + XX(4)*R+XX(5)*R**2+
     *        XX(6)*R**3+XX(7)*R**4+XX(8)*R**5
     *        +XX(9)*R**6+XX(10)*R**7+XX(11)*R**8
      TERM3=0.5D0*
     *   (1.D0+ DTANH(1.D0+XX(12)*R))
      TERM4 = XX(13)/(R**5)+XX(14)/(R**6)

       VSGI = (TERM1*TERM2-TERM3*TERM4)
       RETURN
       END

c***H12 COUPLING 3-D MRCI****************
      FUNCTION POT3DH12(RR,rsmall,T)
C F-H2 system STATE:H12 MRCI 
C units R=A E=cm-1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(35)
        DATA XX/
     *    0.146785045039485578D+02,
     *    0.639375370693218148D+01,
     *    0.151309417676178299D+03,
     *    0.274355254615448914D+02,
     *    0.182771692231812133D+01,
     *    0.106138270858910547D+00,
     *   -0.257751121801726955D-01,
     *    0.178071479524053472D+00,
     *   -0.199256700458829187D+03,
     *   -0.348623967196420921D+02,
     *   -0.189030316773231233D+01,
     *    0.141927268991115241D+00,
     *    0.396013774912161512D-01,
     *   -0.201682930487279699D+00,
     *    0.880168246801674172D+02,
     *    0.146941474850966749D+02,
     *    0.625184715211922315D+00,
     *   -0.166144804878691321D+00,
     *   -0.867688162199622276D-02,
     *    0.763923456040812299D-01,
     *   -0.132195819894863451D+02,
     *   -0.206769189296769618D+01,
     *   -0.628160580822506864D-01,
     *    0.358512231466525388D-01,
     *   -0.105914455880216072D-02,
     *   -0.970636190221498517D-02,
     *    0.799138571452743122D+03,
     *    0.102809848102600995D+03,
     *    0.733001295294727928D+01,
     *    0.399438638172067906D+01,
     *    0.500726073851736242D+01,
     *    0.111641442903222288D+02,
     *    0.263797551918497923D+02,
     *   -0.147978375673764972D+01,
     *   -0.192970042731921265D+01/
      PI = DACOS(-1.D0)
      CT=DCOS(T*PI/180.d0)
      DTERM = 0.D0
      LLEGD=1
      LLEGB=1
      LLEGG=6
      DO 80 L=1,LLEGD
      PLNORMD = DSQRT((2.D0*DFLOAT(L-1)+1.D0)/2.D0)
       DTERM = DTERM + XX(L)*PLNORMD*PLGNDR(L-1,0,CT)
80    CONTINUE
      BTERM = 0.D0
      DO 90 L=1,LLEGB
      PLNORMB = DSQRT((2.D0*DFLOAT(L-1)+1.D0)/2.D0)
       BTERM = BTERM + XX(LLEGD+L)*PLNORMB*PLGNDR(L-1,0,CT)
90    CONTINUE
      GTERM = 0.D0
      DO 100 L=1,LLEGG
      PLNORMG = DSQRT((2.D0*DFLOAT(2*L)+1.D0)/2.D0)
       GTERM = GTERM +
     * (XX(LLEGD+LLEGB+L)+
     * RR*XX(LLEGD+LLEGB+LLEGG+L)+
     * RR*RR*XX(LLEGD+LLEGB+LLEGG+LLEGG+L)+
     * RR*RR*RR*XX(LLEGD+LLEGB+LLEGG+LLEGG+LLEGG+L))*
     * PLNORMG*PLGNDR(2*L,1,CT)
100   CONTINUE
      LLEGC = LLEGD+LLEGB+LLEGG+LLEGG+LLEGG+LLEGG
       RAD3 = RR**3
       RAD4 = RR*RAD3
       RAD5 = RR*RAD4
       RAD6 = RR*RAD5
       RAD7 = RR*RAD6
       RAD8 = RR*RAD7
       RAD9 = RR*RAD8
       RAD10 = RR*RAD9
       RAD11 = RR*RAD10
       RAD12 = RR*RAD11
       RAD13 = RR*RAD12
C VAN DER WAALS PART
       re=0.7408D0
       CTERM5 = XX(LLEGC+1)/(RAD5)*D(5,BTERM,RR)*
     * DSQRT((2.D0*2.D0+1.D0)/2.D0)*PLGNDR(2,1,CT)
       CTERM6 = XX(LLEGC+2)/(RAD6)*D(6,BTERM,RR)*
     * DSQRT((2.D0*4.D0+1.D0)/2.D0)*PLGNDR(4,1,CT)
       CTERM7 = XX(LLEGC+3)/(RAD7)*D(7,BTERM,RR)*
     * DSQRT((2.D0*6.D0+1.D0)/2.D0)*PLGNDR(6,1,CT)
C INTERNAL COORDINATE
       ZVAR = (rsmall-re)/re
       VTERM = XX(LLEGC+4)+XX(LLEGC+5)*ZVAR+
     *         XX(LLEGC+6)*ZVAR**2+XX(LLEGC+7)*ZVAR**3
       VTERMEXP=DEXP(-XX(LLEGC+8)*ZVAR-XX(LLEGC+9)*ZVAR**2)
         CTERM = CTERM5+CTERM6+CTERM7
C FUNCTION VALUE

       POT3DH12 =(GTERM*DEXP(DTERM -RR*BTERM)*VTERMEXP
     *          + CTERM)*VTERM

       RETURN
C LAST CARD IN POT3DH12 FUNCTION
       END

       FUNCTION D(nn,beta,r)
c
c     calculate the damping factor (small R correct)
c
      implicit real*8 (a-h,o-z)
      br=beta*r
      sum=1.0d0
      term=1.0d0
      ncn=nn
      do i=1,ncn
        term=term*br/i
        sum=sum+term
      enddo
      d=1.0d0 - dexp(-br)*sum
c     in case of d --> 0 use
c     d=1.0d0 - dexp(-br)*sum = sum_m=ncn+1^\infty br^m/m!
      if(dabs(d).lt.1.0d-8) then
        d=0.0d0
        do i=ncn+1,1000
          term=term*br/i
          d=d+term
          if(term/d .lt. 1.0d-8) go to 111
        enddo
        write(6,*) 'No convergence in d'
  111 continue
      d=d*dexp(-br)
      endif
      return
      end

      FUNCTION PLGNDR(L,M,x)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.)pause
     *'bad arguments in plgndr'
      pmm=1.
      if(m.gt.0) then
        somx2=dsqrt((1.-x)*(1.+x))
        fact=1.
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.
11      continue
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
         do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndr=pll
        endif
      endif
      return
      END
