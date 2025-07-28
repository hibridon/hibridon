# Restarting Calculations

If you have a long-running calculation which is interrupted by a power failure or other system crash, you may be able to restart by
- Setting the [flag](./Flags.md) [RSFLAG](./param/rsflag.md) = .TRUE.
- and then re-running the job with the original input data, unchanged except for [RSFLAG](./param/rsflag.md)

Similarly, if you finish a calculation, and wish to extend it to include addition partial waves, you can
- Set the [flag](./Flags.md) RSFLAG = .TRUE.
- Increase the value of [`JTOT2`](./param/jtot2.md)
- then re-run the job with the original input data, otherwise unchanged

Restarting will work only if either (or both) of the flags [`WRPART`](./param/wrpart.md) or [`WRXSEC`](./param/wrxsec.md) were .TRUE. in the input to the first job.

## Restarting a calculation: an example

In this example, I have added illustrative comments to the actual output of the run. In addition, the various commands
have been linked to the corresponding help files

The input file used here is [tests/Arn2_test.inp] in the Hibridon directory tree

The potential subroutine is `src/pot/pot_arn2.f`

Initiate execution of your code

```
% hib_arn2_151

 --------------------------------------------------------------------------
           HIBRIDON SCATTERING CODE V 4.1 05/17/97 14:23:43 EDT

     AUTHORS: M. ALEXANDER, D. MANOLOPOULOS, H.-J. WERNER, B. FOLLMEG
 CONTRIBUTORS: D. LEMOINE, P. VOHRALIK, G. COREY, R. JOHNSON, T. ORLIKOWSKI
               A. BERNING, A. DEGLI-ESPOSTI, C. RIST, P. DAGDIGIAN, B. POUILLY
               G. VAN DER SANDEN, M. YANG, F. DE WEERD, S. GREGURICK
 --------------------------------------------------------------------------
```

read in the input file, note that the name is in lower case.
the first character will be converted automatically to upper case)

```
Hibridon> inp=arn2_test.inp
Hibridon> show
     *** Parameters (scattering):
 JTOT1  =  20        JTOT2  =  20        JTOTD  =   5        JLPAR  =   1
 NERG   =   1        NUMAX  =   0        NUMIN  =   0        NUD    =   1
 LSCREEN=  48        IPRINT =   0
 FSTFAC =  15.00     RINCR  =  3.000     RCUT   =  30.00     RENDAI =  25.00
 RENDLD =  8.000     RSTART =  5.600     SPAC   =  .1500     TOLAI  =  1.150
 XMU    =  16.47
 NOUT:  4; JOUT:   0   2   4   6
 INDOUT:    0
     *** 1-SIGMA system parameters:
 NTERM  =   1        VMIN   =   0        VMAX   =   0        JMIN   =   0
 JMAX   =   4
 BROT   =  2.010     DROT   = 0.0000E+00 HROT   = 0.0000E+00 EVIB   = 0.0000E+00
 LAMMIN:    2
 LAMMAX:    2
 MPROJ:     0
     *** Flags:
 AIRYFL= T    BASTST= F    BATCH = F    CHLIST= T    CSFLAG= F    FLAGHF= F
 FLAGSU= F    IHOMO = T    IPOS  = F    LOGDFL= T    NOPRIN= F    NUCROS= F
 PHOTOF= F    PRAIRY= F    PRLOGD= T    PRPART= F    PRSMAT= T    PRT2  = T
 PRXSEC= F    READPT= F    RSFLAG= F    T2TEST= F    TWOMOL= F    WAVEFL= F
 WRPART= F    WRSMAT= F    WRXSEC= F    BOUNDC= F
 ** Maximum Channels:  151; Anisotropic Terms:  80
 ** Energies:      500.000000
 ** Label:           N2-Ar CC PATTENGILL POTENTIAL
 ** Pot name:      PATTENGILL-LABUDDE-BERNSTEIN AR-N2
 ** Input File:  Arn2_test.inp
 ** Output file: Outpt
 ** Jobname:     Job
 ```

Change the intial values of JTOT1, JTOT2, and JTOTD so that the calculation will cover 6 values of JTOT

Also supress printing by setting [NOPRIN](./param/noprin.md) = .TRUE.

```
Hibridon> jtot1=0;jtot2=10;jtotd=2;noprin=t
```

 set the flag [WRXSEC](./param/wrxsec.md) .TRUE., so that the calculation will accumulate the integral cross sections

```
 Hibridon> wrxsec=t
 Hibridon> save
 *** INPUT FILE Arn2_test.inp OVERWRITTEN
 Hibridon> run
 NOPRIN = .TRUE., SO IPRINT SET TO -1
 ** J =    0 JLPAR = 1 FINISHED;  CPU: 00:00:00.02   WALL: 00:00:00.14   DATE:
17-May-97  14:27:55
 ** J =    2 JLPAR = 1 FINISHED;  CPU: 00:00:00.03   WALL: 00:00:00.29   DATE:
17-May-97  14:27:55
 ** J =    4 JLPAR = 1 FINISHED;  CPU: 00:00:00.05   WALL: 00:00:00.46   DATE:
17-May-97  14:27:55
 ** J =    6 JLPAR = 1 FINISHED;  CPU: 00:00:00.08   WALL: 00:00:00.62   DATE:
17-May-97  14:27:55
 ** J =    8 JLPAR = 1 FINISHED;  CPU: 00:00:00.11   WALL: 00:00:00.78   DATE:
17-May-97  14:27:55
 ** J =   10 JLPAR = 1 FINISHED;  CPU: 00:00:00.14   WALL: 00:00:00.95   DATE:
17-May-97  14:27:55
 ===============================================================================
 **** END OF CALCULATION ****
      MAXIMUM NUMBER OF CHANNELS USED WAS:    9
      TIMING:  ELAPSED 00:00:00.95 / CPU 00:00:00.14
      CURRENT DATE:   17-May-97  14:27:55
 ===============================================================================
 Hibridon> exit
```

Now you wish to extend the calculations up to JTOT2 = 20
 you first reinitiate the calculation

```
ibm4> hib_arn2_151

 --------------------------------------------------------------------------
           HIBRIDON SCATTERING CODE V 4.1 05/17/97 14:23:43 EDT

     AUTHORS: M. ALEXANDER, D. MANOLOPOULOS, H.-J. WERNER, B. FOLLMEG
 CONTRIBUTORS: D. LEMOINE, P. VOHRALIK, G. COREY, R. JOHNSON, T. ORLIKOWSKI
               A. BERNING, A. DEGLI-ESPOSTI, C. RIST, P. DAGDIGIAN, B. POUILLY
               G. VAN DER SANDEN, M. YANG, F. DE WEERD, S. GREGURICK
 --------------------------------------------------------------------------
```

and supply the initial input variables, which you saved before initiating the first job

```
 Hibridon> inp=arn2_test.inp
 NOPRIN = .TRUE., SO IPRINT SET TO -1
```
now extend the value of JTOT2, and set the flag RSFLAG = .TRUE.

```
 Hibridon> jtot2=20
 Hibridon> rsflag=t
```

re-run the calculation

```
 Hibridon> run
 NOPRIN = .TRUE., SO IPRINT SET TO -1
```

The output of the original job is read in (from the file Job.sav)

```
 RESTART DATA FOUND FOR JTOT= 10  JLPAR= 1  NU=  0

 ENERGY=   500.000 NLEVOP=  3

 CHECKING INTERPOLATION DATA IN RESTART FILE:

 READ INTEGRAL CROSS SECTIONS FOR J=  8  JP= 1
 READ PARTIAL  CROSS SECTIONS FOR J= 10  JP= 1
 READ PARTIAL  CROSS SECTIONS FOR J=  8  JP= 1
 READ PARTIAL  CROSS SECTIONS FOR J=  6  JP= 1
 READ INTEGRAL CROSS SECTIONS FOR J= 10  JP= 1

 NO ERRORS DETECTED.
 ** CONTINUE CC CALCULATION AT JTOT=   12, JLPAR=    1
 ** J =   12 JLPAR = 1 FINISHED;  CPU: 00:00:00.03   WALL: 00:00:00.23   DATE:
17-May-97  14:28:14
 ** J =   14 JLPAR = 1 FINISHED;  CPU: 00:00:00.05   WALL: 00:00:00.41   DATE:
17-May-97  14:28:14
 ** J =   16 JLPAR = 1 FINISHED;  CPU: 00:00:00.07   WALL: 00:00:00.57   DATE:
17-May-97  14:28:14
 ** J =   18 JLPAR = 1 FINISHED;  CPU: 00:00:00.10   WALL: 00:00:00.73   DATE:
17-May-97  14:28:14
 ** J =   20 JLPAR = 1 FINISHED;  CPU: 00:00:00.13   WALL: 00:00:00.90   DATE:
17-May-97  14:28:14
 ===============================================================================
 **** END OF CALCULATION ****
      MAXIMUM NUMBER OF CHANNELS USED WAS:    9
      TIMING:  ELAPSED 00:00:00.90 / CPU 00:00:00.13
      CURRENT DATE:   17-May-97  14:28:14
 ===============================================================================
```

Now print out the integral cross sections

```
 Hibridon> printc

 INTEGRAL CROSS SECTIONS READ FROM FILE Job1.ics
 WRITTEN:    17-May-97  14:28:14
 LABEL:          N2-Ar CC PATTENGILL POTENTIAL
 POT NAME:  PATTENGILL-LABUDDE-BERNSTEIN AR-N2

 ** INTEGRAL CROSS SECTIONS; IEN= 1 **
    RMU=  16.4700  E=   500.00  JLPAR= 1  JTOT-1=  0  JTOT-2=  20  JTOT-D=  2
 ** CC CALCULATION, JLPAR= 1 **

 ** COLUMN HEADINGS ARE INTIAL STATES, ROW HEADINGS ARE FINAL STATES **

            J=   0        2        4
   J    I | I=   0        0        0

    0    0  3.89E+00 2.22E-01 6.50E-02
    2    0  1.08E+00 2.84E+00 2.76E-01
    4    0  5.38E-01 4.68E-01 3.05E+00
```

To check these, they can be compared with the results of a full calculation extending from JTOT = 0 to 20

```
 Hibridon> jtot1=0
 Hibridon> run
 NOPRIN = .TRUE., SO IPRINT SET TO -1
 ** J =    0 JLPAR = 1 FINISHED;  CPU: 00:00:00.00   WALL: 00:00:00.13   DATE:
17-May-97  14:28:34
 ** J =    2 JLPAR = 1 FINISHED;  CPU: 00:00:00.02   WALL: 00:00:00.29   DATE:
17-May-97  14:28:34
 ** J =    4 JLPAR = 1 FINISHED;  CPU: 00:00:00.05   WALL: 00:00:00.45   DATE:
17-May-97  14:28:34
 ** J =    6 JLPAR = 1 FINISHED;  CPU: 00:00:00.08   WALL: 00:00:00.61   DATE:
17-May-97  14:28:34
 ** J =    8 JLPAR = 1 FINISHED;  CPU: 00:00:00.11   WALL: 00:00:00.78   DATE:
17-May-97  14:28:35
 ** J =   10 JLPAR = 1 FINISHED;  CPU: 00:00:00.13   WALL: 00:00:00.94   DATE:
17-May-97  14:28:35
 ** J =   12 JLPAR = 1 FINISHED;  CPU: 00:00:00.15   WALL: 00:00:01.10   DATE:
17-May-97  14:28:35
 ** J =   14 JLPAR = 1 FINISHED;  CPU: 00:00:00.18   WALL: 00:00:01.27   DATE:
17-May-97  14:28:35
 ** J =   16 JLPAR = 1 FINISHED;  CPU: 00:00:00.21   WALL: 00:00:01.43   DATE:
17-May-97  14:28:35
 ** J =   18 JLPAR = 1 FINISHED;  CPU: 00:00:00.25   WALL: 00:00:01.60   DATE:
17-May-97  14:28:35
 ** J =   20 JLPAR = 1 FINISHED;  CPU: 00:00:00.28   WALL: 00:00:01.76   DATE:
17-May-97  14:28:35
 ===============================================================================
 **** END OF CALCULATION ****
      MAXIMUM NUMBER OF CHANNELS USED WAS:    9
      TIMING:  ELAPSED 00:00:01.76 / CPU 00:00:00.28
      CURRENT DATE:   17-May-97  14:28:35
 ===============================================================================
 Hibridon> printc

 INTEGRAL CROSS SECTIONS READ FROM FILE Job1.ics
 WRITTEN:    17-May-97  14:28:35
 LABEL:          N2-Ar CC PATTENGILL POTENTIAL
 POT NAME:  PATTENGILL-LABUDDE-BERNSTEIN AR-N2

 ** INTEGRAL CROSS SECTIONS; IEN= 1 **
    RMU=  16.4700  E=   500.00  JLPAR= 1  JTOT-1=  0  JTOT-2=  20  JTOT-D=  2
 ** CC CALCULATION, JLPAR= 1 **

 ** COLUMN HEADINGS ARE INTIAL STATES, ROW HEADINGS ARE FINAL STATES **

            J=   0        2        4
   J    I | I=   0        0        0

    0    0  3.89E+00 2.22E-01 6.50E-02
    2    0  1.08E+00 2.84E+00 2.76E-01
    4    0  5.38E-01 4.68E-01 3.05E+00
```

The results are the same!

```
 Hibridon> quit
 %
```
