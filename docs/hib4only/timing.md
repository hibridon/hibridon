# Timing and Benchmarks

The command `hib_timer` initiates two test close-coupled calculations for benchmarking the speed of your platform.  The test calculations are:
- A calculation for collisions of Ar with N<sub>2</sub>, involving 961 channels with 96 sectors in the log-derivative propagation and 30 sectors in the Airy propagation

- A calculation for collisions of Ar with NO, involving 1982 channels with 18 sectors in the log-derivative
propagation and 12 sectors in the Airy propagation

The command `hib_timer` creates the file `.../hib43/timing/time_output` which contains the wall and cpu times for
calculations at an initial value of the total energy and at subsequent vales of the total energy. For both the Airy and Log-derivative propagation algorithms fewer matrix-matrix operations are required at a subsequent energy.  Typically, the computation time is reduced by a factor of 0.66 (2/3)  in the log-derivative propagation and  a factor of 0.2-0.25 in the Airy propagation.

In addition, version 10 of the Intel MKL math library or (to a lesser extent) the latest version of Apple's Accelerate framework take advantage of parallelism to obtain significant time reductions on systems wih more than one core per cpu.  The file `.../hib43/timing/time_output` reports cpu and wall clock times for each calculation.  The wall time is the actual elapsed time for the job, while the cpu time is the time the job would take on a single cpu. The ratio of <i>t<sub>wall</sub></i> to <i>t<sub>cpu</sub></i> is a  measure of this reduction.

Sample times for the two test calculations are listed in the following tables.  It is clear that execution time depends not only on the clock speed of the processor but also on the size of the available L2 cache.

Also, a significant speed up can be obtained from use of the 64-bit, rather than 32-bit, versions of Intel's compiler and MKL libraries.

## Timing (min:sec); 961 Channel Close-Coupling Ar+N<sub>2</sub>[^2]

| OS                  | Compiler                    | Library                   | CPU                       | speed<br>(GHz) | RAM<br>(GB) | L2 Cache<br>(MB)[^b] | Bus Speed<br>(MHz) | Log-Deriv<br>1<sup>st</sup> Energy | Log-Deriv<br>2<sup>nd</sup> Energy | Airy<br>1<sup>st</sup> Energy | Airy<br>2<sup>nd</sup> Energy |
| ------------------- | --------------------------- | ------------------------- | ------------------------- | -------------- | ----------- | -------------------- | ------------------ | ---------------------------------- | ---------------------------------- | ----------------------------- | ----------------------------- |
| OSX 10.5.2          | ifort 10.1<br>32 bit [^c]   | MKL 10.0.2<br>32 bit [^d] | Core 2 Duo                | 2.4            | 4           | 4                    | 800                | 2:58<br>1:59                       | 1:58<br>1:09                       | 1:28<br>0:53                  | 0:30<br>0:15                  |
| OSX 10.5.2          | ifort 10.1<br>64 bit[^e]    | MKL 10.0.2<br>64 bit [^f] | Core 2 Duo                | 2.4            | 4           | 4                    | 800                | 2:25<br>1:31                       | 1:31<br>0:46                       | 1:14<br>0:42                  | 0:25<br>0:14                  |
| OSX 10.5.2          | ifort 10.1<br>64 bit[^e]    | vecLib 1.4.2              | Core 2 Duo                | 2.4            | 4           | 4                    | 800                | 2:54<br>2:28                       | 1:58<br>1:20                       | 2:02<br>1:35                  | 0:38<br>0:24                  |
| OSX 10.3.9          | xlf 8.1                     | vecLib 1.0.3              | G5                        | 2              | 4           | 0.5                  | 1000               | 5:33<br>5:33                       | 3:43<br>3:43                       | 2:47<br>2:47                  | 0:48<br>0:48                  |
| OSX 10.5.2          | ifort 10.1<br>64 bit[^e]    | MKL 10.0.2<br>64 bit[^f]  | Dual Xeon<br>dual-core    | 2.66           | 5           | 4                    | 1333               | 3:04<br>0:48                       | 2:01<br>0:31                       | 1:25<br>0:24                  | 0:28<br>0:07                  |
| OSX 10.6.5          | ifort 12.0.0<br>64 bit[^e]  | MKL 10.3<br>64 bit[^f]    | Core i7<br>dual-core      | 2.66           | 8           | 4                    | 1067               | 1:52<br>0:59                       | 1:16<br>0:38                       | 0:58<br>0:30                  | 0:19<br>0:10                  |
| OSX 10.7.2          | ifort 12.1.0<br>64 bit [^e] | MKL 10.3<br>64 bit [^f]   | Xeon W3680<br>hexa-core   | 3.33           | 24          | 2                    | 1333               | 1:40<br>0:17                       | 1:05<br>0:11                       | 1:02<br>0:11                  | 0:19<br>0:03                  |
| Linux<br>RedHat 4   | ifort 10.1<br>64 bit[^e]    | MKL 10.0.0<br>64 bit[^g]  | Quad Xeon<br>dual-core    | 2.0            | 4           | 1                    | 1333               | 2:42<br>0:57                       | 2:02<br>1:27                       | 1:47<br>0:30                  | 0:36<br>0:10                  |
| Linux<br>RedHat 4   | ifort 10.1<br>64 bit[^e]    | MKL 10.0.0<br>64 bit[^g]  | Eight Xeon<br>quad-core   | 2.33           | 8           | 3                    | 1333               | 4:17<br>0:41                       | 3:18<br>1:19                       | 5:47<br>1:27                  | 0:50<br>0:12                  |
| Linux<br>RedHat 4   | ifort 10.1<br>64 bit[^e]    | MKL 10.0.0<br>64 bit[^g]  | Quad Xeon<br>dual-core    | 2.66           | 8           | 1.5                  | 1333               | 2:16<br>0:43                       | 1:44<br>1:22                       | 1:26<br>0:24                  | 0:31<br>0:14                  |
| Linux<br>RedHat 4   | ifort 10.1<br>64 bit[^e]    | MKL 10.0.0<br>64 bit[^g]  | Quad Xeon<br>single-core  | 3.0            | 4           | 2                    | 800                | 6:36<br>3:20                       | 4:30<br>1:59                       | 4:25<br>1:10                  | 1:35<br>0:25                  |
| Linux<br>RedHat 4   | ifort 10.1<br>64 bit[^e]    | MKL 10.0.0<br>64 bit[^g]  | Dual Xeon<br>single-core  | 3.2            | 4           | 0.5                  | 800                | 2:40<br>2:50                       | 1:51<br>1:546                      | 1:47<br>1:04                  | 0:37<br>0:19                  |
| Linux<br>Arch 3.1.8 | ifort 12.1.2<br>64 bit      | MKL 10.3<br>64 bit        | Core i7 2600<br>quad-core | 3.4            | 8           | 2                    | 1333               | 0:42<br>0:11                       | 0:29<br>0:07                       | 0:37<br>0:10                  | 0:10<br>0:03                  |


- [^a] : The first and second entries correspond, respectively, to the total cpu time and to the elapsed time. The ratio of the total cpu time to the elapsed time is a measure of the gain due to  parallelism in the LAPACK routines in the indicated math libraries.
- [^b] : Cache per processor
- [^c] : `/opt/intel/fc`
- [^d] : `/Library/Frameworks/Intel_MKL.framework/Versions/Current/lib/32`
- [^e] : `/opt/intel/fce`
- [^f] : `/Library/Frameworks/Intel_MKL.framework/Versions/Current/lib/em64t`
- [^g] : `.../intel/mkl/10.0.011/lib/em64t`

## Timing (min:sec); 1982 Channel Close-Coupling Ar+NO [^a]

| OS                  | Compiler               | Library               | CPU                       | speed<br>(GHz) | RAM<br>(GB) | L2 Cache<br>(MB)[^b] | Bus Speed<br>(MHz) | Log-Deriv<br>1<sup>st</sup> Energy | Log-Deriv<br>2<sup>nd</sup> Energy | Airy<br>1<sup>st</sup> Energy | Airy<br>2<sup>nd</sup> Energy |
| ------------------- | ---------------------- | --------------------- | ------------------------- | -------------- | ----------- | -------------------- | ------------------ | ---------------------------------- | ---------------------------------- | ----------------------------- | ----------------------------- |
| OSX 10.5.2          | ifort 10.1 <br>32 bit  | MKL 10.0.2 <br>32 bit | Core 2 Duo                | 2.4            | 4           | 4                    | 800                | 3:25<br>2:28                       | 2:14<br>1:13                       | 5:19<br>3:02                  | 1:50<br>0:60                  |
| OSX 10.5.2          | ifort 10.1 <br>64 bit  | MKL 10.0.2 <br>64 bit | Core 2 Duo                | 2.4            | 4           | 4                    | 800                | 2:39<br>2:03                       | 1:43<br>0:57                       | 4:42<br>2:41                  | 1:31<br>0:49                  |
| OSX 10.5.2          | ifort 10.1 <br>64 bit  | vecLib 1.4.2          | Core 2 Duo                | 2.4            | 4           | 4                    | 800                | 5:18<br>3:08                       | 2:49<br>1:39                       | 6:47<br>4:25                  | 2:27<br>1:33                  |
| OSX 10.3.9          | xlf 8.1                | vecLib                | G5                        | 2              | 4           | 0.5                  | 1000               | 10:38<br>10:38                     | 6:58<br>6:58                       | 8:43<br>8:43                  | 3:05<br>3:05                  |
| OSX 10.5.2          | ifort 10.1             | MKL 10.0.2            | Dual Xeon<br>dual-core    | 2.66           | 5           | 4                    | 1333               | 3:04<br>1:03                       | 1:58<br>0:34                       | 4:37<br>1:29                  | 1:33<br>0:27                  |
| OSX 10.6.5          | ifort 12.0.0           | MKL 10.3              | Core i7<br>dual-core      | 2.66           | 5           | 4                    | 1067               | 2:04<br>1:10                       | 1:19<br>0:41                       | 3:14<br>1:44                  | 1:09<br>0:36                  |
| OSX 10.7.2          | ifort 12.1.0<br>64 bit | MKL 10.3<br>64 bit    | Xeon W3680<br>hexa-core   | 3.33           | 24          | 2                    | 1333               | 2:08<br>0:29                       | 1:23<br>0:16                       | 2:54<br>0:39                  | 1:08<br>0:13                  |
| Linux<br>RedHat 4   | ifort 10.1<br>64 bit   | MKL 10.0.0<br>64 bit  | Quad Xeon<br>dual-core    | 2.0            | 4           | 4                    | 1333               | 3:51<br>3:18                       | 2:34<br>2:54                       | 5:44<br>2:18                  | 2:01<br>1:24                  |
| Linux<br>RedHat 4   | ifort 10.1<br>64 bit   | MKL 10.0.0<br>64 bit  | Eight Xeon<br>quad-core   | 2.33           | 8           | 1.5                  | 1333               | 5:46<br>1:27                       | 3:56<br>1:45                       | 5:44<br>1:04                  | 2:25<br>0:29                  |
| Linux<br>RedHat 4   | ifort 10.1<br>64 bit   | MKL 10.0.0<br>64 bit  | Quad Xeon<br>dual-core    | 2.66           | 8           | 2                    | 1333               | 3:13<br>1:09                       | 2:08<br>1:43                       | 4:40<br>1:29                  | 1:38<br>0:44                  |
| Linux<br>RedHat 4   | ifort 10.1<br>64 bit   | MKL 10.0.0<br>64 bit  | Quad Xeon<br>single-core  | 3.0            | 4           | 0.5                  | 800                | 9:27<br>2:55                       | 6:17<br>2:28                       | 16:07<br>4:30                 | 5:30<br>1:37                  |
| Linux<br>Arch 3.1.8 | ifort 12.1.2<br>64 bit | MKL 10.3<br>64 bit    | Core i7 2600<br>quad-core | 3.4            | 8           | 2                    | 1333               | 1:08<br>0:22                       | 0:44<br>0:12                       | 1:48<br>0:36                  | 0:35<br>0:09                  |

[^a] The first and second entries correspond, respectively, to the total cpu time and to the elapsed 
time. The ratio of the total cpu time to the elapsed time is a measure of the gain due to  parallelism in the LAPACK
routines in the indicated math libraries.
