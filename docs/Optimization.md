## Introduction

In order to carry out scattering calculations on a particular system you must first decide on the values of the parameters which control:
  - The starting point of the integration of the coupled equations. This is set by the parameter [RSTART](RSTART,-RENDLD,-and-RENDAI)

  As the total angular momentum (Jtot) increases, the classical turning point will also increase. As this occurs, the program automatically adjusts [RSTART](RSTART,-RENDLD,-and-RENDAI) to remain a constant distance inside the innermost classical turning point.

  - The accuracy of the numerical integration
    - In the LOGD propagation the accuracy is controlled by the parameter [SPAC](SPAC), which sets the step size.
    - In the AIRY propagation the accuracy is controlled by two parameters:
      - [TOLAI](TOLAI-and-RINCR), which determines how large can be the maximum of the diagonal and off-diagonal terms which are neglected
      - [RINCR](TOLAI-and-RINCR), which sets the point beyond which the step size is allowed to increase

  - The point at which the switch is made from the LOGD to the AIRY propagator. This is set by the parameter [RENDLD](RSTART,-RENDLD,-and-RENDAI).
As [RSTART](RSTART,-RENDLD,-and-RENDAI) increases, at subsequent values of Jtot, the value of [RENDLD](RSTART,-RENDLD,-and-RENDAI) increases accordingly, so that the distance covered by the LOGD propagator remains constant. If, however, you are not using the AIRY propagator ([flag](Flags) [AIRYFL](AIRYFL-and-PRAIRY) = .false.), then the program does not change [RENDLD](RSTART,-RENDLD,-and-RENDAI) at higher values of Jtot.

  - The ending point of the integration. This is set by the parameter [RENDAI](RSTART,-RENDLD,-and-RENDAI).
The number of partial waves needed and the required partial wave increment. This is set by the parameters [JTOT2](JTOT) and [JTOTD]JTOT).

  - The number of quantum states required in the channel expansion

To set these parameters, use the [OPTIMIZE](OPTIMIZE) [command](Commands) and proceed as follows:

## Optimization of RSTART and LODG propagation

First, disable the AIRY propagation by setting [AIRYFL](AIRYFL-and-PRAIRY) = .false.
In this step you will adjust the two parameters that control the LOGD integration. OPTIMIZE checks for convergence of the square of T matrix elements. The only levels checked are those which are explicity designated in the level lists [JOUT](JOUT) and [INDOUT](INDOUT). ⚠️ Make sure these parameters are set so that you have selected a representative set of levels for your problem. Set the total angular momentum [JTOT1](JTOT) to be small. Set [JTOT2](JTOT) = [JTOT1](JTOT), and [JTOTD](JTOT) = 1.

- First set [RSTART](RSTART,-RENDLD,-and-RENDAI) sufficiently small that the integration is starting inside the classically forbidden region. Set [RENDLD](RSTART,-RENDLD,-and-RENDAI) to be roughly 2 bohr outside the innermost classical turning point, which you can estimate using the [command](Commands) [TURN](TURN-and-MINPOT). First optimize the parameter [SPAC](SPAC as follows:
    ```
    OPT,SPAC,{valmax},{valmin},0.7
    ```
    where {valmax} is the initial value, typically 0.5 times the deBroglie wavelength (which can be found using the [command](Commands) [DEBROGLI](DEBROGLI)), and {valmin} is typically 0.05*valmax.

- When you have found a reasonable value of [SPAC](SPAC) save the parameters (using the [command](Commands) [SAVE](SHOW,-SAVE,-and-READ)) and then proceed to adjust the value of [RSTART](RSTART,-RENDLD,-and-RENDAI) as follows:
    ```
    OPT,RSTART,{valmax},{valmin},1,-0.25
    ```
  where now 
    - {valmax} = innermost turning point - 0.2 bohr 
    - {valmin} = {valmax} - 1.5 bohr
    when this is complete, [save](save) the optimized value of [RSTART](RSTART,-RENDLD,-and-RENDAI).


## Optimization of AIRY propagation

Now, you can adjust the parameters that control the AIRY integration.
Set [AIRYFL](AIRYFL-and-PRAIRY) = .true., [TOLAI](TOLAI-and-RINCR) = 1.0, and [RENDAI](RSTART,-RENDLD,-and-RENDAI) equal to a value at which you expect the coupling potential to be asymptotically weak (for problems with long-range potentials that vary only as R-6, RENDAI should be typically 20 - 30 bohr; the value of the potential at any value of R can be ascertained using the command [TESTPOT](TESTPOT)).

- First adjust FSTFAC, which controls the relative size of the AIRY steps to the LOGD steps, as follows:
 ```
  OPT,FSTFAC,{valmax},{valmin},0.8
  ```
with, typically, {valmax} = 20, {valmin} = 2

- Then adjust [RENDAI](RSTART,-RENDLD,-and-RENDAI) as follows:
  ```
  OPT,RENDAI,{valmin},{valmax},1,3
  ```
where, typically, {valmin} = 15, {valmax} = 40 . For long-range potentials (dipole-dipole, charge-quadrupole, charge-dipole coupling these values will have to be much larger).
Then adjust [TOLAI](TOLAI-and-RINCR) as follows:
  ```
  OPT,TOLAI,{valmax},1,1,{step}
  ```
where, typically {valmax} = 1.5 - 2 and {step} = - 0.05. [TOLAI](TOLAI-and-RINCR) controls the rate at which step sizes increase in the AIRY propagation. This is very system dependent; typical values are 1 ¾ TOLAI ¾ 1.2.

## Location of LOGD/AIRY switch

At this point, you need to set [RENDLD](RSTART,-RENDLD,-and-RENDAI), the integration switching point. The goal is to minimize the cpu time spent in propagating the log-derivative matrix from RSTART to RENDAI. Since the LOGD propagator uses a constant step size, its “speed” - how long it takes to propagate the log-derivative matrix over a given distance - is constant. For a given sector width, the AIRY propagator involves significantly more matrix operations and is hence slower than the LOGD propagator. In general the speed ratio is about 6:1 in favor of the LOGD propagator on scalar machines, but drops to less than 3:1 on vector machines with fast matrix multiply library routines. However, in the AIRY propagation the interval widths increase with R, so that the speed becomes faster as R increases. Eventually, the AIRY propagator will become faster than the LOGD propagator, which uses constant step sizes.

The total time will be proportional to
```
NL + 6 NA ($N_L$ + 3 $N_A$, on vector machines),
```

where $N_L$ designates the number of steps in the LOGD propagation and $N_A$ designates the number of intervals in the AIRY propagation. For a given value of the parameter spac, $N_L$ will be directly proportional to [RENDLD](RSTART,-RENDLD,-and-RENDAI). Thus to minimize the total cpu time you need to decrease [RENDLD](RSTART,-RENDLD,-and-RENDAI) (thereby decreasing $N_L$) while at the same time trying to minimize the total ($N_L$ + f * $N_A$), where f=6 for scalar machines but 3 for vector machines). You can accomplish this as follows:

The scattering program prints out $(N_L)$ and $N_A$. First decrease [RENDLD](RSTART,-RENDLD,-and-RENDAI), reoptimize the AIRY parameters, as described above, and see if $N_L$ + f * $N_A$ (or the total cpu time) decreases. If so, continue to decrease [RENDLD](RSTART,-RENDLD,-and-RENDAI). Otherwise, increase [RENDLD](RSTART,-RENDLD,-and-RENDAI).

In general, you don't need to search for a precise minimum in the function ($N_L$ + f * $N_A$) since the number of intervals used in the AIRY propagation will vary substantially at different values of Jtot. For systems with small deBroglie wavelengths (high mass, high collision energy) the minimal cpu time will correspond to using the LOGD propagator over only a short range (0.5 -2 Bohr). For systems with large deBroglie wavelengths (collisions involving H2 or He, low collision energy), you will even want to use the LOGD propagator over a longer range. In some cases you may want to use only the LOGD integrator (set [AIRYFL](AIRYFL-and-PRAIRY) = .false.) , or only the AIRY integrator (set [LOGDFL](LOGDFL-and-PRLOGD) = .false.).


## Determination of number of partial waves

At this point, if you are treating a collision problem, you need to determine how many partial waves are required to optain convergence in the differential or integral cross sections of interest. This is governed by the value of the parameter [JTOT2](JTOT)

- In the case of differential cross sections, use a representative channel basis, but smaller than that necessary for convergence, and a large value of [JTOT2](JTOT)

  Determine the differential cross sections for several transitions of interest, over a coarse angular grid, using the [DIFCRS](DIFCRS) command. Successively increase the parameter jtotend in the call list to [DIFCRS](DIFCRS) until convergence is reached.

- In the case of integral cross sections, use a representative channel basis, but smaller than that necessary for convergence, and a large value of [JTOT2](JTOT), and a moderate value of [JTOTD](JTOT) (JTOTD = 3 - 5).

Determine the partial inelastic (or elastic) cross sections for several transitions of interest, using the [PARTC](PARTC) command. For the output (which will appear at the console as well as in the file {jobname}.psc, you will be able to judge how large a value of [JTOT2](JTOT) is necessary.

The, decrease [JTOTD](JTOT) until you have reached convergence in the integral cross sections.

⚠️ Remember, convergence of integral or differential cross sections to within 1% is usually well within the errors bars of most (if not all) experiments.

⚠️ In the case of [photodissociation](Photodissociation), or collisions with [surfaces](FLAGSU), only one partial wave is involved. Consequently, this optimization step can be skipped.


## Size of channel basis

Once you have decided on the parameters which control the integration and the parameters which control the partial wave sum, you need only decide on the size of the channel basis. The number of channels is, of course, dependnet on the particular system and/or [BASIS](basis) subroutine selected. In many cases, this is controled by the parameter JMAX. By following a similar procedure outlined in step #4 above (optimization of the partial wave sum), one can decide how large a channel basis is necessary to obtain convergence.

Unfortunately, since the cpu requirement of any close-coupled calculation goes up as the thirdpower of the number of channels, the choice of the maximum size of the channel basis may involve significant computational effort.

⚠️ Remember, convergence of integral or differential cross sections to within 1% is usually well within the errors bars of most (if not all) experiments.