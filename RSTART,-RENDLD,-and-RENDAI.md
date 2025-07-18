RSTART is the starting point of the integration. If the flag [LOGDFL](LOGDFL-and-PRLOGD) = .FALSE., then the entire integration, from RSTART to RENDAI, is done using the Airy propagator. Normally, [LOGDFL](LOGDFL-and-PRLOGD) = .TRUE., and the log-derivative propagator is used from RSTART to RENDLD, and then the Airy propagator from RENDLD to RENDAI

RENDLD is the step size (sector width) in the log-derivative integration. If [LOGDFL](LOGDFL-and-PRLOGD) = .false., then RENDLD is set equal to RSTART.

RENDAI is the end point for the integration, the point at which the S matrix is determined. If the [flag](Flags) [AIRYFL](AIRYFL-and-PRAIRY) = .FALSE., then RENDAI is set equal to RENDLD.

The parameter [FSTFAC](FSTFAC) controls the ratio between the step size in the log-derivative integration (given by the parameter [SPAC](SPAC)) and the step size in the Airy integration.

In [bound-state](Bound-States) calculations (flag [BOUNDC](BOUNDC) = .TRUE.), then RSTART and RENDAI define the range of the R-dependent wavefunction.