# Parameters defining integration ranges

## RSTART

RSTART is the starting point of the integration. If the flag [`LOGDFL`](./logdfl.md) = .FALSE., then the entire integration, from RSTART to RENDAI, is done using the Airy propagator. Normally, [`LOGDFL`](./logdfl.md) = .TRUE., and the log-derivative propagator is used from RSTART to RENDLD, and then the Airy propagator from RENDLD to RENDAI

## RENDLD

RENDLD is the step size (sector width) in the log-derivative integration. If [`LOGDFL`](./logdfl.md) = .false., then RENDLD is set equal to RSTART.


## RENDAI

RENDAI is the end point for the integration, the point at which the S matrix is determined. If the [flag](../Flags.md) [`AIRYFL`](./airyfl.md) = .FALSE., then RENDAI is set equal to RENDLD.

## notes
The parameter [`FSTFAC`](../param/fstfac.md) controls the ratio between the step size in the log-derivative integration (given by the parameter [`SPAC`](../param/spac.md)) and the step size in the Airy integration.

In [bound-state](../Bound-States.md) calculations (flag [`BOUNDC`](../param/boundc.md) = .TRUE.), then RSTART and RENDAI define the range of the R-dependent wavefunction.