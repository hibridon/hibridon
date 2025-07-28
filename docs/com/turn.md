# TURN

The command `TURN` invokes a simple search procedure to find the classical turning point for the potential which is defined by the potential subroutine which is linked in the currently executing scattering program. This information is useful in the determination of the starting point for the integration of the coupled equation, [`RSTART`]RSTART,-RENDLD,-and-RENDAI).

The turning point is estimated using the isotropic term in the potential and a centrifugal barrier of
```
JTOT1 * (JTOT1 + 1) / (2 * XMU * R^2),
```
