# LOGDFL

The flag `LOGDFL` should be set
  - .TRUE. if the Logd integrator is to be used
  - .FALSE. if the Logd integrator is not to be used

If `LOGDFL` is .FALSE., then [`RENDLD`](./rendld.md) is set equal to [`RSTART`](./rstart.md).

⚠️ Since the width of the first interval in the Airy propagation is equal to [`SPAC`](../param/spac.md) x [`FSTFAC`](../param/fstfac.md), even if LOGDFL is set .FALSE., SPAC should be given a non-zero value.

# PRLOGD

If the flag `PRLOGD` is set .TRUE., then the log-derivative matrix is printed out after both the Logd and Airy propagation steps. Also the $K$ matrix is printed out prior to the determination of the $S$ matrix.

⚠️ Do not set [`AIRYFL`](./airyfl.md) and [`LOGDFL`](./logdfl.md) both .FALSE.