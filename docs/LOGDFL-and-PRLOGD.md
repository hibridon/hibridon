The flag LOGDFL should be set
  - .TRUE. if the Logd integrator is to be used
  - .FALSE. if the Logd integrator is not to be used

If LOGDFL is .FALSE., then RENDLD is set equal to RSTART.

⚠️ Since the width of the first interval in the Airy propagation is equal to [SPAC](SPAC) x [FSTFAC](FSTFAC), even if LOGDFL is set .FALSE., SPAC should be given a non-zero value.