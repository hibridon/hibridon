# FSTFAC

The parameter FSTFAC is the ratio between the first interval size in the AIRY integration, and the value of [`SPAC`](../param/spac.md), the step-size in the LOGD integration. Specifically,
```
DRNOW (first interval) = FSTFAC * SPAC
```
The default value is FSTFAC=1.