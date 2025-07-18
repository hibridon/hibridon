The key parameters `TOLAI` and `RINCR` control the interval size in the AIRY propagation.

Whenever R < `RINCR`, then no change is allowed in the interval size.

If, however, R > `RINCR` and `TOLAI` is .ge. 1, then the width of the next interval is set to be a constant multiple of the present interval size, namely:
```
DRNOW(n+1) = DRNOW(n) * TOLAI
```

Thus, to ensure constant step sizes in the AIRY propagation, set `TOLAI` = 1. By setting `TOLAI` > 1, the interval sizes will be monotonically increasing, so that you will use larger step sizes at long range, where the potential presumably varies more slowly. Values of `TOLAI` in the range 1.1 - 1.3 are appropriate for in most calculations.
If, however, the control parameter `TOLAI` is set to be < 1, then, again only if R > `RINCR`, the interval sizes are adjusted using crude estimates of the error intoduced by the linear approximation to the W matrix (see M. H. Alexander, J. Chem. Phys. **81**, 4510). These estimates are designated `COFF` and `CDIAG`. The width of the next interval is adjusted according to the following alogrithm:
```
DRNOW(n+1) = DRNOW(n) * (TOLAI/CMAX) ** (1/3)
```
where CMAX = max(`COFF`,`CDIAG`). If CMAX < `TOLAI`, then the interval width is increased, if CMAX > `TOLAI`, the interval width is decreased.

In this case (`TOLAI` < 1) decreasing the value of `TOLAI` should increase the accuracy of the integration, although there is no quantitive relation between the value of `TOLAI` and the desired accuracy in the S matrix.

The default value is `TOLAI`=1, which forces constant interval sizes in the AIRY propagation.

⚠️ Regardless of the sign of `TOLAI`, no change in the interval width is made if R < `RINCR`

Click [here](Parameters) for more information on reading in, changing, or saving the values of any parameter.