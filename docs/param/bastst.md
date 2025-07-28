# FLAG BASTST

Setting the flag `BASTST` = .TRUE. allows the user to check the output of the subroutine <a href=basis.html>BASIS</a> without carrying out a full integration of the coupled equations. 

In addition, running with `BASTST` = .TRUE. allows you to print quickly all the internal states used in the expansion of the scattering (or photodissociation) wavefunction, and their energies.  This is particularly useful in determining the asymptotic correlation of the adiabatic energies;  see the [`EADIAB`](../com/eadiab.md) help file.

If `BASTST` = .TRUE., then following the command [`RUN`](../com/run.md), the program returns control to the user after calling the subroutine [`BASIS`](../basis.md).

---

⚠️ If the FLAG NOPRIN is set .TRUE. the output from `BASTST` might not be complete
