The scattering calculation is controlled by a number of logical flags
To change values of flags, enter the command
```
FLAG1=VAL1,FLAG2=VAL2 etc.
```
different flags can be separated by commas or by semicolons.

The Hibridon input driver is not case sensitive, so flags can be entered in either upper or lower case. Flags may be abbreviated, although you will have to expand the abbreviation the flag if there is an ambiguity.

The default values of the various flags are:

  - [`AIRYFL`](./param/airyfl.md)=T 
  - [`BASTST`](./param/bastst.md)=F
  - [`BATCH`](./com/batch.md)=F
  - [`BOUNDC`](./param/boundc.md)=F
  - [`CHLIST`](./param/chlist.md)=T
  - [`CSFLAG`](./param/csflag.md)=F
  - [`FLAGHF`](./param/flaghf.md)=F
  - [`FLAGSU`](./param/flagsu.md)=F
  - [`IHOMO`](./param/ihomo.md)=T
  - [`IPOS`](./param/ipos.md)=F
  - [`LOGDFL`](./param/logdfl.md)=T
  - [`NOPRIN`](./param/noprin.md)=F
  - [`NUCROS`](./param/nucros.md)=F
  - [`PHOTOF`](./param/photof.md)=F
  - [`PRAIRY`](./param/prairy.md)=F
  - [`PRLOGD`](./param/prlogd.md)=F
  - [`PRPART`](./param/prpart.md)=F
  - [`PRSMAT`](./param/prsmat)=F
  - [`PRT2`](./param/prt2.md)=T
  - [`PRXSEC`](./param/prxsec.md)=T
  - [`READPT`](./param/readpt.md)=F
  - [`RSFLAG`](./param/rsflag.md)=F
  - [`T2TEST`](./param/t2test.md)=F
  - [`TWOMOL`](./param/twomol.md)=F
  - [`WAVEFL`](./param/wavefl.md)=F
  - [`WRPART`](./param/wrpart.md)=F
  - [`WRSMAT`](./param/wrsmat.md)=T
  - [`WRXSEC`](./param/wrxsec.md)=F