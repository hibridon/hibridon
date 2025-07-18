The scattering calculation is controlled by a number of logical flags
To change values of flags, enter the command
```
FLAG1=VAL1,FLAG2=VAL2 etc.
```
different flags can be separated by commas or by semicolons.

The Hibridon input driver is not case sensitive, so flags can be entered in either upper or lower case. Flags may be abbreviated, although you will have to expand the abbreviation the flag if there is an ambiguity.

The default values of the various flags are:

  - [AIRYFL](AIRYFL-and-PRAIRY)=T 
  - [BASTST](BASTST)=F
  - [BATCH](BATCH)=F
  - [BOUNDC](BOUNDC)=F
  - [CHLIST](CHLIST)=T
  - [CSFLAG](CSFLAG)=F
  - [FLAGHF](FLAGHF)=F
  - [FLAGSU](FLAGSU)=F
  - [IHOMO](IHOMO-and-TWOMOL)=T
  - [IPOS](IPOS)=F
  - [LOGDFL](LOGDFL-and-PRLOGD)=T
  - [NOPRIN](NOPRIN)=F
  - [NUCROS](NUCROS)=F
  - [PHOTOF](PHOTOF)=F
  - [PRAIRY](AIRYFL-and-PRAIRY)=F
  - [PRLOGD](LOGDFL-and-PRLOGD)=F
  - [PRPART](PRPART-and-WRPART)=F
  - [PRSMAT](PRSMAT-and-WRSMAT)=F
  - [PRT2](PRT2-and-T2TEST)=T
  - [PRXSEC](PRXSEC-and-WRXSEC)=T
  - [READPT](READPT)=F
  - [RSFLAG](RSFLAG)=F
  - [T2TEST](PRT2-and-T2TEST)=F
  - [TWOMOL](IHOMO-and-TWOMOL)=F
  - [WAVEFL](WAVEFL)=F
  - [WRPART](PRPART-and-WRPART)=F
  - [WRSMAT](PRSMAT-and-WRSMAT)=T
  - [WRXSEC](PRXSEC-and-WRXSEC)=F