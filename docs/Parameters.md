The scattering calculation is controlled by GENERAL (system independent) and [`SYSTEM`](basis) parameters (those which are specific to a particular collision system).

More information on the collision systems which can be treated is given [here](basis).

The system independent parameters are:
- [`ENERG`](./com/energ.md)
- [`FSTFAC`](./param/fstfac.md)
- [`INDOUT`](./com/indout.md)
- [`IPRINT`](./param/iprint.md)
- [`J1J2`](./com/j1j2.md)
- [`JOUT`](./com/jout.md)
- [`JTOT1`](./param/jtot1.md) 
- [`JTOT2`](./param/jtot2.md)
- [`JTOTD`](./param/jtotd.md)
- [`JLPAR`](./param/jlpar.md) 
- [`LSCREEN`](./param/lscreen.md)
- [`NERG`](./com/nerg.md)
- [`NUMAX`](./param/numax.md)
- [`NUMIN`](./param/numin.md) 
- [`NUD`](./param/nud.md)
- [`RCUT`](./param/rcut.md)
- [`READPT`](./param/readpt.md)
- [`RENDAI`](./param/rendai.md)
- [`RENDLD`](./param/rendld.md)
- [`RINCR`](./param/rincr.md) 
- [`RSTART`](./param/rstart.md)
- [`SPAC`](./param/spac.md)
- [`TOLAI`](./param/tolai.md)
- [`XMU`](./param/xmu.md)


To change values of any parameters enter the command
```
PARAMETER1=VAL1,PARAMETER2=VAL2 etc.
```
Different parameters can be separated by commas or by semicolons.

Note that values for the parameters [`ENERG`](./com/energ.md), [`INDOUT`](./com/indout.md), or [`JOUT`](./com/jout.md) must be terminated by a backslash (\) if other PARAMETER=VAL strings follow, e.g.
```
ENERG=1000,2000\AIRYFL=T
```

The Hibridon input parser is not case sensitive, so parameters can be entered in either upper or lower case.

Parameters may be abbreviated, although you will be asked to re-enter the parameter if there is an ambiguity.

The [default](Defaults) parameters are those for the collision of N2 with Ar, which is used as a [test case](tests).