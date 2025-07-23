The scattering calculation is controlled by GENERAL (system independent) and [SYSTEM](Basis-subroutines) parameters (those which are specific to a particular collision system).

More information on the collision systems which can be treated is given [here](Basis-subroutines).

The system independent parameters are:
- [ENERG](ENERG)
- [FSTFAC](FSTFAC)
- [INDOUT](INDOUT)
- [IPRINT](IPRINT)
- [J1J2](J1J2)
- [JOUT](JOUT)
- [JTOT1](JTOT1) 
- [JTOT2](JTOT2)
- [JTOTD](JTOTD)
- [JLPAR](JLPAR) 
- [LSCREEN](LSCREEN)
- [NERG](NERG)
- [NUMAX](NUMAX)
- [NUMIN](NUMIN) 
- [NUD](NUD)
- [RCUT](RCUT)
- [RENDAI](RSTART,-RENDLD,-and-RENDAI)
- [RENDLD](RSTART,-RENDLD,-and-RENDAI)
- [RINCR](RINCR) 
- [RSTART](RSTART,-RENDLD,-and-RENDAI)
- [SPAC](SPAC)
- [TOLAI](TOLAI)
- [XMU](XMU)


To change values of any parameters enter the command
```
PARAMETER1=VAL1,PARAMETER2=VAL2 etc.
```
Different parameters can be separated by commas or by semicolons.

Note that values for the parameters [ENERG](ENERG), [INDOUT](INDOUT), or [JOUT](JOUT) must be terminated by a backslash (\) if other PARAMETER=VAL strings follow, e.g.
```
ENERG=1000,2000\AIRYFL=T
```

The Hibridon input parser is not case sensitive, so parameters can be entered in either upper or lower case.

Parameters may be abbreviated, although you will be asked to re-enter the parameter if there is an ambiguity.

The [default](Defaults) parameters are those for the collision of N2 with Ar, which is used as a [test case](tests).