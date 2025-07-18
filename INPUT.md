All the [parameters](Parameters) needed by the Hibridon code are assigned default values appropriate to the rotationally inelastic collision of Ar with N2. You can either redefine the full set of input parameters by setting them one by one as described in more detail in the [parameters](Parameters) help file, or read in a previously defined set of parameters, as described in more detail below.

Once the general and system parameters have been set, then can be saved using the [SAVE](SHOW,-SAVE,-and-READ) command
```
SAVE,{inputfile}
```
where {inputfile} designates the name of the file in which the parameters will be stored. The variable {inputfile} can be lower or upper case, but is converted internally so that the first character is upper case and all subsequent characters are lower case.
Executing the command
```
SAVE
```
with {inputfile} left blank, will save (overwrite) the current parameters to the current input file. The default value of the variable {inputfile} is Inpt.

To use a preexisting set of parameters in a subsequent calculation, execute the command
```
INPUT={inputfile}
```
The command [READ](SHOW,-SAVE,-and-READ):
```
READ
```
causes all parameters and flags to be reread from the current input file.

The command
```
READ = {inputfile}
```
causes all parameters and flags to be reread from the current input file.
The [READ](SHOW,-SAVE,-and-READ) command is useful if you have changed any parameters or flags and wish to restore a previous set