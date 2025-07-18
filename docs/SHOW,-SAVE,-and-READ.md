The command

`SHOW`

causes the current values of all parameters and flags to be displayed on the screen and written to the main [output](output) file.

---

Changes made to the currently accessed input file can be saved using the command

`SAVE`


Similarly, the current set of parameters and flags can be saved to a new file by the command

`SAVE={filename}`


{filename} can be lower or upper case, but is converted internally so that the first character is upper case and all subsequent characters are lower case.

---

Initially, the default [input](input) file is named Inpt

The `SAVE` command is useful for creating updated input data files.

The command 

`READ`

causes all parameters and flags to be reread from the current input file.

The command

`READ,{filename}`

causes all parameters and flags to be reread from the file {filename}.

This command is useful if you have changed any parameters or flags and wish to restore a previous set