The Hibridon code writes all output, except for the specialized files discussed in the in the [JOB](Files) help entry, to a main output file (fortran logical unit 9, a sequential formatted file). In addition some output is directed to stdout (fortran logical unit 6). The name of the main output file is initally set to Outpt. This name can be changed by resetting the variable OUT, as follows:
```
OUT={outname},
```
where {outname} is a character string (not in quotes). This name is converted internally to lower case, with the first character upper case. 

⚠️ On some machines existing output files will be overwritten unless you give them a new name.