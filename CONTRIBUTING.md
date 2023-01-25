## FORTRAN STANDARD USED
Fortran 2008 (maximum)

## CODING CONVENTIONS USED IN THIS PROJECT
- Implicit none everywhere !
- Indentations: two spaces
- Gotos written as: `goto` (not `go to`)
- `end <block> <block name>` at the end of a subroutine/function/module.
  e.g.:
  ```fortran
  subroutine toto
  end subroutine toto
  ```
- Lowercase for everything except for preprocessing which would be uppercase
- Extended types are suffixed by `_type` (e.g.: `grovec_type`)
- Constants (declared as parameter) prefixed with `k_`
- Enums are prefixed with the name of the enum
  e.g.:
  ```fortran
  enum, bind(c) :: color
    enumerator :: color_red = 4
    enumerator :: color_blue = 9
    enumerator :: color_yellow = 12
  end enum
  
  type(color) :: mycolor = color_blue
  ```
- Use the `intent` statement when possible.
- Prefix module names by `mod_`
- When using a module, precise the variable used:
  ```fortran
  use mod_mymodule, only: myvariable
  ```
- Add a space after a comma

  
## REMOVE OBSOLETE CONSTRUCTS
- Gotos should be removed
- Variable declaration using double colons and following:
  ```fortran
  character(len=22) :: mystring
  integer :: myinteger
  real(8) :: myrealarray(10)
  ```
- Replace operators by modern equivalents when possible (e.g.: `.eq.`->`==`)

- Replace old labeled loops: 
  ```fortran
        do 100 i=istart,ilast,istep
  100    isum = isum + i
  ```
  Should be replaced by:
  ```fortran
  do i = istart, ilast, istep
    isum = isum +1
  end do
  ```
 
