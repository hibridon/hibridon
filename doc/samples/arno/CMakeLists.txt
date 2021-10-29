# this is a minimal cmake example for creating a hibridon executable with a user-defined potential
cmake_minimum_required (VERSION 3.3)
project (myarno)
enable_language (Fortran)

# add hibridon library
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../../.." hibridon)

# declare a new executable using hibridon's add_hibexe cmake function, where:
# - the 1st argument (here myarnoexe) is the name of the resulting executable
# - the 2nd argument (here "pot_arno.F90") is the name of the user provided potential file
# - the 3rd argument (here "kmax") is the size of the t matrix
#    - kmax : for normal cases
#    - kbig : for special cases (only arn2_big test uses it)
add_hibexe(myarnoexe "pot_arno.F90" "kmax")



