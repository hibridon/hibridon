#------------------------------------------------------------------------------
# These test the basis subroutines using the BASTST flag
# For each basis, consistancy of the basis set is tested
# as well as the saving of the basis parameters to an input file
#------------------------------------------------------------------------------


set(TEST_POT_SRC_FILE "pot_bastst.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_KMAX 500)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")

# (BASIS 06 is not supported)
set(supported_bases 01 02 03 04 05 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)

foreach(ibasty IN LISTS supported_bases)
  set(TEST_ID bastst_${ibasty})
  set(TEST_INPUT_FILES  "Bastst_${ibasty}.inp")
  set(TEST_OUTPUT_FILES "Bastst_${ibasty}.out Save_${ibasty}.inp")
  set(TEST_COMMAND_FILE "bastst_${ibasty}.com")
  add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
endforeach()
