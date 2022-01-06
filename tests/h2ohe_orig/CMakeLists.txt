#------------------------------------------------------------------------------
# h2ohe_orig test
#------------------------------------------------------------------------------

set(TEST_ID h2ohe_orig)
set(TEST_POT_SRC_FILE "pot_h2ohe.F90")
set(TEST_POT_DATA_FILES "h2o_coefd.dat h2o_params.dat")
set(TEST_COMMAND_FILE "h2ohe_test.com")
set(TEST_INPUT_FILES "H2ohe_test.inp")
set(TEST_OUTPUT_FILES "H2ohe1.ics")
set(TEST_KMAX 45)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
