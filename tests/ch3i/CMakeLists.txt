#------------------------------------------------------------------------------
# ch3i test
#------------------------------------------------------------------------------

set(TEST_ID ch3i)
set(TEST_POT_SRC_FILE "pot_ch3i.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "ch3i_test.com")
set(TEST_INPUT_FILES "Ch3i.inp")
set(TEST_OUTPUT_FILES "Ch3itest.flx Ch3itest.psi")
set(TEST_KMAX 100)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

#------------------------------------------------------------------------------
# ch3i test with scattering
#------------------------------------------------------------------------------

set(TEST_ID ch3iscat)
set(TEST_POT_SRC_FILE "pot_ch3i.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "ch3iscat_test.com")
set(TEST_INPUT_FILES "Ch3iscat.inp")
set(TEST_OUTPUT_FILES "Ch3iscat.flx")   #  TODO: add Ch3iscat.psi to the checked files
set(TEST_KMAX 100)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")