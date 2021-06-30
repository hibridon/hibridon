#------------------------------------------------------------------------------
# vfit test
#------------------------------------------------------------------------------

set(TEST_ID vfit)
set(TEST_POT_SRC_FILE "pot_vfit.F")
set(TEST_POT_DATA_FILES "FOLLMEG.BIN")
set(TEST_COMMAND_FILE "vfit_test.com")
set(TEST_INPUT_FILES "N2phetest.inp")
set(TEST_OUTPUT_FILES "Vfit_test1.ics")
set(TEST_KMAX 36)
set(TEST_T_MATRIX_SIZE kmax)

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_KMAX}" "${TEST_T_MATRIX_SIZE}")
