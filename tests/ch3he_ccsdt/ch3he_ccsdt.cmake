#------------------------------------------------------------------------------
# ch3he_ccsdt test
#------------------------------------------------------------------------------

set(TEST_ID ch3he_ccsdt)
set(TEST_POT_SRC_FILE "pot_ch3he_ccsdt.F")
set(TEST_POT_DATA_FILES "ch3he_pot.dat")
set(TEST_COMMAND_FILE "ch3he_test.com")
set(TEST_INPUT_FILES "Ch3he_test.inp")
set(TEST_OUTPUT_FILES "Ch3he1.ics")
set(TEST_KMAX 151)
set(TEST_T_MATRIX_SIZE kmax)

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_KMAX}" "${TEST_T_MATRIX_SIZE}")
