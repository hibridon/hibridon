#------------------------------------------------------------------------------
# ch3he test
# CH3(X)--He Dagdigian CCSD(T) PES
#------------------------------------------------------------------------------

set(TEST_ID ch3he)
set(TEST_POT_SRC_FILE "pot_ch3he_ccsdt.F90")
set(TEST_POT_DATA_FILES "ch3he_pot.dat")
set(TEST_COMMAND_FILE "ch3he_test.com")
set(TEST_INPUT_FILES "Ch3he_test.inp")
set(TEST_OUTPUT_FILES "Ch3he1.dcs Ch3he1.ics")
set(TEST_KMAX 150)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")


#------------------------------------------------------------------------------
# ch3he test (Quick version for code coverage)
# CH3(X)--He Dagdigian CCSD(T) PES
#------------------------------------------------------------------------------

set(TEST_ID ch3he_quick)
set(TEST_POT_SRC_FILE "pot_ch3he_ccsdt.F90")
set(TEST_POT_DATA_FILES "ch3he_pot.dat")
set(TEST_COMMAND_FILE "ch3he_test_quick.com")
set(TEST_INPUT_FILES "Ch3he_testq.inp")
set(TEST_OUTPUT_FILES "Ch3heq1.dcs Ch3heq1.ics")
set(TEST_KMAX 150)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
