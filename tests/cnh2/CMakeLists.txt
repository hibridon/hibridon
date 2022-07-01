#------------------------------------------------------------------------------
# cnph2 test
#------------------------------------------------------------------------------

set(TEST_ID cnph2)
set(TEST_POT_SRC_FILE "pot_cnh2_lique.F90")
set(TEST_POT_DATA_FILES "pot_cnh2_lique.dat")
set(TEST_COMMAND_FILE "cnph2.com")
set(TEST_INPUT_FILES "Cnph2.inp")
set(TEST_OUTPUT_FILES "E10p1.xxsc")
set(TEST_KMAX 300)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

#------------------------------------------------------------------------------
# cnoh2 test
#------------------------------------------------------------------------------

set(TEST_ID cnoh2)
set(TEST_POT_SRC_FILE "pot_cnh2_lique.F90")
set(TEST_POT_DATA_FILES "pot_cnh2_lique.dat")
set(TEST_COMMAND_FILE "cnoh2.com")
set(TEST_INPUT_FILES "Cnoh2.inp")
set(TEST_OUTPUT_FILES "E10o1.xxsc")
set(TEST_KMAX 300)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")