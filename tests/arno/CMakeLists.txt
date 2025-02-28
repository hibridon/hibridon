#------------------------------------------------------------------------------
# arno test
#------------------------------------------------------------------------------

set(TEST_ID arno)
set(TEST_POT_SRC_FILE "pot_arno.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "arno_test.com")
set(TEST_INPUT_FILES "Arno_test.inp")
set(TEST_OUTPUT_FILES "Arno_tes1.ics Arno_tes2.ics Arno_tes3.ics")
set(TEST_KMAX 191)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
