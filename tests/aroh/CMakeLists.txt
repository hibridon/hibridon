#------------------------------------------------------------------------------
# arohx_ump4 test
#------------------------------------------------------------------------------

set(TEST_ID aroh)
set(TEST_POT_SRC_FILE "pot_arohx_ump4.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "aroh_jatest.com")
set(TEST_INPUT_FILES "Aroh_jatest.inp")
set(TEST_OUTPUT_FILES "Aroh_new1.ics")
set(TEST_KMAX 301)
set(TEST_T_MATRIX_SIZE kmax)

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")


#------------------------------------------------------------------------------
# arohx_ump4 test (quick version for code coverage)
#------------------------------------------------------------------------------

set(TEST_ID aroh_quick)
set(TEST_POT_SRC_FILE "pot_arohx_ump4.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "aroh_jatest_quick.com")
set(TEST_INPUT_FILES "Aroh_jatestq.inp")
set(TEST_OUTPUT_FILES "Aroh_q1.ics")
set(TEST_KMAX 301)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
