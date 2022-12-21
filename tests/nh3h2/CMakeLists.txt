#------------------------------------------------------------------------------
# nh3h2 test
# NH3--H2 2009 Potential
#------------------------------------------------------------------------------

set(TEST_ID nh3h2)
set(TEST_POT_SRC_FILE "pot_nh3h2_2009.F90")
set(TEST_POT_DATA_FILES "fitvij_bf_62.h2")
set(TEST_COMMAND_FILE "nh3h2_test.com")
set(TEST_INPUT_FILES "Nh3h2_po.inp")
set(TEST_OUTPUT_FILES "Job1.ics")
set(TEST_KMAX 390)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

#------------------------------------------------------------------------------
# issue188 test
# non physical test to detect regressions in sav file mechanism
# based on unpublished h3ohe potential
#------------------------------------------------------------------------------

set(TEST_ID issue188)
set(TEST_POT_SRC_FILE "pot_nh3h2_2009.F90")
set(TEST_POT_DATA_FILES "fitvij_bf_62.h2")
set(TEST_COMMAND_FILE "issue188_test.com")
set(TEST_INPUT_FILES "Issue188.inp")
set(TEST_OUTPUT_FILES "I188_1.ics I188_1.xsc")
set(TEST_KMAX 600)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
