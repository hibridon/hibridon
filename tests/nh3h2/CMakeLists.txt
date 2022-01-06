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

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
