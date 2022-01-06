#------------------------------------------------------------------------------
# nh3he test
#------------------------------------------------------------------------------

set(TEST_ID nh3he)
set(TEST_POT_SRC_FILE "pot_nh3he_wheatley.F90")
set(TEST_POT_DATA_FILES "pot_nh3he_wheatley_ylm.txt")
set(TEST_COMMAND_FILE "nh3he.com")
set(TEST_INPUT_FILES "Nh3he_para.inp")
set(TEST_OUTPUT_FILES "Job1.ics Job1.xsc Job1.pcs")
set(TEST_KMAX 200)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
