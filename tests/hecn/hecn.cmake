#------------------------------------------------------------------------------
# hecn test
#------------------------------------------------------------------------------

set(TEST_ID hecn)
set(TEST_POT_SRC_FILE "pot_hecn_dgels.F")
set(TEST_POT_DATA_FILES "hecn_fitmlv.dat")
set(TEST_COMMAND_FILE "hecnx_hyp.com")
set(TEST_INPUT_FILES "Hecnx.inp")
set(TEST_OUTPUT_FILES "Hecn1.ics")
set(TEST_KMAX 81)
set(TEST_T_MATRIX_SIZE kmax)

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_KMAX}" "${TEST_T_MATRIX_SIZE}")

