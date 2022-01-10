#------------------------------------------------------------------------------
# hecn test
#------------------------------------------------------------------------------

set(TEST_ID hecn)
set(TEST_POT_SRC_FILE "pot_hecn_dgels.F90")
set(TEST_POT_DATA_FILES "hecn_fitmlv.dat")
set(TEST_COMMAND_FILE "hecnx_hyp.com")
set(TEST_INPUT_FILES "Hecnx.inp")
set(TEST_OUTPUT_FILES "Hecn1.ics")
set(TEST_KMAX 81)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

#------------------------------------------------------------------------------
# hecn test (quick version for code coverage)
#------------------------------------------------------------------------------

set(TEST_ID hecn_quick)
set(TEST_POT_SRC_FILE "pot_hecn_dgels.F90")
set(TEST_POT_DATA_FILES "hecn_fitmlv.dat")
set(TEST_COMMAND_FILE "hecn_quick.com")
set(TEST_INPUT_FILES "Hecn_quick.inp")
set(TEST_OUTPUT_FILES "Hecnq1.ics Hecnq1.pcs Hecnq1.hfx")
set(TEST_KMAX 81)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
