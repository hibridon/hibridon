#------------------------------------------------------------------------------
# He--OH(X) RCCSD(T) Potential (Klos and Tobola)
#------------------------------------------------------------------------------

set(TEST_ID heohx)
set(TEST_POT_SRC_FILE "pot_heohx_lmax10.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "heohx_test.com")
set(TEST_INPUT_FILES "Heohx_test.inp")
set(TEST_OUTPUT_FILES "Heohx1.dcs Heohx1.xsc")
set(TEST_KMAX 151)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")


#------------------------------------------------------------------------------
# He--OH(X) RCCSD(T) Potential (Klos and Tobola) (quick version for code coverage)
#------------------------------------------------------------------------------

set(TEST_ID heohx_quick)
set(TEST_POT_SRC_FILE "pot_heohx_lmax10.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "heohx_quick.com")
set(TEST_INPUT_FILES "Heohx_quick.inp")
set(TEST_OUTPUT_FILES "Heohxq1.dcs Heohxq1.xsc")
set(TEST_KMAX 151)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
