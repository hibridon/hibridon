#------------------------------------------------------------------------------
# He--OH(X) RCCSD(T) Potential (Klos and Tobola)
#------------------------------------------------------------------------------

set(TEST_ID heohx)
set(TEST_POT_SRC_FILE "pot_heohx_lmax10.F")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "heohx_test.com")
set(TEST_INPUT_FILES "Heohx_test.inp")
set(TEST_OUTPUT_FILES "Heohx1.dcs Heohx1.xsc")
set(TEST_KMAX 151)
set(TEST_T_MATRIX_SIZE kmax)

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_KMAX}" "${TEST_T_MATRIX_SIZE}")
