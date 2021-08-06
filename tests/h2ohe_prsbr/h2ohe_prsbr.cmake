#------------------------------------------------------------------------------
# H2O--He Pressure Broadening Cross Section
#------------------------------------------------------------------------------

set(TEST_ID h2ohe_prsbr)
set(TEST_POT_SRC_FILE "pot_h2ohe.F")
set(TEST_POT_DATA_FILES "h2o_coefd.dat h2o_params.dat h2o_coefi.dat")
set(TEST_COMMAND_FILE "prsbr.com")
set(TEST_INPUT_FILES "H2ohe_ortho.inp H2ohe_para.inp")
set(TEST_OUTPUT_FILES "H2o11.ppb")
set(TEST_KMAX 200)
set(TEST_T_MATRIX_SIZE kmax)

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_KMAX}" "${TEST_T_MATRIX_SIZE}")
