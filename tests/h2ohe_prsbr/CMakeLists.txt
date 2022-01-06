#------------------------------------------------------------------------------
# H2O--He Pressure Broadening Cross Section
#------------------------------------------------------------------------------

set(TEST_ID h2ohe_prsbr)
set(TEST_POT_SRC_FILE "pot_h2ohe.F90")
set(TEST_POT_DATA_FILES "h2o_coefd.dat h2o_params.dat h2o_coefi.dat")
set(TEST_COMMAND_FILE "prsbr.com")
set(TEST_INPUT_FILES "H2ohe_ortho.inp H2ohe_para.inp")
set(TEST_OUTPUT_FILES "H2o11.ppb")
set(TEST_KMAX 200)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")


#------------------------------------------------------------------------------
# H2O--He Pressure Broadening Cross Section (Quick version)
#------------------------------------------------------------------------------

set(TEST_ID h2ohe_prsbr_quick)
set(TEST_POT_SRC_FILE "pot_h2ohe.F90")
set(TEST_POT_DATA_FILES "h2o_coefd.dat h2o_params.dat h2o_coefi.dat")
set(TEST_COMMAND_FILE "prsbr_quick.com")
set(TEST_INPUT_FILES "H2ohe_orthoq.inp H2ohe_paraq.inp")
set(TEST_OUTPUT_FILES "H2o1q1.ppb")
set(TEST_KMAX 200)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
