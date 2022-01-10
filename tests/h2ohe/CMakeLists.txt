#------------------------------------------------------------------------------
# H2O--He SAPT Potential
#------------------------------------------------------------------------------

set(TEST_ID h2ohe)
set(TEST_POT_SRC_FILE "pot_h2ohe.F90")
set(TEST_POT_DATA_FILES "h2o_coefd.dat h2o_params.dat h2o_coefi.dat")
set(TEST_COMMAND_FILE "h2ohe_test.com")
set(TEST_INPUT_FILES "H2ohe_test.inp")
set(TEST_OUTPUT_FILES "H2ohe1.pcs H2ohe1.ics")
set(TEST_KMAX 200)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

#------------------------------------------------------------------------------
# H2O--He transport cross section
#------------------------------------------------------------------------------

set(TEST_ID h2ohe_trnprt)
set(TEST_POT_SRC_FILE "pot_h2ohe.F90")
set(TEST_POT_DATA_FILES "h2o_coefd.dat h2o_params.dat h2o_coefi.dat")
set(TEST_COMMAND_FILE "h2ohe_trnprt_test.com")
set(TEST_INPUT_FILES "H2ohe_test.inp")
set(TEST_OUTPUT_FILES "H2ohetrn1.ics H2ohetrn1.pcs H2ohetrn1.trn")
set(TEST_KMAX 200)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
