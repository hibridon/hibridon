#------------------------------------------------------------------------------
# OH--H2 j2=1 user-defined basis
#------------------------------------------------------------------------------

set(TEST_ID ohh2_bausr_j1)
set(TEST_POT_SRC_FILE "pot_ohh2_bausr.F")
set(TEST_POT_DATA_FILES "pot_ohh2_bausr_b pot_ohh2_bausr_f")
set(TEST_COMMAND_FILE "ohh2_bausr_j1.com")
set(TEST_INPUT_FILES "Ohh2.inp")
set(TEST_OUTPUT_FILES "Job1.ics Job1.pcs Job1.xxsc")
set(TEST_KMAX 500)
set(TEST_T_MATRIX_SIZE kmax)

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_KMAX}" "${TEST_T_MATRIX_SIZE}")

#------------------------------------------------------------------------------
# OH(X)--He 3D PES using OH--H2 basis routine
#------------------------------------------------------------------------------

set(TEST_ID ohh2_bausr_heohx)
set(TEST_POT_SRC_FILE "pot_ohh2_bausr.F")
set(TEST_POT_DATA_FILES "pot_heohx_3d_b pot_heohx_3d_f")
set(TEST_COMMAND_FILE "ohh2_bausr_heohx.com")
set(TEST_INPUT_FILES "Testpisg.inp")
set(TEST_OUTPUT_FILES "Heohx1.ics Heohx1.pcs Heohx1.xxsc")
set(TEST_KMAX 150)
set(TEST_T_MATRIX_SIZE kmax)

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_KMAX}" "${TEST_T_MATRIX_SIZE}")