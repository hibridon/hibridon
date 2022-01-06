#------------------------------------------------------------------------------
# OH--H2 j2=1 user-defined basis
#------------------------------------------------------------------------------

set(TEST_ID ohh2_bausr_j1)
set(TEST_POT_SRC_FILE "pot_ohh2_bausr.F90")
set(TEST_POT_DATA_FILES "pot_ohh2_bausr_b pot_ohh2_bausr_f")
set(TEST_COMMAND_FILE "ohh2_bausr_j1.com")
set(TEST_INPUT_FILES "Ohh2.inp")
set(TEST_OUTPUT_FILES "Job1.ics Job1.pcs Job1.xxsc")
set(TEST_KMAX 500)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

#------------------------------------------------------------------------------
# OH--H2 j2=1 user-defined basis (quick version for code coverage)
#------------------------------------------------------------------------------

set(TEST_ID ohh2_bausr_j1_quick)
set(TEST_POT_SRC_FILE "pot_ohh2_bausr.F90")
set(TEST_POT_DATA_FILES "pot_ohh2_bausr_b pot_ohh2_bausr_f")
set(TEST_COMMAND_FILE "ohh2_bausr_j1_quick.com")
set(TEST_INPUT_FILES "Ohh2_quick.inp")
set(TEST_OUTPUT_FILES "Quick1.ics Quick1.pcs Quick1.xxsc")
set(TEST_KMAX 500)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

#------------------------------------------------------------------------------
# OH(X)--He 3D PES using OH--H2 basis routine
#------------------------------------------------------------------------------

set(TEST_ID ohh2_bausr_heohx)
set(TEST_POT_SRC_FILE "pot_ohh2_bausr.F90")
set(TEST_POT_DATA_FILES "pot_heohx_3d_b pot_heohx_3d_f")
set(TEST_COMMAND_FILE "ohh2_bausr_heohx.com")
set(TEST_INPUT_FILES "Testpisg.inp")
set(TEST_OUTPUT_FILES "Heohx1.ics Heohx1.pcs Heohx1.xxsc")
set(TEST_KMAX 150)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

#------------------------------------------------------------------------------
# OH(X)--He 3D PES using OH--H2 basis routine (quick version for code coverage)
#------------------------------------------------------------------------------

set(TEST_ID ohh2_bausr_heohx_quick)
set(TEST_POT_SRC_FILE "pot_ohh2_bausr.F90")
set(TEST_POT_DATA_FILES "pot_heohx_3d_b pot_heohx_3d_f")
set(TEST_COMMAND_FILE "ohh2_bausr_heohx_quick.com")
set(TEST_INPUT_FILES "Testpisg_quick.inp")
set(TEST_OUTPUT_FILES "Heohxq1.ics Heohxq1.pcs Heohxq1.xxsc")
set(TEST_KMAX 150)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
