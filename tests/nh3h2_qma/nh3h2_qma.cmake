#------------------------------------------------------------------------------
# NH3--H2 2009 Potential, custom basis routine by Q. Ma
#------------------------------------------------------------------------------

set(TEST_ID nh3h2_qma)
set(TEST_POT_SRC_FILE "pot_nh3h2_qma.F90")
set(TEST_POT_DATA_FILES "pot_nh3h2_2009_fitvij_bf_62")
set(TEST_COMMAND_FILE "nh3h2_test.com")
set(TEST_INPUT_FILES "Nh3h2_po.inp")
set(TEST_OUTPUT_FILES "Job1.ics")
set(TEST_KMAX 390)
set(TEST_T_MATRIX_SIZE kmax)

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_KMAX}" "${TEST_T_MATRIX_SIZE}")