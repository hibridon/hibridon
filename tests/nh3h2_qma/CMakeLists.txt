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
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

set(TEST_ID nh3h2_qma_long)
set(TEST_POT_SRC_FILE "pot_nh3h2_qma.F90")
set(TEST_POT_DATA_FILES "pot_nh3h2_2009_fitvij_bf_62")
set(TEST_COMMAND_FILE "nh3h2_long_test.com")
set(TEST_INPUT_FILES "Nh3h2_po_long.inp")
set(TEST_OUTPUT_FILES "Nh3h2_long1.ics")
set(TEST_KMAX 5000)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "benchmark")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
