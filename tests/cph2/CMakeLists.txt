# HCO--H transport cross sections
# title = C+(2P) + H2

set(TEST_ID cph2)
set(TEST_POT_SRC_FILE "pot_cph2.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "Cph2_test.com")
set(TEST_INPUT_FILES "Cph2_test.inp")
set(TEST_OUTPUT_FILES "Cph21.pcs Cph21.xsc")  #  
set(TEST_KMAX 200)
set(TEST_T_MATRIX_SIZE kmax)

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")