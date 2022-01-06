# CO--H Transport Cross Section

set(TEST_ID hco_trnprt)
set(TEST_POT_SRC_FILE "pot_hco.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "Hco_trnprt.com")
set(TEST_INPUT_FILES "Hco_test.inp")
set(TEST_OUTPUT_FILES "Hcotrn1.ics Hcotrn1.trn Hcotrn1.xsc")
set(TEST_KMAX 301)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

# CO--H Transport Cross Section (QUICK)

set(TEST_ID hco_trnprt_quick)
set(TEST_POT_SRC_FILE "pot_hco.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "Hco_trnprt_quick.com")
set(TEST_INPUT_FILES "Hco_testq.inp")
set(TEST_OUTPUT_FILES "Hcotrnq1.ics Hcotrnq1.trn Hcotrnq1.xsc")
set(TEST_KMAX 301)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
