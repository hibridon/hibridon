#------------------------------------------------------------------------------
# vfit test
#------------------------------------------------------------------------------

set(TEST_ID vfit)
set(TEST_POT_SRC_FILE "pot_vfit.F90")
set(TEST_POT_DATA_FILES "FOLLMEG.BIN")
set(TEST_COMMAND_FILE "vfit_test.com")
set(TEST_INPUT_FILES "N2phetest.inp")
set(TEST_OUTPUT_FILES "Vfit_test1.ics")
set(TEST_KMAX 36)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

#------------------------------------------------------------------------------
# vfit_tenxsc test
# tests tenxsc with all allowed values of iframe.
#------------------------------------------------------------------------------

set(TEST_ID vfit_tenxsc)
set(TEST_POT_SRC_FILE "pot_vfit.F90")
set(TEST_POT_DATA_FILES "FOLLMEG.BIN")
set(TEST_COMMAND_FILE "vfit_tenxsc.com")
set(TEST_INPUT_FILES "N2phetest.inp")
set(TEST_OUTPUT_FILES "")  # todo : add some output files to check, because at the moment the check only tests if it crashes or not
set(TEST_KMAX 36)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
