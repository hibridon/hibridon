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
# vfit_tenxsc_f0 test
# tests tenxsc with iframe = 0.
#------------------------------------------------------------------------------

set(TEST_ID vfit_tenxsc_f0)
set(TEST_POT_SRC_FILE "pot_vfit.F90")
set(TEST_POT_DATA_FILES "FOLLMEG.BIN")
set(TEST_COMMAND_FILE "vfit_tenxsc_f0.com")
set(TEST_INPUT_FILES "N2phetest.inp")
set(TEST_OUTPUT_FILES "Vfit_tenf01.tcb Vfit_tenf01.tcs")
set(TEST_KMAX 36)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

#------------------------------------------------------------------------------
# vfit_tenxsc_f1 test
# tests tenxsc with iframe = 1.
#------------------------------------------------------------------------------

set(TEST_ID vfit_tenxsc_f1)
set(TEST_POT_SRC_FILE "pot_vfit.F90")
set(TEST_POT_DATA_FILES "FOLLMEG.BIN")
set(TEST_COMMAND_FILE "vfit_tenxsc_f1.com")
set(TEST_INPUT_FILES "N2phetest.inp")
set(TEST_OUTPUT_FILES "Vfit_tenf11.tcb Vfit_tenf11.tcs")
set(TEST_KMAX 36)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

#------------------------------------------------------------------------------
# vfit_tenxsc_f2 test
# tests tenxsc with iframe = 2.
#------------------------------------------------------------------------------

set(TEST_ID vfit_tenxsc_f2)
set(TEST_POT_SRC_FILE "pot_vfit.F90")
set(TEST_POT_DATA_FILES "FOLLMEG.BIN")
set(TEST_COMMAND_FILE "vfit_tenxsc_f2.com")
set(TEST_INPUT_FILES "N2phetest.inp")
set(TEST_OUTPUT_FILES "Vfit_tenf21.dch")
set(TEST_KMAX 36)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

#------------------------------------------------------------------------------
# vfit_tenxsc_f3 test
# tests tenxsc with iframe = 3.
#------------------------------------------------------------------------------

set(TEST_ID vfit_tenxsc_f3)
set(TEST_POT_SRC_FILE "pot_vfit.F90")
set(TEST_POT_DATA_FILES "FOLLMEG.BIN")
set(TEST_COMMAND_FILE "vfit_tenxsc_f3.com")
set(TEST_INPUT_FILES "N2phetest.inp")
set(TEST_OUTPUT_FILES "Vfit_tenf31.dcga")
set(TEST_KMAX 36)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
