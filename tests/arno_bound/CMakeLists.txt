#------------------------------------------------------------------------------
# NO--Ar bound states
# See "Setting the parameters for a bound-state calculation" in Hibridon document
# See also J. Phys. Chem. A 112, 9483 (2008)
#------------------------------------------------------------------------------

set(TEST_ID arno_bound)
set(TEST_POT_SRC_FILE "pot_arno_ccsdt.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "arno_bound.com")
set(TEST_INPUT_FILES "Arno.inp")
set(TEST_OUTPUT_FILES "Job.evl Jmax2.evl Jmax3.evl Jmax8.evl Jmax10.evl")
set(TEST_KMAX 42)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
