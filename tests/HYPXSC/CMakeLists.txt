#------------------------------------------------------------------------------
# Tests for hypxsc subroutine (hyperfine cross sections calculations)
#------------------------------------------------------------------------------


set(TEST_POT_SRC_FILE "pot_hypxsc.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_INPUT_FILES "")
set(TEST_KMAX 500)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


set(supported_cases atom-molecule-1spin atom-molecule-2spins molecule-molecule-1spin)

foreach(case IN LISTS supported_cases)
  set(TEST_ID hypxsc_${case})
  set(TEST_INPUT_FILES  "Hypxsc_${case}1.smt")
  set(TEST_OUTPUT_FILES "Hypxsc_${case}1.hfx")
  set(TEST_COMMAND_FILE "hypxsc_${case}.com")
  add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
endforeach()