#------------------------------------------------------------------------------
#  B(2P)--ortho H2 bound states
#------------------------------------------------------------------------------

set(TEST_ID boh2_b)
set(TEST_POT_SRC_FILE "pot_boh2.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "boh2_b.com")
set(TEST_INPUT_FILES "Boh2_b.inp")
set(TEST_OUTPUT_FILES "Boh2_b.evl boh2_b.stdout")
set(TEST_KMAX 120)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")


#------------------------------------------------------------------------------
#  B(2P)--ortho H2 bound states with wavefl
#------------------------------------------------------------------------------

set(TEST_ID boh2_bw)
set(TEST_POT_SRC_FILE "pot_boh2.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "boh2_bw.com")
set(TEST_INPUT_FILES "Boh2_bw.inp")
set(TEST_OUTPUT_FILES "Boh2_bw.evl boh2_bw.stdout")
set(TEST_KMAX 120)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")


#------------------------------------------------------------------------------
# B(2P)--ortho H2 bound states (modified from the original hibtest)
#------------------------------------------------------------------------------

set(TEST_ID boh2_b2)
set(TEST_POT_SRC_FILE "pot_boh2.F90")
set(TEST_POT_DATA_FILES "")
set(TEST_COMMAND_FILE "boh2_b2.com")
set(TEST_INPUT_FILES "Boh2_b2.inp")
set(TEST_OUTPUT_FILES "Boh2_b2a.evl Boh2_b2b.evl boh2_b2.stdout")
set(TEST_KMAX 5)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
