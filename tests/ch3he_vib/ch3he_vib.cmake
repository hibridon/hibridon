#------------------------------------------------------------------------------
# CH3(X)--He Vibrational Relaxation
# 4D CCSD(T) PES, user defined pot/basis routine (ch3he_vib) by Q. Ma
#------------------------------------------------------------------------------

set(TEST_ID ch3he_vib)
set(TEST_POT_SRC_FILE "pot_ch3he_vib.F")
set(TEST_POT_DATA_FILES "pot_ch3he_vib_data pot_ch3he_vib_ylmasym pot_ch3he_vib_ylmsym")
set(TEST_COMMAND_FILE "ch3he_vib_test.com")
set(TEST_INPUT_FILES "Vises.inp")
set(TEST_OUTPUT_FILES "Job1.ics Job1.xxsc")
set(TEST_KMAX 1000)
set(TEST_T_MATRIX_SIZE kmax)

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_KMAX}" "${TEST_T_MATRIX_SIZE}")