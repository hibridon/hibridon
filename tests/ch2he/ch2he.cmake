#------------------------------------------------------------------------------
# CH2(a)--He
#------------------------------------------------------------------------------

set(TEST_ID ch2he_a)
set(TEST_POT_SRC_FILE "pot_ch2he_a190_c20.F")
set(TEST_POT_DATA_FILES "ch2he_a190.dat")
set(TEST_COMMAND_FILE "ch2ahe_stmixtst.com")
set(TEST_INPUT_FILES "Ch2he_para.inp")
set(TEST_OUTPUT_FILES "Ch2_p1.ics Ch2_p1.xsc")
set(TEST_KMAX 1451)
set(TEST_T_MATRIX_SIZE kmax)

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_KMAX}" "${TEST_T_MATRIX_SIZE}")

#------------------------------------------------------------------------------
# CH2(X)--He and spin-resolved cross sections
#------------------------------------------------------------------------------

set(TEST_ID ch2he_x)
set(TEST_POT_SRC_FILE "pot_ch2he_x52_c20_v3.F")
set(TEST_POT_DATA_FILES "ch2he_x52_v3.dat")
set(TEST_COMMAND_FILE "ch2xhe_stmixtst.com")
set(TEST_INPUT_FILES "Ch2x3he_para.inp")
set(TEST_OUTPUT_FILES "Ch2x_p1.xsc Ch2x_p1.ics Ch2x_p1.hfx")
#set(TEST_OUTPUT_FILES "Ch2x_p1.xsc Ch2x_p1.ics Ch2x_p1.hfx ch2he_x.stdout") # stdout files can't be automatically compared at the moment
# todo: make sure that this test works : not sure it actually works because ch2he_x.stdout shows a message about a missing Ch2_p1.smt file... ? And ch2he_x.stdout is missing the end part

set(TEST_KMAX 1451)
set(TEST_T_MATRIX_SIZE kmax)

add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_KMAX}" "${TEST_T_MATRIX_SIZE}")

#------------------------------------------------------------------------------
# CH2(X,a) singlet-triplet mixing
#------------------------------------------------------------------------------

# disabled because ch2he_stmix.stdout simply says something about a missing Ch2_p1.smt file... ?

# set(TEST_ID ch2he_stmix)
# set(TEST_POT_SRC_FILE "pot_ch2he_x52_c20_v3.F")
# set(TEST_POT_DATA_FILES "ch2he_x52_v3.dat")
# set(TEST_COMMAND_FILE "ch2xhe_stmixtst2.com")
# set(TEST_INPUT_FILES "Ch2x3he_para.inp")
# set(TEST_OUTPUT_FILES "ch2he_stmix.stdout")
# set(TEST_KMAX 1451)
# set(TEST_T_MATRIX_SIZE kmax)

# add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_KMAX}" "${TEST_T_MATRIX_SIZE}")