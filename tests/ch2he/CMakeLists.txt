#------------------------------------------------------------------------------
# CH2(a)--He
#------------------------------------------------------------------------------

set(TEST_ID ch2he_a)
set(TEST_POT_SRC_FILE "pot_ch2he_a190_c20.F90")
set(TEST_POT_DATA_FILES "ch2he_a190.dat")
set(TEST_COMMAND_FILE "ch2ahe_stmixtst.com")
set(TEST_INPUT_FILES "Ch2he_para.inp")
set(TEST_OUTPUT_FILES "Ch2_p1.ics Ch2_p1.xsc")
set(TEST_KMAX 1451)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

#------------------------------------------------------------------------------
# CH2(a)--He (quick version for test coverage)
#------------------------------------------------------------------------------

set(TEST_ID ch2he_a_quick)
set(TEST_POT_SRC_FILE "pot_ch2he_a190_c20.F90")
set(TEST_POT_DATA_FILES "ch2he_a190.dat")
set(TEST_COMMAND_FILE "ch2ahe_stmixtst_quick.com")
set(TEST_INPUT_FILES "Ch2he_para.inp")
set(TEST_OUTPUT_FILES "Ch2_pq1.ics Ch2_pq1.xsc")
set(TEST_KMAX 1451)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

#------------------------------------------------------------------------------
# CH2(X)--He and spin-resolved cross sections
#------------------------------------------------------------------------------

set(TEST_ID ch2he_x)
set(TEST_POT_SRC_FILE "pot_ch2he_x52_c20_v3.F90")
set(TEST_POT_DATA_FILES "ch2he_x52_v3.dat")
set(TEST_COMMAND_FILE "ch2xhe_stmixtst.com")
set(TEST_INPUT_FILES "Ch2x3he_para.inp Ch2_p1.smt") # Ch2_p1.smt is an output of ch2he_a test
set(TEST_OUTPUT_FILES "Ch2x_p1.xsc Ch2x_p1.ics Ch2x_p1.hfx")
set(TEST_KMAX 1451)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

#------------------------------------------------------------------------------
# CH2(X)--He and spin-resolved cross sections (quick version for test coverage)
#------------------------------------------------------------------------------

set(TEST_ID ch2he_x_quick)
set(TEST_POT_SRC_FILE "pot_ch2he_x52_c20_v3.F90")
set(TEST_POT_DATA_FILES "ch2he_x52_v3.dat")
set(TEST_COMMAND_FILE "ch2xhe_stmixtst_quick.com")
set(TEST_INPUT_FILES "Ch2x3he_para.inp Ch2_pq1.smt") # Ch2_p1q.smt is an output of ch2he_a quick test
set(TEST_OUTPUT_FILES "Ch2x_pq1.xsc Ch2x_pq1.ics Ch2x_pq1.hfx")
set(TEST_KMAX 1451)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

#------------------------------------------------------------------------------
# CH2(X,a) singlet-triplet mixing
#------------------------------------------------------------------------------

set(TEST_ID ch2he_stmix)
set(TEST_POT_SRC_FILE "pot_ch2he_x52_c20_v3.F90")
set(TEST_POT_DATA_FILES "ch2he_x52_v3.dat")
set(TEST_COMMAND_FILE "ch2xhe_stmixtst2.com")
set(TEST_INPUT_FILES "Ch2x3he_para.inp Ch2_p1.smt Ch2x_p1.smt") # Ch2x_p1.smt is an output of ch2he_x test
set(TEST_OUTPUT_FILES "ch2he_stmix.stdout")
set(TEST_KMAX 1451)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")

# Test disabled because stdout files are not a good way to test that the program runs well.
# We should modify the .com file so hibridon generates a real data output file 
#add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
