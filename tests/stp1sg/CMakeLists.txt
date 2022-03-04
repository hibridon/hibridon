#------------------------------------------------------------------------------
# CH3(X)--He Dagdigian CCSD(T) PES with SYMTOP--1Sigma basis
#------------------------------------------------------------------------------

set(TEST_ID stp1sg_ch3he)
set(TEST_POT_SRC_FILE "pot_stp1sg_qma.F90")
set(TEST_POT_DATA_FILES "pot_ch3he_rigid_stp1sg.txt")
set(TEST_COMMAND_FILE "stp1sg_och3he.com")
set(TEST_INPUT_FILES "Stp_och3he.inp")
set(TEST_OUTPUT_FILES "Stp_och3he1.dcs Stp_och3he1.ics Stp_och3he1.pcs Stp_och3he1.xxsc")
set(TEST_KMAX 1000)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

#------------------------------------------------------------------------------
# CH3(X)--He Dagdigian CCSD(T) PES with SYMTOP--1Sigma basis (quick version for code coverage)
#------------------------------------------------------------------------------

set(TEST_ID stp1sg_ch3he_quick)
set(TEST_POT_SRC_FILE "pot_stp1sg_qma.F90")
set(TEST_POT_DATA_FILES "pot_ch3he_rigid_stp1sg.txt")
set(TEST_COMMAND_FILE "stp1sg_ch3he_quick.com")
set(TEST_INPUT_FILES "Stp_ch3he_quick.inp")
set(TEST_OUTPUT_FILES "Ch3heq1.dcs Ch3heq1.ics Ch3heq1.pcs Ch3heq1.xxsc")
set(TEST_KMAX 1000)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")

#------------------------------------------------------------------------------
# p-NH3--p-H2 test with Symtop-1sigma basis
#------------------------------------------------------------------------------

set(TEST_ID stp1sg_nh3h2)
set(TEST_POT_SRC_FILE "pot_stp1sg_qma.F90")
set(TEST_POT_DATA_FILES "pot_nh3h2_cp.dat")
set(TEST_COMMAND_FILE "stp1sg_pnh3ph2.com")
set(TEST_INPUT_FILES "Stp_pnh3ph2.inp")
set(TEST_OUTPUT_FILES "Stp_pnh3ph21.xxsc Stp_pnh3ph21.ics")
set(TEST_KMAX 1000)
set(TEST_T_MATRIX_SIZE kmax)
set(TEST_LABELS "coverage")


add_hibridon_test("${TEST_ID}" "${TEST_POT_SRC_FILE}" "${TEST_POT_DATA_FILES}" "${TEST_COMMAND_FILE}" "${TEST_INPUT_FILES}" "${TEST_OUTPUT_FILES}" "${TEST_T_MATRIX_SIZE}" "${TEST_LABELS}")
