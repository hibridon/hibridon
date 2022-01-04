# CMake project file for hibridon

# cmake version >= 3.3 is required for IN_LIST
cmake_minimum_required (VERSION 3.3)
project (hibridon VERSION 5.0)
enable_language (Fortran)

find_package( BLAS REQUIRED )
find_package( LAPACK REQUIRED )


add_library(hib
     src/assert.f90
     src/hiamp.F90
     src/hibound.F90
     src/hibrid1.F90
     src/hibrid2.F90
     src/hibrid3.F90
     src/hibrid4.F90
     src/hibrid5.F90
     src/hicommon.F90
     src/hidifcrs.F90
     src/hiflow.F90
     src/hihypxsc.F90
     src/hiiolib_f.F90
     src/himatrix.F90
     src/himem.F90
     src/hinput.F90
     src/hiprsbr.F90
     src/histmix.F90
     src/hisystem.F90
     src/hitensor.F90
     src/hitrnprt.F90
     src/hiutil.F90
     src/hivector.F90
     src/hipotutil.F90
     src/hiversion.F90
     src/bases/hiba01_1sg.F90
     src/bases/hiba02_2sg.F90
     src/bases/hiba03_2pi.F90
     src/bases/hiba04_sgpi.F90
     src/bases/hiba05_pi.F90
     src/bases/hiba06_stp.F90
     src/bases/hiba07_13p.F90
     src/bases/hiba08_2mol.F90
     src/bases/hiba09_stpln.F90
     src/bases/hiba10_22p.F90
     src/bases/hiba11_1del.F90
     src/bases/hiba12_h2p.F90
     src/bases/hiba13_h3p.F90
     src/bases/hiba14_2del.F90
     src/bases/hiba15_diat2p.F90
     src/bases/hiba16_astp.F90
     src/bases/hiba17_ch2x.F90
     src/bases/hiba18_stp1.F90
     src/bases/hiba19_sgpi1.F90
     src/bases/hiba20_2pi1sg.F90
     src/bases/hiba21_stp1sg.F90
     src/bases/hiba22_1d3p.F90
     src/bases/hiba23_3p2s.F90
     src/bases/hiba24_sphtp.F90
     src/bases/hiba25_1sg1sg.F90
     src/bases/hiba26_2sg1sg.F90
     src/bases/hiba27_astp1.F90
     src/bases/hiba28_3sg1sg.F90
     src/bases/hiba29_astp2.F90
     src/bases/hiba30_astp3.F90)
