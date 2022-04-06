# unit tests

file (GLOB unit_test_sources ${hibridon_SOURCE_DIR}/src/unittests/*.F90)
message ("unit_test_sources = ${unit_test_sources}")
#set (unit_tests_list "")
foreach(unit_test_source_file IN LISTS unit_test_sources)
   get_filename_component(unit_test_name ${unit_test_source_file} NAME_WE)
   message ("Adding ${unit_test_name} unit test")
   set(unit_test_exe_name "unit_${unit_test_name}")
   add_executable(${unit_test_exe_name} ${unit_test_source_file})
   target_compile_definitions(${unit_test_exe_name} PUBLIC T_MATRIX_SIZE=${TEST_T_MATRIX_SIZE})
   target_link_libraries(${unit_test_exe_name} hib ${HIBRIDON_LINK_LIBRARIES})
   target_include_directories(${unit_test_exe_name} PRIVATE ${hibridon_SOURCE_DIR}/lib) # to find the included common files
   target_include_directories(${unit_test_exe_name} PRIVATE ${hibridon_BINARY_DIR}/lib) # to find the .mod files
   target_include_directories(${unit_test_exe_name} PUBLIC ${hibridon_SOURCE_DIR}/src)
   target_compile_options(${unit_test_exe_name} PRIVATE ${HIBRIDON_COMPILE_OPTIONS})
   add_test(unit_${unit_test_name}_build
         "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${unit_test_exe_name})
   add_test(NAME unit_${unit_test_name} COMMAND ${unit_test_exe_name})
   SET_TESTS_PROPERTIES(unit_${unit_test_name} PROPERTIES DEPENDS unit_${unit_test_name}_build)
   #list(APPEND unit_tests_list unit_${unit_test_name})
endforeach()
