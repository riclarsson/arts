add_executable(test_lbl_perf test_lbl_perf.cpp)

target_link_libraries(test_lbl_perf PUBLIC lbl rng)

add_executable(test_lbl_ecs test_lbl_ecs.cpp)
target_link_libraries(test_lbl_ecs PRIVATE lbl)
add_test(NAME cpp.fast.lbl.ecs COMMAND test_lbl_ecs)
add_dependencies(check-deps test_lbl_ecs)
