add_library(nonlte STATIC nonlte_atm.cpp)

target_link_libraries(nonlte PUBLIC matpack atm lbl)
target_include_directories(nonlte PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
