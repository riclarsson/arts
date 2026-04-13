add_library(sensor STATIC obsel.cpp sensor_meta_info.cpp)
target_link_libraries(sensor PUBLIC matpack rtepack arts_enum_options)
target_include_directories(sensor PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
