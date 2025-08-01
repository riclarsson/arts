add_library(xml_partfun STATIC xml_io_partfun.cc)
target_link_libraries (xml_partfun PUBLIC matpack xml_io coretypes)
target_include_directories(xml_partfun PUBLIC ${ARTS_SOURCE_DIR}/src)

add_executable(partfun_generator make_auto_partfuns.cc)
target_link_libraries (partfun_generator PRIVATE xml_partfun species)

file(GLOB_RECURSE partfun_xml_files LIST_DIRECTORIES FALSE RELATIVE ${ARTS_PARTFUN_DATA_PATH} "*.xml")
set(partfun_cpp_files ${partfun_xml_files})
list(TRANSFORM partfun_cpp_files REPLACE ".xml" ".cpp")
list(TRANSFORM partfun_xml_files PREPEND "${ARTS_PARTFUN_DATA_PATH}/")

# Generate auto_partfun.h for the ARTS interface
add_custom_command(
    OUTPUT auto_partfun.h ${partfun_cpp_files}
    COMMAND partfun_generator ${partfun_xml_files}
    DEPENDS partfun_generator
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating automatic partition function interface"
)

add_library(auto_partfun STATIC ${partfun_cpp_files} auto_partfun.h)
target_link_libraries(auto_partfun PUBLIC matpack binio xml_partfun)
target_include_directories(auto_partfun PUBLIC ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(auto_partfun PUBLIC ${ARTS_SOURCE_DIR}/src)

if (IPO_SUPPORTED)
    set_property(TARGET auto_partfun PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

add_library(partfun STATIC partfun.cc)
target_link_libraries(partfun PUBLIC auto_partfun species)
