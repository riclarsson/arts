add_library(quantum STATIC
            hitran_species.cc
            quantum.cc
            quantum_term_symbol.cc
            hitran_species_map.h
            hitran_species_info.cc
            )
target_link_libraries(quantum PUBLIC species arts_enum_options matpack binio)
target_include_directories(quantum PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(quantum PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

add_executable(make_molparam_map make_molparam_map.cc hitran_species_info.cc)
target_link_libraries(make_molparam_map PUBLIC species)
target_include_directories(make_molparam_map PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_custom_command(
  OUTPUT hitran_species_map.h
  COMMAND make_molparam_map ${ARTS_CAT_DATA_DIR}
  DEPENDS make_molparam_map
)
