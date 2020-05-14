if (ARTS_GUI)
  cmake_minimum_required(VERSION 3.10)
  project(imgui_interface)

  add_subdirectory(imgui)
  add_subdirectory(implot)
  
  add_executable(imgui_example imgui_example.cpp)
  target_include_directories(imgui_example PUBLIC ${PROJECT_SOURCE_DIR})
  target_link_libraries(imgui_example imgui)

  # Add several errors here
  
  message(STATUS "Added GUI libraries (-DARTS_GUI=0 to deactivate)")
else()
  message(STATUS "Added no GUI libraries (-DARTS_GUI=1 to activate)")
endif()
