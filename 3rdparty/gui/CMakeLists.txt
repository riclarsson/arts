cmake_minimum_required(VERSION 3.6)

if (ARTS_GUI)
  project(imgui_interface)
  
  add_subdirectory(imgui)
  add_subdirectory(implot)
  add_subdirectory(filebrowser)
  
  include_directories(${PROJECT_SOURCE_DIR})
  add_executable(imgui_example imgui_example.cpp)
  target_link_libraries(imgui_example imgui)

  # glfw
  find_package(glfw3 REQUIRED)
  include_directories(${GLFW_INCLUDE_DIRS})
  link_libraries(${GLFW_LIBRARY_DIRS})
  
  # opengl
  set(OpenGL_GL_PREFERENCE GLVND)
  find_package(OpenGL REQUIRED)
  include_directories(${OPENGL_INCLUDE_DIRS})
  
  # glew
  find_package(GLEW REQUIRED)
  include_directories(${GLEW_INCLUDE_DIRS})
  
  target_link_libraries(
        imgui_example
        glfw
        ${OPENGL_LIBRARIES}
        ${GLEW_LIBRARIES}
        ${EXTRA_LIBS}
  )
  
  # Add several errors here
  
  message(STATUS "Added GUI libraries (-DARTS_GUI=0 to deactivate)")
else()
  message(STATUS "Added no GUI libraries (-DARTS_GUI=1 to activate)")
endif()

