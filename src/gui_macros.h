#ifndef GUI_MACROS_H
#define GUI_MACROS_H

/** Begins an endless GUI loop 
 * 
 * Sets:
 *  clear_color: Background color of loop
 * Expects:
 *  window: A working glfw window
 */
#define BeginWhileLoopARTSGUI                       \
  ImVec4 clear_color = ImVec4(0.f, 0.f, 0.f, 0.f);  \
  while (!glfwWindowShouldClose(window))            \
  {                                                 \
    glfwPollEvents();                               \
                                                    \
    ImGui_ImplOpenGL3_NewFrame();                   \
    ImGui_ImplGlfw_NewFrame();                      \
    ImGui::NewFrame() 

/** Ends an endless GUI loop 
 * 
 * Expects:
 *  clear_color: Background color of loop
 *  window: A working glfw window
 */
#define EndWhileLoopARTSGUI                                                   \
    ImGui::Render();                                                          \
    int display_w, display_h;                                                 \
    glfwGetFramebufferSize(window, &display_w, &display_h);                   \
    glViewport(0, 0, display_w, display_h);                                   \
    glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w); \
    glClear(GL_COLOR_BUFFER_BIT);                                             \
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());                   \
    glfwSwapBuffers(window);                                                  \
  }                                                                           \

/** Sets the glsl_version variable and hints for the windows
 * 
 * System dependent...
 * 
 * Sets:
 *  glsl_version: Version of GLSL
 */
#if __APPLE__
  #define Selectglsl_versionARTSGUI                                           \
    const char* glsl_version = "#version 150";                                \
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);                            \
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);                            \
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);            \
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#else
  #define Selectglsl_versionARTSGUI                                           \
    const char* glsl_version = "#version 130";                                \
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);                            \
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
#endif

/** Initialize the ARTSGUI
 * 
 * Should set up all things needed to begin the main loop of drawing things
 * 
 * Sets:
 *  window: A glfw window
 *  glsl_version: Version of GLSL
 */
#define InitializeARTSGUI                                                   \
  glfwSetErrorCallback(glfw_error_callback);                                \
  if (!glfwInit())                                                          \
    throw std::runtime_error("Cannot initialize a window");                 \
                                                                            \
  Selectglsl_versionARTSGUI                                                 \
                                                                            \
  GLFWwindow* window = glfwCreateWindow(1280, 720, "ARTS GUI", NULL, NULL); \
  if (window == NULL)                                                       \
    throw std::runtime_error("Cannot create a window");                     \
  glfwMakeContextCurrent(window);                                           \
  glfwSwapInterval(1);                                                      \
                                                                            \
  if (glewInit() != GLEW_OK)                                                \
    throw std::runtime_error("Cannot initialize OpenGL loader");            \
                                                                            \
  IMGUI_CHECKVERSION();                                                     \
  ImGui::CreateContext();                                                   \
  ImGui::StyleColorsDark();                                                 \
                                                                            \
  ImGui_ImplGlfw_InitForOpenGL(window, true);                               \
  ImGui_ImplOpenGL3_Init(glsl_version)

/** Cleanup for the ARTSGUI
 * 
 * Should clean up all things InitializeARTSGUI
 * 
 * Destroys:
 *  window: The glfw window
 */
#define CleanupARTSGUI          \
  ImGui_ImplOpenGL3_Shutdown(); \
  ImGui_ImplGlfw_Shutdown();    \
  ImGui::DestroyContext();      \
                                \
  glfwDestroyWindow(window);    \
  glfwTerminate()

#endif  // GUI_MACROS_H
