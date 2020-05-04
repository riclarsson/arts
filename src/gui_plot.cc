

// Our includes
#include "gui_plot.h"

static void glfw_error_callback(int error, const char* description)
{
    fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

int main(int, char**)
{
    // Setup window
    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit())
        return 1;

    // Decide GL+GLSL versions
#if __APPLE__
    // GL 3.2 + GLSL 150
    const char* glsl_version = "#version 150";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // Required on Mac
#else
    // GL 3.0 + GLSL 130
    const char* glsl_version = "#version 130";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
    //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // 3.0+ only
#endif

    // Create window with graphics context
    GLFWwindow* window = glfwCreateWindow(1280, 720, "ARTS GUI", NULL, NULL);
    
    int* width = new int(0);
    int* height = new int(0);
    glfwGetWindowSize(window, width, height);
    
    if (window == NULL)
        return 1;
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // Enable vsync

    // Initialize OpenGL loader
#if defined(IMGUI_IMPL_OPENGL_LOADER_GL3W)
    bool err = gl3wInit() != 0;
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLEW)
    bool err = glewInit() != GLEW_OK;
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLAD)
    bool err = gladLoadGL() == 0;
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLBINDING2)
    bool err = false;
    glbinding::Binding::initialize();
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLBINDING3)
    bool err = false;
    glbinding::initialize([](const char* name) { return (glbinding::ProcAddress)glfwGetProcAddress(name); });
#else
    bool err = false; // If you use IMGUI_IMPL_OPENGL_LOADER_CUSTOM, your loader is likely to requires some form of initialization.
#endif
    if (err)
    {
        fprintf(stderr, "Failed to initialize OpenGL loader!\n");
        return 1;
    }

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();

    // Setup Platform/Renderer bindings
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);

    // Load Fonts
    // - If no fonts are loaded, dear imgui will use the default font. You can also load multiple fonts and use ImGui::PushFont()/PopFont() to select them.
    // - AddFontFromFileTTF() will return the ImFont* so you can store it if you need to select the font among multiple.
    // - If the file cannot be loaded, the function will return NULL. Please handle those errors in your application (e.g. use an assertion, or display an error and quit).
    // - The fonts will be rasterized at a given size (w/ oversampling) and stored into a texture when calling ImFontAtlas::Build()/GetTexDataAsXXXX(), which ImGui_ImplXXXX_NewFrame below will call.
    // - Read 'docs/FONTS.txt' for more instructions and details.
    // - Remember that in C/C++ if you want to include a backslash \ in a string literal you need to write a double backslash \\ !
    //io.Fonts->AddFontDefault();
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Roboto-Medium.ttf", 16.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Cousine-Regular.ttf", 15.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/DroidSans.ttf", 16.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/ProggyTiny.ttf", 10.0f);
    //ImFont* font = io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\ArialUni.ttf", 18.0f, NULL, io.Fonts->GetGlyphRangesJapanese());
    //IM_ASSERT(font != NULL);

    // Our states
    bool show_about=false;
    bool show_metrics=false;
    bool show_style=false;
    bool fullscreen=false;
    ImPlotRange copy_range;
    int stored_window_height, stored_window_width, stored_window_xpos, stored_window_ypos;
    ImVec4 clear_color = ImVec4(0.f, 0.f, 0.f, 0.f);
    
    // Test data
    auto x = PlotData("X-axis", 50000);
    auto y = PlotData("Y-axis", 50000);
    Vector tmp(50000);
    for (int i=0; i<50000; i++)
      tmp[i] = (6.28*5*i)/100000;
    x.set(tmp);
    
    for (int i=0; i<50000; i++)
      tmp[i] = std::sin(tmp[i]);
    y.set(tmp);            
    
    Line line("XY-plot", &x, &y);
    
    // Main loop
    while (!glfwWindowShouldClose(window))
    {
        // Poll and handle events (inputs, window resize, etc.)
        // You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui wants to use your inputs.
        // - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main application.
        // - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main application.
        // Generally you may always pass all inputs to dear imgui, and hide them from your application based on those two flags.
        glfwPollEvents();

        // Start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        
        // Menu bar
        if (ImGui::BeginMainMenuBar()) {
          
          // Main menu "Files"
          if (ImGui::BeginMenu("File")) {
            
            if (ImGui::MenuItem("Fullscreen", "F11")) {
              if (not fullscreen) {
                glfwGetWindowSize(window, &stored_window_width, &stored_window_height);
                glfwGetWindowPos(window, &stored_window_xpos, &stored_window_ypos);
                const GLFWvidmode * mode = glfwGetVideoMode(get_current_monitor(window));
                glfwSetWindowMonitor(window, get_current_monitor(window), 0, 0, mode->width, mode->height, 0);
              }
              else
                glfwSetWindowMonitor(window, NULL, stored_window_xpos, stored_window_ypos, stored_window_width, stored_window_height, 0);
              fullscreen = not fullscreen;
            }
            
            ImGui::Separator();
            
            if (ImGui::MenuItem("Quit", "Ctrl+X")) {
              glfwSetWindowShouldClose(window, 1);
            }
            
            // Must end the menu
            ImGui::EndMenu();
          }
          
          // Main menu "Tabs"
          if (ImGui::BeginMenu("Tabs")) {
            
            // Must end the menu
            ImGui::EndMenu();
          }
          
          // Main menu "Help"
          if (ImGui::BeginMenu("Help")) {
            
            ImGui::Checkbox("  About ImGui", &show_about);
            ImGui::Checkbox("  Metrics ImGui", &show_metrics);
            ImGui::Checkbox("  Style ImGui", &show_style);
            
            // Must end the menu
            ImGui::EndMenu();
          }
          
          // Must end the menubar
          ImGui::EndMainMenuBar();
        }
        
        // Menu shortcuts
        if (io.KeyCtrl and ImGui::IsKeyPressed(GLFW_KEY_X)) {
          glfwSetWindowShouldClose(window, 1);
        } else if (ImGui::IsKeyPressed(GLFW_KEY_F11) or (fullscreen and ImGui::IsKeyPressed(GLFW_KEY_ESCAPE))) {
          if (not fullscreen) {
            glfwGetWindowSize(window, &stored_window_width, &stored_window_height);
            glfwGetWindowPos(window, &stored_window_xpos, &stored_window_ypos);
            const GLFWvidmode * mode = glfwGetVideoMode(get_current_monitor(window));
            glfwSetWindowMonitor(window, get_current_monitor(window), 0, 0, mode->width, mode->height, 0);
          }
          else
            glfwSetWindowMonitor(window, NULL, stored_window_xpos, stored_window_ypos, stored_window_width, stored_window_height, 0);
          fullscreen = not fullscreen;
        }
        
        // Plotting defaults:
        ImGui::GetPlotStyle().LineWeight = 4;
        
        //Cursors and sizes
        glfwGetWindowSize(window, width, height);
        ImVec2 pos = ImGui::GetCursorPos();
        ImVec2 size = {float(*width)-2*pos.x, float(*height)-pos.y};
        
        // Show a simple window
        ImGui::SetNextWindowPos(ImGui::GetCursorPos());
        ImGui::SetNextWindowSize(size);
        if (ImGui::Begin("Plot tool 4", NULL, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_MenuBar)) {
          ImGui::SetNextPlotRange(1, 10, 0, 6);
          if (ImGui::BeginPlot("XY-plot", "X-values", "Y-values", {-1, -1}, ImPlotFlags_MousePos | ImPlotFlags_Legend | ImPlotFlags_Highlight | ImPlotFlags_Selection | ImPlotFlags_ContextMenu)) {
            ImGui::Plot(line.name(), line.getter(), (void*)&line, line.size());
            ImGui::EndPlot();
          }
          
          // Menu bar
          if (ImGui::BeginMenuBar()) {
            
            if (ImGui::BeginMenu("Range")) {
              ImPlotRange range = ImGui::GetPlotRange();
              bool change = false;
              change |= ImGui::InputFloat2("X", &range.XMin, "%g", ImGuiInputTextFlags_CharsScientific);
              change |= ImGui::InputFloat2("Y", &range.YMin, "%g", ImGuiInputTextFlags_CharsScientific);
              
              ImGui::Separator();
              
              if (ImGui::MenuItem("Copy")) {
                copy_range = range;
              }
              
              bool enable = not (std::isnan(copy_range.XMax) or std::isnan(copy_range.XMin) or std::isnan(copy_range.YMax) or std::isnan(copy_range.YMin));
              if (ImGui::MenuItem("Paste", NULL, false, enable)) {
                range = copy_range;
                change = true;
              }
              
              if (change) ImGui::SetPlotRange(range);
              
              // Must end the menu
              ImGui::EndMenu();
            }
            
            if (ImGui::BeginMenu("Scale")) {
              int N = int(line.TypeModifier());
              if (ImGui::SliderInt("Linear Average", &N, 1, int(line.maxsize()))) {
                line.Linear(N);
              }
              
              // Must end the menu
              ImGui::EndMenu();
            }
            
            // Must end the menubar
            ImGui::EndMenuBar();
          }
        }
        ImGui::End();
        
        // Show special windows
        if (show_about) {
          ImGui::ShowAboutWindow();
        }
        if (show_metrics) {
          ImGui::ShowMetricsWindow();
        }
        if (show_style) {
          ImGui::ShowStyleEditor();
        }

        // Rendering
        ImGui::Render();
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
    }

    // Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
