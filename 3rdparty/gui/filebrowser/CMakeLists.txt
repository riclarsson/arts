cmake_minimum_required(VERSION 3.6)
add_library(
        imfilebrowser STATIC
        ImGuiFileBrowser.cpp
        ImGuiFileBrowser.h
)
target_include_directories(imfilebrowser PUBLIC
    ../imgui
    ./
)

