cmake_minimum_required(VERSION 3.6)
add_library(
        implot STATIC
            implot.cpp
            implot.h
        )
target_include_directories(implot PUBLIC
    ../imgui
    ./
)
target_compile_definitions(implot PUBLIC ImDrawIdx=unsigned\ int)
