cmake_minimum_required(VERSION 3.14)
project(caviar2 LANGUAGES CXX)

# 1. Build the library
add_library(caviar2
    src/caviar2.cpp
)

target_include_directories(caviar2 PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
    $<INSTALL_INTERFACE:include>
)

# 2. Add test executables from tests/ directory
# Assuming each main.cpp is a separate test executable

file(GLOB TEST_SOURCES tests/*.cpp)

foreach(test_src ${TEST_SOURCES})
    get_filename_component(test_name ${test_src} NAME_WE)  # extract filename without extension
    
    add_executable(${test_name} ${test_src})
    target_link_libraries(${test_name} PRIVATE caviar2)
    target_include_directories(${test_name} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/inc
    )
endforeach()

# Optionally enable testing with CTest
enable_testing()

foreach(test_src ${TEST_SOURCES})
    get_filename_component(test_name ${test_src} NAME_WE)
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()
