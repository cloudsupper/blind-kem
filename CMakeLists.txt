cmake_minimum_required(VERSION 3.10)

# Enable compile commands export for IDE support
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(rlwe-signature 
        VERSION 0.1.0
        DESCRIPTION "RLWE Based KEM Scheme Implementation"
        LANGUAGES CXX)

# Require C++17 for better type support and features
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add compiler warnings
if(MSVC)
    add_compile_options(/W4 /WX)
else()
    add_compile_options(-Wall -Wextra -Werror -pedantic)
endif()

# Option for building tests
option(BUILD_TESTS "Build test suite" ON)

# Find required packages
find_package(OpenMP REQUIRED)
find_package(OpenSSL REQUIRED)

# Add subdirectories
add_subdirectory(src)

# Add tests if enabled
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Create symbolic link to compile_commands.json in source directory
if(CMAKE_EXPORT_COMPILE_COMMANDS)
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_BINARY_DIR}/compile_commands.json
        ${CMAKE_SOURCE_DIR}/compile_commands.json
    )
endif()
