#!/bin/bash

# Create build directory if it doesn't exist
mkdir -p build

# Enter build directory
cd build

# Configure CMake
cmake ..

# Build the project
cmake --build . -j$(nproc)

# Run tests if they were built
if [ -f tests/kem_tests ]; then
    ctest --output-on-failure
fi
