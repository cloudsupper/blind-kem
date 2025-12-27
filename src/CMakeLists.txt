# Create the library
add_library(kem
    kem.cpp
    polynomial.cpp
    ntt.cpp
    sha256.cpp
)

# Add include directories
target_include_directories(kem
    PUBLIC
        ${CMAKE_SOURCE_DIR}/include
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link with OpenMP and OpenSSL
target_link_libraries(kem 
    PRIVATE 
        OpenMP::OpenMP_CXX
        OpenSSL::SSL 
        OpenSSL::Crypto
)
