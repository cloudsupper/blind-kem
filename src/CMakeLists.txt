# Create the library
add_library(rlwe
    rlwe.cpp
    polynomial.cpp
    ntt.cpp
    sha256.cpp
)

# Add include directories
target_include_directories(rlwe
    PUBLIC
        ${CMAKE_SOURCE_DIR}/include
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link with OpenMP and OpenSSL
target_link_libraries(rlwe 
    PRIVATE 
        OpenMP::OpenMP_CXX
        OpenSSL::SSL 
        OpenSSL::Crypto
)

# ---------------------------------------------------------------------------
# Attack / demo binaries
# ---------------------------------------------------------------------------

add_executable(oracle_cca
    ${CMAKE_SOURCE_DIR}/attacks/oracle_cca.cpp
)

target_link_libraries(oracle_cca PRIVATE rlwe)
