# Find or download Google Test
include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG release-1.12.1
)
FetchContent_MakeAvailable(googletest)

# Add test executable
add_executable(kem_tests
    polynomial_test.cpp
    sha256_test.cpp
    ntt_test.cpp
    polynomial_ntt_multiply_test.cpp
)

# Link against Google Test and our library
target_link_libraries(kem_tests
    PRIVATE
        gtest_main
        kem
)

# Add the test to CTest
add_test(NAME kem_tests COMMAND kem_tests)
