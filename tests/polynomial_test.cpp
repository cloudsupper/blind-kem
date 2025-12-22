#include <gtest/gtest.h>
#include <polynomial.h>

class PolynomialTest : public ::testing::Test {
protected:
    // Using parameters for easier testing
    // Z[x]/(x^4 + 1) mod 17
    const size_t n = 2;  // degree will be 2n = 4
    const uint64_t q = 17;
};

TEST_F(PolynomialTest, ToBytes) {
    // Create a polynomial
    Polynomial p(4, 17);
    std::vector<uint64_t> coeffs = {1, 2, 3, 4};
    p.setCoefficients(coeffs);
    
    // Get bytes
    auto bytes = p.toBytes();
    
    // Expected size: sizeof(size_t) for ring_dim + sizeof(uint64_t) for modulus + 
    // 4 * sizeof(uint64_t) for coefficients
    size_t expected_size = sizeof(size_t) + sizeof(uint64_t) + 4 * sizeof(uint64_t);
    EXPECT_EQ(bytes.size(), expected_size);
    
    // Create another identical polynomial
    Polynomial p2(4, 17);
    p2.setCoefficients(coeffs);
    
    // Verify that identical polynomials produce identical byte sequences
    EXPECT_EQ(p.toBytes(), p2.toBytes());
    
    // Modify one coefficient
    std::vector<uint64_t> coeffs2 = {1, 2, 3, 5};
    p2.setCoefficients(coeffs2);
    
    // Verify that different polynomials produce different byte sequences
    EXPECT_NE(p.toBytes(), p2.toBytes());
}

TEST_F(PolynomialTest, Addition) {
    // Create polynomials: f = 1 + 2x + 3x^2 + 4x^3, g = 5 + 6x + 7x^2 + 8x^3
    Polynomial f({1, 2, 3, 4}, q);
    Polynomial g({5, 6, 7, 8}, q);
    
    // Expected: (1+5) + (2+6)x + (3+7)x^2 + (4+8)x^3 mod 17
    // = 6 + 8x + 10x^2 + 12x^3
    Polynomial h = f + g;
    
    EXPECT_EQ(h[0], 6);
    EXPECT_EQ(h[1], 8);
    EXPECT_EQ(h[2], 10);
    EXPECT_EQ(h[3], 12);
}

// Rest of the test file...

TEST_F(PolynomialTest, Subtraction) {
    // f = 1 + 2x + 3x^2 + 4x^3
    // g = 5 + 6x + 7x^2 + 8x^3
    // f - g (mod 17) = (1-5, 2-6, 3-7, 4-8) mod 17 = (13, 13, 13, 13)
    Polynomial f({1, 2, 3, 4}, q);
    Polynomial g({5, 6, 7, 8}, q);

    Polynomial h = f - g;

    EXPECT_EQ(h[0], 13u);
    EXPECT_EQ(h[1], 13u);
    EXPECT_EQ(h[2], 13u);
    EXPECT_EQ(h[3], 13u);
}

TEST_F(PolynomialTest, Negation) {
    // p = 0 + 1x + 16x^2 + 8x^3 over Z_17
    // -p = (0, 16, 1, 9)
    Polynomial p({0, 1, 16, 8}, q);

    Polynomial neg = -p;

    EXPECT_EQ(neg[0], 0u);      // -0 ≡ 0 (mod 17)
    EXPECT_EQ(neg[1], 16u);     // -1 ≡ 16 (mod 17)
    EXPECT_EQ(neg[2], 1u);      // -16 ≡ 1 (mod 17)
    EXPECT_EQ(neg[3], 9u);      // -8 ≡ 9 (mod 17)
}

TEST_F(PolynomialTest, MultiplicationByOneAndZero) {
    // ring_dim = 4 so we work in Z_17[x]/(x^4 + 1)
    Polynomial one({1, 0, 0, 0}, q);   // represents 1
    Polynomial zero({0, 0, 0, 0}, q);  // represents 0
    Polynomial f({3, 5, 7, 9}, q);

    Polynomial f_times_one = f * one;
    Polynomial one_times_f = one * f;
    Polynomial f_times_zero = f * zero;
    Polynomial zero_times_f = zero * f;

    for (size_t i = 0; i < f.degree(); ++i) {
        EXPECT_EQ(f_times_one[i], f[i]);
        EXPECT_EQ(one_times_f[i], f[i]);
        EXPECT_EQ(f_times_zero[i], 0u);
        EXPECT_EQ(zero_times_f[i], 0u);
    }
}

TEST_F(PolynomialTest, MultiplicationWrapAroundX4Plus1) {
    // In Z_17[x]/(x^4 + 1), we have x * x^3 = x^4 ≡ -1 (mod x^4 + 1)
    // which corresponds to (16, 0, 0, 0) in coefficient form
    Polynomial x({0, 1, 0, 0}, q);   // x
    Polynomial x3({0, 0, 0, 1}, q);  // x^3

    Polynomial prod1 = x * x3;
    Polynomial prod2 = x3 * x;

    EXPECT_EQ(prod1[0], 16u);
    EXPECT_EQ(prod1[1], 0u);
    EXPECT_EQ(prod1[2], 0u);
    EXPECT_EQ(prod1[3], 0u);

    EXPECT_EQ(prod2[0], 16u);
    EXPECT_EQ(prod2[1], 0u);
    EXPECT_EQ(prod2[2], 0u);
    EXPECT_EQ(prod2[3], 0u);
}

TEST_F(PolynomialTest, ScalarMultiplication) {
    // f = 1 + 2x + 3x^2 + 4x^3, scalar = 5
    // 5 * f (mod 17) = (5, 10, 15, 3)
    Polynomial f({1, 2, 3, 4}, q);
    uint64_t scalar = 5;

    Polynomial h = f * scalar;

    EXPECT_EQ(h[0], 5u);
    EXPECT_EQ(h[1], 10u);
    EXPECT_EQ(h[2], 15u);
    EXPECT_EQ(h[3], 3u);
}

TEST_F(PolynomialTest, PolySignalRounding) {
    // For q = 17, half_mod = 8. The rounding chooses whichever of 0 or 8
    // is closer on the cyclic group Z_17, with ties going to 0.
    // Coefficients: 0,4,5,7,8,9,13,16 -> 0,0,8,8,8,8,0,0
    Polynomial p({0, 4, 5, 7, 8, 9, 13, 16}, q);

    Polynomial s = p.polySignal();

    std::vector<uint64_t> expected = {0, 0, 8, 8, 8, 8, 0, 0};
    ASSERT_EQ(s.degree(), expected.size());
    for (size_t i = 0; i < expected.size(); ++i) {
        EXPECT_EQ(s[i], expected[i]) << "Mismatch at index " << i;
    }
}

TEST_F(PolynomialTest, SetCoefficientsModReductionAndSizeCheck) {
    Polynomial p(4, q);

    // Values larger than q should be reduced modulo q
    std::vector<uint64_t> coeffs = {q + 1, 2 * q, 0, q - 1};
    p.setCoefficients(coeffs);

    EXPECT_EQ(p[0], 1u);         // (q + 1) mod q
    EXPECT_EQ(p[1], 0u);         // (2q) mod q
    EXPECT_EQ(p[2], 0u);
    EXPECT_EQ(p[3], (q - 1));    // already in range

    // Setting coefficients with the wrong size should throw
    std::vector<uint64_t> wrong_size = {1, 2, 3};
    EXPECT_THROW(p.setCoefficients(wrong_size), std::invalid_argument);
}

TEST_F(PolynomialTest, RingDimensionAndModulusAccessors) {
    Polynomial p({1, 2, 3, 4}, q);

    EXPECT_EQ(p.degree(), 4u);
    EXPECT_EQ(p.getModulus(), q);
}

TEST_F(PolynomialTest, OperationsRequireSameRingAndModulus) {
    Polynomial f({1, 2, 3, 4}, q);
    Polynomial different_dim({1, 2}, q);       // ring_dim = 2
    Polynomial different_mod({1, 2, 3, 4}, q + 1);  // same dim, different modulus

    EXPECT_THROW(f + different_dim, std::invalid_argument);
    EXPECT_THROW(f - different_dim, std::invalid_argument);
    EXPECT_THROW(f * different_dim, std::invalid_argument);

    EXPECT_THROW(f + different_mod, std::invalid_argument);
    EXPECT_THROW(f - different_mod, std::invalid_argument);
    EXPECT_THROW(f * different_mod, std::invalid_argument);
}
