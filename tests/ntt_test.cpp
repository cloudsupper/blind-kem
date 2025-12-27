#include <gtest/gtest.h>

#include <kem.h>
#include <ntt.h>

#include <random>

namespace {

// Helper: check that for given (n, q) the NTT roundtrip on Polynomial
// recovers the original polynomial exactly.
static void check_ntt_roundtrip_for_params(const RLWEParams& params,
                                           std::uint64_t seed) {
    const std::size_t n = params.n;
    const std::uint64_t q = params.q;

    NTT ntt(n, q, /*negacyclic=*/true);

    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<std::uint64_t> dist(0, q - 1);

    // Try several random polynomials plus some simple structured ones.
    auto run_single = [&](const std::vector<std::uint64_t>& coeffs) {
        ASSERT_EQ(coeffs.size(), n);
        Polynomial p(coeffs, q);
        Polynomial original = p;

        ntt.forward(p);
        ntt.inverse(p);

        const auto& out = p.getCoeffs();
        const auto& orig = original.getCoeffs();
        ASSERT_EQ(out.size(), orig.size());
        for (std::size_t i = 0; i < n; ++i) {
            EXPECT_EQ(out[i], orig[i])
                << "Mismatch at index " << i
                << " for n=" << n << ", q=" << q;
        }
    };

    // Zero polynomial
    run_single(std::vector<std::uint64_t>(n, 0));

    // Single 1 at each position (delta basis vectors)
    for (std::size_t pos = 0; pos < std::min<std::size_t>(n, 8); ++pos) {
        std::vector<std::uint64_t> v(n, 0);
        v[pos] = 1;
        run_single(v);
    }

    // Increasing pattern modulo q
    {
        std::vector<std::uint64_t> v(n);
        for (std::size_t i = 0; i < n; ++i) {
            v[i] = static_cast<std::uint64_t>(i) % q;
        }
        run_single(v);
    }

    // A few random polynomials
    for (int t = 0; t < 5; ++t) {
        std::vector<std::uint64_t> v(n);
        for (std::size_t i = 0; i < n; ++i) {
            v[i] = dist(rng);
        }
        run_single(v);
    }
}

} // namespace

TEST(NTTTest, RoundtripAllSecurityLevels) {
    // Test Tiny / Small (insecure, but good for fast tests)
    check_ntt_roundtrip_for_params(
        KEM::getParameterSet(SecurityLevel::TEST_TINY), 0xA1B2C3D4ULL);
    check_ntt_roundtrip_for_params(
        KEM::getParameterSet(SecurityLevel::TEST_SMALL), 0xBEEF1234ULL);

    // Kyber-like 128-bit, NTT-friendly parameters
    check_ntt_roundtrip_for_params(
        KEM::getParameterSet(SecurityLevel::KYBER512), 0x12345678ULL);

    // Moderate and High security levels
    check_ntt_roundtrip_for_params(
        KEM::getParameterSet(SecurityLevel::MODERATE), 0xCAFEBABEULL);
    check_ntt_roundtrip_for_params(
        KEM::getParameterSet(SecurityLevel::HIGH), 0xDEADBEEFULL);
}
