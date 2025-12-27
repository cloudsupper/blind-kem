#include <gtest/gtest.h>

#include <kem.h>
#include <polynomial.h>

#include <random>
#include <vector>

namespace {

// For each RLWE parameter set, compare NTT-based Polynomial::operator*
// against the reference schoolbook implementation for random inputs.
static void check_ntt_multiply_matches_schoolbook(const RLWEParams& params,
                                                  std::uint64_t seed) {
    const std::size_t n = params.n;
    const std::uint64_t q = params.q;

    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<std::uint64_t> dist(0, q - 1);

    auto random_poly = [&](void) {
        std::vector<std::uint64_t> coeffs(n);
        for (std::size_t i = 0; i < n; ++i) {
            coeffs[i] = dist(rng);
        }
        return Polynomial(coeffs, q);
    };

    // Local schoolbook reference (same logic as multiplySchoolbook).
    auto schoolbook = [&](const Polynomial& a, const Polynomial& b) {
        const auto& ac = a.getCoeffs();
        const auto& bc = b.getCoeffs();

        std::vector<std::uint64_t> temp(2 * n, 0);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                std::uint64_t prod =
                    (static_cast<std::uint64_t>(ac[i]) *
                     static_cast<std::uint64_t>(bc[j])) % q;
                temp[i + j] = (temp[i + j] + prod) % q;
            }
        }

        std::vector<std::uint64_t> res(n, 0);
        for (std::size_t i = 0; i < n; ++i) {
            std::uint64_t coeff = temp[i];
            std::size_t higher_degree = i + n;
            while (higher_degree < temp.size()) {
                const std::int64_t diff = static_cast<std::int64_t>(coeff) -
                                          static_cast<std::int64_t>(temp[higher_degree]);
                const std::int64_t m = static_cast<std::int64_t>(q);
                std::int64_t r = diff % m;
                if (r < 0) r += m;
                coeff = static_cast<std::uint64_t>(r);
                higher_degree += n;
            }
            res[i] = coeff;
        }

        return Polynomial(res, q);
    };

    // Try a few random pairs.
    for (int t = 0; t < 5; ++t) {
        Polynomial a = random_poly();
        Polynomial b = random_poly();

        Polynomial ref = schoolbook(a, b);
        Polynomial got = a * b;  // This uses NTT where available.

        ASSERT_EQ(ref.degree(), got.degree());
        const auto& rc = ref.getCoeffs();
        const auto& gc = got.getCoeffs();
        for (std::size_t i = 0; i < n; ++i) {
            EXPECT_EQ(rc[i], gc[i])
                << "Mismatch at index " << i
                << " for n=" << n << ", q=" << q;
        }
    }
}

} // namespace

TEST(PolynomialNTTMultiplyTest, AllSecurityLevelsMatchSchoolbook) {
    using std::pair;
    using std::vector;

    vector<pair<RLWEParams, uint64_t>> configs = {
        {KEM::getParameterSet(SecurityLevel::TEST_TINY), 0x0102030405060708ULL},
        {KEM::getParameterSet(SecurityLevel::TEST_SMALL), 0x1112131415161718ULL},
        {KEM::getParameterSet(SecurityLevel::KYBER512), 0x2122232425262728ULL},
        {KEM::getParameterSet(SecurityLevel::MODERATE), 0x3132333435363738ULL},
        {KEM::getParameterSet(SecurityLevel::HIGH), 0x4142434445464748ULL}
    };

    for (const auto& [params, seed] : configs) {
        check_ntt_multiply_matches_schoolbook(params, seed);
    }
}
