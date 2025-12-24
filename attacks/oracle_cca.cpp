#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "rlwe.h"
#include "polynomial.h"
#include "logging.h"

// -----------------------------------------------------------------------------
// oracle_cca.cpp
//
// CCA-style oracle attack demo against the RLWE blind signature / KEM scheme.
//
// Idea:
//   - The signer holds a secret key polynomial s(x).
//   - The blind signing algorithm, for an input polynomial m(x), computes
//         d(x) = s(x) * m(x) + e(x)
//     where e(x) is fresh noise each time.
//   - If the attacker can choose m(x) = 1, then
//         d(x) = s(x) * 1 + e(x) = s(x) + e(x).
//   - Repeating this query many times and averaging coefficients of d(x),
//     the noise averages out and we recover an estimate of s(x).
//
// This program does exactly that using the existing BlindKEM implementation:
//   1. Instantiate BlindKEM with test parameters (insecure, but fast).
//   2. Generate a fresh key pair, which samples secret s(x).
//   3. Build the polynomial p(x) = 1 in the ring Z_q[x]/(x^n + 1).
//   4. Query blindSign(p) many times to obtain samples d_i(x).
//   5. Convert coefficients to a signed representation in (-q/2, q/2].
//   6. Average coefficients across all queries to estimate s(x).
//   7. Compare the recovered s_hat(x) against the true secret s(x) and
//      report whether the attack was successful and which coefficients
//      deviate (and by how much).
//
// NOTE: This is *not* a secure usage pattern. It is an attack demo
// to explore how malformed / chosen-message queries can leak the
// underlying RLWE secret key via simple statistics.
// -----------------------------------------------------------------------------

// Convert a coefficient in [0, q-1] to a signed representative in (-q/2, q/2].
static long double to_signed(uint64_t c, uint64_t q) {
    if (c > q / 2) {
        return static_cast<long double>(c) - static_cast<long double>(q);
    }
    return static_cast<long double>(c);
}

int main(int argc, char** argv) {
    // Parse number of oracle queries from command line (default: 10,000)
    size_t num_samples = 10000;
    if (argc >= 2) {
        try {
            num_samples = static_cast<size_t>(std::stoull(argv[1]));
        } catch (const std::exception&) {
            std::cerr << "[oracle_cca] Invalid sample count '" << argv[1]
                      << "', using default " << num_samples << "\n";
        }
    }

    // Keep the library logging mostly quiet for the attack demo
    Logger::setOutputStream(std::cout);
    Logger::enable_logging = false;

    std::cout << "[oracle_cca] RLWE CCA-style oracle attack demo" << std::endl;
    std::cout << "[oracle_cca] Using num_samples = " << num_samples << "\n";

    // 1. Instantiate RLWE scheme with insecure TEST_SMALL parameters for speed.
    //    You can switch to KYBER512 to mimic more realistic parameters.
    BlindKEM rlwe(SecurityLevel::MODERATE);
    rlwe.generateKeys();

    RLWEParams params = rlwe.getParameters();
    const size_t n = params.n;
    const uint64_t q = params.q;

    std::cout << "[oracle_cca] Parameters: n=" << n
              << ", q=" << q
              << ", sigma=" << params.sigma << "\n";

    // 2. Build p(x) = 1 in Z_q[x]/(x^n + 1): coefficient 0 is 1, rest are 0.
    Polynomial p(n, q);
    {
        std::vector<uint64_t> coeffs(n, 0);
        coeffs[0] = 1;
        p.setCoefficients(coeffs);
    }

    std::cout << "[oracle_cca] Built oracle query polynomial p(x) = 1" << std::endl;

    // 3. Repeatedly query the blind signing oracle with p(x) = 1.
    //    Each response d(x) = s(x) + e(x) (since s * 1 = s).
    std::vector<long double> sum_coeffs(n, 0.0L);

    for (size_t t = 0; t < num_samples; ++t) {
        Polynomial d = rlwe.blindSign(p);
        const auto& dc = d.getCoeffs();

        for (size_t i = 0; i < n; ++i) {
            sum_coeffs[i] += to_signed(dc[i], q);
        }

        if ((t + 1) % (num_samples / 10 == 0 ? 1 : num_samples / 10) == 0) {
            std::cout << "[oracle_cca] Collected " << (t + 1)
                      << " / " << num_samples << " samples..." << std::endl;
        }
    }

    // 4. Compute averaged coefficients as an estimate s_hat(x) ≈ s(x).
    std::vector<long double> avg_coeffs(n, 0.0L);
    std::vector<long long>  s_hat(n, 0);

    long double inv_T = 1.0L / static_cast<long double>(num_samples);
    for (size_t i = 0; i < n; ++i) {
        avg_coeffs[i] = sum_coeffs[i] * inv_T;
        s_hat[i] = static_cast<long long>(std::llround(avg_coeffs[i]));
    }

    // 5. Obtain the true secret key polynomial and convert to signed reps.
    Polynomial s_true_poly = rlwe.getSecretKeyForTesting();
    const auto& s_true_modq = s_true_poly.getCoeffs();

    std::vector<long double> s_true_signed(n, 0.0L);
    std::vector<long long>  s_true_int(n, 0);
    for (size_t i = 0; i < n; ++i) {
        s_true_signed[i] = to_signed(s_true_modq[i], q);
        s_true_int[i] = static_cast<long long>(std::llround(s_true_signed[i]));
    }

    // 6. Compare recovered coefficients with the true secret.
    bool all_match = true;
    size_t mismatch_count = 0;
    long double max_abs_err = 0.0L;
    size_t max_err_idx = 0;

    std::vector<long long> delta_int(n, 0);

    for (size_t i = 0; i < n; ++i) {
        long long diff = s_hat[i] - s_true_int[i];
        delta_int[i] = diff;
        if (diff != 0) {
            all_match = false;
            ++mismatch_count;
        }

        long double err = std::fabs(static_cast<long double>(avg_coeffs[i] - s_true_signed[i]));
        if (err > max_abs_err) {
            max_abs_err = err;
            max_err_idx = i;
        }
    }

    // 7. Print a summary of the recovered secret coefficients and errors.
    std::cout << "\n[oracle_cca] Estimated secret polynomial coefficients (signed)" << std::endl;
    std::cout << "[oracle_cca] Showing first 32 coefficients (or all if n < 32):" << std::endl;

    size_t limit = std::min<size_t>(32, n);
    for (size_t i = 0; i < limit; ++i) {
        std::cout << "  i=" << i
                  << ": s_hat=" << s_hat[i]
                  << ", s_true=" << s_true_int[i]
                  << ", delta=" << delta_int[i]
                  << " (avg=" << static_cast<double>(avg_coeffs[i])
                  << ", true_signed=" << static_cast<double>(s_true_signed[i])
                  << ")\n";
    }

    if (n > limit) {
        std::cout << "[oracle_cca] (" << (n - limit)
                  << " more coefficients omitted in preview)" << std::endl;
    }

    // 8. Report attack success/failure and list deviating coefficients.
    if (all_match) {
        std::cout << "\n[oracle_cca] ATTACK SUCCESS: recovered secret exactly for all "
                  << n << " coefficients." << std::endl;
    } else {
        std::cout << "\n[oracle_cca] ATTACK PARTIAL: " << mismatch_count
                  << " / " << n << " coefficients deviate from the true secret." << std::endl;
        std::cout << "[oracle_cca] Maximum absolute error "
                  << static_cast<double>(max_abs_err)
                  << " at index i=" << max_err_idx << "." << std::endl;

        std::cout << "[oracle_cca] Listing all deviating coefficients" << std::endl;
        std::cout << "[oracle_cca] Format: i: s_hat, s_true, delta, avg, true_signed" << std::endl;
        for (size_t i = 0; i < n; ++i) {
            if (delta_int[i] != 0) {
                std::cout << "  i=" << i
                          << ": s_hat=" << s_hat[i]
                          << ", s_true=" << s_true_int[i]
                          << ", delta=" << delta_int[i]
                          << ", avg=" << static_cast<double>(avg_coeffs[i])
                          << ", true_signed=" << static_cast<double>(s_true_signed[i])
                          << "\n";
            }
        }
    }

    std::cout << "\n[oracle_cca] Attack complete."
              << " The closer sigma is to 0 and the more samples you take,"
              << " the more accurate this estimate of s(x) becomes.\n";

    return 0;
}
