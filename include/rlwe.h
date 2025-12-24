#ifndef RLWE_H
#define RLWE_H

#include <cmath>
#include <polynomial.h>
#include <vector>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <logging.h>

/**
 * @brief Supported security levels for the RLWE signature scheme.
 *
 * The values correspond to different parameter sets based on
 * NIST and lattice-based cryptography recommendations. Some
 * levels are intended only for testing and are not secure.
 */
enum class SecurityLevel {
    /**
     * @brief Tiny, insecure test parameters.
     *
     * Very small ring dimension and modulus. Only suitable for
     * functional tests and examples, not for any real security.
     */
    TEST_TINY,

    /**
     * @brief Small, insecure test parameters.
     *
     * Slightly larger than TEST_TINY but still far below
     * recommended security margins.
     */
    TEST_SMALL,

    /**
     * @brief Parameters corresponding roughly to NIST Kyber-512.
     *
     * Provides approximately 128 bits of classical security and
     * around 64 bits of quantum security.
     */
    KYBER512,

    /**
     * @brief Moderate security parameter set.
     *
     * Intended to provide ~192 bits of classical security and
     * ~96 bits of quantum security.
     */
    MODERATE,

    /**
     * @brief High security parameter set.
     *
     * Intended to provide ~256 bits of classical security and
     * ~128 bits of quantum security.
     */
    HIGH,
};

/**
 * @brief Describes a concrete RLWE parameter set.
 */
struct RLWEParams {
    /** Ring dimension (polynomial degree bound). */
    size_t n;
    /** Coefficient modulus. */
    uint64_t q;
    /** Gaussian standard deviation for noise sampling. */
    double sigma;
    /** Human-readable name of the parameter set. */
    const char* name;
    /** Estimated classical security level in bits. */
    int classical_bits;
    /** Estimated quantum security level in bits. */
    int quantum_bits;
    /** Whether this parameter set is considered cryptographically secure. */
    bool is_secure;
};

/**
 * @brief RLWE-based blind signature scheme implementation.
 *
 * This class implements a simple RLWE blind signature construction
 * over the ring @f$Z_q[x]/(x^n + 1)@f$. It supports key generation,
 * message blinding, blind signing, and verification.
 */
class BlindKEM {
public:
    /**
     * @brief Construct an RLWE instance with explicit parameters.
     *
     * @param n Ring dimension (must be a power of two).
     * @param q Coefficient modulus.
     * @param sigma Standard deviation of the discrete Gaussian
     *             used for noise sampling. If zero or negative,
     *             a reasonable default is chosen.
     *
     * @throws std::invalid_argument If @p n is not a power of two.
     */
    BlindKEM(size_t n, uint64_t q, double sigma = 0.0);

    /**
     * @brief Construct an RLWE instance from a named security level.
     *
     * This constructor selects a predefined parameter set based on
     * @p level. It is the recommended way to construct an instance
     * for most applications.
     *
     * @param level Desired security level.
     *
     * @throws std::invalid_argument If the derived ring dimension is
     *         not a power of two.
     */
    explicit BlindKEM(SecurityLevel level = SecurityLevel::KYBER512);

    /**
     * @brief Generate a fresh key pair.
     *
     * Samples a uniform public polynomial @f$a@f$, a secret key
     * polynomial @f$s@f$ from a discrete Gaussian, and an error
     * polynomial @f$e@f$. The public key is @f$(a, b = a s + e)@f$.
     */
    void generateKeys();

    /**
     * @brief Compute a blind signature on a blinded message polynomial.
     *
     * @param blindedMessage Blinded message polynomial.
     * @return Blind signature polynomial.
     */
    Polynomial blindSign(const Polynomial& blindedMessage);

    /**
     * @brief Verify a signature on a message.
     *
     * The message is first hashed to a polynomial, and the signature
     * is checked against the expected value derived from the secret
     * key. Verification is performed in a noise-tolerant manner using
     * a binary signal representation.
     *
     * @param secret Message bytes to verify.
     * @param signature Signature polynomial.
     * @return true if the signature is accepted, false otherwise.
     */
    bool verify(const std::vector<uint8_t>& secret,
                const Polynomial& signature);

    /**
     * @brief Retrieve the public key.
     *
     * @return Pair (a, b) representing the public key.
     */
    std::pair<Polynomial, Polynomial> getPublicKey() const {
        return std::make_pair(a, b);
    }

    /**
     * @brief Hash a message to a polynomial with coefficients in {0, q/2}.
     *
     * A counter-based SHA-256 construction is used to fill the
     * coefficients deterministically.
     *
     * @param message Input message bytes.
     * @return Polynomial representation of the hash.
     */
    Polynomial hashToPolynomial(const std::vector<uint8_t>& message);

    /**
     * @brief Compute a blinded message polynomial and blinding factor.
     *
     * Given a message, this method samples a blinding factor and
     * returns the blinded message as well as the factor required to
     * unblind the final signature.
     *
     * @param secret Message to blind.
     * @return Pair (blindedMessage, blindingFactor).
     */
    std::pair<Polynomial, Polynomial> computeBlindedMessage(const std::vector<uint8_t>& secret);

    /**
     * @brief Compute the final signature from a blind signature.
     *
     * @param blindSignature Signature on the blinded message.
     * @param blindingFactor Blinding factor used when creating
     *        the blinded message.
     * @param publicKey Public key component used for unblinding.
     * @return Unblinded signature polynomial.
     */
    Polynomial computeSignature(const Polynomial& blindSignature,
                                const Polynomial& blindingFactor,
                                const Polynomial& publicKey);

    /**
     * @brief Get the current effective parameters.
     *
     * @return Description of the active parameter set.
     */
    RLWEParams getParameters() const;

    /**
     * @brief Retrieve a predefined parameter set for a security level.
     *
     * @param level Named security level.
     * @return Corresponding parameter set.
     */
    static RLWEParams getParameterSet(SecurityLevel level);

    /**
     * @brief Expose the secret key polynomial for testing/attacks.
     *
     * WARNING: This must never be used in production code. It is
     * provided solely to allow experiments like oracle_cca to
     * compare their recovered secret against the ground truth.
     */
    Polynomial getSecretKeyForTesting() const { return s; }

private:
    size_t ring_dim_n;
    uint64_t modulus;
    double gaussian_stddev;

    Polynomial a;
    Polynomial b;
    Polynomial s;

    /**
     * @brief Generate a uniformly random 64-bit integer.
     *
     * Uses a platform-specific secure random source.
     */
    uint64_t getRandomUint64();

    /**
     * @brief Generate a random double from a standard normal distribution.
     */
    double getRandomDouble();

    /**
     * @brief Sample a polynomial with coefficients uniformly in [0, q).
     */
    Polynomial sampleUniform();

    /**
     * @brief Sample a polynomial with coefficients drawn from a
     *        discretized Gaussian distribution.
     *
     * @param stddev Standard deviation of the Gaussian.
     */
    Polynomial sampleGaussian(double stddev);

    /**
     * @brief Encode a message as a polynomial with 0/1 coefficients.
     *
     * The message bits are packed into coefficients starting from
     * the most significant bit of each byte.
     *
     * @param message Message bytes.
     * @return Polynomial with coefficients in {0,1}.
     */
    Polynomial messageToPolynomial(const std::vector<uint8_t>& message);

    /**
     * @brief Log and validate the chosen security parameters.
     */
    void validateSecurityParameters();

    static constexpr double LARGE_THRESHOLD_DIVISOR = 4.0;
    static constexpr size_t MIN_DIFFERENT_COEFFS = 1;

    /**
     * @brief Helper to log a byte sequence in hexadecimal form.
     *
     * @param prefix Descriptive prefix for the log line.
     * @param message Bytes to log.
     */
    void logMessageBytes(const std::string& prefix, const std::vector<uint8_t>& message) {
        std::stringstream ss;
        ss << prefix << " bytes: [";
        for (size_t i = 0; i < message.size(); ++i) {
            if (i > 0) {
                ss << ", ";
            }
            ss << std::hex << std::uppercase << std::setw(2) << std::setfill('0')
               << static_cast<int>(message[i]);
        }
        ss << "]";
        Logger::log(ss.str());
    }
};

#endif // RLWE_H
