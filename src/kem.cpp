#include <polynomial.h>
#include <kem.h>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <limits>
#include <random>
#include <sha256.h>

#if defined(_WIN32)
#include <windows.h>
#include <bcrypt.h>
#elif defined(__APPLE__)
#include <Security/SecRandom.h>
#endif

static void getSecureRandomBytes(uint8_t* buffer, size_t length) {
#if defined(_WIN32)
    BCRYPT_ALG_HANDLE hAlg = NULL;
    NTSTATUS status = BCryptOpenAlgorithmProvider(&hAlg, BCRYPT_RNG_ALGORITHM, NULL, 0);
    if (!BCRYPT_SUCCESS(status)) {
        throw std::runtime_error("Failed to open BCrypt algorithm provider");
    }
    
    status = BCryptGenRandom(hAlg, buffer, static_cast<ULONG>(length), 0);
    BCryptCloseAlgorithmProvider(hAlg, 0);
    
    if (!BCRYPT_SUCCESS(status)) {
        throw std::runtime_error("Failed to generate random bytes using BCrypt");
    }
#elif defined(__APPLE__)
    if (SecRandomCopyBytes(kSecRandomDefault, length, buffer) != 0) {
        throw std::runtime_error("Failed to generate random bytes using SecRandomCopyBytes");
    }
#else
    std::random_device rd("/dev/urandom");
    if (!rd.entropy()) {
        throw std::runtime_error("Failed to access secure random source");
    }
    
    for (size_t i = 0; i < length; i += sizeof(uint32_t)) {
        uint32_t random = rd();
        size_t remaining = std::min(sizeof(uint32_t), length - i);
        std::memcpy(buffer + i, &random, remaining);
    }
#endif
}

uint64_t KEM::getRandomUint64() {
    uint64_t result;
    getSecureRandomBytes(reinterpret_cast<uint8_t*>(&result), sizeof(result));
    return result;
}

double KEM::getRandomDouble() {
    uint64_t r1, r2;
    getSecureRandomBytes(reinterpret_cast<uint8_t*>(&r1), sizeof(r1));
    getSecureRandomBytes(reinterpret_cast<uint8_t*>(&r2), sizeof(r2));
    
    double u1 = static_cast<double>(r1) / std::numeric_limits<uint64_t>::max();
    double u2 = static_cast<double>(r2) / std::numeric_limits<uint64_t>::max();
    
    double radius = std::sqrt(-2 * std::log(u1));
    double theta = 2 * M_PI * u2;
    
    return radius * std::cos(theta);
}

static bool isPowerOfTwo(size_t n) {
    return n != 0 && (n & (n - 1)) == 0;
}

static bool validatePowerOfTwo(size_t n) {
    if (!isPowerOfTwo(n)) {
        return false;
    }

#if defined(_MSC_VER)
    unsigned long index;
    _BitScanReverse(&index, static_cast<unsigned long>(n));
    return (static_cast<size_t>(1) << index) == n;
#elif defined(__GNUC__) || defined(__clang__)
    unsigned long long n_ull = static_cast<unsigned long long>(n);
    return __builtin_popcountll(n_ull) == 1;
#else
    size_t count = 0;
    size_t temp = n;
    while (temp > 0) {
        count += temp & 1;
        temp >>= 1;
    }
    return count == 1;
#endif
}

RLWEParams KEM::getParameterSet(SecurityLevel level) {
    switch (level) {
        case SecurityLevel::TEST_TINY:
            return {8, 7681, 3.0, "TEST_TINY (INSECURE)", 4, 2, false};
        
        case SecurityLevel::TEST_SMALL:
            return {32, 7681, 3.0, "TEST_SMALL (INSECURE)", 16, 8, false};
        
        case SecurityLevel::KYBER512:
            return {256, 7681, 3.0, "KYBER512-like (NTT-friendly)", 128, 64, true};
        
        case SecurityLevel::MODERATE:
            return {512, 12289, 3.2, "MODERATE", 192, 96, true};
        
        case SecurityLevel::HIGH:
            return {1024, 18433, 3.2, "HIGH", 256, 128, true};
        
        default:
            return {256, 7681, 3.0, "KYBER512-like (NTT-friendly)", 128, 64, true};
    }
}

KEM::KEM(size_t n, uint64_t q, double sigma)
    : ring_dim_n(n),
      modulus(q),
      gaussian_stddev(sigma > 0 ? sigma : 3.2),
      a(n, q),
      b(n, q),
      s(n, q)
{
    if (!validatePowerOfTwo(n)) {
        throw std::invalid_argument("n must be a power of 2");
    }

    Logger::log("Created RLWE instance with n=" + std::to_string(n) + 
                ", q=" + std::to_string(q) + ", σ=" + std::to_string(gaussian_stddev));
    
    validateSecurityParameters();
}

KEM::KEM(SecurityLevel level) 
    : ring_dim_n(0),
      modulus(0),
      gaussian_stddev(0),
      a(1, 1),
      b(1, 1),
      s(1, 1)
{
    RLWEParams params = getParameterSet(level);
    
    ring_dim_n = params.n;
    modulus = params.q;
    gaussian_stddev = params.sigma;
    
    a = Polynomial(params.n, params.q);
    b = Polynomial(params.n, params.q);
    s = Polynomial(params.n, params.q);
    
    if (!validatePowerOfTwo(ring_dim_n)) {
        throw std::invalid_argument("n must be a power of 2");
    }
    
    Logger::log("\n" + std::string(70, '='));
    Logger::log("RLWE INSTANCE CREATED");
    Logger::log(std::string(70, '='));
    Logger::log("Security Level: " + std::string(params.name));
    Logger::log("Parameters: n=" + std::to_string(params.n) + 
                ", q=" + std::to_string(params.q) + 
                ", σ=" + std::to_string(params.sigma));
    Logger::log("Estimated Security:");
    Logger::log("  Classical: ~" + std::to_string(params.classical_bits) + " bits");
    Logger::log("  Quantum:   ~" + std::to_string(params.quantum_bits) + " bits");
    
    if (!params.is_secure) {
        Logger::log("\n⚠️  WARNING: INSECURE PARAMETERS ⚠️");
        Logger::log("These parameters provide insufficient security!");
        Logger::log("Only use for testing and development.");
        Logger::log("For production, use SecurityLevel::KYBER512 or higher.");
        Logger::log("⚠️  DO NOT USE IN PRODUCTION ⚠️\n");
    } else {
        Logger::log("\n✓ Parameters meet cryptographic security requirements");
    }
    Logger::log(std::string(70, '=') + "\n");
    
    validateSecurityParameters();
}

RLWEParams KEM::getParameters() const {
    RLWEParams params;
    params.n = ring_dim_n;
    params.q = modulus;
    params.sigma = gaussian_stddev;
    params.name = "Custom";
    
    if (ring_dim_n < 128) {
        params.classical_bits = static_cast<int>(ring_dim_n * 0.5);
        params.quantum_bits = static_cast<int>(ring_dim_n * 0.25);
        params.is_secure = false;
    } else if (ring_dim_n < 256) {
        params.classical_bits = 80;
        params.quantum_bits = 40;
        params.is_secure = false;
    } else {
        params.classical_bits = static_cast<int>(ring_dim_n * 0.6);
        params.quantum_bits = static_cast<int>(ring_dim_n * 0.3);
        params.is_secure = (ring_dim_n >= 256);
    }
    
    return params;
}

void KEM::validateSecurityParameters() {
    double alpha = gaussian_stddev / modulus;
    
    Logger::log("\nValidating security parameters...");
    Logger::log("Ring dimension (n):     " + std::to_string(ring_dim_n));
    Logger::log("Modulus (q):            " + std::to_string(modulus));
    Logger::log("Gaussian σ:             " + std::to_string(gaussian_stddev));
    Logger::log("Noise ratio (α = σ/q):  " + std::to_string(alpha));
    
    if (ring_dim_n < 256) {
        Logger::log("⚠️  WARNING: Ring dimension n=" + std::to_string(ring_dim_n) + 
                   " is below recommended minimum of 256");
        Logger::log("   Current security: ~" + std::to_string(ring_dim_n * 0.5) + 
                   " bits (INSECURE)");
        Logger::log("   Recommended: n >= 256 for production use");
    }
    
    if (alpha > 0.01) {
        Logger::log("⚠️  WARNING: Large noise ratio α=" + std::to_string(alpha) + 
                   " may affect correctness");
    }
    
    if (validatePowerOfTwo(ring_dim_n)) {
        Logger::log("✓ Ring dimension is a power of 2 (required for efficiency)");
    }
    
    Logger::log("Parameter validation complete.\n");
}

void KEM::generateKeys() {
    Logger::log("\nGenerating keys...");
    a = sampleUniform();
    s = sampleGaussian(gaussian_stddev);
    
    Logger::log("Sampling gaussian polynomial e with σ=" + std::to_string(gaussian_stddev));
    Polynomial e = sampleGaussian(gaussian_stddev);
    
    Logger::log("Computing b = a*s + e");
    b = a * s + e;
    
    Logger::log("Public key a: " + a.toString());
    Logger::log("Public key b: " + b.toString());  
    Logger::log("Secret key s: " + s.toString());
}

Polynomial KEM::sampleUniform() {
    std::vector<uint64_t> coeffs(ring_dim_n);
    
    for (size_t i = 0; i < ring_dim_n; i++) {
        coeffs[i] = getRandomUint64() % modulus;
    }
    
    return Polynomial(coeffs, modulus);
}

Polynomial KEM::sampleGaussian(double stddev) {
    std::vector<uint64_t> coeffs(ring_dim_n);
    
    for (size_t i = 0; i < ring_dim_n; i++) {
        double sample = getRandomDouble() * stddev;
        int64_t rounded = static_cast<int64_t>(std::round(sample));
        
        if (rounded < 0) {
            rounded += modulus;
        }
        
        coeffs[i] = rounded;
    }
    
    return Polynomial(coeffs, modulus);
}

Polynomial KEM::messageToPolynomial(const std::vector<uint8_t>& message) {
    std::vector<uint64_t> coeffs(ring_dim_n, 0);
    
    size_t coeff_idx = 0;
    for (size_t byte_idx = 0; byte_idx < message.size() && coeff_idx < ring_dim_n; byte_idx++) {
        uint8_t byte = message[byte_idx];
        for (int j = 7; j >= 0 && coeff_idx < ring_dim_n; j--) {
            coeffs[coeff_idx++] = (byte >> j) & 1;
        }
    }
        
    return Polynomial(coeffs, modulus);
}

Polynomial KEM::hashToPolynomial(const std::vector<uint8_t>& message) {
    Logger::log("\nConverting message to polynomial using counter-based hashing");
    logMessageBytes("Input message", message);
    
    std::vector<uint64_t> coeffs(ring_dim_n, 0);
        
    size_t coeff_idx = 0;
    uint32_t counter = 0;
    
    while (coeff_idx < ring_dim_n) {
        std::vector<uint8_t> block;
        block.reserve(message.size() + sizeof(counter));
        const uint8_t* counter_bytes = reinterpret_cast<const uint8_t*>(&counter);
        block.insert(block.end(), counter_bytes, counter_bytes + sizeof(counter));
        block.insert(block.end(), message.begin(), message.end());
        Logger::log("Block " + std::to_string(counter) + " content:");
        logMessageBytes("  ", block);
        std::vector<uint8_t> hash = SHA256::hash(block);
        std::stringstream ss;
        ss << "Block " << counter << " hash: ";
        for (uint8_t b : hash) {
            ss << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(b);
        }
        Logger::log(ss.str());

        for (size_t byte_idx = 0; coeff_idx < ring_dim_n && byte_idx < hash.size(); byte_idx++) {
            for (int bit = 7; bit >= 0 && coeff_idx < ring_dim_n; bit--) {
                bool bit_value = (hash[byte_idx] >> bit) & 1;
                coeffs[coeff_idx++] = bit_value ? (modulus / 2) : 0;
            }
        }
        
        counter++;
    }
    
    Logger::log("Final polynomial coefficients:");
    std::stringstream result_ss;
    for (size_t i = 0; i < coeffs.size(); i++) {
        if (i > 0) result_ss << ", ";
        result_ss << coeffs[i];
    }
    Logger::log(result_ss.str());
    
    return Polynomial(coeffs, modulus);
}