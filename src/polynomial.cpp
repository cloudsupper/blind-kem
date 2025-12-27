#include <polynomial.h>
#include <stdexcept>
#include <ntt.h>

Polynomial Polynomial::polySignal() const {
    Polynomial result(ring_dim, modulus);
    uint64_t half_mod = modulus / 2;
    
    for (size_t i = 0; i < ring_dim; i++) {
        uint64_t coeff = coeffs[i];
        uint64_t dist_to_zero = std::min(coeff, modulus - coeff);
        uint64_t dist_to_half = std::min(
            (coeff >= half_mod) ? coeff - half_mod : half_mod - coeff,
            (coeff >= half_mod) ? modulus - coeff + half_mod : modulus - half_mod + coeff
        );
        
        result[i] = (dist_to_zero <= dist_to_half) ? 0 : half_mod;
    }
    
    Logger::log("Rounded polynomial coefficients to binary signal");
    return result;
}

Polynomial Polynomial::operator+(const Polynomial& other) const {
    if (ring_dim != other.ring_dim || modulus != other.modulus) {
        throw std::invalid_argument("Polynomials must be in the same ring");
    }

    Logger::log("Adding polynomials:\n  " + toString() + "\n  " + other.toString());

    Polynomial result(ring_dim, modulus);
    for (size_t i = 0; i < ring_dim; i++) {
        result[i] = (coeffs[i] + other.coeffs[i]) % modulus;
    }

    Logger::log("Addition result:\n  " + result.toString());
    return result;
}

Polynomial Polynomial::operator-(const Polynomial& other) const {
    if (ring_dim != other.ring_dim || modulus != other.modulus) {
        throw std::invalid_argument("Polynomials must be in the same ring");
    }

    Logger::log("Subtracting polynomials:\n  " + toString() + "\n  " + other.toString());

    Polynomial result(ring_dim, modulus);
    for (size_t i = 0; i < ring_dim; i++) {
        result[i] = mod(static_cast<int64_t>(coeffs[i]) - 
                       static_cast<int64_t>(other.coeffs[i]), modulus);
    }

    Logger::log("Subtraction result:\n  " + result.toString());
    return result;
}

Polynomial Polynomial::operator-() const {
    Logger::log("Negating polynomial:\n  " + toString());

    Polynomial result(ring_dim, modulus);
    for (size_t i = 0; i < ring_dim; i++) {
        result[i] = (coeffs[i] == 0) ? 0 : modulus - coeffs[i];
    }

    Logger::log("Negation result:\n  " + result.toString());
    return result;
}

Polynomial Polynomial::operator*(const Polynomial& other) const {
    if (ring_dim != other.ring_dim || modulus != other.modulus) {
        throw std::invalid_argument("Polynomials must be in the same ring");
    }

    Logger::log("Multiplying polynomials (NTT-accelerated where available):\n  " +
                toString() + "\n  " + other.toString());


    NTT ntt(ring_dim, modulus, /*negacyclic=*/true);

    std::vector<std::uint64_t> a_vec = coeffs;
    std::vector<std::uint64_t> b_vec = other.coeffs;

    ntt.forward(a_vec);
    ntt.forward(b_vec);

    for (std::size_t i = 0; i < ring_dim; ++i) {
        a_vec[i] = (a_vec[i] * b_vec[i]) % modulus;
    }

    ntt.inverse(a_vec);

    Polynomial result(ring_dim, modulus);
    result.setCoefficients(a_vec);

    Logger::log("NTT-based multiplication result:\n  " + result.toString());
    return result;
}

Polynomial Polynomial::operator*(uint64_t scalar) const {
    Logger::log("Multiplying polynomial by scalar " + std::to_string(scalar) + ":\n  " + toString());

    Polynomial result(ring_dim, modulus);
    for (size_t i = 0; i < ring_dim; i++) {
        result[i] = (coeffs[i] * scalar) % modulus;
    }

    Logger::log("Scalar multiplication result:\n  " + result.toString());
    return result;
}

void Polynomial::setCoefficients(const std::vector<uint64_t>& new_coeffs) {
    if (new_coeffs.size() != ring_dim) {
        throw std::invalid_argument("New coefficient vector size must match polynomial ring dimension");
    }
    coeffs = new_coeffs;
    for (auto& c : coeffs) {
        c = mod(c, modulus);
    }
    Logger::log("Updated polynomial coefficients to: " + Logger::vectorToString(coeffs));
}

std::string Polynomial::toString() const {
    std::stringstream ss;
    ss << "Polynomial(dim=" << ring_dim << ", q=" << modulus << "): ";
    ss << Logger::vectorToString(coeffs);
    return ss.str();
}

std::vector<uint8_t> Polynomial::toBytes() const {
    std::vector<uint8_t> bytes;
    bytes.reserve(sizeof(size_t) + sizeof(uint64_t) + coeffs.size() * sizeof(uint64_t));

    const uint8_t* dim_bytes = reinterpret_cast<const uint8_t*>(&ring_dim);
    bytes.insert(bytes.end(), dim_bytes, dim_bytes + sizeof(size_t));

    const uint8_t* mod_bytes = reinterpret_cast<const uint8_t*>(&modulus);
    bytes.insert(bytes.end(), mod_bytes, mod_bytes + sizeof(uint64_t));

    for (const uint64_t& coeff : coeffs) {
        const uint8_t* coeff_bytes = reinterpret_cast<const uint8_t*>(&coeff);
        bytes.insert(bytes.end(), coeff_bytes, coeff_bytes + sizeof(uint64_t));
    }

    return bytes;
}