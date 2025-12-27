#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <sstream>
#include <logging.h>

/**
 * @brief Represents a polynomial in the quotient ring Z_q[x]/(x^n + 1).
 *
 * The polynomial is stored as a vector of coefficients
 * \f$(c_0, c_1, \ldots, c_{n-1})\f$ corresponding to the element
 * \f$c_0 + c_1 x + \cdots + c_{n-1} x^{n-1}\f$ in the ring
 * \f$\mathbb{Z}_q[x]/(x^n + 1)\f$.
 *
 * All coefficients are maintained in the canonical range
 * \f$[0, q-1]\f$ and arithmetic is performed modulo both \f$q\f$
 * and the polynomial modulus \f$x^n + 1\f$ where applicable.
 */
class Polynomial {
public:
    /**
     * @brief Construct a zero polynomial in Z_q[x]/(x^n + 1).
     *
     * @param n Ring dimension (number of coefficients).
     * @param q Coefficient modulus.
     */
    Polynomial(size_t n, uint64_t q)
        : coeffs(n, 0), ring_dim(n), modulus(q) {
        Logger::log("Created zero polynomial of degree " + std::to_string(n - 1) +
                    " with modulus " + std::to_string(q));
    }

    /**
     * @brief Construct a polynomial from a vector of coefficients.
     *
     * The coefficients are interpreted as elements of Z_q and reduced
     * modulo @p q when set via setCoefficients().
     *
     * @param coefficients Coefficient vector in ascending degree order.
     * @param q Coefficient modulus.
     */
    Polynomial(const std::vector<uint64_t>& coefficients, uint64_t q)
        : coeffs(coefficients), ring_dim(coefficients.size()), modulus(q) {
        Logger::log("Created polynomial from coefficients: " +
                    Logger::vectorToString(coefficients) +
                    " with modulus " + std::to_string(q));
    }

    /**
     * @brief Access a coefficient by index.
     *
     * @param idx Zero-based coefficient index.
     * @return Reference to the coefficient at @p idx.
     *
     * @note No bounds checking is performed.
     */
    uint64_t& operator[](size_t idx) {
        return coeffs[idx];
    }

    /**
     * @brief Access a coefficient by index (const overload).
     *
     * @param idx Zero-based coefficient index.
     * @return Const reference to the coefficient at @p idx.
     *
     * @note No bounds checking is performed.
     */
    const uint64_t& operator[](size_t idx) const {
        return coeffs[idx];
    }

    /**
     * @brief Get the ring dimension.
     *
     * This is equal to the number of coefficients and the degree bound
     * @f$n@f$ in @f$Z_q[x]/(x^n + 1)@f$.
     *
     * @return Ring dimension.
     */
    size_t degree() const {
        return ring_dim;
    }

    /**
     * @brief Get the coefficient modulus.
     *
     * @return Modulus @f$q@f$.
     */
    uint64_t getModulus() const {
        return modulus;
    }

    /**
     * @brief Add two polynomials coefficient-wise modulo the common modulus.
     *
     * @param other Polynomial to add.
     * @return Sum @f*this + other@f$ modulo @f$q@f$.
     *
     * @throws std::invalid_argument If the ring dimension or modulus does not match.
     */
    Polynomial operator+(const Polynomial& other) const;

    /**
     * @brief Subtract another polynomial coefficient-wise modulo the common modulus.
     *
     * @param other Polynomial to subtract.
     * @return Difference @f$this - other@f$ modulo @f$q@f$.
     *
     * @throws std::invalid_argument If the ring dimension or modulus does not match.
     */
    Polynomial operator-(const Polynomial& other) const;

    /**
     * @brief Negate the polynomial modulo the coefficient modulus.
     *
     * @return Polynomial whose coefficients are @f$-c_i \bmod q@f$.
     */
    Polynomial operator-() const;

    /**
     * @brief Multiply two polynomials in Z_q[x]/(x^n + 1).
     *
     * The product is reduced modulo both @f$q@f$ and the polynomial
     * modulus @f$x^n + 1@f$.
     *
     * @param other Polynomial to multiply by.
     * @return Product @f$this \cdot other@f$ in the quotient ring.
     *
     * @throws std::invalid_argument If the ring dimension or modulus does not match.
     */
    Polynomial operator*(const Polynomial& other) const;

    /**
     * @brief Multiply the polynomial by a scalar modulo the coefficient modulus.
     *
     * @param scalar Multiplier in Z_q.
     * @return Scaled polynomial with each coefficient multiplied by @p scalar
     *         modulo @f$q@f$.
     */
    Polynomial operator*(uint64_t scalar) const;

    /**
     * @brief Get a const reference to the internal coefficient vector.
     *
     * @return Coefficient vector in ascending degree order.
     */
    const std::vector<uint64_t>& getCoeffs() const {
        return coeffs;
    }

    /**
     * @brief Map coefficients to a binary "signal" representation.
     *
     * Each coefficient is rounded to either @f$0@f$ or @f$q/2@f$ depending on
     * which value is closer in the cyclic group @f$\mathbb{Z}_q@f$.
     *
     * This is used by the signature verification procedure to obtain a
     * coarse, noise-tolerant representation of a polynomial.
     *
     * @return Polynomial with coefficients in the set {0, q/2}.
     */
    Polynomial polySignal() const;

    /**
     * @brief Replace the polynomial coefficients.
     *
     * The input vector must match the current ring dimension. Each value
     * is reduced modulo the coefficient modulus @f$q@f$.
     *
     * @param new_coeffs New coefficient vector.
     *
     * @throws std::invalid_argument If @p new_coeffs has a different size
     *         than the ring dimension.
     */
    void setCoefficients(const std::vector<uint64_t>& new_coeffs);

    /**
     * @brief Convert the polynomial to a human-readable string.
     *
     * The string contains the ring dimension, modulus, and the coefficient
     * vector, primarily for logging and debugging.
     *
     * @return Descriptive string representation.
     */
    std::string toString() const;

    /**
     * @brief Serialize the polynomial to a byte vector.
     *
     * The byte encoding has the following layout:
     * - ring dimension (size_t)
     * - modulus (uint64_t)
     * - coefficients (uint64_t[ring_dim])
     *
     * All values are written in the native endianness of the host.
     *
     * @return Contiguous byte representation of the polynomial.
     */
    std::vector<uint8_t> toBytes() const;

private:
    /**
     * @brief Coefficient storage in ascending degree order.
     */
    std::vector<uint64_t> coeffs;

    /**
     * @brief Polynomial ring dimension (number of coefficients).
     */
    size_t ring_dim;

    /**
     * @brief Coefficient modulus.
     */
    uint64_t modulus;

    /**
     * @brief Reduce an integer modulo a positive modulus.
     *
     * @param x Value to reduce.
     * @param m Positive modulus.
     * @return Result in the canonical range @f$[0, m-1]@f$.
     */
    static uint64_t mod(int64_t x, uint64_t m) {
        int64_t r = x % static_cast<int64_t>(m);
        return r < 0 ? r + m : r;
    }
};

#endif // POLYNOMIAL_H
