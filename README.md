# RLWE Blind DH Exchange

This project implements a Ring Learning With Errors (RLWE) based blind exchange protocol, designed as an alternative to the current BDHKE (Blind Diffie-Hellman Key Exchange) used in Cashu. This implementation replaces elliptic curve primitives with lattice-based cryptography primitives that are believed to be quantum-resistant.

## ⚠️ Security Notice

**Default parameters are now set to KYBER512 (NIST standard)** providing ~128 bits of classical security.

The previous test parameters (n=8, n=32) are **INSECURE** and only suitable for development/testing. They are still available via `SecurityLevel::TEST_TINY` and `SecurityLevel::TEST_SMALL` for fast unit tests.

### Security Levels

| Level | n | q | σ | Classical | Quantum | Status |
|-------|---|---|---|-----------|---------|---------|
| **KYBER512** (default) | 256 | 3329 | 1.6 | ~128 bits | ~64 bits | ✅ SECURE |
| **MODERATE** | 512 | 12289 | 3.2 | ~192 bits | ~96 bits | ✅ SECURE |
| **HIGH** | 1024 | 16384 | 3.2 | ~256 bits | ~128 bits | ✅ SECURE |
| TEST_TINY | 8 | 7681 | 3.0 | ~4 bits | ~2 bits | ⚠️ INSECURE |
| TEST_SMALL | 32 | 7681 | 3.0 | ~16 bits | ~8 bits | ⚠️ INSECURE |

See [SECURITY_ANALYSIS.md](SECURITY_ANALYSIS.md) for detailed security analysis.

## Quick Start

### Using Secure Defaults (Recommended)

```cpp
#include <rlwe.h>

// Create instance with KYBER512 parameters (default, secure)
RLWESignature rlwe;  // or explicitly: RLWESignature rlwe(SecurityLevel::KYBER512);

// Generate keys
rlwe.generateKeys();
auto [a, b] = rlwe.getPublicKey();

// Client: Blind a secret
std::vector<uint8_t> secret = {0xDE, 0xAD, 0xBE, 0xEF};
auto [blindedMessage, blindingFactor] = rlwe.computeBlindedMessage(secret);

// Server: Sign blinded message
Polynomial blindSignature = rlwe.blindSign(blindedMessage);

// Client: Unblind signature
Polynomial signature = rlwe.computeSignature(blindSignature, blindingFactor, b);

// Server: Verify
bool verified = rlwe.verify(secret, signature);  // Should be true
```

### Choosing Security Level

```cpp
// For different security levels:
RLWESignature rlwe_kyber(SecurityLevel::KYBER512);    // 128-bit security (NIST standard)
RLWESignature rlwe_moderate(SecurityLevel::MODERATE);  // 192-bit security
RLWESignature rlwe_high(SecurityLevel::HIGH);          // 256-bit security

// For testing only (INSECURE):
RLWESignature rlwe_test(SecurityLevel::TEST_TINY);     // Fast but insecure
```

### Custom Parameters

```cpp
// Advanced: Use custom parameters (not recommended unless you know what you're doing)
RLWESignature rlwe(512, 12289, 3.2);  // n, q, sigma
```

## Overview

The blind exchange protocol works in the polynomial ring $R = Z[x]/(x^n + 1)$, where $n$ is a power of 2 and replaces key Cashu operations with RLWE equivalents:

### Key Transformations

- **Point-Scalar Multiplication → Noisy Polynomial Multiplication**
  - Instead of ECC point-scalar multiplication, we use RLWE's fundamental operation: `a * s + e`
  - Where:
    - `a` and `s` are polynomials
    - `e` is a small "error" polynomial sampled from a Gaussian distribution
    
- **Point Addition → Polynomial Addition**
  - ECC point addition is replaced by simple polynomial addition in the ring
  - Operations are performed modulo both the ring polynomial $(x^n + 1)$ and a prime modulus q

### Key Function Replacements

- **hashToCurve → hashToPoly**
  - Cashu's `hashToCurve` function is replaced with `hashToPoly`
  - `hashToPoly` takes a secret message and maps it to a polynomial using SHA256
  - The resulting polynomial maintains the security properties needed for the blind exchange

### Security Basis

The security of this implementation relies on the Ring Learning With Errors problem, which is considered to be:
- Quantum-resistant (unlike ECC-based systems)
- Shown to be as hard as solving worst-case lattice problems

## Warning

This is an experimental implementation meant for research and learning purposes. It should not be used in production environments without thorough security review and analysis.
