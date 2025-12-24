#include <bits/stdc++.h>
using namespace std;

static uint64_t modMul(uint64_t a, uint64_t b, uint64_t m) {
    return (unsigned __int128)a * b % m;
}

static uint64_t modPow(uint64_t base, uint64_t exp, uint64_t m) {
    uint64_t res = 1;
    while (exp) {
        if (exp & 1) res = modMul(res, base, m);
        base = modMul(base, base, m);
        exp >>= 1;
    }
    return res;
}

static uint64_t modInverse(uint64_t a, uint64_t m) {
    long long t = 0, new_t = 1;
    long long r = (long long)m;
    long long new_r = (long long)(a % m);
    while (new_r != 0) {
        long long q = r / new_r;
        long long tmp_t = t - q * new_t;
        t = new_t;
        new_t = tmp_t;
        long long tmp_r = r - q * new_r;
        r = new_r;
        new_r = tmp_r;
    }
    if (r > 1) throw runtime_error("no inverse");
    if (t < 0) t += (long long)m;
    return (uint64_t)t;
}

static uint64_t findPsi(uint64_t q, uint64_t n) {
    uint64_t k = 2 * n;
    uint64_t order_factor = (q - 1) / k;
    for (uint64_t g = 2; g < q; ++g) {
        uint64_t cand = modPow(g, order_factor, q); // order divides k
        if (modPow(cand, k, q) != 1) continue;
        // ensure no smaller power gives 1, and cand^{k/2} = -1
        bool ok = true;
        uint64_t tmp_k = k;
        while (tmp_k % 2 == 0 && tmp_k > 1) {
            tmp_k /= 2;
            if (modPow(cand, tmp_k, q) == 1) { ok = false; break; }
        }
        if (!ok) continue;
        if (modPow(cand, k/2, q) != q-1) continue;
        return cand;
    }
    return 0;
}

int main() {
    vector<pair<uint64_t,uint64_t>> params = {
        {8, 7681},
        {32, 7681},
        {256, 7681},
        {512, 12289},
        {1024, 18433},
    };

    cout << "#ifndef NTT_TABLES_H\n#define NTT_TABLES_H\n\n";
    cout << "#include <cstddef>\n#include <cstdint>\n\n";
    cout << "namespace ntt_tables {\n\n";
    cout << "struct PsiTables {\\n    std::size_t n;\\n    std::uint64_t q;\\n    std::uint64_t psi;\\n    std::uint64_t psi_inv;\\n    const std::uint64_t* twist;\\n    const std::uint64_t* twist_inv;\\n};\\n\\n";

    for (auto [n,q] : params) {
        uint64_t k = 2*n;
        if ((q-1) % k != 0) {
            cerr << "(q-1) not divisible by 2n for n=" << n << " q=" << q << "\n";
            continue;
        }
        uint64_t psi = findPsi(q, n);
        if (!psi) {
            cerr << "Failed to find psi for n=" << n << " q=" << q << "\n";
            continue;
        }
        uint64_t psi_inv = modInverse(psi, q);

        string tag = to_string(n) + string("_") + to_string(q);

        cout << "inline constexpr std::uint64_t psi_" << tag << " = " << psi << "ULL;\n";
        cout << "inline constexpr std::uint64_t psi_inv_" << tag << " = " << psi_inv << "ULL;\n";

        // twist = psi^{2i+1}
        cout << "inline constexpr std::uint64_t twist_" << tag << "[" << n << "] = {";
        for (uint64_t i = 0; i < n; ++i) {
            uint64_t e = 2*i + 1;
            uint64_t v = modPow(psi, e, q);
            cout << v << "ULL";
            if (i + 1 != n) cout << ",";
        }
        cout << "};\n";

        // twist_inv = psi^{-(2i+1)} = (psi_inv)^{2i+1}
        cout << "inline constexpr std::uint64_t twist_inv_" << tag << "[" << n << "] = {";
        for (uint64_t i = 0; i < n; ++i) {
            uint64_t e = 2*i + 1;
            uint64_t v = modPow(psi_inv, e, q);
            cout << v << "ULL";
            if (i + 1 != n) cout << ",";
        }
        cout << "};\n";

        cout << "inline constexpr PsiTables tables_" << tag << "{" << n << ", " << q
             << ", psi_" << tag << ", psi_inv_" << tag
             << ", twist_" << tag << ", twist_inv_" << tag << "};\n\n";
    }

    cout << "inline constexpr const PsiTables* getPsiTables(std::size_t n, std::uint64_t q) {\n";
    for (auto [n,q] : params) {
        string tag = to_string(n) + string("_") + to_string(q);
        cout << "    if (n == " << n << " && q == " << q << ") return &tables_" << tag << ";\n";
    }
    cout << "    return nullptr;\n}" << "\n\n";

    cout << "} // namespace ntt_tables\n\n#endif // NTT_TABLES_H\n";
}
