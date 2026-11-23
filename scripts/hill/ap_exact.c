/* The weaker problem of HILL_NOTE.md, section 5, searched without floating point.

   w = (p^2+q^2)^2 / (4pq(p^2-q^2)) for coprime p > q >= 1 of different parity, p <= P.  Question: are there three distinct
   values with w1 + w3 = 2 w2 ?

   Every value is stored as a reduced fraction a/b and as its residue r = a * b^(-1) modulo the prime M = 2^61 - 1 (b < M,
   so b is invertible).  Reduction modulo M is a ring homomorphism on the fractions with denominator prime to M, so
       w1 + w3 = 2 w2   implies   r1 + r3 = 2 r2  (mod M).
   The program checks that the residues of the distinct values are distinct, puts the residues 2 r (mod M) in a hash table,
   and looks up r1 + r3 for every pair.  So no solution can be missed; a hit is a candidate, printed with its fractions (they can then be
   added by hand or with any exact arithmetic) and tested again modulo two other primes.  No candidate has occurred.
   Three made-up fractions in arithmetic progression are added to the list; the program must find them.

   gcc -O2 -fopenmp -o ap_exact ap_exact.c ;  ./ap_exact P                                                             */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned __int128 u128; typedef uint64_t u64;
static const u64 M = (1ULL << 61) - 1, M2 = 4611686018427387847ULL, M3 = 4611686018427387751ULL;   /* three primes below 2^62 */
static u64 mulmod(u64 x, u64 y, u64 m) { return (u64)((u128)x * y % m); }
static u64 powmod(u64 x, u64 e, u64 m) { u64 r = 1; x %= m; while (e) { if (e & 1) r = mulmod(r, x, m); x = mulmod(x, x, m); e >>= 1; } return r; }
static u64 res(u64 a, u64 b, u64 m) { return mulmod(a % m, powmod(b % m, m - 2, m), m); }
static u64 gcd(u64 a, u64 b) { while (b) { u64 t = a % b; a = b; b = t; } return a; }
typedef struct { u64 a, b, r; int p, q; } val;
static int cmpfrac(const void *x, const void *y) { const val *u = x, *v = y; u128 l = (u128)u->a * v->b, r = (u128)v->a * u->b; return l < r ? -1 : l > r; }
static int cmpres(const void *x, const void *y) { u64 u = *(const u64 *)x, v = *(const u64 *)y; return u < v ? -1 : u > v; }
int main(int argc, char **argv) {
    if (argc < 2) return 2;
    int P = atoi(argv[1]); size_t cap = (size_t)P * P / 4 + 16, n = 0; val *v = malloc(cap * sizeof(val));
    for (int p = 2; p <= P; p++) for (int q = 1; q < p; q++) if (((p - q) & 1) && gcd(p, q) == 1) {
        u64 y = (u64)p * p + (u64)q * q, a = y * y, b = 4 * (u64)p * q * ((u64)p * p - (u64)q * q), g = gcd(a, b);
        v[n].a = a / g; v[n].b = b / g; v[n].p = p; v[n].q = q; n++; }
    u64 pa = 1000003, pb = 7919, da = 104729, db = 15485863;                 /* planted: pa/pb + k da/db, k = 0, 1, 2 */
    for (int k = 0; k < 3; k++) { u64 a = pa * db + k * da * pb, b = pb * db, g = gcd(a, b); v[n].a = a / g; v[n].b = b / g; v[n].p = 0; v[n].q = k; n++; }
    qsort(v, n, sizeof(val), cmpfrac);
    size_t m = 0; for (size_t i = 0; i < n; i++) if (m == 0 || cmpfrac(&v[m - 1], &v[i]) != 0) v[m++] = v[i];      /* exact duplicates removed */
    n = m; u64 *rs = malloc(n * sizeof(u64));
    for (size_t i = 0; i < n; i++) { v[i].r = res(v[i].a, v[i].b, M); rs[i] = v[i].r; }
    qsort(rs, n, sizeof(u64), cmpres); for (size_t i = 1; i < n; i++) if (rs[i] == rs[i - 1]) { printf("two distinct values with the same residue: stop\n"); return 1; }
    int hb = 1; while (((size_t)1 << hb) < 4 * n) hb++;
    size_t hs = (size_t)1 << hb; u64 *hk = malloc(hs * sizeof(u64)); uint32_t *hv = malloc(hs * sizeof(uint32_t));
    for (size_t i = 0; i < hs; i++) hk[i] = ~0ULL;
    int fb = 27; unsigned char *bits = calloc((size_t)1 << (fb - 3), 1);                                                   /* a first filter of 2^27 bits */
    for (size_t i = 0; i < n; i++) { u64 t = 2 * v[i].r; if (t >= M) t -= M; size_t h = (size_t)((t * 0x9E3779B97F4A7C15ULL) >> (64 - hb));
        while (hk[h] != ~0ULL) h = (h + 1) & (hs - 1);
        hk[h] = t; hv[h] = (uint32_t)i;
        size_t f = (size_t)((t * 0xD6E8FEB86659FD93ULL) >> (64 - fb)); bits[f >> 3] |= 1 << (f & 7); }
    long cand = 0, planted = 0;
    #pragma omp parallel for schedule(dynamic, 256) reduction(+:cand, planted)
    for (size_t i = 0; i < n; i++) { u64 ri = v[i].r;
        for (size_t k = i + 1; k < n; k++) { u64 s = ri + v[k].r; if (s >= M) s -= M;
            size_t f = (size_t)((s * 0xD6E8FEB86659FD93ULL) >> (64 - fb)); if (!(bits[f >> 3] >> (f & 7) & 1)) continue;
            size_t h = (size_t)((s * 0x9E3779B97F4A7C15ULL) >> (64 - hb));
            while (hk[h] != ~0ULL && hk[h] != s) h = (h + 1) & (hs - 1);
            if (hk[h] != s) continue;
            size_t j = hv[h]; int more = 1;
            { u64 x = (res(v[i].a, v[i].b, M2) + res(v[k].a, v[k].b, M2)) % M2, y = 2 * res(v[j].a, v[j].b, M2) % M2; if (x != y) more = 0; }
            { u64 x = (res(v[i].a, v[i].b, M3) + res(v[k].a, v[k].b, M3)) % M3, y = 2 * res(v[j].a, v[j].b, M3) % M3; if (x != y) more = 0; }
            if (v[i].p == 0 && v[j].p == 0 && v[k].p == 0) { planted += more; continue; }
            #pragma omp critical
            { cand++; printf("candidate%s: (%d,%d) %llu/%llu  (%d,%d) %llu/%llu  (%d,%d) %llu/%llu\n", more ? " (also modulo two other primes)" : " (not modulo the other primes)",
                v[i].p, v[i].q, (unsigned long long)v[i].a, (unsigned long long)v[i].b, v[j].p, v[j].q, (unsigned long long)v[j].a, (unsigned long long)v[j].b,
                v[k].p, v[k].q, (unsigned long long)v[k].a, (unsigned long long)v[k].b); } } }
    printf("shapes with q < p <= %d: %zu distinct values of w, %zu pairs, each looked up exactly modulo 2^61 - 1\n", P, n - 3, (n - 3) * (n - 4) / 2);
    printf("self-test, the planted progression is found: %s\n", planted == 1 ? "yes" : "NO");
    printf("candidates (to be confirmed in exact arithmetic): %ld\n", cand);
    return planted == 1 ? 0 : 1; }
