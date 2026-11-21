// Five-form system of the G1 family modulo an odd prime p:
//   F0 = m^2 b^2 + n^2 a^2, F8 = n^2 b^2 + m^2 a^2, G4 = 2(m^2+n^2)(a^2+b^2),
//   GL = G4 - 8 m n a b, GU = G4 + 8 m n a b     (G4, GL, GU are 4 F4, 4 L, 4 U)
// For every projective ratio (m:n) mod p, is there (a:b) with all five values squares (0 allowed)?
// Prints the ratios with no such (a:b).  usage: classify LIMIT
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
int main(int argc, char **argv) {
  int LIM = atoi(argv[1]);
  char *sieve = calloc(LIM + 1, 1);
  char *sq = malloc(LIM + 1);
  long total = 0;
  for (int p = 3; p <= LIM; p++) {
    if (sieve[p]) continue;
    for (long q = (long)p * p; q <= LIM; q += p) sieve[q] = 1;
    if (p % 2 == 0) continue;
    memset(sq, 0, p);
    for (u64 x = 0; x < (u64)p; x++) sq[(x * x) % p] = 1;
    int cnt = 0;
    // ratios: (m:1) for m = 0..p-1, and (1:0) coded as m = p
    for (int r = 0; r <= p; r++) {
      u64 m = (r < p) ? r : 1, n = (r < p) ? 1 : 0;
      u64 m2 = m * m % p, n2 = n * n % p, s2 = 2 * (m2 + n2) % p, mn8 = 8 * m % p * n % p;
      int found = 0;
      // (a:b) = (t:1), t = 0..p-1, then (1:0)
      for (int t = 0; t <= p && !found; t++) {
        u64 a = (t < p) ? t : 1, b = (t < p) ? 1 : 0;
        u64 a2 = a * a % p, b2 = b * b % p, ab = a * b % p;
        u64 f0 = (m2 * b2 + n2 * a2) % p; if (!sq[f0]) continue;
        u64 f8 = (n2 * b2 + m2 * a2) % p; if (!sq[f8]) continue;
        u64 g4 = s2 * ((a2 + b2) % p) % p; if (!sq[g4]) continue;
        u64 x = mn8 * ab % p;
        if (!sq[(g4 + p - x) % p]) continue;
        if (!sq[(g4 + x) % p]) continue;
        found = 1;
      }
      if (!found) { printf("%d %d\n", p, r); cnt++; }
    }
    total += cnt;
  }
  fprintf(stderr, "total obstructed classes: %ld\n", total);
  return 0;
}
