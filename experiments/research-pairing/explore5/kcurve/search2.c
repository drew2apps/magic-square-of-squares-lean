// Sieve search for K=Q(i) points on E: Y^2 = (X-e1)(X-e2)(X-e3) through the descent form
//   general form: curve y^2 = x (x^2 + A1 x + A0), A0 = F1*F2, x = b U^2/W^2,
//   N := b (b^2 U^4 + A1 b U^2 W^2 + A0 W^4) must be a square in Z[i]
// U = a + b' i, W = c + d i,  0 <= a <= H, |b'| <= H, 1 <= c <= H, 0 <= d <= H (W up to units, U up to sign).
// A candidate is printed if N is a square or zero modulo every sieving place (degree-one primes of Z[i]); exact check is done elsewhere.
// usage: search H tag bre bim A1re A1im F1re F1im F2re F2im
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll; typedef __int128 i128;
#define MAXP 80
static int P[MAXP], R[MAXP], np = 0;
static unsigned char *T[MAXP], *okW[MAXP]; static int Tinf[MAXP];
static unsigned char *qr[MAXP];
static ll md(i128 x, ll p) { ll r = (ll)(x % p); return r < 0 ? r + p : r; }
static ll pw(ll a, ll e, ll p) { ll r = 1; a %= p; while (e) { if (e & 1) r = r*a % p; a = a*a % p; e >>= 1; } return r; }
int main(int argc, char **argv) {
  int H = atoi(argv[1]); const char *tag = argv[2];
  ll v[8]; for (int k = 0; k < 8; k++) v[k] = atoll(argv[3 + k]);
  for (int p = 13; np < MAXP - 1 && p < 700; p++) {
    if (p % 4 != 1) continue; int ok = 1; for (int q = 2; q*q <= p; q++) if (p % q == 0) ok = 0; if (!ok) continue;
    int r = 0; for (int x = 2; x < p; x++) if ((x*x + 1) % p == 0) { r = x; break; }
    for (int sgn = 0; sgn < 2; sgn++) { P[np] = p; R[np] = sgn ? p - r : r; np++; }
  }
  for (int k = 0; k < np; k++) {
    int p = P[k]; ll r = R[k];
    qr[k] = calloc(p, 1); for (int x = 0; x < p; x++) qr[k][(ll)x*x % p] = 1;
    ll b = md((i128)v[0] + (i128)v[1]*r, p), A1 = md((i128)v[2] + (i128)v[3]*r, p), F1 = md((i128)v[4] + (i128)v[5]*r, p), F2 = md((i128)v[6] + (i128)v[7]*r, p), A0 = F1*F2 % p;
    T[k] = malloc(p); okW[k] = malloc(p);
    for (int u = 0; u < p; u++) { ll u2 = (ll)u*u % p; ll bu = b*u2 % p; ll N = b*((bu*bu + A1*bu + A0) % p) % p; T[k][u] = qr[k][N]; }
    Tinf[k] = qr[k][b];
  }
  ll cnt = 0;
  for (int c = 1; c <= H; c++) for (int d = 0; d <= H; d++) {
    for (int k = 0; k < np; k++) {
      int p = P[k]; ll w = md((ll)c + (ll)d*R[k], p);
      if (w == 0) { memset(okW[k], Tinf[k], p); okW[k][0] = 1; }
      else { ll wi = pw(w, p - 2, p); for (int u = 0; u < p; u++) okW[k][u] = T[k][u*wi % p]; }
    }
    int p0 = P[0], r0 = R[0], p1 = P[1], r1 = R[1], p2 = P[2], r2 = R[2];
    for (int a = 0; a <= H; a++) {
      int u0 = (int)md((ll)a - (ll)H*r0, p0), u1 = (int)md((ll)a - (ll)H*r1, p1), u2 = (int)md((ll)a - (ll)H*r2, p2);
      for (int bb = -H; bb <= H; bb++) {
        if (okW[0][u0] && okW[1][u1] && okW[2][u2]) {
          int k; for (k = 3; k < np; k++) { int u = (int)md((ll)a + (ll)bb*R[k], P[k]); if (!okW[k][u]) break; }
          if (k == np && !(a == 0 && bb == 0)) { printf("HIT %s %d %d %d %d\n", tag, a, bb, c, d); fflush(stdout); cnt++; }
        }
        u0 += r0; if (u0 >= p0) u0 -= p0; u1 += r1; if (u1 >= p1) u1 -= p1; u2 += r2; if (u2 >= p2) u2 -= p2;
      }
    }
  }
  printf("DONE %s H=%d places=%d hits=%lld\n", tag, H, np, cnt);
  return 0;
}
