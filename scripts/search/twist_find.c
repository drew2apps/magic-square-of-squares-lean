// Non-degenerate t0 = p/q for the ratio m:n with F0*F8 and GL*GU both squares:
//   (m^2 q^2 + n^2 p^2)(n^2 q^2 + m^2 p^2) = square,   (S T)^2 - (4 m n p q)^2 = square,  S = m^2+n^2, T = p^2+q^2.
// Then the four points above t0, -t0, 1/t0, -1/t0 are rational on the twist d y^2 = F0 GL, d = squarefree part of F0 GL (t0).
// Output: one line "m n p q" per solution with p != q (none is known).  Result of the run to 400: twist_search.txt.
// usage: twist_find MAXM MAXP [m_lo m_hi]
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef unsigned __int128 u128; typedef unsigned long long u64;
static u64 gcd(u64 a, u64 b){ while(b){ u64 t=a%b; a=b; b=t;} return a; }
static int issq(u128 v){ if(v==0) return 1; long double r = sqrtl((long double)v); u128 s=(u128)r; while(s*s>v) s--; while((s+1)*(s+1)<=v) s++; return s*s==v; }
int main(int argc, char**argv){
  u64 MAXM=(u64)atol(argv[1]), MAXP=(u64)atol(argv[2]); u64 mlo = argc>3?(u64)atol(argv[3]):2, mhi = argc>4?(u64)atol(argv[4]):MAXM;
  static unsigned char q64[64], q63[63], q65[65];
  for(int i=0;i<64;i++) q64[(i*i)%64]=1;
  for(int i=0;i<63;i++) q63[(i*i)%63]=1;
  for(int i=0;i<65;i++) q65[(i*i)%65]=1;
  for(u64 m=mlo;m<=mhi;m++) for(u64 n=1;n<m;n++){ if(gcd(m,n)!=1) continue; u64 S=m*m+n*n;
    for(u64 p=1;p<=MAXP;p++) for(u64 q=1;q<=MAXP;q++){ if(p==q||gcd(p,q)!=1) continue;
      u64 T=p*p+q*q; u128 a=(u128)S*T, b=(u128)4*m*n*p*q; u128 lo=a-b, hi=a+b;
      // filter mod 64, 63, 65 on the product lo*hi
      u64 l64=(u64)(lo%64)*(u64)(hi%64)%64; if(!q64[l64]) continue;
      u64 l63=(u64)(lo%63)*(u64)(hi%63)%63; if(!q63[l63]) continue;
      u64 l65=(u64)(lo%65)*(u64)(hi%65)%65; if(!q65[l65]) continue;
      // lo*hi may exceed 128 bits? S,T < 2*MAX^2 each: a < 4e12 for MAX=1000, product < 1.6e25 fits
      if(!issq(lo*hi)) continue;
      u128 x=(u128)(m*q)*(m*q)+(u128)(n*p)*(n*p), y=(u128)(n*q)*(n*q)+(u128)(m*p)*(m*p);
      if(!issq(x*y)) continue;
      printf("%llu %llu %llu %llu\n", m,n,p,q); fflush(stdout);
    } }
  return 0; }
