// Pairs of ratios (m : n), (m' : n') of height at most H, and the necessary conditions of G1Centre.lean.
//
// A primitive 3x3 magic square of squares with b + d, h + f, f - h > 0 is given by two ratios in lowest terms,
// m : n = (b + d) : (h + f) and m' : n' = (b + d) : (f - h), through
//   b = mm' + nn', d = mm' - nn', f = nm' + mn', h = nm' - mn', 2 e^2 = (m^2 + n^2)(m'^2 + n'^2),
//   a^2 = n^2 m'^2 + m^2 n'^2, i^2 = m^2 m'^2 + n^2 n'^2, c^2 = e^2 - 2mnm'n', g^2 = e^2 + 2mnm'n'.
// This program counts the ordered pairs with 1 <= n < m <= H coprime and 1 <= m', n' <= H coprime, and how
// many pass, cumulatively:
//   F1  2 (m^2 + n^2)(m'^2 + n'^2) is a square (IsMSS.two_mul_norms_isSquare; the centre is a square);
//   F2  every prime dividing m^2 + n^2 to an odd power is 1 mod 8 (IsMSS.norm_isSquare_or_prime_dvd_center);
//   F3  neither ratio lies in one of the 484 classes of modp_classes.txt (G1.fullTable_valid);
//   F4  the four corner entries are squares as well;
//   F5  the nine entries are distinct and nonzero.
// F1 is found by grouping the ratios by the odd squarefree part of m^2 + n^2 and its parity, not by testing
// every pair. Usage: ratio_pairs H modp_classes.txt   (H <= 6000)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef unsigned long long u64; typedef __int128 i128;
typedef struct { u64 c; unsigned m, n; unsigned char even, f2, f3; } R;
static unsigned *spf; static int ncls; static unsigned clsp[600], clsr[600];
static unsigned gcdu(unsigned a, unsigned b){ while(b){unsigned t=a%b;a=b;b=t;} return a; }
static int cmp(const void*x,const void*y){ const R*a=x,*b=y; if(a->c!=b->c) return a->c<b->c?-1:1; return (int)a->even-(int)b->even; }
static int issq(i128 v, u64 *r){ if(v<0) return 0; u64 s=(u64)sqrtl((long double)v); while((i128)s*s>v) s--; while((i128)(s+1)*(s+1)<=v) s++; *r=s; return (i128)s*s==v; }
static void fill(R *r, unsigned m, unsigned n){
  u64 S=(u64)m*m+(u64)n*n; r->m=m; r->n=n; r->even=!(S&1); r->f2=1; r->f3=1;
  u64 x=S, c=1; while(!(x&1)) x>>=1;
  while(x>1){ unsigned p=spf[x]; int e=0; while(x%p==0){x/=p;e++;} if(e&1){ c*=p; if(p%8!=1) r->f2=0; } }
  r->c=c;
  for(int k=0;k<ncls;k++){ unsigned p=clsp[k]; if(n%p==0) continue; if(m%p==(u64)clsr[k]*(n%p)%p){ r->f3=0; break; } }
}
int main(int argc,char**argv){
  if(argc!=3){ fprintf(stderr,"usage: ratio_pairs H modp_classes.txt\n"); return 2; }
  unsigned H=atoi(argv[1]); if(H<2||H>6000) return 2;
  FILE*fp=fopen(argv[2],"r"); if(!fp) return 2; char line[256];
  while(fgets(line,sizeof line,fp)){ unsigned p,r; if(sscanf(line,"%u %u",&p,&r)==2 && ncls<600){ clsp[ncls]=p; clsr[ncls]=r; ncls++; } }
  fclose(fp); if(ncls!=484){ fprintf(stderr,"expected 484 classes, read %d\n",ncls); return 3; }
  u64 M=2ULL*H*H+1; spf=calloc(M+1,sizeof(unsigned)); if(!spf) return 3;
  for(u64 i=2;i<=M;i++) if(!spf[i]) for(u64 j=i;j<=M;j+=i) if(!spf[j]) spf[j]=(unsigned)i;
  u64 cap=(u64)H*H; R*A=malloc(cap*sizeof(R)), *B=malloc(cap*sizeof(R)); if(!A||!B) return 3;
  u64 na=0, nb=0;
  for(unsigned m=1;m<=H;m++) for(unsigned n=1;n<=H;n++){ if(gcdu(m,n)!=1) continue; fill(&B[nb++],m,n); if(n<m) fill(&A[na++],m,n); }
  qsort(B,nb,sizeof(R),cmp);
  u64 f1=0,f2=0,f3=0,f4=0,f5=0, triv=0;
  for(u64 i=0;i<na;i++){
    R key=A[i]; key.even=!A[i].even;                                   // the other number has the other parity
    u64 lo=0,hi=nb; while(lo<hi){ u64 mid=(lo+hi)/2; if(cmp(&B[mid],&key)<0) lo=mid+1; else hi=mid; }
    for(u64 j=lo;j<nb&&B[j].c==key.c&&B[j].even==key.even;j++){
      f1++; if(A[i].c==1) triv++;
      if(!A[i].f2) continue;
      f2++;
      if(!A[i].f3||!B[j].f3) continue;
      f3++;
      u64 m=A[i].m,n=A[i].n,mp=B[j].m,np=B[j].n, e,a,ii,c,g;
      i128 S=(i128)(m*m+n*n)*(mp*mp+np*np);                              // = 2 e^2
      if(!issq(S/2,&e)||(S&1)){ fprintf(stderr,"F1 inconsistent\n"); return 4; }
      i128 q=(i128)2*m*n*mp*np;
      if(!issq((i128)n*n*mp*mp+(i128)m*m*np*np,&a)) continue;
      if(!issq((i128)m*m*mp*mp+(i128)n*n*np*np,&ii)) continue;
      if(!issq((i128)e*e-q,&c)) continue;
      if(!issq((i128)e*e+q,&g)) continue;
      f4++;
      long long b=(long long)(m*mp+n*np), d=(long long)(m*mp)-(long long)(n*np), f=(long long)(n*mp+m*np), h=(long long)(n*mp)-(long long)(m*np);
      u64 v[9]={a,(u64)llabs(b),c,(u64)llabs(d),e,(u64)llabs(f),g,(u64)llabs(h),ii}; int ok=1;
      for(int x=0;x<9;x++){ if(v[x]==0) ok=0; for(int y=0;y<x;y++) if(v[x]==v[y]) ok=0; }
      if(ok){ f5++; printf("SOLUTION m=%llu n=%llu m'=%llu n'=%llu\n",m,n,mp,np); }
    }
  }
  printf("H=%u ratios: first=%llu second=%llu pairs=%llu | F1=%llu (of these with both numbers a square or twice a square: %llu) F2=%llu F3=%llu F4=%llu F5=%llu\n",
         H,na,nb,na*nb,f1,triv,f2,f3,f4,f5);
  return 0;
}
