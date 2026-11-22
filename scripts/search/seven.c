// 3x3 magic squares with centre x^2 and at least 7 square entries (distinct entries), L <= x <= N.
// Usage: seven L N   (memory is proportional to N - L; x up to 10^9 fits in 64-bit arithmetic)
// Entries: E, E±s, E±t, E±(s+t), E±(s-t), E = x^2. A "full line" has deviation d with E-d and E+d both squares:
// d = 2*leg1*leg2 for a Pythagorean triangle with hypotenuse x. Seven squares with a square centre need two
// full lines; the other two deviations are then determined up to four cases, and need two more squares.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
typedef unsigned long long u64; typedef long long i64;
static u64 gcd(u64 a, u64 b){ while(b){u64 t=a%b;a=b;b=t;} return a; }
static int issq(i64 v){ if(v<0) return 0; u64 r=(u64)sqrtl((long double)v); while(r*r>(u64)v) r--; while((r+1)*(r+1)<=(u64)v) r++; return r*r==(u64)v; }
int main(int argc,char**argv){
  u64 L=strtoull(argv[1],0,10), N=strtoull(argv[2],0,10); u64 W=N-L+1;
  unsigned *cnt=calloc(W+2,sizeof(unsigned));
  for(u64 m=2;m*m<=N;m++) for(u64 n=1;n<m&&m*m+n*n<=N;n++){ if(((m-n)&1)==0||gcd(m,n)!=1) continue; u64 h=m*m+n*n; for(u64 x=((L+h-1)/h)*h;x<=N;x+=h) cnt[x-L]++; }
  u64 *start=malloc((W+2)*sizeof(u64)); u64 tot=0; for(u64 x=0;x<W;x++){start[x]=tot;tot+=cnt[x];} start[W]=tot;
  u64 *D=malloc((tot+1)*sizeof(u64)); unsigned *fill=calloc(W+2,sizeof(unsigned));
  for(u64 m=2;m*m<=N;m++) for(u64 n=1;n<m&&m*m+n*n<=N;n++){ if(((m-n)&1)==0||gcd(m,n)!=1) continue; u64 h=m*m+n*n,A=m*m-n*n,B=2*m*n; for(u64 k=(L+h-1)/h;k*h<=N;k++){ u64 x=k*h; D[start[x-L]+fill[x-L]++]=2*k*k*A*B; } }
  fprintf(stderr,"representations: %llu\n",tot);
  u64 pairs=0,hits=0;
  for(u64 x=(L?L:1);x<=N;x++){ unsigned r=cnt[x-L]; if(r<2) continue; i64 E=(i64)(x*x); u64*d=D+start[x-L];
    for(unsigned i=0;i<r;i++) for(unsigned j=i+1;j<r;j++){ pairs++; i64 a=d[i],b=d[j]; if(a==b) continue; if(a>b){i64 t=a;a=b;b=t;}
      i64 c[4][2]={{a+b,b-a},{b-a,llabs(2*a-b)},{b-a,2*b-a},{((a+b)&1)?0:(a+b)/2,((a+b)&1)?0:(b-a)/2}};
      for(int k=0;k<4;k++){ i64 r1=c[k][0],r2=c[k][1]; if(r1<=0||r2<=0||r1==r2||r1==a||r1==b||r2==a||r2==b) continue;
        int s=issq(E+r1)+issq(E-r1)+issq(E+r2)+issq(E-r2);
        if(s>=2 && E-r1>0 && E-r2>0){ hits++; printf("x=%llu full=(%lld,%lld) rest=(%lld,%lld) squares=%d\n",x,a,b,r1,r2,5+s); fflush(stdout);} } } }
  fprintf(stderr,"%llu<=x<=%llu pairs=%llu hits=%llu\n",L,N,pairs,hits); return 0; }
