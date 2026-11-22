// 3x3 magic squares with centre x^2 (L <= x <= N) and at least 7 square entries among nine distinct positive entries.
//
// Entries: E, E±s, E±t, E±(s+t), E±(s-t) with E = x^2 and s > t > 0: four lines through the centre, with
// deviations s, t, s+t, s-t. A line is "full" if E-d and E+d are both squares; then d = 2*leg1*leg2 for a
// Pythagorean triangle with hypotenuse x (E-d = (u-v)^2, E+d = (u+v)^2, u^2+v^2 = x^2), and conversely.
// Seven squares with a square centre means at least 6 squares on the four lines, hence at least two full lines.
// For two full lines with deviations a < b, the other two deviations are, according to which two of the four
// lines they are:
//   {t, s}     -> a+b, b-a            {s, s+t} or {t, s+t} -> b-a, |2a-b|      {s-t, s} -> b-a, 2b-a
//   {s-t, s+t} -> (a+b)/2, (b-a)/2    {t, s-t} -> a+b, 2a+b  or  a+b, a+2b  (a = t or a = s-t)
// All six cases are tested for every pair of triangles. Every hit is compared with the multiple of Bremner's
// square (x = 425k, deviations 97104k^2, 138600k^2, 41496k^2, 180096k^2); other hits are printed.
// Arithmetic: d <= x^2 <= 10^18, so every quantity is below 5*10^18 < 2^63. Usage: seven L N, N <= 10^9.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef unsigned long long u64; typedef long long i64;
static u64 gcd(u64 a, u64 b){ while(b){u64 t=a%b;a=b;b=t;} return a; }
static int issq(i64 v){ if(v<0) return 0; u64 r=(u64)sqrtl((long double)v); while(r*r>(u64)v) r--; while((r+1)*(r+1)<=(u64)v) r++; return r*r==(u64)v; }
int main(int argc,char**argv){
  if(argc!=3){ fprintf(stderr,"usage: seven L N   (1 <= L <= N <= 1000000000)\n"); return 2; }
  u64 L=strtoull(argv[1],0,10), N=strtoull(argv[2],0,10);
  if(L<1||N<L||N>1000000000ULL){ fprintf(stderr,"need 1 <= L <= N <= 10^9\n"); return 2; }
  u64 W=N-L+1;
  unsigned *cnt=calloc(W+2,sizeof(unsigned)); u64 *start=malloc((W+2)*sizeof(u64)); unsigned *fill=calloc(W+2,sizeof(unsigned));
  if(!cnt||!start||!fill){ fprintf(stderr,"out of memory\n"); return 3; }
  for(u64 m=2;m*m<=N;m++) for(u64 n=1;n<m&&m*m+n*n<=N;n++){ if(((m-n)&1)==0||gcd(m,n)!=1) continue; u64 h=m*m+n*n; for(u64 x=((L+h-1)/h)*h;x<=N;x+=h) cnt[x-L]++; }
  u64 tot=0; for(u64 x=0;x<W;x++){start[x]=tot;tot+=cnt[x];} start[W]=tot;
  u64 *D=malloc((tot+1)*sizeof(u64)); if(!D){ fprintf(stderr,"out of memory\n"); return 3; }
  for(u64 m=2;m*m<=N;m++) for(u64 n=1;n<m&&m*m+n*n<=N;n++){ if(((m-n)&1)==0||gcd(m,n)!=1) continue; u64 h=m*m+n*n,A=m*m-n*n,B=2*m*n; for(u64 k=(L+h-1)/h;k*h<=N;k++){ u64 x=k*h; D[start[x-L]+fill[x-L]++]=2*k*k*A*B; } }
  u64 pairs=0,bremner=0,other=0,by[10]={0};
  for(u64 x=L;x<=N;x++){ unsigned r=cnt[x-L]; if(r<2) continue; i64 E=(i64)(x*x); u64*d=D+start[x-L];
    for(unsigned i=0;i<r;i++) for(unsigned j=i+1;j<r;j++){ pairs++; i64 a=d[i],b=d[j]; if(a==b) continue; if(a>b){i64 t=a;a=b;b=t;}
      int ev=((a+b)&1)==0;
      i64 c[6][2]={{a+b,b-a},{b-a,llabs(2*a-b)},{b-a,2*b-a},{ev?(a+b)/2:0,ev?(b-a)/2:0},{a+b,2*a+b},{a+b,a+2*b}};
      for(int k=0;k<6;k++){ i64 r1=c[k][0],r2=c[k][1]; if(r1<=0||r2<=0||r1==r2||r1==a||r1==b||r2==a||r2==b) continue;
        if(E-r1<=0||E-r2<=0) continue;
        int s=issq(E+r1)+issq(E-r1)+issq(E+r2)+issq(E-r2);
        if(s>=2){ by[5+s]++; i64 q=(i64)(x/425); i64 lo=r1<r2?r1:r2, hi=r1<r2?r2:r1;
          if(x%425==0 && s==2 && a==97104*q*q && b==138600*q*q && lo==41496*q*q && hi==180096*q*q) bremner++;
          else { other++; printf("OTHER x=%llu full=(%lld,%lld) rest=(%lld,%lld) squares=%d\n",x,a,b,r1,r2,5+s); fflush(stdout);} } } } }
  printf("%llu<=x<=%llu triangles=%llu pairs=%llu hits: 7 squares=%llu 8 squares=%llu 9 squares=%llu; multiples of Bremner's square=%llu other=%llu\n",L,N,tot,pairs,by[7],by[8],by[9],bremner,other);
  return 0; }
