// Independent check of seven.c on a small range, with a different reduction.
// For every x <= X: H = all d > 0 with x^2 + d a square or x^2 - d a positive square. A magic square
// E, E±s, E±t, E±(s+t), E±(s-t) (E = x^2, s > t > 0, all entries positive) with at least 6 squares among the
// 8 other entries has at least two of its four lines with deviation in H. So: for every ordered pair (d1, d2)
// in H and every ordered pair of lines, solve the two linear equations for (s, t) and count the squares.
// No table of cases is used. Prints "x s t squares" for every solution found (possibly several times).
#include <stdio.h>
#include <stdlib.h>
int main(int argc,char**argv){
  if(argc!=2){ fprintf(stderr,"usage: seven_brute X   (X <= 3000)\n"); return 2; }
  long X=atol(argv[1]); if(X<1||X>3000){ fprintf(stderr,"need 1 <= X <= 3000\n"); return 2; }
  long M=3*X*X+10; char *sq=calloc(M,1); for(long r=1;r*r<M;r++) sq[r*r]=1;       /* 0 is not counted: entries are positive */
  long *H=malloc(4*X*sizeof(long)+64);
  static const int al[4]={1,0,1,1}, be[4]={0,1,1,-1};                             /* lines: s, t, s+t, s-t */
  for(long x=1;x<=X;x++){ long E=x*x; long nh=0;
    for(long d=1;d<2*E;d++) if(sq[E+d]||(d<E&&sq[E-d])) H[nh++]=d;
    for(long i=0;i<nh;i++) for(long j=0;j<nh;j++){ if(i==j) continue; long d1=H[i],d2=H[j];
      for(int p=0;p<4;p++) for(int q=0;q<4;q++){ if(p==q) continue;
        long det=al[p]*be[q]-al[q]*be[p]; long ns=d1*be[q]-d2*be[p], nt=al[p]*d2-al[q]*d1;
        if(ns%det||nt%det) continue;
        long s=ns/det,t=nt/det;
        if(!(s>t&&t>0)||s+t>=E) continue;
        int c=sq[E+s]+sq[E-s]+sq[E+t]+sq[E-t]+sq[E+s+t]+sq[E-s-t]+sq[E+s-t]+sq[E-s+t];
        if(c>=6) printf("%ld %ld %ld %d\n",x,s,t,c+1); } } }
  return 0; }
