\\ Cyclic 4-isogenies with K-rational kernel generator in the 2-power K-isogeny class of the (627,364) curve, K = Q(i).
\\ For psi: A -> B = A/<Q>, Q in A(K) of order 4, find T' in B(K) of order 4 with B/<T'> = A (so ker psihat = <T'>).
\\ Output models Y^2 = X^3 + c2 X^2 + c4 X + c6 (X = 4x, Y = 8y + 4(a1 x + a3)) and the points in those models.
read("../isog_627_364.gp");
bnf = bnfinit(y^2+1,1);
sh(C) = { my(b2 = C.b2, b4 = C.b4, b6 = C.b6); [lift(b2), lift(8*b4), lift(16*b6)] };
shp(C, P) = { my(X = 4*P[1], Y = 8*P[2] + 4*(C.a1*P[1] + C.a3)); [lift(X), lift(Y)] };
ord4(C) = { my(T = elltors(C), G = T[3], pts = List());
  if(#G == 1, for(k=0,3, my(P = ellmul(C, G[1], k)); if(ellorder(C,P)==4, listput(pts,P))),
     for(k=0,3, for(l=0,1, my(P = elladd(C, ellmul(C,G[1],k), ellmul(C,G[2],l))); if(ellorder(C,P)==4, listput(pts,P)))));
  Vec(pts) };
out = List();
{for(j=1,#Curves,
  my(A = ellinit(Curves[j], bnf), P4 = ord4(A), done = List());
  for(k=1,#P4, my(Q = P4[k], sub = [Q, ellneg(A,Q)]);
     if(#select(z->z==Q || z==ellneg(A,Q), Vec(done)), next);
     listput(done, Q); listput(done, ellneg(A,Q));
     my(iso = ellisogeny(A, Q, 1), B = ellinit(iso, bnf), jB = B.j, tgt = -1, TB);
     \\ which curve is B
     for(m=1,#Curves, if(ellinit(Curves[m],bnf).j == jB, tgt = m));
     \\ find T' on B with B/<T'> = A
     my(PB = ord4(B), Tp = 0);
     for(t=1,#PB, my(C2 = ellinit(ellisogeny(B, PB[t], 1), bnf)); if(C2.j == A.j, Tp = PB[t]; break));
     if(Tp == 0, error("no T'"));
     print("A = E", j, "  Q = ", lift(Q), "  B = E", tgt, "  T' = ", lift(Tp));
     listput(out, [j, tgt, sh(A), shp(A, Q), sh(B), shp(B, Tp)]))) ;}
write("isog4_data.gp", "D4 = ", Vec(out), ";");
print(#out, " cyclic 4-isogenies with rational kernel generator");
