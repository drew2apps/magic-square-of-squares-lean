\\ gen4(Curves, outfile): all cyclic 4-isogenies psi: A -> B with ker psi = <Q>, Q in A(K) of order 4, A in Curves,
\\ with T' in B(K) generating ker psihat; models Y^2 = X^3 + c2 X^2 + c4 X + c6 (X = 4x, Y = 8y + 4(a1 x + a3)).
bnf = bnfinit(y^2+1,1);
sh(C) = { [lift(C.b2), lift(8*C.b4), lift(16*C.b6)] };
shp(C, P) = { [lift(4*P[1]), lift(8*P[2] + 4*(C.a1*P[1] + C.a3))] };
alltors(C) = { my(T = elltors(C), G = T[3], pts = List(), n1 = T[2][1], n2 = if(#T[2]>1, T[2][2], 1));
  for(k=0,n1-1, for(l=0,n2-1, my(P = ellmul(C,G[1],k)); if(n2>1, P = elladd(C, P, ellmul(C,G[2],l))); if(P != [0], listput(pts, P))));
  Vec(pts) };
ord4(C) = select(P->ellorder(C,P)==4, alltors(C));
gen4(Curves, outfile) = {
  my(out = List());
  for(j=1,#Curves,
    my(A = ellinit(Curves[j], bnf), P4 = ord4(A), done = List());
    for(k=1,#P4, my(Q = P4[k]);
       if(#select(z->z==Q, Vec(done)), next);
       listput(done, Q); listput(done, ellneg(A,Q));
       my(B = ellinit(ellisogeny(A, Q, 1), bnf), tgt = -1);
       for(m=1,#Curves, if(ellinit(Curves[m],bnf).j == B.j, tgt = m));
       my(PB = ord4(B), Tp = 0);
       for(t=1,#PB, my(C2 = ellinit(ellisogeny(B, PB[t], 1), bnf)); if(C2.j == A.j, Tp = PB[t]; break));
       if(Tp == 0, error("no T'"));
       listput(out, [j, tgt, sh(A), shp(A, Q), sh(B), shp(B, Tp), apply(P->shp(A,P), alltors(A)), apply(P->shp(B,P), alltors(B))])));
  write(outfile, "D4 = ", Vec(out), ";"); #out };
