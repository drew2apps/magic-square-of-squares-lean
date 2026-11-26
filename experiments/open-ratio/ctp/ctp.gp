lsqrt_fail = 0;
UNITRED = 0;  \\ 1 turns on the unit reduction in conicpt (slow: it computes a class group of the quartic field)
tr(x) = if(type(x) == "t_PADIC", truncate(x), x);   \\ p-adic to rational; leaves exact rationals alone
\\ Cassels-Tate pairing on Sel^2 for E: y^2=(x-e1)(x-e2)(x-e3) over a number field K (variable y), full 2-torsion.
\\ Fisher, "On binary quartics and the Cassels-Tate pairing" (arXiv 2208.14977), Theorem 3.1, with L = K^3 split.
\\ EXPERIMENT, written by an AI model (Claude). Not reviewed.

\\ ---------- local square classes via Hilbert symbols ----------
\\ a "place" is [pr] for a prime ideal, or [0] for the real place of Q
isre(pl) = type(pl[1]) == "t_INT";
hil(nf, a, b, pl) = if(isre(pl), if(simplify(a)<0 && simplify(b)<0, -1, 1), nfhilbert(nf, a, b, pl[1]));
expdim(nf, pl) = { if(isre(pl), return(1)); my(pr=pl[1]); if(pr.p==2, 2 + pr.e*pr.f, 2); }
\\ spanning set of K_v*/K_v*^2
lsqbasis(nf, pl) = { my(d = expdim(nf, pl), C = List(), B = List(), M, p);
  if(isre(pl), return([-1]));
  p = pl[1].p;
  listput(C, p); listput(C, nfbasistoalg(nf, pl[1].gen[2]));
  for(a=-7, 7, for(b=-7*(poldegree(nf.pol)>1), 7*(poldegree(nf.pol)>1), if(a||b, listput(C, a + b*y))));
  for(k=1, 30, listput(C, random(2*p^3)+1 - p^3 + (poldegree(nf.pol)>1)*(random(2*p^3)-p^3)*y));
  C = select(c -> c != 0, Vec(C));
  \\ greedy: add candidates whose Hilbert vector against C increases rank
  M = matrix(0, #C);
  for(i=1, #C, my(r = vector(#C, j, (1-hil(nf, C[i], C[j], pl))/2)); my(M2 = concat(M, r));
     if(matrank(Mod(M2,2)) > matrank(Mod(M,2)), M = M2; listput(B, C[i])); if(#B==d, break));
  if(#B != d, error("lsqbasis: only ", #B, " of ", d));
  Vec(B); }
lvec(nf, u, B, pl) = vector(#B, j, (1 - hil(nf, u, B[j], pl))/2);
lsq(nf, u, B, pl) = vecmax(lvec(nf, u, B, pl)) == 0;

\\ ---------- curve data ----------
fx(e, t) = (t-e[1])*(t-e[2])*(t-e[3]);
\\ places: primes dividing 2*prod(e_i-e_j), plus real place if K=Q
badplaces(nf, e) = { my(P = Set([2]), L = List(), N);
  for(i=1,3, for(j=i+1,3, N = abs(norm(Mod(e[i]-e[j], nf.pol))); P = setunion(P, Set(factor(numerator(N))[,1]~)); P = setunion(P, Set(factor(denominator(N))[,1]~))));
  for(k=1, #P, my(D = idealprimedec(nf, P[k])); for(t=1, #D, listput(L, [D[t]])));
  if(poldegree(nf.pol)==1, listput(L, [0]));
  Vec(L); }

\\ Kummer image (x-e1, x-e2) of a sampled local point
kum(e, t) = { if(t==e[1], return([(e[1]-e[2])*(e[1]-e[3]), e[1]-e[2]]));
  if(t==e[2], return([e[2]-e[1], (e[2]-e[1])*(e[2]-e[3])])); [t-e[1], t-e[2]]; }
randK(nf, R) = random(2*R+1)-R + (poldegree(nf.pol)>1)*(random(2*R+1)-R)*y;
\\ local image of E(K_v)/2E(K_v) in F_2^(2 dim)
locimg(nf, e, pl, B) = { my(d = #B, tgt = expdim(nf, pl), rows = List(), M, pi, v, t, w);
  listput(rows, concat(lvec(nf, kum(e,e[1])[1], B, pl), lvec(nf, kum(e,e[1])[2], B, pl)));
  listput(rows, concat(lvec(nf, kum(e,e[2])[1], B, pl), lvec(nf, kum(e,e[2])[2], B, pl)));
  M = Mat(Vec(rows)~)~; M = matrix(#rows, 2*d, i, j, rows[i][j]);
  pi = if(isre(pl), 1, nfbasistoalg(nf, pl[1].gen[2])); if(!isre(pl) && pl[1].e==1, pi = pl[1].p);
  for(k=1, 20000, if(matrank(Mod(M,2)) == tgt, break);
     my(md = k%4);
     if(md==0, t = e[random(3)+1] + pi^random(12)*randK(nf, 60),
        md==1, t = randK(nf, 10^4),
        md==2, t = randK(nf, 300)/pi^(2*(random(4)+1)),
        t = pi^random(8)*randK(nf, 300));
     if(red(fx(e,t))==0, next);
     if(!lsq(nf, fx(e,t), B, pl), next);
     my(kk = kum(e,t)); w = concat(lvec(nf, kk[1], B, pl), lvec(nf, kk[2], B, pl));
     my(M2 = matconcat([M; w])); if(matrank(Mod(M2,2)) > matrank(Mod(M,2)), M = M2));
  if(matrank(Mod(M,2)) != tgt, error("locimg: rank ", matrank(Mod(M,2)), " of ", tgt, " at ", pl));
  M; }

\\ 2-Selmer group: returns list of pairs [d1, d2]
selmer2(bnf, e) = { my(nf = bnf.nf, PL = badplaces(nf, e), S = List(), gens, g, rows = List(), K, out = List());
  for(k=1, #PL, if(!isre(PL[k]), listput(S, PL[k][1])));
  my(su = bnfsunit(bnf, Vec(S)));
  gens = concat([nfbasistoalg(nf, bnf.tu[2])], concat(apply(u->nfbasistoalg(nf,u), bnf.fu), apply(u->nfbasistoalg(nf,u), su[1])));
  gens = apply(u -> lift(Mod(u, nf.pol)), gens);
  g = #gens;
  for(k=1, #PL, my(pl = PL[k], B = lsqbasis(nf, pl), d, W, Mv, Ann);
     d = #B; W = locimg(nf, e, pl, B);
     Ann = lift(matker(Mod(W,2)~)~);   \\ rows a with a.w = 0 for w in image... (left annihilator of row space)
     Ann = lift(matker(Mod(W,2)))~;   \\ kernel of W: vectors v with W*v=0 i.e. orthogonal to all rows
     \\ map: s in F2^(2g) -> (sum s_i vec(gens_i), sum s_{g+i} vec(gens_i))
     my(GV = matrix(d, g, a, b, lvec(nf, gens[b], B, pl)[a]));
     my(Mmap = matconcat([GV, matrix(d,g); matrix(d,g), GV]));
     my(C = Ann * Mmap); for(r=1, #C~, listput(rows, C[r,])));
  my(A = matrix(#rows, 2*g, i, j, rows[i][j]));
  K = lift(matker(Mod(A,2)));
  for(c=1, #K, my(v = K[,c]); listput(out, [prod(i=1,g, gens[i]^v[i]), prod(i=1,g, gens[i]^v[g+i])]));
  [Vec(out), gens, PL]; }

\\ ---------- 2-coverings as binary quartics ----------
\\ solve A X^2 + B Y^2 + C W^2 = 0 over K, return [X,Y,W]
conicpt(bnf, A, B, C) = { my(nf = bnf.nf, r, s, t, T, a);
  if(poldegree(nf.pol)==1, my(v = qfsolve(matdiagonal([A,B,C]))); if(type(v)!="t_COL", error("conic over Q")); return(v~));
  r = red(-B/A); s = red(-C/A);
  t = nfroots(nf, x^2 - r);
  if(#t, t = lift(t[1]); return([(1+s)/2, (s-1)/(2*t), 1]));
  if(type(r) != "t_POL" && type(s) != "t_POL", my(v = qfsolve(matdiagonal([1, -r, -s]))); if(type(v) == "t_COL", return([v[1]/v[3], v[2]/v[3], 1]~)));
  my(D = denominator(content(r)), rr = r*D^2);   \\ integral model: X^2 - rr Y'^2 = s with Y = D Y'
  T = rnfisnorminit(nf.pol, x^2 - rr);
  a = iferr(rnfisnorm(T, s, 100), E, rnfisnorm(T, s, 0));
  if(a[2] != 1, error("conicpt: no solution found"));
  a = lift(a[1]); a = if(!UNITRED, [red(polcoef(a,0,x)), red(polcoef(a,1,x))], iferr(unitreduce(nf, rr, T, red(polcoef(a,0,x)), red(polcoef(a,1,x))), E, [red(polcoef(a,0,x)), red(polcoef(a,1,x))]));
  [a[1], D*a[2], 1]; }
\\ alpha = X + Y w (w^2 = rr) has N_{L/K}(alpha) = s; multiply by the power of a norm-1 unit of L that balances |alpha| at the
\\ two complex places of L, which removes the huge unit factor rnfisnorm tends to leave
unitreduce(nf, rr, T, X, Y) = { my(rnf, Lb, fu, NE, m, eta, al, rw, l1, l2, lam, n, M);
  if(poldegree(nf.pol) == 1, return([X, Y]));
  rnf = rnfinit(nf, x^2 - rr); Lb = bnfinit(rnf.polabs, 1);
  if(#Lb.fu == 0, return([X, Y]));
  fu = lift(rnfeltabstorel(rnf, Lb.fu[1]));
  M = Mod(fu, x^2 - rr); NE = red(lift(M * subst(fu, x, -x)));
  m = 1; while(red(NE^m) != 1, m++; if(m > 4, return([X, Y])));
  eta = M^m;
  rw = sqrt(subst(lift(Mod(rr, y^2+1)), y, I));
  my(av(c, sg) = abs(subst(subst(liftall(c), x, sg*rw), y, I)));
  al = Mod(X + Y*x, x^2 - rr);
  if(av(al, 1) < 1e-30 || av(al, -1) < 1e-30 || av(eta, 1) < 1e-30, return([X, Y]));
  l1 = log(av(al, 1)); l2 = log(av(al, -1)); lam = log(av(eta, 1));
  if(abs(lam) < 1e-6, return([X, Y]));
  n = round((l2 - l1)/(2*lam));
  al = liftall(al * eta^n);
  [red(polcoef(al,0,x)), red(polcoef(al,1,x))]; }
\\ quartic g (Pol in x, z=1 dehomogenized, degree <= 4) of the 2-covering x - e_i = d_i z_i^2
quart(bnf, e, d1, d2) = { my(d3 = d1*d2, c2 = e[2]-e[1], c3 = e[3]-e[1], P0, U, V, R, F, B, P, g);
  P0 = conicpt(bnf, d1, -d2, -c2);
  \\ pick u, v completing P0 to a basis
  my(E3 = [[1,0,0],[0,1,0],[0,0,1]], ok = 0);
  for(i=1,3, for(j=i+1,3, if(!ok && matdet(Mat([P0~, E3[i]~, E3[j]~]))!=0, U = E3[i]; V = E3[j]; ok = 1)));
  R = x*U + V;
  P = (d1*R[1]^2 - d2*R[2]^2 - c2*R[3]^2)*P0 - 2*(d1*P0[1]*R[1] - d2*P0[2]*R[2] - c2*P0[3]*R[3])*R;
  g = d3*(d1*P[1]^2 - c3*P[3]^2);
  g = redp(g);
  g; }
qcoef(g) = vector(5, k, polcoef(g, 5-k, x));   \\ [a,b,c,d,e]
qIJ(g) = { my(q = qcoef(g), a=q[1], b=q[2], c=q[3], d=q[4], ee=q[5]);
  [12*a*ee - 3*b*d + c^2, 72*a*c*ee - 27*a*d^2 - 27*b^2*ee + 9*b*c*d - 2*c^3]; }
qhess(g) = { my(q = qcoef(g), a=q[1], b=q[2], c=q[3], d=q[4], ee=q[5]);
  (3*b^2-8*a*c)*x^4 + 4*(b*c-6*a*d)*x^3 + 2*(2*c^2-24*a*ee-3*b*d)*x^2 + 4*(c*d-6*b*ee)*x + (3*d^2-8*c*ee); }
red(z) = lift(Mod(z, nf_pol_glob));
redp(g) = if(type(g)!="t_POL" || variable(g)!=x, red(g), Pol(apply(c->red(c), Vec(g)), x));

\\ ---------- the pairing (Fisher Thm 3.1, L split) ----------
nfsqrt(nf, v) = { my(r = nfroots(nf, x^2 - v)); if(#r==0, return(0)); lift(r[1]); }
evq(g, X, Z) = { my(q = qcoef(g)); q[1]*X^4 + q[2]*X^3*Z + q[3]*X^2*Z^2 + q[4]*X*Z^3 + q[5]*Z^4; }
evq2(h, X, Z) = polcoef(h,2,x)*X^2 + polcoef(h,1,x)*X*Z + polcoef(h,0,x)*Z^2;
\\ scale quartic g to have invariants [I0,J0]; returns 0 if impossible
scaleto(nf, g, IJ0) = { my(IJ = qIJ(g), u2 = red(IJ[2]*IJ0[1]/(IJ0[2]*IJ[1])));
  if(nfsqrt(nf, u2) == 0, error("scaleto: u^2 not a square"));
  redp(g / u2); }
\\ H for (g, phi): G = (4 phi g + hess)/3, H quadratic with H^2 = G(1,0) G, H(1,0) = G(1,0)
Hform(g, phi, zconv) = { my(G = redp((4*phi*g + qhess(g))/(if(zconv==3, 1, 3))), z0 = polcoef(G,4,x), q = redp(z0*G), h2, h1, h0, H);
  h2 = z0; h1 = red(polcoef(q,3,x)/(2*h2)); h0 = red((polcoef(q,2,x) - h1^2)/(2*h2));
  H = h2*x^2 + h1*x + h0;
  if(redp(H^2 - q) != 0, error("Hform: not a square"));
  [z0, H]; }
\\ local point on y^2 = g1 at place pl with gamma(x) != 0
locpt(nf, g1, gam, pl, B, nt) = { my(pi, t, lst = List());
  pi = if(isre(pl), 1, nfbasistoalg(nf, pl[1].gen[2])); if(!isre(pl) && pl[1].e==1, pi = pl[1].p);
  pi = lift(pi);
  if(isre(pl), my(rr = select(r -> abs(imag(r)) < 1e-20, polroots(g1)), cand = List());
     rr = vecsort(real(rr));
     for(i=1, #rr-1, listput(cand, bestappr((rr[i]+rr[i+1])/2, 10^6)));
     if(#rr, listput(cand, floor(rr[1])-1); listput(cand, ceil(rr[#rr])+1));
     foreach(cand, X, my(gv = simplify(red(evq(g1, X, 1))), cv = simplify(red(evq2(gam, X, 1)))); if(gv > 0 && cv != 0, listput(lst, [X,1]))));
  for(k=1, 40000, if(#lst >= nt, break);
    my(md = k%4, X, Z = 1);
    if(md==0, X = randK(nf, 30), md==1, X = randK(nf, 10^4), md==2, X = 1; Z = pi^random(6)*randK(nf, 50), X = pi^random(8)*randK(nf,300));
    if(Z == 0, next);
    my(gv = simplify(red(evq(g1, X, Z))), cv = simplify(red(evq2(gam, X, Z))));
    if(gv == 0 || cv == 0, next);
    if(isre(pl), if(gv > 0, listput(lst, [X,Z])); next);
    if(lsq(nf, gv, B, pl), listput(lst, [X,Z])));
  Vec(lst); }
\\ <[g1],[g2]> with g3 for the third class; returns [values over sign choices of m, per-place data]
\\ add the places dividing the norms of the given elements (numerators and denominators)
augment(nf, PL, els) = { my(P = Set(), L = List(PL), have = Set());
  for(k=1, #PL, if(!isre(PL[k]), have = setunion(have, Set([PL[k][1].p]))));
  for(k=1, #els, my(c = els[k], N); if(c == 0, next); N = abs(norm(Mod(c, nf.pol)));
     P = setunion(P, Set(factor(numerator(N), 10^9)[,1]~)); P = setunion(P, Set(factor(denominator(N), 10^9)[,1]~)));
  for(k=1, #P, if(P[k] > 1 && !setsearch(have, P[k]), if(!isprime(P[k]), print("  warning: unfactored ", P[k])); my(D = idealprimedec(nf, P[k])); for(t=1,#D, listput(L,[D[t]]))));
  Vec(L); }
ctpair(bnf, g1, g2, g3, PL, opt) = { my(nf = bnf.nf, IJ = qIJ(g1), ph, res = List(), zc = if(opt, opt, 1));
  ph = apply(lift, nfroots(nf, x^3 - 3*IJ[1]*x + IJ[2]));
  if(#ph != 3, error("ctpair: L not split"));
  my(D = vector(3, j, prod(k=1,3, if(k==j, 1, ph[j]-ph[k]))));
  my(Z = vector(3, k, vector(3, j, Hform([g1,g2,g3][k], ph[j], zc))));
  my(m0 = vector(3, j, nfsqrt(nf, red(Z[1][j][1]*Z[2][j][1]*Z[3][j][1]))));
  if(vecmin(apply(v->v!=0, m0)) == 0, return(["no-m", m0]));
  forvec(s = vector(3, i, [0,1]),
    my(m = vector(3, j, (-1)^s[j]*m0[j]));
    my(gam = redp(sum(j=1,3, red(Z[2][j][1]*Z[3][j][1]/m[j]) * Z[1][j][2] / D[j])));
    my(a2 = polcoef(g2, 4, x), tot = 0, per = List());
    if(a2 == 0, error("g2(1,0)=0"));
    my(PLx = augment(nf, PL, concat([a2, red(4*IJ[1]^3 - IJ[2]^2), polcoef(gam,2,x), polcoef(gam,1,x), polcoef(gam,0,x)], qcoef(g1))));
    for(k=1, #PLx, my(pl = PLx[k], B = lsqbasis(nf, pl), pts = locpt(nf, g1, gam, pl, B, 3), vals);
      if(#pts == 0, error("no local point at ", pl));
      vals = vector(#pts, i, my(c = red(evq2(gam, pts[i][1], pts[i][2]))); (1 - hil(nf, a2, c, pl))/2);
      if(vecmin(vals) != vecmax(vals), listput(per, ["inconsistent", pl[1], vals]));
      tot += vals[1]; listput(per, vals[1]));
    listput(res, [s, tot % 2, Vec(per)]));
  Vec(res); }
\\ g(alpha x + beta z, gamma x + delta z), det must be 1 to keep I, J
qtrans(g, al, be, ga, de) = { my(q = qcoef(g)); redp(sum(k=0, 4, q[k+1]*(al*x+be)^(4-k)*(ga*x+de)^k)); }
\\ random det-1 change of variables so that z(g) != 0 on every component
qgood(nf, g) = { my(IJ = qIJ(g), ph = apply(lift, nfroots(nf, x^3 - 3*IJ[1]*x + IJ[2])), h = g);
  for(t=0, 50, my(ok = polcoef(h,4,x) != 0);
     if(ok, for(j=1, #ph, if(red(4*ph[j]*polcoef(h,4,x) + polcoef(qhess(h),4,x)) == 0, ok = 0)));
     if(ok, return(h));
     my(k = random(21)-10); h = qtrans(g, 1, k, random(7)-3, 0); h = qtrans(g, k, -1, 1, 0); h = qtrans(qtrans(g, 1, random(11)-5, 0, 1), 0, -1, 1, random(9)-4));
  error("qgood failed"); }

\\ ================= Cassels' formula (Shukla-Stoll arXiv 2302.01640, eq. (3.10) and Thm 5.6) =================
\\ <a,a'> = prod_v prod_{i=1..3} ( L_i(q_v), b'_i )_v ,  e_i in K.
\\ H_i = (b_j G_j^2 - b_k G_k^2)/(e_j - e_k) + T^2, (i,j,k) cyclic; q_i = (G*_j : G*_k : 1) a K-point on H_i = 0;
\\ L_i = 2 (b_j G*_j G_j - b_k G*_k G_k)/(e_j - e_k) + 2 T;  q_v = (w_1 : w_2 : w_3 : 1) with x(Q) - e_i = b_i w_i^2, Q in E(K_v).
\\ local field contexts: [1, p, N, iota] for K_v = Q_p (K = Q, or p split in Q(i)); [2, p, N] for K_v = Q_p(i); [0] real.
lctx(nf, pl, N) = { if(isre(pl), return([0])); my(pr = pl[1], p = pr.p);
  if(poldegree(nf.pol) == 1, return([1, p, N, 0]));
  if(pr.f == 1 && pr.e == 1, my(r = polrootsmod(x^2+1, p), i0 = 0);
     for(k=1, #r, if(nfeltval(nf, y - lift(r[k]), pr) > 0, i0 = lift(r[k])));
     my(io = sqrt(-1 + O(p^N))); if(tr(io) % p != i0 % p, io = -io); return([1, p, N, io]));
  [2, p, N]; }
emb(c, ctx) = { my(t0 = ctx[1]); c = lift(c); if(type(c)=="t_POL" || type(c)=="t_RFRAC", c = red(c));
  if(t0 == 0, return(if(type(c)=="t_POL", polcoef(c,0,y), c)*1.));
  if(t0 == 1, return(if(type(c)=="t_POL", polcoef(c,0,y) + polcoef(c,1,y)*ctx[4], c) + O(ctx[2]^ctx[3])));
  Mod(if(type(c)=="t_POL", polcoef(c,0,y) + polcoef(c,1,y)*t, c) + O(ctx[2]^ctx[3]), t^2+1); }
back(z, ctx) = { if(ctx[1] == 1, return(tr(z))); my(w = lift(z)); tr(polcoef(w,0,t)) + tr(polcoef(w,1,t))*y; }
\\ square root in Q_2(i): strip (1+i)^v, find a root mod 8 by search, then Newton
lsqrt2(z, ctx) = { my(N = ctx[3], w = lift(z), a = tr(polcoef(w,0,t)), b = tr(polcoef(w,1,t)), vn, z1, a1, b1, w0 = 0, P = N + 24, u, h);
  if(a == 0 && b == 0, error("lsqrt2: zero"));
  vn = valuation(a^2 + b^2, 2); if(vn % 2, error("lsqrt2: odd valuation"));
  z1 = lift(Mod(a + b*t, t^2+1) / Mod(2*t, t^2+1)^(vn/2)); a1 = polcoef(z1,0,t); b1 = polcoef(z1,1,t);
  for(A=0, 7, for(B=0, 7, if(w0 == 0, my(r = lift(Mod((A + B*t)^2, t^2+1)) - z1);
     if(valuation(polcoef(r,0,t), 2) >= 3 && valuation(polcoef(r,1,t), 2) >= 3, w0 = A + B*t))));
  if(w0 == 0, error("lsqrt2: not a square"));
  u = Mod(polcoef(w0,0,t) + O(2^P) + (polcoef(w0,1,t) + O(2^P))*t, t^2+1);
  h = Mod(a1 + O(2^P) + (b1 + O(2^P))*t, t^2+1);
  for(k=1, ceil(log(P)/log(2)) + 3, u = (u + h/u)/2);
  u * Mod(1 + t, t^2+1)^(vn/2); }
lsqrt(z, ctx) = { my(t0 = ctx[1]);
  if(t0 == 0, return(sqrt(z)));
  if(t0 == 1, return(sqrt(z)));
  my(w = lift(z), a = polcoef(w,0,t), b = polcoef(w,1,t), Nr, u, v);
  if(ctx[2] == 2, return(lsqrt2(z, ctx)));
  if(b == 0 && issquare(a), return(Mod(sqrt(a), t^2+1)));
  if(b == 0, return(Mod(sqrt(-a)*t, t^2+1)));
  \\ u in Q_p with u^4 - a u^2 - b^2/4 = 0, v = b/(2u); found with polrootspadic on the rational truncations
  my(p = ctx[2], N = ctx[3], ar = tr(a), br = tr(b), R = polrootspadic(x^4 - ar*x^2 - br^2/4, p, N));
  if(#R == 0, error("lsqrt: not a square in Q_p(i)"));
  for(k=1, #R, u = R[k]; if(u == 0, next); v = (br + O(p^N))/(2*u);
    my(c = Mod(u + v*t, t^2+1), d = lift(c^2 - z));
    if(valuation(polcoef(d,0,t), p) >= N/2 && valuation(polcoef(d,1,t), p) >= N/2, return(c)));
  lsqrt_fail++; if(lsqrt_fail <= 5, print("  lsqrt fail p=", p, " val ", [valuation(a,p), valuation(b,p)], " #R ", #R)); error("lsqrt: no root validated"); }
\\ K-point (X, Y) with A X^2 + B Y^2 + 1 = 0
conicaff(bnf, A, B) = { my(P = conicpt(bnf, A, B, 1), F, Bl, R, Q);
  if(P[3] != 0, return([red(P[1]/P[3]), red(P[2]/P[3])]));
  for(k=1, 50, R = [randK(bnf.nf, 9), randK(bnf.nf, 9), 1];
    Q = (A*R[1]^2 + B*R[2]^2 + R[3]^2)*P - 2*(A*P[1]*R[1] + B*P[2]*R[2] + P[3]*R[3])*R;
    if(red(Q[3]) != 0, return([red(Q[1]/Q[3]), red(Q[2]/Q[3])])));
  error("conicaff"); }
\\ global data for class b = [b1,b2,b3]
casdata(bnf, e, b) = { my(G = vector(3), cyc = [[1,2,3],[2,3,1],[3,1,2]]);
  for(t=1, 3, my(i = cyc[t][1], j = cyc[t][2], k = cyc[t][3], s = e[j] - e[k]);
     G[i] = conicaff(bnf, red(b[j]/s), red(-b[k]/s)));   \\ G[i] = [G*_j, G*_k]
  G; }
\\ local term at place pl for x0 with Kummer image b at pl
caslocal(nf, e, b, bp, G, pl, x0, N) = { my(ctx = lctx(nf, pl, N), w = vector(3), cyc = [[1,2,3],[2,3,1],[3,1,2]], tot = 0);
  for(i=1, 3, w[i] = iferr(lsqrt(emb(red((x0 - e[i])/b[i]), ctx), ctx), E, 0); if(w[i] === 0, return(-1)));
  for(t=1, 3, my(i = cyc[t][1], j = cyc[t][2], k = cyc[t][3], s = e[j] - e[k], L);
     L = 2*(emb(b[j]*G[i][1], ctx)*w[j] - emb(b[k]*G[i][2], ctx)*w[k])/emb(s, ctx) + 2;
     if(ctx[1] == 0, tot += (L < 0 && simplify(bp[i]) < 0),
        my(Lg = back(L, ctx)); if(Lg == 0, return(-1)); tot += (1 - nfhilbert(nf, Lg, bp[i], pl[1]))/2));
  tot % 2; }
\\ the pairing; PL: places to use. Returns [value, per-place values]
cassels(bnf, e, b, bp, PL0, N) = { my(nf = bnf.nf, G = casdata(bnf, e, b), per = List(), tot = 0, PL);
  PL = augment(nf, PL0, concat(concat(b, bp), concat(G[1], concat(G[2], G[3]))));
  for(k=1, #PL, my(pl = PL[k], B = lsqbasis(nf, pl), vals = List(), x0, cnt = 0);
    for(tr=1, 60000, if(#vals >= 3, break);
      my(md = tr % 4, pi = if(isre(pl), 1, if(pl[1].e==1, pl[1].p, lift(nfbasistoalg(nf, pl[1].gen[2])))));
      x0 = if(md==0, e[random(3)+1] + pi^random(10)*randK(nf, 60), md==1, randK(nf, 10^4), md==2, randK(nf, 300)/pi^(2*(random(4)+1)), pi^random(8)*randK(nf,300));
      if(red(fx(e, x0)) == 0, next);
      if(vecmax(concat(vector(3, i, lvec(nf, red((x0 - e[i])/b[i]), B, pl)))) != 0, next);
      my(v = caslocal(nf, e, b, bp, G, pl, x0, N)); if(v < 0, next); listput(vals, v));
    if(#vals == 0, error("cassels: no local point at ", pl));
    if(vecmin(Vec(vals)) != vecmax(Vec(vals)), print("  INCONSISTENT at ", if(isre(pl), "inf", pl[1].p), " ", Vec(vals)));
    tot += vals[1]; listput(per, vals[1]));
  [tot % 2, Vec(per), G, apply(pl->if(isre(pl),0,[pl[1].p, pl[1].f]), PL)]; }

\\ ---------- pairing without factoring the conic points ----------
\\ At a good place P outside PL (P odd, beta, beta', e_i - e_j units), L_i = -Q_i(Gamma - w)/(e_j - e_k), and the local term is
\\ sum_i v_P(Z_i) [beta'_i non-square at P] = sum_i (Z_i, beta'_i)_P, Z_i a generator of the common denominator ideal of Gamma_i.
\\ Hilbert reciprocity then gives the sum over all P outside PL as sum_{v in PL} (Z_i, beta'_i)_v, so no factoring is needed.
dengen(bnf, g) = { my(nf = bnf.nf, D = 1);
  foreach(g, c, my(d = denominator(content(lift(c))), a = red(c*d), Jd = idealdiv(nf, d, idealadd(nf, d, a))); D = idealintersect(nf, D, Jd));
  if(poldegree(nf.pol) == 1, return(D[1,1]));
  if(nf.pol != y^2 + 1, error("dengen: only Q and Q(i)"));
  D = idealhnf(nf, D); ggcd(D[1,1], D[1,2] + D[2,2]*y); }
\\ gcd in Z[i] (elements as polynomials in y), Euclid with rounded quotients
ggcd(a, b) = { a = red(a); b = red(b);
  while(b != 0, my(q = red(a * red(conj_i(b))) / norm(Mod(b, y^2+1)), qr = round(polcoef(q,0,y)) + round(polcoef(q,1,y))*y, r = red(a - qr*b)); a = b; b = r);
  a; }
conj_i(z) = subst(lift(Mod(z, y^2+1)), y, -y);
cassels2(bnf, e, b, bp, PL0, N) = { my(nf = bnf.nf, G = casdata(bnf, e, b), per = List(), tot = 0, PL, Z, corr = 0);
  PL = augment(nf, PL0, concat(b, bp));
  for(k=1, #PL, my(pl = PL[k], B = lsqbasis(nf, pl), vals = List(), x0);
    for(tr=1, 60000, if(#vals >= 3, break);
      my(md = tr % 4, pi = if(isre(pl), 1, if(pl[1].e==1, pl[1].p, lift(nfbasistoalg(nf, pl[1].gen[2])))));
      x0 = if(md==0, e[random(3)+1] + pi^random(10)*randK(nf, 60), md==1, randK(nf, 10^4), md==2, randK(nf, 300)/pi^(2*(random(4)+1)), pi^random(8)*randK(nf,300));
      if(red(fx(e, x0)) == 0, next);
      if(vecmax(concat(vector(3, i, lvec(nf, red((x0 - e[i])/b[i]), B, pl)))) != 0, next);
      my(v = caslocal(nf, e, b, bp, G, pl, x0, N)); if(v < 0, next); listput(vals, v));
    if(#vals == 0, error("cassels2: no local point at ", pl));
    if(vecmin(Vec(vals)) != vecmax(Vec(vals)), print("  INCONSISTENT at ", if(isre(pl), "inf", pl[1].p), " ", Vec(vals)));
    tot += vals[1]; listput(per, vals[1]));
  Z = vector(3, i, dengen(bnf, G[i]));
  for(i=1, 3, for(k=1, #PL, corr += (1 - hil(nf, Z[i], bp[i], PL[k]))/2));
  [(tot + corr) % 2, Vec(per), corr % 2, Z]; }
