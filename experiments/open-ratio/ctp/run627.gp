\\ Cassels-Tate pairing on Sel^2 of E1 and E4 for (627,364) over Q(i).  EXPERIMENT; run: gp -q -s 4G run627.gp
read("/mnt/project-files/mss-experiments/open-ratio/ctp/ctp.gp");
PLc(nf, e, extra) = { my(L = List(badplaces(nf, e)), have = Set(apply(pl -> if(isre(pl), 0, pl[1].p), badplaces(nf, e))));
  foreach(extra, p, if(!setsearch(have, p), my(D = idealprimedec(nf, p)); foreach(D, pr, listput(L, [pr])))); Vec(L); }
casmat(bnf, e) = { my(r = selmer2(bnf, e), S = r[1], n = #S, M = matrix(n, n), PL = PLc(bnf.nf, e, [3,5,7]));
  print("  Selmer elements ", n, ": ", S);
  for(k=1, n, for(l=1, n, my(b = [S[k][1], S[k][2], S[k][1]*S[k][2]], bp = [S[l][1], S[l][2], S[l][1]*S[l][2]], c);
     my(vals = List());
     for(rep = 1, 2, for(att = 1, 4, c = iferr(cassels2(bnf, e, b, bp, PL, 60), E, print("  retry: ", E); 0); if(c != 0, break)); if(c == 0, error("gave up"));
        listput(vals, c[1]); print("  <", k, ",", l, "> run ", rep, " = ", c[1], "  per place ", c[2], " correction ", c[3]));
     if(vals[1] != vals[2], print("  DISAGREE at <", k, ",", l, ">")); M[k,l] = vals[1]));
  [n, M, matrank(Mod(M,2)), M == M~ && vecmax(vector(n,i,M[i,i]))==0, S]; }
bk = bnfinit(y^2+1, 1); nf_pol_glob = bk.nf.pol;
E1 = [-57020325625 - 39655832216*y, 216970115625 - 39655832216*y, -159949790000 + 79311664432*y];
E4 = [-176651052841 - 288326974936*y, 97339388409 + 128377184936*y, 79311664432 + 159949790000*y];
{foreach([["E4", E4], ["E1", E1]], c, my(t = getabstime(), r);
  print("=== ", c[1], " e = ", c[2]);
  r = casmat(bk, c[2]);
  print(c[1], ": dim Sel2 = ", r[1], ", pairing matrix ", r[2], ", rank ", r[3], ", alternating ", r[4], ", ", (getabstime()-t)\1000, " s"));}
