\\ E4 only: upper-triangle entries of the Cassels-Tate pairing on Sel^2 of E4 for (627,364) over Q(i), each entry written
\\ to run627b.log as soon as it is computed.  EXPERIMENT; run: gp -q -s 4G run627b.gp
read("/mnt/project-files/mss-experiments/open-ratio/ctp/ctp.gp");
PLc(nf, e, extra) = { my(L = List(badplaces(nf, e)), have = Set(apply(pl -> if(isre(pl), 0, pl[1].p), badplaces(nf, e))));
  foreach(extra, p, if(!setsearch(have, p), my(D = idealprimedec(nf, p)); foreach(D, pr, listput(L, [pr])))); Vec(L); }
LOG = "/mnt/project-files/mss-experiments/open-ratio/ctp/run627b.log";
bk = bnfinit(y^2+1, 1); nf_pol_glob = bk.nf.pol;
E4 = [-176651052841 - 288326974936*y, 97339388409 + 128377184936*y, 79311664432 + 159949790000*y];
{ my(r = selmer2(bk, E4), S = r[1], n = #S, PL = PLc(bk.nf, E4, [3,5,7]), ord = [[1,2],[2,4],[2,3],[3,4],[1,3],[1,4]]);
  write(LOG, "E4 Selmer elements ", n, ": ", S);
  \\ which elements are torsion images: the image of E(K)[2]-related torsion is computed by selmer2's second output if present
  write(LOG, "selmer2 extra output: ", if(#r > 1, r[2..#r], "none"));
  foreach(ord, kl, my(k = kl[1], l = kl[2], b = [S[k][1], S[k][2], S[k][1]*S[k][2]], bp = [S[l][1], S[l][2], S[l][1]*S[l][2]], c, vals = List());
    for(rep = 1, 2, for(att = 1, 4, c = iferr(cassels2(bk, E4, b, bp, PL, 60), E, write(LOG, "  retry: ", E); 0); if(c != 0, break));
      if(c == 0, write(LOG, "  <", k, ",", l, "> gave up"); break);
      listput(vals, c[1]); write(LOG, "  <", k, ",", l, "> run ", rep, " = ", c[1], "  per place ", c[2], " correction ", c[3]));
    if(#vals == 2 && vals[1] != vals[2], write(LOG, "  DISAGREE at <", k, ",", l, ">"))); }
write(LOG, "done");
