\\ Validation over Q(i) on curves base-changed from Q: pairing rank vs dim Sel2 - 2 - (rank E(Q) + rank E^-1(Q)).
\\ Run: gp -q -s 3G val_Qi.gp   Output: val_Qi.log
read("/mnt/project-files/mss-experiments/open-ratio/ctp/ctp.gp");
LOG = "/mnt/project-files/mss-experiments/open-ratio/ctp/validation/val_Qi.log";
PLc(nf, e, extra) = { my(L = List(badplaces(nf, e)), have = Set(apply(pl -> if(isre(pl), 0, pl[1].p), badplaces(nf, e))));
  foreach(extra, p, if(!setsearch(have, p), my(D = idealprimedec(nf, p)); foreach(D, pr, listput(L, [pr])))); Vec(L); }
casmat(bnf, e) = { my(r = selmer2(bnf, e), S = r[1], n = #S, M = matrix(n, n), PL = PLc(bnf.nf, e, [3,5,7]));
  for(k=1, n, for(l=1, n, my(b = [S[k][1], S[k][2], S[k][1]*S[k][2]], bp = [S[l][1], S[l][2], S[l][1]*S[l][2]]);
     M[k,l] = cassels2(bnf, e, b, bp, PL, 60)[1]));
  [n, M, matrank(Mod(M,2)), M == M~ && vecmax(vector(n,i,M[i,i]))==0, S]; }
bk = bnfinit(y^2+1, 1); nf_pol_glob = bk.nf.pol;
write(LOG, "over Q(i), curves from Q: [e] | dimSel2-2 | pairing rank | alternating | rank E(Q)+rank E^-1(Q) [lo,hi] | Sha2 lower bounds over Q for E, E^-1");
{foreach([[0,195,-376],[0,333,-238],[0,360,-41],[0,145,-109],[0,240,-34],[0,34,-34],[-7,3,11],[0,41,-41]], e,
  my(E = ellinit([0, -(e[1]+e[2]+e[3]), 0, e[1]*e[2]+e[1]*e[3]+e[2]*e[3], -e[1]*e[2]*e[3]]), Et = elltwist(E, -4), c = casmat(bk, e), r1 = ellrank(E), r2 = ellrank(Et));
  write(LOG, e, " | ", c[1]-2, " | ", c[3], " | ", c[4], " | [", r1[1]+r2[1], ",", r1[2]+r2[2], "] | ", r1[3], " ", r2[3]); write(LOG, "   matrix ", c[2]));}

write(LOG, "done");
