\\ Validation over Q: direct pairing (cassels) vs reciprocity shortcut (cassels2), and pairing rank vs dim Sel2 - 2 - rank.
\\ Run: gp -q -s 2G val_Q.gp   Output: val_Q.log
read("/mnt/project-files/mss-experiments/open-ratio/ctp/ctp.gp");
LOG = "/mnt/project-files/mss-experiments/open-ratio/ctp/validation/val_Q.log";
PLc(nf, e, extra) = { my(L = List(badplaces(nf, e)), have = Set(apply(pl -> if(isre(pl), 0, pl[1].p), badplaces(nf, e))));
  foreach(extra, p, if(!setsearch(have, p), my(D = idealprimedec(nf, p)); foreach(D, pr, listput(L, [pr])))); Vec(L); }
bq = bnfinit(y, 1); nf_pol_glob = bq.nf.pol;
write(LOG, "[e] | dim Sel2 | direct == shortcut | pairing rank | alternating | ellrank [lo,hi,Sha2 lower] | dimSel2-2-rank");
{foreach([[0,195,-376],[0,360,-41],[0,240,-283],[0,145,-109],[0,333,-238],[0,240,-34],[0,34,-34],[-7,3,11],[0,41,-41]], e,
  my(r = selmer2(bq, e), S = r[1], n = #S, M1 = matrix(n,n), M2 = matrix(n,n), PL = PLc(bq.nf, e, [3,5,7]),
     E = ellinit([0, -(e[1]+e[2]+e[3]), 0, e[1]*e[2]+e[1]*e[3]+e[2]*e[3], -e[1]*e[2]*e[3]]), rk = ellrank(E));
  for(k=1,n, for(l=1,n, my(b = [S[k][1], S[k][2], S[k][1]*S[k][2]], bp = [S[l][1], S[l][2], S[l][1]*S[l][2]]);
    M1[k,l] = cassels(bq, e, b, bp, PL, 60)[1]; M2[k,l] = cassels2(bq, e, b, bp, PL, 60)[1]));
  write(LOG, e, " | ", n, " | ", M1 == M2, " | ", matrank(Mod(M2,2)), " | ", M2 == M2~ && vecmax(vector(n,i,M2[i,i]))==0, " | ", rk[1..3], " | ", if(rk[1]==rk[2], n-2-rk[1], "?"));
  write(LOG, "   matrix ", M2));}
write(LOG, "done");
