read("gen4.gp");
\\ test curves over Q (Tate normal form, Z/4): rank over K = rank E(Q) + rank E^(-1)(Q)
{for(t = 2, 40, my(E = ellinit([1,-t,-t,0,0]), Et = elltwist(E, -4), r1 = ellrank(E), r2 = ellrank(ellinit(Et)));
   if(r1[1]==r1[2] && r2[1]==r2[2], print("t=", t, " rank E(Q)=", r1[1], " rank E^-1(Q)=", r2[1], "  rank E(K)=", r1[1]+r2[1], "  cond ", ellglobalred(E)[1], " ", ellglobalred(ellinit(Et))[1])));}
