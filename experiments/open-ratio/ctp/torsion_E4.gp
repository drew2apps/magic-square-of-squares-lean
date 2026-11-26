nf = nfinit(y^2+1); R(z) = lift(Mod(z, y^2+1));
e = [-176651052841 - 288326974936*y, 97339388409 + 128377184936*y, 79311664432 + 159949790000*y];
S = [[991*y + 991, y^2 + y], [-4*y^3 + 8*y^2 + 5*y, -2*y^2 - 9*y - 10], [6710*y^4 + 94299*y^3 + 367434*y^2 + 259200*y, 1764730*y^4 + 24800637*y^3 + 96635142*y^2 + 68169600*y], [95732*y - 621357, 25177516*y - 163416891]];
S = apply(v -> apply(R, v), S);
PR = List(); forprime(p = 1000, 200000, if(p%4==1 && #PR < 300, listput(PR, p)));
issq(a) = { my(z = R(a)); foreach(PR, p, my(r = sqrt(Mod(-1,p))); foreach([r, -r], t, my(v = subst(z, y, t)); if(v != 0 && !issquare(v), return(0)))); 1; }
sq(a) = { my(r = nfroots(nf, x^2 - R(a))); if(#r, lift(r[1]), 0); }
img(x0) = if(x0 == e[1], [R((e[1]-e[2])*(e[1]-e[3])), R(e[1]-e[2])], x0 == e[2], [R(e[2]-e[1]), R((e[2]-e[1])*(e[2]-e[3]))], [R(x0-e[1]), R(x0-e[2])]);
gens = List([img(e[1]), img(e[2])]);
{for(i=1,3, my(j = (i%3)+1, k = ((i+1)%3)+1, a = sq(e[i]-e[j]), b = sq(e[i]-e[k]));
  print("T", i, ": ei-ej square ", a!=0, ", ei-ek square ", b!=0);
  if(a && b, my(x0 = R(e[i] + a*b)); listput(gens, img(x0)); print("  4-torsion x = ", x0)));}
{foreach(gens, v, print("image ", v);
  forvec(c = vector(4, i, [0,1]), my(w = v); for(i=1,4, if(c[i], w = [w[1]*S[i][1], w[2]*S[i][2]])); if(issq(w[1]) && issq(w[2]), print("  = Selmer combination ", c))));}
