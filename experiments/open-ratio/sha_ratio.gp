\\ Isogeny invariance of the BSD quotient over K = Q(i) across the 2-power isogeny class.
\\ If rank E_j(K) = 0 for all j (they are isogenous, so all or none) and Sha finite, then
\\ |Sha(E_j)| * Omega_j * c_j / |tors_j|^2 is the same for all j (Cassels; K has class number 1,
\\ so global minimal models exist and the unit scaling of omega has absolute value 1).
\\ Omega_j here = 2*area of the period lattice at the complex place (constant factors cancel in ratios).
read("isog_627_364.gp");
bnf = bnfinit(y^2+1,1);
default(realprecision, 60);
res = vector(#Curves);
{for(j=1,#Curves,
  my(E = ellinit(Curves[j], bnf), M = ellminimalmodel(E), gr, c, t, EC, w, area);
  gr = ellglobalred(M); c = gr[3];
  t = elltors(M)[1];
  EC = ellinit(apply(z->subst(lift(z), y, I), M[1..5]));
  w = EC.omega; area = abs(imag(conj(w[1])*w[2]));
  res[j] = [t, c, 2*area];
  print(j, ": tors ", t, "  Tamagawa product ", c, "  Omega ", 2*area, "  disc norm ", idealnorm(bnf, M.disc));
);}
\\ T_j = tors^2/(Omega*c);  |Sha_j| proportional to T_j
T = vector(#Curves, j, res[j][1]^2/(res[j][3]*res[j][2]));
print("ratios |Sha_j|/|Sha_1| = ", vector(#Curves, j, T[j]/T[1]));
print("log2 of ratios = ", vector(#Curves, j, log(T[j]/T[1])/log(2)));
