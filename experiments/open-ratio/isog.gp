\\ The K-isogeny graph of E (2-power isogenies with K-rational kernels), and roots of the 2-division cubic for each curve.
read("curve_627_364.gp");
bnf = bnfinit(y^2+1,1); E = ellinit(Ecoef, bnf);
T = elltors(E); print("torsion ", T[1..2]); print("gens ", lift(T[3]));
\\ cubic roots
roots2(Ee) = { my(f = x^3 + Ee.a2*x^2 + Ee.a4*x + Ee.a6); \\ assumes a1=a3=0
  nfroots(bnf, lift(f)) };
print("2-division roots of E: ", lift(roots2(E)));
\\ isogenous curves by 2-torsion points
curves = List([E]); 
q = [E]; seen = [E.j];
{while(#q, my(C = q[1]); q = q[2..#q]; 
  my(F = lift(subst(x^3 + C.a2*x^2 + C.a4*x + C.a6, y, y)), R = nfroots(bnf, F));
  for(k=1,#R, my(iso = ellisogeny(C, [R[k], 0], 1), C2 = ellinit(iso, bnf), C2m = ellminimalmodel(C2), j2 = C2m.j);
     if(#select(z->z==j2, seen)==0, seen = concat(seen, [j2]); q = concat(q, [C2m]); listput(curves, C2m))));}
print("number of curves in the 2-power K-isogeny class found: ", #curves);
for(k=1,#curves, my(C=curves[k]); print(k, ": a-invs ", lift(C[1..5]), "  tors ", elltors(C)[1..2], "  #2-div roots in K ", #nfroots(bnf, lift(x^3 + C.a2*x^2 + C.a4*x + C.a6))));
write("isog_627_364.gp", "Curves = ", apply(C->lift(C[1..5]), Vec(curves)), ";");
