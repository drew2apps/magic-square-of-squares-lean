\\ classK(coeffs): the 2-power K-isogeny class (K = Q(i)) by BFS over K-rational 2-torsion points.  needs bnf from gen4.gp
classK(c) = {
  my(E = ellminimalmodel(ellinit(c, bnf)), curves = List([E]), q = [E], seen = [E.j]);
  while(#q, my(C = q[1]); q = q[2..#q];
    my(T = select(P->ellorder(C,P)==2, alltors(C)));
    for(k=1,#T, my(C2 = ellminimalmodel(ellinit(ellisogeny(C, T[k], 1), bnf)));
       if(#select(z->z==C2.j, seen)==0, seen = concat(seen,[C2.j]); q = concat(q,[C2]); listput(curves, C2))));
  apply(C->lift(C[1..5]), Vec(curves)) };
