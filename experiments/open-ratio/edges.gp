\\ For every 2-isogeny edge in the 2-power K-isogeny class: integral model y^2 = x(x^2 + a x + b) over Z[i] with the kernel at 0.
read("isog_627_364.gp"); bnf=bnfinit(y^2+1,1);
toab(C) = { my(b2=C.b2, b4=C.b4, b6=C.b6, f = x^3 + b2*x^2 + 8*b4*x + 16*b6, R = nfroots(bnf, lift(f)), out=List());
  for(k=1,#R, my(g = lift(subst(f, x, x + R[k]))); \\ g = x^3 + a x^2 + b x
     listput(out, [R[k], polcoef(g,2), polcoef(g,1)]) ); Vec(out) };
{for(k=1,#Curves, my(C=ellinit(Curves[k],bnf), L=toab(C)); for(j=1,#L, my(a=Mod(L[j][2],y^2+1), b=Mod(L[j][3],y^2+1), iso=ellinit([0,-2*a,0,a^2-4*b,0],bnf), jj = ellinit(iso,bnf).j, tgt=0);
   for(l=1,#Curves, if(ellinit(Curves[l],bnf).j == jj, tgt=l));
   print("EDGE ", k, " ", j, " ", tgt, " ", [real(subst(lift(a),y,I)), imag(subst(lift(a),y,I)), real(subst(lift(b),y,I)), imag(subst(lift(b),y,I))])))}
