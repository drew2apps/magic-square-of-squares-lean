import cypari2
pari = cypari2.Pari(); pari.allocatemem(4*10**8)
# E: Y^2 = z(z-1)(z+1)(z-2).  z = 2 + 1/t ... simpler: z = 1/t:  (t^2 Y)^2 = (1-t)(1+t)(1-2t) = 2t^3 - t^2 - 2t + 1 ; X = 2t, W = 2 t^2 Y:  W^2 = X^3 - X^2 - 4X + 4 = (X-1)(X-2)(X+2)
print(pari('nf = nfinit(y^2+1); E = ellinit([0,-1,0,-4,4], nf); T = elltors(E); print("torsion over Q(i): ", T[1], " ", T[2]); T[1]'))
print(pari('EQ = ellinit([0,-1,0,-4,4]); print("over Q: ", elltors(EQ)); print("rank ", ellrank(EQ)); Em = ellinit([0,1,0,-4,-4]); print("twist -1 over Q: ", elltors(Em), " rank ", ellrank(Em)); 0'))
# all torsion points over Q(i): enumerate group generated
print(pari('''G = T[3]; ord = T[2]; pts = List();
for(a=0, ord[1]-1, for(b=0, if(#ord>1, ord[2]-1, 0), P = elladd(E, ellmul(E, G[1], a), if(#ord>1, ellmul(E, G[2], b), [0])); listput(pts, P)));
for(k=1, #pts, P = pts[k]; if(P == [0], print("O  -> X = oo, t = oo, z = 0"), X = lift(nfbasistoalg(nf, P[1])); print("X = ", X, "   z = 2/X = ", if(X==0, "oo", lift(Mod(2,y^2+1)/Mod(X,y^2+1))))));
#pts'''))
