read("isog_627_364.gp"); bnf=bnfinit(y^2+1,1);
for(k=1,#Curves, my(C=ellinit(Curves[k],bnf)); if(C.a1==0 && C.a3==0, print(k, " ", lift(nfroots(bnf, lift(x^3+C.a2*x^2+C.a4*x+C.a6))))));
