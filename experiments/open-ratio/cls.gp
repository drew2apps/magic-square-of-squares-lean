K=nfinit(y^2+1);
f(z)=idealfactor(K,z);
print(factor(1095961765000)); print(factor(1201132190341915225000000));
print(f(2047701738184*y - 3586214119488));
print(factor(norm(Mod(2047701738184*y - 3586214119488,y^2+1))));
