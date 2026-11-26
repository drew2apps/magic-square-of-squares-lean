read("gen4.gp"); read("classK.gp");
Cl = classK([0, 1, 0, 2657389687335443426815780334514559, -46004163125583148363653912915452814645634607488641]);
print("K-class size ", #Cl);
for(k=1,#Cl, my(C=ellinit(Cl[k],bnf)); print(k, " tors ", elltors(C)[2]));
write("class609_K.gp", "Curves = ", Cl, ";");
print("4-isogenies: ", gen4(Cl, "e609_4data.gp"));
