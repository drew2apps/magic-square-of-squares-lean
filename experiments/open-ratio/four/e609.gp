E0 = ellinit(ellminimalmodel(ellinit([0, 58077506951413592707/100, 0, 1150906168103852495278998427225573449283/10000, 7723823874569200669382970534520117704028821195707960451009/1000000])));
print("E0 min model ", E0[1..5]); print("cond ", ellglobalred(E0)[1]); print("tors Q ", elltors(E0)[1]);
print("ellrank E0 ", ellrank(E0)[1..3]);
Et = ellinit(ellminimalmodel(ellinit(elltwist(E0, -4)))); print("twist min ", Et[1..5]);
print("ellrank E0^-1 ", ellrank(Et)[1..3]);
print("isogeny class over Q: ", #ellisomat(E0)[1]);
