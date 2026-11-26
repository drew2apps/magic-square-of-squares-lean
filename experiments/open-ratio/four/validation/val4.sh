# Validation of the cyclic 4-isogeny descent over Q(i) on test curves of known rank (base change from Q, Tate normal form, Z/4).
# Output: val4.log in this folder.
cd "$(dirname "$0")/.."
L=validation/val4.log; : > $L
echo "== known ranks: rank E(K) = rank E(Q) + rank E^-1(Q) (PARI ellrank, lower = upper)" >> $L
gp -q -s 1G test_t.gp >> $L 2>&1
for t in 2 3 7 15 25; do echo "== t=$t: run4.py test_t$t.gp" >> $L; timeout 3000 python3 run4.py test_t$t.gp "" 1 >> $L 2>&1; done
echo done >> $L
