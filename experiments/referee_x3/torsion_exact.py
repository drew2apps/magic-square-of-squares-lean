"""x3: torsion_E4.gp with the heuristic square test (300 split primes) replaced by an exact one (nfroots of x^2 - a over
Q(i)); data copied from the committed ctp/torsion_E4.gp.  Run through cypari2 (PARI 2.17.4)."""
import cypari2
pari = cypari2.Pari(); pari.allocatemem(10**9)
src = open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', 'open-ratio', 'ctp', 'torsion_E4.gp')).read()
lines = [l for l in src.split('\n') if l.strip()]
# keep the data lines and helper definitions; replace the final loop by an exact test
for l in lines[:9]:
    pari(l.replace('{', '').replace('}', '') if l.startswith('{') else l)
pari('isq(a) = (R(a) == 0) || #nfroots(nf, x^2 - R(a)) > 0')
pari('for(i=1,3, my(j = (i%3)+1, k = ((i+1)%3)+1, a = sq(e[i]-e[j]), b = sq(e[i]-e[k])); print("T", i, ": ei-ej square ", a!=0, ", ei-ek square ", b!=0); if(a && b, my(x0 = R(e[i] + a*b)); listput(gens, img(x0)); print("  4-torsion x = ", x0)))')
pari('foreach(gens, v, print("image ", v); forvec(c = vector(4, i, [0,1]), my(w = v); for(i=1,4, if(c[i], w = [w[1]*S[i][1], w[2]*S[i][2]])); my(h = issq(w[1]) && issq(w[2]), x = isq(w[1]) && isq(w[2])); if(h || x, print("  = Selmer combination ", c, "  heuristic ", h, "  exact ", x))))')
