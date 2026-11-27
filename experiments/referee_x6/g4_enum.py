"""Review x6 (run 6), written by an AI model; not read by a person.  Does g4_search.points list every small combination of
the generators (not only the generators)?  Prints only.  usage: python3 -I g4_enum.py   (about 1 minute; PARI-dependent)"""
import sys, os, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'roadmap-g', 'g4'))
import g4_search as S
pari = S.pari
tested = missing = 0
for m, n in [(12, 5), (40, 9), (195, 28), (440, 279), (7735, 4968), (56, 33), (84, 13)]:
    for name, E, to_t, Xs in S.curve_list(m, n):
        rk = pari.ellrank(E, 1); lo, hi = int(rk[0]), int(rk[1]); gens = list(rk[3])
        if not gens or len(gens) != lo or lo > 3: continue
        gens = list(pari.ellsaturation(E, gens, 100))
        hmax = max(float(pari.ellheight(E, P)) for P in gens); H = 64 * hmax
        pts, vecs = S.points(E, gens, H)
        have = set(str(P) for P in pts)
        r = len(gens)
        for v in itertools.product(range(-4, 5), repeat=r):
            Q = pari('[0]')
            for a, g in zip(v, gens): Q = pari.elladd(E, Q, pari.ellmul(E, g, a))
            h = float(pari.ellheight(E, Q)) if Q != pari('[0]') else 0.0
            if h <= 0.999 * H:
                tested += 1; missing += str(Q) not in have
print('combinations sum a_i g_i, |a_i| <= 4, of height <= 0.999 H on the curves searched: %d; not listed: %d' % (tested, missing))
