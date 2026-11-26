"""Full 2-Selmer bound over Q(i) for the two curves with full 2-torsion in the 2-power K-isogeny class of the (627,364)
bound-2 curve (L1 L3 M1 M2), using the repository's kdescent.py unchanged (read-only import)."""
import sys, time
sys.path.insert(0, "/tmp/claude-0/-home-user/9ec7a8d3-42dc-53c8-8518-1a663ac809f1/scratchpad/repo/scripts/second_descent")
import kdescent as kd
curves = {
 1: [(-57020325625, -39655832216), (216970115625, -39655832216), (-159949790000, 79311664432)],
 4: [(-176651052841, -288326974936), (97339388409, 128377184936), (79311664432, 159949790000)],
}
for k, e in curves.items():
    for seed in (1, 2):
        t0 = time.time(); b = kd.selmer_bound(e, seed=seed, verbose=True)
        print("curve", k, "seed", seed, "dim Sel - 2 =", b, " (%.0fs)" % (time.time() - t0), flush=True)
