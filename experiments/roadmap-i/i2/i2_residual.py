"""I2, the residual tuples (EXPERIMENT, written by an AI model; not read by a person; nothing here changes the census).
For each ratio of a list that keeps a class and needs a good prime: the tuples (d0, d8, dL, dU) allowed at every bad
prime (the repository's exact local test, every bad prime; the archimedean place of Q(i) is complex and imposes
nothing), and for each such tuple the first good prime (in increasing order, to 5000) at which it is not allowed.
The killing prime is the largest of these: the first good prime at which no residual tuple is left.  A tuple that
no good prime <= 5000 kills is printed with kill null (a survivor candidate; none expected).
Output: one JSON object per line, ratios in the order of the list:
  {"m": m, "n": n, "residual": [[d0, d8, dL, dU], ...], "kill": [p, ...], "killing_prime": P}
(residual sorted; kill[j] is the prime that kills residual[j]).  This is the machine-readable file asked for in the
roadmap I addendum, for testing whether a Cassels-Tate pairing on a genus-1 quotient removes the same tuples.
usage: i2_residual.py LIST i N [--write]     worker i of N; output i2_residual_<i>of<N>.jsonl (.part while running)
       i2_residual.py --merge N [--write]     the N workers' files in the order of the list -> i2_residual.jsonl
       i2_residual.py --negative [--write]    forged records refused, the true one accepted -> i2_residual_negative.txt
Folders: G1_DIR, G2_DIR (the repository's experiments/roadmap-g), H2_DIR, H4_DIR (roadmap H's h2, h4)."""
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__))
H2 = os.environ.get('H2_DIR', next((d for d in (os.path.join(HERE, '..', '..', 'roadmap-h', 'h2'), os.path.join(HERE, '..', 'h2')) if os.path.isdir(d)), os.path.join(HERE, '..', '..', 'roadmap-h', 'h2'))); sys.path.insert(0, H2); sys.path.insert(0, HERE)
argv = sys.argv; sys.argv = [argv[0]]
import h2_hunt as HH
import i2_shape as I2
sys.argv = argv
P = HH.P; X = HH.X; K = HH.K; H4 = HH.H4
LIST = 'i2_needgood_list.txt'

def residual(m, n):
    forms, sup, known = K.TD.magic(m, n)
    divs = K.divisors(sup)
    bad = sorted(HH.primes_of(2*m*n*(m*m - n*n), m*m + n*n, m**8 + 14*m**4*n**4 + n**8))
    for q in bad: assert P.isprime(q)
    locs = []
    for p in bad:
        L, al = X.allowed_at(forms, p)
        for kt in known: assert tuple(L.dlabel(d) for d in kt) in al, ('known point rejected', m, n, p)
        if len(al) < len(L.labels()) ** len(forms): locs.append((p, L, al))
    S = H4.chunked(divs, locs, False)
    kill = [None] * len(S); left = set(range(len(S)))
    for q in P.primes([3, 5000]):
        q = int(q)
        if not left: break
        if q in bad: continue
        L, al = X.allowed_at(forms, q)
        for j in sorted(left):
            if tuple(L.dlabel(d) for d in S[j]) not in al: kill[j] = q; left.discard(j)
    Pk = 0 if not S else None if left else max(kill)        # 0: the bad primes alone leave nothing
    return dict(m=m, n=n, residual=[list(t) for t in S], kill=kill, killing_prime=Pk)

def line(r): return json.dumps(r, separators=(', ', ': '))

def verify(r):
    """True when the record agrees with a recomputation (residual set, each tuple's killing prime, the killing prime)."""
    return residual(r['m'], r['n']) == r

def negative():
    out = ['I2 residual tuples, negative and positive controls. Written by an AI model; not read by a person.']
    m, n = 15001, 2081                                       # the first ratio of the list (|S| = 2, killed at 29)
    r = residual(m, n)
    out.append('true record for (%d,%d): %d residual tuples, killing prime %s: %s' % (m, n, len(r['residual']), r['killing_prime'], 'accepted' if verify(r) else 'REFUSED'))
    assert verify(r)
    f = json.loads(json.dumps(r)); f['residual'][0][3] *= 3   # a tuple outside the set left by the bad primes
    out.append('forged residual tuple %s: %s' % (f['residual'][0], 'refused' if not verify(f) else 'ACCEPTED')); assert not verify(f)
    f = json.loads(json.dumps(r)); f['killing_prime'] = 23; f['kill'] = [min(k, 23) for k in f['kill']]
    out.append('forged killing prime 23: %s' % ('refused' if not verify(f) else 'ACCEPTED')); assert not verify(f)
    f = json.loads(json.dumps(r)); f['residual'] = f['residual'][:1]; f['kill'] = f['kill'][:1]; f['killing_prime'] = f['kill'][0]
    out.append('forged record with a residual tuple left out: %s' % ('refused' if not verify(f) else 'ACCEPTED')); assert not verify(f)
    return out

if __name__ == '__main__':
    if '--negative' in argv:
        s = '\n'.join(negative()) + '\n'; f = os.path.join(HERE, 'i2_residual_negative.txt'); part = None
    elif '--merge' in argv:
        N = int([x for x in argv[1:] if not x.startswith('--')][0])
        files = [[l for l in open(os.path.join(HERE, 'i2_residual_%dof%d.jsonl' % (i, N)))] for i in range(N)]
        rows = []; j = 0
        while any(files):
            rows.append(files[j % N].pop(0)); j += 1
        s = ''.join(rows); f = os.path.join(HERE, 'i2_residual.jsonl'); part = None
    else:
        a = [x for x in argv[1:] if not x.startswith('--')]
        lst, i, N = a[0], int(a[1]), int(a[2])
        rs = [r for j, r in enumerate(I2.read(os.path.join(HERE, lst))) if j % N == i]
        f = os.path.join(HERE, 'i2_residual_%dof%d.jsonl' % (i, N)); part = f + '.part'
        done = [l.rstrip('\n') for l in open(part)] if '--write' in argv and os.path.exists(part) else []
        pf = open(part, 'a') if '--write' in argv else None
        out = []
        for j, (m, n) in enumerate(rs):
            if j < len(done): out.append(done[j]); continue
            try: s = line(residual(m, n))
            except MemoryError: s = line(dict(m=m, n=n, status='over the per-process memory limit, not done'))
            out.append(s)
            if pf: pf.write(s + '\n'); pf.flush()
        s = '\n'.join(out) + '\n'
    if '--write' in argv:
        open(f, 'w').write(s)
        if part and os.path.exists(part): os.remove(part)
    else:
        sys.stdout.write(s)
        if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS from ' + os.path.basename(f))
