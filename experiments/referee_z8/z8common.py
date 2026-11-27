"""Shared helpers for the review z8 of roadmap H (written by an AI model; not read by a person; experiment, not counted).

Paths: the repository is MSS_REPO, or else ../.. (this folder is meant to sit at experiments/referee_z8); G1's folder is
G1_DIR or else <repo>/experiments/roadmap-g/g1; H1's and H3's are H1_DIR, H3_DIR or else <repo>/experiments/roadmap-h/h1,
h3.  emit(lines, name): without --write print the lines; with --write also write <name> in this folder.
load_phantoms(): the 24 ratios of H1 (the scan's 13, (1109503,597671), the partners) with H1's known and extra tuples,
read from H1's committed h1_scan_k*.txt and h1_partners.txt (tuples in tuple_descent's order d0, d8, dL, dU)."""
import os, sys, re, ast, glob
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.environ.get('MSS_REPO', os.path.join(HERE, '..', '..')))
G1 = os.environ.get('G1_DIR', os.path.join(REPO, 'experiments', 'roadmap-g', 'g1'))
H1 = os.environ.get('H1_DIR', os.path.join(REPO, 'experiments', 'roadmap-h', 'h1'))
H3 = os.environ.get('H3_DIR', os.path.join(REPO, 'experiments', 'roadmap-h', 'h3'))
os.environ.setdefault('MSS_REPO', REPO)
WRITE = '--write' in sys.argv
ARGS = [a for a in sys.argv[1:] if a != '--write']

def emit(lines, name):
    s = '\n'.join(lines) + '\n'
    sys.stdout.write(s); sys.stdout.flush()
    if WRITE: open(os.path.join(HERE, name), 'w').write(s)

def isprime(p): return p > 1 and all(p % q for q in range(2, int(p ** 0.5) + 1))

def good(m, n, p):
    return all(x % p for x in (m, n, m - n, m + n, m * m + n * n, m ** 8 + 14 * m ** 4 * n ** 4 + n ** 8))

CONTROL = {(1109503, 597671): ([(1109503, 597671, 2, 1), (1195342, 2219006, 2, 1)],
           [(1109503, 597671, 2, 959), (1109503, 597671, 1918, 1), (1195342, 2219006, 2, 959), (1195342, 2219006, 1918, 1)])}

def load_phantoms():
    out = dict(CONTROL)
    for f in sorted(glob.glob(os.path.join(H1, 'h1_scan_k*.txt'))):
        for l in open(f):
            mm = re.match(r'\((\d+),(\d+)\) PHANTOM: .*known (\[.*\]); extra (\[.*\])', l.strip())
            if mm: out[(int(mm.group(1)), int(mm.group(2)))] = (ast.literal_eval(mm.group(3)), ast.literal_eval(mm.group(4)))
    for l in open(os.path.join(H1, 'h1_partners.txt')):
        mm = re.match(r'\((\d+),(\d+)\) -> partner \((\d+),(\d+)\).*: PHANTOM; .*known (\[.*\]); extra (\[.*\])', l.strip())
        if mm:
            b = (int(mm.group(3)), int(mm.group(4)))
            if b not in out: out[b] = (ast.literal_eval(mm.group(5)), ast.literal_eval(mm.group(6)))
    assert len(out) == 24, len(out)
    return out

# the 7 partner pairs whose extra tuples survive every prime (full_descent.txt), as 14 ratios
TRUE_PHANTOMS = [(204724, 97893), (302617, 106831), (413004, 148603), (561607, 264401), (756492, 255635), (1012127, 500857),
                 (853587, 255916), (1109503, 597671), (1992873, 999664), (2992537, 993209), (189089400, 156063271),
                 (345152671, 33026129), (231641256, 72371617), (304012873, 159269639)]

def pari():
    import cypari2
    P = cypari2.Pari()
    sys.stdout.flush(); fd = os.dup(1); nul = os.open(os.devnull, os.O_WRONLY); os.dup2(nul, 1)
    try: P.allocatemem(2 * 10 ** 9)          # PARI prints its new stack size on stdout: kept out of the output
    finally: os.dup2(fd, 1); os.close(fd); os.close(nul)
    return P

def fac_proved(P, x):
    """the prime factors of |x| (PARI's factor), each proved prime by PARI's isprime"""
    x = abs(int(x))
    if x <= 1: return []
    f = P.factor(x); ps = [int(q) for q in f[0]]
    for q in ps: assert P.isprime(q), ('not proved prime', q)
    return sorted(ps)
