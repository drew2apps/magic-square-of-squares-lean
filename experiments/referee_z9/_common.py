"""Review z9 of roadmap H (run 7), written by an AI model; not read by a person.  Shared helpers for the scripts of this
folder: paths (the repository root is ../.., or MSS_REPO), output collected and printed at the end, PARI's messages on
stdout sent to stderr, and the comparison with the committed .txt (printed to stderr: '# matches NAME.txt')."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get('MSS_REPO', os.path.join(HERE, '..', '..'))
EXP = os.path.join(REPO, 'experiments')
H2 = os.path.join(EXP, 'roadmap-h', 'h2'); H4 = os.path.join(EXP, 'roadmap-h', 'h4'); H6 = os.path.join(EXP, 'roadmap-h', 'h6')
G1 = os.path.join(EXP, 'roadmap-g', 'g1'); G2 = os.path.join(EXP, 'roadmap-g', 'g2'); G3 = os.path.join(EXP, 'roadmap-g', 'g3')
STRUCT = os.path.join(REPO, 'scripts', 'cover_descent', 'structure')

class Out:
    """lines are kept and printed by finish(); while the computation runs, file descriptor 1 goes to stderr, so that
    PARI's own messages ('PARI stack size set to ...') never reach the output."""
    def __init__(self, name):
        self.name = name; self.lines = []
        sys.stdout.flush(); self.fd = os.dup(1); os.dup2(2, 1)
    def __call__(self, s=''): self.lines.append(s)
    def finish(self):
        sys.stdout.flush(); os.dup2(self.fd, 1); os.close(self.fd)
        s = '\n'.join(self.lines) + '\n'; f = os.path.join(HERE, self.name)
        if '--write' in sys.argv: open(f, 'w').write(s)
        sys.stdout.write(s); sys.stdout.flush()
        if '--write' not in sys.argv:
            ok = os.path.exists(f) and open(f).read() == s
            sys.stderr.write(('# matches %s\n' if ok else '# DIFFERS from %s\n') % self.name)
            sys.exit(0 if ok else 1)
