"""Paths and the output helper of referee_z11 (review agent z11; written by an AI model; not read by a person).
Paths relative to this folder: repository root ../.. (or MSS_REPO), t10's folder ../conductor-t10 (or T10_DIR),
z7's folder ../referee_z7 (or Z7_DIR); each falls back to MSS_REPO/experiments/<name>."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.environ.get('MSS_REPO', os.path.join(HERE, '..', '..')))
def _dir(env, name):
    if env in os.environ: return os.path.abspath(os.environ[env])
    d = os.path.join(HERE, '..', name)
    return os.path.abspath(d if os.path.isdir(d) else os.path.join(REPO, 'experiments', name))
T10 = _dir('T10_DIR', 'conductor-t10')
Z7 = _dir('Z7_DIR', 'referee_z7')
G1 = os.path.join(REPO, 'experiments', 'roadmap-g', 'g1')
G2 = os.path.join(REPO, 'experiments', 'roadmap-g', 'g2')
H1 = os.path.join(REPO, 'experiments', 'roadmap-h', 'h1')
for p in (HERE, Z7, os.path.join(REPO, 'scripts', 'cover_descent'), G1):
    if p not in sys.path: sys.path.insert(0, p)

class Out:
    """print to stdout and keep the lines; finish() compares with NAME.txt (on stderr) or writes it with --write"""
    def __init__(s, name):
        s.name = name; s.lines = []
    def __call__(s, line=''):
        s.lines.append(line); print(line, flush=True)
    def finish(s):
        path = os.path.join(HERE, s.name)
        text = '\n'.join(s.lines) + '\n'
        if '--write' in sys.argv:
            open(path, 'w').write(text); sys.stderr.write('# wrote %s\n' % s.name); return
        if os.path.exists(path) and open(path).read() == text:
            sys.stderr.write('# matches %s\n' % s.name)
        else:
            sys.stderr.write('# DIFFERS from %s\n' % s.name); sys.exit(1)
