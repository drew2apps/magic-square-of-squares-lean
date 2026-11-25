#!/usr/bin/env python3
"""MANIFEST.sha256: the SHA-256 of every file of the repository, so that a copy can be checked without network access.

    python3 scripts/verify_manifest.py            verify: every listed file has the listed hash, and no file is unlisted
    python3 scripts/verify_manifest.py --write    rewrite MANIFEST.sha256 from the files present

The manifest lists every file under the root except itself, the directories .git and .lake, and Python caches.  With
git available, the list of files is checked against `git ls-files` as well.  CI runs the verification on every push, so
a result file cannot change without the manifest changing in the same commit.  The manifest says that the files are
the ones committed; it does not say that they are right."""
import hashlib, os, subprocess, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NAME = 'MANIFEST.sha256'
SKIP_DIRS = {'.git', '.lake', '__pycache__'}

def files():
    out = []
    for d, sub, fs in os.walk(ROOT):
        sub[:] = sorted(x for x in sub if x not in SKIP_DIRS)
        for f in sorted(fs):
            rel = os.path.relpath(os.path.join(d, f), ROOT).replace(os.sep, '/')
            if rel != NAME and rel != '.git' and not rel.endswith('.pyc'): out.append(rel)   # .git is a file in a worktree
    return sorted(out)

def sha(rel):
    h = hashlib.sha256()
    with open(os.path.join(ROOT, rel), 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''): h.update(block)
    return h.hexdigest()

if __name__ == '__main__':
    present = files()
    if len(sys.argv) > 1 and sys.argv[1] == '--write':
        with open(os.path.join(ROOT, NAME), 'w', newline='\n') as fh:
            for rel in present: fh.write('%s  %s\n' % (sha(rel), rel))
        print('wrote %s: %d files' % (NAME, len(present))); sys.exit(0)
    listed = {}
    for line in open(os.path.join(ROOT, NAME), encoding='utf-8'):
        h, rel = line.rstrip('\n').split('  ', 1); listed[rel] = h
    bad = [rel for rel in listed if rel in present and sha(rel) != listed[rel]]
    missing = sorted(set(listed) - set(present)); unlisted = sorted(set(present) - set(listed))
    try:
        tracked = subprocess.run(['git', '-C', ROOT, 'ls-files'], capture_output=True, text=True, check=True).stdout.split('\n')
        tracked = sorted(t for t in tracked if t and t != NAME)
        untracked_listed = sorted(set(listed) - set(tracked)); git_note = '; git: %d tracked files' % len(tracked)
    except Exception:
        untracked_listed = []; git_note = '; git not available'
    for name, lst in (('hash differs', bad), ('listed but missing', missing), ('present but not listed', unlisted), ('listed but not tracked by git', untracked_listed)):
        for rel in lst: print('%s: %s' % (name, rel))
    ok = not (bad or missing or unlisted or untracked_listed)
    print('MANIFEST %s: %d files listed%s' % ('OK' if ok else 'FAILED', len(listed), git_note))
    sys.exit(0 if ok else 1)
