"""PARI port of scripts/second_descent/wide.py's one(): first descent (must leave the 2 expected tuples), second descent,
then the partner if still open.  EXPERIMENT, written by an AI model (Claude).  usage: wide_pari.py pairs.json K J  (does pairs[K::J])"""
import sys, json, time
sys.path.insert(0, '/mnt/project-files/mss-experiments/census-wide')
import sd_pari as sd
def one(m, n):
    c, differ = sd.first(m, n); out = dict(pair=[m, n], first=c, differ=differ)
    if c != 2 or not differ: out['status'] = 'first descent: not the expected 2 tuples'; return out
    r = sd.close(m, n); out['status'] = r['status']; out['by'] = r.get('by')
    if r['status'] == 'closed':
        cc = r['curves'][-1]; out['rank_bounds'] = cc['rank_bounds']; out['rational_t'] = cc['rational_t']
    return out
if __name__ == '__main__':
    P = json.load(open(sys.argv[1])); k, j = int(sys.argv[2]), int(sys.argv[3])
    for m, n in P[k::j]:
        t0 = time.time()
        try:
            o = one(m, n)
            if o['status'] == 'open':
                a, b = m + n, m - n
                if a % 2 == 0: a, b = a//2, b//2
                po = one(a, b); o['partner'] = po
                if po['status'] == 'closed': o['status'] = 'closed by partner'
        except Exception as ex:
            o = dict(pair=[m, n], status='error %s %s' % (type(ex).__name__, str(ex)[:80]))
        o['sec'] = round(time.time() - t0, 1)
        print(json.dumps(o, default=str), flush=True)
