set -e
cd "$(dirname "$0")"
python3 ../census2000/gen_pairs.py 9000 11000 pairs_9000_11000.json > gen11k.log
python3 -c "
import json,math
P=json.load(open('pairs_9000_11000.json'))['open']
def sq(x): r=math.isqrt(x); return r*r==x
D=[p for p in P if sq(p[0]**2+p[1]**2) or ((p[0]**2+p[1]**2)%2==0 and sq((p[0]**2+p[1]**2)//2))]
json.dump(D,open('degen_9000_11000.json','w')); print(len(P),len(D))
" >> gen11k.log
for k in 0 1; do nohup python3 wide_pari.py degen_9000_11000.json $k 2 > s11_$k.jsonl 2>&1 & done
wait
