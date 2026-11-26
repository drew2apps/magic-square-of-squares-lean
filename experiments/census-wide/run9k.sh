set -e
cd /mnt/project-files/mss-experiments/census-wide
python3 ../census2000/gen_pairs.py 7000 9000 pairs_7000_9000.json > gen9k.log
python3 -c "
import json,math
P=json.load(open('pairs_7000_9000.json'))['open']
def sq(x): r=math.isqrt(x); return r*r==x
D=[p for p in P if sq(p[0]**2+p[1]**2) or ((p[0]**2+p[1]**2)%2==0 and sq((p[0]**2+p[1]**2)//2))]
json.dump(D,open('degen_7000_9000.json','w')); print(len(P),len(D))
" >> gen9k.log
for k in 0 1 2; do nohup python3 wide_pari.py degen_7000_9000.json $k 3 > s9_$k.jsonl 2>&1 & done
wait
