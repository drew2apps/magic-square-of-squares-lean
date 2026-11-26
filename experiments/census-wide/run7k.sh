set -e
cd /mnt/project-files/mss-experiments/census-wide
python3 ../census2000/gen_pairs.py 5000 7000 pairs_5000_7000.json > gen7k.log
python3 -c "
import json,math
P=json.load(open('pairs_5000_7000.json'))['open']
def sq(x): r=math.isqrt(x); return r*r==x
D=[p for p in P if sq(p[0]**2+p[1]**2) or ((p[0]**2+p[1]**2)%2==0 and sq((p[0]**2+p[1]**2)//2))]
json.dump(D,open('degen_5000_7000.json','w')); print(len(P),len(D))
" >> gen7k.log
for k in 0 1 2 3; do nohup python3 wide_pari.py degen_5000_7000.json $k 4 > s7_$k.jsonl 2>&1 & done
wait
