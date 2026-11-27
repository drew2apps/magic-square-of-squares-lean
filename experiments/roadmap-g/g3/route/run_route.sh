#!/bin/bash
# G3: the rank-1 route (../../../census-wide/route/route_pari.py, unchanged) on ratios left after Demjanenko-Manin whose
# '+13' class d gives E1 rank bounds [1,1] (scope_pari.py).  usage: bash run_route.sh m,n d   (writes log_<m>.txt)
cd "$(dirname "$0")"; R=../../../census-wide/route
m=${1%%,*}; ( cd $R && nice python3 route_pari.py $1 $2 ) > log_$m.txt 2>&1
