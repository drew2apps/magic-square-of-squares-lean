#!/bin/bash
# usage: bash run_h6_rich.sh m n   (Richelot descent on the genus-2 pieces; times.txt)
cd "$(dirname "$0")"
s=$(date +%s); python3 h6_richelot.py $1 $2 --write 2> err_rich_$1_$2.txt; echo "h6_richelot $1 $2 exit=$? $(( $(date +%s) - s )) s" >> times.txt
