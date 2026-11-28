#!/bin/bash
H=$1
cat tasks_iso.txt | xargs -P 2 -L 1 sh -c './search2 '$H' "$0" "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8"' > isohits_H$H.txt
