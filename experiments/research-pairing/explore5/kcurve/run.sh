#!/bin/bash
# run.sh H : all tasks, two at a time
H=$1
cat tasks.txt | xargs -P 2 -L 1 sh -c './search '$H' "$0" "$1" "$2" "$3" "$4" "$5" "$6"' > hits_H$H.txt
