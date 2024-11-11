#!/bin/bash

N1=$1
N2=$2
echo "N1 = "$N1
echo "N2 = "$N2

for ((i = ${N1}; i <= ${N2}; i++)); do
    RUN=Z-${i}
    echo "Run ${RUN}"
    ./run.sh ${RUN}
done

exit
