#!/bin/sh

echo 'Running with HOME = '${HOME}

CODEBASE=$HOME/work/gpProcess

RUN=$1
echo 'Found RUN = '${RUN}
cd $HOME/work/GPRuns/Run-${RUN}

mkdir PP
cd PP

cp -p ../pair.dat.gz .
cp -p ../lumi.ee.out.gz lumi-${RUN}.ee.out.gz
gunzip lumi-${RUN}.ee.out.gz
${CODEBASE}/delLL.sh lumi-${RUN}.ee.out

gunzip pairs.dat.gz
mv pairs.dat pairs-${RUN}.dat

sort -k 1n -k2gr -o pairs-${RUN}_sorted.dat pairs-${RUN}.dat

gunzip pairs0.dat.gz
mv pairs0.dat pairs0-${RUN}.dat

awk 'NR % 2 == 1 { printf "%s ", $0; next } { print }' pairs0-${RUN}.dat >cpairs0-${RUN}.dat
awk 'NR % 2 == 1 { printf "%s ", $0; next } { print }' pairs-${RUN}_sorted.dat >cpairs-${RUN}_sorted.dat
paste -d' ' cpairs0-${RUN}.dat cpairs-${RUN}_sorted.dat >qcombpairs-${RUN}.dat

cat ${CODEBASE}/Header-32field.txt qcombpairs-${RUN}.dat >qcombpairs-${RUN}.csv

module load root

root -l -b -q 'dfplot.C('\"${RUN}\"')';

exit
