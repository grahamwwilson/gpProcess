#!/bin/sh
#
# Process Guinea-PIG files for specified run on HPC.
#

echo 'Running with HOME = '${HOME}

CODEBASE=$HOME/work/gpProcess

RUN=$1
echo 'Found RUN = '${RUN}
cd $HOME/work/GPRuns/Run-${RUN}

mkdir PP
cd PP

# Prepares "pairs" files for analayzing the scattered Bhabhas after EM deflection
cp -p ../pair*.dat.gz .
gunzip pairs.dat.gz
mv pairs.dat pairs-${RUN}.dat

sort -k 1n -k2gr -o pairs-${RUN}_sorted.dat pairs-${RUN}.dat

gunzip pairs0.dat.gz
mv pairs0.dat pairs0-${RUN}.dat

awk 'NR % 2 == 1 { printf "%s ", $0; next } { print }' pairs0-${RUN}.dat >cpairs0-${RUN}.dat
awk 'NR % 2 == 1 { printf "%s ", $0; next } { print }' pairs-${RUN}_sorted.dat >cpairs-${RUN}_sorted.dat
paste -d' ' cpairs0-${RUN}.dat cpairs-${RUN}_sorted.dat >qcombpairs-${RUN}.dat
cat ${CODEBASE}/Header-32field.txt qcombpairs-${RUN}.dat >qcombpairs-${RUN}.csv
rm *.dat

module load root

cp ${CODEBASE}/dfplot.C .
root -l -b -q 'dfplot.C('\"${RUN}\"')';
rm qcombpairs-${RUN}.csv

cp -p ../lumi.ee.out.gz lumi-${RUN}.ee.out.gz
gunzip lumi-${RUN}.ee.out.gz
${CODEBASE}/delLL.sh lumi-${RUN}.ee.out
cat ${CODEBASE}/LumiHeader-17field.txt lumi-${RUN}.ee.out >lumiee-${RUN}.csv
rm lumi-${RUN}.ee.out
cp ${CODEBASE}/dlplot.C .
root -l -b -q 'dlplot.C('\"${RUN}\"')';

rm lumiee-${RUN}.csv

#Clean up any remaining intermediate files ?
rm *.C

exit
