#!/bin/sh
#
#
# Process Guinea-PIG files for specified run on HPC.
# This has four main parts as documented below.
#

module load root
module list
root --version

echo 'Running with HOME = '${HOME}

CODEBASE=$HOME/work/gpProcess

RUN=$1
echo 'Found RUN = '${RUN}
cd $HOME/work/GPRuns/Run-${RUN}

if [ -d PP ]; then
   echo "PP directory already exists .."
else   
   mkdir PP
fi
cd PP

#----------------------------------------------------------------------
# STEP 0: Parse the results file for just the binned luminosity numbers
#----------------------------------------------------------------------
${CODEBASE}/lumiprofile.sh ${RUN}

# -----------------------------------------------------------------------------
# STEP 1
# Prepare "pairs" files for analyzing the scattered Bhabhas after EM deflection
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

cp ${CODEBASE}/dfplot.C .
root -l -b -q 'dfplot.C('\"${RUN}\"')';
rm qcombpairs-${RUN}.csv

# -----------------------------------------------------------------------------
# STEP 2
# Prepare lumi root file

cp -p ../lumi.ee.out.gz lumi-${RUN}.ee.out.gz
gunzip lumi-${RUN}.ee.out.gz
${CODEBASE}/delLL.sh lumi-${RUN}.ee.out
cat ${CODEBASE}/LumiHeader-17field.txt lumi-${RUN}.ee.out >lumiee-${RUN}.csv
rm lumi-${RUN}.ee.out
cp ${CODEBASE}/dlplot.C .
root -l -b -q 'dlplot.C('\"${RUN}\"')';

rm lumiee-${RUN}.csv

# Clean up any remaining intermediate files ?
rm *.C

# -----------------------------------------------------------------------------
# STEP 3
# Run example ReadandDerive code
${CODEBASE}/ReadandDerive -h
${CODEBASE}/ReadandDerive -n 10 >RandD-${RUN}.out
mv EMD-Analysis.root EMD-${RUN}-Analysis.root

# -----------------------------------------------------------------------------
# STEP 4
# Run current eepairs file analysis
cp ${CODEBASE}/Ana.* .
cp ${CODEBASE}/macro.C .
# Run root in "batch mode" using macro.C to process the eepairs root file 
root -l -b -q ${RUN}.root macro.C
mv Ana.root Ana-${RUN}.root
# Clean up
rm Ana_*
# Keep Ana.C and Ana.h so that they can be checked
#rm *.C
#rm *.h

module unload root

exit
