#!/bin/sh
#
# Read and histogram the binned luminosity numbers
#

RUN=$1
echo ${RUN}
pwd

CODEBASE=$HOME/work/gpProcess

#---------------------------------------------------------------
# STEP 0. Parse the file for just the binned luminosity numbers
#---------------------------------------------------------------

cp ../GPResults.out.gz .
gunzip GPResults.out.gz

# Try to identify the line numbers associated with 
#step   lumi total   lumi peak
# and
#------------- other information  ------------

grep -n "step   lumi total   lumi peak" GPResults.out
grep -n "step   lumi total   lumi peak" GPResults.out | cut -d: -f1

FLINE="$(grep -n "step   lumi total   lumi peak" GPResults.out | cut -d: -f1)"
echo 'FLINE '$FLINE

tail -n +${FLINE} GPResults.out > temp && mv temp GPResults.out
tail -n +2 GPResults.out > temp && mv temp GPResults.out

grep -n "other information" GPResults.out
grep -n "other information" GPResults.out | cut -d: -f1
LLINE="$(grep -n "other information" GPResults.out | cut -d: -f1)"
echo 'LLINE '$LLINE

head -n ${LLINE} GPResults.out > temp && mv temp GPResults.out

# Now we want to delete the last 3 lines

head -n -3 GPResults.out > temp && mv temp GPLumi-Summary.txt
rm GPResults.out

${CODEBASE}/CheckLumiProfile -h
${CODEBASE}/CheckLumiProfile >CheckLumi.out

mv GPLumi.root GPLumi-${RUN}.root

exit
