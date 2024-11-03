#!/bin/sh

FN=$1
echo 'Examining file '${FN}

LLINE="$(tail -n 1 ${FN} | wc -c)"
echo $LLINE

if [ $LLINE -eq 1 ]; then
   echo 'Found empty last line in file'
   echo 'Will proceed with deleting it'
   echo 'Deleting last line of file named '${FN}
   echo 'Current file length: '
   wc -l ${FN}
   echo 'Last 3 lines of current file (with $ signifying EOL)'
   tail -3 ${FN} | cat -vet

   head -n -1 ${FN} >temp.txt
   mv temp.txt ${FN}

   echo 'New file length'
   wc -l ${FN}
   echo 'Last 2 lines of current file (with $ signifying EOL)'
   tail -2 ${FN} | cat -vet
else
   echo 'Last line is not empty'
   echo 'Leaving file as is'
fi

exit
