#!/bin/sh
#
# To compile ReadandDerive.cpp  
# ./cl.sh 
# Then the executable can be executed using ./filename

module load root/6.32.2
module list

target=ReadandDerive
echo 'Compiling (with ROOT libraries): '${target}.cpp

g++ -g -o ${target} ${target}.cpp `root-config --cflags --glibs`

module unload root/6.32.2

exit
