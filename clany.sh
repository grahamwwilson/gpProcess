#!/bin/sh
#
# To compile filename.cpp  
# ./clany.sh filename
# Then the executable can be executed using ./filename
# Here we don't use ROOT so don't make it depend on it.

module load root

target=$1
echo 'Compiling (with ROOT libraries): '${target}.cpp

g++ -g -o ${target} ${target}.cpp `root-config --cflags --glibs`

exit
