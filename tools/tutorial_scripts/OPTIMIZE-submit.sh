#!/bin/bash
#
EXECUTABLE=$EXCITINGROOT/bin/exciting_mpismp

label=`ls -d *_??`
for dirn in $label ; do
    cd $dirn
    cp -f $dirn.xml input.xml
    echo
    echo 'SCF calculation of "'$dirn'" starts --------------------------------'
    time mpirun -n 4 $EXECUTABLE | tee output.screen
    cd ../
done
echo 
