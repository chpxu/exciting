#!/bin/bash
#
function run_exciting() {

mpirun -n 8 "$1/bin/exciting_smp"
}

# EXECUTABLE=

label=`ls -d *_??`
for dirn in $label ; do
    cd $dirn
    cp -f $dirn.xml input.xml
    echo
    echo 'SCF calculation of "'$dirn'" starts --------------------------------'
    # time $EXECUTABLE | tee output.screen
    time mpirun $EXCITINGROOT/bin/exciting_purempi | tee output.screen
    cd ../
done
echo 
