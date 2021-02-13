#!/bin/bash

Lx=50
Ly=51
NR=1
ND=1
B=10
N1=10
N2=20
S=3
W=0.1
energies="0.05 0.1 0.2 0.3"
input1="--geometry $Lx $Ly --nrandom $NR --ndisorder $ND --mult $B --seed $S --anderson $W --nmoments $N2"
input2="--geometry $Lx $Ly --nrandom $NR --ndisorder $ND --mult $B --seed $S --anderson $W --nmoments $N1"
input3="--readfrom restart2.h5 --moremoments $N1"

print_flags="--energies $energies --print_to_cout"

restart_flags1="--save_restart_to restart1.h5"
restart_flags2="--save_restart_to restart2.h5"
restart_flags3="--save_restart_to restart2.h5"

log_flags="--log_status"

options1="$input1 $print_flags $log_flags $restart_flags1"
options2="$input2 $print_flags $log_flags $restart_flags2"
options3="$input3 $print_flags $log_flags $restart_flags3"

VERBOSE=1 STRIDE=9999999 make -e -j6
./kestrel $options1 > /dev/null
./kestrel $options2 > /dev/null
./kestrel $options3 > /dev/null

h5dump -d /MU restart1.h5 | tail -n +6 | head -n +4
h5dump -d /MU restart2.h5 | tail -n +6 | head -n +4
