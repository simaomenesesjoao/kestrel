#!/bin/bash

set -e
Lx=1000
Ly=1001
NR=1
ND=1
B=100
N1=200
N2=300
N3=500
S=3
W=0.1
energies="0.05 0.1 0.2 0.3"
input1="--geometry $Lx $Ly --nrandom $NR --ndisorder $ND --mult $B --seed $S --anderson $W --nmoments $N1"
input2="--geometry $Lx $Ly --nrandom $NR --ndisorder $ND --mult $B --seed $S --anderson $W --nmoments $N2"
input3="--geometry $Lx $Ly --nrandom $NR --ndisorder $ND --mult $B --seed $S --anderson $W --nmoments $N3"
#input="--readfrom restart.h5 --moremoments 50"
print_flags="--energies $energies --print_to_cout"
restart_flags="--save_restart_to restart.h5"
log_flags="--log_status"

options1="$input1 $print_flags $log_flags $restart_flags"
options2="$input2 $print_flags $log_flags $restart_flags"
options3="$input3 $print_flags $log_flags $restart_flags"

VERBOSE=1 STRIDE=9999999 make -e -j6 
./kestrel $options1 &
./kestrel $options2 &
./kestrel $options3 &
