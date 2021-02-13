#!/bin/bash

set -e
Lx=200
Ly=200
nx=2
ny=2
NR=1
ND=1
B=0
N=500
S=3
W=0.40
#energies="0.01 0.02 0.03 0.05 0.1 0.2 0.3"
energies=1000
input="--thread_division $nx $ny --geometry $Lx $Ly --nrandom $NR --ndisorder $ND --mult $B --seed $S --anderson $W --nmoments $N"
#input="--readfrom restart.h5 --moremoments 50"
#print_flags="--energies $energies --print_to_cout"
print_flags="--NEnergies $energies --print_to_file" 
#restart_flags="--save_restart_to restart.h5"
log_flags="--log_status"

options="$input $print_flags $log_flags $restart_flags"

VERBOSE=4 STRIDE=9999999 make -e -j6 
./kestrel $options
