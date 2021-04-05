#!/bin/bash
# Test different geometries of the same model

set -e
Ly=100
Lx=101
nx=2
ny=2
NR=1
ND=1
B=10
N=65536
S=3
W=0.00
C=0.00
energies=200000
filename="dos_N${N}_W0.000000_C0.000000_B${B}.dat"

print_flags="--NEnergies $energies --print_to_file" 
log_flags="--log_status"

cd ../..
MODEL=1 VERBOSE=1 make -e -j6 
rm src/*.o
mv kestrel tests/test8
cd tests/test8

input="--thread_division $nx $ny --geometry $Lx $Ly --nrandom $NR --ndisorder $ND --mult $B --seed $S --anderson $W --vacancies $C --nmoments $N"
options="$input $print_flags $log_flags $restart_flags"
./kestrel $options
