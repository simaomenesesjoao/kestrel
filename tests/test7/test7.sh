#!/bin/bash
# Test graphene

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
energies=100000

print_flags="--NEnergies $energies --print_to_file" 
log_flags="--log_status"

cd ../..
MODEL=3 VERBOSE=1 STRIDE=9999999 make -e -j6 
rm src/*.o
mv kestrel tests/test7
cd tests/test7

# Geometry
input="--thread_division $nx $ny --geometry $Lx $Ly --nrandom $NR --ndisorder $ND --mult $B --seed $S --anderson $W --vacancies $C --nmoments $N"
options="$input $print_flags $log_flags $restart_flags"
./kestrel $options
