#!/bin/bash
# Test different models of the same Hamiltonian

set -e
Ly=10
nx=1
ny=1
NR=1
ND=1
B=1
N=4096
S=3
W=0.00
C=0.00
energies=16384
filename="dos_N4096_W0.000000_C0.000000_B${B}.dat"

print_flags="--NEnergies $energies --print_to_file" 
log_flags="--log_status"

# Model 1
Lx=10
input="--thread_division $nx $ny --geometry $Lx $Ly --nrandom $NR --ndisorder $ND --mult $B --seed $S --anderson $W --vacancies $C --nmoments $N"
options="$input $print_flags $log_flags $restart_flags"

cd ..
MODEL=1 VERBOSE=1 STRIDE=9999999 make -e -j6 
rm src/*.o
mv kestrel tests
cd tests
./kestrel $options
mv $filename dos_sq1.dat

# Model 2
Lx=5
input="--thread_division $nx $ny --geometry $Lx $Ly --nrandom $NR --ndisorder $ND --mult $B --seed $S --anderson $W --vacancies $C --nmoments $N"
options="$input $print_flags $log_flags $restart_flags"

cd ..
MODEL=2 VERBOSE=1 STRIDE=9999999 make -e -j6 
rm src/*.o
mv kestrel tests
cd tests
./kestrel $options
mv $filename dos_sq2.dat
