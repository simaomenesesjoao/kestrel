#!/bin/bash
# Test different geometries of the same model

set -e
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

cd ../..
MODEL=1 VERBOSE=1 STRIDE=9999999 make -e -j6 
rm src/*.o
mv kestrel tests/test5
cd tests/test5

# 210=21x10=3x7x5x2
# Geometry
Ly=14
Lx=15
input="--thread_division $nx $ny --geometry $Lx $Ly --nrandom $NR --ndisorder $ND --mult $B --seed $S --anderson $W --vacancies $C --nmoments $N"
options="$input $print_flags $log_flags $restart_flags"
./kestrel $options
mv $filename dos_sq1_geo1.dat


# Model 2
Ly=10
Lx=21
input="--thread_division $nx $ny --geometry $Lx $Ly --nrandom $NR --ndisorder $ND --mult $B --seed $S --anderson $W --vacancies $C --nmoments $N"
options="$input $print_flags $log_flags $restart_flags"

./kestrel $options
mv $filename dos_sq1_geo2.dat
