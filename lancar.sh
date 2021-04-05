#!/bin/bash
#source /users/js2892/virtualenv_python/bin/activate

Ly=4000
Lx=4001
nx=2
ny=2
NR=1
ND=1
B=10
N=2048
S=3
W=0.50
C=0.00
energies=16384


# This script requires the variables Lx,Ly,N,W,C,D,COMP,SE
if [ -z "$Lx" ]; then echo "system size Lx is not defined"; exit; fi
if [ -z "$Ly" ]; then echo "system size Ly is not defined"; exit; fi
if [ -z "$N" ]; then echo "number of cheb moments N is not defined"; exit; fi
if [ -z "$C" ]; then echo "vacancy concentration C is not defined"; exit; fi
if [ -z "$W" ]; then echo "Anderson strength W is not defined"; exit; fi
if [ -z "$B" ]; then echo "Magnetic field B is not defined"; exit; fi
if [ -z "$D" ]; then echo "divisions D is not defined"; exit; fi
if [ -z "$COMP" ]; then echo "compilation flag COMP is not defined"; exit; fi
if [ -z "$SE" ]; then echo "seed SE is not defined"; exit; fi

THREADS=$((D*D))

# Check if SLURM_JOBID is defined
if [ -z "$SLURM_JOBID" ]; then 
    a=$(cat ~/inc.txt); b=$(echo $((a+1))); echo $b > ~/inc.txt
    SLURM_JOBID=$b
fi

script_name=config.py
launchdir=${HOME}/projects_sync/edelstein_aires/
kitedir=${HOME}/projects_sync/codes/kite_code/kite_ss/
workdir=$launchdir/simulations/test3/$SLURM_JOBID
mkdir -p $workdir

if [ $HOME == "/homes/cfp-smjoao" ]; then 
  module unload hdf/5.1.10.1
  module load hdf/5.1.10.1
  module load gcc/9.3.1
fi

# With compilation
if [ $COMP -eq 1 ]; then
    cd $workdir

    # Build required directory structure
    mkdir kite
    ln -s $kitedir $workdir/kite/link
    mkdir kite/build
    mkdir kite/tools
    mkdir kite/tools/build

    cp -r $kitedir/Src                  $workdir/kite/
    cp -r $kitedir/CMakeLists.txt       $workdir/kite/
    cp -r $kitedir/tools/src            $workdir/kite/tools/
    cp -r $kitedir/tools/CMakeLists.txt $workdir/kite/tools/
    chmod -R 755 $workdir/kite

    # Compilation step
    cd $workdir/kite/build
    cmake .. && make -j$THREADS
    cd $workdir/kite/tools/build
    cmake .. && make -j$THREADS

    # Move the compiled codes to the directory needed
    mv $workdir/kite/build/KITEx $workdir/kite
    mv $workdir/kite/tools/build/KITE-tools $workdir/kite

    # Remove excess code
    rm -r $workdir/kite/Src
    rm -r $workdir/kite/build
    rm -r $workdir/kite/CMakeLists.txt
    rm -r $workdir/kite/tools

# Without compilation
else
    mkdir $workdir/kite
    cp $kitedir/build/KITEx $workdir/kite/
    cp $kitedir/tools/build/KITE-tools $workdir/kite/
    chmod -R 755 $workdir/kite


fi


cd $workdir
cp $launchdir/*.py $workdir
mkdir params
echo $INFO > $workdir/params/info.txt

echo "building conf script"
echo E$E S$S SE$SE L$L N$N C$C R$R D$D COMP$COMP > $workdir/params/params.txt
python $workdir/$script_name $L $N $R $C $D $E $S $flag
SEED=$SE $workdir/kite/KITEx config.h5

halfM=$((N/2))
NE=8192
#max_E=4.0
$workdir/kite/KITE-tools config.h5 --DOS -E $NE -M $halfM -N dos_N${halfM}.dat
$workdir/kite/KITE-tools config.h5 --DOS -E $NE -N dos_N${N}.dat

rm kite/KITEx
rm kite/KITE-tools
