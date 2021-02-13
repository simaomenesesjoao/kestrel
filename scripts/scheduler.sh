#!/bin/bash


dir=/home/simao/Sync/article_magnetic/program/kestrel/scripts/lockfile
source_dir=/home/simao/Sync/article_magnetic/program/kestrel
workdir=/home/simao/Sync/article_magnetic/program/kestrel/scripts
exec 100>$dir

while true; do
    flock 100
    # Read from the task list to see what is the next task
    cd $workdir
    next=$(python find_next.py)

    # If all tasks are done, then exit
    if [[ "$next" == "done" ]]; then
        echo DONE
        break
    fi

    # Fetch all the parameters needed for the next simulation, if
    # there are still tasks to do 
    Lx=$(echo $next | awk '{print $1}' | awk -F= '{print $2}')
    Ly=$(echo $next | awk '{print $2}' | awk -F= '{print $2}')
    S=$(echo $next | awk '{print $3}' | awk -F= '{print $2}')
    W=$(echo $next | awk '{print $4}' | awk -F= '{print $2}')
    nx=$(echo $next | awk '{print $5}' | awk -F= '{print $2}')
    ny=$(echo $next | awk '{print $6}' | awk -F= '{print $2}')
    B=$(echo $next | awk '{print $7}' | awk -F= '{print $2}')
    NR=$(echo $next | awk '{print $9}' | awk -F= '{print $2}')
    ND=$(echo $next | awk '{print $10}' | awk -F= '{print $2}')
    N=$(echo $next | awk '{print $11}' | awk -F= '{print $2}')
    dN=$(echo $next | awk '{print $12}' | awk -F= '{print $2}')
    Na=$(echo $next | awk '{print $13}' | awk -F= '{print $2}')
    stt=$(echo $next | awk '{print $14}' | awk -F= '{print $2}')

    ener=$(echo $next | awk '{print $8}' | awk -F= '{print $2}')
    energies=$(echo $ener | awk -F[ '{print $2}' | awk -F] '{print $1}')

    # Put the 'running' flag in that task so no other program
    # tries to complete it
    python update_tasks.py "B=$B" "status=R"
    flock -u 100

    # Do the task
    cur_dir=$workdir/simulations/$B
    mkdir -p $cur_dir
    cd $cur_dir

    cp -r $source_dir/src .
    cp $source_dir/Makefile .



    VERBOSE=1 STRIDE=9999 make -j4
    input="--thread_division $nx $ny --geometry $Lx $Ly --nrandom $NR --ndisorder $ND --mult $B --seed $S --anderson $W --nmoments $N"
    print_flags="--energies $energies --print_to_file"
    ./kestrel $input $print_flags > log.txt

    flock 100
    cd $workdir
    python update_tasks.py "B=$B" "Na=$N"
    python update_tasks.py "B=$B" "status=C"
    flock -u 100

done

