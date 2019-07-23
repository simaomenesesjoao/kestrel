#include <Eigen/Core>
#include <fstream>
#include <chrono>
#include <iostream>
#include <complex>

#include "general.hpp"
#include "kpm_vector.hpp"
#include "hamiltonian.hpp"
#include "cheb.hpp"

chebinator::chebinator(){
    ;
};

void chebinator::set_hamiltonian(hamiltonian *new_H){
    H = new_H;
    Lx = H->Lx;
    Ly = H->Ly;
    Norb = H->Norb;

    KPM0.set_geometry(Lx, Ly, Norb);
    KPM1.set_geometry(Lx, Ly, Norb);
    KPM_initial.set_geometry(Lx, Ly, Norb);

}

void chebinator::cheb_iteration(unsigned n, unsigned ndis, unsigned nrand){
    debug_message("Entered cheb_iteration\n");


    // set some useless variables that are only going to be needed
    // for displaying the output with system signals
    current_ndis = ndis;
    current_moments = n;
    current_nrand = nrand;



    start = std::chrono::high_resolution_clock::now();



    unsigned moments = n;
    max_iter = moments*nrand*ndis;
    mu = Eigen::Array<TR, -1, -1>::Zero(moments, 1);
    mu(0) = 1.0;

    // Implement a more efficient way to compute averages
    for(unsigned d = 0; d < ndis; d++){
        debug_message("disorder realization " + std::to_string(d) + "\n")

        // new realization of disorder
        H->set_anderson();
        for(unsigned r = 0; r < nrand; r++){
            debug_message("random realization " + std::to_string(r) + "\n")

            KPM0.random_uniform();
            KPM0.fill_ghosts();
            KPM1.zero();
            KPM_initial = KPM0;
            KPM_initial.empty_ghosts();


            for(unsigned i = 0; i < moments-1; i++){
                debug_message("cheb iteration " + std::to_string(i) + "\n")

                H->cheb(KPM1, KPM0, i);
                swap_pointers(&KPM1, &KPM0);
                KPM0.fill_ghosts();
                mu(i+1) += (KPM_initial*KPM0)/nrand/ndis;

                current_iter = i + r*moments + d*moments*nrand;
                current = std::chrono::high_resolution_clock::now();
            }
        }
    }
    debug_message("Left cheb_iteration\n");
}

void chebinator::set_seed(int seed2){
    seed = seed2;
}

void chebinator::set_mult(int mult2){
    mult = mult2;
}

void chebinator::calc_finish(){
    debug_message("Entered chebinator::calc_finis()\n");
    /* Prints information about the current state of the program to the user's
     * home directory. This information consists of all the system's parameters
     * as well as the progress to being finished and the estimated time remaining
     */

    // Calculate the progress and the time remaining
    std::chrono::duration<double> time_span;
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(current - start);
    double way_done = current_iter/double(max_iter);
    double avg_iter_time = time_span.count()/current_iter;

    // Append to a file in the user's home directory
    const char *homedir;
    homedir = getenv("HOME");
    std::ofstream file;
    file.open(std::string(homedir) + "/kestrel_info.txt", std::ofstream::app);

    // Output the information about the program
    file << "Lx,Ly=" << Lx << "," << Ly;
    file << " N=" << current_moments;
    file << " ND=" << current_ndis;
    file << " NR=" << current_nrand;
    file << " W=" << H->W*SCALE;
    file << " S=" << seed;
    file << " B=" << mult;
    file << " Progress: " << way_done*100 << "%.";
    file << " Time left: " << avg_iter_time*max_iter*(1.0-way_done) << "s\n";
    file.close();

    debug_message("Left chebinator::calc_finis()\n");
}

double chebinator::get_estimate(double time, unsigned moments, unsigned num_disorder, unsigned num_random){
    return time*moments*num_disorder*num_random;
}

void chebinator::set_geometry(unsigned lx, unsigned ly, unsigned norb){
    Lx = lx;
    Ly = ly;
    Norb = norb;
    //LxG = Lx + 2;
    //LyG = Ly + 2;
    //Npoints = Lx*Ly*Norb;

    // Now that we know the number of orbitals, allocate
    // the memory needed for all the KPM vectors
    KPM0.set_geometry(Lx, Ly, norb);
    KPM1.set_geometry(Lx, Ly, norb);
    KPM_initial.set_geometry(Lx, Ly, norb);

}
