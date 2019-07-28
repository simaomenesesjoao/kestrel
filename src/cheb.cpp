#include <Eigen/Core>
#include <fstream>
#include <chrono>
#include <iostream>
#include <complex>

#include "general.hpp"
#include "kpm_vector.hpp"
#include "hamiltonian.hpp"
#include "cheb.hpp"

#include <H5Cpp.h>
#include <H5Group.h>
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"

chebinator::chebinator(){
    debug_message("Entered chebinator constructor\n");
    ;
    debug_message("Left chebinator constructor\n");
};

//chebinator::~chebinator(){
    //debug_message("Entered chebinator destructor\n");
    //;
    //debug_message("Left chebinator destructor\n");
//};

void chebinator::set_hamiltonian(hamiltonian *new_H){
    H = new_H;
    Lx = H->Lx;
    Ly = H->Ly;
    Norb = H->Norb;

    KPM0.set_geometry(Lx, Ly, Norb);
    KPM1.set_geometry(Lx, Ly, Norb);
    KPM_initial.set_geometry(Lx, Ly, Norb);

}


void chebinator::cheb_iteration_restart(unsigned n, KPM_vector &KPM0_restart, KPM_vector &KPM1_restart, Eigen::Array<TR, -1, -1> mu_restart){
//void chebinator::cheb_iteration_restart(unsigned n, Eigen::Array<TR, -1, -1> mu_restart){
    debug_message("Entered cheb_iteration\n");


    // set some useless variables that are only going to be needed
    // for displaying the output with system signals
    current_ndis = 1;
    current_nrand = 1;



    start = std::chrono::high_resolution_clock::now();



    unsigned moments = n;
    mu = Eigen::Array<TR, -1, -1>::Zero(moments, 1);
    mu(0) = 1.0;


    // These two steps are needed so that the disorder realization
    // and the initial random vectors are the same
    H->set_anderson();
    KPM0.random_uniform();
    KPM_initial = KPM0;
    KPM_initial.empty_ghosts();

    //load(name);
    unsigned previous_moments = mu_restart.size();
    max_iter = moments + previous_moments;

    KPM0 = KPM0_restart;
    KPM1 = KPM1_restart;
    mu = Eigen::Array<TR, -1, -1>::Zero(previous_moments + n, 1);
    mu.block(0,0,previous_moments,1) = mu_restart;

    //Eigen::Array<TR, -1, -1> mu_new = mu;
    //std::cout << "previous_moments: " << previous_moments << "\n";
    current_moments = max_iter;

    //std::cout << "mu:\n" << mu << "\n";
    //std::cout << "KPM0[0]" << KPM0.KPM[0] << "\n";
    // only after the seed has been used, load the vectors



    for(unsigned i = previous_moments; i < previous_moments + moments; i++){
        debug_message("cheb iteration " + std::to_string(i) + "\n");
        //std::cout << "iteration: " << i << "\n";
        //std::cout << "Before H\n";
        //std::cout << "KPM0.0:\n" << KPM0.KPM[0] << "\n";
        //std::cout << "KPM1.0:\n" << KPM1.KPM[0] << "\n";

        H->cheb(KPM1, KPM0, i);
        //std::cout << "after H\n";
        //std::cout << "KPM0.0:\n" << KPM0.KPM[0] << "\n";
        //std::cout << "KPM1.0:\n" << KPM1.KPM[0] << "\n";

        swap_pointers(&KPM1, &KPM0);
        KPM0.fill_ghosts();
        mu(i) = KPM_initial*KPM0;

        current_iter = i;
        current = std::chrono::high_resolution_clock::now();
    }
    debug_message("Left cheb_iteration\n");
    //std::cout << "mu:\n" << mu << "\n";
    //std::cout << "KPM0[0]" << KPM0.KPM[0] << "\n";
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
        debug_message("disorder realization " + std::to_string(d) + "\n");

        // new realization of disorder
        H->set_anderson();
        for(unsigned r = 0; r < nrand; r++){
            debug_message("random realization " + std::to_string(r) + "\n");

            // use the same seed to generate the same initial vector
            KPM0.random_uniform();
            KPM1.zero();
            KPM0.fill_ghosts();
            KPM_initial = KPM0;
            KPM_initial.empty_ghosts();
            //std::cout << "KPM_initial:\n" << KPM_initial.KPM[0] << "\n";

            for(unsigned i = 0; i < moments-1; i++){
                debug_message("cheb iteration " + std::to_string(i) + "\n");
                //std::cout << "iteration: " << i << "\n";

                //std::cout << "before H\n";
                //std::cout << "KPM0.0:\n" << KPM0.KPM[0] << "\n";
                //std::cout << "KPM1.0:\n" << KPM1.KPM[0] << "\n";
                H->cheb(KPM1, KPM0, i);
                //std::cout << "after H\n";
                //std::cout << "KPM0.0:\n" << KPM0.KPM[0] << "\n";
                //std::cout << "KPM1.0:\n" << KPM1.KPM[0] << "\n";

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

void chebinator::save(std::string name){
    debug_message("Entered chebinator::save\n");

    H5::H5File * file1 = new H5::H5File(name, H5F_ACC_TRUNC);

    // Write information about the Hamiltonian
    debug_message("Writing information about the Hamiltonian\n");
    write_hdf5(&(H->Lx), file1, "/Lx");
    write_hdf5(&(H->Ly), file1, "/Ly");
    write_hdf5(&(H->Norb), file1, "/Norb");
    write_hdf5(&(H->W), file1, "/anderson_w");
    write_hdf5(&mult, file1, "/mult");

    // Save information about the Chebyshev object
    debug_message("Writing information about the Chebyshev iteration\n");
    write_hdf5(&current_ndis, file1, "/num_disorder");
    write_hdf5(&current_nrand, file1, "/num_random");
    write_hdf5(&current_moments, file1, "/num_moments");
    write_hdf5(&seed, file1, "/seed");

    // Save current information about the chebyshev iterator
    for(unsigned n = 0; n < H->Norb; n++){
        write_hdf5(KPM0.KPM[n], file1, "/KPM0_"+std::to_string(n));
        write_hdf5(KPM1.KPM[n], file1, "/KPM1_"+std::to_string(n));
    }
    write_hdf5(mu, file1, "/MU");
    file1->close();
    delete file1;

    debug_message("Left chebinator::save\n");
}

void chebinator::load(std::string name){
    debug_message("Entered chebinator::load\n");

    H5::H5File * file = new H5::H5File(name, H5F_ACC_RDWR);

    Eigen::Array<TR, -1, -1> mu1;
    unsigned moms;

    get_hdf5(&moms, file, (char *) "/num_moments");
    mu1 = Eigen::Array<TR, -1, -1>::Zero(moms, 1);
    get_hdf5(mu1.data(), file, (char *) "/MU");

    mu = mu1;

    for(unsigned n = 0; n < H->Norb; n++){
        std::string name1 = "/KPM0_" + std::to_string(n);
        std::string name2 = "/KPM1_" + std::to_string(n);
        get_hdf5(KPM0.KPM[n].data(), file, name1);
        get_hdf5(KPM1.KPM[n].data(), file, name2);
    }
    file->close();
    delete file;

    debug_message("Left chebinator::load\n");
}

void chebinator::calc_finish(){
    debug_message("Entered chebinator::calc_finish()\n");
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
    const char *id = getenv("PBS_JOBID");
    if(id!=NULL){
        file << "id=" << std::string(id);
    }
    file << " Lx,Ly=" << Lx << "," << Ly;
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
