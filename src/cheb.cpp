#include <Eigen/Core>
#include <fstream>
#include <chrono>
#include <iostream>
#include <complex>

#include "general.hpp"
#include "kpm_vector.hpp"
#include "aux.hpp"
#include "hamiltonian.hpp"
#include "cheb.hpp"

#include <H5Cpp.h>
#include <H5Group.h>
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"


chebinator::chebinator(){
    debug_message("Entered chebinator constructor\n");
    initialized = false;

    debug_message("Left chebinator constructor\n");
};

chebinator::~chebinator(){
    debug_message("Entered chebinator destructor\n");
    debug_message("Left chebinator destructor\n");
};

void chebinator::set_hamiltonian(hamiltonian *new_H){
    debug_message("Entered chebinator set hamiltonian\n");
    H = new_H;
    Lx = H->Lx;
    Ly = H->Ly;
    
    n_threads_x = H->N_threads_x;
    n_threads_y = H->N_threads_y;

    Norb = H->Norb;

    n_threads = omp_get_num_threads();
    //std::cout << "num threads: " << n_threads << "\n" << std::flush;

    thread_id = omp_get_thread_num();
    KPM0.thread_id = thread_id;
    KPM1.thread_id = thread_id;
    KPM_initial.thread_id = thread_id;

    KPM0.NTx = n_threads_x;
    KPM0.NTy = n_threads_y;
    KPM0.Npoints_lattice = H->Lattice_Lx*H->Lattice_Ly*Norb;

    KPM1.NTx = n_threads_x;
    KPM1.NTy = n_threads_y;
    KPM1.Npoints_lattice = H->Lattice_Lx*H->Lattice_Ly*Norb;

    KPM_initial.NTx = n_threads_x;
    KPM_initial.NTy = n_threads_y;
    KPM_initial.Npoints_lattice = H->Lattice_Lx*H->Lattice_Ly*Norb;

    KPM0.corners = corners;
    KPM1.corners = corners;
    KPM_initial.corners = corners;

    KPM0.sides = sides;
    KPM1.sides = sides;
    KPM_initial.sides = sides;

#pragma omp barrier

    KPM0.set_geometry(Lx, Ly, Norb);
    KPM1.set_geometry(Lx, Ly, Norb);
    KPM_initial.set_geometry(Lx, Ly, Norb);
    debug_message("Left chebinator set hamiltonian\n");

}



void chebinator::cheb_iteration_restart(unsigned n, KPM_vector &KPM0_restart, KPM_vector &KPM1_restart, Eigen::Array<TR, -1, -1> mu_restart){
//void chebinator::cheb_iteration_restart(unsigned n, Eigen::Array<TR, -1, -1> mu_restart){
    debug_message("Entered cheb_iteration\n");


    // set some useless variables that are only going to be needed
    // for displaying the output with system signals
    current_ndis = 1;
    current_nrand = 1;



    start = std::chrono::high_resolution_clock::now();



    unsigned previous_moments = mu_restart.size();
    unsigned moments = n;
    unsigned total_moments = previous_moments + moments;
    current_moments = total_moments;
    max_iter = moments;


    // These two steps are needed so that the disorder realization
    // and the initial random vectors are the same
    H->set_anderson();
    KPM0.random_uniform();
    KPM_initial = KPM0;
    KPM_initial.empty_ghosts();

    //load(name);

    KPM0 = KPM0_restart;
    KPM1 = KPM1_restart;
    mu = Eigen::Array<TR, -1, -1>::Zero(total_moments, 1);
    mu.block(0,0,previous_moments,1) = mu_restart;


    // only after the seed has been used, load the vectors



    std::chrono::duration<double> time_span;
    for(unsigned i = previous_moments; i < total_moments; i++){
        debug_message("cheb iteration " + std::to_string(i) + "\n");
        H->cheb(KPM1, KPM0, i);

        swap_pointers(&KPM1, &KPM0);
        KPM0.fill_ghosts();
        mu(i) = KPM_initial*KPM0;

        current_iter = i - previous_moments;
        current = std::chrono::high_resolution_clock::now();

        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(current - start);
        //avg_iter_time = time_span.count()/(1+current_iter);

        //vars->avg_time = time_span.count()/(1+current_iter); 
        //vars->iter = current_iter;
        //vars->update_status();
    }
    debug_message("Left cheb_iteration\n");
}

void chebinator::cheb_iteration(unsigned n, unsigned ndis, unsigned nrand){ 
    debug_message("Entered cheb_iteration\n");

    //std::chrono::duration<double> time_span;

    // set some useless variables that are only going to be needed
    // for displaying the output with system signals
    current_ndis = ndis;
    current_moments = n;
    current_nrand = nrand;


    //vars->max_iter = n*nrand*ndis;
    //sprintf(vars->buffer, "boas");


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
            //KPM0.site(0, 1, 1);

            unsigned r1,r2;
            debug_message("A vacancies\n");
            for(unsigned i = 0; i < H->NvacA; i++){
                r1 = H->vacanciesA(i,0);
                r2 = H->vacanciesA(i,1);
                debug_message(i); debug_message(" ");
                debug_message(r1); debug_message(" ");
                debug_message(r2) debug_message("\n");;
                KPM0.KPM[0](r1+1,r2+1) = 0;
            }
            debug_message("B vacancies\n");
            for(unsigned i = 0; i < H->NvacB; i++){
                r1 = H->vacanciesB(i,0);
                r2 = H->vacanciesB(i,1);
                debug_message(i); debug_message(" ");
                debug_message(r1); debug_message(" ");
                debug_message(r2); debug_message("\n");
                KPM0.KPM[1](r1+1,r2+1) = 0;
            }
            debug_message("Finished putting vacancies");
            static unsigned totVac = 0;
#pragma omp barrier
#pragma omp critical
            {
                totVac += H->NvacA + H->NvacB;
            }
#pragma omp barrier
            // Fix normalization of the vectors
#pragma omp critical
            {
                KPM0.KPM[0] *= sqrt(KPM0.Npoints_lattice)/sqrt(KPM0.Npoints_lattice-totVac);
                KPM0.KPM[1] *= sqrt(KPM0.Npoints_lattice)/sqrt(KPM0.Npoints_lattice-totVac);

            }

            static double tot_norm = 0;
#pragma omp barrier
#pragma omp critical
            {
                tot_norm += KPM0.KPM[0].block(1,1,H->Ly,H->Lx).abs2().sum();
                tot_norm += KPM0.KPM[1].block(1,1,H->Ly,H->Lx).abs2().sum();
                //std::cout << "element:" << KPM0.KPM[0](1,1) << "\n";
            }
#pragma omp barrier
//#pragma omp master
            //{
                    //std::cout << "NORM:" << tot_norm << "\n" << std::flush;
                    //std::cout << "tot vac: " << totVac << "\n" << std::flush;
            //}
//#pragma omp barrier

            KPM1.zero();
            //KPM0.zero();
            //KPM0.KPM[0] += 1;
            //KPM0.KPM[1] += 1;
            KPM0.fill_ghosts();
            KPM_initial = KPM0;
            KPM_initial.empty_ghosts();
            //std::cout << "KPM_initial:\n" << KPM_initial.KPM[0] << "\n";

            //start = std::chrono::high_resolution_clock::now();
            //double time_count = 0.0;
            for(unsigned i = 0; i < moments-1; i++){
                debug_message("cheb iteration " + std::to_string(i) + "\n");
                H->cheb(KPM1, KPM0, i);

                swap_pointers(&KPM1, &KPM0);
                KPM0.fill_ghosts();
                mu(i+1) += (KPM_initial*KPM0)/nrand/ndis;

                //current = std::chrono::high_resolution_clock::now();
                //time_span = std::chrono::duration_cast<std::chrono::duration<double>>(current - start);
                current_iter = i + r*moments + d*moments*nrand;
                //vars->avg_time = time_span.count()/(1+current_iter); 
                //vars->iter = current_iter;
                //vars->update_status();

            }
#pragma omp barrier
#pragma omp critical
            {
            //std::cout << time_span.count() << "\n" << std::flush;
            }
#pragma omp barrier
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
    double anderson_w = H->W*SCALE;
    write_hdf5(&anderson_w, file1, "/anderson_w");
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


double chebinator::get_estimate(double time, unsigned moments, unsigned num_disorder, unsigned num_random){
    return time*moments*num_disorder*num_random;
}

void chebinator::set_geometry(unsigned lx, unsigned ly, unsigned norb){
    Lx = lx;
    Ly = ly;
    Norb = norb;

    // Now that we know the number of orbitals, allocate
    // the memory needed for all the KPM vectors
    KPM0.set_geometry(Lx, Ly, norb);
    KPM1.set_geometry(Lx, Ly, norb);
    KPM_initial.set_geometry(Lx, Ly, norb);

}
