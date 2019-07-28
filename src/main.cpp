#include <iostream>
#include <functional>
#include <csignal>
#include <exception>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <Eigen/Dense>
#include <H5Cpp.h>
#include <H5Group.h>

#include "general.hpp"
#include "kpm_vector.hpp"
#include "aux.hpp"
#include "hamiltonian.hpp"
#include "cheb.hpp"
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"

//
// https://stackoverflow.com/questions/11468414/using-auto-and-lambda-to-handle-signal
//


namespace {
    std::function<void(int)> shutdown_handler;
    void signal_handler(int signal) { shutdown_handler(signal); }
} 


int main(int argc, char **argv){

    // Parse input from the command line
    parameters P;
    parse_input(argc, argv, &P);

    // Rescale the relevant quantities
    P.anderson_W = P.anderson_W/SCALE;
    P.energies = P.energies/SCALE;

    // Calculate the compatible magnetic fields and
    // the minimum magnetic flux allowed
    int min_flux;
    int M12, M21;
    extended_euclidean(P.Lx, P.Ly, &M12, &M21, &min_flux);
    Eigen::Matrix<int, 2, 2> gauge_matrix;
    gauge_matrix(0,0) = 0;
    gauge_matrix(1,1) = 0;
    gauge_matrix(0,1) = -M12*P.mult;
    gauge_matrix(1,0) = M21*P.mult;
    
    // Print the parameters gotten from the command line
    print_compilation_info();
    print_ham_info(P);
    print_magnetic_info(P, M12, M21, min_flux);
    print_cheb_info(P);
    print_output_info(P);
    print_restart_info(P);

    // set the global seed
    // NOTE: should have separate seed for disorder and
    //       random realisation of KPM vector
    if(P.seed != -1){
        std::srand(P.seed);
    } else {
        std::srand((unsigned int) time(0));
    }


    // Set the Hamiltonian object and initialize it with
    // the values obtained from the command line
    hamiltonian H;
    H.set_geometry(P.Lx, P.Ly);
    H.set_regular_hoppings();
    H.set_anderson_W(P.anderson_W);
    H.set_peierls(gauge_matrix);
    double time = H.time_H();


    // Chebyshev object to iterate KPM vectors with the
    // Hamiltonian object
    chebinator C;
    C.set_hamiltonian(&H);

    // Set the external signal handler to use a method from the 
    // Chebyshev object
    std::signal(SIGUSR1, signal_handler);
    shutdown_handler = [&](int signal) { C.calc_finish(); };

    // There two are not required for the computation, but are used 
    // to obtain output information about the program
    C.set_seed(P.seed);
    C.set_mult(P.mult);

    // Estimate the time for completion
    double estimate = C.get_estimate(time, P.nmoments, P.ndisorder, P.nrandom);
    verbose1("Estimated time for completion: " << estimate << " seconds.\n");

    if(P.need_read){
        KPM_vector KPM0, KPM1;
        KPM0.set_geometry(H.Lx, H.Ly, H.Norb);
        KPM1.set_geometry(H.Lx, H.Ly, H.Norb);
        Eigen::Array<TR, -1, -1> mu1;
        load(P.filename_read, &C.KPM0, &C.KPM1, &mu1);
        C.cheb_iteration_restart(P.moremoments, KPM0, KPM1, mu1);
    } else {
        C.cheb_iteration(P.nmoments, P.ndisorder, P.nrandom);
    }

    if(P.need_write){
        C.save(P.filename_write);
    }




    Eigen::Array<TR, -1, -1> mu(P.nmoments,1);
    mu = C.mu;
    verbose2("Finished Chebychev iterations\n");




    if(P.output_energies){
        // shorten the mu matrix if so desired (to assess convergence)
        unsigned moments_trunc = P.nmoments/2;
        Eigen::Array<TR, -1, -1> mu_trunc(moments_trunc, 1);
        mu_trunc = mu.block(0, 0, moments_trunc, 1);

        // Calculate the density of states
        Eigen::Array<TR, -1, -1> dos_jack, dos_jack_trunc;
        dos_jack            = calc_dos(mu,       P.energies, "jackson");
        dos_jack_trunc      = calc_dos(mu_trunc, P.energies, "jackson");

        P.energies *= SCALE;

        verbose2("Finished calculating DoS\n");
        for(unsigned i = 0; i < dos_jack.size(); i++){
            std::string metadata;
            metadata  = "Lx:" + std::to_string(P.Lx);
            metadata += " Ly:" + std::to_string(P.Ly);
            metadata += " Bmult:" + std::to_string(P.mult);
            metadata += " nrandom:" + std::to_string(P.nrandom);
            metadata += " ndisorder:" + std::to_string(P.ndisorder);
            metadata += " seed:" + std::to_string(P.seed);


            std::cout << metadata << " N:" << P.nmoments << " en:" << P.energies(i) << " dos:" << dos_jack(i) << "\n";
            std::cout << metadata << " N:" << moments_trunc << " en:" << P.energies(i) << " dos:" << dos_jack_trunc(i) << "\n";
        }
        
         //Save to file
        std::ofstream file;
        //file.open("dos_green1.dat");
        unsigned N_energies = P.energies.size();

        file.open("dos_W" + std::to_string(P.anderson_W*SCALE) + "_B" + std::to_string(P.mult) + ".dat");
        for(unsigned i = 0; i < N_energies; i++)
            file << P.energies(i) << " " << dos_jack(i) << "\n";
        file.close();

        file.open("dos_trunc_W" + std::to_string(P.anderson_W*SCALE) + "_B" + std::to_string(P.mult) + ".dat");
        for(unsigned i = 0; i < N_energies; i++)
            file << P.energies(i) << " " << dos_jack_trunc(i) << "\n";
        file.close();
    }
    return 0;
}

