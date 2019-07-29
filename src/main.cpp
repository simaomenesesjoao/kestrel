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

void print_dos(parameters P, unsigned N_lists, unsigned *nmoments_list, Eigen::Array<TR, -1, -1> *dos);

namespace {
    std::function<void(int)> shutdown_handler;
    void signal_handler(int signal) { shutdown_handler(signal); }
} 



int main(int argc, char **argv){
    std::signal(SIGUSR1, signal_handler);
    shutdown_handler = [&](int signal) {  };

    // Parse input from the command line
    parameters P;
    parse_input(argc, argv, &P);

    // Rescale the relevant quantities
    //P.anderson_W = P.anderson_W/SCALE;
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
    H.set_anderson_W(P.anderson_W/SCALE);
    H.set_peierls(gauge_matrix);


    // Chebyshev object to iterate KPM vectors with the
    // Hamiltonian object
    chebinator C;
    C.set_hamiltonian(&H);

    // Set the external signal handler to use a method from the 
    // Chebyshev object
    shutdown_handler = [&](int signal) { C.calc_finish(); };

    // There two are not required for the computation, but are used 
    // to obtain output information about the program
    C.set_seed(P.seed);
    C.set_mult(P.mult);

    // Estimate the time for completion
    double time = H.time_H();
    unsigned moments_to_process = P.nmoments;
    if(P.need_read) moments_to_process = P.moremoments;
    double estimate = C.get_estimate(time, moments_to_process, P.ndisorder, P.nrandom);
    verbose1("Estimated time for completion: " << estimate << " seconds.\n");

    if(P.need_read){
        verbose2("Need read\n");
        KPM_vector KPM0, KPM1;
        KPM0.set_geometry(H.Lx, H.Ly, H.Norb);
        KPM1.set_geometry(H.Lx, H.Ly, H.Norb);
        Eigen::Array<TR, -1, -1> mu1;
        load(P.filename_read, &KPM0, &KPM1, &mu1);

        C.cheb_iteration_restart(P.moremoments, KPM0, KPM1, mu1);
    } else {
        verbose2("Doesn't need read\n");
        C.cheb_iteration(P.nmoments, P.ndisorder, P.nrandom);
    }

    if(P.need_write){
        C.save(P.filename_write);
    }




    Eigen::Array<TR, -1, -1> mu;
    mu = C.mu;
    verbose2("Finished Chebychev iterations\n");



    if(P.need_print){
        const unsigned N_lists = 2;
        Eigen::Array<TR, -1, -1> *dos_list, *mu_list;
        unsigned *nmoments_list;
        dos_list = new Eigen::Array<TR, -1, -1>[N_lists];
        mu_list  = new Eigen::Array<TR, -1, -1>[N_lists];
        nmoments_list = new unsigned[N_lists];

        unsigned total_moments = P.nmoments;
        if(P.need_read) total_moments += P.moremoments;

        nmoments_list[0] = total_moments;
        nmoments_list[1] = total_moments/2;
        mu_list[0] = mu.block(0,0,total_moments,1);
        mu_list[1] = mu.block(0,0,total_moments/2,1);

        for(unsigned i = 0; i < N_lists; i++){
            dos_list[i] = calc_dos(mu_list[i], P.energies, "jackson");
        }
        delete []mu_list;
        verbose2("Finished calculating DoS\n");


        print_dos(P, N_lists, nmoments_list, dos_list);
        delete []dos_list;
        delete []nmoments_list;
    }
    return 0;
}


void print_dos(parameters P, unsigned N_lists, unsigned *nmoments_list, Eigen::Array<TR, -1, -1> *dos){
    unsigned N_energies = dos[0].size(); // They all have the same number of points

    if(P.need_print_to_cout){
        std::cout << "________Density of states________\n";
        for(unsigned i = 0; i < N_energies; i++){
            std::string metadata;
            metadata  = "Lx:" + std::to_string(P.Lx);
            metadata += " Ly:" + std::to_string(P.Ly);
            metadata += " Bmult:" + std::to_string(P.mult);
            metadata += " nrandom:" + std::to_string(P.nrandom);
            metadata += " ndisorder:" + std::to_string(P.ndisorder);
            metadata += " W:" + std::to_string(P.anderson_W);
            metadata += " seed:" + std::to_string(P.seed);

            for(unsigned j = 0; j < N_lists; j++){
                std::cout << metadata << " N:" << nmoments_list[j] << " en:" << P.energies(i)*SCALE << " dos:" << dos[j](i) << "\n";
            }
        }
    }
    
    if(P.need_print_to_file){

         //Save to file
        std::ofstream file;

        for(unsigned j = 0; j < N_lists; j++){
            std::string name;
            name = "dos_N" + std::to_string(nmoments_list[j]) + "_W" + std::to_string(P.anderson_W) + "_B" + std::to_string(P.mult) + ".dat";
            file.open(name);
            for(unsigned i = 0; i < N_energies; i++)
                file << P.energies(i)*SCALE << " " << dos[j](i) << "\n";
            file.close();
        }

    }

    if(!P.need_print) std::cout << "Nothing to print.\n";
}
