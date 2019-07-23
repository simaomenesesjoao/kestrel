#include <iostream>
#include <functional>
#include <csignal>
#include <exception>
#include <chrono>
#include <algorithm>
#include <fstream>
#include <Eigen/Dense>

#include "general.hpp"
#include "kpm_vector.hpp"
#include "aux.hpp"
#include "hamiltonian.hpp"
#include "cheb.hpp"

//
// https://stackoverflow.com/questions/11468414/using-auto-and-lambda-to-handle-signal
//


namespace {
    std::function<void(int)> shutdown_handler;
    void signal_handler(int signal) { shutdown_handler(signal); }
} 


int main(int argc, char **argv){

    unsigned Lx, Ly, num_random, moments, mult, num_disorder;
    double W;
    int seed = -1;

    Eigen::Array<TR, -1, -1> energies;
    parse_input(argc, argv, &Lx, &Ly, &num_random, &moments, &mult, &seed, &energies, &W, &num_disorder);
    W = W/SCALE;
    energies = energies/SCALE;

    energies = Eigen::Array<TR, -1, 1>::LinSpaced(50000, -0.99, 0.99);

    
    int min_flux;
    int M12, M21;
    extended_euclidean(Lx, Ly, &M12, &M21, &min_flux);
    


    // Print the parameters gotten from the command line
    verbose1("_____Compilation parameters_____\n");
    verbose1("STRIDE: " << STRIDE << "\n");
    verbose1("VERBOSE: " << VERBOSE << "\n");

    verbose1("______System parameters_______\n");
    verbose1("Physical system:     Graphene\n");
    verbose1("Lx:                  " + std::to_string(Lx)           + "\n");
    verbose1("Ly:                  " + std::to_string(Ly)           + "\n");
    verbose1("num random:          " + std::to_string(num_random)   + "\n");
    verbose1("num disorder:        " + std::to_string(num_disorder) + "\n");
    verbose1("num moments:         " + std::to_string(moments)      + "\n");
    verbose1("mult:                " + std::to_string(mult)         + "\n");
    verbose1("seed:                " + std::to_string(seed)         + "\n");
    verbose1("anderson:            " + std::to_string(W*SCALE)      + "\n");
    if(energies.size() > 10){
        verbose1("list of energies:    ");
        verbose1(energies.topRows(2).transpose()*SCALE);
        verbose1(" ... (" + std::to_string(energies.size() - 4) + " more) ... ");
        verbose1(energies.bottomRows(2).transpose()*SCALE);
        verbose1("\n");
    } else {
        verbose1("list of energies:    " << energies.transpose()*SCALE     << "\n");
    }


    verbose1("\nInformation about the magnetic field:\n");
    verbose1("Minimum flux:        " + std::to_string(min_flux) + "\n");
    verbose1("Flux:                " + std::to_string(min_flux) + "\n");
    verbose1("Gauge matrix M(0,1): " + std::to_string(-M12*int(mult)) + "\n");
    verbose1("Gauge matrix M(1,0): " + std::to_string(M21*int(mult)) + "\n");


    if(seed != -1){
        std::srand(seed);
    } else {
        std::srand((unsigned int) time(0));
    }




    Eigen::Matrix<int, 2, 2> gauge_matrix;
    gauge_matrix(0,0) = 0;
    gauge_matrix(1,1) = 0;
    gauge_matrix(0,1) = -M12*int(mult);
    gauge_matrix(1,0) = M21*int(mult);

    //std::cout << "gauge:\n" << gauge_matrix << "\n";
    hamiltonian H;
    H.set_geometry(Lx, Ly);
    H.set_regular_hoppings();
    H.set_anderson_W(W);
    H.set_peierls(gauge_matrix);
    //double time = H.time_H();


    chebinator C;
    std::signal(SIGUSR1, signal_handler);
    shutdown_handler = [&](int signal) { C.calc_finish(); };
    C.set_hamiltonian(&H);

    // There two are not required for the computation, but are used 
    // to obtain output information about the program
    C.set_seed(seed);
    C.set_mult(mult);
    //double estimate = C.get_estimate(time, moments, num_disorder, num_random);
    //verbose1("Estimated time for completion: " << estimate << "\n");
    C.cheb_iteration(moments, num_disorder, num_random);




    Eigen::Array<TR, -1, -1> mu(moments,1);
    mu = C.mu;
    verbose2("Finished Chebychev iterations\n");


    // shorten the mu matrix if so desired (to assess convergence)
    unsigned moments_trunc = moments/2;
    Eigen::Array<TR, -1, -1> mu_trunc(moments_trunc, 1);
    mu_trunc = mu.block(0, 0, moments_trunc, 1);

    // Calculate the density of states
    Eigen::Array<TR, -1, -1> dos_green1, dos_green1_trunc, dos_green2, dos_green2_trunc, dos_jack, dos_jack_trunc;


    // These values are in units of t
    //double S1 = 0.003/SCALE; 
    //double S2 = 0.015/SCALE;

    //dos_green1          = calc_dos(mu,       energies, "green " + std::to_string(S1));
    //dos_green1_trunc    = calc_dos(mu_trunc, energies, "green " + std::to_string(S1));
    //dos_green2          = calc_dos(mu,       energies, "green " + std::to_string(S2));
    //dos_green2_trunc    = calc_dos(mu_trunc, energies, "green " + std::to_string(S2));
    dos_jack            = calc_dos(mu,       energies, "jackson");
    dos_jack_trunc      = calc_dos(mu_trunc, energies, "jackson");

    energies *= SCALE;

    //std::cout << "dos: " << dos << "\n";
    for(unsigned i = 0; i < dos_jack.size(); i++){
        std::string metadata;
        metadata  = "Lx:" + std::to_string(Lx);
        metadata += " Ly:" + std::to_string(Ly);
        metadata += " Bmult:" + std::to_string(mult);
        metadata += " nrandom:" + std::to_string(num_random);
        metadata += " ndisorder:" + std::to_string(num_disorder);
        metadata += " seed:" + std::to_string(seed);


        std::cout << metadata << " N:" << moments << " en:" << energies(i) << " dos:" << dos_jack(i) << "\n";
        std::cout << metadata << " N:" << moments_trunc << " en:" << energies(i) << " dos:" << dos_jack_trunc(i) << "\n";
    }
    
     //Save to file
    std::ofstream file;
    //file.open("dos_green1.dat");
    unsigned N_energies = energies.size();
    //for(unsigned i = 0; i < N_energies; i++)
        //file << energies(i) << " " << dos_green1(i) << "\n";
    //file.close();

    //file.open("dos_green1_trunc.dat");
    //for(unsigned i = 0; i < N_energies; i++)
        //file << energies(i) << " " << dos_green1_trunc(i) << "\n";
    //file.close();

    //file.open("dos_green2.dat");
    //for(unsigned i = 0; i < N_energies; i++)
        //file << energies(i) << " " << dos_green2(i) << "\n";
    //file.close();

    //file.open("dos_green2_trunc.dat");
    //for(unsigned i = 0; i < N_energies; i++)
        //file << energies(i) << " " << dos_green2_trunc(i) << "\n";
    //file.close();

    file.open("dos_W" + std::to_string(W*SCALE) + "_B" + std::to_string(mult) + ".dat");
    for(unsigned i = 0; i < N_energies; i++)
        file << energies(i) << " " << dos_jack(i) << "\n";
    file.close();

    file.open("dos_trunc_W" + std::to_string(W*SCALE) + "_B" + std::to_string(mult) + ".dat");
    for(unsigned i = 0; i < N_energies; i++)
        file << energies(i) << " " << dos_jack_trunc(i) << "\n";
    file.close();
    return 0;
}

