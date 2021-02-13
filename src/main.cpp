#include <iostream>
#include <functional>
#include <csignal>
#include <exception>
#include <chrono>
#include <thread>
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

int calculation(int argc, char **argv){

    // Parse input from the command line
    // P is a class containing the parsed inputs
    parameters P;
    parse_input(argc, argv, &P);

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
    print_log_info(P);
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





    unsigned n_threads;
    n_threads = P.nx*P.ny;
    unsigned Norb = 2;
    omp_set_num_threads(n_threads);

    // Parameters that are going to be shared among all threads: the sides and 
    // corners of the lattice sections
    T ***corners;
    Eigen::Array<T, -1, -1> ***sides;
    corners = new T**[n_threads];
    sides   = new Eigen::Array<T, -1, -1>**[n_threads];
    for(unsigned t = 0; t < n_threads; t++){
        corners[t] = new T*[Norb];
        sides[t] = new Eigen::Array<T, -1, -1>*[Norb];

        for(unsigned orb = 0; orb < Norb; orb++){
            corners[t][orb] = new T[4];
            sides[t][orb] = new Eigen::Array<T, -1, -1>[4];
        }
    }


    Eigen::Array<TR, -1, -1> Global_MU;
    Global_MU = Eigen::Array<TR, -1, -1>::Zero(P.nmoments, 1);


#pragma omp parallel shared(corners, sides, Global_MU) firstprivate(gauge_matrix, P, Norb, n_threads)
    {
    unsigned id; 
    unsigned size_x, size_y;
    unsigned lx, ly;
    unsigned tx, ty;

    id = omp_get_thread_num();
    tx = id%P.nx;
    ty = id/P.nx;
    //std::cout << "tx,ty:" << tx << " " << ty << "\n" << std::flush;

    size_y = P.Ly/double(P.ny) + 0.99;
    size_x = P.Lx/double(P.nx) + 0.99;

    ly = std::min(size_y, P.Ly - ty*size_y);
    lx = std::min(size_x, P.Lx - tx*size_x);

    // Set the Hamiltonian object and initialize it with
    // the values obtained from the command line
    hamiltonian H;

    H.Lattice_Lx = P.Lx;
    H.Lattice_Ly = P.Ly;
    H.N_threads_x = P.nx;
    H.N_threads_y = P.ny;
    
    H.set_geometry(lx, ly);
    H.set_regular_hoppings();
    H.set_anderson_W(P.anderson_W/SCALE);
    H.set_anderson();
    H.set_peierls(gauge_matrix);

#pragma omp barrier


    chebinator C;
    C.corners = corners;
    C.sides = sides;
    C.set_hamiltonian(&H);

    C.set_seed(P.seed*n_threads + id);
    C.cheb_iteration(P.nmoments, 1, 1);
#pragma omp barrier
#pragma omp critical
    {
    Global_MU += C.mu.real()/n_threads;
    }
#pragma omp barrier
    } // End parallelization


    //if(P.need_read){
        //verbose2("Need read\n");
        //KPM_vector KPM0, KPM1;
        //KPM0.set_geometry(H.Lx, H.Ly, H.Norb);
        //KPM1.set_geometry(H.Lx, H.Ly, H.Norb);
        //Eigen::Array<TR, -1, -1> mu1;
        //load(P.filename_read, &KPM0, &KPM1, &mu1);

        //C.cheb_iteration_restart(P.moremoments, KPM0, KPM1, mu1);
    //} else {
        //verbose2("Doesn't need read\n");
        //C.cheb_iteration(P.nmoments, P.ndisorder, P.nrandom);
    //}

    //if(P.need_write){
        //C.save(P.filename_write);
    //}




    //Eigen::Array<TR, -1, -1> mu;
    //mu = C.mu;
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
        mu_list[0] = Global_MU.block(0,0,total_moments,1);
        mu_list[1] = Global_MU.block(0,0,total_moments/2,1);

        for(unsigned i = 0; i < N_lists; i++){
            dos_list[i] = calc_dos(mu_list[i], P.energies/SCALE, "jackson");
        }
        delete []mu_list;
        verbose2("Finished calculating DoS\n");


        print_dos(P, N_lists, nmoments_list, dos_list);
        delete []dos_list;
        delete []nmoments_list;
    }


    for(unsigned t = 0; t < n_threads; t++){
        for(unsigned orb = 0; orb < Norb; orb++){
            delete []corners[t][orb];
            delete []sides[t][orb];
        }
        delete []corners[t];
        delete []sides[t];
    }
    delete []corners;
    delete []sides;
    return 0;
}



int main(int argc, char **argv){

    calculation(argc, argv);
    return 0;
}


