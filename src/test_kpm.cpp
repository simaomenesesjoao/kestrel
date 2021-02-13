#include <iostream>
#include <omp.h>
#include <Eigen/Dense>
#include <H5Cpp.h>
#include <H5Group.h>
#include <chrono>

#include "general.hpp"
#include "kpm_vector.hpp"
#include "aux.hpp"
#include "hamiltonian.hpp"
#include "cheb.hpp"
#include "ComplexTraits.hpp"
#include "myHDF5.hpp"

#define NUM_MOMENTS 2000

int main(){

    unsigned nx, ny, n_threads;
    nx = NX_T;
    ny = NY_T;
    n_threads = nx*ny;
    unsigned Norb = 2;


    omp_set_num_threads(n_threads);
    
    // (in order) upper left (ul), upper right (ur), down left (dl), down right (dr)
    // (in order) left, right, up, down
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
    Eigen::Array<TR, -1, -1> MU;

    MU = Eigen::Array<TR, -1, -1>::Zero(NUM_MOMENTS, 1);

#pragma omp parallel shared(corners, sides, MU)
    {
    unsigned id = omp_get_thread_num();
    unsigned Lx, Ly;
    unsigned mult;
    mult = 1;
    Lx = 15;
    Ly = 17;
    unsigned size_y = Ly/double(ny) + 0.9;
    unsigned size_x = Lx/double(nx) + 0.9;
    unsigned lx, ly, Norb;
    unsigned tx, ty;
    Norb = 2;
    tx = id%nx;
    ty = id/nx;
    ly = std::min(size_y, Ly - ty*size_y);
    lx = std::min(size_x, Lx - tx*size_x);
//#pragma omp barrier
//#pragma omp critical
    //{
        //std::cout << "sizes: " << size_x << " " << size_y << "\n" << std::flush;
        //std::cout << "lx,ly: " << lx << " " << ly << "\n" << std::flush;
    //std::cout << "id:" << id << " tx,ty:" << tx << " " << ty << "\n" << std::flush;
    //}

#pragma omp barrier

    Eigen::Matrix<int, 2, 2> gauge_matrix;
    int min_flux;
    int M12, M21;
    extended_euclidean(Lx, Ly, &M12, &M21, &min_flux);
    gauge_matrix(0,0) = 0;
    gauge_matrix(1,1) = 0;
    gauge_matrix(0,1) = -M12*mult;
    gauge_matrix(1,0) = M21*mult;
//#pragma omp barrier
//#pragma omp critical
    //{
    //std::cout << min_flux << " \n" << gauge_matrix << "\n\n"<< std::flush;
    //}
//#pragma omp barrier

    hamiltonian H;
    double W = 0.0;
    H.set_geometry(lx, ly);
    H.Lattice_Lx = Lx;
    H.Lattice_Ly = Ly;
    H.N_threads_x = nx;
    H.N_threads_y = ny;

    H.set_regular_hoppings();
    H.set_anderson_W(W);
    H.set_anderson();
    H.set_peierls(gauge_matrix);

    chebinator C;
    C.n_threads_x = NX_T;
    C.n_threads_y = NY_T;
    C.corners = corners;
    C.sides = sides;
    C.set_hamiltonian(&H);

    // There two are not required for the computation, but are used 
    // to obtain output information about the program
    C.set_seed(1*n_threads + id);
    C.set_mult(1);
    C.cheb_iteration(NUM_MOMENTS, 1, 1);
#pragma omp barrier
#pragma omp critical
    {
    MU += C.mu.real()/n_threads;
    }
#pragma omp barrier
    }

    Eigen::Array<TR, -1, -1> dos_list;
    Eigen::Array<TR, -1, -1> energies;
    //std::cout << "MU:\n" << MU << "\n";
    double lim = 0.99;
    unsigned N_Energies = 5000;
    energies = Eigen::Array<TR, -1, 1>::LinSpaced(N_Energies, -lim, lim);
    dos_list = calc_dos(MU, energies, "jackson");
    std::cout << "##################\n";
    for(unsigned i = 0; i < N_Energies; i++){
        std::cout << energies(i)*SCALE << " " << dos_list(i) << "\n";
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
