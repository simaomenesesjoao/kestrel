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
    Eigen::Array<int, 2, 2> gauge_matrix;
    gauge_matrix(0,0) = 0;
    gauge_matrix(1,1) = 0;
    gauge_matrix(0,1) = M12*P.mult;
    gauge_matrix(1,0) = -M21*P.mult; //don't forget this sign has to be minus
    
    // Print the parameters gotten from the command line
    print_compilation_info();
    print_ham_info(P);
    print_magnetic_info(P, M12, M21, min_flux);
    print_cheb_info(P);
    print_log_info(P);
    print_output_info(P);
    print_restart_info(P);




    


    //
    //
    //
    // Information about the system
    unsigned dim   = 2;
    unsigned Norb  = 1;
    unsigned Nhops = 4;
    double scale = 4.1;

    // Information about the primitive vectors
    Eigen::Array<TR, 2, 1> a1, a2;
    a1 << 1.0, 0.0;
    a2 << 0.0, 1.0;



    // Information about the hoppings
    Eigen::Array<int, -1, 2> origin;
    Eigen::Array<unsigned, 1,-1> to_matrix, from_matrix;
    from_matrix = Eigen::Array<unsigned, 1, -1>::Zero(1,Nhops);
    to_matrix   = Eigen::Array<unsigned, 1, -1>::Zero(1,Nhops);
    origin      = Eigen::Array<int, -1, -1>::Zero(Nhops,dim);

    origin << 1,0,
              -1,0,
              0,0,
              0,-1;

    from_matrix << 0,0,0,0;
    to_matrix   << 0,0,0,0;


    // Vectors between two carbons inside the same unit cell,
    // in units of the lattice vectors
    Eigen::Array<TR, 2, 1> delta;
    delta << 0.0, 0.0;

    // positions for the right hand side
    Eigen::Array<double, -1, -1> orb_pos;
    orb_pos = Eigen::Array<double, -1,-1>::Zero(Norb,dim);
    orb_pos.row(0) = delta*0.0;
    //orb_pos.row(1) = delta


















    // set the global seed
    // NOTE: should have separate seed for disorder and
    //       random realisation of KPM vector
    if(P.seed != -1){
        std::srand(P.seed);
    } else {
        std::srand((unsigned int) time(0));
    }

    // Vacancies 
    //double conc = 0.01; //1%
    double conc = P.conc;
    unsigned Nvac = conc*P.Ly*P.Lx*Norb;
    unsigned Nsites = P.Ly*P.Lx*Norb;
    if(Nvac > Nsites){
        std::cout << "There are more vacancies than sites\n";
        exit(1);
    }


    // Create a list of vacancies, no superpositions
    unsigned *occupancies;
    unsigned *vac_sites;

    occupancies = new unsigned[Nsites];
    vac_sites = new unsigned[Nvac];

    for(unsigned i = 0; i < Nsites; i++){ occupancies[i] = 0; }


    unsigned total = 0;
    unsigned prop;
    while(total < Nvac){
        prop = rand() % Nsites;
        if(occupancies[prop] == 0){
            // Accept proposal
            vac_sites[total] = prop;
            occupancies[prop] = 1;
            total++;
        }
    }
    delete[] occupancies;

    // print
    //for(unsigned i = 0; i < Nvac; i++){ std::cout << vac_sites[i] << " "; } std::cout << "\n";

    // Check
    unsigned a, b;
    bool double_vac = false;
    for(unsigned i = 0; i < Nvac; i++){
        a = vac_sites[i];
        for(unsigned j = i+1; j < Nvac; j++){
            b = vac_sites[j];
            if(b==a){
                double_vac = true;
            }
        }
    }
    if(double_vac){
        std::cout << "Two vacancies exist in the same place\n";
        exit(1);
    }










    unsigned n_threads;
    n_threads = P.nx*P.ny;
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


#pragma omp parallel shared(Nvac, corners, sides, vac_sites, Global_MU) firstprivate(gauge_matrix, P, Norb, n_threads)
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


    unsigned NvacA = 0;
    unsigned NvacB = 0;
    Eigen::Array<unsigned, -1, -1> vacAtemp, vacBtemp, vacA, vacB;
    vacAtemp = Eigen::Array<unsigned, -1, -1>::Zero(Nvac,dim);
    vacBtemp = Eigen::Array<unsigned, -1, -1>::Zero(Nvac,dim);
    unsigned orb, vacpos, r1, r2;
#pragma omp critical
    {
        //std::cout << "thread: " << omp_get_thread_num() << "\n";
        for(unsigned i = 0; i < Nvac; i++){
            vacpos = vac_sites[i]; // vacpos = Lx*Ly * o + Lx*y + x
            r1 = vacpos % P.Lx;
            r2 = ((vacpos - r1)/P.Lx ) % P.Ly;
            orb = (vacpos - r1 - r2*P.Lx)/(P.Lx*P.Ly);

            unsigned r1_tx, r2_ty;
            // Which thread should it go to?
            r1_tx = r1 / size_x;
            r2_ty = r2 / size_y;

            // coordinates local
            unsigned r1_loc, r2_loc;
            r1_loc = r1 - r1_tx*size_x;
            r2_loc = r2 - r2_ty*size_y;

            //std::cout << "vac:" << vacpos << "    r1,r2,o: " << r1 << " " << r2 << " " << orb << "    tx,ty:" << r1_tx << "," << r2_ty << "    loc:" << r1_loc << "," << r2_loc << "\n";
            if(r1_tx == tx && r2_ty == ty){
                if(orb == 0){
                    vacAtemp(NvacA,0) = r1_loc;
                    vacAtemp(NvacA,1) = r2_loc;
                    NvacA++;
                }
                if(orb == 1){
                    vacBtemp(NvacB,0) = r1_loc;
                    vacBtemp(NvacB,1) = r2_loc;
                    NvacB++;
                }
            }
        }

    vacA = Eigen::Array<unsigned, -1, -1>::Zero(NvacA,2);
    vacB = Eigen::Array<unsigned, -1, -1>::Zero(NvacB,2);
    for(unsigned i = 0; i < NvacA; i++){ vacA.row(i) = vacAtemp.row(i); }
    for(unsigned i = 0; i < NvacB; i++){ vacB.row(i) = vacBtemp.row(i); }
    }
#pragma omp barrier



    // Set the Hamiltonian object and initialize it with
    // the values obtained from the command line
    hamiltonian H(Norb, Nhops);

    H.Lattice_Lx = P.Lx;
    H.Lattice_Ly = P.Ly;
    H.N_threads_x = P.nx;
    H.N_threads_y = P.ny;
    
    H.set_geometry(lx, ly);
    H.set_regular_hoppings();
    H.scale = scale;

    H.set_primitive2(a1,a2);
    H.set_origin_to_from(origin, to_matrix, from_matrix);
    H.set_orbpos(orb_pos);

    H.set_anderson_W(P.anderson_W/SCALE);
    H.set_anderson();
    H.set_vacancies(vacA, vacB);
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
    Global_MU += C.mu.real();
    }
#pragma omp barrier
    } // End parallelization
    Global_MU(0) = 1;
    //std::cout << "Global Mu\n" << Global_MU << "\n";


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

    // Clean up the vacancy list
    delete[] vac_sites;
    return 0;
}



int main(int argc, char **argv){

    calculation(argc, argv);
    return 0;
}


