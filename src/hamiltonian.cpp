#include <iostream>
#include <chrono>
#include <exception>
#include <algorithm>
#include <fstream>
#include <Eigen/Dense>

#include "general.hpp"

#include "kpm_vector.hpp"
#include "aux.hpp"
#include "hamiltonian.hpp"



void hamiltonian::set_peierls(Eigen::Matrix<int, 2, 2> gauge_matrix){
    Eigen::Matrix<TR, -1, -1> A(2,2);
    A(0,0) = gauge_matrix(0,0)*1.0/Lattice_Lx;
    A(1,1) = gauge_matrix(1,1)*1.0/Lattice_Ly;
    A(0,1) = gauge_matrix(0,1)*1.0/Lattice_Lx;
    A(1,0) = gauge_matrix(1,0)*1.0/Lattice_Ly;

    // Lattice vectors
    Eigen::Matrix<TR, 2, 1> a1, a2;
    a1(0) = sqrt(3.0); 
    a1(1) = 0.0;
    a2(0) = sqrt(3.0)/2.0;
    a2(1) = 3.0/2.0;

    // Lattice vectors in units of the lattice vectors
    Eigen::Matrix<TR, 2, 1> a1c, a2c;
    a1c(0) = 1.0; 
    a1c(1) = 0.0;
    a2c(0) = 0.0;
    a2c(1) = 1.0;

    // Vectors between two carbons inside the same unit cell,
    // in units of the lattice vectors
    Eigen::Matrix<TR, 2, 1> delta;
    delta = - 1.0/3.0 * a1c + 2.0/3.0 * a2c;
    //delta = 2.0/3.0 * a1c - 1.0/3.0 * a2c;
    //delta = delta*0.0;

    // neighbour distance for each hopping
    std::complex<TR> im = std::complex<TR>(0, 1.0);
    Eigen::Matrix<TR, 2, 1> d[6];



    
    /*      A     A              A       A
     *       \   /            t5  \     /  t4
     *        \ /                  \   /
     *         B                    >B<              B   
     *         |                     ^               |  t0 
     *         |                     |  t3           v
     *         A                     A             _ A _
     *        / \                             t1   /| |\  t2            ;      
     *       /   \                                /     \
     *      B     B                              B       B 
     */


    d[0] = -delta;
    d[1] = -delta + a2c;
    d[2] = -delta + a2c - a1c;
    d[3] = delta;
    d[4] = delta - a2c;
    d[5] = delta - a2c + a1c;

    // positions for the right hand side
    Eigen::Matrix<TR, 2, 1> pos1[6], pos0[6];
    pos0[0] = delta;
    pos0[1] = delta;
    pos0[2] = delta;
    pos0[3] = delta*0.0;
    pos0[4] = delta*0.0;
    pos0[5] = delta*0.0;

    pos1[0] = delta*0.0;
    pos1[1] = delta*0.0;
    pos1[2] = delta*0.0;
    pos1[3] = delta;
    pos1[4] = delta;
    pos1[5] = delta;


    Eigen::Matrix<TR, 2, 1> origin[6];
    origin[0] = Eigen::Matrix<TR, 2, 1>(0, 0);
    origin[1] = Eigen::Matrix<TR, 2, 1>(0, -1);
    origin[2] = Eigen::Matrix<TR, 2, 1>(1, -1);
    origin[3] = Eigen::Matrix<TR, 2, 1>(0, 0);
    origin[4] = Eigen::Matrix<TR, 2, 1>(0, 1);
    origin[5] = Eigen::Matrix<TR, 2, 1>(-1, 1);

    unsigned id;
    unsigned tx, ty;
    unsigned offset_x, offset_y;
#pragma omp barrier
#pragma omp critical
    {
    id = omp_get_thread_num(); // id = ty*N_threads_x + tx
    tx = id%N_threads_x;
    ty = id/N_threads_x;
    unsigned size_y = Lattice_Ly/double(N_threads_y) + 0.9;
    unsigned size_x = Lattice_Lx/double(N_threads_x) + 0.9;

    offset_x = tx*size_x;
    offset_y = std::max(int(Lattice_Ly)- int((ty + 1)*size_y), 0);
    }

#pragma omp barrier

    double t = 1.0/SCALE;
    for(unsigned orb = 0; orb < 6; orb++){
        for(unsigned j = 0; j < Ly; j++){
            for(unsigned i = 0; i < Lx; i++){


                Eigen::Matrix<TR, 2, 1> r, ri, rf, r_med, r_dif;
                r(0) = offset_x + i; 
                r(1) = offset_y + Ly - 1 - j;

                r += origin[orb];
                
                rf = r + pos0[orb] + d[orb];
                ri = r + pos0[orb];

                r_med = (rf + ri)/2.0;
                r_dif = rf - ri;

                double phase1 = r_med.transpose()*A*r_dif;
                double phase2 = -rf.transpose()*A*pos1[orb];
                double phase3 = -ri.transpose()*A*pos0[orb];

                hoppings[orb](j, i) = exp(2.0*M_PI*im*(phase1 + phase2 - phase3))*t;
            }
        }
    }

}

void hamiltonian::set_geometry(unsigned lx, unsigned ly){
    Lx = lx;
    Ly = ly;
}


hamiltonian::hamiltonian(){
    Norb = 2;
    N_hoppings = 6;
    hoppings = new Eigen::Array<T, -1, -1>[N_hoppings];
    anderson = new Eigen::Array<T, -1, -1>[Norb];

    is_anderson_set = false;
}

hamiltonian::~hamiltonian(){
    debug_message("Entered hamiltonian destructor\n");
    delete []hoppings;
    delete []anderson;
    debug_message("Left hamiltonian destructor\n");
}

void hamiltonian::set_regular_hoppings(){

    for(unsigned i = 0; i < N_hoppings; i++){
        hoppings[i] = Array::Zero(Ly, Lx) + 1.0/SCALE;
    }
}

void hamiltonian::set_anderson_W(double w_anderson){
    W = w_anderson;
}

void hamiltonian::set_anderson(){
    debug_message("Entered hamiltonian::set_anderson\n");
        for(unsigned j = 0; j < Norb; j++){
            anderson[j] = Eigen::Array<TR, -1, -1>::Random(Ly, Lx)*W/2.0; // random numbers between -W/2 and W/2
        }
    is_anderson_set = true;
    debug_message("Left hamiltonian::set_anderson\n");
}

void hamiltonian::set_vacancies(Eigen::Array<unsigned, -1, -1> vacA, Eigen::Array<unsigned, -1, -1> vacB){
    debug_message("Entered hamiltonian::set_vacancies\n");
    vacanciesA = vacA;
    vacanciesB = vacB;
    NvacA = vacanciesA.rows();
    NvacB = vacanciesB.rows();

//#pragma omp critical
    //{
    //unsigned id = omp_get_thread_num();
    //std::cout << "id:" << id << "\n";
    //std::cout << vacanciesA << "\n" << vacanciesB << "\n";
    //}
    //debug_message("Left hamiltonian::set_vacancies\n");
}

void hamiltonian::H(KPM_vector &KPM_new, KPM_vector &KPM_old, unsigned f){
    debug_message("Entered hamiltonian::H\n");
    // Writes the result of H*KPM_old into KPM_new by adding. It does not replace the values
    // This is intentional but may cause some problems if KPM_new is not initialized to zero

    unsigned j = 0;
    unsigned i = 0;
    unsigned final_i = Lx;
    unsigned final_j = Ly;

//#pragma omp barrier
//#pragma omp critical
    //{

    // Anderson disorder for each orbital
    KPM_new.KPM[0].block(1+j, 1+i, final_j, final_i) += KPM_old.KPM[0].block(1+j, 1+i, final_j, final_i)*anderson[0].block(j, i, final_j, final_i)*f;
    KPM_new.KPM[1].block(1+j, 1+i, final_j, final_i) += KPM_old.KPM[1].block(1+j, 1+i, final_j, final_i)*anderson[1].block(j, i, final_j, final_i)*f;

    // Regular hoppings
    KPM_new.KPM[0].block(1+j, 1+i, final_j, final_i) += KPM_old.KPM[1].block(1+j, 1+i, final_j, final_i)*hoppings[0].block(j, i, final_j, final_i)*f;
    KPM_new.KPM[0].block(1+j, 1+i, final_j, final_i) += KPM_old.KPM[1].block(2+j, 1+i, final_j, final_i)*hoppings[1].block(j, i, final_j, final_i)*f;
    KPM_new.KPM[0].block(1+j, 1+i, final_j, final_i) += KPM_old.KPM[1].block(2+j, 2+i, final_j, final_i)*hoppings[2].block(j, i, final_j, final_i)*f;

    KPM_new.KPM[1].block(1+j, 1+i, final_j, final_i) += KPM_old.KPM[0].block(1+j, 1+i, final_j, final_i)*hoppings[3].block(j, i, final_j, final_i)*f;
    KPM_new.KPM[1].block(1+j, 1+i, final_j, final_i) += KPM_old.KPM[0].block(0+j, 1+i, final_j, final_i)*hoppings[4].block(j, i, final_j, final_i)*f;
    KPM_new.KPM[1].block(1+j, 1+i, final_j, final_i) += KPM_old.KPM[0].block(0+j, 0+i, final_j, final_i)*hoppings[5].block(j, i, final_j, final_i)*f;

    // Vacancies
    debug_message("Before vacancies in Hamiltonian\n");
    debug_message("thread number:");
    debug_message(omp_get_thread_num());
    debug_message("\n");
    unsigned r1,r2;
    //std::cout << "size: " << KPM_new.KPM[0].rows() << "," << KPM_new.KPM[0].cols() << "\n";
    //std::cout << "size: " << KPM_new.KPM[1].rows() << "," << KPM_new.KPM[1].cols() << "\n";
    debug_message("A vacancies\n");
    for(unsigned i = 0; i < NvacA; i++){

        r1 = vacanciesA(i,0);
        r2 = vacanciesA(i,1);
        //std::cout << "r1,r2: " << r1 << "," << r2 << "\n" << std::flush;
        KPM_new.KPM[0](r1+1,r2+1) = 0;
    }
    debug_message("B vacancies\n");
    for(unsigned i = 0; i < NvacB; i++){
        r1 = vacanciesB(i,0);
        r2 = vacanciesB(i,1);
        //std::cout << "r1,r2: " << r1 << "," << r2 << "\n" << std::flush;
        KPM_new.KPM[1](r1+1,r2+1) = 0;
    }
    debug_message("After vacancies in Hamiltonian\n");
    //
    debug_message("Left hamiltonian::H\n");
    //}
//#pragma omp barrier
}



void hamiltonian::cheb(KPM_vector &KPM_new, KPM_vector &KPM_old, unsigned it_num){
    debug_message("Entered hamiltonian::cheb\n");
    
    // First Chebyshev iteration
    // this if branch requires KPM_new to be empty because it ADDS to KPM_new
    if(it_num == 0){
        H(KPM_new, KPM_old, 1);
    }

    if(it_num != 0){
        KPM_new.KPM[0] = -KPM_new.KPM[0];
        KPM_new.KPM[1] = -KPM_new.KPM[1];
        H(KPM_new, KPM_old, 2);
    }

    debug_message("Left hamiltonian::cheb\n");
}


