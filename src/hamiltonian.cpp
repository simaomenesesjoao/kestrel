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

void hamiltonian::set_origin_to_from(Eigen::Array<int, -1,2> origin_temp, Eigen::Array<unsigned, 1,-1> to_temp, Eigen::Array<unsigned, 1, -1> from_temp){
    debug_message("Entered hamiltonian::set_origin_to_from\n");
    origin = origin_temp*1;
    to_matrix = to_temp*1;
    from_matrix = from_temp*1;
    debug_message("Left hamiltonian::set_origin_to_from\n");
}

void hamiltonian::set_orbpos(Eigen::Array<double, -1, -1> orbpos_temp){
    debug_message("Entered hamiltonian::set_orbpos\n");

    orb_pos = orbpos_temp*1.0;

    debug_message("Left hamiltonian::set_orbpos\n");
}

void hamiltonian::set_primitive2(Eigen::Array<double, 2,1> a1_temp, Eigen::Array<double, 2,1> a2_temp){
    debug_message("Entered hamiltonian::set_primitive2\n");
    
    // Lattice vectors
    a1 = a1_temp*1.0;
    a2 = a2_temp*1.0;
    Vcell = abs(a1[0]*a2[1] - a1[1]*a2[0]);

    debug_message("Left hamiltonian::set_primitive2\n");
}

void hamiltonian::set_peierls(Eigen::Array<int, 2, 2> gauge_matrix){
    debug_message("Entered hamiltonian::set_peierls\n");
    Eigen::Array<TR, -1, -1> A(2,2);
    A(0,0) = gauge_matrix(0,0)*1.0/Lattice_Lx;
    A(1,1) = gauge_matrix(1,1)*1.0/Lattice_Ly;
    A(0,1) = gauge_matrix(0,1)*1.0/Lattice_Lx;
    A(1,0) = gauge_matrix(1,0)*1.0/Lattice_Ly;

    double flux = A(0,1) - A(1,0); //magnetic flux per unit cell, in units of the flux quantum


    // neighbour distance for each hopping
    std::complex<TR> im = std::complex<TR>(0, 1.0);



    
    /*      A     A              A       A
     *       \   /            t5  \     /  t4
     *        \ /                  \   /
     *         B                    >B<              B   
     *         |                     ^               |  t0 
     *         |                     |  t3           v
     *         A                     A             _ A _
     *        / \                             t1   /| |\  t2        
     *       /   \                                /     \
     *      B     B                              B       B 
     */

    //



    unsigned id;
    unsigned tx, ty;
    unsigned offset_x, offset_y;
#pragma omp barrier
#pragma omp critical
    {
    id = omp_get_thread_num(); // id = ty*N_threads_x + tx
    tx = id%N_threads_x;
    ty = id/N_threads_x;
    unsigned size_x = Lattice_Lx/double(N_threads_x) + 0.9;
    unsigned size_y = Lattice_Ly/double(N_threads_y) + 0.9;

    offset_x = tx*size_x;
    offset_y = ty*size_y;
    }

//#pragma omp critical
#pragma omp barrier
    //{

    int di, dj;
    double phase1;
    double phase2;
    double phase;
    for(unsigned hop = 0; hop < N_hoppings; hop++){
        for(unsigned i = 0; i < Lx; i++){
            for(unsigned j = 0; j < Ly; j++){

                di = origin(hop,0);
                dj = origin(hop,1);

                double area1, area2, area3, area4, area;
                Eigen::Matrix<TR, 2,1> da, db, de;
                de = di*a1 + dj*a2;

                unsigned to, from;
                to   = to_matrix[hop];
                from = from_matrix[hop];

                da = orb_pos(to  ,0)*a1 + orb_pos(to  ,1)*a2;
                db = orb_pos(from,0)*a1 + orb_pos(from,1)*a2;
                area1 = (db[0] + da[0])*de[1];
                area2 = -de[0]*(db[1] + da[1]);
                area3 = da[0]*db[1];
                area4 = - db[1]*da[0];
                area = 0.5*(area1 + area2 + area3 + area4);

                phase1  = (1.0*(i + offset_x) + di/2.0)*dj*A(0,1) + (1.0*(j + offset_y) + dj/2.0)*di*A(1,0);
                phase2 = flux*area/Vcell;
                phase = phase1 + phase2;

                //std::cout << "hop" << hop << " " << i << "," << j << " " << "di,dj:" << di << "," << dj<<" " << phase1 << " " << phase2 << " area1,2,3,4:" << area1 << " " << area2 << " " << area3 << " " << area4 << "  da,db" << da[0] << "," << da[1] << " " << db[0] << "," << db[1] << "\n";
                //
                //std::cout << "thread:" << id << " hop" << hop << " " << i + offset_x << "," << j + offset_y << " " << "di,dj:" << di << "," << dj<<" phases_reg,inv:" << phase1 << " " << phase2 << " Area:" << area << "  da,db" << da[0] << "," << da[1] << " " << db[0] << "," << db[1] << "\n";


                hoppings[hop](i, j) *= exp(2.0*M_PI*im*phase);
            }
        }
    }
    //}
//#pragma omp barrier

    debug_message("Left hamiltonian::set_peierls\n");
}

void hamiltonian::set_geometry(unsigned lx, unsigned ly){
    Lx = lx;
    Ly = ly;
}


hamiltonian::hamiltonian(unsigned Norb_temp, unsigned N_hoppings_temp){
    std::cout << "N_hoppings: " << N_hoppings_temp << "\n";
    Norb = Norb_temp;
    N_hoppings = N_hoppings_temp;
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
        hoppings[i] = Array::Zero(Lx, Ly) + 1.0/SCALE;
    }
}

void hamiltonian::set_anderson_W(double w_anderson){
    W = w_anderson;
}

void hamiltonian::set_anderson(){
    debug_message("Entered hamiltonian::set_anderson\n");
        for(unsigned j = 0; j < Norb; j++){
            anderson[j] = Eigen::Array<TR, -1, -1>::Random(Lx, Ly)*W/2.0; // random numbers between -W/2 and W/2
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
    debug_message("Left hamiltonian::set_vacancies\n");
}

void hamiltonian::H(KPM_vector &KPM_new, KPM_vector &KPM_old, unsigned f){
    debug_message("Entered hamiltonian::H\n");
    // Writes the result of H*KPM_old into KPM_new by adding. It does not replace the values
    // This is intentional but may cause some problems if KPM_new is not initialized to zero

//#pragma omp barrier
//#pragma omp critical
    //{

    // Anderson disorder for each orbital
    debug_message("Adding Anderson disorder\n");
    for(unsigned i = 0; i < Norb; i++){
        KPM_new.KPM[i].block(1, 1, Lx, Ly) += KPM_old.KPM[i].block(1, 1, Lx, Ly)*anderson[i]*f;
    }



    int d1, d2;
    unsigned to, from;

    for(unsigned hop = 0; hop < N_hoppings; hop++){
        to = to_matrix[hop];
        from = from_matrix[hop];
        d1 = origin(hop,0);
        d2 = origin(hop,1);

        //std::cout << "hop:" << hop << " to,from:" << to << "," << from << " d1,d2:" << d1 << "," << d2 << "\n" <<std::flush;

        KPM_new.KPM[to].block(1,1,Lx,Ly) += KPM_old.KPM[from].block(1+d1,1+d2,Lx,Ly)*hoppings[hop]*f;
    }

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
        for(unsigned orb = 0; orb < Norb; orb++){
            KPM_new.KPM[orb] = -KPM_new.KPM[orb];
        }
        H(KPM_new, KPM_old, 2);
    }

    debug_message("Left hamiltonian::cheb\n");
}


