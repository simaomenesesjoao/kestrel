#include <iostream>
#include <omp.h>
#include <Eigen/Dense>

#include "general.hpp"
#include "kpm_vector.hpp"

TR KPM_vector::operator*(const KPM_vector& other){
    // Make sure to empty ghosts first
    TR sum = 0;
    for(unsigned orb = 0; orb < Norb; orb++){
#if nis_complex == 1
        sum += (KPM[orb].conjugate()*other.KPM[orb]).sum().real();
#endif
#if nis_complex == 0
        sum += (KPM[orb]*other.KPM[orb]).sum();
#endif

    }
    return sum;
}

KPM_vector& KPM_vector::operator=(const KPM_vector& other){
    bool c1 = other.Norb==Norb;
    bool c2 = other.Lx==Lx;
    bool c3 = other.Ly==Ly;
    if(!(c1 && c2 && c3)){
        std::cout << "Cannot copy. Aborting.\n";
        exit(1);
    }
    for(unsigned orb = 0; orb < Norb; orb++){
        KPM[orb] = other.KPM[orb];
    }



    return *this;
};

KPM_vector::KPM_vector(){}

KPM_vector::KPM_vector(unsigned lx, unsigned ly, unsigned norb){
    debug_message("Entered KPM_vector constructor.\n");
    Lx = lx;
    Ly = ly;
    Norb = norb;
    LxG = Lx + 2;
    LyG = Ly + 2;
    Npoints = Lx*Ly*Norb;

    // Now that we know the number of orbitals, allocate
    // the memory needed for all the KPM vectors
    KPM = new Eigen::Array<T, -1, -1>[Norb];

    // Initialize each of the KPM vectors to the right dimensions
    for(unsigned i = 0; i < Norb; i++){
        KPM[i] = Eigen::Array<T, -1, -1>::Zero(LxG, LyG);
    }
    debug_message("Left KPM_vector constructor.\n");
}

KPM_vector::~KPM_vector(){
    debug_message("Entered KPM_vector destructor.\n");
    delete []KPM;
    debug_message("Left KPM_vector destructor.\n");
}

void KPM_vector::zero(){
    for(unsigned i = 0; i < Norb; i++){
        KPM[i].setZero();
    }
}

void KPM_vector::site(unsigned orb, unsigned x, unsigned y){
    // Initialize all the elements of KPM to zero
    // except for one single site which is one

    for(unsigned i = 0; i < Norb; i++){
        KPM[i].setZero();
    }

    if(thread_id == 0){
        KPM[orb](x,y) = 1;
    }
#pragma omp barrier
}


void KPM_vector::random_uniform(){
    debug_message("Entered random_uniform\n");
    // Initializes the KPM vector to a random 
    // vector of real (or complex) numbers
    // All orbitals get set to random numbers

    // Real case
    if(!nis_complex){
        for(unsigned i = 0; i < Norb; i++){
            KPM[i].setRandom();
            KPM[i]/=sqrt(Npoints_lattice/3.0);
        }
    }


    // complex case
    std::complex<double> im = std::complex<double>(0.0, 1.0);
    if(nis_complex){
        for(unsigned o = 0; o < Norb; o++){
            for(unsigned i = 0; i < LxG; i++){
                for(unsigned j = 0; j < LyG; j++){
                    KPM[o](i,j) = exp(rand()*2.0*M_PI/RAND_MAX*im)/sqrt(Npoints_lattice);
                }
            }
        }
    }
    debug_message("Left random_uniform\n");

}

void KPM_vector::fill_ghosts(){
    debug_message("Entered KPM_vectors::fill_ghosts\n");
    //unsigned cols, rows;
    unsigned t_up, t_down, t_left, t_right;
    unsigned t_ur, t_ul, t_lr, t_ll;

    // save to memory
    for(unsigned orb = 0; orb < Norb; orb++){


        //      Lattice
        //   # # # # # # # 
        //   # o o u o o #    ^y
        //   # l o o o r #    |
        //   # o o d o o #    | 
        //   # # # # # # #    ---> x

        // Ghosts on the sides
        // (in order) Down, Up, Left, Right
        sides[thread_id][orb][0] = KPM[orb].block(    1,     1, LxG-2,     1);   // Down
        sides[thread_id][orb][1] = KPM[orb].block(    1, LyG-2, LxG-2,     1);   // Up
        sides[thread_id][orb][2] = KPM[orb].block(    1,     1,     1, LyG-2);   // Left
        sides[thread_id][orb][3] = KPM[orb].block(LxG-2,     1,     1, LyG-2);   // Right

        // Ghosts in the corners
        // (in order) Lower left (ll), Upper left (ul), Lower right (lr), Upper right (ur)
        corners[thread_id][orb][0] = KPM[orb](    1,     1); // Lower left
        corners[thread_id][orb][1] = KPM[orb](    1, LyG-2); // Upper left
        corners[thread_id][orb][2] = KPM[orb](LxG-2,     1); // Lower right
        corners[thread_id][orb][3] = KPM[orb](LxG-2, LyG-2); // Upper right
    }

    // Make sure all threads have finished writing into memory
#pragma omp barrier

    unsigned tx, ty;
    for(unsigned orb = 0; orb < Norb; orb++){
        tx = thread_id%NTx;
        ty = thread_id/NTx;
        // tx + NTx * ty
        //
        // thread disposition in space
        // 3 4 5
        // 0 1 2 

        // Ghosts on the sides
        t_up    = tx                 + ((ty + 1)%NTy)*NTx;
        t_down  = tx                 + ((ty - 1 + NTy)%NTy)*NTx;
        t_left  = (tx - 1 + NTx)%NTx + ty*NTx;
        t_right = (tx + 1)%NTx       + ty*NTx;

        KPM[orb].block(    1,    0, LxG-2,     1) = sides[t_down ][orb][1]; // lower part of KPM connects to upper part (1) of thread below
        KPM[orb].block(    1,LyG-1, LxG-2,     1) = sides[t_up   ][orb][0]; // upper part of KPM connects to lower part (0) of thread above
        KPM[orb].block(    0,    1,     1, LyG-2) = sides[t_left ][orb][3]; // left part of KPM connects to right part (3) of thread to the left
        KPM[orb].block(LxG-1,    1,     1, LyG-2) = sides[t_right][orb][2]; // right part of KPM connects to left part (2) of thread to the right

        // Ghosts in the corners
        t_ll = (tx - 1 + NTx)%NTx + ((ty - 1 + NTy)%NTy)*NTx; // ll - lower left
        t_lr = (tx + 1)%NTx       + ((ty - 1 + NTy)%NTy)*NTx; // lr - lower right
        t_ul = (tx - 1 + NTx)%NTx + ((ty + 1)%NTy)*NTx;       // ul - upper left
        t_ur = (tx + 1)%NTx       + ((ty + 1)%NTy)*NTx;       // ur - upper right

        KPM[orb](    0,     0) = corners[t_ll][orb][3]; // lower left part of KPM connects to upper right (3) part of thread to lower left
        KPM[orb](    0, LyG-1) = corners[t_ul][orb][2]; // upper left part of KPM connects to lower right (2) part of thread to upper left
        KPM[orb](LxG-1,     0) = corners[t_lr][orb][1]; // lower right part of KPM connects to upper left (1) part of thread to lower right
        KPM[orb](LxG-1, LyG-1) = corners[t_ur][orb][0]; // upper right part of KPM connects to lower left (0) part of thread to upper right

    }
#pragma omp barrier
    debug_message("Left KPM_vectors::fill_ghosts\n");
}

void KPM_vector::empty_ghosts(){
    debug_message("Entered KPM_vectors::empty_ghosts\n");

    for(unsigned orb = 0; orb < Norb; orb++){
        // Ghosts on the sides
        KPM[orb].block(    0,    1,    1,LyG-2) = Eigen::Array<T, -1, -1>::Zero(    1,LyG-2);
        KPM[orb].block(LxG-1,    1,    1,LyG-2) = Eigen::Array<T, -1, -1>::Zero(    1,LyG-2);
        KPM[orb].block(    1,    0,LxG-2,    1) = Eigen::Array<T, -1, -1>::Zero(LxG-2,    1);
        KPM[orb].block(    1,LyG-1,LxG-2,    1) = Eigen::Array<T, -1, -1>::Zero(LxG-2,    1);

        // Ghosts in the corners
        KPM[orb](    0,     0) = 0;
        KPM[orb](LxG-1, LyG-1) = 0;
        KPM[orb](    0, LyG-1) = 0;
        KPM[orb](LxG-1,     0) = 0;
    }
    debug_message("Left KPM_vectors::empty_ghosts\n");
}



void swap_pointers(KPM_vector *KPM1, KPM_vector *KPM0){
    debug_message("Entered KPM_vector::swap_pointers.\n");
    unsigned Norb = KPM0->Norb;
    Eigen::Array<T, -1, -1> *KPM_temp;
    KPM_temp = new Eigen::Array<T, -1, -1>[Norb];

    // Swap the pointers for each of the orbitals between
    // each of the two KPM vectors
    for(unsigned orb = 0; orb < Norb; orb++){
        KPM_temp[orb] = KPM1->KPM[orb];
        KPM1->KPM[orb] = KPM0->KPM[orb];
        KPM0->KPM[orb] = KPM_temp[orb];
    }

    debug_message("Left KPM_vector::swap_pointers.\n");
    delete []KPM_temp;
}


void KPM_vector::set_geometry(unsigned lx, unsigned ly, unsigned norb){
    debug_message("Entered KPM_vector set_geometry.\n");
    Lx = lx;
    Ly = ly;
    Norb = norb;
    LxG = Lx + 2;
    LyG = Ly + 2;
    Npoints = Lx*Ly*Norb;
    //std::cout << "thraedid: " << thread_id << "\n" << std::flush;

    // Now that we know the number of orbitals, allocate
    // the memory needed for all the KPM vectors
    KPM = new Eigen::Array<T, -1, -1>[Norb];

#pragma omp barrier
#pragma omp critical
    {
    // Initialize each of the KPM vectors to the right dimensions
    for(unsigned orb = 0; orb < Norb; orb++){
        KPM[orb] = Eigen::Array<T, -1, -1>::Zero(LxG, LyG);
        corners[thread_id][orb][0] = 0;
        corners[thread_id][orb][1] = 0;
        corners[thread_id][orb][2] = 0;
        corners[thread_id][orb][3] = 0;

        sides[thread_id][orb][0] = Eigen::Array<T, -1, -1>::Zero(Lx, 1);
        sides[thread_id][orb][1] = Eigen::Array<T, -1, -1>::Zero(Lx, 1);
        sides[thread_id][orb][2] = Eigen::Array<T, -1, -1>::Zero(1, Ly);
        sides[thread_id][orb][3] = Eigen::Array<T, -1, -1>::Zero(1, Ly);
    }
    }
#pragma omp barrier

    debug_message("Left KPM_vector set_geometry.\n");
}
