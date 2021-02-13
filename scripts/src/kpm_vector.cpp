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
        KPM[i] = Eigen::Array<T, -1, -1>::Zero(LyG, LxG);
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
        KPM[orb](y,x) = 1;
    }
#pragma omp barrier
}


void KPM_vector::random_uniform(){
    // Initializes the KPM vector to a random 
    // vector of real (or complex) numbers
    // All orbitals get set to random numbers
    for(unsigned i = 0; i < Norb; i++){
        KPM[i].setRandom();
        KPM[i]/=sqrt(Npoints_lattice/3.0);
        if(nis_complex){
            KPM[i]/=sqrt(2.0);
        }
    }
}

void KPM_vector::fill_ghosts(){
    //unsigned cols, rows;
    unsigned t, t_up, t_down, t_left, t_right;
    unsigned t_ur, t_ul, t_dr, t_dl;
    unsigned cols, rows;


    cols = LxG;
    rows = LyG;

    // save to memory
    for(unsigned orb = 0; orb < Norb; orb++){

        // Ghosts on the sides
        // (in order) left, right, up, down
        sides[thread_id][orb][0] = KPM[orb].block(1,1,rows-2,1);        // Left
        sides[thread_id][orb][1] = KPM[orb].block(1,cols-2,rows-2,1);   // Right
        sides[thread_id][orb][2] = KPM[orb].block(1,1,1,cols-2);        // Up
        sides[thread_id][orb][3] = KPM[orb].block(rows-2,1,1,cols-2);   // Down

        // Ghosts in the corners
        // (in order) upper left (ul), upper right (ur), down left (dl), down right (dr)
        corners[thread_id][orb][0] = KPM[orb](1, 1);            // Upper left
        corners[thread_id][orb][1] = KPM[orb](1, cols-2);       // Upper right
        corners[thread_id][orb][2] = KPM[orb](rows-2, 1);       // Lower left
        corners[thread_id][orb][3] = KPM[orb](rows-2, cols-2);  // Lower right
    }

    // Make sure all threads have finished writing into memory
#pragma omp barrier

    unsigned tx, ty;
    for(unsigned orb = 0; orb < Norb; orb++){
        tx = thread_id%NTx;
        ty = thread_id/NTx;

        // Ghosts on the sides. Note that my notion of increasing y (up) is different
        // from the one printed in the shell (down)
        t_up    = tx + ((ty - 1 + NTy)%NTy)*NTx;
        t_down  = tx + ((ty + 1)%NTy)*NTx;
        t_left  = (tx - 1 + NTx)%NTx + ty*NTx;
        t_right = (tx + 1)%NTx + ty*NTx;

        KPM[orb].block(1,0,rows-2,1)        = sides[t_left][orb][1];    // left
        KPM[orb].block(1,cols-1,rows-2,1)   = sides[t_right][orb][0];   // right
        KPM[orb].block(0,1,1,cols-2)        = sides[t_up][orb][3];      // up
        KPM[orb].block(rows-1,1,1,cols-2)   = sides[t_down][orb][2];    // down

        t_ul = (tx - 1 + NTx)%NTx + ((ty - 1 + NTy)%NTy)*NTx;   // ul - upper left
        t_ur = (tx + 1)%NTx       + ((ty - 1 + NTy)%NTy)*NTx;   // ur - upper right
        t_dl = (tx - 1 + NTx)%NTx + ((ty + 1)%NTy)*NTx;         // dl - down left
        t_dr = (tx + 1)%NTx       + ((ty + 1)%NTy)*NTx;         // dr - down right

        // Ghosts in the corners
        KPM[orb](0, 0)            = corners[t_ul][orb][3];
        KPM[orb](rows-1, cols-1)  = corners[t_dr][orb][0];
        KPM[orb](0, cols-1)       = corners[t_ur][orb][2];
        KPM[orb](rows-1, 0)       = corners[t_dl][orb][1];

    }
#pragma omp barrier
}

void KPM_vector::empty_ghosts(){

    unsigned cols, rows;

    cols = LxG;
    rows = LyG;

    for(unsigned orb = 0; orb < Norb; orb++){
        // Ghosts on the sides
        KPM[orb].block(0,1,1,cols-2) = Eigen::Array<T, -1, -1>::Zero(1,cols-2);
        KPM[orb].block(rows-1,1,1,cols-2) = Eigen::Array<T, -1, -1>::Zero(1,cols-2);
        KPM[orb].block(1,0,rows-2,1) = Eigen::Array<T, -1, -1>::Zero(rows-2,1);
        KPM[orb].block(1,cols-1,rows-2,1) = Eigen::Array<T, -1, -1>::Zero(rows-2,1);

        // Ghosts in the corners
        KPM[orb](0, 0) = 0;
        KPM[orb](rows-1, cols-1) = 0;
        KPM[orb](0, cols-1) = 0;
        KPM[orb](rows-1, 0) = 0;
    }
}



void swap_pointers(KPM_vector *KPM1, KPM_vector *KPM0){
    unsigned Norb = KPM0->Norb;
    Eigen::Array<T, -1, -1> *KPM_temp;
    KPM_temp = new Eigen::Array<T, -1, -1>[Norb];

    // Swap the pointers for each of the orbitals between
    // each of the two KPM vectors
    for(unsigned orb = 0; orb < 2; orb++){
        KPM_temp[orb] = KPM1->KPM[orb];
        KPM1->KPM[orb] = KPM0->KPM[orb];
        KPM0->KPM[orb] = KPM_temp[orb];
    }

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
        KPM[orb] = Eigen::Array<T, -1, -1>::Zero(LyG, LxG);
        corners[thread_id][orb][0] = 0;
        corners[thread_id][orb][1] = 0;
        corners[thread_id][orb][2] = 0;
        corners[thread_id][orb][3] = 0;

        sides[thread_id][orb][0] = Eigen::Array<T, -1, -1>::Zero(Ly, 1);
        sides[thread_id][orb][1] = Eigen::Array<T, -1, -1>::Zero(Ly, 1);
        sides[thread_id][orb][2] = Eigen::Array<T, -1, -1>::Zero(1, Lx);
        sides[thread_id][orb][3] = Eigen::Array<T, -1, -1>::Zero(1, Lx);
    }
    }
        //if(thread_id == 2)
            //corners[thread_id][0][0] = 99;
#pragma omp barrier

    debug_message("Left KPM_vector set_geometry.\n");
}
