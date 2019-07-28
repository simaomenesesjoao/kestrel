#include <iostream>
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
        std::cout << "Cannot copy.\n";
    
    }
    for(unsigned orb = 0; orb < Norb; orb++){
        KPM[orb] = other.KPM[orb];
    }



    return *this;
};

KPM_vector::KPM_vector(){}

KPM_vector::KPM_vector(unsigned lx, unsigned ly, unsigned norb){
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

void KPM_vector::random_uniform(){
    // Initializes the KPM vector to a random 
    // vector of real (or complex) numbers
    // All orbitals get set to random numbers
    for(unsigned i = 0; i < Norb; i++){
        KPM[i].setRandom();
        KPM[i]/=sqrt(Npoints/3.0);
        if(nis_complex){
            KPM[i]/=sqrt(2.0);
        }
    }
}

void KPM_vector::fill_ghosts(){
    unsigned cols, rows;

    cols = LxG;
    rows = LyG;

    for(unsigned orb = 0; orb < Norb; orb++){
        // Ghosts on the sides
        KPM[orb].block(0,1,1,cols-2) = KPM[orb].block(rows-2,1,1,cols-2);
        KPM[orb].block(rows-1,1,1,cols-2) = KPM[orb].block(1,1,1,cols-2);
        KPM[orb].block(1,0,rows-2,1) = KPM[orb].block(1,cols-2,rows-2,1);
        KPM[orb].block(1,cols-1,rows-2,1) = KPM[orb].block(1,1,rows-2,1);

        // Ghosts in the corners
        KPM[orb](0, 0) = KPM[orb](rows-2, cols-2);
        KPM[orb](rows-1, cols-1) = KPM[orb](1, 1);
        KPM[orb](0, cols-1) = KPM[orb](rows-2, 1);
        KPM[orb](rows-1, 0) = KPM[orb](1, cols-2);

    }

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

}
