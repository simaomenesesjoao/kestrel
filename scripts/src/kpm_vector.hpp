
class KPM_vector{
    public:
        unsigned Lx, Ly;            // System dimensions
        unsigned LxG, LyG;          // System dimensions with ghosts
        unsigned Norb;              // Number of orbitals
        unsigned Npoints;           // Number of points in this subdomain

        unsigned Npoints_lattice;   // Number of points in the whole lattice
        unsigned NTx, NTy, thread_id;
        Eigen::Array<T, -1, -1> *KPM;   // KPM vector
        T ***corners;
        Eigen::Array<T, -1, -1> ***sides;

        KPM_vector& operator=(const KPM_vector& other);

        KPM_vector(unsigned, unsigned, unsigned);
        void set_geometry(unsigned, unsigned, unsigned);
        KPM_vector();
        ~KPM_vector();
        void random_uniform();
        void site(unsigned, unsigned, unsigned);
        void zero();
        void fill_ghosts();
        void empty_ghosts();
        TR operator*(const KPM_vector& other);
};


void swap_pointers(KPM_vector *KPM1, KPM_vector *KPM0);
