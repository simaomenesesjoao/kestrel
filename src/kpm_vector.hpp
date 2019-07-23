
class KPM_vector{
    public:
        unsigned Lx, Ly;        // System dimensions
        unsigned LxG, LyG;      // System dimensions with ghosts
        unsigned Norb;          // Number of orbitals
        unsigned Npoints;
        Eigen::Array<T, -1, -1> *KPM;   // KPM vector

        KPM_vector& operator=(const KPM_vector& other);

        KPM_vector(unsigned, unsigned, unsigned);
        void set_geometry(unsigned, unsigned, unsigned);
        KPM_vector();
        ~KPM_vector();
        void random_uniform();
        void zero();
        void fill_ghosts();
        void empty_ghosts();
        TR operator*(const KPM_vector& other);
};


void swap_pointers(KPM_vector *KPM1, KPM_vector *KPM0);
