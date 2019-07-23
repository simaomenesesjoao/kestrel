
class chebinator{
    //using namespace std::chrono;
    private:
        unsigned Lx, Ly, Norb;
        KPM_vector KPM0, KPM1, KPM_initial;
        hamiltonian *H;
        std::chrono::high_resolution_clock::time_point start;
        std::chrono::high_resolution_clock::time_point current;
        unsigned current_moments, current_nrand, current_ndis;
        int seed;
        int mult;
    public:
        unsigned long int current_iter;
        unsigned long int max_iter;
        Eigen::Array<TR, -1, -1> mu;
        chebinator();
        void cheb_iteration(unsigned, unsigned, unsigned);
        void set_hamiltonian(hamiltonian*);
        void set_geometry(unsigned, unsigned, unsigned);
        void set_seed(int);
        void set_mult(int);

        double get_estimate(double, unsigned, unsigned, unsigned);
        void calc_finish();
};
