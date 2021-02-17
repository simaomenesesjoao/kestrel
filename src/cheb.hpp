
class chebinator{
    //using namespace std::chrono;
    private:
        unsigned Lx, Ly, Norb;
        hamiltonian *H;
        std::chrono::high_resolution_clock::time_point start;
        std::chrono::high_resolution_clock::time_point current;
        unsigned current_moments, current_nrand, current_ndis;
        int seed;
        int mult;
        bool initialized;

    public:
        T ***corners;
        Eigen::Array<T, -1, -1> ***sides;
        unsigned thread_id, n_threads, n_threads_x, n_threads_y;
        KPM_vector KPM0, KPM1, KPM_initial;
        //variables *vars;
        unsigned long int current_iter;
        unsigned long int max_iter;
        Eigen::Array<TR, -1, -1> mu;
        chebinator();
        ~chebinator();

        parameters P;

        bool need_log_status;
        std::string filename_status;

        void cheb_iteration(unsigned, unsigned, unsigned);
        void cheb_iteration_restart(unsigned, KPM_vector&, KPM_vector&, Eigen::Array<TR, -1, -1>);
        void set_hamiltonian(hamiltonian*);
        void set_geometry(unsigned, unsigned, unsigned);
        void set_seed(int);
        void set_mult(int);

        double get_estimate(double, unsigned, unsigned, unsigned);
        void save(std::string);
};
