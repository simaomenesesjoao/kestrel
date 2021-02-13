void H(KPM_vector &KPM_new, KPM_vector &KPM_old, Array *hoppings, Eigen::Array<T, -1, -1> *anderson, unsigned Lx, unsigned Ly, unsigned f);
void cheb(KPM_vector &KPM_new, KPM_vector &KPM_old, Array *hoppings, Eigen::Array<T, -1, -1> *anderson, unsigned Lx, unsigned Ly, unsigned it_num);

T peierls(unsigned Lx, unsigned Ly, Array *hoppings, unsigned mult);


class hamiltonian{
    //private:


    public:
        unsigned Lx, Ly, Norb;
        unsigned N_hoppings;
        Array* hoppings;
        Array* anderson;
        double W;
        bool is_anderson_set;

        unsigned Lattice_Lx, Lattice_Ly;
        unsigned N_threads_x, N_threads_y;

        hamiltonian();
        ~hamiltonian();
        void set_geometry(unsigned, unsigned);
        void H(KPM_vector&, KPM_vector &, unsigned);
        double time_H();
        void cheb(KPM_vector&, KPM_vector &, unsigned);

        void set_anderson();
        void set_anderson_W(double);

        void set_regular_hoppings();
        void set_peierls(Eigen::Matrix<int, 2, 2>);



};

