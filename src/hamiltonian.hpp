void H(KPM_vector &KPM_new, KPM_vector &KPM_old, Array *hoppings, Eigen::Array<T, -1, -1> *anderson, unsigned Lx, unsigned Ly, unsigned f);
void cheb(KPM_vector &KPM_new, KPM_vector &KPM_old, Array *hoppings, Eigen::Array<T, -1, -1> *anderson, unsigned Lx, unsigned Ly, unsigned it_num);

T peierls(unsigned Lx, unsigned Ly, Array *hoppings, unsigned mult);


class hamiltonian{
    private:
        Eigen::Array<double, 2, 1> a1, a2;
        Eigen::Array<unsigned, 1, -1> from_matrix, to_matrix;
        Eigen::Array<int, -1, -1> origin;
        Eigen::Array<double, -1, -1> orb_pos; // each row is the position of an orbital in units of the primitive vectors


    public:
        unsigned Lx, Ly, Norb;
        unsigned N_hoppings;
        double scale;
        double Vcell; //Volume (area, in 2D) of the unit cell)
        Array* hoppings;
        Array* anderson;
        Eigen::Array<unsigned, -1, -1> vacanciesA, vacanciesB;
        unsigned NvacA, NvacB;
        double W;
        bool is_anderson_set;

        unsigned Lattice_Lx, Lattice_Ly;
        unsigned N_threads_x, N_threads_y;

        hamiltonian(unsigned, unsigned);
        ~hamiltonian();
        void set_geometry(unsigned, unsigned);
        void H(KPM_vector&, KPM_vector &, unsigned);
        void cheb(KPM_vector&, KPM_vector &, unsigned);

        void set_anderson();
        void set_anderson_W(double);
        void set_vacancies(Eigen::Array<unsigned, -1, -1>, Eigen::Array<unsigned, -1, -1>);

        void set_regular_hoppings();
        void set_peierls(Eigen::Array<int, 2, 2>);
        void set_primitive2(Eigen::Array<double, 2,1>, Eigen::Array<double, 2,1>);
        void set_origin_to_from(Eigen::Array<int, -1,2>, Eigen::Array<unsigned, 1,-1>, Eigen::Array<unsigned, 1, -1>);
        void set_orbpos(Eigen::Array<double, -1, -1>);



};

