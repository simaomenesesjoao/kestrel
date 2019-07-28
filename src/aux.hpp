
double jackson(unsigned n, unsigned N);
Eigen::Array<TR, -1, -1> calc_dos(Eigen::Array<TR, -1, -1> mu, Eigen::Array<TR, -1, -1> energies, std::string mode);
void extended_euclidean(int a, int b, int *x, int *y, int *gcd);

struct parameters{
    unsigned Lx, Ly;
    unsigned nrandom, ndisorder, nmoments;
    int mult, seed;
    double anderson_W;
    Eigen::Array<TR, -1, -1> energies;
    unsigned NEnergies;
    bool found_NEnergies;
    bool found_energies;

    bool output_energies;
    bool need_write;

    std::string saveto, readfrom;
    bool need_moremoments;
    unsigned moremoments;
    bool need_read;
    
    std::string filename_read;
    std::string filename_write;

};
void parse_input(int argc, char **argv, parameters*);
void load(std::string, KPM_vector *, KPM_vector *, Eigen::Array<TR, -1, -1> *);
void print_ham_info(parameters P);
void print_cheb_info(parameters P);
void print_compilation_info();
void print_magnetic_info(parameters, int, int, int);
void print_output_info(parameters);
void print_restart_info(parameters);
