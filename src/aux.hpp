
double jackson(unsigned n, unsigned N);
Eigen::Array<TR, -1, -1> calc_dos(Eigen::Array<TR, -1, -1> mu, Eigen::Array<TR, -1, -1> energies, std::string mode);
void extended_euclidean(int a, int b, int *x, int *y, int *gcd);

class variables{

    public:
        double avg_time;
        unsigned iter;
        unsigned max_iter;
        bool has_initialized;
        bool has_finished;

        unsigned Lx, Ly;
        unsigned nrandom, ndisorder, nmoments;
        int mult, seed;
        double anderson_W;
        std::string status;

        variables(){
            has_initialized = false;
            has_finished = false;
            status="";
        }

        void update_status(){
            std::string init = "0";
            std::string fin = "0";
            if(has_initialized)
                init = "1";
            if(has_finished)
                fin = "1";
            //std::cout << "before string\n" << std::flush;
            status = init + " " + fin + " running " 
                + "Lx=" + std::to_string(Lx)         + " "
                + "Ly=" + std::to_string(Ly)         + " "
                + "NR=" + std::to_string(nrandom)    + " "
                + "ND=" + std::to_string(ndisorder)  + " "
                + "N="  + std::to_string(nmoments)   + " "
                + "B="  + std::to_string(mult)       + " "
                + "S="  + std::to_string(seed)       + " "
                + "W="  + std::to_string(anderson_W) + " "
                + "i="  + std::to_string(iter)       + " " 
                + "M="  + std::to_string(max_iter)   + " " 
                + "T="  + std::to_string(avg_time);
        };
        void update_status(std::string custom){
            status = custom;
        }

};
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

    bool need_print_to_file, need_print_to_cout, need_print;
    bool need_log_status;

    std::string filename_status;

};
void parse_input(int argc, char **argv, parameters*);
void load(std::string, KPM_vector *, KPM_vector *, Eigen::Array<TR, -1, -1> *);
void print_ham_info(parameters P);
void print_cheb_info(parameters P);
void print_compilation_info();
void print_magnetic_info(parameters, int, int, int);
void print_output_info(parameters);
void print_restart_info(parameters);
void print_log_info(parameters);
void print_dos(parameters P, unsigned N_lists, unsigned *nmoments_list, Eigen::Array<TR, -1, -1> *dos);
