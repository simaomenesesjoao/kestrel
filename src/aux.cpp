#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include <H5Cpp.h>
#include <H5Group.h>

#include "general.hpp"
#include "kpm_vector.hpp"
#include "aux.hpp"

#include "ComplexTraits.hpp"
#include "myHDF5.hpp"
#include "messages.hpp"

void print_dos(parameters P, unsigned N_lists, unsigned *nmoments_list, Eigen::Array<TR, -1, -1> *dos){
    unsigned N_energies = dos[0].size(); // They all have the same number of points

    if(P.need_print_to_cout){
        std::cout << "________Density of states________\n";
        for(unsigned i = 0; i < N_energies; i++){
            std::string metadata;
            metadata  = "Lx:" + std::to_string(P.Lx);
            metadata += " Ly:" + std::to_string(P.Ly);
            metadata += " Bmult:" + std::to_string(P.mult);
            metadata += " nrandom:" + std::to_string(P.nrandom);
            metadata += " ndisorder:" + std::to_string(P.ndisorder);
            metadata += " W:" + std::to_string(P.anderson_W);
            metadata += " C:" + std::to_string(P.conc);
            metadata += " seed:" + std::to_string(P.seed);

            for(unsigned j = 0; j < N_lists; j++){
                std::cout << metadata << " N:" << nmoments_list[j] << " en:" << P.energies(i) << " dos:" << dos[j](i)/SCALE << "\n";
            }
        }
    }
    
    if(P.need_print_to_file){

         //Save to file
        std::ofstream file;

        for(unsigned j = 0; j < N_lists; j++){
            std::string name;
            name = "dos_N" + std::to_string(nmoments_list[j]) + "_W" + std::to_string(P.anderson_W) + 
                "_C" + std::to_string(P.conc) + "_B" + std::to_string(P.mult) + ".dat";
            file.open(name);
            for(unsigned i = 0; i < N_energies; i++)
                file << P.energies(i) << " " << dos[j](i) << "\n";
            file.close();
        }

    }

    if(!P.need_print) std::cout << "Nothing to print.\n";
}

Eigen::Array<TR, -1, -1> calc_dos(Eigen::Array<TR, -1, -1> mu, Eigen::Array<TR, -1, -1> energies, std::string mode){
    debug_message("Entered calc_dos\n");
    //std::cout << energies << "\n";

    unsigned N_energies = energies.size();
    unsigned moments = mu.size();
    if(moments <= 0){
        std::cout << "Invalid number of chebyshev moments (" << moments << ")"
            " inside the calculation of the density of states. Aborting.\n";
        exit(1);
    }
    if(N_energies <=0){
        std::cout << "Invalid number of energies (" << N_energies << ") inside the "
            "calculation of the density of states. Aborting.\n";
        exit(1);
    }
    std::complex<TR> im(0.0, 1.0);
    Eigen::Array<TR, -1, -1> dos(N_energies, 1);


    if(mode == "jackson"){
        verbose2("Entered mode: jackson\n");
        if(moments > 3){
            debug_message("Number of moments inside DoS is > 3.\n");
            Eigen::Array<TR, -1, -1> temp[2], sq_en(N_energies, 1);
            sq_en = Eigen::rsqrt(1.0 - energies*energies)*2.0/M_PI;
            temp[0] = Eigen::Array<TR, -1, -1>::Zero(N_energies, 1) + 1;
            temp[1] = energies;
            //temp[0] = temp[1]*2*energies - temp[0]

            dos  = temp[0]*mu(0)/2.0*jackson(0, moments);
            dos += temp[1]*mu(1)*jackson(1, moments);
            for(unsigned i = 2; i < moments; i++){
                temp[i%2] = 2*energies*temp[(i+1)%2] - temp[i%2];
                dos += temp[i%2]*mu(i)*jackson(i, moments);
            }
            dos = dos*sq_en;
        } else {
            debug_message("Number of moments inside DoS is <= 3 ("<<moments<<").\n");
            Eigen::Array<TR, -1, -1> temp[3], sq_en(N_energies, 1);
            sq_en = Eigen::rsqrt(1.0 - energies*energies)*2.0/M_PI;

            temp[0] = Eigen::Array<TR, -1, -1>::Zero(N_energies, 1) + 1;
            temp[1] = energies;
            temp[2] = 2*energies*energies - 1;


            dos.setZero();
            for(unsigned i = 0; i < moments; i++){
                dos += temp[i]*mu(i)*jackson(i, moments);
            }
            dos = dos*sq_en;
        }
    } else {
        double parameter;
        std::string mode1;
        try{
            mode1 = mode.substr(0, mode.find(" "));
            parameter = std::stod(mode.substr(mode.find(" ") + 1, mode.length()));
        }
        catch(int e){
            std::cout << "Unkown mode for the density of states. Aborting.\n";
            exit(1);
        }

        if(mode1 == "green"){
            verbose2("Entered mode: " << mode1 << " with parameter " << parameter << "\n");
            TR lambda = parameter;
            Eigen::Array<std::complex<TR>, -1, -1> green(N_energies, 1), complex_sqrt(N_energies, 1);
            complex_sqrt = -Eigen::rsqrt(1.0 - (energies + im*lambda)*(energies + im*lambda))*(-2.0*im)/M_PI;
            green = Eigen::exp(-Eigen::acos(energies + im*lambda)*im);

            dos = (complex_sqrt*mu(0)/2.0).imag();
            for(unsigned i = 1; i < moments; i++){
                complex_sqrt = complex_sqrt*green;
                dos += complex_sqrt.imag()*mu(i);
            }
        } else {
            std::cout << "Unkown mode for the density of states. Aborting.\n";
            exit(1);
        }
    }

    debug_message("Left DoS");
    return dos/SCALE;

}

template <typename T>
bool is_numeric (std::string & str){
    //auto result = T();
    //auto i = std::istringstream(str);

    //i >> result;
    //i >> std::ws;

    //return !i.fail() && i.eof();
    return true;
}

void parse_input(int argc, char **argv, parameters *P){ 
    debug_message("Entered parse_input");

    // Convert everything to a string to make it easier to process
    std::string *inputs;
    inputs = new std::string[argc-1];
    for(int i = 0; i < argc-1; i++){
        inputs[i] = std::string(argv[i+1]);
    }

    std::string *p_threads;
    std::string *p_geometry, *p_mult, *p_anderson, *p_conc;      // Hamiltonian
    std::string *p_nrandom, *p_nmoments, *p_numdis, *p_seed;    // Cheb
    std::string *p_NEnergies, *p_energies;              // output
    std::string *p_readfrom, *p_saveto;                 // restart
    std::string *p_help;                                // help
    std::string *p_moremoments;
    std::string *p_print_to_cout, *p_print_to_file;
    std::string *p_log_status;

    // find the pointers to each of the quantities
    verbose2("finding pointers");
    p_threads     = std::find(inputs, inputs+argc-1, "--thread_division");

    p_geometry    = std::find(inputs, inputs+argc-1, "--geometry");
    p_mult        = std::find(inputs, inputs+argc-1, "--mult");
    p_anderson    = std::find(inputs, inputs+argc-1, "--anderson");
    p_conc        = std::find(inputs, inputs+argc-1, "--vacancies");

    p_nrandom     = std::find(inputs, inputs+argc-1, "--nrandom");
    p_numdis      = std::find(inputs, inputs+argc-1, "--ndisorder");
    p_nmoments    = std::find(inputs, inputs+argc-1, "--nmoments");
    p_seed        = std::find(inputs, inputs+argc-1, "--seed");

    p_NEnergies   = std::find(inputs, inputs+argc-1, "--NEnergies");
    p_energies    = std::find(inputs, inputs+argc-1, "--energies");

    p_readfrom    = std::find(inputs, inputs+argc-1, "--readfrom");
    p_moremoments = std::find(inputs, inputs+argc-1, "--moremoments");
    p_saveto      = std::find(inputs, inputs+argc-1, "--save_restart_to");
    p_help        = std::find(inputs, inputs+argc-1, "--help");

    p_print_to_cout = std::find(inputs, inputs+argc-1, "--print_to_cout");
    p_print_to_file = std::find(inputs, inputs+argc-1, "--print_to_file");

    p_log_status    = std::find(inputs, inputs+argc-1, "--log_status");





    // Check which of the possible quantities has been defined 
    // through the command line
    verbose2("checking existence of command line instructions");
    bool has_nmoments, has_ndisorder, has_nrandom, has_seed;
    bool has_geometry, has_anderson, has_mult, has_conc;
    bool has_thread;
    bool need_read, need_write, need_help;
    bool found_NEnergies, found_energies;
    bool need_moremoments;
    bool need_print_to_file, need_print_to_cout, need_print;
    bool need_log_status;

    std::string *p_final = inputs + argc - 1;

    has_nmoments     = p_nmoments    != p_final;
    has_ndisorder    = p_numdis      != p_final;
    has_nrandom      = p_nrandom     != p_final;
    has_seed         = p_seed        != p_final;

    found_NEnergies  = p_NEnergies   != p_final;
    found_energies   = p_energies    != p_final;

    has_geometry     = p_geometry    != p_final;
    has_thread       = p_threads     != p_final;
    has_mult         = p_mult        != p_final;
    has_anderson     = p_anderson    != p_final;
    has_conc         = p_conc        != p_final;
 
    need_read        = p_readfrom    != p_final;
    need_write       = p_saveto      != p_final;
    need_help        = p_help        != p_final;
    need_moremoments = p_moremoments != p_final;

    need_print_to_file = p_print_to_file != p_final;
    need_print_to_cout = p_print_to_cout != p_final;
    need_print         = need_print_to_file || need_print_to_cout;


    P->need_print_to_file = need_print_to_file;
    P->need_print_to_cout = need_print_to_cout;
    P->need_print = need_print;

    need_log_status    = p_log_status != p_final;
    P->need_log_status = need_log_status;
    P->filename_status = "status";


    // Check if the input from the command line is complete
    verbose2("Checking if input is complete\n");
    bool has_any_cheb, has_any_ham, has_any_cl_input;
    bool has_all_cheb, has_all_ham, has_all_cl_input;
    bool output_energies;

    has_any_cheb = has_nmoments || has_ndisorder || has_nrandom || has_seed;
    has_any_ham  = has_geometry || has_mult || has_anderson || has_thread || has_conc;
    has_any_cl_input = has_any_cheb || has_any_ham;

    has_all_cheb = has_nmoments && has_ndisorder && has_nrandom && has_seed;
    has_all_ham  = has_geometry && has_mult && has_anderson && has_conc;
    has_all_cl_input = has_all_ham && has_all_cheb;

    output_energies = found_NEnergies || found_energies;

    if(need_help){
        std::cout << help_message;
        exit(0);
    }

    if(need_print && !output_energies){
        std::cout << "Printing was requested but no energies were specified. "
            "Please specify --energies or --NEnergies when using --print_to_cout "
            "or --print_to_file. Aborting.\n";
        exit(1);
    }
    if(!need_print && output_energies){
        std::cout << "Printing was not requested but energies were specified. "
            "Please specify --print_to_file or --print_to_cout when using "
            "--energies or --NEnergies. Aborting.\n";
        exit(1);
    }

    // Cannot ask to use the quantities specified by the command 
    // line (calculate new) if you're also fetching these quantities
    // from the input file
    if(has_any_cl_input && need_read){
        std::cout << "The input data must be inserted either through a "
            "configuration file through --readfrom or from the command line "
            "using the needed parameters. Both cannot be specified "
            "simultaneously. Aborting.\n";
        exit(1);
    }

    // If no input data is specified
    if(!has_any_cl_input && !need_read){
        std::cout << "The input data must be inserted either through a "
            "configuration file through --readfrom or from the command line "
            "using the needed parameters. At least one must be specified. "
            "Aborting.\n";
        exit(1);
    }

    // Not output at all would be given. This would be a useless 
    // waste of computational resources and should give an error
    if(!output_energies && !need_write){
        std::cout << "No output is being requested. Please specify either "
            "--NEnergies, --energies or --save_restart_to. Aborting.\n";
        exit(1);
    }

    if(found_energies && found_NEnergies){
        std::cout << "Both the number of energies (NEnergies) and the energies have ";
        std::cout << "been specified. Only one of them may be specified. Aborting.\n";
        exit(1);
    }

    // Not all the needed information has been provided
    if(has_any_cl_input && !has_all_cl_input && !need_read){
        std::cout << "The input is being specified through the command line. "
            "Please provide the complete input. Aborting.\n";
        exit(1);
    }
    
    if(!need_read && need_moremoments){
        std::cout << "The option --moremoments should only be specified "
            "when reading from a configuration file with information about "
            "a previously saved system state (KPM0 and KPM1). Aborting.\n";
        exit(1);
    }


    // Find the values of each of the quantities
    verbose2("Finding the values of each of the quantities.\n");
    if(!need_read){
        try{
            std::string s_nmoments, s_nrandom, s_numdis, s_mult;
            std::string s_Lx, s_Ly, s_seed, s_anderson, s_conc;
            std::string s_nx, s_ny;

            s_nmoments  = *(p_nmoments+1);
            s_nrandom   = *(p_nrandom+1 );
            s_numdis    = *(p_numdis+1  );
            s_mult      = *(p_mult+1    );                
            s_Lx        = *(p_geometry+1);
            s_Ly        = *(p_geometry+2);
            s_nx        = *(p_threads+1 );
            s_ny        = *(p_threads+2 );
            s_seed      = *(p_seed+1);
            s_anderson  = *(p_anderson+1);
            s_conc      = *(p_conc+1);

            //std::cout << "SNX: " << s_nx << " " << s_ny << "\n" << std::flush;
            P->nmoments       = atoi((s_nmoments     ).c_str());
            P->nrandom        = atoi((s_nrandom      ).c_str());
            P->ndisorder      = atoi((s_numdis       ).c_str());
            P->mult           = atoi((s_mult         ).c_str());
            P->Lx             = atoi((s_Lx           ).c_str());
            P->Ly             = atoi((s_Ly           ).c_str());
            P->nx             = atoi((s_nx           ).c_str());
            P->ny             = atoi((s_ny           ).c_str());
            P->seed           = atoi((s_seed         ).c_str());
            P->anderson_W     = atof((s_anderson     ).c_str());
            P->conc           = atof((s_conc         ).c_str());

            //std::cout << "ints: " << P->nx << " " << P->ny << "\n" << std::flush;
            bool b_nmoments, b_nrandom, b_numdis, b_mult;
            bool b_Lx, b_Ly, b_seed, b_anderson, b_conc;
            bool b_nx, b_ny;

            b_nmoments  = is_numeric<int>(s_nmoments); 
            b_nrandom   = is_numeric<int>(s_nrandom);   
            b_numdis    = is_numeric<int>(s_numdis);  
            b_mult      = is_numeric<int>(s_mult);  
            b_Lx        = is_numeric<int>(s_Lx);  
            b_Ly        = is_numeric<int>(s_Ly);  
            b_nx        = is_numeric<int>(s_nx);  
            b_ny        = is_numeric<int>(s_ny);  
            b_seed      = is_numeric<int>(s_seed);  
            b_anderson  = is_numeric<double>(s_anderson);  
            b_conc      = is_numeric<double>(s_conc);  

            bool all_quantities_correct = b_nmoments && b_nrandom && b_numdis && b_mult && b_Lx && b_Ly && b_seed && b_anderson && b_nx && b_ny && b_conc;


            verbose2("Checking if all the quantities are correctly formatted.");
            if(!all_quantities_correct){
                std::cout << "Some input values have incorrect formatting: ";
                if(!b_nmoments)  std::cout << "nmoments ";
                if(!b_nrandom)   std::cout << "nrandom ";
                if(!b_numdis)    std::cout << "ndisorder ";
                if(!b_mult)      std::cout << "mult ";
                if(!b_Lx)        std::cout << "geometry ";
                if(!b_Ly)        std::cout << "geometry ";
                if(!b_nx)        std::cout << "threads ";
                if(!b_ny)        std::cout << "threads ";
                if(!b_seed)      std::cout << "seed ";
                if(!b_anderson)  std::cout << "anderson ";
                if(!b_conc)      std::cout << "concentration ";
                std::cout << ". Aborting.\n";
                exit(1);
            }



        } catch (int e){
            std::cout << "An exception occured while reading the input "
                "from the command line. Make sure all the quantities "
                " are in the correct format. Aborting.\n";
            exit(1);
        }
    }


    P->need_read = need_read;
    verbose2("Checking if need to read from restart file\n");
    if(need_read){

        // number of extra chebyshev moments to be calculated
        if(need_moremoments){
            std::string s_moremoments;
            s_moremoments  = *(p_moremoments+1);
            P->moremoments = atoi(s_moremoments.c_str());
            P->need_moremoments = need_moremoments;
        }

        std::string filename_read;
        filename_read  = *(p_readfrom+1);
        P->filename_read = filename_read;


        std::string name =  (*(p_readfrom+1)).c_str();

        H5::H5File * file = new H5::H5File(name, H5F_ACC_RDWR);


        // Fetch information about the Hamiltonian
        get_hdf5(&P->Lx, file, (char *) "/Lx");
        get_hdf5(&P->Ly, file, (char *) "/Ly");
        get_hdf5(&P->anderson_W, file, (char *) "/anderson_w");
        get_hdf5(&P->conc, file, (char *) "/concentration");
        get_hdf5(&P->mult, file, (char *) "/mult");

        // Fetch information about the Chebyshev iteration
        get_hdf5(&P->seed, file, (char *) "/seed");
        P->nrandom = 1;
        P->ndisorder = 1;

        get_hdf5(&P->nmoments, file, (char *) "/num_moments");

        file->close();
        delete file;

        //std::cout << "must read from";
    }



    
    
    
     
    // Part regarding the output
    verbose2("output\n");


    P->output_energies  = output_energies;
    P->need_write      = need_write;


    if(need_write){
        std::string filename_write;
        filename_write = *(p_saveto+1);
        P->filename_write = filename_write;
    }




    // the number of energies has been specified
    if(output_energies){
        verbose2("The number of energies has been specified.\n");
        const double lim = 0.99;
        if(!found_energies && found_NEnergies){
            P->NEnergies = atoi((*(p_NEnergies+1)).c_str());
            P->energies  = Eigen::Array<TR, -1, 1>::LinSpaced(P->NEnergies, -lim, lim)*SCALE;
        }

        // Specific values of the energies have been specified
        verbose2("Specific energies have been specified.\n");
        if(found_energies && !found_NEnergies){
            debug_message("found energies\n");
            
            unsigned n_energies = 0;
            for(unsigned i = 0; i < 999; i++){
                std::string *current = p_energies + i + 1;

                // Check if the end has been reached
                if(current == inputs + argc -1){
                    break;
                }

                if(current == p_geometry or 
                   current == p_threads or
                   current == p_mult or
                   current == p_anderson or
                   current == p_conc or
                   current == p_nrandom or
                   current == p_nmoments or
                   current == p_numdis or
                   current == p_seed or
                   current == p_NEnergies or
                   current == p_energies or
                   current == p_help or
                   current == p_readfrom or
                   current == p_saveto or
                   current == p_moremoments or
                   current == p_print_to_cout or
                   current == p_print_to_file){
                    break;
                }
                n_energies++;

            }

            P->energies = Eigen::Array<TR, -1, -1>::Zero(n_energies, 1);
            for(unsigned i = 0; i < n_energies; i++){
                (P->energies)(i) = atof((*(p_energies + i +1)).c_str());
            }
        }
    }
    debug_message("Left aux parse_input.\n");
    delete []inputs;
}

void extended_euclidean(int a, int b, int *x, int *y, int *gcd){
    debug_message("Entered extended_euclidean.\n");
    // Uses the extended euclidean algorithm to calculate
    // The gauge matrix. The diophantine equation we need
    // to solve is
    //
    //          a*x + b*y = gcd(a,b)
    //
    // We use the algorithm to compute the gcd, and then use 
    // the gcd to compute x and y
    

    int s = 0;
    int old_s = 1;
    int t = 1;
    int old_t = 0;
    int r = a;
    int old_r = b;
    int quotient;
    int temp;

    unsigned num_iterations = 0;
    while(r != 0){
        quotient = old_r/r;
        temp  = old_r;
        old_r = r;
        r     = temp - quotient*r;

        temp  = old_s;
        old_s = s;
        s     = temp - quotient*s;

        temp  = old_t;
        old_t = t;
        t     = temp - quotient*t;
        num_iterations++;
        if(num_iterations > 9999){
            std::cout << "Extended euclidean algorithm is not converging. Check for errors. Aborting.\n";
            exit(1);
        }
    }

    *gcd = old_r;
    *x = old_t;
    *y = old_s;


    std::cout << "a,b: " << a << "," << b << " x,y: " << *x << "," << *y << " gcd: " << *gcd << "\n";
    debug_message("Left extended_euclidean.\n");
}

void save(KPM_vector *KPM0, KPM_vector *KPM1){
    static Eigen::Array<T, -1, -1> kpmA, kpmB;
    int Lx = KPM0->Lattice_Lx;
    int Ly = KPM0->Lattice_Ly;
#pragma omp master
    {
    kpmA = Eigen::Array<T, -1,- 1>::Zero(Ly, Lx);
    kpmB = Eigen::Array<T, -1,- 1>::Zero(Ly, Lx);
    //kpm(1,1) = 3;
    //H5::H5File * file1 = new H5::H5File(name, H5F_ACC_TRUNC);

    //for(unsigned n = 0; n < H->Norb; n++){
        //write_hdf5(KPM0.KPM[n], file1, "/KPM0_"+std::to_string(n));
        //write_hdf5(KPM1.KPM[n], file1, "/KPM1_"+std::to_string(n));
    }
#pragma omp barrier
#pragma omp critical
    {
        std::cout << "inside critical " << omp_get_thread_num() << "\n" << std::flush;
        unsigned tx, ty;
        unsigned id;
        id = KPM0->thread_id;
        tx = id%KPM0->NTx;
        ty = id/KPM0->NTx;

        unsigned size_x, size_y;
        size_y = Ly/double(KPM0->NTy) + 0.99;
        size_x = Lx/double(KPM0->NTx) + 0.99;

        unsigned i, j;
        i = ty*size_y;
        j = tx*size_x;
        std::cout << "i,j " << i << " " << j << " " << KPM0->Ly << " " << KPM0->Lx << "\n" << std::flush;

        kpmA.block(i,j,KPM0->Ly,KPM0->Lx) = KPM0->KPM[0].block(1,1,KPM0->Ly,KPM0->Lx);
        kpmB.block(i,j,KPM0->Ly,KPM0->Lx) = KPM0->KPM[1].block(1,1,KPM0->Ly,KPM0->Lx);

        std::cout << "PARTS\n";
        std::cout << "A:\n" << KPM0->KPM[0] << "\n";
        std::cout << "B:\n" << KPM0->KPM[1] << "\n";

        std::cout << "TOTAL\n";
        std::cout << "A:\n" << kpmA << "\n";
        std::cout << "B:\n" << kpmB << "\n";
    }

    //file1->close();
    //delete file1;
#pragma omp barrier
    verbose2("Left save");
}

void load(std::string name, KPM_vector *KPM0, KPM_vector *KPM1, Eigen::Array<TR, -1, -1> *mu){
    debug_message("Entered load\n");

    H5::H5File * file = new H5::H5File(name, H5F_ACC_RDWR);


    verbose2("Reading number of moments from restart file.\n");
    unsigned moms;
    get_hdf5(&moms, file, (char *) "/num_moments");
    verbose2("number of moments: " << moms);

    verbose2("Reading MU array from restart file.\n");
    Eigen::Array<TR, -1, -1> mu1;
    //std::cout << "before allocation\n" << std::flush; 
    mu1 = Eigen::Array<TR, -1, -1>::Zero(moms, 1);
    //std::cout << "before get_hdf\n" << std::flush; 
    get_hdf5(mu1.data(), file, (char *) "/MU");

    //std::cout << "before copy\n" << std::flush; 
    *mu = mu1;
    //std::cout << "after copy\n" << std::flush; 


    verbose2("Reading KPM vectors from restart file.\n");
    unsigned Norb = 2;
    for(unsigned n = 0; n < Norb; n++){
        std::string name1 = "/KPM0_" + std::to_string(n);
        std::string name2 = "/KPM1_" + std::to_string(n);
        get_hdf5(KPM0->KPM[n].data(), file, name1);
        get_hdf5(KPM1->KPM[n].data(), file, name2);
    }


    file->close();
    delete file;

    debug_message("Left load\n");
};

double jackson(unsigned n, unsigned N){
    double arg = M_PI/(N + 1);
    double factor = (N - n + 1)*cos(arg*n) + sin(arg*n)/tan(arg);
    return factor/(N + 1);
}


void print_ham_info(parameters P){
    verbose1("______Hamiltonian Parameters______\n");
    verbose1("Physical system:     Graphene\n");
    verbose1("Lx:                  " + std::to_string(P.Lx)           + "\n");
    verbose1("Ly:                  " + std::to_string(P.Ly)           + "\n");
    verbose1("mult:                " + std::to_string(P.mult)         + "\n");
    verbose1("seed:                " + std::to_string(P.seed)         + "\n");
    verbose1("anderson:            " + std::to_string(P.anderson_W)   + "\n");
    verbose1("vacancies:           " + std::to_string(P.conc)         + "\n");
}
void print_cheb_info(parameters P){
    verbose1("______Cheb iteration parameters_______\n");
    verbose1("num random:          " + std::to_string(P.nrandom)   + "\n");
    verbose1("num disorder:        " + std::to_string(P.ndisorder) + "\n");
#if VERBOSE >= 1
    unsigned moments_to_calculate = P.nmoments;
#endif
    if(P.need_read){
#if VERBOSE >= 1
    moments_to_calculate = P.moremoments;
#endif
        verbose1("num moments to calc: " + std::to_string(moments_to_calculate) + "\n");
        verbose1("total num moments:   " + std::to_string(moments_to_calculate + P.nmoments) + "\n");
    } else {
        verbose1("num moments:         " + std::to_string(moments_to_calculate) + "\n");
    }
}

void print_compilation_info(){
    verbose1("_______Compilation parameters_______\n");
    verbose1("VERBOSE: " << VERBOSE << "\n");
    verbose1("SCALE: " << SCALE << "\n");
}

void print_magnetic_info(parameters P, int M12, int M21, int min_flux){
    verbose1("___________Magnetic info___________\n");
    verbose1("Minimum flux:            " + std::to_string(min_flux) + "\n");
    verbose1("Flux:                    " + std::to_string(min_flux) + "\n");
    verbose1("Gauge matrix M(0,1)=M12: " + std::to_string(M12*int(P.mult)) + "\n");
    verbose1("Gauge matrix M(1,0)=M21: " + std::to_string(-M21*int(P.mult)) + "\n");
}

void print_output_info(parameters P){
    verbose1("____________Output info____________\n");
    if(P.energies.size() > 10){
        verbose1("Calculating DoS at energies:    ");
        verbose1(P.energies.topRows(2).transpose()*SCALE);
        verbose1(" ... (" + std::to_string(P.energies.size() - 4) + " more) ... ");
        verbose1(P.energies.bottomRows(2).transpose()*SCALE);
        verbose1("\n");
    } else {
        verbose1("Calculating DoS at energies:    " << P.energies.transpose()*SCALE     << "\n");
    }
}


void print_log_info(parameters P){
    verbose1("____________Log info____________\n");
    if(P.need_log_status){
        verbose1("Log being saved to directory: " << P.filename_status << "\n");
    } else {
        verbose1("No logging information requested.\n");
    }

}

void print_restart_info(parameters P){
    verbose1("____________Restart info____________\n");
    if(P.need_read){
        verbose1("Reading restart file " << P.filename_read << ".\n");
        if(P.need_moremoments){
            verbose1("Calculating " << P.moremoments << " more moments.\n");
        }
    } else {
        verbose1("Reading from restart file not requested.\n");
    }
    if(P.need_write){
        verbose1("Saving to restart file " << P.filename_write << "\n");
    } else {
        verbose1("Saving to restart file not requested.\n");
    }

}



