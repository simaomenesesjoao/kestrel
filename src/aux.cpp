#include <iostream>
#include <Eigen/Dense>
#include "general.hpp"
#include "aux.hpp"



Eigen::Array<TR, -1, -1> calc_dos(Eigen::Array<TR, -1, -1> mu, Eigen::Array<TR, -1, -1> energies, std::string mode){
    debug_message("Entered calc_dos\n");

    unsigned N_energies = energies.size();
    unsigned moments = mu.size();
    std::complex<TR> im(0.0, 1.0);
    Eigen::Array<TR, -1, -1> dos(N_energies, 1);


    if(mode == "jackson"){
        verbose2("Entered mode: jackson\n");
        Eigen::Array<TR, -1, -1> temp[2], sq_en(N_energies, 1);
        sq_en = Eigen::rsqrt(1.0 - energies*energies)*2.0/M_PI;
        temp[0] = Eigen::Array<TR, -1, -1>::Zero(N_energies, 1) + 1;
        temp[1] = energies;

        dos  = temp[0]*mu(0)/2.0*jackson(0, moments);
        dos += temp[1]*mu(1)*jackson(1, moments);
        for(unsigned i = 1; i < moments; i++){
            temp[i%2] = 2*energies*temp[(i+1)%2] - temp[i%2];
            dos += temp[i%2]*mu(i)*jackson(i, moments);
        }
        dos = dos*sq_en;
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

    return dos/SCALE;

}

void parse_input(int argc, char **argv, unsigned *Lx, unsigned *Ly, unsigned *nrandom, unsigned *nmoments, unsigned *mult, int *seed, Eigen::Array<TR, -1, -1> *energies, double *W, unsigned *num_disorder){

    // Convert everything to a string to make it easier to process
    std::string *inputs;
    inputs = new std::string[argc-1];
    for(int i = 0; i < argc-1; i++){
        inputs[i] = std::string(argv[i+1]);
    }

    std::string *p_geometry, *p_mult, *p_nrandom, *p_nmoments, *p_energies, *p_seed, *p_anderson, *p_num_dis;

    // find the pointers to each of the quantities
    p_geometry  = std::find(inputs, inputs+argc-1, "--geometry");
    p_mult      = std::find(inputs, inputs+argc-1, "--mult");
    p_nrandom   = std::find(inputs, inputs+argc-1, "--nrandom");
    p_num_dis   = std::find(inputs, inputs+argc-1, "--ndisorder");
    p_nmoments  = std::find(inputs, inputs+argc-1, "--nmoments");
    p_energies  = std::find(inputs, inputs+argc-1, "--energies");
    p_seed      = std::find(inputs, inputs+argc-1, "--seed");
    p_anderson  = std::find(inputs, inputs+argc-1, "--anderson");


    // Find the values of each of the quantities
    // NOTE: should make sure the values exist at all
    *nmoments       = atoi(( *(p_nmoments+1)  ).c_str());
    *nrandom        = atoi(( *(p_nrandom+1 )  ).c_str());
    *num_disorder   = atoi(( *(p_num_dis+1 )  ).c_str());
    *mult           = atoi(( *(p_mult+1    )  ).c_str());
    *Lx             = atoi(( *(p_geometry+1)  ).c_str());
    *Ly             = atoi(( *(p_geometry+2)  ).c_str());
    *seed           = atoi(( *(p_seed+1)      ).c_str());
    *W              = atof(( *(p_anderson+1)  ).c_str());

    // Find the energies
    unsigned n_energies = 0;
    for(unsigned i = 0; i < 999; i++){
        //std::cout << "i: " << i << " ";
        std::string *current = p_energies + i + 1;
        //std::cout << *current << "\n";

        // Check if the end has been reached
        if(current == inputs + argc -1){
            break;
        }

        if(current == p_geometry or 
           current == p_nrandom or
           current == p_nmoments or
           current == p_mult){
            break;
        }
        n_energies++;

    }

    *energies = Eigen::Array<TR, -1, -1>::Zero(n_energies, 1);
    for(unsigned i = 0; i < n_energies; i++){
        (*energies)(i) = atof((*(p_energies + i +1)).c_str());
    }
}

void extended_euclidean(int a, int b, int *x, int *y, int *gcd){
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


}

double jackson(unsigned n, unsigned N){
    double arg = M_PI/(N + 1);
    double factor = (N - n + 1)*cos(arg*n) + sin(arg*n)/tan(arg);
    return factor/(N + 1);
}
