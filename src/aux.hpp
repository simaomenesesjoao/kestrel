
double jackson(unsigned n, unsigned N);
Eigen::Array<TR, -1, -1> calc_dos(Eigen::Array<TR, -1, -1> mu, Eigen::Array<TR, -1, -1> energies, std::string mode);
void parse_input(int argc, char **argv, unsigned *Lx, unsigned *Ly, unsigned *nrandom, unsigned *nmoments, unsigned *mult, int *seed, Eigen::Array<TR, -1, -1> *energies, double *W, unsigned *num_disorder);
void extended_euclidean(int a, int b, int *x, int *y, int *gcd);
