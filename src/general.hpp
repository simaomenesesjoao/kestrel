
#if VERBOSE >= 3
#define debug_message(VAR) std::cout<<VAR<<std::flush;
#define verbose3(VAR) std::cout<<VAR<<std::flush;
#endif

#if VERBOSE < 3
#define debug_message(VAR)
#define verbose3(VAR)
#endif



#if VERBOSE >= 2
#define verbose2(VAR) std::cout <<VAR<<std::flush;
#endif

#if VERBOSE < 2
#define verbose2(VAR)
#endif

#if VERBOSE >= 1
#define verbose1(VAR) std::cout <<VAR<<std::flush;
#endif

#if VERBOSE < 1
#define verbose1(VAR)
#endif


#define SCALE 3.5

typedef double TR;
#if nis_complex == 1
typedef std::complex<double> T;
#endif
#if nis_complex == 0
typedef double T;
#endif


typedef Eigen::Array<T, -1, -1> Array;
typedef Eigen::Matrix<T, -1, -1> Matrix;
