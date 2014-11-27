#include	<stdlib.h>
#include	<iostream>
#include	<sstream>
#include	<fstream>
#include	<cmath>
#include	<immintrin.h>
#include    <mpi.h>

#ifdef USE_LIKWID
extern "C" {
#include <likwid.h>
}
#endif

#include	"Timer.h"

static const double k = 2*M_PI;

/**
  Converts a string to an arbitrary type. >> operator must be defined for the target type.
  @param string string which should be converted
  @return converted string
 **/
template<typename T>
T StringTo(const std::string& string){
	T valor;

	std::stringstream stream(string);
	stream >> valor;
	return valor;
}

constexpr	double f(const double x, const double y){
	return 4.0*M_PI*M_PI*sin(2.0*M_PI*x)*sinh(2.0*M_PI*y);
}

constexpr	double border(const double x, const double y){
	return sin(2.0*M_PI*x)*sinh(2.0*M_PI*y);
}

int main(int argc, char **argv) {
	///******************************************************
	///********************** INPUT *************************
	///******************************************************
    if (argc != 5) {
		std::cout << "Invalid number of arguments!" << std::endl;
		std::cout << "./rgbs nx ny c" << std::endl;
		exit(EXIT_FAILURE);
	}

    // MPI initialization
    MPI_Init( &argc, &argv );

    // Variable loading
    int nx = 0;
	int	ny = 0;
	int	c = 0;
    double eps = 0.0;
	nx = StringTo<int>(argv[1]);
	ny = StringTo<int>(argv[2]);
	c = StringTo<int>(argv[3]);
    eps = StringTo<double>(argv[4]);
	
	std::cout << "nx," << nx << std::endl;
	std::cout << "ny," << ny << std::endl;
	std::cout << "c," << c <<std::endl;
    std::cout << "eps,," << eps <<std::endl;

	///******************************************************
	///********************** CALCULATION *******************
	///******************************************************
	double time = 0;
	
#ifdef USE_LIKWID
	likwid_markerInit();
	likwid_markerStartRegion("dummy");
#endif

	siwir::Timer	timer;
	
	double	hx = 2.0/nx;
	double	hy = 1.0/ny;
	
	double	invHx2 = 1.0/hx/hx;
	double	invHy2 = 1.0/hy/hy;
	
	double	preF = 1.0 / (2*invHx2 + 2*invHy2 + k*k);
    std::cout << preF << std::endl;

    int size(0); // The total number of processes
    int rank(0); // The rank/number of this process

    // Determining the number of CPUs and the rank of this process
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    // 'Hello World' output on process 0
    if( rank == 0 )
       std::cout << "Hello World!" << std::endl;

    // Output of the process rank for each process
    std::cout << "I am CPU " << rank << " of " << size << " CPUs" << std::endl;
	
	time = timer.elapsed();
	std::cout << "time," << time << std::endl;

#ifdef USE_LIKWID
	likwid_markerStopRegion("dummy");
	likwid_markerClose();
#endif

	///******************************************************
	///********************** OUTPUT ************************
	///******************************************************
	
    // MPI finalization
    MPI_Finalize();
};
