#include	<stdlib.h>
#include	<iostream>
#include	<sstream>
#include	<fstream>
#include	<cmath>
#include	<immintrin.h>
#include    <mpi.h>

#include	"Timer.h"

#include	"grid.hpp"
#include	"util.hpp"

int main(int argc, char **argv) {
	///******************************************************
	///********************** INPUT *************************
	///******************************************************
	Params params(argc, argv);
		
	params.bx      = params.nx;
	params.by      = params.ny;
	params.offsetX = 0;
	params.offsetY = 0;
	
	Grid	u(params.bx, params.by);

	///******************************************************
	///********************** CALCULATION *******************
	///******************************************************
	double time = 0;

	siwir::Timer	timer;

	time = timer.elapsed();
	std::cout << "time," << time << std::endl;

	///******************************************************
	///********************** OUTPUT ************************
	///******************************************************
};
