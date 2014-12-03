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
	
	// Initialization of MPI
	// ----------------------------------------------------------------
	MPI_Init( &argc, &argv );
	// ----------------------------------------------------------------
	
	Params params(argc, argv);

	// Determining the number of CPUs and the rank of this process
	// ----------------------------------------------------------------
	MPI_Comm_size( MPI_COMM_WORLD, &params.size );
	MPI_Comm_rank( MPI_COMM_WORLD, &params.rank );
	// ----------------------------------------------------------------
	
	params.subdivideGrid();

	// The new MPI communicator for the Cartesian topology
	MPI_Comm cartcomm( MPI_COMM_NULL );

	// Creating the Cartesian topology:
	//  - Creating a 2D grid with 2 processes in x- and 2 processes in y-direction
	//  - Determining the coordinates of the processes
	//  - Determining the neighbors in UP, DOWN, LEFT and RIGHT direction
	// ----------------------------------------------------------------
	MPI_Cart_create( MPI_COMM_WORLD, 2, params.dims, params.periods, params.reorder, &cartcomm );
	MPI_Comm_rank( cartcomm, &params.cartrank );
	MPI_Cart_coords( cartcomm, params.cartrank, 2, params.coords );
	MPI_Cart_shift( cartcomm, 0, 1, &params.nbrs[Params::LEFT], &params.nbrs[Params::RIGHT] );
	MPI_Cart_shift( cartcomm, 1, 1, &params.nbrs[Params::DOWN], &params.nbrs[Params::UP] );
	// ----------------------------------------------------------------
	
	params.createBlock();
	
	Grid	lookupF(params.bx, params.by);
	
	for (int y = 0; y < params.by; ++y){
		for (int x = 0; x < params.bx; ++x){
			lookupF(x, y) = f(params.getXCoord(x, y), params.getYCoord(x,y));
		}
	}
	
	Grid	u(params.bx, params.by);
	for (int x = 0; x < params.bx; ++x){
		u(x, params.by) = border(params.getXCoord(x, params.by), params.getYCoord(x,params.by));
	}
	
	Grid	r(params.bx, params.by);
	Grid	d(params.bx, params.by);
	Grid	z(params.bx, params.by);
	
	double	delta0 = 0.0;

	///******************************************************
	///********************** CALCULATION *******************
	///******************************************************
	double time = 0;

	siwir::Timer	timer;
	
	for (int y = 0; y < params.by; ++y){
		for (int x = 0; x < params.bx; ++x){
			r(x, y) = lookupF(x, y) + params.invHx2 * ( u(x - 1, y) + u(x + 1, y) ) + params.invHy2 * ( u(x, y - 1) + u(x, y + 1) ) - params.preF * u(x, y);
			d(x, y) = r(x,y);
			delta0 += r(x,y) * r(x,y);
		}
	}
	
	//if (delta0 < params.eps2) {...}
	
	for (int c = 0; c < params.c; ++c){
		double	alpha = 0.0;
		for (int y = 0; y < params.by; ++y){
			for (int x = 0; x < params.bx; ++x){
				z(x, y) = - params.invHx2 * ( d(x - 1, y) + d(x + 1, y) ) - params.invHy2 * ( d(x, y - 1) + d(x, y + 1) ) + params.preF * d(x, y);
				alpha += z(x, y) * d(x, y);
			}
		}
		alpha = delta0 / alpha;
		double	delta1 = 0.0;
		for (int y = 0; y < params.by; ++y){
			for (int x = 0; x < params.bx; ++x){
				u(x, y) += alpha * d(x, y);
				r(x, y) -= alpha * z(x, y);
				delta1 += r(x, y) * r(x, y);
			}
		}
		if (delta1 < params.eps2) break;
		double beta = delta1 / delta0;
		for (int y = 0; y < params.by; ++y){
			for (int x = 0; x < params.bx; ++x){
				d(x, y) = r(x, y) + beta * d(x,y);
			}
		}
		delta0 = delta1;
		//std::cout << delta0 << std::endl;
	}

	MPI_Barrier( MPI_COMM_WORLD );
	time = timer.elapsed();
	if (params.rank == 0){
		std::cout << "time," << time << std::endl;
	}

	///******************************************************
	///********************** OUTPUT ************************
	///******************************************************
	/*
	std::ofstream	fOut("data/solution.txt");
	for (int y = -params.offsetY; y < params.by + params.offsetY; ++y) {
		for (int x = -params.offsetX; x < params.bx + params.offsetX; ++x) {
			fOut << params.getXCoord(x, y) << "\t" << params.getYCoord(x, y) << "\t" << u(x, y) << std::endl;
		}
		fOut << std::endl;
	}
	fOut.close();
	*/
	
	// MPI finalizations
	// ----------------------------------------------------------------
	MPI_Finalize();
	// ----------------------------------------------------------------

	return 0;
};
