#include	<stdlib.h>
#include	<cstdio>
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
	
	params.printConfiguration();
	
	params.subdivideGrid();

	// Creating the Cartesian topology:
	//  - Creating a 2D grid with 2 processes in x- and 2 processes in y-direction
	//  - Determining the coordinates of the processes
	//  - Determining the neighbors in UP, DOWN, LEFT and RIGHT direction
	// ----------------------------------------------------------------
	MPI_Cart_create( MPI_COMM_WORLD, 2, params.dims, params.periods, params.reorder, &params.cartcomm );
	MPI_Comm_rank( params.cartcomm, &params.cartrank );
	MPI_Cart_coords( params.cartcomm, params.cartrank, 2, params.coords );
	MPI_Cart_shift( params.cartcomm, 0, 1, &params.nbrs[Params::LEFT], &params.nbrs[Params::RIGHT] );
	MPI_Cart_shift( params.cartcomm, 1, 1, &params.nbrs[Params::DOWN], &params.nbrs[Params::UP] );
	params.nbrs[Params::UPLEFT]    = params.getNeighbour(Params::UPLEFT);
	params.nbrs[Params::UPRIGHT]   = params.getNeighbour(Params::UPRIGHT);
	params.nbrs[Params::DOWNLEFT]  = params.getNeighbour(Params::DOWNLEFT);
	params.nbrs[Params::DOWNRIGHT] = params.getNeighbour(Params::DOWNRIGHT);
	//if (params.rank==1)
		//std::cout << MPI_PROC_NULL << "\t" << params.nbrs[Params::UPLEFT] << "\t" << params.nbrs[Params::UPRIGHT] << "\t" << params.nbrs[Params::DOWNLEFT] << "\t" << params.nbrs[Params::DOWNRIGHT];
	// ----------------------------------------------------------------
	
	params.createBlock();
	
	Grid	lookupF(params.bx, params.by, params.numGhostLayers);
	
	for (int y = 0; y < params.by; ++y){
		for (int x = 0; x < params.bx; ++x){
			lookupF(x, y) = f(params.getXCoord(x, y), params.getYCoord(x,y));
		}
	}
	
	Grid	u(params.bx, params.by, params.numGhostLayers);
	if (params.dims[1]-1 == params.coords[1]){
		//std::cout << "border," << params.rank << std::endl;
		for (int x = 1 - params.numGhostLayers; x < params.bx + params.numGhostLayers - 1; ++x){
			u(x, params.by) = border(params.getXCoord(x, params.by), params.getYCoord(x,params.by));
		}
	}
	
	Grid	r(params.bx, params.by, params.numGhostLayers);
	Grid	d(params.bx, params.by, params.numGhostLayers);
	Grid	z(params.bx, params.by, params.numGhostLayers);
	
	// Creating a new derived data type
	MPI_Datatype	verticalBorderType;   
	MPI_Type_vector( params.by, 1, u.ld, MPI_DOUBLE, &verticalBorderType );
	MPI_Type_commit( &verticalBorderType );
	
	// Creating a new derived data type
	MPI_Datatype	horizontalBorderType;   
	MPI_Type_vector( params.numGhostLayers, params.bx + 2*params.numGhostLayers - 2, u.ld, MPI_DOUBLE, &horizontalBorderType );
	MPI_Type_commit( &horizontalBorderType );
	
	// Creating a new derived data type
	MPI_Datatype	cornerBorderType;   
	MPI_Type_vector( params.numGhostLayers, params.numGhostLayers, u.ld, MPI_DOUBLE, &cornerBorderType );
	MPI_Type_commit( &cornerBorderType );
	
	MPI_Request reqs[16];
	MPI_Status stats[16];
	
	double	localDelta0 = 0.0;
	double	delta0      = 0.0;
	double	delta1      = 0.0;
	double	alpha       = 0.0;

	///******************************************************
	///********************** CALCULATION *******************
	///******************************************************
	double time         = 0.0;

	siwir::Timer	timer;
	
	for (int y = 0; y < params.by; ++y){
		for (int x = 0; x < params.bx; ++x){
			r(x, y) = lookupF(x, y) + params.invHx2 * ( u(x - 1, y) + u(x + 1, y) ) + params.invHy2 * ( u(x, y - 1) + u(x, y + 1) ) - params.preF * u(x, y);
			d(x, y) = r(x,y);
			localDelta0 += r(x,y) * r(x,y);
		}
	}
	//delta0 = 0;
	MPI_Allreduce(&localDelta0, &delta0, 1, MPI_DOUBLE, MPI_SUM, params.cartcomm);
	//if (params.rank==0) std::cout << "delta0, " << delta0 << std::endl;
	
	//if (delta0 < params.eps2) {...}
	
	MPI_Isend( &d(0, 0), 1, verticalBorderType, params.nbrs[Params::LEFT], 0, params.cartcomm, &reqs[0]   );
	MPI_Isend( &d(params.bx-1, 0), 1, verticalBorderType, params.nbrs[Params::RIGHT], 0, params.cartcomm, &reqs[1]   );
	MPI_Isend( &d(0, 0), 1, horizontalBorderType, params.nbrs[Params::DOWN], 0, params.cartcomm, &reqs[2]   );
	MPI_Isend( &d(0, params.by-params.numGhostLayers), 1, horizontalBorderType, params.nbrs[Params::UP], 0, params.cartcomm, &reqs[3]   );
	
	MPI_Isend( &d(0, params.by - params.numGhostLayers), 1, cornerBorderType, params.nbrs[Params::UPLEFT], 0, params.cartcomm, &reqs[8]   );
	MPI_Isend( &d(params.bx - params.numGhostLayers, params.by - params.numGhostLayers), 1, cornerBorderType, params.nbrs[Params::UPRIGHT], 0, params.cartcomm, &reqs[9]   );
	MPI_Isend( &d(0, 0), 1, cornerBorderType, params.nbrs[Params::DOWNLEFT], 0, params.cartcomm, &reqs[10]   );
	MPI_Isend( &d(params.bx - params.numGhostLayers, 0), 1, cornerBorderType, params.nbrs[Params::DOWNRIGHT], 0, params.cartcomm, &reqs[11]   );
	
	MPI_Irecv( &d(-1, 0), 1, verticalBorderType, params.nbrs[Params::LEFT], 0, params.cartcomm, &reqs[4] );
	MPI_Irecv( &d(params.bx, 0), 1, verticalBorderType, params.nbrs[Params::RIGHT], 0, params.cartcomm, &reqs[5] );
	MPI_Irecv( &d(0, -params.numGhostLayers), params.bx, horizontalBorderType, params.nbrs[Params::DOWN], 0, params.cartcomm, &reqs[6] );
	MPI_Irecv( &d(0, params.by), params.bx, horizontalBorderType, params.nbrs[Params::UP], 0, params.cartcomm, &reqs[7] );
	
	MPI_Irecv( &d(-params.numGhostLayers, params.by), 1, cornerBorderType, params.nbrs[Params::UPLEFT], 0, params.cartcomm, &reqs[12] );
	MPI_Irecv( &d(params.bx, params.by), 1, cornerBorderType, params.nbrs[Params::UPRIGHT], 0, params.cartcomm, &reqs[13] );
	MPI_Irecv( &d(-params.numGhostLayers, -params.numGhostLayers), 1, cornerBorderType, params.nbrs[Params::DOWNLEFT], 0, params.cartcomm, &reqs[14] );
	MPI_Irecv( &d(params.bx, -params.numGhostLayers), 1, cornerBorderType, params.nbrs[Params::DOWNRIGHT], 0, params.cartcomm, &reqs[15] );
	
	MPI_Waitall( 16, reqs, stats );
	
	int c = 0;
	
	while (c < params.c){
		int	c2 = c + params.numGhostLayers;
		if (c >= params.c) c = params.c;
		double	cInnerLoop = 0;
		while (c < c2){
			//if (params.rank == 0) std::cout << c << "\t" << c2 << std::endl;
			for (int y = 1 - params.numGhostLayers + cInnerLoop; y < params.by + params.numGhostLayers - 1 - cInnerLoop; ++y){
				for (int x = 1 - params.numGhostLayers + cInnerLoop; x < params.bx + params.numGhostLayers - 1 - cInnerLoop; ++x){
					z(x, y) = - params.invHx2 * ( d(x - 1, y) + d(x + 1, y) ) - params.invHy2 * ( d(x, y - 1) + d(x, y + 1) ) + params.preF * d(x, y);
				}
			}
			
			double	localAlpha = 0.0;
			for (int y = 0; y < params.by; ++y){
				for (int x = 0; x < params.bx; ++x){
					localAlpha += z(x, y) * d(x, y);
				}
			}
			MPI_Allreduce(&localAlpha, &alpha, 1, MPI_DOUBLE, MPI_SUM, params.cartcomm);
			
			//if (params.rank==0) std::cout << "alpha, " << alpha << std::endl;
			alpha = delta0 / alpha;
			
			for (int y = 1 - params.numGhostLayers + cInnerLoop; y < params.by + params.numGhostLayers - 1 - cInnerLoop; ++y){
				for (int x = 1 - params.numGhostLayers + cInnerLoop; x < params.bx + params.numGhostLayers - 1 - cInnerLoop; ++x){
					u(x, y) += alpha * d(x, y);
					r(x, y) -= alpha * z(x, y);
				}
			}
			
			double	localDelta1 = 0.0;
			for (int y = 0; y < params.by; ++y){
				for (int x = 0; x < params.bx; ++x){
					localDelta1 += r(x, y) * r(x, y);
				}
			}
			MPI_Allreduce(&localDelta1, &delta1, 1, MPI_DOUBLE, MPI_SUM, params.cartcomm);
			
			//if (params.rank==0) std::cout << "delta1, " << delta1 << std::endl;
			if (delta1 < params.eps2) break;
			double beta = delta1 / delta0;
			for (int y = 1 - params.numGhostLayers + cInnerLoop; y < params.by + params.numGhostLayers - 1 - cInnerLoop; ++y){
				for (int x = 1 - params.numGhostLayers + cInnerLoop; x < params.bx + params.numGhostLayers - 1 - cInnerLoop; ++x){
					d(x, y) = r(x, y) + beta * d(x,y);
				}
			}
			delta0 = delta1;
			++cInnerLoop;
			++c;
		}
		
		MPI_Isend( &d(0, 0), 1, verticalBorderType, params.nbrs[Params::LEFT], 0, params.cartcomm, &reqs[0]   );
		MPI_Isend( &d(params.bx-1, 0), 1, verticalBorderType, params.nbrs[Params::RIGHT], 0, params.cartcomm, &reqs[1]   );
		MPI_Isend( &d(0, 0), 1, horizontalBorderType, params.nbrs[Params::DOWN], 0, params.cartcomm, &reqs[2]   );
		MPI_Isend( &d(0, params.by-params.numGhostLayers), 1, horizontalBorderType, params.nbrs[Params::UP], 0, params.cartcomm, &reqs[3]   );
		
		MPI_Isend( &d(0, params.by - params.numGhostLayers), 1, cornerBorderType, params.nbrs[Params::UPLEFT], 0, params.cartcomm, &reqs[8]   );
		MPI_Isend( &d(params.bx - params.numGhostLayers, params.by - params.numGhostLayers), 1, cornerBorderType, params.nbrs[Params::UPRIGHT], 0, params.cartcomm, &reqs[9]   );
		MPI_Isend( &d(0, 0), 1, cornerBorderType, params.nbrs[Params::DOWNLEFT], 0, params.cartcomm, &reqs[10]   );
		MPI_Isend( &d(params.bx - params.numGhostLayers, 0), 1, cornerBorderType, params.nbrs[Params::DOWNRIGHT], 0, params.cartcomm, &reqs[11]   );
		
		MPI_Irecv( &d(-1, 0), 1, verticalBorderType, params.nbrs[Params::LEFT], 0, params.cartcomm, &reqs[4] );
		MPI_Irecv( &d(params.bx, 0), 1, verticalBorderType, params.nbrs[Params::RIGHT], 0, params.cartcomm, &reqs[5] );
		MPI_Irecv( &d(0, -params.numGhostLayers), params.bx, horizontalBorderType, params.nbrs[Params::DOWN], 0, params.cartcomm, &reqs[6] );
		MPI_Irecv( &d(0, params.by), params.bx, horizontalBorderType, params.nbrs[Params::UP], 0, params.cartcomm, &reqs[7] );
		
		MPI_Irecv( &d(-params.numGhostLayers, params.by), 1, cornerBorderType, params.nbrs[Params::UPLEFT], 0, params.cartcomm, &reqs[12] );
		MPI_Irecv( &d(params.bx, params.by), 1, cornerBorderType, params.nbrs[Params::UPRIGHT], 0, params.cartcomm, &reqs[13] );
		MPI_Irecv( &d(-params.numGhostLayers, -params.numGhostLayers), 1, cornerBorderType, params.nbrs[Params::DOWNLEFT], 0, params.cartcomm, &reqs[14] );
		MPI_Irecv( &d(params.bx, -params.numGhostLayers), 1, cornerBorderType, params.nbrs[Params::DOWNRIGHT], 0, params.cartcomm, &reqs[15] );
		
		MPI_Waitall( 16, reqs, stats );

		//std::cout << delta0 << std::endl;
	}

	MPI_Barrier( MPI_COMM_WORLD );
	time = timer.elapsed();
	if (params.rank == 0){
		std::cout << "residuum," << sqrt(delta0 / ( params.nx - 2 ) / ( params.ny - 2 )) << std::endl;
		std::cout << "time," << time << std::endl;
	}

	///******************************************************
	///********************** OUTPUT ************************
	///******************************************************
	
	MPI_Isend( &u(0, 0), 1, verticalBorderType, params.nbrs[Params::LEFT], 0, params.cartcomm, &reqs[0]   );
	MPI_Isend( &u(params.bx-1, 0), 1, verticalBorderType, params.nbrs[Params::RIGHT], 0, params.cartcomm, &reqs[1]   );
	MPI_Isend( &u(0, 0), params.bx, MPI_DOUBLE, params.nbrs[Params::DOWN], 0, params.cartcomm, &reqs[2]   );
	MPI_Isend( &u(0, params.by-1), params.bx, MPI_DOUBLE, params.nbrs[Params::UP], 0, params.cartcomm, &reqs[3]   );
	
	MPI_Irecv( &u(-1, 0), 1, verticalBorderType, params.nbrs[Params::LEFT], 0, params.cartcomm, &reqs[4] );
	MPI_Irecv( &u(params.bx, 0), 1, verticalBorderType, params.nbrs[Params::RIGHT], 0, params.cartcomm, &reqs[5] );
	MPI_Irecv( &u(0, -1), params.bx, MPI_DOUBLE, params.nbrs[Params::DOWN], 0, params.cartcomm, &reqs[6] );
	MPI_Irecv( &u(0, params.by), params.bx, MPI_DOUBLE, params.nbrs[Params::UP], 0, params.cartcomm, &reqs[7] );
	
	MPI_Waitall( 8, reqs, stats );
	
	localDelta0 = 0;
	for (int y = 0; y < params.by; ++y){
		for (int x = 0; x < params.bx; ++x){
			double temp = lookupF(x, y) + params.invHx2 * ( u(x - 1, y) + u(x + 1, y) ) + params.invHy2 * ( u(x, y - 1) + u(x, y + 1) ) - params.preF * u(x, y);
			localDelta0 += temp * temp;
		}
	}
	delta0 = 0;
	MPI_Allreduce(&localDelta0, &delta0, 1, MPI_DOUBLE, MPI_SUM, params.cartcomm);
	if (params.rank == 0){
		std::cout << "residuum2," << sqrt(delta0 / ( params.nx - 2 ) / ( params.ny - 2 )) << std::endl;
	}
	
	if (params.output){
		if (params.rank == 0) remove("data/solution.txt");
		MPI_Barrier( params.cartcomm );
		/// run through all rows
		for (int y = 0; y < params.ny; ++y){
			/// run through grid rows
			for (int dimY = 0; dimY < params.dims[1]; ++dimY){
				/// test if row is inside grid row
				if ((y>=(params.offsetY+params.getBottomBorderOffset())) && (y<(params.offsetY+params.getTopBorderOffset()+params.by))){
					/// run through grid cols
					for (int dimX = 0; dimX < params.dims[0]; ++dimX){
						if ((dimX == params.coords[0]) && (dimY == params.coords[1])){
							std::fstream	fOut("data/solution.txt", std::fstream::out | std::fstream::app);
							for (int x = 0 + params.getLeftBorderOffset(); x < params.bx + params.getRightBorderOffset(); ++x) {
								//std::cout << y << "\t" << dimY << "\t" << params.offsetY << "\t" << dimX << "\t" << x << std::endl;
								//std::cout << "block," << params.rank << "\t" << params.coords[0] << "\t" << params.coords[1] << "\t" << params.offsetX << "\t" << params.offsetY << "\t" << x << "\t" << y << std::endl;
								fOut << params.getXCoord(x, y - params.offsetY) << "\t" << params.getYCoord(x, y - params.offsetY) << "\t" << u(x, y - params.offsetY) << std::endl;
							}
							if (params.coords[0] == params.dims[0]-1) fOut << std::endl;
							fOut.close();						
						}
						MPI_Barrier( params.cartcomm );
					}
				} else {	
					for (int dimX = 0; dimX < params.dims[0]; ++dimX){
						MPI_Barrier( params.cartcomm );
					}
				}
			}		
		}
	}
	
	if (params.debugOutput){
		std::stringstream	ss;
		ss << "data/grid_" << params.coords[0] << "x" << params.coords[1] << ".txt";
		std::ofstream	fOut(ss.str());
		for (int y = -1; y < params.by + 1; ++y) {
			for (int x = -1; x < params.bx + 1; ++x) {
				fOut << params.getXCoord(x, y) << "\t" << params.getYCoord(x, y) << "\t" << u(x, y) << std::endl;
			}
			fOut << std::endl;
		}
		fOut.close();
	}
	
	// Freeing derived data type
	MPI_Type_free( &verticalBorderType );	
	MPI_Type_free( &horizontalBorderType );
	MPI_Type_free( &cornerBorderType );
	
	// MPI finalizations
	// ----------------------------------------------------------------
	MPI_Finalize();
	// ----------------------------------------------------------------

	return 0;
};
