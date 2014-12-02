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
		
	params.bx      = params.nx - 2;
	params.by      = params.ny - 2;
	params.offsetX = 1;
	params.offsetY = 1;
	
	
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
			delta0 = r(x,y) * r(x,y);
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
		std::cout << delta0 << std::endl;
	}

	time = timer.elapsed();
	std::cout << "time," << time << std::endl;

	///******************************************************
	///********************** OUTPUT ************************
	///******************************************************
	
	std::ofstream	fOut("data/solution.txt");
	for (int y = -params.offsetY; y < params.by + params.offsetY; ++y) {
		for (int x = -params.offsetX; x < params.bx + params.offsetX; ++x) {
			fOut << params.getXCoord(x, y) << "\t" << params.getYCoord(x, y) << "\t" << u(x, y) << std::endl;
		}
		fOut << std::endl;
	}
	fOut.close();
};
