#ifndef	UTIL_HPP_INCLUDED
#define	UTIL_HPP_INCLUDED

#include	<stdlib.h>
#include	<iostream>
#include	<sstream>
#include	<fstream>
#include	<cmath>

/**
  Converts a string to an arbitrary type. >> operator must be defined for the target type.
  @param string string which should be converted
  @return converted string
  **/
template<typename T>
T StringTo(const std::string& string) {
	T valor;

	std::stringstream stream(string);
	stream >> valor;
	return valor;
}

constexpr	double f(const double x, const double y) {
	return 4.0*M_PI*M_PI*sin(2.0*M_PI*x)*sinh(2.0*M_PI*y);
}

constexpr	double border(const double x, const double y) {
	return sin(2.0*M_PI*x)*sinh(2.0*M_PI*y);
}

class	Params{
public:
	static constexpr double k = 2.0 * M_PI;
	
	int nx      = 0;
	int	ny      = 0;
	int	c       = 0;
	double eps  = 0.0;
	double eps2 = 0.0;

	double	hx     = 0.0;
	double	hy     = 0.0;

	double	invHx2 = 0.0;
	double	invHy2 = 0.0;

	double	preF   = 0.0;
	
	//size of block
	int	bx	= 0;
	int	by	= 0;
	int offsetX = 0;
	int offsetY = 0;
		
	Params(int argc, char **argv){
		if (argc != 5) {
			std::cout << "Invalid number of arguments!" << std::endl;
			std::cout << "./rgbs nx ny c eps" << std::endl;
			exit(EXIT_FAILURE);
		}
		nx          = StringTo<int>(argv[1]);
		ny          = StringTo<int>(argv[2]);
		c           = StringTo<int>(argv[3]);
		eps         = StringTo<double>(argv[4]);
		eps2 		= copysign(eps * eps, eps);

		// output configuration parameters
		std::cout << "nx," << nx << std::endl;
		std::cout << "ny," << ny << std::endl;
		std::cout << "c," << c << std::endl;
		std::cout << "eps," << eps << std::endl;
		std::cout << "eps2," << eps2 << std::endl;
		
		// calculate global parameters
		hx     = 2.0 / nx;
		hy     = 1.0 / ny;

		invHx2 = 1.0 / hx / hx;
		invHy2 = 1.0 / hy / hy;

		preF   = ( 2 * invHx2 + 2 * invHy2 + k*k );
		
		//make real number of points
		++nx;
		++ny;
	}
	
	inline
	double	getXCoord(const int indexX, const int indexY){
		(void)(indexY);
		return	hx * (indexX + offsetX);
	}
	
	inline
	double	getYCoord(const int indexX, const int indexY){
		(void)(indexX);
		return	hy * (indexY + offsetY);
	}

};


	
#endif