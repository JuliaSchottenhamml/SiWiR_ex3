#ifndef FdGrid_H
#define FdGrid_H
 
#endif // FdGrid_H


#include <iostream>
#include "immintrin.h"
#include <memory>
#define LD 64
#define domxl 0
#define domxh 2
#define domyl 0
#define domyh 1


using namespace std;


class FdGrid{

    int dimM, dimN;
    double hx,hy;
    int totdim;

    //__declspec(align(16,sizeof(double)))/////
    double * matarray, *bx;

public:

    FdGrid()
    {
        dimM = 0;
        dimN = 0;
        matarray = new double[0];

    }

    FdGrid(int m, int n)
    {
        dimM = m;
        dimN = n;
        hx = (double)((domxh-domxl)/(n+1.0));
        hy = (double)((domyh-domyl)/(m+1.0));
        
        totdim = (dimM+LD) * (dimN+LD);
        matarray = new double[totdim];
       for (int i=0; i < totdim ; i++)
       {
         matarray [i] = 0;
       }

    }

    FdGrid(const int m, const int n, const double * data)
    {
        dimM = n;
        dimN = m;

        hx = (domxh-domxl)/(n+1);
        hy = (domyh-domyl)/(m+1);
        
        std::cout << "..in FGrid Hx and Hy.." << hx << "  " << hy << "\n";

        int newdimN = dimN + LD;
        int newdimM = dimM + LD;
        matarray = new double[(dimM+LD) * (dimN+LD)];
        bx = new double[dimM];

        for(int i=0; i<newdimM; i++)
        {
            for(int j=0; j<newdimN; j++)
            {
                if (i< dimM && j < dimN)
                matarray [i*newdimN +j] = data[i*(dimN+LD) +j];
                else
                matarray [i*newdimN +j] = 0;

            }

        }
    }

    inline void initialize()
    {
        for (int i=0; i < totdim ; i++)
        {
            /*if(i>= totdim - dimN)
            {
                double x = (dimN  - i +totdim)*hx;
                double y = domyh;

                matarray [i]= border(x,y);
            }
           else*/
                  matarray [i] = 0;
        }
    }

    int totalGridPoints()
    {
        return dimN*dimM;
    }

    inline double fxy(const int i, const int j){

        double x = j*hx;
        double y = i*hy;

        return 4.0*M_PI*M_PI*sin(2.0*M_PI*x)*sinh(2.0*M_PI*y);
    }

    inline  double border(const double x, const double y){
        return sin(2.0*M_PI*x)*sinh(2.0*M_PI*y);
    }


    void storeDataA(int m, int n, double data)
    {
        matarray [m * (dimN + LD) + n] = data;
    }

    double getDataA(int m, int n)
    {
        return matarray[m*(dimN+LD)+n];
    }

    double * getDataAdd(int m, int n)
    {
        return &matarray[m*(dimN+LD)+n];
    }


    void storeDataB(int m ,double data)
    {
        bx[m] = data;
    }

    double getDataB(int m)
    {
        return bx[m];
    }

   /* void storeAVXData(int m, int n, __m256d data)
    {

        _mm256_store_pd(&matarray[m*(dimN+LD)+n], data);

    }

    __m256d getAVXData(int m, int n)
    {
        __m256d m1 = _mm256_load_pd(&matarray[m*dimN+n]);

        return m1;
    }*/

    int getDimM()
    {
        return dimM;
    }

    int getDimN()
    {
        return dimN;
    }
    
    
    double getHx()
    {
        return hx;
    }
    

    double getHy()
    {
        return hy;
    }
    
    
    ~FdGrid()
    {
        delete matarray;
    }

};
