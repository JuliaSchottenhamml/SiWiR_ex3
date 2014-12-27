#include <iostream>
#include <numeric>
#include <fstream>
#include <algorithm>
#include "FdGrid.h"
#include "immintrin.h"
#include <memory>
#include <vector>
#include <mpi.h>
#include "GridCaptain.h"
#include	"Timer.h"
#include	<cmath>

# define k = 2.0*M_PI

#define ERRLIMIT 0;

inline double fxy(const int x, const int y){

         return 4.0*M_PI*M_PI*sin(2.0*M_PI*x)*sinh(2.0*M_PI*y);
    }

inline double border(const double x, const double y){
        return sin(2.0*M_PI*x)*sinh(2.0*M_PI*y);
    }

inline double compute2norm(std::vector<double> vec)
{

    __m256d a, r;

    r = _mm256_setzero_ps();

    for(int i = 0 ; i<vec.size(); i+=4)
    {
        a = _mm256_load_pd(vec+i);
        r = _mm256_add_ps(_mm256_mul_ps(a,a),r);

    }

    return r[0]+r[1]+r[2]+r[3];

}

inline std::vector matMult(FdGrid& fgrid, vector<double> vec, GridCaptain& gcap,const double alpha, const double bita, const double gamma)
{

   int size(0); // The total number of processes
   int rank(0); // The rank/number of this process (within MPI_COMM_WORLD)
   int proc = gcap.proc;
   int veclen = vec.size();
   vector<double> fresult(veclen,0);
   int rec_cnt = new int[proc];
   int rec_disp = new int[proc];
   
   // Initialization of MPI
   // ----------------------------------------------------------------
   MPI_Init();
   // ----------------------------------------------------------------

   // Determining the number of CPUs and the rank of this process
   // ----------------------------------------------------------------
   MPI_Comm_size( MPI_COMM_WORLD, &size );
   MPI_Comm_rank( MPI_COMM_WORLD, &rank );
   // ----------------------------------------------------------------   

   //int MPI_Cart_shift(MPI_COMM_WORLD,0,1,rank,rank+2);
 
    int * dim = new int [2];
   if (rank == 0)
   {  

    dim[0]=fgrid.getDimM();
    dim[1]=fgrid.getDimN();
   
    MPI_Bcast(dim,1,MPI_INT,0,MPI_COMM_WORLD);
    
   }
       int eval=0, wval = 0, sval = 0, nval = 0, cval = 0;
    int bleny =  dim[1];    
    
    int blenx = gcap.worksheet[rank*3+1];
    int sx = gcap.worksheet[rank*3+0];
    rec_disp[rank] =  gcap.worksheet[rank*3+2];
    int sy = 0;
    int ex = sx+blenx-1;
    int ey = bleny;
    int sz=blenx*bleny;
    array[rank]=sz;
    double * result = new double[sz];
    int gridno = 0;

    double *ss = new double[blenx];
    double *ns = new double[blenx];
            double gama1 = 0;
            double gama2 = 0;
            double beta1 = 0;
            double beta2 = 0;

   
    double *sr = new double[blenx];
    double *nr = new double[blenx];
                 
    for ( int i=0; i <blenx;i++)
    {
        sr[i]=0;
        er[i]=0;
    }

    ss=fgrid.getDataAdd(sx,sy);
    ns=fgrid.getDataAdd(ex,sy);    

    int rn = MPI_Cart_shift(MPI_COMM_WORLD,0,1,rank,rank+1 );
    int rs = MPI_Cart_shift(MPI_COMM_WORLD,0,-1,rank,rank-1);

    for(int i=sx; i<blenx ; i++)
    {
        int l = (i-sx)%bleny;
        for(int j=sy; j<bleny ; j++)
        {
            
            gama1 = 0;
            gama2 = 0;
            beta1 = 0;
            beta2 = 0;
            gridno= i*dim[1] + j;            
            int k = (j-sy)%blenx;
            if(rn != MPI_NULL_PROC)
             gama1 = gama;
                        
            if(rs != MPI_NULL_PROC)
             gama2 = gama;
                                                           
            if(j!=sy)
                beta1=beta;
            
            if(j != bleny-1)
                beta2=beta;
            
            int pm = 0;
            int ppm =0;
            int fm =0;
            int ffm =0;  
             
            int cm = alpha*vec[gridno];
            if(gridno-1>=0)
            pm = beta1*vec[gridno-1];
            if(gridno-3>=0)
            ppm = gama1*vec[gridno-3];
            if(gridno+1<sz)
            fm = beta2*vec[gridno+1];
            if(gridno+3<sz)
            ffm = gama2*vec[gridno+3];
            result[gridno]=cm+pm+ppm+fm+ffm;

        }
    }

    MPI_Allgatherv(result,sz, MPI_DOUBLE, fresult, rec_cnt,rec_disp, MPI_DOUBLE,MPI_COMM_WORLD );
    MPI_Finalize();
    
    return fresult;
}


inline std::vector<double> cal_fVec(FdGrid& fgrid, GridCaptain& gcap)
{

   int size(0); // The total number of processes
   int rank(0); // The rank/number of this process (within MPI_COMM_WORLD)
   int proc = gcap.proc;
   int veclen = vec.size();
   vector<double> fresult(veclen,0);
   int rec_cnt = new int[proc];
   int rec_disp = new int[proc];
   
   // Initialization of MPI
   // ----------------------------------------------------------------
   MPI_Init();
   // ----------------------------------------------------------------

   // Determining the number of CPUs and the rank of this process
   // ----------------------------------------------------------------
   MPI_Comm_size( MPI_COMM_WORLD, &size );
   MPI_Comm_rank( MPI_COMM_WORLD, &rank );
   // ----------------------------------------------------------------   

   //int MPI_Cart_shift(MPI_COMM_WORLD,0,1,rank,rank+2);
   
    int * dim = new int [2];
   if (rank == 0)
   {  

    dim[0]=fgrid.getDimM();
    dim[1]=fgrid.getDimN();
   
    MPI_Bcast(dim,1,MPI_INT,0,MPI_COMM_WORLD);
    
   }
       int eval=0, wval = 0, sval = 0, nval = 0, cval = 0;
    int bleny =  dim[1];    
    
    int blenx = gcap.worksheet[rank*3+1];
    int sx = gcap.worksheet[rank*3+0];
    rec_disp[rank] =  gcap.worksheet[rank*3+2];
    int sy = 0;
    int ex = sx+blenx-1;
    int ey = bleny;
    int sz=blenx*bleny;
    array[rank]=sz;
    double * result = new double[sz];
    int gridno = 0;

    double *ss = new double[blenx];
    double *ns = new double[blenx];
    
    double *sr = new double[blenx];
    double *nr = new double[blenx];
                 
    for ( int i=0; i <blenx;i++)
    {
        sr[i]=0;
        er[i]=0;
    }

    ss=fgrid.getDataAdd(sx,sy);
    ns=fgrid.getDataAdd(ex,sy);
    
    int rn = MPI_Cart_shift(MPI_COMM_WORLD,0,1,rank,rank+1 );
    int rs = MPI_Cart_shift(MPI_COMM_WORLD,0,-1,rank,rank-1);

    /* MPI_ISend(ss,1,MPI_double,rs, MPI_COMM_WORLD);
     MPI_ISend(ns,1,MPI_double,rn, MPI_COMM_WORLD);

     MPI_recv(sr,1,MPI_double,rs, MPI_COMM_WORLD);
     MPI_recv(nr,1,MPI_double,rn, MPI_COMM_WORLD);*/
    int x = 0;
    int y = 0;

    for(int i=sx; i<blenx ; i++)
    {
        for(int j=sy; j<bleny ; j++)
        {
            
            gama1 = 0;
            gama2 = 0;
            beta1 = 0;
            beta2 = 0;
            gridno= i*dim[1] + j;            
            int k = (j-sy)%blenx;
            x = (((gridno-1)%dimN)+1)*hx;
            y = (((gridno-1)/dimN)+1)*hy;
            int f = fxy(x,y);
                                  
            if(rs == MPI_NULL_PROC)
            {
             gama2 = gama*border(x,y);
             result[gridno] = f-gama2;
            }
            else 
                result[gridno] = f;
        }
    }

    MPI_Allgatherv(result,sz, MPI_DOUBLE, fresult, rec_cnt,rec_disp, MPI_DOUBLE,MPI_COMM_WORLD );
    MPI_Finalize();
    
    return fresult; 
}


inline double * callCG(FdGrid& fgrid, int const iter, int const proc, int const err)
{

    std::vector<double> Xvec (fgrid.totalGridPoints(),0);
    std::vector<double> Rvec (fgrid.totalGridPoints(),0);
    std::vector<double> Fvec (fgrid.totalGridPoints(),0);
    std::vector<double> Tvec (fgrid.totalGridPoints(),0);
    std::vector<double> Tmpvec (fgrid.totalGridPoints(),0);
    double alpha = 0;
    double hx = fgrid.hx;
    double hy = fgrid.hx;
    
    double alfa=0;
        double bita=0;
            double gama=0;
            
            bita = 1/hx/hx;
            gama = 1/hy/hy;
            alfa = -(2/gama+ 2/bita + k * k);

    GridCaptain gcap = new GridCaptain(fgrid);
    
    Tvec = initMatMult(fgrid,Xvec,gcap,alfa, bita, gama);
    
    Fvec = cal_fVec(fgrid,gcap);

    std::transform (Fvec.begin(), Fvec.end(), TVec.begin(), Rvec.begin(),  std::minus<double>());

    double dt0 = std::inner_product(Rvec.begin(), Rvec.end(), Rvec.begin(),0);

    double resd = compute2norm(Rvec);

    if(resd < err)
        return Xvec;

      std::vector<double> Dvec (Rvec);

    for(int i = 0 ; i<iter; i++)

    {
        Tvec = matMult(fgrid,Dvec,gcap,alfa, bita, gama);

        double dt = std::inner_product(Dvec.begin(), Dvec.end(), Tvec.begin(),0);

        alpha = dt0 / dt;

        std::transform (Dvec.begin(), Dvec.end(), Tmpvec.begin(),  std::multiplies<double>(),alpha);

        std::transform (Xvec.begin(), Xvec.end(), Tmpvec.begin(), Tvec.begin(),   std::plus<double>());

        std::transform (Tvec.begin(), Tvec.end(), Tmpvec.begin(),  std::multiplies<double>(),alpha);

        std::transform (Rvec.begin(), Rvec.end(), Tmpvec.begin(), Rvec.begin(),  std::minus<double>());

        double dt1 = std::inner_product(Rvec.begin(), Rvec.end(), Rvec.begin(),0);

        resd = compute2norm(Rvec);

        if(resd < err)
            return Xvec;

        double beta = dt1/dt0;

        std::transform (Dvec.begin(), Dvec.end(), Tmpvec.begin(),  std::multiplies<double>(),beta);

        std::transform (Rvec.begin(), Rvec.end(), Tmpvec.begin(), TVec.begin(),   std::plus<double>());

        dt0 = dt1;

    }

    return Xvec;

}


int main(int argc, char** argv)
{

    int nx = 0;
    int ny = 0;
    int iter = 0;
    int proc = 0;
    int error=0;

    if (argc != 6)
    {
        std::cout << "invalid number of argument.. Program exiting..";
        exit(EXIT_FAILURE);
    }

    nx = atoi(argv[5]);
    ny = atoi(argv[6]);
    iter = atoi(argv[7]);
    error = atoi(argv[8]);
    proc = atoi(argv[3]);

    int nnx = nx-1;
    int nny = ny-1;

    int totdim = nnx*nny;

    FdGrid *fGrid = new FdGrid (nnx,nny);
    
    std::cout << "nx," << nx << std::endl;
	std::cout << "ny," << ny << std::endl;
	std::cout << "c," << iter <<std::endl;
	
    ///******************************************************
	///********************** CALCULATION *******************
	///******************************************************
	double time = 0;
	
#ifdef USE_LIKWID
	likwid_markerInit();
	likwid_markerStartRegion("dummy");
#endif

	siwir::Timer	timer;

    std::vector<double> xsol = callCG(fGrid,iter,error,proc);
    
    	time = timer.elapsed();
	std::cout << "time," << time << std::endl;

#ifdef USE_LIKWID
	likwid_markerStopRegion("dummy");
	likwid_markerClose();
#endif

    for (int i= 0; i< totdim; i++ )
        std::cout << xsol[i] << ' ';

    return 0;

}
