#include <iostream>
#include <numeric>
#include <fstream>
#include <algorithm>
#include "immintrin.h"
#include <memory>
#include <vector>
#include <mpi.h>
#include "GridCaptain.h"
#include	"Timer.h"
#include	<cmath>
#include <functional>

# define k 2.0*M_PI

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

    r = _mm256_setzero_pd();

    for(int i = 0 ; i< (int) vec.size(); i+=4)
    {
        a = _mm256_load_pd( &vec[i]);
        r = _mm256_add_pd(_mm256_mul_pd(a,a),r);
    }

    return (double) r[0]+r[1]+r[2]+r[3];

}

inline std::vector<double> matMult( std::vector<double> vec,GridCaptain gcap,const double alpha, const double beta, const double gama,
 const int * dim ,int size, int rank)
{  
   int veclen = (int)vec.size();
   vector<double> fresult(veclen,0);
   
   std::cout << "1 " << "\n";
   
   
    int *rec_cnt = new int[size];
   int *rec_disp = new int[size];
    int bleny =  dim[1];  
   //  int proc = size;  
    
    int blenx = gcap.worksheet[rank*3+1];
    int sx = gcap.worksheet[rank*3+0];
    rec_disp[rank] =  gcap.worksheet[rank*3+2];
    int sy = 0;
    //int ex = sx+blenx-1;
   // int ey = bleny;
    int sz=blenx*bleny;
    double * result = new double[sz];
    int gridno = 0;


            double gama1 = 0;
            double gama2 = 0;
            double beta1 = 0;
            double beta2 = 0;


    int rank1 =  rank+1;  
    int rank2 =  rank-1; 

    int rn = MPI_Cart_shift(MPI_COMM_WORLD,0,1,&rank,&rank1 );
    int rs = MPI_Cart_shift(MPI_COMM_WORLD,0,-1,&rank,&rank2);

    for(int i=sx; i<blenx ; i++)
    {
       // int l = (i-sx)%bleny;
        for(int j=sy; j<bleny ; j++)
        {
            
            gama1 = 0;
            gama2 = 0;
            beta1 = 0;
            beta2 = 0;
            gridno= i*dim[1] + j;            
            //int kl = (j-sy)%blenx;
            if(rn != MPI_PROC_NULL)
             gama1 = gama;
                        
            if(rs != MPI_PROC_NULL)
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

    MPI_Allgatherv(result,sz, MPI_DOUBLE, (void*)&fresult, rec_cnt,rec_disp, MPI_DOUBLE,MPI_COMM_WORLD );
   // MPI_Finalize();
    
    if(rank == 0)
    return fresult;
    else
    return 0;
}


inline std::vector<double> cal_fVec(GridCaptain gcap, const double gama, const int * dim, int size, int rank, int tgrdpoint, double hx, double hy)
{ 
   vector<double> fresult(tgrdpoint,0);

    int *rec_cnt = new int[size];
   int *rec_disp = new int[size];

    int bleny =  dim[1];    
    
    int blenx = gcap.worksheet[rank*3+1];
    int sx = gcap.worksheet[rank*3+0];
    rec_disp[rank] =  gcap.worksheet[rank*3+2];
    int sy = 0;
    int sz=blenx*bleny;
    double * result = new double[sz];
    int gridno = 0;
    


    int rank2 =  rank-1; 

    int rs = MPI_Cart_shift(MPI_COMM_WORLD,0,-1,&rank,&rank2);
   
     double gama2 = 0;

    int x = 0;
    int y = 0;

    for(int i=sx; i<blenx ; i++)
    {
        for(int j=sy; j<bleny ; j++)
        {
            
            gama2 = 0;
            gridno= i*dim[1] + j;            
            //int k = (j-sy)%blenx;
            x = (((gridno-1)%dim[1])+1)*hx;
            y = (((gridno-1)/dim[1])+1)*hy;
            int f = fxy(x,y);
                                  
            if(rs == MPI_PROC_NULL)
            {
             gama2 = gama*border(x,y);
             result[gridno] = f-gama2;
            }
            else 
                result[gridno] = f;
        }
    }

    MPI_Allgatherv(result,sz, MPI_DOUBLE,  (void*)&fresult, rec_cnt,rec_disp, MPI_DOUBLE,MPI_COMM_WORLD );
    //MPI_Finalize();
    if(rank == 0)
    return fresult;
    else
    return 0;
}


int main(int argc, char** argv)
{

    int size(0); // The total number of processes
   int rank(0); // The rank/number of this process (within MPI_COMM_WORLD)
  // int proc = gcap.proc;
   //int veclen = (int)vec.size();
   
   int nx = 0;
    int ny = 0;
    int iter = 0;
    //int proc = 0;
    int error=0;
 
  for (int i = 0; i<argc; i++)
        std::cout << argv[i] << "\n";
        
    if (argc != 5)
    {
        for (int i = 0; i<5; i++)
        std::cout << argv[i] << "\n";
        
        std::cout << argc << "= invalid number of argument.. Program exiting..";
        exit(EXIT_FAILURE);
        
        
    }

    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
    iter = atoi(argv[3]);
    error = atoi(argv[4]);
   // proc = atoi(argv[3]);

    int nnx = nx-1;
    int nny = ny-1;

    int totdim = nnx*nny;
    
     //std::cout << "111 " << "\n";

    FdGrid* fgrid = new FdGrid (nnx,nny);
    
    std::cout << "nx," << nx << std::endl;
	std::cout << "ny," << ny << std::endl;
	std::cout << "c," << iter <<std::endl;
	
	int gridpoint = fgrid->totalGridPoints();
	double hx = fgrid->getHx();
    double hy = fgrid->getHy();
   	
    ///******************************************************
	///********************** CALCULATION *******************
	///******************************************************
	double time = 0;
	double dt0 = 0.0;
	double resd = 0.0;
#ifdef USE_LIKWID
	likwid_markerInit();
	likwid_markerStartRegion("dummy");
#endif

	siwir::Timer	timer;

    //std::vector<double> xsol = callCG(*fGrid,iter,error);
        
    int len = 0;

    if(gridpoint%4 != 0)
    len = gridpoint + (4-gridpoint%4);
    std::vector<double> Xvec (len,0);
    std::vector<double> Rvec (len,0);
    std::vector<double> Fvec (len,0);
    std::vector<double> Tvec (len,0);
    std::vector<double> Tmpvec (len,0);
    int * dim = new int [2];         
    GridCaptain* gcap = NULL;
    //double * result=NULL;
    double alpha = 0;
      
    double alfa=0;
    double bita=0;
    double gama=0;
            
            bita = 1/hx/hx;
            gama = 1/hy/hy;
            alfa = -(2/gama+ 2/bita + k * k);
            
 // Initialization of MPI
   // ----------------------------------------------------------------
   MPI_Init(&argc, &argv);
   // ----------------------------------------------------------------

   // Determining the number of CPUs and the rank of this process
   // ----------------------------------------------------------------
   MPI_Comm_size( MPI_COMM_WORLD, &size );
   MPI_Comm_rank( MPI_COMM_WORLD, &rank );
   // ----------------------------------------------------------------   

   if (rank == 0)
   {  
    gcap = new GridCaptain(size,*fgrid);
    dim[0]=fgrid->getDimM();
    dim[1]=fgrid->getDimN();   
    MPI_Bcast(dim,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(gcap,1,MPI_INT,0,MPI_COMM_WORLD);
    // std::cout << "3 " << "\n";
   }
     
    //GridCaptain* gcap = new GridCaptain(proc,fgrid);
    
    Tvec = matMult(Xvec,*gcap, alfa, bita,gama, dim, size, rank);    
       
    Fvec = cal_fVec(*gcap,gama, size, rank, gridpoint,hx ,hy);
    
    if(rank==0)
    {

    std::transform (Fvec.begin(), Fvec.end(), Tvec.begin(), Rvec.begin(),  std::minus<double>());

    dt0 = std::inner_product(Rvec.begin(), Rvec.end(), Rvec.begin(),0);

    resd = compute2norm(Rvec);
      std::vector<double> Dvec (Rvec);
    }

    if(resd > error)
     {            
        for(int i = 0 ; i<iter; i++)
       {
        
        MPI_Bcast(Dvec,1,MPI_INT,0,MPI_COMM_WORLD);
        
        Tvec = matMult(Dvec,*gcap, alfa, bita, gama, size, rank);
        
        if(rank == 0)
        {
        double dt = std::inner_product(Dvec.begin(), Dvec.end(), Tvec.begin(),0);

        alpha = dt0 / dt;
              
        for(int j=0; j< (int)Dvec.size();j+=4)
        {
            Tmpvec[j] = alpha * Dvec[j];
             Tmpvec[j+1] = alpha * Dvec[j+1];
              Tmpvec[j+2] = alpha * Dvec[j+2];
               Tmpvec[j+3] = alpha * Dvec[j+3]; 
            
        }     

        std::transform (Xvec.begin(), Xvec.end(), Tmpvec.begin(), Tvec.begin(),   std::plus<double>());
        
        for(int j=0; j< (int)Tvec.size();j+=4)
        {
           Tmpvec[j] = alpha * Tvec[j];
             Tmpvec[j+1] = alpha * Tvec[j+1];
              Tmpvec[j+2] = alpha * Tvec[j+2];
               Tmpvec[j+3] = alpha * Tvec[j+3]; 
            
        }      
        
        std::transform (Rvec.begin(), Rvec.end(), Tmpvec.begin(), Rvec.begin(),  std::minus<double>());

        double dt1 = std::inner_product(Rvec.begin(), Rvec.end(), Rvec.begin(),0);

        resd = compute2norm(Rvec);

        if(resd < error)
            break;

        double beta = dt1/dt0;

        //std::transform (Dvec.begin(), Dvec.end(), Tmpvec.begin(),  std::multiplies<double>(),beta);
        
         for(int j=0; j< (int)Dvec.size();j+=4)
        {
             Tmpvec[j] = beta * Dvec[j];
             Tmpvec[j+1] = beta * Dvec[j+1];
              Tmpvec[j+2] = beta * Dvec[j+2];
               Tmpvec[j+3] = beta * Dvec[j+3]; 
        }      
        
        std::transform (Rvec.begin(), Rvec.end(), Tmpvec.begin(), Tvec.begin(),   std::plus<double>());

        dt0 = dt1;
        }
    }

}   
    
    MPI_Finalize();
    	time = timer.elapsed();
	std::cout << "time," << time << std::endl;

#ifdef USE_LIKWID
	likwid_markerStopRegion("dummy");
	likwid_markerClose();
#endif

    for (int i= 0; i< totdim; i++ )
        std::cout << Xvec[i] << ' ';

    return 0;

}
