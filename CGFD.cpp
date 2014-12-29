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

inline double compute2norm(double * vec)
{

    __m256d a, r;

    r = _mm256_setzero_pd();

    for(int i = 0 ; i< (int)sizeof(vec); i+=4)
    {
        a = _mm256_load_pd( &vec[i]);
        r = _mm256_add_pd(_mm256_mul_pd(a,a),r);
    }

    return (double) r[0]+r[1]+r[2]+r[3];

}

inline double compute2normVec(vector<double> vec)
{

    __m256d a, r;

    r = _mm256_setzero_pd();

    for(int i = 0 ; i< (int)vec.size(); i+=4)
    {
        a = _mm256_load_pd( &vec[i]);
        r = _mm256_add_pd(_mm256_mul_pd(a,a),r);
    }

    return (double) r[0]+r[1]+r[2]+r[3];

}

inline double * matMult( std::vector<double> vec,int blenx,int bleny,int sx,const double alpha, const double beta, const double gama,
   int destn, int dests)
{  
   
    int sy = 0;
    int sz=blenx*bleny;
    double * result = new double[sz];
    int gridno = 0;
    double gama1 = 0;
    double gama2 = 0;
    double beta1 = 0;
    double beta2 = 0;

  for(int i=sx; i<blenx ; i++)
    {
       // int l = (i-sx)%bleny;
        for(int j=sy; j<bleny ; j++)
        {
            gama1 = 0;
            gama2 = 0;
            beta1 = 0;
            beta2 = 0;
            gridno= i*bleny + j;            
            //int kl = (j-sy)%blenx;
            if(destn != -1)
             gama1 = gama;
                        
            if(dests != -1)
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

    return result;
    
}

inline double * cal_fVec(int blenx,int bleny ,int sx,const double gama,  double hx, double hy, int dests)
{ 
  // vector<double> fresult(tgrdpoint,0);

    //int bleny =  dim[1];    
    
       
    int sy = 0;
    int sz=blenx*bleny;
    double * result = new double[sz];
    int gridno = 0;
     
     double gama2 = 0;

    int x = 0;
    int y = 0;

    for(int i=sx; i<blenx ; i++)
    {
        for(int j=sy; j<bleny ; j++)
        {
            
            gama2 = 0;
            gridno= i*bleny + j;            
            //int k = (j-sy)%blenx;
            x = (((gridno-1)%bleny)+1)*hx;
            y = (((gridno-1)/bleny)+1)*hy;
            int f = fxy(x,y);
                                  
            if(dests == -1)
            {
             gama2 = gama*border(x,y);
             result[gridno] = f-gama2;
            }
            else 
                result[gridno] = f;
        }
    }

       return result;
   
}


int main(int argc, char** argv)
{

            
    if (argc < 5)
    {
        for (int i = 0; i<5; i++)
        std::cout << argv[i] << "\n";
        
        std::cout << argc << "= invalid number of argument.. Program exiting..";
        exit(EXIT_FAILURE);
    
    }
    
    int size(0); // The total number of processes
    int rank(0); // The rank/number of this process (within MPI_COMM_WORLD)
    int nx = 0;
    int ny = 0;
    int iter = 0;
    int error=0;
    FdGrid* fgrid;
   	double time = 0;
	double  dt0[1];
 
   
	//double  iresd;
#ifdef USE_LIKWID
	likwid_markerInit();
	likwid_markerStartRegion("dummy");
#endif

	siwir::Timer	timer;

    //std::vector<double> xsol = callCG(*fGrid,iter,error);
       // Initialization of MPI
   // ----------------------------------------------------------------
   MPI_Init(&argc, &argv);
   // ----------------------------------------------------------------

   // Determining the number of CPUs and the rank of this process
   // ----------------------------------------------------------------
   MPI_Comm_size( MPI_COMM_WORLD, &size );
   MPI_Comm_rank( MPI_COMM_WORLD, &rank );
   // ----------------------------------------------------------------   

         MPI_Request request;
         MPI_Status status;  
    

    double alfa=0;
    double bita=0;
    double gama=0;
    int gridpoint = 0;
    double hx = 0.0, hy=0.0;
    int len = 0;
    int dests=0, destn=0, blenx=0, sx =0; 
 
    int * dim = new int [2]; 
        
    GridCaptain* gcap = NULL ;
    
    double alpha = 0;
      int nnx =0, nny=0;
 
   if (rank == 0)
   { 
    std::cout << "3 " << "\n";
    
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
    nnx = nx-1;
    nny = ny-1;
    iter = atoi(argv[3]);
    error = atoi(argv[4]);
     fgrid = new FdGrid (nnx,nny); 
     gridpoint = fgrid->totalGridPoints();
	 hx = fgrid->getHx();
     hy = fgrid->getHy();
    bita = 1/hx/hx;
    gama = 1/hy/hy;
    alfa = -(2/gama+ 2/bita + k * k);
    gcap = new GridCaptain(size,*fgrid);
    
    for(int t=0;t<size;t++)
    {
    blenx = gcap->worksheet[t*3+1];
    sx = gcap->worksheet[t*3+0];
    MPI_Isend(&blenx,1,MPI_INT,t,t+100,MPI_COMM_WORLD,&request);
    MPI_Isend(&sx,1,MPI_INT,t,t+100,MPI_COMM_WORLD,&request);
    }
    dim[0]=fgrid->getDimM();
    dim[1]=fgrid->getDimN();
    MPI_Bcast(&nnx,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&nny,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&nx,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&ny,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&iter,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&error,1,MPI_INT,0,MPI_COMM_WORLD);   
    MPI_Bcast(dim,2,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(gcap,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&gridpoint,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&hx,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&hy,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&alfa,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&bita,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&gama,1,MPI_INT,0,MPI_COMM_WORLD);
       
      std::cout << "nx," << nx << std::endl;
	std::cout << "ny," << ny << std::endl;
	std::cout << "c," << iter <<std::endl;
   }
   
    int totdim = nnx*nny;
    MPI_Recv(&blenx,1, MPI_INT,0, rank+100, MPI_COMM_WORLD,&status);
    MPI_Recv(&sx,1, MPI_INT,0, rank+100, MPI_COMM_WORLD,&status);
    
   MPI_Barrier(MPI_COMM_WORLD);
   
   if(gridpoint%4 != 0)
    len = gridpoint + (4-gridpoint%4);
         
    std::vector<double> Xvec (len,0);
    std::vector<double> Rvec (len,0);
    std::vector<double> Fvec (len,0);
    std::vector<double> Tvec (len,0);
    std::vector<double> Tmpvec (len,0); 
 
    int bleny =  dim[1];  
           
    int sz=blenx*bleny;
    
    double * tresult = new double[sz];
    double * fresult = new double[sz];
    double * mresult = new double[sz];
    double * nresult = new double[sz];
    double resd =0.0;
   
    if(rank == 0)
    destn = -1;
    else 
    destn = rank -1;
    
    if(rank == size)
    dests = -1;
    else 
    destn = rank +1;    
    
    tresult = matMult(Xvec,blenx,bleny,sx, alfa, bita,gama, destn,dests);        
     
    fresult = cal_fVec(blenx,bleny,sx,gama, hx ,hy,dests);
    std::cout << "3222 " << sizeof(mresult) << " " << sizeof(fresult) << " "<< sizeof(tresult) << "\n";     
    for(int i = 0; i< (int)sizeof(tresult); i++)
    {
        mresult[i] = fresult[i]-tresult[i];
        std::cout << "\n" << rank << " " << mresult[i];
    } 
    std::cout << "3333 " << "\n"; 
     double resdlocal=0.0; 
      std::cout << "1### " << resdlocal;
    //  std::cout << "2### " <<  ;
    resdlocal = compute2norm(mresult);
    std::cout << "2### " << resdlocal ;
    MPI_Reduce(&resdlocal, dt0,1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
    std::cout << "3### " << "\n";
    MPI_Isend(mresult,(int)sizeof(mresult), MPI_DOUBLE, 0, rank, MPI_COMM_WORLD,&request);
    std::cout << "4### " << "\n";
    int jn=0;
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
    {
     for( int i=0; i< size; i++)
    {
      MPI_Recv(nresult,(int)sizeof(mresult), MPI_DOUBLE, i, i, MPI_COMM_WORLD,&status);
      std::cout << "5### " << "\n";
      for(int l=0; l< (int)sizeof(nresult);l++)
          Rvec[++jn]= nresult[l];
    }         
    }     
    
    if(*dt0 > error)
     {       
        std::vector<double> Dvec (Rvec);     
        for(int i = 0 ; i<iter; i++)
       {
        
        MPI_Bcast((void*)&Dvec,1,MPI_INT,0,MPI_COMM_WORLD);
        
        tresult = matMult(Dvec, blenx,bleny,sx, alfa, bita, gama, destn,dests);
        
       // MPI_Allgatherv((void*)tresult,rec_cnt[rank], MPI_DOUBLE, (void*)&Tvec, rec_cnt,rec_disp, MPI_DOUBLE,MPI_COMM_WORLD );
        MPI_Isend(tresult,(int)sizeof(tresult), MPI_DOUBLE, 0, rank+10, MPI_COMM_WORLD,&request);
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank == 0)
        {
            int jk = 0;
            for( int km=0; km < size; km++)
            {
              MPI_Recv(nresult,(int)sizeof(tresult), MPI_DOUBLE,km, km+10, MPI_COMM_WORLD,&status);
             // std::cout << "5### " << "\n";
              for(int l=0; l< (int)sizeof(nresult);l++)
                  Rvec[++jk]= nresult[l];
             }  
        double dt = std::inner_product(Dvec.begin(), Dvec.end(), Tvec.begin(),0);

        alpha = *dt0 / dt;
              
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

        resd = compute2normVec(Rvec);

        if(resd < error)
            break;

        double beta = dt1/(*dt0);

        //std::transform (Dvec.begin(), Dvec.end(), Tmpvec.begin(),  std::multiplies<double>(),beta);
        
         for(int j=0; j< (int)Dvec.size();j+=4)
        {
             Tmpvec[j] = beta * Dvec[j];
             Tmpvec[j+1] = beta * Dvec[j+1];
              Tmpvec[j+2] = beta * Dvec[j+2];
               Tmpvec[j+3] = beta * Dvec[j+3]; 
        }      
        
        std::transform (Rvec.begin(), Rvec.end(), Tmpvec.begin(), Tvec.begin(),   std::plus<double>());

        *dt0 = dt1;
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
