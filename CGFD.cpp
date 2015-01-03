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

inline double fxy(const double x, const double y){

         return 4.0*M_PI*M_PI*sin(2.0*M_PI*x)*sinh(2.0*M_PI*y);
    }

inline double border(const double x, const double y){
        return sin(2.0*M_PI*x)*sinh(2.0*M_PI*y);
    }

inline double * matMult( std::vector<double> vec,int blenx,int bleny,int sx,const double alpha, const double beta, const double gama,
   int destn, int dests)
{  
    
    //std::cout << " in cal matmult ";
   
    int sy = 0;
    int len=0;
    int sz=blenx*bleny;
    if(sz%4 != 0)
    len = sz + (4-sz%4);
    else
    len = sz;
    double * result = new double[len];
    int gridno = 0;
    double gama1 = 0.0;
    double gama2 = 0.0;
    double beta1 = 0.0;
    double beta2 = 0.0;
    for(int h=0;h<len;h+=4)
    {
      result[h]=0.0;  
      result[h+1]=0.0;
      result[h+2]=0.0;
      result[h+3]=0.0;
    }
    int l =0;

  for(int i=sx; i<sx+blenx ; i++)
      {
       // int l = (i-sx)%bleny;
        for(int j=sy; j<bleny ; j++)
        {
            gama1 = 0.0;
            gama2 = 0.0;
            beta1 = 0.0;
            beta2 = 0.0;
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
            result[l++]=cm+pm+ppm+fm+ffm;

        }
    }

    return result;
    
}

inline double * cal_fVec(int blenx,int bleny ,int sx,const double gama,  double hx, double hy, int dests)
{ 
  // vector<double> fresult(tgrdpoint,0);

    //int bleny =  dim[1];    
    //std::cout << rank << " in cal fvec ";
       
    int sy = 0;
    int len=0;
    int sz=blenx*bleny;
    if(sz%4 != 0)
    len = sz + (4-sz%4);
    else
    len = sz;
    double * result = new double[len];
    
    int gridno = 0;
     
    double gama2 = 0.0;

    double x = 0.0;
    double y = 0.0;

     for(int h=0;h<len;h+=4)
    {
      result[h]=0.0;  
      result[h+1]=0.0;
      result[h+2]=0.0;
      result[h+3]=0.0;
    }

    int l = 0;
    for(int i=sx; i< sx+blenx ; i++)
    {
        for(int j=sy; j<bleny ; j++)
        {
            
            gama2 = 0.0;
            gridno= i*bleny + j;            
            //int k = (j-sy)%blenx;
            x = (((gridno-1)%bleny)+1)*hx;
            y = (((gridno-1)/bleny)+1)*hy;
            //std::cout << "x, y" << x << " " <<y;
            double f = fxy(x,y);
                        //std::cout << rank << " x, y, f " << x << " " <<y << " " << f << "\n";                      
            if(dests == -1)
            {
             gama2 = gama*border(x,y);
             result[l++] = f-gama2;
            }
            else 
                result[l++] = f;
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
    

    double alfa=0.0;
    double bita=0.0;
    double gama=0.0;
    int gridpoint = 0;
    double hx = 0.0, hy=0.0;
    int len = 0;
    int dests=0, destn=0, blenx=0, sx =0;
     int nx = 0;
    int ny = 0;
    int iter = 0;
    double error=0.0;
    FdGrid* fgrid;
   	double time = 0;
	double  dt0[1]; 
 
    int * dim = new int [2]; 
        
    GridCaptain* gcap = NULL ;
    
    double alpha = 0.0;
      int nnx =0, nny=0;
 
   if (rank == 0)
   { 
    //std::cout << "3 " << "\n";
    
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
    nnx = nx-1;
    nny = ny-1;
    iter = atoi(argv[3]);
    error = atof(argv[4]);
     fgrid = new FdGrid (nnx,nny); 
     gridpoint = fgrid->totalGridPoints();
	 hx = fgrid->getHx();
     hy = fgrid->getHy();
    bita = 1/hx/hx;
    gama = 1/hy/hy;
    alfa = (-1.0)*(2.0*gama+ 2.0*bita + k * k);
    //std::cout << "..Rank 0.." << atof(argv[4]) << "  " << error << "\n";
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
   
       
   std::cout << "nx," << nx << std::endl;
    std::cout << "ny," << ny << std::endl;
	std::cout << "c," << iter <<std::endl;
   }
      
       MPI_Bcast(&nnx,1,MPI_INT,0,MPI_COMM_WORLD);
       MPI_Bcast(&nny,1,MPI_INT,0,MPI_COMM_WORLD);
       MPI_Bcast(&nx,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&ny,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&iter,1,MPI_INT,0,MPI_COMM_WORLD);
       MPI_Bcast(&error,1,MPI_DOUBLE,0,MPI_COMM_WORLD);   
       MPI_Bcast(dim,2,MPI_INT,0,MPI_COMM_WORLD);
       MPI_Bcast(&gridpoint,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&hx,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(&hy,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
       
    
    MPI_Bcast(&alfa,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&bita,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&gama,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Recv(&blenx,1, MPI_INT,0, rank+100, MPI_COMM_WORLD,&status);
    MPI_Recv(&sx,1, MPI_INT,0, rank+100, MPI_COMM_WORLD,&status);
   // std::cout << rank << " gridpoint =  " << gridpoint << "\n"; 
    MPI_Barrier(MPI_COMM_WORLD);
    int totdim = nnx*nny;
     
  /*  std::cout << rank << " nnx = " << nnx << " "; 
    std::cout << rank << " nny = " << nny << " "; 
     std::cout << rank << " nx = " << nx << " "; 
      std::cout << rank << " ny = " << ny << " "; 
      std::cout << rank << " iter = " << iter << " ";  
       std::cout << rank << " error = " << error << " ";              
        std::cout << rank << " gridpoint = " << gridpoint << " "; 
         std::cout << rank << " hx = " << hx << " "; 
         std::cout << rank << " hy = " << hy << " "; 
         std::cout << rank << " alpha = " << alfa << " "; 
         std::cout << rank << " beta = " << bita << " "; 
         std::cout << rank << " gama = " << gama << " "; 
          std::cout << rank << " blenx = " << blenx << " "; 
         std::cout << rank << " sx = " << sx << " ";*/
     
     
   if(gridpoint%4 != 0)
    len = gridpoint + (4-gridpoint%4);
   else 
   len = gridpoint;
   

    std::vector<double> Xvec (len,0);
    std::vector<double> Rvec (len,0);
    std::vector<double> Fvec (len,0);
    std::vector<double> Tvec (len,0);
    std::vector<double> Tmpvec (len,0); 
 
    int bleny =  dim[1];
    int broke=0;  
           
    int sz=blenx*bleny;
    
    if(sz%4 != 0)
    len = sz + (4-sz%4);
    else
    len = sz;
    
    double * tresult = new double[len];
    double * fresult = new double[len];
    double * mresult = new double[len];
    double * nresult = new double[len];
              
    //double resd =0.0;
   
    if(rank == 0)
    destn = -1;
    else 
    destn = rank -1;
    
    if(rank == size-1)
    dests = -1;
    else 
    dests = rank +1;    
    
       tresult = matMult(Xvec,blenx,bleny,sx, alfa, bita,gama, destn,dests);        
       fresult = cal_fVec(blenx,bleny,sx,gama, hx ,hy,dests);
   // std::cout << "\n" << rank << " " << blenx << " " << bleny << " " << sx << " " << gama << " " << hx << " " << hy << " " << dests;
      for(int i = 0; i< (int)sizeof(tresult); i++)
    {
       // std::cout << "\n" << rank << " " << fresult[i];
        mresult[i] = fresult[i]-tresult[i];
        std::cout << "\n" << rank << " " << (int)sizeof(tresult) << " " << fresult[i] << " " << tresult[i] << " " << mresult[i];
    } 
     
     double resdlocal=0.0; 
     
     
     
     for(int i = 0 ; i< (int)sizeof(mresult); i+=4)
    {   
        resdlocal += mresult[i] * mresult[i];
        resdlocal += mresult[i+1] * mresult[i+1];
        resdlocal += mresult[i+2] * mresult[i+2];
        resdlocal += mresult[i+3] * mresult[i+3];
    }
    
    std::cout << "\n %%%%%%%%%%%%%%%%%%%  resedual=  "<< << rank << " " << resdlocal;
    
    MPI_Allreduce(&resdlocal, dt0,1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
    
    std::cout << "\n %%%%%%%%%%%%%%%%%%%  resedual=  " << *dt0 ;
    
    MPI_Isend(mresult,(int)sizeof(mresult), MPI_DOUBLE, 0, rank, MPI_COMM_WORLD,&request);
  
    if(rank==0)
    {
     int jn=0;
     for( int i=0; i< size; i++)
    {
      MPI_Recv(nresult,(int)sizeof(mresult), MPI_DOUBLE, i, i, MPI_COMM_WORLD,&status);
    
      for(int l=0; l< (int)sizeof(nresult);l++)
          Rvec[++jn]= nresult[l];
          
          MPI_Isend(&Rvec[0],(int)Rvec.size(),MPI_DOUBLE,i,i*10,MPI_COMM_WORLD,&request);
    
    }
  } 
  
    
     MPI_Recv(&Rvec[0],(int)Rvec.size(),MPI_DOUBLE,0,rank*10,MPI_COMM_WORLD,&status);
        
  
  
    if(abs(*dt0) > abs(error))
     {    
        std::vector<double> Dvec (Rvec); 
      
        for(int i = 0 ; i<iter; i++)
       {

        if(i>0 && rank!=0)
        {
           MPI_Recv(&broke,1,MPI_INT,0,rank*13,MPI_COMM_WORLD,&status);
           if(broke == 1)
           break;
           std::cout << rank << " &&&&@@@@@@@@@@@@########## " << "\n"; 
           MPI_Recv(&Dvec[0],(int)Dvec.size(),MPI_DOUBLE,0,rank*11,MPI_COMM_WORLD,&status);
        }

        tresult = matMult(Dvec, blenx,bleny,sx, alfa, bita, gama, destn,dests);
        MPI_Isend(tresult,(int)sizeof(tresult), MPI_DOUBLE, 0, rank*19, MPI_COMM_WORLD,&request);
        
        if(rank == 0)
        {       
            //std::cout << rank << " &&&&@@@@@@@@@@@@########## " << i << "\n";      
            int jk = 0;
            for( int km=0; km < size; km++)
            {
            //   std::cout << rank << " " << "33" <<std::endl; 
              MPI_Recv(nresult,(int)sizeof(tresult), MPI_DOUBLE,km, km*19, MPI_COMM_WORLD,&status);
             // std::cout << "5### " << "\n";
              for(int l=0; l< (int)sizeof(nresult);l++)
                  Tvec[++jk]= nresult[l];
            }  
            //std::cout << "5### " << "\n";
            double dt = std::inner_product(Dvec.begin(), Dvec.end(), Tvec.begin(),0);
           
            //double dt =1;

            alpha = *dt0 / dt;
             //std::cout << rank << " 6### " << "\n";   
        for(int j=0; j< (int)Dvec.size();j+=4)
        {
            Tmpvec[j] = alpha * Dvec[j];
             Tmpvec[j+1] = alpha * Dvec[j+1];
              Tmpvec[j+2] = alpha * Dvec[j+2];
               Tmpvec[j+3] = alpha * Dvec[j+3]; 
        }     
        //std::cout << rank <<  "7### " << "\n";
        std::transform (Xvec.begin(), Xvec.end(), Tmpvec.begin(), Xvec.begin(),   std::plus<double>());
         // std::cout << rank << "8### " << "\n";
        for(int j=0; j< (int)Tvec.size();j+=4)
        {
           Tmpvec[j] = alpha * Tvec[j];
             Tmpvec[j+1] = alpha * Tvec[j+1];
              Tmpvec[j+2] = alpha * Tvec[j+2];
               Tmpvec[j+3] = alpha * Tvec[j+3]; 
            
        }      
        
        std::transform (Rvec.begin(), Rvec.end(), Tmpvec.begin(), Rvec.begin(),  std::minus<double>());

        double dt1 = std::inner_product(Rvec.begin(), Rvec.end(), Rvec.begin(),0);

        if(abs(dt1) < abs(error))
        {
            broke = 1;

         for(int jb=1; jb< size;jb++)
        {
          MPI_Isend(&broke,1,MPI_INT,jb,jb*13,MPI_COMM_WORLD,&request);   
        } 
            break;
        }

        double beta = dt1/(*dt0);
        
         for(int j=0; j< (int)Dvec.size();j+=4)
        {
             Tmpvec[j] = beta * Dvec[j];
             Tmpvec[j+1] = beta * Dvec[j+1];
              Tmpvec[j+2] = beta * Dvec[j+2];
               Tmpvec[j+3] = beta * Dvec[j+3]; 
        }      
        
        std::transform (Rvec.begin(), Rvec.end(), Tmpvec.begin(), Dvec.begin(),   std::plus<double>());
        *dt0 = dt1;
        for(int j=1; j< size;j++)
        {
          MPI_Isend(&broke,1,MPI_INT,j,j*13,MPI_COMM_WORLD,&request);   
          MPI_Isend(&Dvec[0],(int)Dvec.size(),MPI_DOUBLE,j,j*11,MPI_COMM_WORLD,&request);  
        }   
        
        }
    }

}       
    
    if(rank == 0)
    {
        for (int i= 0; i< totdim; i++ )
        std::cout << Xvec[i] << ' ';
    	time = timer.elapsed();
	std::cout << rank << " time," << time << std::endl;
	
    }
    MPI_Finalize();
#ifdef USE_LIKWID
	likwid_markerStopRegion("dummy");
	likwid_markerClose();
#endif

    

    return 0;

}
