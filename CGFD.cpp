#include <iostream>
#include <numeric>
#include <fstream>
#include <algorithm>
#include "emmintrin.h"
#include <memory>
#include <vector>
#include <mpi.h>
#include	"Timer.h"
#include	<cmath>
#include <functional>
#include "immintrin.h"
#include <memory>
#define LD 64
#define domxl 0.0
#define domxh 2.0
#define domyl 0.0
#define domyh 1.0
#define BLOCKSIZE 4

#define LD1 3

# define k 2.0*M_PI

#define ERRLIMIT 0;
int iterat = 0;

inline double fxy(const double x, const double y){

         return 4.0*M_PI*M_PI*sin(2.0*M_PI*x)*sinh(2.0*M_PI*y);
    }

inline double border(const double x, const double y){
        return sin(2.0*M_PI*x)*sinh(2.0*M_PI*y);
    }
    
    
  inline double * gcapt(double * worksheet, int const procn, int dimM, int dimN)
   {
        
                      
       int m =0;
        int n =0;
        int p =0;
        int s = 0;
        int ver =0;
    
       //worksheet = new double[proc*LD1];

       if(procn%4==0)
           ver = 4;
       else if (procn%5==0)
           ver = 5;
       else if (procn%6==0)
           ver = 6;
       else if (procn%7==0)
           ver = 7;
       int remM = dimM%ver;
       int bn = dimM/ver;
        
        
       for (int i = 0; i < ver; i++)
         {
                    int q = i*LD1;
                //worksheet[q+1]=dimN;
                worksheet[q]=n;
                worksheet[q+2]=s;
                //worksheet[q+3]=0;
                if(i<remM)
                {
                  m=bn+1;  
                  worksheet[q+1]=m; 
                  n += m;
                  p= m*dimN;
                } 
                else
                {
                worksheet[q+1]=bn;  
                n += bn;
                p = bn*dimN;
                }                
                s += p;                 
               //worksheet[q+5]=dimN;
            }   
       return worksheet;
   }

inline double * matMult( double* vec,int blenx,int bleny,int sx,const double alpha, const double beta, const double gama,
   /*int destn, int dests, */int len, int startpnt, double ev, double wv, double nv, double sv)
{  
    
    
    int le=sx+blenx;
  	//double time = 0;
    siwir::Timer	timer;
     double * result = new double[len];
    int gridno = 0;
    //double gama1 = 0.0;
    //double gama2 = 0.0;
    int index=0;
   
    int l =0;
    gridno= sx*bleny;
    index = gridno-startpnt;
    __m128d a,b,c,d,e,f,g;
   
    for(int i=sx; i<le ; i++)
      {
        for(int j=0; j<bleny ; j++,gridno++, index++)
        {
            //int pm = 0;
            //int ppm =0;
            //int fm =0;
            //int ffm =0;  
            //time = timer.elapsed();
            //if(iterat >75 && iterat < 78 && i==sx && j==0)
            //std::cout <<  " time, 1:" << time << " " << vec[index]; 
             a[0]=alpha;
             a[1]= vec[index];;
             b[0]= beta;
             b[1] =beta;
             c[0]=gama;
             c[1]=gama;
             //b[0] = vec[index];
              
            if(j!=0 && i!=sx)
               d[0]=vec[index-1];
               
            if(j!=0 && i==sx)
               d[0]=wv;
            if(j==0)
               d[0] = 0.0;
               
            if(j != bleny-1 && i!=le)
                d[1] = vec[index+1];
            if(j != bleny-1 && i==le)
                d[1] = ev;
            if(j == bleny-1)
               d[1] =0.0; 
            
            if(index-3>=0)
              e[0] = vec[index-3];
              else
              e[0] = nv;
              
              if(index+3>=0)
              e[1] = vec[index+3];
              else
              e[1] = sv;
         
            f = _mm_mul_pd(b,d);
            g = _mm_mul_pd(c,e);  
              
      //      time = timer.elapsed();
            //if(iterat >75 && iterat < 78 && i==sx && j==0)
            //std::cout <<  " time, 2:" << time << "\n";          
            
            e = _mm_hadd_pd(f,g);
                
            result[l++]=e[0]+e[1]+a[0]*a[1];

        }
    }

    return result;
    
}

inline double * cal_fVec(int blenx,int bleny ,int sx,const double gama,  double hx, double hy, int dests, int len)
{ 
  // vector<double> fresult(tgrdpoint,0);

    //int bleny =  dim[1];    
    //std::cout << rank << " in cal fvec ";
       
    //int sy = 0;
   // int len=0;
    //int sz=blenx*bleny;
   // int abc = sz%4;
    /*if(abc != 0)
    len = sz + (4-sz%4);
    else
    len = sz;*/
    double * result = new double[len];
    
    int gridno = 0;
     
    double gama2 = 0.0;

    double x = 0.0;
    double y = 0.0;


    int l = 0;
    for(int i=sx; i< sx+blenx ; i++)
    {
        for(int j=0; j<bleny ; j++)
        {
            
            gama2 = 0.0;
            gridno= i*bleny + j;            
            //int k = (j-sy)%blenx;
            x = (((gridno-1)%bleny)+1)*hx;
            y = (((gridno-1)/bleny)+1)*hy;
            //std::cout << "x, y" << x << " " <<y;
            double f = fxy(x,y);
           // std::cout << " x, y, f " << x << " " <<y << " " << f << "\n";                      
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
    

    int gridpoint = 0;
    //int len = 0;
    int dests=0, blenx=0, sx =0;
    int nx = 0;
    int ny = 0;
    int iter = 0;
    //int broke=0;  
    int nnx =0, nny=0;
    double alpha = 0.0;
    double resdlocal=0.0; 
    double error=0.0;
    double* worksheet = new double[size*LD1];;
   	double time = 0;
	double dt0[1]={0.0};   
    double alfa=0.0;
    double bita=0.0;
    double gama=0.0;
     double dt1 = 0.0;
    double hx = 0.0, hy=0.0; 
    int startpnt =0;
    double dt = 0.0; 
    double *start = new double[2];
    double *end = new double[2];
    
   double ev=0.0,wv=0.0,sv=0.0,nv=0.0;
    MPI_Datatype columntype;   
    MPI_Type_vector( 2, 1, 2, MPI_DOUBLE, &columntype );
    MPI_Type_commit( &columntype );  
      
    //*dt0=0.0;   
    nx = atoi(argv[1]);
    ny = atoi(argv[2]);
    nnx = nx-1;
    nny = ny-1;
    iter = atoi(argv[3]);
    error = atof(argv[4]);
    gridpoint = nnx*nny;
    hx = (double)((domxh-domxl)/(double)(nx));
    hy = (double)((domyh-domyl)/(double)(ny));
    bita = 1/(hx*hx);
    gama = 1/(hy*hy);
    alfa = (-1.0)*(2.0*gama+ 2.0*bita + k * k);


    if (rank == 0)
   {
    worksheet = gcapt(worksheet, size, nny, nnx); 
    
    for(int t=0;t<size;t++)
    {
    sx = worksheet[t*3];
    blenx = worksheet[t*3+1];
    startpnt = worksheet[t*3+2];
    
    
    MPI_Isend(&blenx,1,MPI_INT,t,t+100,MPI_COMM_WORLD,&request);
    MPI_Isend(&sx,1,MPI_INT,t,t+110,MPI_COMM_WORLD,&request);
    MPI_Isend(&startpnt,1,MPI_INT,t,t+120,MPI_COMM_WORLD,&request);
    }
           
    }
      
    MPI_Recv(&blenx,1, MPI_INT,0, rank+100, MPI_COMM_WORLD,&status);
    MPI_Recv(&sx,1, MPI_INT,0, rank+110, MPI_COMM_WORLD,&status);
    MPI_Recv(&startpnt,1, MPI_INT,0, rank+120, MPI_COMM_WORLD,&status);
  
    std::cout << rank << " sx = " << sx << " " << status.MPI_ERROR; 
    std::cout << rank << " blenx = " << nnx << " " << status.MPI_ERROR; 
     std::cout << rank << " startpoint = " << startpnt << " " << status.MPI_ERROR; 
  
  /*std::cout << rank << " nnx = " << nnx << " "; 
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
   
   int abc = 0;
     
   //int bleny =  nnx;
              
    int sz=nnx*blenx;
    abc = sz%4;
    
    if(abc != 0)
    sz = sz + (4-abc);
    else 
    sz = sz + 4;
        
    double * tresult = new double[sz];
    double * fresult = new double[sz];
    //double * mresult = new double[sz];
   // double * nresult = new double[sz];
    double * Xvec = new double[sz];
    double * Rvec = new double[sz];
    double * Dvec = new double[sz]; 
    double * Fvec = new double[gridpoint];
    
    for(int i=0;i<sz;i+=4)
    {
        Xvec[i]=0.0;
        Xvec[i+1]=0.0;
        Xvec[i+2]=0.0;
        Xvec[i+3]=0.0;
    }
    //double * Tvec = new double[sz];
    //double * Tmpvec = new double[sz];
              
    //double resd =0.0;
   
   /* if(rank == 0)
    destn = -1;
    else 
    destn = rank -1;*/  
    
    if(rank == size-1)
    dests = -1;
     
    
       tresult = matMult(Xvec,blenx,nnx,sx, alfa, bita,gama,/*destn,dests,*/sz,startpnt,0.0,0.0,0.0,0.0);        
       fresult = cal_fVec(blenx,nnx,sx,gama, hx ,hy,dests,sz);
       
       
   // std::cout << "\n" << rank << " " << iter << " " << blenx << " " << nnx << " " << sx << " " << gama << " " << hx << " " << hy << " " << startpnt;
      for(int i = 0; i< sz; i+=4)
    {
       // std::cout << "\n" << rank << " " << fresult[i];
        Rvec[i] = fresult[i]-tresult[i];
        Rvec[i+1] = fresult[i+1]-tresult[i+1];
        Rvec[i+2] = fresult[i+2]-tresult[i+2];
        Rvec[i+3] = fresult[i+3]-tresult[i+3];
        Dvec[i]=Rvec[i];
        Dvec[i+1]=Rvec[i+1];
        Dvec[i+2]=Rvec[i+2];
        Dvec[i+3]=Rvec[i+3];
        resdlocal += Rvec[i] * Rvec[i];
        resdlocal += Rvec[i+1] * Rvec[i+1];
        resdlocal += Rvec[i+2] * Rvec[i+2];
        resdlocal += Rvec[i+3] * Rvec[i+3];
       //std::cout << "\n" << rank << " " << (int)sizeof(tresult) << " " << fresult[i] << " " << tresult[i] << " " << mresult[i];
    } 
     
   
   // std::cout << "\n %%%%%%%%%%%%%%%%%%%  resedual=  " <<  rank << " " << resdlocal;
    
    MPI_Allreduce(&resdlocal, dt0,1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
    
    //std::cout << "\n %%%%%%%%%%%%%%%%%%%  resedual=  " << *dt0 ;
    //MPI_Barrier(MPI_COMM_WORLD);
   
    if(fabs(sqrt(*dt0)) > fabs(error))
     {   
      //std::cout << "\n %% rank = " << rank << "before for "<< size ; 
      int ik=0;
      while(ik < iter)
       {
            iterat = ik;
        //std::cout << "\n %% rank = " << rank << "iteration number= " << ik << "\n";
        ev=0.0;
        wv=0.0;
        sv=0.0;
        nv=0.0;
        
        int gn=0;
        int hn =0;
        if(rank!=0)
        gn=rank-1;
        else 
        gn = size-1;
        if(rank!=size-1)
        hn=rank+1;
        else 
        hn = 0;
        
         MPI_Isend(&Dvec[0],1,columntype, gn, gn+130, MPI_COMM_WORLD,&request);
         
         MPI_Recv(end,1, columntype,hn, rank+130, MPI_COMM_WORLD,&status);
         if(rank != size-1)
         {
         ev = end[0];
         sv = end[1];
        }
        
         MPI_Isend(&Dvec[sz-3],1,columntype, hn, hn+140, MPI_COMM_WORLD,&request);
         MPI_Recv(start,1, columntype,gn, rank+140, MPI_COMM_WORLD,&status);
        if(rank!=0)
        {
          wv = start[0];
          nv = start[1];
        }
          //MPI_Barrier(MPI_COMM_WORLD);
//          time = timer.elapsed();
//           std::cout <<  " time, 1:" << time;
         tresult = matMult(Dvec, blenx,nnx,sx, alfa, bita, gama,sz,startpnt,ev,wv,nv,sv);
        // time = timer.elapsed();
//         std::cout << " 2: " << time ;
         
          for( int km=0; km < sz; km+=4)
            {
                    dt += Dvec[km]*tresult[km];
                     dt += Dvec[km+1]*tresult[km+1];
                      dt += Dvec[km+2]*tresult[km+2];
                       dt += Dvec[km+3]*tresult[km+3];
            }
          time = timer.elapsed();
         //std::cout << " 3: " << time << "\n";
         double dt3 = 0.0;  
         MPI_Allreduce(&dt, &dt3,1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
          // time = timer.elapsed();
	   //std::cout <<  " time, 4:" << time;
         alpha = *dt0 / dt3;
         // time = timer.elapsed();
//	    std::cout << " 5: " << time << "\n";
        
         dt = 0.0;
         
        for(int j=0; j< sz;j+=4)
        {
            Xvec[j] += alpha * Dvec[j];
            Xvec[j+1] += alpha * Dvec[j+1];
             Xvec[j+2] += alpha * Dvec[j+2];
               Xvec[j+3] += alpha * Dvec[j+3]; 
               Rvec[j] -= alpha * tresult[j];
             Rvec[j+1] -= alpha * tresult[j+1];
              Rvec[j+2] -= alpha * tresult[j+2];
               Rvec[j+3] -= alpha * tresult[j+3]; 
               dt += Rvec[j]*Rvec[j];
              dt += Rvec[j+1]*Rvec[j+1];
                 dt += Rvec[j+2]*Rvec[j+2];
                dt += Rvec[j+3]*Rvec[j+3];
        }
                    
        MPI_Allreduce(&dt, &dt1,1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);

        if(fabs(sqrt(dt1)) < fabs(error))
        {
            std::cout << rank << "\n residue I am here to break you loop " << fabs(dt1) << " " << fabs(error) << "\n";
            break;
        }

        double beta = dt1/(*dt0);        
        
         for(int j=0; j< sz;j+=4)
        {
            Dvec[j]   = Rvec[j]-beta * Dvec[j];
            Dvec[j+1] = Rvec[j+1]-beta * Dvec[j+1];
            Dvec[j+2] = Rvec[j+2]-beta * Dvec[j+2];
            Dvec[j+3] = Rvec[j+3]-beta * Dvec[j+3]; 
        }      
        
         *dt0 = dt1;
         ik++;
   
        }                
    }       
     //std::cout << "\n %%%%%%%%%%%%%%%%%%%  at end " ;
     MPI_Isend(&sz,1,MPI_INT, 0, rank+49, MPI_COMM_WORLD,&request); 
     MPI_Isend(&Xvec[0],sz,MPI_DOUBLE, 0, rank+39, MPI_COMM_WORLD,&request); 
   
    
   // std::cout << "\n %%%%%%%%%%%%%%%%%%%  at end " ;
    
    if(rank == 0)
    {     
       int nsz =0; 
       int y=0;
      for(int j=0; j<size;j++)
      {
        y+=nsz;
        MPI_Recv(&nsz,1, MPI_INT,j, j+49, MPI_COMM_WORLD,&status);     
        MPI_Recv(&Fvec[y],nsz, MPI_DOUBLE,j, j+39, MPI_COMM_WORLD,&status); 
      }
    time = timer.elapsed();
	std::cout << rank << " time," << time << std::endl;
	double residual =0.0;
	residual = sqrt(dt1);
    std::cout << "finale residuam:  " << residual << "\n Writing data to data/solution.txt" ;
	
	std::ofstream	fOut("data/solution.txt");

	for (int i= 0; i< gridpoint; i++ )
        {
        //if(i%nnx == 0 )
        //std::cout << "\n"; 
        //residual +=  Fvec[i]*Fvec[i];
         
        fOut << Fvec[i] << "\n";
        
        }
        fOut << std::endl;
        
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Type_free( &columntype );
    MPI_Finalize();
#ifdef USE_LIKWID
	likwid_markerStopRegion("dummy");
	likwid_markerClose();
#endif

    

    return 0;

}
