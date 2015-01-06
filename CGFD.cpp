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
#include <stdlib.h>
#define LD 64
#define domxl 0.0
#define domxh 2.0
#define domyl 0.0
#define domyh 1.0
#define BLOCKSIZE 4

#define LD1 3

# define k 2.0*M_PI

#define ERRLIMIT 0;
//int iterat = 0;

inline double fxy(const double x, const double y){

         return 4.0*M_PI*M_PI*sin(2.0*M_PI*x)*sinh(2.0*M_PI*y);
    }

inline double border(const double x, const double y){
        return sin(2.0*M_PI*x)*sinh(2.0*M_PI*y);
    }
    
    
  inline double * gcapt(double * worksheet, int const procn, int dimM/*, int dimN*/)
   {
        
                      
       int m =0;
        int n =0;
      //  int p =0;
       // int s = 0;
        int ver =0;
    
         ver = procn;
       int remM = dimM%ver;
       int bn = dimM/ver;
        
        
       for (int i = 0; i < ver; i++)
         {
                    int q = i*LD1;
                //worksheet[q+1]=dimN;
                worksheet[q]=n;
                //worksheet[q+2]=s;
                //worksheet[q+3]=0;
                if(i<remM)
                {
                  m=bn+1;  
                  worksheet[q+1]=m; 
                  n += m;
                  //p= m*dimN;
                } 
                else
                {
                worksheet[q+1]=bn;  
                n += bn;
                //p = bn*dimN;
                }                
                //s += p;                 
               //worksheet[q+5]=dimN;
            }   
       return worksheet;
   }

inline double * matMult( double* vec,int blenx,int bleny,int sx, double alpha,  double beta,  double gama,
   /*int destn, int dests, */int len/*, int startpnt*/, double * start, double *end)
{  
    
     double * result = (double *) calloc(len, sizeof(double));
     int index=0;
     int le=sx+blenx;
       
    __m128d a,b,c,d,e,f,g;
   
    for(int i=sx; i<le ; i++)
      {
        for(int j=0; j<bleny ; j++)
        {            
             a[0]=alpha;
             a[1]= vec[index];
             b[0]= beta;
             b[1] =beta;
             c[0]=gama;
             c[1]=gama;
                          
            if(j > 0)
               d[0]=vec[index-1];            
            if(j==0)
               d[0] = 0.0;
               
            if(j < bleny-1)
                d[1] = vec[index+1];
            if(j == bleny-1)
               d[1] =0.0; 
            
              if(i==sx)
              e[0] = start[j];
              else
              e[0] = vec[index-bleny];
              
              
              if(i==le-1)
              e[1] = end[j];
              else
              e[1] = vec[index+bleny];
              
            f = _mm_mul_pd(b,d);
            g = _mm_mul_pd(c,e);  
              
      //      time = timer.elapsed();
            //if(iterat >75 && iterat < 78 && i==sx && j==0)
            //std::cout <<  " time, 2:" << time << "\n";          
            
            e = _mm_add_pd(f,g);
                
            result[index++]=e[0]+e[1]+a[0]*a[1];

        }
    }
    //std:: cout << "####";
 //     for(int i=0;i<len;i+=2)
//    {
//        std:: cout << " " << result[i] ;
//        std:: cout << " " << result[i+1] ;
//    }

    return result;
    
}

inline double * cal_fVec(int blenx,int bleny ,int sx,const double gama,  double hx, double hy,int dests, int len)
{ 
  
     double * result =(double *) calloc(len, sizeof(double));
    
    int gridno = 0;
     
    //double gama2 = 0.0;

    double x = 0.0;
    double y = 0.0;

    int le = sx+blenx;
    int l = 0;
    for(int i=sx; i< le ; i++)
    {
        for(int j=0; j<bleny ; j++)
        {
            gridno = i*bleny + j;            
            x = (((gridno)%bleny)+1)*hx;
            y = (((gridno)/bleny)+1)*hy;
            double f = fxy(x,y);
            if(dests == -1 && i == le-1)
            {
              result[l++] = f-(gama*border(x,domyh));
            }
            else 
                result[l++] = f;
        }
    }

       return result;
   
}


int main(int argc, char** argv)
{

    
    //std::cout << argc << "= invalid number of argument.. Program exiting..&&&&&&&&&&&&";
            
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
	double dt0=0.0;   
    double alfa=0.0;
    double bita=0.0;
    double gama=0.0;
    double dt1 = 0.0;
    double hx = 0.0, hy=0.0; 
    //Sint startpnt =0;
    double dt = 0.0; 

    
  // double ev=0.0,wv=0.0,sv=0.0,nv=0.0;
   // MPI_Datatype columntype;   
    //MPI_Type_vector( 2, 1, 2, MPI_DOUBLE, &columntype );
    //MPI_Type_commit( &columntype );  
    
   // std::cout << rank << "111 ";
      
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
    bita = -(1.0)/(hx*hx);
    gama = -(1.0)/(hy*hy);
    alfa = ( k * k - (2.0)*gama - (2.0)*bita );
     //std::cout << rank << "222 ";
    if (rank == 0)
   {
    worksheet = gcapt(worksheet, size, nny); 
    
    for(int t=0;t<size;t++)
    {
    sx = worksheet[t*3];
    blenx = worksheet[t*3+1];
    //startpnt = worksheet[t*3+2];
    
    
    MPI_Isend(&blenx,1,MPI_INT,t,t+100,MPI_COMM_WORLD,&request);
    MPI_Isend(&sx,1,MPI_INT,t,t+110,MPI_COMM_WORLD,&request);
   // MPI_Isend(&startpnt,1,MPI_INT,t,t+120,MPI_COMM_WORLD,&request);
    }
           
    }
     
     //std::cout << rank << "333 ";  
    MPI_Recv(&blenx,1, MPI_INT,0, rank+100, MPI_COMM_WORLD,&status);
    MPI_Recv(&sx,1, MPI_INT,0, rank+110, MPI_COMM_WORLD,&status);
   // MPI_Recv(&startpnt,1, MPI_INT,0, rank+120, MPI_COMM_WORLD,&status);
  
    std::cout << rank << " sx = " << sx << " " << status.MPI_ERROR; 
    std::cout << rank << " blenx = " << blenx << " " << status.MPI_ERROR; 
    
     // std::cout << rank << " startpoint = " << startpnt << " " << status.MPI_ERROR; 
  
 
   int abc = 0;
    int len=0; 
   //int bleny =  nnx;
              
    int sz=nnx*blenx;
    abc = sz%2;
    
    if(abc != 0)
    len = sz + 1;
    else 
    len = sz + 2;
        
    double * tresult = (double *) calloc(len, sizeof(double));
    double * fresult = (double *) calloc(len, sizeof(double));
    
   // double * nresult = new double[sz];
    double * Xvec = (double *) calloc(len, sizeof(double));
    double * Rvec = (double *) calloc(len, sizeof(double));
    double * Dvec = (double *) calloc(len, sizeof(double)) ;
    double * Fvec = (double *) calloc(gridpoint, sizeof(double));
  //  for(int i=0;i<len;i+=2)
//    {
//        std:: cout << " " << Dvec[i] ;
//        std:: cout << " " << Dvec[i+1] ;
//    }
       
    
 //   for(int i=0;i<len;i+=2)
//    {
//        Xvec[i]=0.0;
//        Xvec[i+1]=0.0;
//        Rvec[i]=0.0;
//        Rvec[i+1]=0.0;
//        Dvec[i]=0.0;
//        Dvec[i+1]=0.0;
//    }
    
    if(rank == size-1)
    dests = -1;
   
    double *start = (double *) calloc(nnx, sizeof(double));
    double *end = (double *) calloc(nnx, sizeof(double));
     

       tresult = matMult(Xvec,blenx,nnx,sx, alfa, bita,gama,/*destn,dests,*/len,start,end);        
       fresult = cal_fVec(blenx,nnx,sx,gama, hx ,hy,dests,len);
       
 //      std:: cout << "@@@@@@";
//      for(int i=0;i<len;i+=2)
//    {
//        std:: cout << " " << tresult[i] ;
//        std:: cout << " " << tresult[i+1] ;
//    }
    
  //   std:: cout << "$$$$$$$";
//      for(int i=0;i<len;i+=2)
//    {
//        std:: cout << " " << fresult[i] ;
//        std:: cout << " " << fresult[i+1] ;
//    }

       
       
    __m128d a,b,c,d,e,f,g,hh,ii,jj;   
       
   // std::cout << "\n" << rank << " " << iter << " " << blenx << " " << nnx << " " << sx << " " << gama << " " << hx << " " << hy << " " << startpnt;
      for(int i = 0; i< len; i+=2)
    {
        //std::cout << "\n" << rank << " " << fresult[i];
        a = _mm_load_pd(&fresult[i]);
        b = _mm_load_pd(&tresult[i]);
        c = _mm_sub_pd(a,b);
        _mm_store_pd (&Rvec[i], c);
        _mm_store_pd (&Dvec[i], c);
        d = _mm_mul_pd(c,c);
        resdlocal += d[0] + d[1];
     } 
    
    free(tresult);
    free(fresult);
     
   
    //std::cout << "\n %%%%%%%%%%%%%%%%%%%  resedual=  " <<  rank << " " << resdlocal;
    
    MPI_Allreduce(&resdlocal, &dt0,1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
    
    //std::cout << "\n %%%%%%%%%%%%%%%%%%%  resedual=  " << *dt0 ;
    //MPI_Barrier(MPI_COMM_WORLD);
   
    if(fabs(sqrt(dt0)) > error)
     {   
      //std::cout << "\n %% rank = " << rank << "before for "<< size ; 
      int ik=0;
      while(ik < iter)
       {
            //iterat = ik;
        //std::cout << "\n %% rank = " << rank << "iteration number= " << ik << "\n";
                
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
//       std::cout << "\n ";
       // for(int s=0;s<len;s++)
//         {
//              std::cout << rank << " " << Dvec[s] ;
//         }
//        
//                 std::cout << "\n ";
         MPI_Isend(&Dvec[0],nnx,MPI_DOUBLE, gn, gn+130, MPI_COMM_WORLD,&request);
         
         MPI_Recv(end,nnx, MPI_DOUBLE,hn, rank+130, MPI_COMM_WORLD,&status);
         
             
        
         MPI_Isend(&Dvec[sz-nnx],nnx,MPI_DOUBLE, hn, hn+140, MPI_COMM_WORLD,&request);
         MPI_Recv(start,nnx, MPI_DOUBLE,gn, rank+140, MPI_COMM_WORLD,&status);
        
        
        
        for( int r=0;r<nnx;r++)
         {
         if(rank==0)
        {       
         start[r]=0.0;
        }
         if(rank == size-1)
         {               
         end[r]=0.0;
         }
        }
        
     //   std::cout << "\n" << rank << " ghost layer end sent \n";
//         for(int s=0;s<nnx;s++)
//         {
//              std::cout << rank << " " << Dvec[s] ;
//         }
//         
//        std::cout << "\n" << rank << " ghost layer end  recvd \n";
//        
//        for(int s=0;s<nnx;s++)
//         {
//              std::cout << rank << " " << end[s] ;
//         }
//         
//         std::cout << "\n" << rank << "  ghost layer start sent\n";
//          for(int s=0;s<nnx;s++)
//         {
//              std::cout << rank << " " << Dvec[sz-nnx+s] ;
//         }
//        
//        std::cout << "\n" << rank << "  ghost layer start recvd\n";
//        for(int s=0;s<nnx;s++)
//         {
//              std::cout << rank << " " << start[s] ;
//         }
//        
         //std::cout <<  "\n" << rank << " sv nv ev wv " << sv << nv << ev << wv;
        
          //MPI_Barrier(MPI_COMM_WORLD);
          //time = timer.elapsed();
           //std::cout <<  " time, 1:" << time;
         double * mresult = matMult(Dvec, blenx,nnx,sx, alfa, bita, gama,len,start,end);
        // mresult = matMult(Dvec, blenx,nnx,sx, alfa, bita, gama,len,start,end);
       
       
      dt = 0.0;
          for( int km=0; km < len; km+=2)
            {
                
                a = _mm_load_pd(&Dvec[km]);
                b = _mm_load_pd(&mresult[km]);
                c = _mm_mul_pd(a,b);
                dt += c[0]+c[1];
                    
            }
        
         double dt3 = 0.0;  
         MPI_Allreduce(&dt, &dt3,1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
      
         alpha = dt0 / dt3;
     
       
         dt = 0.0;
         
        for(int j=0; j< len;j+=2)
        {
            
            a = _mm_load_pd(&Dvec[j]);
            b = _mm_load_pd(&mresult[j]);
            c[0] = alpha;
            c[1] = alpha;            
            d = _mm_load_pd(&Xvec[j]);
            e = _mm_load_pd(&Rvec[j]);
            
            f = _mm_mul_pd(c,a);
            g = _mm_mul_pd(c,b);
            hh = _mm_add_pd(d,f);
            ii = _mm_sub_pd(e,g);
            jj = _mm_mul_pd(ii,ii);
            
            _mm_store_pd (&Xvec[j], hh);
            _mm_store_pd (&Rvec[j], ii);
            
               dt+= jj[0]+jj[1];

        }
                    
        MPI_Allreduce(&dt, &dt1,1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);

        if(fabs(sqrt(dt1)) <= fabs(error))
        {
            std::cout << rank << "\n residue I am here to break you loop " << fabs(dt1) << " " << fabs(error) << "\n";
            break;
        }

        double beta = dt1/(dt0);        
        
         // std::cout << " 7: " <<  beta  <<  "\n";
         for(int j=0; j< len;j+=2)
        {
             a = _mm_load_pd(&Dvec[j]);
            b = _mm_load_pd(&Rvec[j]);
            c[0] = beta;
            c[1] = beta;
           // std::cout << " 6: " <<  c[0] << " " << beta << " " << c[1] <<"\n";            
            f = _mm_mul_pd(c,a);
            g = _mm_add_pd(b,f);    
            _mm_store_pd (&Dvec[j], g);
           
        }      
        // std::cout << " 7: " <<  "\n";
         dt0 = dt1;
         ik++;
         free(mresult);
   
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
    std::cout << "final residuum:  " << residual << " " << dt1 << "\n Writing data to data/solution.txt" ;
	
	std::ofstream	fOut("data/solution.txt");

    double xx =0.0, yy=0.0;
	for (int i= 0; i< gridpoint; i++ )
        {
       
       // gridno = i*bleny + j;            
            xx = ((i%nnx)+1)*hx;
            yy = ((i/nnx)+1)*hy;
         
        fOut << xx << " " << yy << " " << Fvec[i] << "\n";
        
        }
        fOut << std::endl;
        
    }
    
    
    MPI_Barrier(MPI_COMM_WORLD);
    free(Dvec);
    free(Rvec);
    free(Xvec);
    free(Fvec);
   //:wq
  // MPI_Type_free( &columntype );
    MPI_Finalize();
#ifdef USE_LIKWID
	likwid_markerStopRegion("dummy");
	likwid_markerClose();
#endif

    

    return 0;

}
