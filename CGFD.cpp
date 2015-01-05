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
//int iterat = 0;

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
   /*int destn, int dests, */int len/*, int startpnt*/, double * start, double *end)
{  
    
    
    int le=sx+blenx;
  	//double time = 0;
   // siwir::Timer	timer;
     double * result = new double[len];
     int index=0;
   
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
                          
            if(index > 0)
               d[0]=vec[index-1];            
            if(j==0)
               d[0] = 0.0;
            if(index == 0)
               d[0]=start[2];  
               
            if(index < len-1)
                d[1] = vec[index+1];
            if(j == bleny-1)
               d[1] =0.0; 
            if(index == len-1)
               d[1] = end[0];
            
             if(index-1 < 0)
              e[0] = start[0];
              else if(index-2 < 0)
              e[0] = start[1];
              else if(index-3 < 0)
              e[0] = start[2];
              else
              e[0] = vec[index-3];
              
              
              if(index +1 ==len)
              e[1] = end[2];
              else if(index +2 ==len)
              e[1] = end[1];
               else if(index+3 ==len)
              e[1] = end[0];
              else
              e[1] = vec[index+3];
            
            f = _mm_mul_pd(b,d);
            g = _mm_mul_pd(c,e);  
              
      //      time = timer.elapsed();
            //if(iterat >75 && iterat < 78 && i==sx && j==0)
            //std::cout <<  " time, 2:" << time << "\n";          
            
            e = _mm_hadd_pd(f,g);
                
            result[index++]=e[0]+e[1]+a[0]*a[1];

        }
    }

    return result;
    
}

inline double * cal_fVec(int blenx,int bleny ,int sx,const double gama,  double hx, /*double hy,*/int dests, int len)
{ 
  
     double * result = new double[len];
    
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
            
            //gama2 = 0.0;
            gridno = i*bleny + j;            
            //int k = (j-sy)%blenx;
            x = (((gridno)%bleny)+1)*hx;
            //y = (((gridno)/bleny)+1)*hy;
            //std::cout << "x, y" << x << " " <<y;
            double f = fxy(x,y);
           // std::cout << " x, y, f " << x << " " <<y << " " << f << "\n";                      
            if(dests == -1 && i == le-1)
            {
             //gama2 = 
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
    double *start = new double[3];
    double *end = new double[3];
    
  // double ev=0.0,wv=0.0,sv=0.0,nv=0.0;
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
    bita = -(1.0)/(hx*hx);
    gama = -(1.0)/(hy*hy);
    alfa = ((-2.0)*gama+ (-2.0)*bita + k * k);


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
    int len=0; 
   //int bleny =  nnx;
              
    int sz=nnx*blenx;
    abc = sz%2;
    
    if(abc != 0)
    len = sz + 1;
    else 
    len = sz + 2;
        
    double * tresult = new double[len];
    double * fresult = new double[len];
    
   // double * nresult = new double[sz];
    double * Xvec = new double[len];
    double * Rvec = new double[len];
    double * Dvec = new double[len]; 
    double * Fvec = new double[gridpoint];
    
    
    
    for(int i=0;i<len;i+=2)
    {
        Xvec[i]=0.0;
        Xvec[i+1]=0.0;
        Rvec[i]=0.0;
        Rvec[i+1]=0.0;
        Dvec[i]=0.0;
        Dvec[i+1]=0.0;
        //tresult[i]=0.0;
//        tresult[i+1]=0.0;
//        fresult[i]=0.0;
//        fresult[i+1]=0.0;        
    }
    
    if(rank == size-1)
    dests = -1;
     
      end[0]=0.0;
         end[1]=0.0;
         end[2]=0.0;
         
          start[0]=0.0;
          start[1]=0.0;
          start[2]=0.0;
         
       tresult = matMult(Xvec,blenx,nnx,sx, alfa, bita,gama,/*destn,dests,*/sz,start,end);        
       fresult = cal_fVec(blenx,nnx,sx,gama, hx ,dests,sz);
       
       if(abc != 0)
       {
        tresult[sz]=0.0;
        fresult[sz]=0.0;
      }
        else
       {
        tresult[sz]=0.0;
        tresult[sz+1]=0.0;
        fresult[sz]=0.0;
        fresult[sz+1]=0.0;
       }
       
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
        /*a = _mm_load_pd(&fresult[i+2]);
        b = _mm_load_pd(&tresult[i+2]);
        c = _mm_sub_pd(a,b);
      _mm_store_sd (&Rvec[i+2], c);
        _mm_store_sd (&Dvec[i+2], c);
         d = _mm_mul_pd(c,c);
        resdlocal += d[0] + d[1];*/
       //std::cout << "\n" << rank << "I am here 1 " ;
    } 
    
    delete[] tresult;
    delete[] fresult;
     
   
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
        
         MPI_Isend(&Dvec[0],3,MPI_DOUBLE, gn, gn+130, MPI_COMM_WORLD,&request);
         
         MPI_Recv(end,3, MPI_DOUBLE,hn, rank+130, MPI_COMM_WORLD,&status);
         if(rank == size-1)
         {
         end[0]=0.0;
         end[1]=0.0;
         end[2]=0.0;
        }
        
         MPI_Isend(&Dvec[sz-3],3,MPI_DOUBLE, hn, hn+140, MPI_COMM_WORLD,&request);
         MPI_Recv(start,3, MPI_DOUBLE,gn, rank+140, MPI_COMM_WORLD,&status);
        if(rank==0)
        {
          start[0]=0.0;
          start[1]=0.0;
          start[2]=0.0;
        }
         //std::cout <<  "\n" << rank << " sv nv ev wv " << sv << nv << ev << wv;
        
          //MPI_Barrier(MPI_COMM_WORLD);
          //time = timer.elapsed();
           //std::cout <<  " time, 1:" << time;
         double * mresult = new double[len];
         mresult = matMult(Dvec, blenx,nnx,sx, alfa, bita, gama,sz,start,end);
        if(abc != 0)
       {
        mresult[sz]=0.0;
       }
        else
       {
        mresult[sz]=0.0;
        mresult[sz+1]=0.0;
       }
       
                 dt = 0.0;
          for( int km=0; km < len; km+=2)
            {
                
                a = _mm_load_pd(&Dvec[km]);
                b = _mm_load_pd(&mresult[km]);
                c = _mm_mul_pd(a,b);
                dt += c[0]+c[1];
                    
            }
          //time = timer.elapsed();
         //std::cout << " 3: " << time << "\n";
         double dt3 = 0.0;  
         MPI_Allreduce(&dt, &dt3,1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
          // time = timer.elapsed();
	   //std::cout <<  " time, 4:" << time;
	     std::cout << rank<<  " 5: rank " <<  dt3<<"\n";
         alpha = *dt0 / dt3;
         // time = timer.elapsed();
//	    std::cout << " 5: " << time << "\n";
       
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
            
            /*  a = _mm_load_pd(&Dvec[j+2]);
            b = _mm_load_pd(&tresult[j+2]);
            //c[0] = alpha;
            //c[1] = alpha;            
            d = _mm_load_pd(&Xvec[j+2]);
            e = _mm_load_pd(&Rvec[j+2]);
            
            f = _mm_mul_pd(c,a);
            g = _mm_mul_pd(c,b);
            hh = _mm_add_pd(d,f);
            ii = _mm_sub_pd(e,g);
            jj = _mm_mul_pd(hh,hh);
            
             _mm_store_sd (&Xvec[j+2], hh);
            _mm_store_sd (&Rvec[j+2], ii);
            dt+= jj[0]+jj[1];*/
        }
                    
        MPI_Allreduce(&dt, &dt1,1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);

        if(fabs(sqrt(dt1)) <= fabs(error))
        {
            std::cout << rank << "\n residue I am here to break you loop " << fabs(dt1) << " " << fabs(error) << "\n";
            break;
        }

        double beta = dt1/(*dt0);        
        
          std::cout << " 7: " <<  beta  <<  "\n";
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
            //std::cout << " 7: " <<  g[0] << " " << f[0] << " " << f[1] << " " << g[1] <<"\n";
           /*
            a = _mm_load_pd(&Dvec[j+2]);
            b = _mm_load_pd(&Rvec[j+2]);
            f = _mm_mul_pd(c,a);
            g = _mm_add_pd(b,f);
                       
            _mm_store_sd (&Dvec[j+2], g); */
        }      
        // std::cout << " 7: " <<  "\n";
         *dt0 = dt1;
         ik++;
         delete[] mresult;
   
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
    std::cout << "finale residuam:  " << residual << " " << dt1 << "\n Writing data to data/solution.txt" ;
	
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
