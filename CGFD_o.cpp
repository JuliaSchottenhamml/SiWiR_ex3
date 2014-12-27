#include <iostream>
#include <numeric>
#include <fstream>
#include "FDGRID.h"
#include "immintrin.h"
#include <memory>
#include <vector>
#include <mpi.h>
#include <GridCaptain.h>

#define ERRLIMIT 0;

inline double compute2norm(std::vector vec)
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


inline std::vector matMult(FDGRID& fgrid, vector<double> vec, GridCaptain& gcap)
{

   int size(0); // The total number of processes
   int rank(0); // The rank/number of this process (within MPI_COMM_WORLD)
   int proc = gcap.proc;
   int sz=vec.size();
   double * result = new double[sz];
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

   if (rank == 0)
   {
   int * dim = new int [2];

    dim[0]=fgrid.getDimM();
    dim[1]=fgrid.getDimN();
   
    MPI_Bcast(dim,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(result,1,MPI_INT,0,MPI_COMM_WORLD);
   }
    int sx = gcap.worksheet[2];
    int sy = gcap.worksheet[3];
    int ex = gcap.worksheet[4];
    int ey = gcap.worksheet[5];
    /*int en = gcap.worksheet[4];
    int wn = gcap.worksheet[5];
    int nn = gcap.worksheet[6];
    int sn = gcap.worksheet[7];*/
    int eval=0, wval = 0, sval = 0, nval = 0, cval = 0;
    int blenx = gcap.worksheet[0];
    int bleny = gcap.worksheet[1];

    int gridno = 0;

    double *ss = new double[blenx];
    double *ns = new double[blenx];

    double *es = new double[blenx];

    double *ws = new double[blenx];

    double *sr = new double[blenx];
    double *nr = new double[blenx];

    double *er = new double[blenx];

    double *wr = new double[blenx];
    
    for ( int i=0; i <blenx;i++)
    {
        sr[i]=0;
        er[i]=0;
        wr[i]=0;
        nr[i]=0;
    }

    ss=fgrid.getDataAdd(sx,sy);
    ns=fgrid.getDataAdd(ex,sy);
    es=fgrid.getDataAdd(sx,ey);
    ws=fgrid.getDataAdd(sx,sy);

    MPI_Datatype columntp;
    MPI_Type_vector(bleny,1,blenx, MPI_DOUBLE, &columntp);
    MPI_Type_commit( &columntp );

    int re = MPI_Cart_shift(MPI_COMM_WORLD,1,1,rank, rank+2);
    int rw = MPI_Cart_shift(MPI_COMM_WORLD,1,-1, rank, rank-2);
    int rn = MPI_Cart_shift(MPI_COMM_WORLD,0,1,rank,rank+2*gcap.hor );
    int rs = MPI_Cart_shift(MPI_COMM_WORLD,0,-1,rank,rank-2*gcap.hor);

     MPI_ISend(es,1,columntp,re, MPI_COMM_WORLD);
     MPI_ISend(ws,1,columntp,rw, MPI_COMM_WORLD);
     MPI_ISend(ss,1,MPI_double,rs, MPI_COMM_WORLD);
     MPI_ISend(ns,1,MPI_double,rn, MPI_COMM_WORLD);

     MPI_recv((er,1,columntp,re, MPI_COMM_WORLD);
     MPI_recv(wr,1,columntp,rw, MPI_COMM_WORLD);
     MPI_recv(sr,1,MPI_double,rs, MPI_COMM_WORLD);
     MPI_recv(nr,1,MPI_double,rn, MPI_COMM_WORLD);

    for(int i=sx; i<blenx ; i++)
    {
        int l = (i-sx)%bleny;
        for(int j=sy; j<bleny ; j++)
        {
            
            gridno= i*dim[1] + j;
            
            cval = fgrid.getDataA(i,j);
            int k = (j-sy)%blenx;
                                   
            if(i==sx)
            {
                  nval == nr[k];
            }
            else
            {
                nval=fgrid.getDataA(i-dimN,j);
            }    
            
            if(j==sy)
            {
             wval == wr[l];   
            }
            else
            {
                wval=fgrid.getDataA(i,j-1);
            }    
            
            if(i == blenx-1)
            {
                sval = sr[k];
            }
            else
            {
                sval=fgrid.getDataA(i+dimN,j);
            }
           
            if(j == bleny-1)
            {
                eval = er[l];
            }
            else
            {
                eval=fgrid.getDataA(i,j+1);
            }      
            
            int pm = 0;
            int ppm =0;
            int fm =0;
            int ffm =0;   
            int cm = cval*vec[gridno];
            if(gridno-1>=0)
            pm = wval*vec[gridno-1];
            if(gridno-3>=0)
            ppm = nval*vec[gridno-3];
            if(gridno+1<sz)
            fm = eval*vec[gridno+1];
            if(gridno+3<sz)
            ffm = sval*vec[gridno+3];
            result[gridno]=cm+pm+ppm+fm+ffm;

        }
    }




    MPI_Finalize();



}

inline double * callCG(FDGRID& fgrid, int const iter, int const proc, int const err)
{

    std::vector<double> Xvec (fgrid.totalGridPoints(),0);
    std::vector<double> Rvec (fgrid.totalGridPoints(),0);
    std::vector<double> Fvec (fgrid.totalGridPoints(),0);
    std::vector<double> Tvec;
     std::vector<double> Tmpvec;
    double alpha = 0;

    GridCaptain gcap = new GridCaptain(fgrid);
    
    TVec = matMult(fgrid,Xvec,gcap);
    
    cal_fVec(fgrid);

    std::transform (Fvec.begin(), Fvec.end(), TVec.begin(), Rvec.begin(),  std::minus<double>());

    double dt0 = std::inner_product(Rvec.begin(), Rvec.end(), Rvec.begin(),0);

    double resd = compute2norm(Rvec);

    if(resd < err)
        return Xvec;

      std::vector<double> Dvec (Rvec);

    for(int i = 0 ; i<iter; i++)

    {
        TVec = matMult(fgrid,Dvec,gcap);

        double dt = std::inner_product(Dvec.begin(), Dvec.end(), Tvec.begin(),0);

        alpha = dt0 / dt;

        std::transform (Dvec.begin(), Dvec.end(), Tmpvec.begin(),  std::multiplies<double>(),alpha);

        std::transform (Xvec.begin(), Xvec.end(), Tmpvec.begin(), TVec.begin(),   std::plus<double>());

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

    return XVec;

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
        std::cout << 'invalid number of argument.. Program exiting..';
        exit(EXIT_FAILURE);
    }

    nx = StringTo<int>(argv[5]);
    ny = StringTo<int>(argv[6]);
    iter = StringTo<int>(argv[7]);
    error = StringTo<int>(argv[8]);
    proc = StringTo<int>(argv[3]);

    int nnx = nx+1;
    int nny = ny+1;

    int totdim = (nx-1)*(ny-1);

    FDGRID fGrid = new FDGRID (nnx,nny);

    std::vector* xsol = callCG(fGrid,iter,error,proc);

    for (int i= 0; i< totdim; i++ )
        std::cout << xsol[i] << ' ';


    return 0;

}
