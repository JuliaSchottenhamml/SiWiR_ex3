#ifndef GRIDCAPTAIN_H
#define GRIDCAPTAIN_H

#endif // GRIDCAPTAIN_H

#include <iostream>
#include "immintrin.h"
#include "FDGRID.h"
#include <memory>

#define BLOCKSIZE 4

#define LD 6

class GridCaptain{


    FDGRID fdgrid;



public:

    double * worksheet;

    int hor,ver;
    int m,n;

    int const proc, dimM , dimN, blocklenX, blocklenY;

   GridCaptain(int const procn, FDGRID& fdg)
   {
       proc = procn;


       fdgrid = fdg;
       dimM= fdg.getDimM();
       dimN= fdg.getDimN();
       blocksz = BLOCKSIZE;

       worksheet = new double[proc*LD];
       m = 0;
       n= 0;

       if(procn%4==0)
           hor = 4;
       else if (procn%5==0)
           hor = 5;
       else if (procn%6==0)
           hor = 6;
       else if (procn%7==0)
           hor = 7;
        
        
        ver = procn/hor;
        if(dimN%hor == 0)
        blocklenX = dimN/hor;
        else
        blocklenX = dimN/hor+1;
        
        blocklenY = dimM/ver;
        assignWork();
   }

   inline void assignWork()
   {

       int remN = dimN%hor;
       int remM = dimM%ver;
       int bN = dimN/hor;
       int bM = dimM/ver;              
       int m =0;
              int n =0;
              int p =0;
       for (int i = 0; i < ver; i++)
         {
            int k = (i+1)*hor;
            int l = i*hor;
            n = 0;
           for (int j = 0; j < hor; j++)
            {
                int q = i*hor*LD+j*LD;
                worksheet[q]=p;
                worksheet[q+2]=i+m;
                worksheet[q+3]=n;
                worksheet[q+4]=worksheet[q+2]+p-1;
            

                /*if(i*ver+j+1 < k)
                    worksheet[i*proc+6]=i*ver+j+1;
                else
                    worksheet[i*proc+6]=-1;

                if(i*ver+j-1 < l)
                    worksheet[i*proc+7]=i*ver+j-1;
                else
                    worksheet[i*proc+7]=i*ver+j-1;

                if(i*ver+j+hor > dimM)
                    worksheet[i*proc+8]=i*ver+j+hor;
                else
                    worksheet[i*proc+8]=-1;

                    worksheet[i*proc+9]=i*ver+j-hor;*/                   
                    
                if(j<remN)
                {
                  worksheet[q+1]=bn+1;  
                  n + = bn+1;
                } 
                else
                {
                worksheet[q+1]=bn;  
                n + = bn;
                }
               worksheet[q+5]=n-1;
            }
            
             if(i<remM)
                {
                  p =bm+1;
                  m + = p; 
                } 
                else
                {
                  p =bM; 
                  m + = p; 
                } 
         }
   }
