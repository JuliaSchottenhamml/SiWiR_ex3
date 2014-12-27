#ifndef GRIDCAPTAIN_H
#define GRIDCAPTAIN_H

#endif // GRIDCAPTAIN_H

#include <iostream>
#include "immintrin.h"
#include "FDGRID.h"
#include <memory>

#define BLOCKSIZE 4

#define LD 3

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
           ver = 4;
       else if (procn%5==0)
           ver = 5;
       else if (procn%6==0)
           ver = 6;
       else if (procn%7==0)
           ver = 7;
        
       assignWork();
   }

   inline void assignWork()
   {     
       int remM = dimM%ver;
       int bn = dimM/ver;              
       int m =0;
        int n =0;
        int p =0;
        int s = 0;
       for (int i = 0; i < ver; i++)
         {
           
                int q = i*LD;
                //worksheet[q+1]=dimN;
                worksheet[q]=n;
                //worksheet[q+3]=0;
                if(i<remM)
                {
                  m=bn+1;  
                  worksheet[q+1]=m; 
                  n + = m;
                  p= m*dimN;
                } 
                else
                {
                worksheet[q+1]=bn;  
                n + = bn;
                p = m*dimN;
                }
                worksheet[q+2]=s;
                s + = p;                 
               //worksheet[q+5]=dimN;
            }           
            
         }
   }
