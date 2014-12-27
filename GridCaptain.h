#ifndef GRIDCAPTAIN_H
#define GRIDCAPTAIN_H

#endif // GRIDCAPTAIN_H

#include <iostream>
#include "immintrin.h"
#include "FdGrid.h"
#include <memory>

#define BLOCKSIZE 4

#define LD1 3

class GridCaptain{


    FdGrid fdgrid;



public:

    double * worksheet;

    int hor,ver;
    //int m,n;

    int proc, dimM , dimN, blocklenX, blocklenY;

   GridCaptain(int const procn, FdGrid fdg)
   {
       proc = procn;


       fdgrid = fdg;
       dimM= fdg.getDimM();
       dimN= fdg.getDimN();
       
       worksheet = new double[proc*LD1];

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
           
                int q = i*LD1;
                //worksheet[q+1]=dimN;
                worksheet[q]=n;
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
                p = m*dimN;
                }
                worksheet[q+2]=s;
                s += p;                 
               //worksheet[q+5]=dimN;
            }           
            
         }
   };
