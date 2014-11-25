//*************************************************************************************************
// MPI_Vector.cpp
//
// ThIS MPI example demonstrates the setup of a new MPI data type "vector". This new data type
// allows "vector"-messages that may be columns of a matrix.
//*************************************************************************************************


// Includes
#include <mpi.h>
#include <iostream>


// Global constant
const int SIZE = 4;


int main( int argc, char* argv[] )
{
   int size(0); // The total number of processes
   int rank(0); // The rank/number of this process (within MPI_COMM_WORLD)

   // Initialization of MPI
   // ----------------------------------------------------------------
   MPI_Init( &argc, &argv );
   // ----------------------------------------------------------------

   // Determining the number of CPUs and the rank of this process
   // ----------------------------------------------------------------
   MPI_Comm_size( MPI_COMM_WORLD, &size );
   MPI_Comm_rank( MPI_COMM_WORLD, &rank );
   // ----------------------------------------------------------------   
   
   // Definition and initialization of the data field a
   double a[ SIZE * SIZE ] = {  1.0,  2.0,  3.0,  4.0,
                                5.0,  6.0,  7.0,  8.0,
                                9.0, 10.0, 11.0, 12.0,
                               13.0, 14.0, 15.0, 16.0 };

   // Definition of the receive buffer b
   double b[SIZE] = { 0.0, 0.0, 0.0, 0.0 };
   
   // Creating a new derived data type
   MPI_Datatype columntype;   
   MPI_Type_vector( SIZE, 1, SIZE, MPI_DOUBLE, &columntype );
   MPI_Type_commit( &columntype );

   // Let the CPU 0 send a column of 'a' to each of the other 3 CPUs,
   // which receive this column in the receive array 'b'
   if( size == SIZE )
   {
      if( rank == 0 )
      {
         b[0] = a[0];
         b[1] = a[4];
         b[2] = a[8];
         b[3] = a[12];

         // ----------------------------------------------------------------
         for( int i = 1; i < 4; ++i )
            MPI_Send( &a[i], 1, columntype, i, 0, MPI_COMM_WORLD );
         // ----------------------------------------------------------------            
      }
      else
      {
         // ----------------------------------------------------------------
         MPI_Status status;
         MPI_Recv( b, SIZE, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status );
         // ----------------------------------------------------------------
      }

      // Output
      std::cout << "rank = " << rank << ": b = ( " << b[0] << ", " << b[1] << ", " << b[2] << ", " << b[3] << " )\n";
   }
   else
   {
      std::cout << "You must specify " << SIZE << " processes!\nTerminating..." << std::endl;
   }
   
   // Freeing derived data type
   MPI_Type_free( &columntype );

   // MPI finalizations
   // ----------------------------------------------------------------
   MPI_Finalize();
   // ----------------------------------------------------------------

   return 0;
}
