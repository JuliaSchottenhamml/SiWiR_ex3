//*************************************************************************************************
// MPI_Parameters.cpp
//
// This MPI example demonstrates how the command line parameters are used in a MPI parallel
// program.
//*************************************************************************************************


// Includes
#include <mpi.h>
#include <iostream>


int main( int argc, char* argv[] )
{
   int size(0); // The total number of processes
   int rank(0); // The rank/number of this process

   // MPI initialization
   MPI_Init( &argc, &argv );

   // Determining the number of CPUs and the rank of this process
   MPI_Comm_size( MPI_COMM_WORLD, &size );
   MPI_Comm_rank( MPI_COMM_WORLD, &rank );

   // Output of the command line arguments (only on process 0!)
   if( rank == 0 ) {
      for( int i = 0; i < argc; ++i )
         std::cout << "Argument " << i << ": " << argv[i] << "\n";
      std::cout << std::endl;;
   }

   // MPI finalization
   MPI_Finalize();

   return 0;
}
