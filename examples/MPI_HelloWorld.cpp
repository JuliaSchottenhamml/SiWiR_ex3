//*************************************************************************************************
// MPI_HelloWorld.cpp
//
// This program demonstrates the synchronization behavior of parallel MPI processes.
// Only process 0 prints the "Hello World" message, whereas all processes print their own rank.
// Run the example several times and note the different ordering of the output messages.
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

   // 'Hello World' output on process 0
   if( rank == 0 )
      std::cout << "Hello World!" << std::endl;

   // Output of the process rank for each process
   std::cout << "I am CPU " << rank << " of " << size << " CPUs" << std::endl;

   // MPI finalization
   MPI_Finalize();

   return 0;
}
