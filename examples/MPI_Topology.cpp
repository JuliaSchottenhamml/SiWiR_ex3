//*************************************************************************************************
// MPI_Topology.cpp
//
// This MPI example demonstrates the setup of a MPI topology. The topology will consist of
// four processes arranged in a 2x2 fashion. One possible process arrangement might look like this:
//
//      | 0  1 |
//      | 2  3 |
//
// After the setup of the MPI topology, each process is supposed to send his own rank to the
// neighboring processes. The neighbors will receive the messages and compare the received
// ranks to the expected ranks.
//*************************************************************************************************


// Includes
#include <mpi.h>
#include <iostream>


// Global constants
const int SIZE  = 4;
const int UP    = 0;
const int DOWN  = 1;
const int LEFT  = 2;
const int RIGHT = 3;


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

   // Creating the Cartesian topology and sending the own rank
   // to the neighboring CPUs
   if( size == SIZE )
   {
      // Definition and initialization of the dimension array 'dims'
      // In this case, we will construct a 2D-grid with 2 processes in x-
      // and 2 processes in y-direction
      int dims[2] = {2,2};
      
      // Definition and initialization of the behavior of the Cartesian topology
      int periods[2] = {0,0}; // not periodic!
      const int reorder = 1;  // allow reordering of process ranks
   
      // The new MPI communicator for the Cartesian topology
      MPI_Comm cartcomm( MPI_COMM_NULL );
      
      // The rank/number of this process (within the Cartesian topology)
      int cartrank(0);
      
      // Coordinates of this process within the Cartesian topology
      int coords[2] = {0,0};
      
      // Definition of the neighbors array 'nrbs'
      // Will contain the ranks of the neighbors in the order UP, DOWN, LEFT and RIGHT
      int nbrs[4] = {0,0,0,0};
   
      // Creating the Cartesian topology:
      //  - Creating a 2D grid with 2 processes in x- and 2 processes in y-direction
      //  - Determining the coordinates of the processes
      //  - Determining the neighbors in UP, DOWN, LEFT and RIGHT direction
      // ----------------------------------------------------------------
      MPI_Cart_create( MPI_COMM_WORLD, 2, dims, periods, reorder, &cartcomm );
      MPI_Comm_rank( cartcomm, &cartrank );
      MPI_Cart_coords( cartcomm, cartrank, 2, coords );
      MPI_Cart_shift( cartcomm, 0, 1, &nbrs[LEFT], &nbrs[RIGHT] );
      MPI_Cart_shift( cartcomm, 1, 1, &nbrs[DOWN], &nbrs[UP] );
      // ----------------------------------------------------------------

      // Setting the message tag (offers the possibility to give
      // each message an explicit number)
      const int tag = 1; // here: same tag for each message

      // Definition and initialization of the receive buffer 'inbuf', where
      // the received messages are stored and which will in this case
      // contain the ranks of the neighboring CPUs
      int inbuf[4] = { MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL };
      
      // Sending the own rank to all four neighboring CPUs
      // and receiving the rank of all four neighbors
      // ----------------------------------------------------------------
      MPI_Request reqs[8];
      MPI_Status stats[8];
      
      for( int i = 0; i < 4; ++i )
      {         
         MPI_Isend( &cartrank, 1, MPI_INT, nbrs[i], tag, cartcomm, &reqs[i]   );
         MPI_Irecv( &inbuf[i], 1, MPI_INT, nbrs[i], tag, cartcomm, &reqs[i+4] );
      }

      MPI_Waitall( 8, reqs, stats );
      // ----------------------------------------------------------------
 
      // output
      for( int i = 0; i < 4; ++i )
      {
         if( rank == i )
         {
            std::cout << "------------------------------------------------------------\n"
                      << "rank (MPI_COMM_WORLD):             " << rank << "\n"
                      << "rank (Cartesian topology):         " << cartrank << "\n"
                      << "Cartesian coordinates:             ( " << coords[0] << ", " << coords[1] << " )\n"
                      << "neighbors (x-direction, expected): " << nbrs[LEFT] << " (left), " << nbrs[RIGHT] << " (right)\n"
                      << "neighbors (y-direction, expected): " << nbrs[DOWN] << " (down), " << nbrs[UP] << " (up)\n"
                      << "neighbors (x-direction, received): " << inbuf[LEFT] << " (left), " << inbuf[RIGHT] << " (right)\n"
                      << "neighbors (y-direction, received): " << inbuf[DOWN] << " (down), " << inbuf[UP] << " (up)"                      
                      << std::endl;
         }
         MPI_Barrier( MPI_COMM_WORLD );
      }
      
      // Everyone sends its Cartesian coordinates and Cartesian rank to process 0
      // so that process 0 can output the entire process topology with Cartesian ranks
      
      int topology[4] = {0,0,0,0};
      
      for( int i = 0; i < SIZE; ++i )
      {
         if( rank == 0 )
         {
            if( i == 0 )
            {
               topology[ coords[1] * 2 + coords[0] ] = cartrank;
            }
            else
            {
               MPI_Status status;
               int othercoords[2] = {0,0};
               MPI_Recv( othercoords,                                      2, MPI_INT, i, 0, MPI_COMM_WORLD, &status );
               MPI_Recv( &topology[ othercoords[1] * 2 + othercoords[0] ], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status );
            }
         }
         else if( rank == i )
         {
            MPI_Send( coords,    2, MPI_INT, 0, 0, MPI_COMM_WORLD );
            MPI_Send( &cartrank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD );
         }
      }
      
      if( rank == 0 )
      {
         std::cout << "------------------------------------------------------------\n"
                   << "topology (Cartesian ranks are shown):\n\n"
                   << "       | " << topology[2] << " " << topology[3] << " |\n"
                   << "  y ^  | " << topology[0] << " " << topology[1] << " |\n"
                   << "    |\n"
                   << "    ---> x\n"
                   << std::endl;
      }
   }
   else
   {
      std::cout << "You must specify " << SIZE << " processes!\nTerminating..." << std::endl;
   }

   // MPI finalizations
   // ----------------------------------------------------------------
   MPI_Finalize();
   // ----------------------------------------------------------------
   
   return 0;
}
