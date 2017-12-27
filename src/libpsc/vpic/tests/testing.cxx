
#include "testing.h"

#include "psc_vpic_bits.h"

#include "src/vpic/vpic.h"

MPI_Comm psc_comm_world;
int psc_world_rank;
int psc_world_size;

void testing_init(int *argc, char ***argv)
{
  boot_checkpt(argc, argv);
  serial.boot(argc, argv);
  thread.boot(argc, argv);
    
  boot_mp(argc, argv);
  //MPI_Init(argc, argv);

  MPI_Comm_dup(MPI_COMM_WORLD, &psc_comm_world);
  MPI_Comm_rank(psc_comm_world, &psc_world_rank);
  MPI_Comm_size(psc_comm_world, &psc_world_size);
}

