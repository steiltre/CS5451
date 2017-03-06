#include <mpi.h> /* Remember to include the MPI header */
#include <stdio.h>

int main(
    int argc,
    char ** argv)
{
  int num_procs, myid, name_len;
  char proc_name[MPI_MAX_PROCESSOR_NAME];

  /* Initialize MPI */
  MPI_Init(&argc, &argv);

  /* Obtain the number of processes */
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  /* Obtain the process id */
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  /* Obtain name of machine process is executing on */
  MPI_Get_processor_name(proc_name, &name_len);

  printf("Hello World from processor %d out of %d, executing on %s\n",
      myid, num_procs, proc_name);

  /* Last call to MPI (REQUIRED) */
  MPI_Finalize();

  return 0;
}
