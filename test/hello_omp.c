#include <omp.h> /* Remember to include the OpenMP Header */
#include <stdio.h>

/* Define the number of threads
 *  * This can be done two ways:
 *   * 1: Define the number of threads within the program by using the
 *    * function call "omp_set_num_threads()" as is done in this program
 *     * 2: Set the environment variable OMP_NUM_THREADS, e.g.
 *      * > setenv OMP_NUM_THREADS 15
 *       * (To use this method, you must remove the function call
 *        * omp_set_num_threads() from below)
 *         */
#define NUM_THREADS 15

int main(
    int argc,
    char ** argv)
{
    int myid;

    /* set the number of threads to use */
    omp_set_num_threads(NUM_THREADS);

    /* Define a parallel segment within the code
     * A total of NUM_THREADS threads will be created and execute
     * the following code in parallel
     * myid is declared as a private variable (each thread will have
     * access to a unique memory address for myid, versus a public variable where
     * each thread would point to the same memory address for myid)
     */
    #pragma omp parallel private(myid)
    {
      myid = omp_get_thread_num();

      printf("Hello World from thread %d\n", myid);
    }

    return 0;
}
