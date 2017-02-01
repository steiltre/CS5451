#include <pthread.h> /* remember to include the pthread header */
#include <stdio.h>
#include <stdlib.h>

/* Define the number of threads to create */
#define NUM_THREADS 20

static void * hello(
            void * id)
{
      int myid;
        myid = * ((int *) id);
          printf("Hello World from thread %d\n", myid);

            /* terminate the thread */
            pthread_exit(NULL);
}


int main(
            int argc,
                char ** argv)
{
      int worker, val;
        int * th_ids;
          pthread_t threads[NUM_THREADS];

            th_ids = (int *) malloc(NUM_THREADS * sizeof(int));

              /* Create the threads and have them start executing the function hello() */
              for(worker = 0; worker < NUM_THREADS; worker++) {
                      th_ids[worker] = worker;
                          val = pthread_create(&threads[worker], NULL, hello,
                                          (void *) &th_ids[worker]);

                              /* Check for successful thread creation */
                              if(val != 0) {
                                        printf("Error creating thread #%d, val = %d\n", worker, val);
                                              exit(-1);
                                                  }
                                }

                free(th_ids);

                  pthread_exit(NULL);
}
