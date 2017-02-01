#include <stdio.h>
#include <pthread.h>

void *PrintHello(void * t)
{
    int i;
    i = (int)t;
    printf("Hello world from thread %d \n", t);
}

int main(int argc, char *argv[])
{
    //int num_threads = (int)argv[0];
    int num_threads = 4;
    int rc;
    int t;

    pthread_t threads[num_threads];

    for (t=0; t<num_threads; t++)
    {
        printf( "In main: creating thread %ld\n", t);
        rc = pthread_create(&threads[t], NULL, PrintHello, (void *)t);
        if (rc)
        {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
        }
    }
    pthread_exit(NULL);
}
