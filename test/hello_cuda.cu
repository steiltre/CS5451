/* Borrowed from http://computer-graphics.se/hello-world-for-cuda.html
 * This program takes the string "Hello ", prints it, then passes it to CUDA
 * with an array * of offsets. Then the offsets are added in parallel to
 * produce the string "World!"
 * By Ingemar Ragnemalm 2010
 */

#include <stdlib.h>
#include <stdio.h>

int const N = 16;
int const blocksize = 16;

__global__
void hello(
    char * const a,
    int const * const b)
{
  a[threadIdx.x] += b[threadIdx.x];
}

int main(
    int argc,
    char ** argv)
{
  char a[N] = "Hello \0\0\0\0\0\0";
  int b[N] = {15, 10, 6, 0, -11, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  char * ad;
  int * bd;
  int const csize = N * sizeof(char);
  int const isize = N * sizeof(int);

  printf("%s", a);

  /* Allocate space for a and b on the device */
  cudaMalloc((void**) &ad, csize );
  cudaMalloc((void**) &bd, isize );

  /* copy the contents of a and b to the device */
  cudaMemcpy(ad, a, csize, cudaMemcpyHostToDevice);
  cudaMemcpy(bd, b, isize, cudaMemcpyHostToDevice);

  dim3 dimBlock(blocksize, 1);
  dim3 dimGrid(1, 1);
  /* call the CUDA kernel */
  hello<<<dimGrid, dimBlock>>>(ad, bd);

  /* copy the result back to the host */
  cudaMemcpy(a, ad, csize, cudaMemcpyDeviceToHost);
  cudaFree(ad);

  printf("%s\n", a);

  return 0;
}
