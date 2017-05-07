
#define THREADSPERBLOCK 512
#define WARPLENGTH 32

#include <stdio.h>
#include <math.h>

extern "C"
{
#include "image.h"
#include "stencil.h"
}

__device__ __constant__ float d_stencil[9];

__global__ void cuda_apply_stencil(
        float * input,
        float * output,
        int width,
        int height)
{
    extern __shared__ float tile[];  /* Tile needed by block */

    float temp;
    temp = 0;

    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    /* Fetch tile values from global memory and place in shared memory (need extra column and row around boundary for applying stencil) */
    if (row > 0 && col > 0 && row < height+1 && col < width+1) {
        tile[ threadIdx.y * (blockDim.x+2) + threadIdx.x ] = input[ (row-1)*width + (col-1) ];
    }

    if (threadIdx.x < 2 && row > 0 && row < height+1 && col + blockDim.x < width) {
        tile[ threadIdx.y * (blockDim.x+2) + threadIdx.x + blockDim.x ] = input[ (row-1)*width + col + blockDim.x - 1 ];
    }

    if (threadIdx.y < 2 && col > 0 && col < width+1 && row + blockDim.y < height) {
        tile[ (threadIdx.y + blockDim.y) * (blockDim.x+2) + threadIdx.x ] = input[ (row + blockDim.y - 1)*width + (col-1) ];
    }

    if (threadIdx.x < 2 && threadIdx.y < 2 && col + blockDim.x < width && row + blockDim.y < height) {
        tile[ (threadIdx.y + blockDim.y) * (blockDim.x+2) + threadIdx.x + blockDim.x ] = input[ (row + blockDim.y - 1)*width + col + blockDim.x - 1 ];
    }

    __syncthreads();

    if (row >= 0 && col >= 0 && row < height && col < width) {
        if (row == 0 || col == 0 || row == height-1 || col == width-1) {
            temp = tile[ (threadIdx.y+1) * (blockDim.x+2) + threadIdx.x+1 ];
            //temp = input[ row * width + col ];
        }
        else {
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    temp += d_stencil[i*3 + j] * tile[ (threadIdx.y+i)*(blockDim.x+2) + threadIdx.x+j ];
                }
            }
        }

        output[ row*width + col ] = temp;
    }

}

image_t * stencil_cuda(
    image_t const * const input,
    float stencil[3][3],
    int const num_times)
{
  image_t * output = image_alloc(input->width, input->height);

  int arr_size = input->width*input->height;

  float *d_redout, *d_greenout, *d_blueout, *d_redin, *d_greenin, *d_bluein;

  cudaMalloc( (void**) &d_redout, 3*arr_size * sizeof(*d_redout) );
  //cudaMalloc( (void**) &d_greenout, arr_size );
  //cudaMalloc( (void**) &d_blueout, arr_size );
  cudaMalloc( (void**) &d_redin, 3*arr_size  * sizeof(*d_redin) );
  //cudaMalloc( (void**) &d_greenin, arr_size );
  //cudaMalloc( (void**) &d_bluein, arr_size );

  d_greenout = d_redout + arr_size;
  d_blueout = d_redout + 2*arr_size;
  d_greenin = d_redin + arr_size;
  d_bluein = d_redin + 2*arr_size;

  cudaMemcpyToSymbol(d_stencil, stencil, 9*sizeof(float));

  cudaMemcpy(d_redin, input->red, 3*arr_size * sizeof(*d_redin) , cudaMemcpyHostToDevice);
  //cudaMemcpy(d_greenin, input->green, arr_size, cudaMemcpyHostToDevice);
  //cudaMemcpy(d_bluein, input->blue, arr_size, cudaMemcpyHostToDevice);

  /* Determine dimensions of blocks and grid */
  /* Want blocks that are roughly square but also have length that is a multiple of WARPLENGTH for coalescing */
  int block_len =  ( (int) ceil( sqrt( float(THREADSPERBLOCK) ) / WARPLENGTH ) ) * WARPLENGTH;
  dim3 block( block_len, THREADSPERBLOCK / block_len );
  dim3 grid( ceil( input->width/float(block.x) ), ceil( input->height/float(block.y) ) );

  for (int i=0; i < num_times; ++i) {
      /* Apply stencil to each channel separately. */
      cuda_apply_stencil<<<grid, block, (block.x+2)*(block.y+2)*sizeof(float)>>>(d_redin, d_redout, input->width, input->height);
      cuda_apply_stencil<<<grid, block, (block.x+2)*(block.y+2)*sizeof(float)>>>(d_greenin, d_greenout, input->width, input->height);
      cuda_apply_stencil<<<grid, block, (block.x+2)*(block.y+2)*sizeof(float)>>>(d_bluein, d_blueout, input->width, input->height);
  }

  cudaMemcpy(output->red, d_redout, 3 * arr_size * sizeof(*d_redout), cudaMemcpyDeviceToHost);
  //cudaMemcpy(output->green, d_greenout, arr_size, cudaMemcpyDeviceToHost);
  //cudaMemcpy(output->blue, d_blueout, arr_size, cudaMemcpyDeviceToHost);

  cudaFree(d_redout);
  //cudaFree(d_greenout);
  //cudaFree(d_blueout);
  cudaFree(d_redin);
  //cudaFree(d_greenin);
  //cudaFree(d_bluein);

  return output;
}


