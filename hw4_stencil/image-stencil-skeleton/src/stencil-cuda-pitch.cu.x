
#define THREADSPERBLOCK 1024
#define WARPLENGTH 32

#include <stdio.h>
#include <math.h>
#include "cuda_profiler_api.h"

extern "C"
{
#include "image.h"
#include "stencil.h"
}

__device__ __constant__ float d_stencil[9];

__global__ void cuda_apply_stencil(
        float * input,
        float * output,
        size_t pitch,
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
        tile[ threadIdx.y * (blockDim.x+2) + threadIdx.x ] = *(float*) ( (char*)input + (row-1)*pitch + (col-1)*sizeof(float) ) ;
    }

    if (threadIdx.x < 2 && row > 0 && row < height+1 && col + blockDim.x < width) {
        tile[ threadIdx.y * (blockDim.x+2) + threadIdx.x + blockDim.x ] = *(float*) ( (char*)input + (row-1)*pitch + (col+blockDim.x-1)*sizeof(float) ) ;
    }

    if (threadIdx.y < 2 && col > 0 && col < width+1 && row + blockDim.y < height) {
        tile[ (threadIdx.y + blockDim.y) * (blockDim.x+2) + threadIdx.x ] = *(float*) ( (char*)input + (row+blockDim.y-1)*pitch + (col-1)*sizeof(float) );
    }

    if (threadIdx.x < 2 && threadIdx.y < 2 && col + blockDim.x < width && row + blockDim.y < height) {
        tile[ (threadIdx.y + blockDim.y) * (blockDim.x+2) + threadIdx.x + blockDim.x ] = *(float*) ( (char*)input + (row+blockDim.y-1)*pitch + (col+blockDim.x-1)*sizeof(float) );
    }

    __syncthreads();

    if (row >= 0 && col >= 0 && row < height && col < width) {
        if (row == 0 || col == 0 || row == height-1 || col == width-1) {
            //temp = tile[ (threadIdx.y+1) * (blockDim.x+2) + threadIdx.x+1 ];
            temp = *(float*) ( (char*)input + row*pitch + col*sizeof(float) );
        }
        else {
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    temp += d_stencil[i*3 + j] * tile[ (threadIdx.y+i)*(blockDim.x+2) + threadIdx.x+j ];
                }
            }
        }

        *(float*) ( (char*)output + row*pitch + col*sizeof(float) ) = temp;
    }

}

image_t * stencil_cuda(
    image_t const * const input,
    float stencil[3][3],
    int const num_times)
{
  image_t * output = image_alloc(input->width, input->height);

  cudaProfilerStart();

  //int arr_size = input->width*input->height * sizeof(float);

  float *d_redout, *d_greenout, *d_blueout, *d_redin, *d_greenin, *d_bluein;

  size_t pitch;  /* Pitch used for padding rows of array for coalescing memory accesses */

  /*
  int width = 35;
  int height = 38;
  float *test = (float*) malloc( width*height*sizeof(float) );
  float *testout = (float*) malloc( width*height*sizeof(float) );
  for (int i=0; i<height; i++) {
      for (int j=0; j<width; j++) {
          test[ i * width + j ] = (float) i*width+j;
      }
  }

  float *d_test;
  float *d_testout;
  cudaMalloc( (void**) &d_test, width*height*sizeof(float) );
  cudaMalloc( (void**) &d_testout, width*height*sizeof(float) );
  cudaMemcpy(d_test, test, height*width*sizeof(float), cudaMemcpyHostToDevice);
  */

  cudaMallocPitch( (void**) &d_redout, &pitch, input->width*sizeof(float), input->height );
  cudaMallocPitch( (void**) &d_greenout, &pitch, input->width*sizeof(float), input->height );
  cudaMallocPitch( (void**) &d_blueout, &pitch, input->width*sizeof(float), input->height );
  cudaMallocPitch( (void**) &d_redin, &pitch, input->width*sizeof(float), input->height );
  cudaMallocPitch( (void**) &d_greenin, &pitch, input->width*sizeof(float), input->height );
  cudaMallocPitch( (void**) &d_bluein, &pitch, input->width*sizeof(float), input->height );

  /*
  cudaMallocPitch( (void**) &d_greenout, arr_size );
  cudaMallocPitch( (void**) &d_blueout, arr_size );
  cudaMallocPitch( (void**) &d_redin, arr_size );
  cudaMallocPitch( (void**) &d_greenin, arr_size );
  cudaMallocPitch( (void**) &d_bluein, arr_size );
  */

  cudaMemcpyToSymbol(d_stencil, stencil, 9*sizeof(float));

  cudaMemcpy2D( d_redin, pitch, input->red, input->width*sizeof(float), input->width*sizeof(float), input->height, cudaMemcpyHostToDevice);
  cudaMemcpy2D( d_greenin, pitch, input->green, input->width*sizeof(float), input->width*sizeof(float), input->height, cudaMemcpyHostToDevice);
  cudaMemcpy2D( d_bluein, pitch, input->blue, input->width*sizeof(float), input->width*sizeof(float), input->height, cudaMemcpyHostToDevice);

  /*
  cudaMemcpy(d_redin, input->red, arr_size, cudaMemcpyHostToDevice);
  cudaMemcpy(d_greenin, input->green, arr_size, cudaMemcpyHostToDevice);
  cudaMemcpy(d_bluein, input->blue, arr_size, cudaMemcpyHostToDevice);
  */

  /* Determine dimensions of blocks and grid */
  /* Want blocks that are roughly square but also have length that is a multiple of WARPLENGTH for coalescing */
  int block_len =  ( (int) ceil( sqrt( float(THREADSPERBLOCK) ) / WARPLENGTH ) ) * WARPLENGTH;
  dim3 block( block_len, THREADSPERBLOCK / block_len );
  dim3 grid( ceil( input->width/float(block.x) ), ceil( input->height/float(block.y) ) );

  /*
  dim3 grid2( ceil( width/float(block.x) ), ceil( height/float(block.y) ) );

  cuda_apply_stencil<<<grid2, block, (block.x+2) * (block.y+2)*sizeof(float)>>>(d_test, d_testout, width, height);
  cudaMemcpy(testout, d_testout, height*width*sizeof(float), cudaMemcpyDeviceToHost);

  for (int i=0; i<height; i++) {
      for (int j=0; j<width; j++) {
            if ( test[i*width+j] != testout[i*width+j] )
                printf("Mismatch at: (%i %i) original: %0.03f output: %0.03f\n", i, j, test[i*width+j], testout[i*width+j] );
      }
  }
  */


  for (int i=0; i < num_times; ++i) {
      /* Apply stencil to each channel separately. */
      cuda_apply_stencil<<<grid, block, (block.x+2)*(block.y+2)*sizeof(float)>>>(d_redin, d_redout, pitch, input->width, input->height);
      cuda_apply_stencil<<<grid, block, (block.x+2)*(block.y+2)*sizeof(float)>>>(d_greenin, d_greenout, pitch, input->width, input->height);
      cuda_apply_stencil<<<grid, block, (block.x+2)*(block.y+2)*sizeof(float)>>>(d_bluein, d_blueout, pitch, input->width, input->height);
  }

  cudaMemcpy2D( output->red, output->width*sizeof(float), d_redout, pitch, output->width*sizeof(float), output->height, cudaMemcpyDeviceToHost);
  cudaMemcpy2D( output->green, output->width*sizeof(float), d_greenout, pitch, output->width*sizeof(float), output->height, cudaMemcpyDeviceToHost);
  cudaMemcpy2D( output->blue, output->width*sizeof(float), d_blueout, pitch, output->width*sizeof(float), output->height, cudaMemcpyDeviceToHost);

  /*
  cudaMemcpy(output->red, d_redout, arr_size, cudaMemcpyDeviceToHost);
  cudaMemcpy(output->green, d_greenout, arr_size, cudaMemcpyDeviceToHost);
  cudaMemcpy(output->blue, d_blueout, arr_size, cudaMemcpyDeviceToHost);
  */

  cudaDeviceSynchronize();
  cudaProfilerStop();

  return output;
}


