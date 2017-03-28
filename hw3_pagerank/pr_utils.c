#include "pr_utils.h"


int GetNumChunks(
        int const chunk_size,
        int const num_items)
{
    int num_chunks = num_items / chunk_size;
    if (num_items % chunk_size > 0)
    {
        ++num_chunks;
    }
    return num_chunks;
}


int GetChunkSize(
        int const chunk_id,
        int const chunk_size,
        int const num_items)
{
    int stop = chunk_size;
    if(stop + (chunk_id * chunk_size) > num_items) {
        stop = num_items % chunk_size;
    }
    return stop;
}

pr_int binary_search(
    pr_int * arr,
    pr_int arr_size,
    pr_int val)
{
  pr_int start, end, mid;

  start = 0;
  end = arr_size-1;
  do {
    if (start > end)
      return -1;

    mid = (start + end) / 2;

    if (arr[mid] < val)
      start = mid+1;
    else if (arr[mid] > val)
      end = mid-1;
    else if (arr[mid] == val)
      return mid;
  } while(1);
}
