#include <cstdio>
#include <cstdlib>
#include <vector>

__global__ void bucket_sort (int *key. int *bucket, int n){
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i>=n) return;
  
  atomicAdd(bucket + key[i], 1);
  __syncthreads();

  for (int j=0, k=0, j<=i; k++){
    key[i] = k;
    j += bucket[k];
  }
}

int main() {
  int n = 50;
  int range = 5;
  //std::vector<int> key(n);
  int* key;
  cudaMallocManaged(&key, n*sizeof(int));
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");

  int* bucket;
  cudaMallocManaged(&bucket, range*sizeof(int));
  int M =32;

  bucket_sort<<<(n+M-1)/M, M, range>>>(key, bucket, n);
  cudaDeviceSynchronize();

/*
  std::vector<int> bucket(range); 
  for (int i=0; i<range; i++) {
    bucket[i] = 0;
  }
  for (int i=0; i<n; i++) {
    bucket[key[i]]++;
  }
  for (int i=0, j=0; i<range; i++) {
    for (; bucket[i]>0; bucket[i]--) {
      key[j++] = i;
    }
  }
*/

  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
}
