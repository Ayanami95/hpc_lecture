#include <cstdio>
#include <cstdlib>
#include <vector>

__global__ void bucket_sort (int *key, int *bucket, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i>=n) return;
  
    atomicAdd(bucket + key[i], i);
    __syncthreads();

    for (int j = 0, k = 0; j <= i; k++) {
        key[i] = k;
        j += bucket[k];
    }
}

int main() {
    int n = 50;
    int range = 5;

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



    for (int i=0; i<n; i++) {
        printf("%d ",key[i]);
    }
    printf("\n");

    cudaFree(key);
    cudaFree(bucket);
}
