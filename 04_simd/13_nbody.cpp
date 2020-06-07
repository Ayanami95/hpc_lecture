#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <immintrin.h>

int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
  }
<<<<<<< HEAD
  __m256 fxvec = _mm256_set1_ps(0);
  __m256 fyvec = _mm256_set1_ps(0);
  __m256 xvec = _mm256_load_ps(x);
  __m256 yvec = _mm256_load_ps(y);
  __m256 mvec = _mm256_load_ps(m);
  
  __m256 zero = _mm256_set1_ps(0);
  __m256 one = _mm256_set1_ps(1);
  
  for(int i=0; i<N; i++){
    __m256 xi = _mm256_set1_ps(x[i]);
    __m256 yi = _mm256_set1_ps(y[i]);
    __m256 rx = _mm256_sub_ps(xi, xvec);
    __m256 ry = _mm256_sub_ps(yi, yvec);
    __m256 rx2 = _mm256_mul_ps(rx, rx);
    __m256 ry2 = _mm256_mul_ps(ry, ry);
    __m256 r2 = _mm256_add_ps(rx2, ry2);
	  
    // To avoid the situation when r2 = 0
    __m256 mask = _mm256_cmp_ps(r2, zero, _CMP_GT_OQ);
    r2 = _mm256_belndv_ps(one, r2, mask);
	  
    __m256 r = _mm256_rsqrt_ps(r2);
    r2 = _mm256_mul_ps(r, r);
    __m256 r3 = _mm256_mul_ps(r2, r);
	  
    __m256 temp = _mm256_mul_ps(rx, mvec);
    temp = _mm256_mul_ps(temp, r3);
    fxvec = _mm256_sub_ps(fxvec, temp);
    temp = _mm256_mul_ps(ry, mvec);
    temp = _mm256_mul_ps(temp, r3);
    fyvec = _mm256_sub_ps(fyvec, temp);
  }

  _mm256_store_ps(fx, fxvec);
  _mm256_store_ps(fy, fyvec);

  for(int i=0; i<N; i++)
    printf("%d %g %g\n",i,fx[i],fy[i]);

  /*for(int i=0; i<N; i++) {
=======
  __m256 zero = _mm256_setzero_ps();
  for(int i=0; i<N; i+=8) {
    __m256 xi = _mm256_load_ps(x+i);
    __m256 yi = _mm256_load_ps(y+i);
    __m256 fxi = zero;
    __m256 fyi = zero;
>>>>>>> f743798ff25f63cf544466b630c34b35525ca76f
    for(int j=0; j<N; j++) {
      __m256 dx = _mm256_set1_ps(x[j]);
      __m256 dy = _mm256_set1_ps(y[j]);
      __m256 mj = _mm256_set1_ps(m[j]);
      __m256 r2 = zero;
      dx = _mm256_sub_ps(xi, dx);
      dy = _mm256_sub_ps(yi, dy);
      r2 = _mm256_fmadd_ps(dx, dx, r2);
      r2 = _mm256_fmadd_ps(dy, dy, r2);
      __m256 mask = _mm256_cmp_ps(r2, zero, _CMP_GT_OQ);
      __m256 invR = _mm256_rsqrt_ps(r2);
      invR = _mm256_blendv_ps(zero, invR, mask);
      mj = _mm256_mul_ps(mj, invR);
      invR = _mm256_mul_ps(invR, invR);
      mj = _mm256_mul_ps(mj, invR);
      fxi = _mm256_fmadd_ps(dx, mj, fxi);
      fyi = _mm256_fmadd_ps(dy, mj, fyi);
    }
<<<<<<< HEAD
    printf("%d %g %g\n",i,fx[i],fy[i]);
  }*/
=======
    _mm256_store_ps(fx+i, fxi);
    _mm256_store_ps(fy+i, fyi);
  }
  for(int i=0; i<N; i++)
    printf("%d %g %g\n",i,fx[i],fy[i]);
>>>>>>> f743798ff25f63cf544466b630c34b35525ca76f
}
