#include <cstdio>
#include <cstdlib>
#include <cmath>

int main() {
  const int N = 8;
  double x[N], y[N], m[N], fx[N], fy[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
  }
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
    for(int j=0; j<N; j++) {
      if(i != j) {
        double rx = x[i] - x[j];
        double ry = y[i] - y[j];
        double r = std::sqrt(rx * rx + ry * ry);
        fx[i] -= rx * m[j] / (r * r * r);
        fy[i] -= ry * m[j] / (r * r * r);
      }
    }
    printf("%d %g %g\n",i,fx[i],fy[i]);
  }*/
}
