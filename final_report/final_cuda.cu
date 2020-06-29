#include <iostream>
#include <fstream>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "device_launch_parameters.h"

using namespace std;
const int nx = 41;
const int ny = 41;

__global__ void cal_b(float *b, float rho, float dt, float *u, float *v, float dx, float dy) {
	int m = blockIdx.x;
	int n = threadIdx.x;

	if (0 <= m&&m < nx - 2 && 0 <= n&& n < ny - 2) {
		float dudx = (u[(m + 1) * nx + n + 2] - u[(m + 1) * nx + n]) / (2.0 * dx);
		float dvdy = (v[(m + 2) * nx + n + 1] - v[(m)*nx + n + 1]) / (2.0 * dy);
		float dudy = (u[(m + 2) * nx + n + 1] - u[(m)*nx + n + 1]) / (2.0 * dy);
		float dvdx = (v[(m + 1) * nx + n + 2] - v[(m + 1) * nx + n]) / (2.0 * dx);
		b[(m + 1) * nx + n + 1] = rho * 1 / dt * (dudx + dvdy) - dudx * dudx - 2 * dudy * dvdx - dvdy * dvdy;
	}
}

__global__ void cal_p(float *p, float dx, float dy, float *b) {
	int m = blockIdx.x;
	int n = threadIdx.x;

	if (0 <= m&&m < nx - 2 && 0 <= n&& n < ny - 2) {
		p[(m + 1)* nx + n + 1] = ((p[(m + 1) * nx + n + 2] + p[(m + 1)* nx + n]) * dy * dy
			+ (p[(m + 2) * nx + n + 1] + p[(m)* nx + n + 1]) * dx * dx
			- b[(m + 1)* nx + n + 1] * dx * dx * dy * dy)
			/ (2 * (dx * dx + dy * dy));
	}

	if (0<= m && m < nx) {
		p[(m)*nx + ny - 1] = p[(m)*nx + ny - 2];
		p[(m)*nx + 0] = p[(m)*nx + 1];
	}
	if (0<= n && n < ny) {
		p[(0) * nx + n] = p[(1) * nx + n];
		p[(nx - 1) * nx + n] = 0;
	}
}
__global__ void cal_uv(float *u, float *v,float*p, float dx, float dy,float rho,float dt, float nu) {
	int m = blockIdx.x;
	int n = threadIdx.x;
	if (0 <= m&&m < nx - 2 && 0 <= n&& n < ny - 2) {
			u[(m + 1) * nx + n + 1] = (u[(m + 1) * nx + n + 1]
				- u[(m + 1) * nx + n + 1] * dt / dx * (u[(m + 1) * nx + n + 1] - u[(m + 1) * nx + n])
				- v[(m + 1) * nx + n + 1] * dt / dy * (u[(m + 1) * nx + n + 1] - u[(m)*nx + n + 1])
				- dt / (2 * rho * dx) * (p[(m + 1) * nx + n + 2] - p[(m + 1) * nx + n])
				+ nu * (dt / (dx * dx) * (u[(m + 1) * nx + n + 2] - 2 * u[(m + 1) * nx + n + 1] + u[(m + 1) * nx + n])
					+ dt / (dy * dy) * (u[(m + 2) * nx + n + 1] - 2 * u[(m + 1) * nx + n + 1] + u[(m)*nx + n + 1])));

			v[(m + 1) * nx + n + 1] = (v[(m + 1) * nx + n + 1]
				- u[(m + 1) * nx + n + 1] * dt / dx * (v[(m + 1) * nx + n + 1] - v[(m + 1) * nx + n])
				- v[(m + 1) * nx + n + 1] * dt / dy * (v[(m + 1) * nx + n + 1] - v[(m)*nx + n + 1])
				- dt / (2 * rho * dx) * (p[(m + 2) * nx + n + 1] - p[(m + 0) * nx + n + 1])
				+ nu * (dt / (dx * dx) * (v[(m + 1) * nx + n + 2] - 2 * v[(m + 1) * nx + n + 1] + v[(m + 1) * nx + n])
					+ dt / (dy * dy) * (v[(m + 2) * nx + n + 1] - 2 * v[(m + 1) * nx + n + 1] + v[(m)*nx + n + 1])));
		
	}
	//  u[(0,:] =  v[(0,:] = v[(-1,:] = 0  u[(-1,:] =1
	if (0<= n&&n < ny) {
		u[(0) * nx + n] = 0;
		u[(nx - 1) * nx + n] = 1;
		v[(0) * nx + n] = 0;
		v[(nx - 1) * nx + n] = 0;
	}

	//  u[(:,0] = u[(:,-1] = v[(:,0] = v[(:,-1] = 0
	if (0 <= m&&m < nx) {
		u[(m)*nx + 0] = 0;
		u[(m)*nx + ny - 1] = 0;
		v[(m)*nx + 0] = 0;
		v[(m)*nx + ny - 1] = 0;
	}
}


void cavity_flow(int nt, float *u, float *v, float dt, float dx, float dy, float *p, float rho, float nu, int nit) {

	float b[nx*ny] = { 0 };

	float*d_b, *d_u, *d_v,*d_p;
	int size = nx*ny * sizeof(float);
	cudaMalloc((void**)&d_b, size);
	cudaMalloc((void**)&d_u, size);
	cudaMalloc((void**)&d_v, size);
	cudaMalloc((void**)&d_p, size);
	dim3 dimBlock(ny);
	dim3 dimGrid(nx);

	cudaMemcpy(d_b, b, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_u, u, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_v, v, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_p, p, size, cudaMemcpyHostToDevice);


	for (int i = 0; i < nt; i++) {
		cal_b << <dimGrid, dimBlock >> > (d_b, rho, dt, d_u,d_v, dx, dy);
		for (int k = 0; k < nit; k++) {
			cal_p << <dimGrid, dimBlock >> > (d_p, dx, dy, d_b);
		}
		cal_uv << <dimGrid, dimBlock >> > (d_u, d_v,d_p, dx, dy, rho, dt, nu);
	}

	cudaMemcpy(p, d_p, size, cudaMemcpyDeviceToHost);
	cudaMemcpy(u, d_u, size, cudaMemcpyDeviceToHost);
	cudaMemcpy(v, d_v, size, cudaMemcpyDeviceToHost);

	cudaFree(d_b);
	cudaFree(d_u);
	cudaFree(d_v);
	cudaFree(d_p);
}

void save_data(float u[], float v[], float p[]) {
	ofstream fs("cavity_cpp_results.txt");

	// u
	fs << "u ";
	for (int m = 0; m < nx; m++) {
		for (int n = 0; n < ny; n++) {
			fs << u[m*nx + n] << "\t";
		}
		fs << "\n";
	}
	fs << "\nv";

	for (int m = 0; m < nx; m++) {
		for (int n = 0; n < ny; n++) {
			fs << v[m * nx + n] << "\t";
		}
		fs << "\n";
	}
	fs << "\np";

	for (int m = 0; m < nx; m++) {
		for (int n = 0; n < ny; n++) {
			fs << p[m * nx + n] << "\t";
		}
		fs << "\n";
	}
	fs << "\n";

	fs.close();
}

int main()
{
	time_t start, end;
	int nt = 700;
	int nit = 50;

	float dx = 2.0 / (nx - 1.0);
	float dy = 2.0 / (ny - 1.0);

	float rho = 1;
	float nu = 0.1;
	float dt = 0.001;

	float u[nx*ny] = { 0 };
	float v[nx*ny] = { 0 };
	float p[nx*ny] = { 0 };
	cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu, nit);
	save_data(u, v, p);

}