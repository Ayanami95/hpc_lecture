#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

using namespace std;

const int nx=41;   
const int ny=41;  
const int nt=100;  
const int nit=50;  
const int size = nx * ny;

const double dx=2.0/(nx-1);
const double dy=2.0/(ny-1);

const double rho=1;
const double nu=0.1;
const double dt=0.001;

void build_up_b(double b[], double u[],double v[]){
    for(int i = 1; i < nx - 1; i++){
        for(int j = 1; j < ny - 1; j++){
            b[j*nx+i] = (rho * (1.0/dt * ((u[j*nx + i+1] - u[j*nx + i-1]) 
            / (2 * dx) + (v[(j+1)*nx+i] - v[(j-1)*nx+i]) / (2 * dy)) - ((u[j*nx+i+1] - u[j*nx+i-1]) 
            / (2*dx)) * ((u[j*nx+i+1] - u[j*nx+i-1]) / (2*dx)) - 2 * ((u[(j+1)*nx+i] - u[(j-1)*nx+i]) 
            / (2*dy) * (v[j*nx+i+1] - v[j*nx + i-1]) / (2*dx)) - ((v[(j+1)*nx+i] - v[(j-1)*nx+i]) 
            / (2*dy)) * ((v[(j+1)*nx+i] - v[(j-1)*nx+i]) / (2*dy)) ));
        }
    }
}

void pressure_poisson(double p[],double pn[], double b[]){
   int i, j;
    for(int k = 0; k < nit ; k++){
        for(i = 0; i < nx ; i++){
            for(j = 0; j < ny ; j++){
                pn[j*nx+i] = p[j*nx+i];
            }
        }
        for(i = 1; i < nx - 1 ; i++){
            for(j = 1; j < ny - 1; j++){
                p[j*nx+i] = (((pn[j*nx+i+1] + pn[j*nx+i-1]) * dy * dy + (pn[(j+1)*nx+i] + pn[(j-1)*nx+i]) * dx * dx)
                / (2 * (dx * dx + dy * dy)) - dx * dx * dy * dy * b[j*nx+i] * rho / (2 * (dx *dx + dy * dy)));
            }
        }
        for(int m = 1 ; m < ny - 1 ; m++){
            p[m*nx+nx-1] = p[m*nx+nx-2]; 
            p[m*nx+0]    = p[m*nx+1];    
        }
        for(int n = 0 ; n < nx ; n++){
            p[0*nx+n] = p[1*nx+n]; 
            p[(nx-1)*nx+n] = 0.0;  
        }
    }
}

void cavity_flow(double u[], double un[], double v[], double vn[], double p[], double pn[], double b[]){
    for(int i = 0; i < nx; i++){
        for(int j = 0; j < ny; j++){
            un[j*nx+i] = u[j*nx+i];
            vn[j*nx+i] = v[j*nx+i];
        }
    }
    
    build_up_b(b, u, v);
    pressure_poisson(p, pn, b);
    
    for(int i = 1; i < nx - 1; i++){
        for(int j = 1; j < ny - 1; j++){
            u[j*nx+i] = (un[j*nx+i] - un[j*nx+i] * dt / dx * (un[j*nx+i] - un[j*nx+i-1]) - vn[j*nx+i] * dt 
            / dy * (un[j*nx+i] - un[(j-1)*nx+i]) - dt / (2 * rho * dx) * (p[j*nx+i+1] - p[j*nx+i-1]) 
            + nu * (dt / (dx * dx) * (un[j*nx+i+1] - 2*un[j*nx+i] + un[j*nx+i-1]) + dt / (dy * dy) * (un[(j+1)*nx+i] 
            - 2*un[j*nx+i] + un[(j-1)*nx+i])));

            v[j*nx+i] = ( vn[j*nx+i] - un[j*nx+i] * dt / dx * (vn[j*nx+i] - vn[j*nx+i-1]) - vn[j*nx+i] * dt 
            / dy * (vn[j*nx+i] - vn[(j-1)*nx+i]) - dt / (2 * rho * dy) * (p[(j+1)*nx+i] - p[(j-1)*nx+i]) 
            + nu * (dt / (dx * dx) * (vn[j*nx+i+1] - 2*vn[j*nx+i] + vn[j*nx+i-1]) + dt / (dy * dy) 
            * (vn[(j+1)*nx+i] - 2*vn[j*nx+i] + vn[(j-1)*nx+i])));
        }
    }
    
    for(int i = 0 ; i < nx; i++){
        u[0*nx+i] = 0.0;
        u[(nx-1)*nx+i] = 1.0;
        v[0*nx+i] = 0.0;
        v[(nx-1)*nx+i] = 0.0;
    }
    
    for(int j = 1; j < ny - 1; j++){
        u[j*nx+0] = 0.0;
        u[j*nx+nx-1] = 0.0;
        v[j*nx+0] = 0.0;
        v[j*nx+nx-1] = 0.0;
    }
}

void writeFile(double u[], double v[], double p[]) {
    ofstream fs("cavity_cpp_results_04.txt");

    fs << "u ";
    for (int i = 0; i < size; i++) {
        fs << u[i] << " ";
    }
    fs << "\nv";

    for (int i = 0; i < size; i++) {
        fs << v[i] << " ";
    }
    fs << "\np";

    for (int i = 0; i < size; i++) {
        fs << p[i] << " ";
    }
    fs << "\n";

    fs.close();
}

int main(void){
    double u[size] = {0.0};
    double un[size] = {0.0};
    double v[size] = {0.0};
    double vn[size] = {0.0};
    double p[size] = {0.0};
    double pn[size] = {0.0};
    double b[size] = {0.0};

    for(int i = 0; i < nx; i++){
        u[(nx-1)*nx + i] = 1.0;
        un[(nx-1)*nx + i] = 1.0;
    }

  //in(u,un,v,vn,p,pn,b,nx,ny);
    for(int i=1;i<nt+1;i++){
        cavity_flow(u,un,v,vn,p,pn,b); 
    }

    writeFile(u,v,p);

//Write in CSV
  //ofstream outFile;
  //outFile.open("data.csv",ios::out);
  //for(int i=0;i<nx*ny;i++){
  //  outFile<<*(u+i)<<','<<*(v+i)<<','<<*(p+i)<<endl;
  //  }
  //outFile.close();
  //cout<<"Finish"<<endl;

  //free(u);
  //free(un);
  //free(v);
  //free(vn);
  //free(p);
  //free(pn);
  //free(b);
return 0;
}