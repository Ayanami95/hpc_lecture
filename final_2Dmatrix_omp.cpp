#include <iostream>
#include <fstream>
#include<vector>
using namespace std;
typedef vector<vector<float>> Matrix;


Matrix build_up_b(Matrix b, float rho, float dt, Matrix u, Matrix v, float dx, float dy) {
    int row = b.size();
    int col = b[0].size();

    #pragma omp parallel for
    for (int m = 0; m < row - 2; m++) {
        for (int n = 0; n < col - 2; n++) {
            float dudx = (u[m + 1][n + 2] - u[m + 1][n]) / (2.0 * dx);
            float dvdy = (v[m + 2][n + 1] - v[m][n + 1]) / (2.0 * dy);
            float dudy = (u[m + 2][n + 1] - u[m][n + 1]) / (2.0 * dy);
            float dvdx = (v[m + 1][n + 2] - v[m + 1][n]) / (2.0 * dx);

            b[m + 1][n + 1] = rho * 1 / dt * (dudx + dvdy) - dudx * dudx - 2 * dudy * dvdx - dvdy * dvdy;
        }
    }
    return b;
}
Matrix pressure_poisson(Matrix p, float dx, float dy, Matrix b, int nit) {
    int row = p.size();
    int col = p[0].size();
    Matrix pn(row, vector<float>(col));
    for (int k = 0; k < nit; k++) {
        pn = p;

        #pragma omp parallel for
        for (int m = 0; m < row - 2; m++) {
            for (int n = 0; n < col - 2; n++) {
                //$$p_{i,j}^{n}=\frac{(p_{i+1,j}^{n}+p_{i-1,j}^{n})\Delta y^2+(p_{i,j+1}^{n}+p_{i,j-1}^{n})\Delta x^2-b_{i,j}^{n}\Delta x^2\Delta y^2}{2(\Delta x^2+\Delta y^2)}$$
                p[m + 1][n + 1] = ((pn[m + 1][n + 2] + pn[m + 1][n]) * dy * dy
                    + (pn[m + 2][n + 1] + pn[m][n + 1]) * dx * dx
                    - b[m + 1][n + 1] * dx * dx * dy * dy)
                    / (2 * (dx * dx + dy * dy));
            }
        }

        //dp/dx = 0 at x = 2 and dp/dx = 0 at x = 0
        #pragma omp parallel for
        for (int m = 0; m < row; m++) {
            p[m][col - 1] = p[m][col - 2];
            p[m][0] = p[m][1];
        }

        //dp/dy = 0 at y = 0 and p = 0 at y = 2
        #pragma omp parallel for
        for (int n = 0; n < col; n++) {
            p[0][n] = p[1][n];
            p[row - 1][n] = 0;
        }

    }
    return p;
}
void cavity_flow(int nt, Matrix& u, Matrix& v, float dt, float dx, float dy, Matrix& p, float rho, float nu, int nit) {
    int row = u.size();
    int col = u[0].size();
    Matrix b(row, vector<float>(col, 0));
    Matrix un(row, vector<float>(col));
    Matrix vn(row, vector<float>(col));

    for (int i = 0; i < nt; i++) {
        un = u;
        vn = v;
        b = build_up_b(b, rho, dt, u, v, dx, dy);
        p = pressure_poisson(p, dx, dy, b, nit);

        #pragma omp parallel for
        for (int m = 0; m < row - 2; m++) {
            for (int n = 0; n < col - 2; n++) {
                u[m + 1][n + 1] = (un[m + 1][n + 1]
                    - un[m + 1][n + 1] * dt / dx * (un[m + 1][n + 1] - un[m + 1][n])
                    - vn[m + 1][n + 1] * dt / dy * (un[m + 1][n + 1] - un[m][n + 1])
                    - dt / (2 * rho * dx) * (p[m + 1][n + 2] - p[m + 1][n])
                    + nu * (dt / (dx * dx) * (un[m + 1][n + 2] - 2 * un[m + 1][n + 1] + un[m + 1][n])
                        + dt / (dy * dy) * (un[m + 2][n + 1] - 2 * un[m + 1][n + 1] + un[m][n + 1])));

                v[m + 1][n + 1] = (vn[m + 1][n + 1]
                    - un[m + 1][n + 1] * dt / dx * (vn[m + 1][n + 1] - vn[m + 1][n])
                    - vn[m + 1][n + 1] * dt / dy * (vn[m + 1][n + 1] - vn[m][n + 1])
                    - dt / (2 * rho * dx) * (p[m + 2][n + 1] - p[m + 0][n+1])
                    + nu * (dt / (dx * dx) * (vn[m + 1][n + 2] - 2 * vn[m + 1][n + 1] + vn[m + 1][n])
                        + dt / (dy * dy) * (vn[m + 2][n + 1] - 2 * vn[m + 1][n + 1] + vn[m][n + 1])));
            }
        }

        //  u[0,:] =  v[0,:] = v[-1,:] = 0  u[-1,:] =1
        #pragma omp parallel for
        for (int n = 0; n < col; n++) {
            u[0][n] = 0;
            u[row - 1][n] = 1;
            v[0][n] = 0;
            v[row - 1][n] = 0;
        }

        //  u[:,0] = u[:,-1] = v[:,0] = v[:,-1] = 0
        #pragma omp parallel for
        for (int m = 0; m < row; m++) {
            u[m][0] = 0;
            u[m][col - 1] = 0;
            v[m][0] = 0;
            v[m][col - 1] = 0;
        }

    }
}

void writeFile(Matrix& u, Matrix& v, Matrix& p) {
    ofstream fs("cavity_cpp_results_final.txt");
    int row = u.size();
    int col = u[0].size();
    // u
    fs << "u ";
    for (int m = 0; m < row; m++) {
        for (int n = 0; n < col; n++) {
            fs << u[m][n] << " ";
        }
        fs << "\n";
    }
    fs << "\nv";

    for (int m = 0; m < row; m++) {
        for (int n = 0; n < col; n++) {
            fs << v[m][n] << " ";
        }
        fs << "\n";
    }
    fs << "\np";

    for (int m = 0; m < row; m++) {
        for (int n = 0; n < col; n++) {
            fs << p[m][n] << " ";
        }
        fs << "\n";
    }
    fs << "\n";

    fs.close();
}

int main()
{
    int nx = 41;
    int ny = 41;
    int nt = 100;
    int nit = 50;

    float c = 1;
    float dx = 2.0 / (nx - 1.0);
    float dy = 2.0 / (ny - 1.0);

    float rho = 1;
    float nu = 0.1;
    float dt = 0.001;

    Matrix u(ny, vector<float>(nx, 0));
    Matrix v(ny, vector<float>(nx, 0));
    Matrix p(ny, vector<float>(nx, 0));

    cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu, nit);
    writeFile(u, v, p);
}

