#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

// 定义二维波动方程的求解器
class WaveEquationSolver {
public:
    int nx, ny;          // 网格点数
    int nt;              // 时间步数
    double dx, dy, dt;   // 空间步长和时间步长
    double c;            // 波速
    vector<vector<double>> u_prev, u_curr, u_next; // 二维向量,前一时刻、当前时刻、下一时刻的解

    WaveEquationSolver(int nx, int ny, int nt, double dx, double dy, double dt, double c)
        : nx(nx), ny(ny), nt(nt), dx(dx), dy(dy), dt(dt), c(c) {
        u_prev.resize(nx, vector<double>(ny, 0.0));
        u_curr.resize(nx, vector<double>(ny, 0.0));
        u_next.resize(nx, vector<double>(ny, 0.0));// 对类的成员变量进行初始化,为三个解的存储矩阵分配空间，并将其初始值设为 0
    }
   
    void initialize(double (*initial_condition)(double, double)) {
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                double x = i * dx;
                double y = j * dy;
                u_curr[i][j] = initial_condition(x, y);// 根据给定的初始条件函数 initial_condition，初始化当前时刻（u_curr）的解
            }
        }
        u_prev = u_curr;// 前一时刻（u_prev）与当前时刻（u_curr）的解相同
    }

    void apply_boundary_conditions() {
        for (int i = 0; i < nx; ++i) {
            u_next[i][0] = u_next[i][ny - 1] = 0.0; // 上下边界，u_next 表示新解，即下一时间步的波动分布，边界条件必须施加在 u_next 上
        }
        for (int j = 0; j < ny; ++j) {
            u_next[0][j] = u_next[nx - 1][j] = 0.0; // 左右边界
        }
    }

    void solve() {
        double c2 = c * c;
        double r1 = c2 * dt * dt / (dx * dx);
        double r2 = c2 * dt * dt / (dy * dy);

        for (int n = 0; n < nt; ++n) {
            for (int i = 1; i < nx - 1; ++i) {
                for (int j = 1; j < ny - 1; ++j) {
                    u_next[i][j] = 2 * u_curr[i][j] - u_prev[i][j] +
                                   r1 * (u_curr[i + 1][j] - 2 * u_curr[i][j] + u_curr[i - 1][j]) +
                                   r2 * (u_curr[i][j + 1] - 2 * u_curr[i][j] + u_curr[i][j - 1]);
                }
            }
            apply_boundary_conditions();// 边界点的值设置为 0
            u_prev = u_curr;
            u_curr = u_next;
        }
    }

    void save_result(const string &filename) {
        ofstream file(filename);
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                file << i * dx << " " << j * dy << " " << u_curr[i][j] << "\n";
            }
            file << "\n";
        }
        file.close();
    }
};

// 初始条件函数
double initial_condition(double x, double y) {
    double sigma = 0.1;
    return exp( -1 * ((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)) / (2 * sigma * sigma)) * (1 - 0.5 * (x * x + y * y)/(sigma * sigma));
}

// 主函数
int main() {
    int nx = 201, ny = 201;     // 空间网格点数
    double dx = 0.005, dy = 0.005; // 空间步长
    double dt = 0.001;           // 时间步长
    double c = 1.0;              // 波速
    int nt = 500;                // 时间步数


    WaveEquationSolver solver(nx, ny, nt, dx, dy, dt, c);
    
    cout << "The time point is: " << nt << endl;
    solver.initialize(initial_condition);
    solver.solve();
    solver.save_result("wave_solution.dat");

    cout << "Wave equation solution saved to 'wave_solution.dat'." << endl;
    return 0;
}
