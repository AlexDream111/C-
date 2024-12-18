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
    vector<vector<double>> u_prev, u_curr, u_next; // 前一时刻、当前时刻、下一时刻的解
    vector<vector<double>> c_field; // 波速场

    // 构造函数
    WaveEquationSolver(int nx, int ny, int nt, double dx, double dy, double dt)
        : nx(nx), ny(ny), nt(nt), dx(dx), dy(dy), dt(dt) {
        u_prev.resize(nx, vector<double>(ny, 0.0));
        u_curr.resize(nx, vector<double>(ny, 0.0));
        u_next.resize(nx, vector<double>(ny, 0.0));
        c_field.resize(nx, vector<double>(ny, 1.0)); // 默认波速为 1.0
    }

    // 初始化波速场
    void initialize_c_field(double (*c_function)(double, double)) {
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                double x = i * dx;
                double y = j * dy;
                c_field[i][j] = c_function(x, y);
            }
        }
    }

    // 初始化初始条件
    void initialize(double (*initial_condition)(double, double)) {
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                double x = i * dx;
                double y = j * dy;
                u_curr[i][j] = initial_condition(x, y);
            }
        }
        u_prev = u_curr;
    }

    // 边界条件
    void apply_boundary_conditions() {
        for (int i = 0; i < nx; ++i) {
            u_next[i][0] = u_next[i][ny - 1] = 0.0; // 上下边界
        }
        for (int j = 0; j < ny; ++j) {
            u_next[0][j] = u_next[nx - 1][j] = 0.0; // 左右边界
        }
    }

    // 求解波动方程
    void solve() {
        for (int n = 0; n < nt; ++n) {
            for (int i = 1; i < nx - 1; ++i) {
                for (int j = 1; j < ny - 1; ++j) {
                    double c2 = c_field[i][j] * c_field[i][j];
                    double r1 = c2 * dt * dt / (dx * dx);
                    double r2 = c2 * dt * dt / (dy * dy);

                    u_next[i][j] = 2 * u_curr[i][j] - u_prev[i][j] +
                                   r1 * (u_curr[i + 1][j] - 2 * u_curr[i][j] + u_curr[i - 1][j]) +
                                   r2 * (u_curr[i][j + 1] - 2 * u_curr[i][j] + u_curr[i][j - 1]);
                }
            }
            apply_boundary_conditions();
            u_prev = u_curr;
            u_curr = u_next;
        }
    }

    // 保存结果到文件
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

    // 保存波速场
    void save_c_field(const string &filename) {
        ofstream file(filename);
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                file << i * dx << " " << j * dy << " " << c_field[i][j] << "\n";
            }
            file << "\n";
        }
        file.close();
    }
};

// 初始条件函数
double initial_condition(double x, double y) {
    double sigma = 0.1;
    return exp(-((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)) / (2 * sigma * sigma));
}

// 波速分布函数
double variable_wave_speed(double x, double y) {
    return 1.0 + 0.5 * sin(2 * M_PI * x) * cos(2 * M_PI * y); // 波速在 [0.5, 1.5] 范围内变化
}

// 主函数
int main() {
    int nx = 201, ny = 201;     // 空间网格点数
    double dx = 0.005, dy = 0.005; // 空间步长
    double dt = 0.001;           // 时间步长
    int nt = 500;                // 时间步数

    WaveEquationSolver solver(nx, ny, nt, dx, dy, dt);

    // 初始化初始条件
    solver.initialize(initial_condition);

    // 初始化波速场
    solver.initialize_c_field(variable_wave_speed);

    // 保存波速场
    solver.save_c_field("wave_speed_field.dat");

    // 求解波动方程
    solver.solve();

    // 保存波动方程结果
    solver.save_result("wave_solution_variable_c.dat");

    cout << "Wave speed field saved to 'wave_speed_field.dat'." << endl;
    cout << "Wave equation solution saved to 'wave_solution_variable_c.dat'." << endl;

    return 0;
}
