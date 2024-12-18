import numpy as np
import matplotlib.pyplot as plt

# 可视化波速场
def visualize_wave_speed(file_name, output_name):
    data = np.loadtxt(file_name)
    x = data[:, 0]
    y = data[:, 1]
    c = data[:, 2]
    nx = int(np.sqrt(len(x)))

    # 检查数据是否可以重塑为正方形网格
    if nx * nx != len(x):
        raise ValueError(f"数据点数量 {len(x)} 无法形成正方形网格！请检查 {file_name} 文件格式。")

    X, Y = x.reshape(nx, nx), y.reshape(nx, nx)
    C = c.reshape(nx, nx)

    plt.figure(figsize=(8, 6))
    plt.contourf(X, Y, C, levels=50, cmap='plasma')
    plt.colorbar(label="Wave Speed (c)")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Wave Speed Field")
    plt.savefig(output_name)
    print(f"波速场图像已保存为 {output_name}")


# 可视化波动解
def visualize_wave_solution(file_name, output_name):
    data = np.loadtxt(file_name)
    x = data[:, 0]
    y = data[:, 1]
    u = data[:, 2]
    nx = int(np.sqrt(len(x)))

    # 检查数据是否可以重塑为正方形网格
    if nx * nx != len(x):
        raise ValueError(f"数据点数量 {len(x)} 无法形成正方形网格！请检查 {file_name} 文件格式。")

    X, Y = x.reshape(nx, nx), y.reshape(nx, nx)
    U = u.reshape(nx, nx)

    plt.figure(figsize=(8, 6))
    plt.contourf(X, Y, U, levels=50, cmap='viridis')
    plt.colorbar(label="Wave Amplitude")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Wave Equation Solution")
    plt.savefig(output_name)
    print(f"波动解图像已保存为 {output_name}")


# 主程序
if __name__ == "__main__":
    # 可视化波速场
    visualize_wave_speed("wave_speed_field.dat", "wave_speed_field.png")

    # 可视化波动解
    visualize_wave_solution("wave_solution_variable_c.dat", "wave_solution_variable_c.png")
