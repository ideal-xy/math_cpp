#include <iostream>
#include <vector>
#include <Accelerate/Accelerate.h>
#include <omp.h>
#include <stdexcept>

// 定义被积函数（向量化版本）
void f_vectorized(const double* x, double* result, vDSP_Length n) {
    // 示例：f(x) = x^2
    // 可以替换为其他函数，或使用 vDSP 加速复杂函数
    vDSP_vsqD(x, 1, result, 1, n); // 向量化的平方运算
}

// Riemann 积分（OpenMP + Accelerate）
double riemann_integral_hybrid(double a, double b, int n) {
    if (n <= 0) throw std::invalid_argument("n must be positive");
    
    double dx = (b - a) / n;
    double sum = 0.0;
    int num_threads = omp_get_max_threads();
    int chunk_size = n / num_threads + 1; // 确保每个线程有足够的工作

    #pragma omp parallel reduction(+:sum)
    {
        std::vector<double> x(chunk_size), fx(chunk_size);
        int thread_id = omp_get_thread_num();
        int start = thread_id * chunk_size;
        int end = std::min(start + chunk_size, n); // 防止越界
        int local_n = end - start;

        if (local_n > 0) {
            // 生成中点 x_i
            for (int i = 0; i < local_n; ++i) {
                x[i] = a + (start + i + 0.5) * dx;
            }

            // 计算 f(x)
            f_vectorized(x.data(), fx.data(), local_n);

            // 使用 Accelerate 的 vDSP 求和
            double local_sum = 0.0;
            vDSP_sveD(fx.data(), 1, &local_sum, local_n);
            sum += local_sum * dx;
        }
    }

    return sum;
}

int main() {
    double a = 0.0; // 积分下限
    double b = 1.0; // 积分上限
    int n = 10000000; // 划分数量

    try {
        double start_time = omp_get_wtime();
        double result = riemann_integral_hybrid(a, b, n);
        double end_time = omp_get_wtime();

        std::cout << "Riemann 积分结果（f(x) = x^2, 区间[" << a << ", " << b << "]）："
                  << result << std::endl;
        std::cout << "计算耗时：" << (end_time - start_time) * 1000 << " 毫秒" << std::endl;
        std::cout << "精确值：1/3 ≈ 0.333333" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "错误：" << e.what() << std::endl;
        return 1;
    }

    return 0;
}