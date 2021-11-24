#include <iostream>
#include <cmath>
#include <vector>
#include <ctime>
#include <limits>
#include <iomanip>
#include <fstream>

/*
 *  Вариант 1. Метод фильтрации: среднее арифметическое, метрика близости: Евклидова
 */

class stochastic_filtering{
public:
    int K = 100;
    double x_min = 0;
    double x_max = M_PI;
    double noise_amplitude = 0.5;
    double a = noise_amplitude / 2;
    int r = 3;
    int L = 10;
    double P = 0.95;
    double e = 0.01;
    int M = static_cast<int>((r - 1) / 2);
    int N = static_cast<int>(round(log(1 - P) / log(1 - e / (x_max - x_min))));
//    double h, w, d, J;
    std::vector<double> alpha;
    std::vector<double> noise;
    double noisy_signal(const int &k_iter){
        return signal(k_iter) + noise[k_iter];
    }
    double signal(const int &x){
        return sin(x_k(x)) + 0.5;
    }
    double x_k(const int &k_iter){
        return x_min + k_iter * (x_max - x_min) / K;
    }
    double filtered_signal(const int &k_iter){
        double result = 0;
        for (int i = k_iter - M; i <= k_iter + M; ++i) {
//            if (k_iter >= M && k_iter <= K - M)
            if (i >= 0 && k_iter <= K)
                result += noisy_signal(i) * alpha[i + M - k_iter];
        }
//        for (int i = k_iter - M; i <= k_iter + M; ++i) {
//            int j = (i + K + 1) % (K + 1);
//            result += noisy_signal(j) * alpha[i + M - k_iter];
//        }
        return result;
    }
    double rnd(const double &first, const double &second){
        if (first == 0 && second == 0){
            return 0;
        }
        int precision = 1000000;
        return static_cast<double>((std::rand() % static_cast<int>((second - first) * precision)) + static_cast<int>(first * precision)) / precision;
    }
    void init_alpha(){ //std::vector<double> &recv_alpha
        for (int i = 0; i < r; ++i) {
            alpha[i] = 0;
        }
        alpha[M] = rnd(0, 1);
        double sum = 0;
        for (int m = 2; m <= M; ++m) {
            for (int s = m; s < r - m; ++s) {
                sum += alpha[s];
            }
            alpha[m - 1] = alpha[r - m] = 0.5 * rnd(0, 1 - sum);
        }
        sum = 0;
        for (int s = 1; s < r - 1; ++s) {
            sum += alpha[s];
        }
        alpha[0] = alpha[r - 1] = 0.5 * (1 - sum);
        sum = 0; // нормализация
        for (int i = 0; i < r; ++i) {
            sum += alpha[i];
        }
        for (int i = 0; i < r; ++i) {
            alpha[i] /= sum;
        }
    }
    double generate_w(){
        double result = 0;
        for (int k = 1; k <= K; ++k) {
//            result += std::pow(filtered_signal(k) - filtered_signal(k - 1), 2);
            result += (filtered_signal(k) - filtered_signal(k - 1)) * (filtered_signal(k) - filtered_signal(k - 1));
        }
        result = sqrt(result);
        return result;
    }
    double generate_d(){
        double result = 0;
        for (int k = 0; k <= K; ++k) {
            result += std::pow(filtered_signal(k) - noisy_signal(k), 2);
        }
        result /= K;
        result = sqrt(result);
        return result;
    }
    void print_line(){
        std::cout << "+-----+--------+";
        for (int i = 0; i < r * 8; ++i) {
            std::cout << '-';
        }
        std::cout << "+--------+--------+\n";
    }
    explicit stochastic_filtering(const std::string& name, const int &r = 3): r(r){
        alpha.resize(r);
        noise.resize(K + 1);
        for (auto &x : noise) {
            x = rnd(-a, a);
        }
        std::vector<double> h(L + 1);
        std::vector<double> w(L + 1);
        std::vector<double> d(L + 1);
        std::vector<double> dist(L + 1);
        std::vector<double> J(L + 1, std::numeric_limits<double>::max());
        std::vector<std::vector<double>> best_alpha(L + 1);
        for (int i = 0; i <= L; ++i) {
            double temp_h = static_cast<double>(i) / L;
            double temp_J = 0, temp_w = 0, temp_d = 0;
//            std::vector<double> temp_alpha(r);
            for (int j = 0; j < N; ++j) {
                init_alpha();
                temp_w = generate_w();
                temp_d = generate_d();
                temp_J = temp_h * temp_w + (1 - temp_h) * temp_d;
                if (temp_J < J[i]){
                    J[i] = temp_J;
                    w[i] = temp_w;
                    d[i] = temp_d;
                    h[i] = temp_h;
                    best_alpha[i] = alpha;
                    dist[i] = sqrt(std::pow(w[i], 2) + std::pow(d[i], 2));
                }
            }
        }
        print_line(); //вывод таблицы
        std::cout << "|  h  |  dist  |";
        for (int i = 0; i < (8 * r - 5) / 2; ++i) {
            std::cout << ' ';
        }
        std::cout << "alpha";
        for (int i = 0; i < (8 * r - 5) / 2; ++i) {
            std::cout << ' ';
        }
        std::cout << " |    w   |    d   |\n";
        print_line();
        for (int i = 0; i <= L; ++i){
            std::cout << "| " << std::fixed << std::setprecision(1) << h[i] << " | "
            << std::setprecision(4) << dist[i] << " | "
            << std::setprecision(4) << best_alpha[i][0];
            for (int j = 1; j < r; ++j) {
                std::cout << ", " << std::fixed << std::setprecision(4) << best_alpha[i][j];
            }
            std::cout << " | " << std::fixed << std::setprecision(4) << w[i]
            << " | " << std::fixed << std::setprecision(4) << d[i] << " |\n";
        }
        print_line();
        size_t best = std::distance(dist.begin(), std::min_element(dist.begin(), dist.end()));
        std::cout << "\n+-----+--------+--------+--------+\n";
        std::cout << "|  h* |    J   |    w   |    d   |\n";
        std::cout << "+-----+--------+--------+--------+\n";
        std::cout << "| " << std::fixed << std::setprecision(1) << h[best] << " | "
        << std::setprecision(4) << J[best] << " | "
        << std::setprecision(4) << w[best] << " | "
        << std::setprecision(4) << d[best] << " |\n";
        std::cout << "+-----+--------+--------+--------+\n\n\n";
        alpha = best_alpha[best];
        std::ofstream out;
        out.open("signals_" + name +".txt");
        for (int i = 0; i <= K; ++i) {
            out << x_k(i) << ' ' << signal(i) << ' '
                << noisy_signal(i) << ' '
                << filtered_signal(i) << '\n';
        }
        out.close();
        out.open("w_d_" + name + ".txt");
        for (int i = 0; i < L + 1; ++i) {
            out << w[i] << ' ' << d[i] << '\n';
        }
        out.close();
    }

};

int main() {
    std::srand(std::time(nullptr));
    stochastic_filtering a("1");
    stochastic_filtering b("2", 5);
//    for (int i = 0; i < 20; ++i) {
//        std::cout << a.rnd(-1.1, 1.1) << '\n';
//    }
    return 0;
}
