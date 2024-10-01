#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <limits>
#include <ctime>
#include <fstream>
#include <iomanip>

using namespace std;

class DataGenerator
{
public:
    DataGenerator(double mean, double stddev)
        : gen(static_cast<unsigned>(std::time(0))), dist(mean, stddev) {}

    double generate()
    {
        return dist(gen);
    }

    std::vector<double> generateData(int size)
    {
        std::vector<double> data(size);
        for (int i = 0; i < size; ++i)
        {
            data[i] = dist(gen);
        }
        return data;
    }

private:
    std::mt19937 gen;
    std::normal_distribution<> dist;
};

class SmoothingSpline
{
public:
    SmoothingSpline(const std::vector<double> &y, double p, double mean)
        : y_(y), p_(p), mean_(mean), n_(y.size()), g_(n_, 0.0) {}
    void fit()
    {
        std::vector<std::vector<double>> A(n_, std::vector<double>(n_, 0.0));
        std::vector<double> b(n_, 0.0);
        for (int i = 0; i < n_; ++i)
        {
            for (int j = 0; j < n_; ++j)
            {
                A[i][j] = basisFunction(i, j);
            }
            b[i] = (1 - p_) * (y_[i] - mean_);
        }
        A[0][0] += p_;
        A[n_ - 1][n_ - 1] += p_;

        luSolve(A, b, g_);
    }
    double evaluate(double x) const
    {
        double result = mean_;
        int xi = static_cast<int>(x);
        if (xi >= n_ - 1)
            xi = n_ - 1;
        for (int i = 0; i < n_; ++i)
        {
            result += g_[i] * basisFunction(xi, i);
        }
        if (xi >= 0)
        {
            return fabs(result);
        }
        else
        {
            return result;
        }
    }

private:
    const std::vector<double> &y_;
    double p_;
    double mean_;
    int n_;
    std::vector<double> g_;
    double basisFunction(int i, int j) const
    {
        if (i == j)
            return 1.0;
        if (abs(i - j) == 1)
            return p_;
        return 0.0;
    }
    void luSolve(std::vector<std::vector<double>> &A, std::vector<double> &b, std::vector<double> &x)
    {
        int n = A.size();
        std::vector<int> pi(n);
        for (int i = 0; i < n; ++i)
            pi[i] = i;
        for (int k = 0; k < n; ++k)
        {
            double maxVal = 0.0;
            int maxRow = k;
            for (int i = k; i < n; ++i)
            {
                if (fabs(A[i][k]) > maxVal)
                {
                    maxVal = fabs(A[i][k]);
                    maxRow = i;
                }
            }
            std::swap(A[k], A[maxRow]);
            std::swap(b[k], b[maxRow]);
            std::swap(pi[k], pi[maxRow]);
            for (int i = k + 1; i < n; ++i)
            {
                double factor = A[i][k] / A[k][k];
                for (int j = k + 1; j < n; ++j)
                {
                    A[i][j] -= factor * A[k][j];
                }
                A[i][k] = factor;
                b[i] -= factor * b[k];
            }
        }
        for (int i = n - 1; i >= 0; --i)
        {
            double sum = 0.0;
            for (int j = i + 1; j < n; ++j)
            {
                sum += A[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / A[i][i];
        }
    }
};

int main()
{
    int n, action;
    double mean;
    double stddev;
    cout << "Введите количество наблюдений n = ";
    cin >> n;

    cout << "Введите математическое ожидание mean = ";
    cin >> mean;

    cout << "Введите среднее квадратичное отклонение stddev = ";
    cin >> stddev;

    DataGenerator generator(mean, stddev);
    vector<double> y(n);

    std::cout << "Введите действие:" << endl;
    std::cout << "1 - Ввод с файла" << endl;
    std::cout << "2 - Генерация и запись в файл" << endl;
    cin >> action;
    std::ifstream inputFile("inputed_values.txt");
    std::ofstream out1File("outputed_values.txt");
    out1File << std::fixed << std::setprecision(5);
    switch (action)
    {
    case 1:
        for (int i = 0; i < n; ++i)
        {
            inputFile >> y[i];
        }
        break;
    case 2:
        for (int i = 0; i < n; ++i)
        {
            y[i] = generator.generate();
        }
        for (int i = 0; i < n; ++i)
        {
            out1File << y[i] << endl;
        }
        break;
    }
    double p;
    cout << "Введите параметр сглаживания p = ";
    cin >> p;
    SmoothingSpline spline(y, p, mean);
    spline.fit();
    std::ofstream outFile("spline_results.txt");
    outFile << std::fixed << std::setprecision(5);
    for (int yi = 0; yi < n; yi++)
    {
        outFile << spline.evaluate(y[yi]) << endl;
    }
    cout << "Все нормально вывелось в файл. Работа завершена.";
    return 0;
}
