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
        : gen(static_cast<unsigned>(time(0))), dist(mean, stddev) {}

    double generate()
    {
        return dist(gen);
    }

    vector<double> generateData(int size)
    {
        vector<double> data(size);
        for (int i = 0; i < size; ++i)
        {
            data[i] = dist(gen);
        }
        return data;
    }

private:
    mt19937 gen;
    normal_distribution<> dist;
};

class SmoothingSpline {
public:
    SmoothingSpline(const vector<double>& y, double p, double mean)
        : y_(y), p_(p), mean_(mean), n_(y.size()), g_(n_, 0.0) {
        fit();
    }

    double evaluate(int i) const {
        if (i < 0) return g_.front();
        if (i >= n_) return g_.back();
        return g_[i];
    }

    void writeToFile(const string& filename) {
        ofstream outFile(filename);
        outFile << fixed << setprecision(5);
        for (int i = 0; i < n_; ++i) {
            outFile << g_[i] << "\n";
        }
        outFile.close();
    }

private:
    const vector<double>& y_;  
    double p_;                     
    double mean_;                
    int n_;                        
    vector<double> g_;         

    void fit() {
        if (p_ == 0) {
            g_ = y_;
            return;
        }
        vector<double> a(n_ - 1, -p_ / 2); 
        vector<double> b(n_, 1 + p_);   
        vector<double> c(n_ - 1, -p_ / 2); 
        vector<double> d(n_, 0.0);  
        for (int i = 0; i < n_; ++i) {
            d[i] = (1 - p_) * y_[i] + p_ * mean_;
        }
        thomasAlgorithm(a, b, c, d, g_);
    }
    void thomasAlgorithm(const vector<double>& a, vector<double>& b, 
                         const vector<double>& c, const vector<double>& d, 
                         vector<double>& x) {
        int n = b.size();
        vector<double> c_prime(n, 0.0);
        vector<double> d_prime(n, 0.0);
        c_prime[0] = c[0] / b[0];
        d_prime[0] = d[0] / b[0];
        for (int i = 1; i < n - 1; ++i) {
            double m = 1.0 / (b[i] - a[i - 1] * c_prime[i - 1]);
            c_prime[i] = c[i] * m;
            d_prime[i] = (d[i] - a[i - 1] * d_prime[i - 1]) * m;
        }

        d_prime[n - 1] = (d[n - 1] - a[n - 2] * d_prime[n - 2]) / 
                         (b[n - 1] - a[n - 2] * c_prime[n - 2]);
        x[n - 1] = d_prime[n - 1];
        for (int i = n - 2; i >= 0; --i) {
            x[i] = d_prime[i] - c_prime[i] * x[i + 1];
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

    cout << "Введите действие:" << endl;
    cout << "1 - Ввод с файла" << endl;
    cout << "2 - Генерация и запись в файл" << endl;
    cin >> action;
    ifstream inputFile("inputed_values.txt");
    ofstream out1File("outputed_values.txt");
    out1File << fixed << setprecision(5);
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
    
    spline.writeToFile("spline_results.txt");

    cout << "Все нормально вывелось в файл. Работа завершена.";
    return 0;
}
