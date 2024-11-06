#include <iostream>
#include <vector>
#include <complex>
#include <iomanip>
#include <functional>
#include <fstream>
class Fourier_Transform
{
private:
    std::vector<std::complex<double>> z, z_calculated, z_restored;
    int N;
    double sign(const double &x);
    double phase(std::complex<double> &z);
    void calculate_dft();
    void calculate_fft();
    void calculate_idft();
public:
    void clear_noise();
    void perform_dft();
    void perform_fft();
    void perform_idft();
    void print();
    void save_to_files();
Fourier_Transform(std::function<std::complex<double>(double)> f, const double &A, const double &B, const double &omega, const double &phi, const double &count);
};
