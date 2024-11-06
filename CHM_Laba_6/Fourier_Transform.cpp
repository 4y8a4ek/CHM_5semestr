#include "./Fourier_Transform.h"
#define eps 1e-7
double Fourier_Transform::sign(const double &x)
{
    if (x < 0)
    {
        return -1.0;
    }
    else
    {
        return 1.0;
    }
}
double Fourier_Transform::phase(std::complex<double> &z)
{
    double result;
    if (z.real() > 0.0)
    {
        return std::atan(z.imag() / z.real());
    }
    if (z.real() == 0.0)
    {
        return sign(z.imag()) * M_PI / 2.0;
    }
    if (z.imag() >= 0.0 && z.real() < 0)
    {
        return M_PI + std::atan(z.imag() / z.real());
    }
    if (z.imag() < 0.0 && z.real() < 0)
    {
        return std::atan(z.imag() / z.real()) - M_PI;
    }
}

void Fourier_Transform::calculate_dft()
{
    std::complex<double> current = {0.0, 0.0};
    for (int m = 0; m < N; m++)
    {
        for (int n = 0; n < N; n++)
        {
            double angle = -2 * M_PI * m * n / N;
            std::complex<double> exp_part(std::cos(angle), std::sin(angle));
            current += z[n] * exp_part;
        }
        z_calculated[m] = current;
        current = {0, 0};
    }
}

void Fourier_Transform::calculate_fft()
{
    std::complex<double> u;
    std::complex<double> v;
    for (int i = 0; i < N / 2.0; i++)
    {
        u = {0.0, 0.0};
        v = {0.0, 0.0};
        for (int j = 0; j < N / 2.0; j++)
        {
            double angle = -2 * M_PI * j * i / (N / 2.0);
            std::complex<double> exp_part(std::cos(angle), std::sin(angle));
            u += z[j * 2] * exp_part;
            v += z[j * 2 + 1] * exp_part;
        }
        double angle = -2 * M_PI * i / float(N);
        std::complex<double> exp_part(std::cos(angle), std::sin(angle));
        z_calculated[i] = u + exp_part * v;
        z_calculated[N / 2.0 + i] = u - exp_part * v;
    }
}

void Fourier_Transform::calculate_idft(){
    for(int n = 0; n < N; n++){
        std::complex<double> current = {0.0,0.0};
        for(int m = 0; m < N; m++){
            double angle = 2 * M_PI * m * n / N;
            std::complex<double> exp_part(std::cos(angle), std::sin(angle));
            current += z_calculated[m]*exp_part;
        }
        z_restored[n] = current / double(N);
    }
}

void Fourier_Transform::clear_noise(){
    double max_amplitude = -1e10;
    for (int i = 0; i < N; i++)
    {
        if(sqrt(pow(z_calculated[i].real(), 2) + pow(z_calculated[i].imag(), 2)) > max_amplitude && sqrt(pow(z_calculated[i].real(), 2) + pow(z_calculated[i].imag(), 2)) > eps){
            max_amplitude = sqrt(pow(z_calculated[i].real(), 2) + pow(z_calculated[i].imag(), 2));
        }
    }
    for (int j = 0; j < N; j++)
    {
        if(sqrt(pow(z_calculated[j].real(), 2) + pow(z_calculated[j].imag(), 2)) < max_amplitude-1 && sqrt(pow(z_calculated[j].real(), 2) + pow(z_calculated[j].imag(), 2)) > eps){
            z_calculated[j] = {0.0, 0.0};
        }
    }    
}

void Fourier_Transform::perform_dft()
{
    calculate_dft();
}

void Fourier_Transform::perform_fft()
{
    calculate_fft();
}

void Fourier_Transform::perform_idft(){
    calculate_idft();
}

void Fourier_Transform::print()
{
    for (int m = 0; m < N; m++)
    {
        if (sqrt(pow(z_calculated[m].real(), 2) + pow(z_calculated[m].imag(), 2)) > eps)
        {
            std::cout << m << "  |  ";
            printf("%e", z[m].real());
            std::cout << "  |  ";
            printf("%e", z_calculated[m].real());
            std::cout << "  |  ";
            printf("%e", z_calculated[m].imag());
            std::cout << "  |  ";
            printf("%e", sqrt(pow(z_calculated[m].real(), 2) + pow(z_calculated[m].imag(), 2))); // Амплитуда
            std::cout << "  |  ";
            printf("%e", phase(z_calculated[m])); // Фаза
            std::cout << std::endl;
        }
    }
}

void Fourier_Transform::save_to_files(){
    std::ofstream file_one("z.txt");
    std::ofstream file_two("z_restore.txt");
    for ( auto& number : z) {
            file_one << number.real() << "\n";
        }
    for ( auto& value : z_restored) {
            file_two << value.real() << "\n";
        }
}

Fourier_Transform::Fourier_Transform(std::function<std::complex<double>(double)> f, const double &A, const double &B, const double &omega, const double &phi, const double &count)
{
    N = count;
    z.resize(N);
    z_calculated.resize(N);
    z_restored.resize(N);
    for (int i = 0; i < N; i++)
    {
        z[i] = f(i);
    }
}
