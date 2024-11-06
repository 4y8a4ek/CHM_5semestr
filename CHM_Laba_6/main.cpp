#include <iostream>
#include <cmath>
#include <chrono>
#include "./Fourier_Transform.h"
double A, B, omega, phi;
int count = 512;
void input_param(){
    std::cout << "Введите значение A =";
    std::cin >> A;
    std::cout << "Введите значение B =";
    std::cin >> B;
    std::cout << "Введите значение omega =";
    std::cin >> omega;
    std::cout << "Введите значение phi = π/";
    std::cin >> phi;
    phi = M_PI/phi;
    return;
}
int main()
{
    input_param();
    auto f1 = [](double x)
    {   
        return std::complex<double>(A + B * cos(2*M_PI* omega*x/count + phi));
    };

    auto f2 = [](double x)
    {   
        return std::complex<double>(cos(2*M_PI*x/count) + 0.01*cos(2*M_PI*omega*x/count));
    };
    Fourier_Transform transformer(f1, A, B, omega, phi, count);
    Fourier_Transform trans_2(f2, A, B, omega, phi, count);
    std::cout << "DFT" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    transformer.perform_dft();
    transformer.print();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Time taken by DFT: " << duration.count() << " seconds" << std::endl;

    std::cout << std::endl << "FFT" << std::endl;
    auto start_fft = std::chrono::high_resolution_clock::now();
    transformer.perform_fft();
    transformer.print();
    auto end_fft = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_fft = end_fft - start_fft;
    std::cout << "Time taken by FFT: " << duration_fft.count() << " seconds" << std::endl;

    trans_2.perform_fft();
    trans_2.print();
    std::cout << std::endl;
    trans_2.clear_noise();
    trans_2.print();
    trans_2.perform_idft();
    trans_2.save_to_files();
}
