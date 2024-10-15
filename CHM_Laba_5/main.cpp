#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <chrono>
#include <iomanip>
#include "integrity.h"

using namespace std;

int main()
{

    auto f = [](double x)
    {
        return - 2 / pow(x, 2);
    };

    double b;
    cout << "Введите верхний предел интегрирования (b): ";
    cin >> b;

    double a;
    cout << "Введите нижний предел интегрирования (a): ";
    cin >> a;

    int n;
    cout << "Введите начальное количество интервалов: ";
    cin >> n;

    Integrator integrator;
    int action;
    cout << "Введите необходимый способ интегрирования:" << endl
         << "1 - Метод прямоугольников" << endl
         << "2 - Метод трапеций" << endl
         << "3 - Метод Симпсона(порабол)" << endl;
    cin >> action;
    IntegrationMethod standartMethod = static_cast<IntegrationMethod>(action);
    integrator.performIntegration(standartMethod, f, a, b, n);
    cout << "Введите метод квадратуры Гаусса:" << endl
         << "1 - Гаусс-2" << endl
         << "2 - Гаусс-3" << endl
         << "3 - Гаусс-4" << endl;
    cin >> action;
    IntegrateGauss gaussMethod = static_cast<IntegrateGauss>(action);
    cout << "Метод квадратуры Гаусса :" << endl;
    integrator.performGaussIntegration(gaussMethod, f, a, b, n);

    

    return 0;
}
