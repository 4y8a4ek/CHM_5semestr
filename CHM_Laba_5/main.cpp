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
        return -2.0 / pow(x, 2);
    };

    double a;
    cout << "Введите верхний предел интегрирования (a): ";
    cin >> a;

    double b;
    cout << "Введите нижний предел интегрирования (b): ";
    cin >> b;

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
    vector<pair<double, int>> resultRect;
    vector<pair<double, int>> resultTri;
    vector<pair<double, int>> resultParab;
    switch (action)
    {
    case 1:
        resultRect = integrator.adaptiveIntegration(Integrator::integrateRectangles, f, a, b, n);
        cout << "Метод прямоугольников:" << endl;
        for (int i = 0; i < 3; i++)
        {
            cout << " Результат = " << resultRect[i].first << ", n = " << resultRect[i].second << endl;
        }
        break;
    case 2:
        resultTri = integrator.adaptiveIntegration(Integrator::integrateTriangles, f, a, b, n);
        cout << "Метод трапеций:" << endl;
        for (int i = 0; i < 3; i++)
        {
            cout << " Результат = " << resultTri[i].first << ", n = " << resultTri[i].second << endl;
        }
        break;
    case 3:
        resultParab = integrator.adaptiveIntegration(Integrator::integrateParabolas, f, a, b, n);
        cout << "Метод парабол (Симпсона):" << endl;
        for (int i = 0; i < 3; i++)
        {
            cout << " Результат = " << resultParab[i].first << ", n = " << resultParab[i].second << endl;
        }
        break;
    }

    vector<vector<double>> gaussWeights = {{1.0, 1.0}, {(5.0 / 9.0), (8.0 / 9.0), (5.0 / 9.0)}, {(18 - sqrt(36)) / 36, (18 + sqrt(36)) / 36, (18 + sqrt(36)) / 36, (18 - sqrt(36)) / 36}};
    vector<vector<double>> gaussNodes = {{-1 / sqrt(3), 1 / sqrt(3)}, {-sqrt((3.0 / 5.0)), 0.0 , sqrt((3.0 / 5.0))}, {-sqrt((3 + 2 * sqrt(6 / 5)) / 7), -sqrt((3 - 2 * sqrt(6 / 5)) / 7), sqrt((3 - 2 * sqrt(6 / 5)) / 7), sqrt((3 + 2 * sqrt(6 / 5)) / 7)}};
    cout << "Введите метод квадратуры Гаусса:" << endl
         << "1 - Гаусс-2" << endl
         << "2 - Гаусс-3" << endl
         << "3 - Гаусс-4" << endl;
    cin >> action;
    double resultGauss;
    cout << "Квадратура Гаусса:" << endl;
    switch (action)
    {
    case 1:
        for (int i = 0; i < 3; i++)
        {
            resultGauss = integrator.integrateGauss(f, a, b, n, gaussWeights[0], gaussNodes[0]);
            cout << " Результат = " << resultGauss << ", количество сегментов n = " << n << endl;
            n *= 2;
        }
        break;
    case 2:
        for (int i = 0; i < 3; i++)
        {
            resultGauss = integrator.integrateGauss(f, a, b, n, gaussWeights[1], gaussNodes[1]);
            cout << " Результат = " << resultGauss << ", количество сегментов n = " << n << endl;
            n *= 2;
        }
        break;
    case 3:
        for (int i = 0; i < 3; i++)
        {
            resultGauss = integrator.integrateGauss(f, a, b, n, gaussWeights[2], gaussNodes[2]);
            cout << " Результат = " << resultGauss << ", количество сегментов n = " << n << endl;
            n *= 2;
        }
        break;
    }
    return 0;
}
