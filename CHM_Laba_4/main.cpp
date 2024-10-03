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

    double epsilon;
    cout << "Введите точность вычисления интеграла (epsilon): ";
    cin >> epsilon;

    auto startRect = chrono::high_resolution_clock::now();
    pair<double, int> resultRect = adaptiveIntegration(integrateRectangles, f, a, b, epsilon);
    auto endRect = chrono::high_resolution_clock::now();
    chrono::duration<double, micro> timeRect = endRect - startRect;
    cout << "Метод прямоугольников: результат = " << resultRect.first << ", n = " << resultRect.second << ", время = " << timeRect.count() << " мс" << endl;

    auto startTri = chrono::high_resolution_clock::now();
    pair<double, int> resultTri = adaptiveIntegration(integrateTriangles, f, a, b, epsilon);
    auto endTri = chrono::high_resolution_clock::now();
    chrono::duration<double, micro> timeTri = endTri - startTri;
    cout << "Метод трапеций: результат = " << resultTri.first << ", n = " << resultTri.second << ", время = " << timeTri.count() << " мс" << endl;

    auto startParab = chrono::high_resolution_clock::now();
    pair<double, int> resultParab = adaptiveIntegration(integrateParabolas, f, a, b, epsilon);
    auto endParab = chrono::high_resolution_clock::now();
    chrono::duration<double, micro> timeParab = endParab - startParab;
    cout << "Метод парабол (Симпсона): результат = " << resultParab.first << ", n = " << resultParab.second << ", время = " << timeParab.count() << " мс" << endl;

    auto startGauss = chrono::high_resolution_clock::now();
    pair<double, int> resultGauss = adaptiveIntegration(integrateGauss, f, a, b, epsilon);
    auto endGauss = chrono::high_resolution_clock::now();
    chrono::duration<double, micro> timeGauss = endGauss - startGauss;
    cout << "Квадратура Гаусса: результат = " << resultGauss.first << ", количество сегментов n = " << resultGauss.second << ", время = " << timeGauss.count() << " мс" << endl;

    auto func = [](double x)
    {
        return x - 1;
    };

    int n;
    cout << "Введите количество узлов интегрирования N: ";
    cin >> n;

    double result = integrate_with_hermite(n, func);

    cout << "Результат = " << std::setprecision(20) << result << endl;
    return 0;
}
