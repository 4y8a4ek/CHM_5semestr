#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include "integrity.h"

using namespace std;

vector<vector<double>> gaussWeights = {{1.0, 1.0}, {(5.0 / 9.0), (8.0 / 9.0), (5.0 / 9.0)}, {(18 - sqrt(36)) / 36, (18 + sqrt(36)) / 36, (18 + sqrt(36)) / 36, (18 - sqrt(36)) / 36}};
vector<vector<double>> gaussNodes = {{-1 / sqrt(3), 1 / sqrt(3)}, {-sqrt((3.0 / 5.0)), 0.0, sqrt((3.0 / 5.0))}, {-sqrt((3 + 2 * sqrt(6 / 5)) / 7), -sqrt((3 - 2 * sqrt(6 / 5)) / 7), sqrt((3 - 2 * sqrt(6 / 5)) / 7), sqrt((3 + 2 * sqrt(6 / 5)) / 7)}};

double Integrator::gaussianQuadrature(const function<double(double)> &f, const double &a, const double &b, const vector<double> &gaussWeights, const vector<double> &gaussNodes)
{
    double sum = 0.0;
    double x;
    for (size_t i = 0; i < gaussNodes.size(); ++i)
    {
        x = ((b - a) / 2) * gaussNodes[i] + (a + b) / 2;
        sum += gaussWeights[i] * f(x);
    }
    return sum * (b - a) / 2;
}

double Integrator::integrateGauss(const function<double(double)> &f, const double &a, const double &b, const int &segments, const vector<double> &gaussWeights, const vector<double> &gaussNodes)
{
    double segmentWidth = (b - a) / segments;
    double totalResult = 0.0;
    double segmentA, segmentB;
    for (int i = 0; i < segments; ++i)
    {
        segmentA = a + i * segmentWidth;
        segmentB = segmentA + segmentWidth;
        totalResult += gaussianQuadrature(f, segmentA, segmentB, gaussWeights, gaussNodes); // gaussianQuadrature для одного сегмента
    }

    return totalResult;
}

double Integrator::integrateRectangles(const function<double(double)> &f, const double &a, const double &b, int &n)
{
    double h = (b - a) / n;
    double sum = 0.0;
    double x;
    for (int i = 0; i < n; ++i)
    {
        x = a + i * h;
        sum += f(x) * h;
    }
    return sum;
}

double Integrator::integrateTriangles(const function<double(double)> &f, const double &a, const double &b, int &n)
{
    double h = (b - a) / n;
    double sum = 0.0;
    double x1, x2;
    for (int i = 0; i < n; ++i)
    {
        x1 = a + i * h;
        x2 = a + (i + 1) * h;
        sum += (f(x1) + f(x2)) * h / 2;
    }
    return sum;
}
double Integrator::integrateParabolas(const function<double(double)> &f, const double &a, const double &b, int &n)
{
    if (n % 2 != 0)
        n++;
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    double x;
    for (int i = 1; i < n; ++i)
    {
        x = a + i * h;
        if (i % 2 == 0)
        {
            sum += 2 * f(x);
        }
        else
        {
            sum += 4 * f(x);
        }
    }

    return sum * h / 3;
}

vector<pair<double, int>> Integrator::adaptiveIntegration(double (Integrator::*method)(const std::function<double(double)> &, const double &, const double &, int &), const std::function<double(double)> &f, const double &a, const double &b, int n)
{
    vector<pair<double, int>> result;
    double current;
    for (int i = 0; i < 4; i++)
    {
        current = (this->*method)(f, a, b, n);
        result.push_back(make_pair(current, n));
        n *= 2;
    }
    return result;
}

void Integrator::performIntegration(const IntegrationMethod &method, const function<double(double)> &f, const double &a, const double &b, const int &n)
{
    Integrator integrator;
    std::vector<std::pair<double, int>> results;

    switch (method)
    {
    case IntegrationMethod::Rectangles:
        results = integrator.adaptiveIntegration(Integrator::integrateRectangles, f, a, b, n);
        cout << "Метод прямоугольников:" << endl;
        break;

    case IntegrationMethod::Triangles:
        results = integrator.adaptiveIntegration(Integrator::integrateTriangles, f, a, b, n);
        cout << "Метод трапеций:" << endl;
        break;

    case IntegrationMethod::Parabolas:
        results = integrator.adaptiveIntegration(Integrator::integrateParabolas, f, a, b, n);
        cout << "Метод парабол (Симпсона):" << endl;
        break;
    }
    for (const auto &result : results)
    {
        cout << " Результат = " << result.first << ", n = " << result.second << endl;
    }
}
void Integrator::performGaussIntegration(const IntegrateGauss &method, const function<double(double)> &f, const double &a, const double &b, const int &n)
{
    int curr_n = n;
    Integrator integrator;
    vector<pair<double, int>> results;
    for (int i = 0; i < 4; i++)
    {
        switch (method)
        {
        case IntegrateGauss::Gauss_two_dots:
            results.push_back(make_pair(integrator.integrateGauss(f, a, b, n, gaussWeights[0], gaussNodes[0]), curr_n));
            break;
        case IntegrateGauss::Gauss_three_dots:
            results.push_back(make_pair(integrator.integrateGauss(f, a, b, n, gaussWeights[1], gaussNodes[1]), curr_n));
            break;

        case IntegrateGauss::Gauss_four_dots:
            results.push_back(make_pair(integrator.integrateGauss(f, a, b, n, gaussWeights[2], gaussNodes[2]), curr_n));
            break;
        }
        curr_n *= 2;
    }
    for (const auto &result : results)
    {
        cout << " Результат = " << result.first << ", n = " << result.second << endl;
    }
}