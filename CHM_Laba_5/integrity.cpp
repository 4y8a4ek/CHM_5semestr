#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include "integrity.h"

using namespace std;


double Integrator::gaussianQuadrature(const function<double(double)> &f, double &a, double &b, vector<double> &gaussWeights,vector<double> &gaussNodes)
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

double Integrator::integrateGauss(const function<double(double)> &f, double &a, double &b, int &segments, vector<double> &gaussWeights,vector<double> &gaussNodes)
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

double Integrator::integrateRectangles(const function<double(double)> &f, double &a, double &b, int &n)
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

double Integrator::integrateTriangles(const function<double(double)> &f, double &a, double &b, int &n)
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
double Integrator::integrateParabolas(const function<double(double)> &f, double &a, double &b, int &n)
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

vector<pair<double, int>> Integrator::adaptiveIntegration(double (Integrator::*method)(const std::function<double(double)> &, double &, double &, int &), const std::function<double(double)> &f, double &a, double &b, int n)
{
    vector<pair<double, int>> result;
    double current;
    for(int i = 0; i < 3; i++){
        current = (this->*method)(f, a, b, n);
        result.push_back(make_pair(current, n));
        n *= 2;
    }
    return result;
}
