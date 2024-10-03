#include <iostream>
#include <cmath>
#include <functional>
#include <vector>

using namespace std;

double integrateRectangles(const function<double(double)> &f, double &a, double &b, int &n)
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

double integrateTriangles(const function<double(double)> &f, double &a, double &b, int &n)
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

double integrateParabolas(const function<double(double)> &f, double &a, double &b, int &n)
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

vector<double> gaussWeights = {1.0, 1.0};
vector<double> gaussNodes = {-1 / sqrt(3), 1 / sqrt(3)};

double gaussianQuadrature(const function<double(double)> &f, double &a, double &b)
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

double integrateGauss(const function<double(double)> &f, double &a, double &b, int &segments)
{
    double segmentWidth = (b - a) / segments;
    double totalResult = 0.0;
    double segmentA, segmentB;
    for (int i = 0; i < segments; ++i)
    {
        segmentA = a + i * segmentWidth;
        segmentB = segmentA + segmentWidth;
        totalResult += gaussianQuadrature(f, segmentA, segmentB); // gaussianQuadrature для одного сегмента
    }

    return totalResult;
}

pair<double, int> adaptiveIntegration(const function<double(const function<double(double)> &, double &, double &, int &)> &method, const function<double(double)> &f, double &a, double &b, double &epsilon)
{
    int n = 2;
    double current = method(f, a, b, n);
    double previous = 0.0;

    do
    {
        previous = current;
        n *= 2;
        current = method(f, a, b, n);
    } while (abs(current - previous) > epsilon);

    return {current, n};
}

double hermite_polynomial(int n, double &x)
{
    if (n == 0)
        return 1.0;
    if (n == 1)
        return 2.0 * x;
    return 2.0 * x * hermite_polynomial(n - 1, x) - 2.0 * (n - 1) * hermite_polynomial(n - 2, x);
}

double hermite_polynomial_derivative(int &n, double &x)
{
    return 2.0 * n * hermite_polynomial(n - 1, x);
}

vector<double> find_hermite_roots(int &n)
{
    vector<double> roots(n);
    const double tolerance = 1e-10;
    for (int i = 0; i < n; ++i)
    {
        roots[i] = cos(M_PI * (i + 0.75) / (n + 0.5));
    }
    for (int i = 0; i < n; ++i)
    {
        double x = roots[i];
        double diff;
        do
        {
            double Hn = hermite_polynomial(n, x);
            double Hn_prime = hermite_polynomial_derivative(n, x);
            diff = Hn / Hn_prime;
            x = x - diff;
        } while (abs(diff) > tolerance);
        roots[i] = x;
    }
    for (int i = 0; i < n / 2; ++i)
    {
        roots[n - i - 1] = -roots[i];
    }
    return roots;
}

vector<double> compute_weights(int &n, const vector<double> &roots)
{
    vector<double> weights(n);
    for (int i = 0; i < n; ++i)
    {
        double xi = roots[i];
        double Hn_prime = hermite_polynomial_derivative(n, xi);
        weights[i] = 2.0 / (Hn_prime * Hn_prime);
    }
    double sum_weights = 0.0;
    for (double weight : weights)
    {
        sum_weights += weight;
    }
    for (double &weight : weights)
    {
        weight /= sum_weights;
    }
    return weights;
}

double integrate_with_hermite(int &n, const function<double(double)> &func)
{
    vector<double> roots = find_hermite_roots(n);
    vector<double> weights = compute_weights(n, roots);
    double integral = 0.0;
    for (int i = 0; i < n; ++i)
    {
        integral += weights[i] * func(roots[i]);
    }
    return integral * sqrt(M_PI);
}