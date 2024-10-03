#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include "integrity.cpp"

using namespace std;

double integrateRectangles(const function<double(double)> &f, double &a, double &b, int &n);
double integrateTriangles(const function<double(double)> &f, double &a, double &b, int &n);
double integrateParabolas(const function<double(double)> &f, double &a, double &b, int &n);
double gaussianQuadrature(const function<double(double)> &f, double &a, double &b);
double integrateGauss(const function<double(double)> &f, double &a, double &b, int &segments);
pair<double, int> adaptiveIntegration(const function<double(const function<double(double)> &, double &, double &, int &)> &method, const function<double(double)> &f, double &a, double &b, double &epsilon);
double hermite_polynomial(int n, double &x);
double hermite_polynomial_derivative(int &n, double &x);
vector<double> find_hermite_roots(int &n);
vector<double> compute_weights(int &n, const vector<double> &roots);
double integrate_with_hermite(int &n, const function<double(double)> &func);