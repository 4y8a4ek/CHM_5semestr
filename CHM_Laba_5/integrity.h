#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <iomanip>

enum class IntegrationMethod
{
    Rectangles = 1,
    Triangles = 2,
    Parabolas = 3
};

enum class IntegrateGauss
{
    Gauss_two_dots = 1,
    Gauss_three_dots = 2,
    Gauss_four_dots = 3
};

using namespace std;

class Integrator
{
private:
    double gaussianQuadrature(const function<double(double)> &f, const double &a, const double &b, const vector<double> &gaussWeights, const vector<double> &gaussNodes);
    double integrateRectangles(const function<double(double)> &f, const double &a, const double &b, int &n);
    double integrateTriangles(const function<double(double)> &f, const double &a, const double &b, int &n);
    double integrateParabolas(const function<double(double)> &f, const double &a, const double &b, int &n);
    double integrateGauss(const function<double(double)> &f, const double &a, const double &b, const int &segments, const vector<double> &gaussWeights, const vector<double> &gaussNodes);
    vector<pair<double, int>> adaptiveIntegration(double (Integrator::*method)(const std::function<double(double)> &, const double &, const double &, int &), const std::function<double(double)> &f, const double &a, const double &b, int n);
public:
    void performIntegration(const IntegrationMethod &method, const function<double(double)> &f, const double &a, const double &b, const int &n);
    void performGaussIntegration(const IntegrateGauss &method, const function<double(double)> &f, const double &a, const double &b, const int &n);
};
