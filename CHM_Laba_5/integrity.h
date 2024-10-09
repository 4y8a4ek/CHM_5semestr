#include <iostream>
#include <cmath>
#include <functional>
#include <vector>

using namespace std;

class Integrator
{
private:
    double gaussianQuadrature(const function<double(double)> &f, double &a, double &b, vector<double> &gaussWeights,vector<double> &gaussNodes);

public:
    double integrateRectangles(const function<double(double)> &f, double &a, double &b, int &n);
    double integrateTriangles(const function<double(double)> &f, double &a, double &b, int &n);
    double integrateParabolas(const function<double(double)> &f, double &a, double &b, int &n);
    double integrateGauss(const function<double(double)> &f, double &a, double &b, int &segments, vector<double> &gaussWeights,vector<double> &gaussNodes);
    vector<pair<double, int>> adaptiveIntegration( double(Integrator::*method)(const std::function<double(double)> &, double &, double &, int &), const std::function<double(double)> &f, double &a, double &b, int n);

};
